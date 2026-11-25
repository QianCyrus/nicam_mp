// =============================================================================
// bench_kokkos_full.cpp
//
// 100% Pure Portable Kokkos 4.x Benchmark & Validation Suite for NICAM NSW6
// (ZERO raw CUDA __global__ / __constant__ / __ballot_sync / CUDA Graphs)
//
// Best-practice portable Kokkos idioms used across ALL modes:
//   - Device-resident 0D Kokkos::View descriptors (dP, dS, dW) keep every
//     functor <= 64 bytes -> fast async launches, zero constant-memory sync!
//   - Explicit ExecSpace instance on Kokkos::deep_copy / RangePolicy / TeamPolicy
//     -> zero unnecessary host-device fences between stages!
//   - Portable Kokkos::TeamPolicy + PerTeam scratch compaction (256x fewer
//     global atomics; portable across CUDA, HIP, SYCL, and OpenMP).
//
// Same full-scale NICAM workload (10 regions x 16,641 cols x 80 levels =
// 13.31M cells, dt = 30 s, strict IEEE FP64) as bench_paper_ci.cu:
//   Mode 0 (Kokkos-BB)   : Ren-style unfused per-region Kokkos, original NSW6
//   Mode 1 (Kokkos-ENG)  : Idiomatic whole-domain Kokkos, original NSW6 math
//   Mode 2 (Kokkos-MATH) : SAME Kokkos framework + our 3 mathematical theorems
//   Mode 3 (Kokkos-MONO) : Mode 2 + O(K) monotone backward-envelope projection
// =============================================================================
#include <Kokkos_Core.hpp>
#include <Kokkos_Timer.hpp>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <utility>
#include <vector>

#include "nsw6_math_core.cuh"

using ExecSpace  = Kokkos::DefaultExecutionSpace;
using MemSpace   = ExecSpace::memory_space;
using ViewD1     = Kokkos::View<double*, MemSpace>;
using ViewU32    = Kokkos::View<uint32_t*, MemSpace>;
using ViewI32    = Kokkos::View<int*, MemSpace>;
using TeamPolicy = Kokkos::TeamPolicy<ExecSpace>;
using TeamMember = TeamPolicy::member_type;
using ScratchU32 = Kokkos::View<uint32_t*, ExecSpace::scratch_memory_space,
                                Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

struct StatePtrs {
    double *rhog, *rhogvx, *rhogvy, *rhogvz, *rhogw, *rhoge, *rhogq;
    double *rho, *pre, *tem, *vx, *vy, *vz, *w, *q;
    double *precip, *precip_lh, *precip_rhoein, *precip_rhophi, *precip_rhokin, *precip_energy;
    double *gdcfrc, *gdclw;
    const double *z, *zh, *phi, *gam2, *gam2h, *gsgam2, *gsgam2h;
    const double *ix, *iy, *iz, *jx, *jy, *jz, *unccn, *factor_n_aie;
    const double *rdgz, *rdgzh, *afact, *bfact, *cfact, *dfact;
};

struct WorkPtrs {
    double *Vt, *kv, *zdis_m, *rceff, *rctop;
    uint32_t *wet_cells, *wet_cols, *col_spmask, *wet_count, *wet_col_count;
    int *col_ktop;
};

KOKKOS_INLINE_FUNCTION uint32_t pack_cell(int ij, int k, int l) {
    return (uint32_t)ij | ((uint32_t)k << 18) | ((uint32_t)l << 26);
}
KOKKOS_INLINE_FUNCTION uint32_t pack_col(int ij, int l) {
    return (uint32_t)ij | ((uint32_t)l << 18);
}

static NSW6Params build_params(int ijdim, int kdim, int opt_mode, int compute_rceff) {
    NSW6Params P = {};
    P.ijdim = ijdim; P.kdim = kdim; P.kmin0 = 1; P.kmax0 = kdim - 2; P.nqmax = 6;
    P.i_qv = 0; P.i_qc = 1; P.i_qr = 2; P.i_qi = 3; P.i_qs = 4; P.i_qg = 5;
    P.opt_mode = opt_mode; P.ein_is_simple = 1; P.prcip_kin2ein = 1; P.precip_4water = 1;
    P.opt_indir = 0; P.opt_sclw = 0; P.opt_explicit_icegen = 1; P.compute_rceff = compute_rceff;
    P.is_top_tem = 1; P.is_top_epl = 0; P.is_btm_tem = 1; P.is_btm_epl = 0;
    P.mp_div_num = 1; P.collect_counters = 0;
    P.dt = 30.0; P.PI = 3.1415926535897932; P.EPS = 2.2204460492503131e-16;
    P.UNDEF = -9.9999e30; P.GRAV = 9.80665; P.Rdry = 287.05; P.Rvap = 461.51;
    P.CVdry = 1004.69 - 287.05; P.CL = 4217.8;
    P.LHV0 = 2.50084e6; P.LHS0 = 2.8345e6; P.LHF0 = P.LHS0 - P.LHV0;
    P.LHV = P.LHV0; P.LHS = P.LHS0; P.LHF = P.LHF0;
    P.PSAT0 = 610.7; P.TEM00 = 273.15; P.PRE00 = 101325.0; P.rho_w = 1000.0;
    P.CPovR_liq = 0.0; P.CVovR_liq = 0.0; P.LovR_liq = P.LHV / P.Rvap;
    P.CPovR_ice = 0.0; P.CVovR_ice = 0.0; P.LovR_ice = P.LHS / P.Rvap;
    for (int i = 0; i < 6; ++i) P.CVW[i] = P.CVdry;
    P.W_VADV1D_MIN = -40.0; P.dtemp_criteria = 1.0e-8;
    P.dens00 = 1.28; P.N0r = 8.0e6; P.N0s = 3.0e6; P.N0g = 4.0e8;
    P.Ar = P.PI * P.rho_w / 6.0; P.As = 0.069; P.Ag = P.PI * 400.0 / 6.0;
    P.Br = 3.0; P.Bs = 2.0; P.Bg = 3.0; P.Cr = 130.0; P.Cs = 4.84;
    P.Cg = std::sqrt((4.0 * 400.0 * P.GRAV) / (3.0 * P.dens00 * 0.6));
    P.Dr = 0.50; P.Ds = 0.25; P.Dg = 0.50;
    P.GAM = 1.0; P.GAM_2 = 1.0; P.GAM_3 = 2.0;
    P.GAM_1br = std::tgamma(1.0 + P.Br); P.GAM_2br = std::tgamma(2.0 + P.Br); P.GAM_3br = std::tgamma(3.0 + P.Br);
    P.GAM_3dr = std::tgamma(3.0 + P.Dr); P.GAM_6dr = std::tgamma(6.0 + P.Dr);
    P.GAM_1brdr = std::tgamma(1.0 + P.Br + P.Dr); P.GAM_5dr_h = std::tgamma(0.5 * (5.0 + P.Dr));
    P.GAM_1bs = std::tgamma(1.0 + P.Bs); P.GAM_2bs = std::tgamma(2.0 + P.Bs); P.GAM_3bs = std::tgamma(3.0 + P.Bs);
    P.GAM_3ds = std::tgamma(3.0 + P.Ds); P.GAM_1bsds = std::tgamma(1.0 + P.Bs + P.Ds);
    P.GAM_5ds_h = std::tgamma(0.5 * (5.0 + P.Ds));
    P.GAM_1bg = std::tgamma(1.0 + P.Bg); P.GAM_3dg = std::tgamma(3.0 + P.Dg);
    P.GAM_1bgdg = std::tgamma(1.0 + P.Bg + P.Dg); P.GAM_5dg_h = std::tgamma(0.5 * (5.0 + P.Dg));
    P.Eiw = 1.0; P.Erw = 1.0; P.Esw = 1.0; P.Egw = 1.0; P.Eri = 1.0; P.Esi = 0.25;
    P.Egi = 0.0; P.Esr = 1.0; P.Egr = 1.0; P.Egs = 0.0;
    P.gamma_sacr = 2.5e-2; P.gamma_gacs = 9.0e-2; P.mi = 4.19e-13; P.Nc_def = 50.0;
    P.beta_saut = 6.0e-3; P.gamma_saut = 6.0e-2; P.beta_gaut = 0.0; P.gamma_gaut = 9.0e-2;
    P.qicrt_saut = 0.0; P.qscrt_gaut = 6.0e-4;
    P.Ka0 = 2.428e-2; P.dKa_dT = 7.47e-5; P.Kd0 = 2.222e-5; P.dKd_dT = 1.37e-7;
    P.nu0 = 1.718e-5; P.dnu_dT = 5.28e-8;
    P.f1r = 0.78; P.f2r = 0.27; P.f1s = 0.65; P.f2s = 0.39; P.f1g = 0.78; P.f2g = 0.27;
    P.A_frz = 0.66; P.B_frz = 100.0; P.mi40 = 2.46e-10; P.mi50 = 4.80e-10; P.vti50 = 1.0; P.Ri50 = 5.0e-5;
    P.Di_max = 500.0e-6; P.Di_a = 11.9;
    P.sw_kk2000 = 0.0; P.sw_roh2014 = 1.0; P.sw_expice = 1.0; P.sw_sr2020 = 0.0;
    P.sw_phillips2007 = 0.0; P.sw_beta_saut = 0.0; P.sw_constVti = 0.0; P.constVti = 0.0;
    P.Ds_cutoff = 0.0; P.Dg_cutoff = 0.0; P.CRIT_SCLW_W = 0.01;
    double ca[10] = { 5.065339, -0.062659, -3.032362, 0.029469, -0.000285,
                      0.31255,   0.000204,  0.003199, 0.0,      -0.015952 };
    double cb[10] = { 0.476221, -0.015896,  0.165977, 0.007468, -0.000141,
                      0.060366,  0.000079,  0.000594, 0.0,      -0.003577 };
    for (int i = 0; i < 10; ++i) { P.coef_a[i] = ca[i]; P.coef_b[i] = cb[i]; }
    P.ln10 = std::log(10.0);
    double Dc = 0.146 - 5.964e-2 * std::log(P.Nc_def / 2000.0);
    double dg = 0.1444 / (Dc * Dc);
    P.rceff_coef_dgam = std::tgamma(dg + 2.0/3.0) / std::tgamma(dg) * std::pow(dg, -2.0/3.0);
    P.rceff_coef_xf   = 3.0 / (4.0 * P.PI * P.rho_w);
    return P;
}

struct HostState {
    int ijdim, kdim, lall, nqmax;
    std::vector<double> gz, gzh, dgz, rdgz, rdgzh, afact, bfact, cfact, dfact;
    std::vector<double> z, zh, phi, gam2, gam2h, gsgam2, gsgam2h;
    std::vector<double> ix, iy, iz, jx, jy, jz, unccn, factor_n_aie;
    std::vector<double> rhog, rhogvx, rhogvy, rhogvz, rhogw, rhoge, rhogq;
    std::vector<double> rho, pre, tem, vx, vy, vz, w, q;
    std::vector<double> precip, precip_lh, precip_rhoein, precip_rhophi, precip_rhokin, precip_energy;
    std::vector<double> gdcfrc, gdclw;

    void init(int ij_, int k_, int l_, const NSW6Params& P) {
        ijdim = ij_; kdim = k_; lall = l_; nqmax = 6;
        size_t n3 = (size_t)ijdim * kdim * lall, n4 = n3 * nqmax, n2 = (size_t)ijdim * lall;
        gz.resize(kdim); gzh.resize(kdim); dgz.resize(kdim); rdgz.resize(kdim); rdgzh.resize(kdim);
        afact.resize(kdim); bfact.resize(kdim); cfact.resize(kdim); dfact.resize(kdim);
        gzh[0] = -33.0; gzh[1] = 0.0;
        for (int k = 1; k < kdim - 1; ++k) gzh[k + 1] = gzh[k] + 66.0 + 1350.0 * std::pow((double)(k - 1) / (double)(kdim - 3), 2.2);
        for (int k = 1; k < kdim - 1; ++k) gz[k] = 0.5 * (gzh[k] + gzh[k + 1]);
        gz[0] = -gz[1]; gz[kdim - 1] = 2.0 * gzh[kdim - 1] - gz[kdim - 2];
        for (int k = 0; k < kdim; ++k) {
            dgz[k] = (k < kdim - 1) ? (gzh[k + 1] - gzh[k]) : (gzh[kdim - 1] - gzh[kdim - 2]);
            rdgz[k] = 1.0 / dgz[k]; rdgzh[k] = 1.0 / ((k > 0) ? (gz[k] - gz[k - 1]) : (gz[1] - gz[0]));
            afact[k] = (k > 0) ? ((gzh[k] - gz[k - 1]) / (gz[k] - gz[k - 1])) : 0.5; bfact[k] = 1.0 - afact[k];
            cfact[k] = (k < kdim - 1) ? ((gz[k] - gzh[k]) / (gzh[k + 1] - gzh[k])) : 0.5; dfact[k] = 1.0 - cfact[k];
        }
        z.resize(n3); zh.resize(n3); phi.resize(n3); gam2.resize(n3); gam2h.resize(n3); gsgam2.resize(n3); gsgam2h.resize(n3);
        ix.resize(n2); iy.resize(n2); iz.resize(n2); jx.resize(n2); jy.resize(n2); jz.resize(n2);
        unccn.resize(n3); factor_n_aie.resize(n2);
        rhog.resize(n3); rhogvx.resize(n3); rhogvy.resize(n3); rhogvz.resize(n3); rhogw.resize(n3); rhoge.resize(n3); rhogq.resize(n4);
        rho.resize(n3); pre.resize(n3); tem.resize(n3); vx.resize(n3); vy.resize(n3); vz.resize(n3); w.resize(n3); q.resize(n4);
        precip.resize(2 * n2, 0.0); precip_lh.resize(n2, 0.0); precip_rhoein.resize(n2, 0.0);
        precip_rhophi.resize(n2, 0.0); precip_rhokin.resize(n2, 0.0); precip_energy.resize(n2, 0.0);
        gdcfrc.resize(n3, 0.0); gdclw.resize(n3, 0.0);

        double ztop = gzh[kdim - 1]; size_t s3 = (size_t)ijdim * kdim;
        for (int ll = 0; ll < lall; ++ll) {
            for (int ij = 0; ij < ijdim; ++ij) {
                int patch = (ij / 64 + ll * 17) % 10;
                double lon = 2.0 * P.PI * (double)(ij + ll * ijdim) / (double)(ijdim * lall);
                double lat = -1.3 + 2.6 * (double)((ij * 13 + ll * 7) % ijdim) / (double)ijdim;
                size_t i2 = (size_t)ll * ijdim + ij;
                ix[i2] = -std::sin(lon); iy[i2] = std::cos(lon); iz[i2] = 0.0;
                jx[i2] = -std::sin(lat) * std::cos(lon); jy[i2] = -std::sin(lat) * std::sin(lon); jz[i2] = std::cos(lat);
                factor_n_aie[i2] = 0.8 + 0.4 * (double)(ij % 7) / 6.0;
                double topo = 800.0 * (0.5 + 0.5 * std::sin(3.0 * lon) * std::cos(2.0 * lat));
                for (int k = 0; k < kdim; ++k) {
                    size_t i3 = (size_t)ll * s3 + (size_t)k * ijdim + ij;
                    size_t i4 = (size_t)ll * nqmax * s3 + (size_t)k * ijdim + ij;
                    double zk = topo + (1.0 - topo / ztop) * gz[k], zhk = topo + (1.0 - topo / ztop) * gzh[k];
                    z[i3] = zk; zh[i3] = zhk; phi[i3] = P.GRAV * zk;
                    gam2[i3] = 1.0 + 0.002 * zk / ztop; gam2h[i3] = 1.0 + 0.002 * zhk / ztop;
                    gsgam2[i3] = (1.0 - topo / ztop) * gam2[i3]; gsgam2h[i3] = (1.0 - topo / ztop) * gam2h[i3];
                    unccn[i3] = 1.0e8;
                    double alt = std::max(zk, 0.0);
                    double Tk = (patch == 7 || patch == 8) ? std::max(200.0, 282.0 - 6.5e-3 * alt) : std::max(195.0, 299.0 - 6.5e-3 * alt);
                    double pk = 101325.0 * std::exp(-alt / 8200.0), rk = pk / (P.Rdry * Tk);
                    double qsl = P.PSAT0 * std::exp(P.LovR_liq * (1.0/P.TEM00 - 1.0/Tk)) / (rk * P.Rvap * Tk);
                    double qsi = P.PSAT0 * std::exp(P.LovR_ice * (1.0/P.TEM00 - 1.0/Tk)) / (rk * P.Rvap * Tk);
                    double qvk = 0.4 * std::min(qsl, qsi), qck = 0.0, qrk = 0.0, qik = 0.0, qsk = 0.0, qgk = 0.0;
                    if (k >= 1 && k <= kdim - 2) {
                        if (patch == 6 && alt < 12000.0) {
                            qvk = 1.12 * qsl;
                            if (Tk > 265.0) qck = 1.6e-3;
                            if (Tk > 270.0) qrk = 3.8e-3 * std::exp(-alt / 5000.0);
                            if (Tk < 275.0) qgk = 2.2e-3 * std::exp(-std::fabs(alt - 6500.0) / 2500.0);
                            if (Tk < 273.0) qsk = 1.1e-3 * std::exp(-std::fabs(alt - 8000.0) / 3000.0);
                            if (Tk < 260.0) qik = 3.5e-4;
                        } else if (patch == 7 && alt < 10000.0) {
                            qvk = 0.96 * std::min(qsl, qsi);
                            if (Tk < 275.0 && Tk > 240.0) qsk = 1.4e-3;
                            if (Tk > 271.0) qrk = 7.0e-4;
                            if (Tk < 270.0) qck = 2.5e-4;
                        } else if (patch == 8 && alt > 6500.0 && alt < 14500.0) {
                            qvk = 1.18 * qsi; qik = 2.2e-4;
                            if (Tk > 240.0) { qck = 6.0e-4; qsk = 4.0e-4; }
                        } else if (patch == 9) {
                            if (k >= 11 && k <= 21) { qrk = 5.5e-3; qgk = 3.5e-3; qsk = 1.6e-3; qvk = 1.05 * qsl; qck = 1.0e-3; }
                            else if (k >= 3 && k <= 7) { qrk = 8.0e-5; qvk = 0.9 * qsl; }
                        }
                    }
                    double qdk = 1.0 - (qvk + qck + qrk + qik + qsk + qgk);
                    rho[i3] = rk; tem[i3] = Tk; pre[i3] = rk * (qdk * P.Rdry + qvk * P.Rvap) * Tk;
                    vx[i3] = 14.0 * std::cos(lon) + 4.0 * std::sin(alt / 3000.0);
                    vy[i3] = -11.0 * std::sin(lon) + 3.0 * std::cos(alt / 3000.0);
                    vz[i3] = 2.5 * std::sin(lat);
                    w[i3]  = (k == 0 || k == kdim - 1) ? 0.0 : ((patch == 6 || patch == 9) ? 3.0 * std::sin(P.PI * alt / 15000.0) : 0.15);
                    rhog[i3] = rk * gsgam2[i3];
                    rhogvx[i3] = rhog[i3] * vx[i3]; rhogvy[i3] = rhog[i3] * vy[i3]; rhogvz[i3] = rhog[i3] * vz[i3];
                    rhoge[i3] = rhog[i3] * Tk * P.CVdry;
                    q[i4 + 0*s3] = qvk; q[i4 + 1*s3] = qck; q[i4 + 2*s3] = qrk;
                    q[i4 + 3*s3] = qik; q[i4 + 4*s3] = qsk; q[i4 + 5*s3] = qgk;
                    for (int nq = 0; nq < 6; ++nq) rhogq[i4 + nq*s3] = rhog[i3] * q[i4 + nq*s3];
                }
                for (int k = 0; k < kdim; ++k) {
                    size_t i3 = (size_t)ll * s3 + (size_t)k * ijdim + ij;
                    double rhogh_k = (k == 0)
                        ? ((afact[1] * rho[i3 + ijdim] + bfact[1] * rho[i3]) * gsgam2h[i3 + ijdim])
                        : ((afact[k] * rho[i3] + bfact[k] * rho[i3 - ijdim]) * gsgam2h[i3]);
                    rhogw[i3] = rhogh_k * w[i3];
                }
            }
        }
    }
};

struct KokkosState {
    ViewD1 z, zh, phi, gam2, gam2h, gsgam2, gsgam2h;
    ViewD1 ix, iy, iz, jx, jy, jz, unccn, factor_n_aie;
    ViewD1 rdgz, rdgzh, afact, bfact, cfact, dfact;
    ViewD1 rhog, rhogvx, rhogvy, rhogvz, rhogw, rhoge, rhogq;
    ViewD1 rho, pre, tem, vx, vy, vz, w, q;
    ViewD1 init_rhog, init_rhogvx, init_rhogvy, init_rhogvz, init_rhogw, init_rhoge, init_rhogq;
    ViewD1 init_rho, init_pre, init_tem, init_vx, init_vy, init_vz, init_w, init_q;
    ViewD1 precip, precip_lh, precip_rhoein, precip_rhophi, precip_rhokin, precip_energy;
    ViewD1 gdcfrc, gdclw;
    Kokkos::View<StatePtrs, MemSpace> d_S;
    const StatePtrs* dp;

    static ViewD1 upload(const std::string& name, const std::vector<double>& v) {
        ViewD1 d(name, v.size());
        auto h = Kokkos::create_mirror_view(d);
        std::copy(v.begin(), v.end(), h.data());
        Kokkos::deep_copy(d, h);
        return d;
    }

    void init(const HostState& H) {
        z = upload("z", H.z); zh = upload("zh", H.zh); phi = upload("phi", H.phi);
        gam2 = upload("gam2", H.gam2); gam2h = upload("gam2h", H.gam2h);
        gsgam2 = upload("gsgam2", H.gsgam2); gsgam2h = upload("gsgam2h", H.gsgam2h);
        ix = upload("ix", H.ix); iy = upload("iy", H.iy); iz = upload("iz", H.iz);
        jx = upload("jx", H.jx); jy = upload("jy", H.jy); jz = upload("jz", H.jz);
        unccn = upload("unccn", H.unccn); factor_n_aie = upload("faie", H.factor_n_aie);
        rdgz = upload("rdgz", H.rdgz); rdgzh = upload("rdgzh", H.rdgzh);
        afact = upload("afact", H.afact); bfact = upload("bfact", H.bfact);
        cfact = upload("cfact", H.cfact); dfact = upload("dfact", H.dfact);

        init_rhog = upload("irhog", H.rhog); init_rhogvx = upload("irhogvx", H.rhogvx);
        init_rhogvy = upload("irhogvy", H.rhogvy); init_rhogvz = upload("irhogvz", H.rhogvz);
        init_rhogw = upload("irhogw", H.rhogw); init_rhoge = upload("irhoge", H.rhoge);
        init_rhogq = upload("irhogq", H.rhogq);
        init_rho = upload("irho", H.rho); init_pre = upload("ipre", H.pre); init_tem = upload("item", H.tem);
        init_vx = upload("ivx", H.vx); init_vy = upload("ivy", H.vy); init_vz = upload("ivz", H.vz);
        init_w = upload("iw", H.w); init_q = upload("iq", H.q);

        rhog = upload("rhog", H.rhog); rhogvx = upload("rhogvx", H.rhogvx);
        rhogvy = upload("rhogvy", H.rhogvy); rhogvz = upload("rhogvz", H.rhogvz);
        rhogw = upload("rhogw", H.rhogw); rhoge = upload("rhoge", H.rhoge);
        rhogq = upload("rhogq", H.rhogq);
        rho = upload("rho", H.rho); pre = upload("pre", H.pre); tem = upload("tem", H.tem);
        vx = upload("vx", H.vx); vy = upload("vy", H.vy); vz = upload("vz", H.vz);
        w = upload("w", H.w); q = upload("q", H.q);

        precip = upload("precip", H.precip); precip_lh = upload("precip_lh", H.precip_lh);
        precip_rhoein = upload("precip_re", H.precip_rhoein); precip_rhophi = upload("precip_rp", H.precip_rhophi);
        precip_rhokin = upload("precip_rk", H.precip_rhokin); precip_energy = upload("precip_en", H.precip_energy);
        gdcfrc = upload("gdcfrc", H.gdcfrc); gdclw = upload("gdclw", H.gdclw);

        StatePtrs sp = {
            rhog.data(), rhogvx.data(), rhogvy.data(), rhogvz.data(), rhogw.data(), rhoge.data(), rhogq.data(),
            rho.data(), pre.data(), tem.data(), vx.data(), vy.data(), vz.data(), w.data(), q.data(),
            precip.data(), precip_lh.data(), precip_rhoein.data(), precip_rhophi.data(), precip_rhokin.data(), precip_energy.data(),
            gdcfrc.data(), gdclw.data(),
            z.data(), zh.data(), phi.data(), gam2.data(), gam2h.data(), gsgam2.data(), gsgam2h.data(),
            ix.data(), iy.data(), iz.data(), jx.data(), jy.data(), jz.data(), unccn.data(), factor_n_aie.data(),
            rdgz.data(), rdgzh.data(), afact.data(), bfact.data(), cfact.data(), dfact.data()
        };
        d_S = Kokkos::View<StatePtrs, MemSpace>("d_S");
        auto h_S = Kokkos::create_mirror_view(d_S);
        h_S() = sp;
        Kokkos::deep_copy(d_S, h_S);
        dp = d_S.data();
    }

    void restore(const ExecSpace& exec) {
        Kokkos::deep_copy(exec, rhog, init_rhog); Kokkos::deep_copy(exec, rhogvx, init_rhogvx);
        Kokkos::deep_copy(exec, rhogvy, init_rhogvy); Kokkos::deep_copy(exec, rhogvz, init_rhogvz);
        Kokkos::deep_copy(exec, rhogw, init_rhogw); Kokkos::deep_copy(exec, rhoge, init_rhoge);
        Kokkos::deep_copy(exec, rhogq, init_rhogq);
        Kokkos::deep_copy(exec, rho, init_rho); Kokkos::deep_copy(exec, pre, init_pre); Kokkos::deep_copy(exec, tem, init_tem);
        Kokkos::deep_copy(exec, vx, init_vx); Kokkos::deep_copy(exec, vy, init_vy); Kokkos::deep_copy(exec, vz, init_vz);
        Kokkos::deep_copy(exec, w, init_w); Kokkos::deep_copy(exec, q, init_q);
        exec.fence();
    }
};

struct KokkosWorkspace {
    ViewD1 Vt, kv, zdis_m, rceff, rctop;
    ViewU32 wet_cells, wet_cols, col_spmask, wet_count, wet_col_count;
    ViewI32 col_ktop;
    Kokkos::View<WorkPtrs, MemSpace> d_W;
    const WorkPtrs* dp;

    void init(int ijdim, int kdim, int lall) {
        size_t n3 = (size_t)ijdim * kdim * lall;
        size_t n2 = (size_t)ijdim * lall;
        Vt            = ViewD1("Vt", 4 * n3);
        kv            = ViewD1("kv", n3);
        zdis_m        = ViewD1("zdis_m", 4 * n3);
        rceff         = ViewD1("rceff", n3);
        rctop         = ViewD1("rctop", n2);
        wet_cells     = ViewU32("wet_cells", n3);
        wet_cols      = ViewU32("wet_cols", n2);
        col_spmask    = ViewU32("col_spmask", n2);
        col_ktop      = ViewI32("col_ktop", n2);
        wet_count     = ViewU32("wet_count", 1);
        wet_col_count = ViewU32("wet_col_count", 1);
        WorkPtrs wp = { Vt.data(), kv.data(), zdis_m.data(), rceff.data(), rctop.data(),
                        wet_cells.data(), wet_cols.data(), col_spmask.data(),
                        wet_count.data(), wet_col_count.data(), col_ktop.data() };
        d_W = Kokkos::View<WorkPtrs, MemSpace>("d_W");
        auto h_W = Kokkos::create_mirror_view(d_W);
        h_W() = wp;
        Kokkos::deep_copy(d_W, h_W);
        dp = d_W.data();
    }
};

// =============================================================================
// Shared Pure-Kokkos Helpers (tiny <= 40-byte functor captures!)
// =============================================================================
static void kokkos_init_kv(const ExecSpace& exec, const NSW6Params* dP, int ijdim, int kdim, int lall,
                           const StatePtrs* dS, const WorkPtrs* dW) {
    const size_t s3 = (size_t)ijdim * kdim, total = (size_t)lall * s3;
    Kokkos::parallel_for("init_kv", Kokkos::RangePolicy<ExecSpace>(exec, 0, total), KOKKOS_LAMBDA(const size_t gid) {
        const NSW6Params& P = *dP; const StatePtrs& S = *dS; const WorkPtrs& W = *dW;
        const int l = (int)(gid / s3);
        const size_t idx = gid - (size_t)l * s3;
        const int k = (int)(idx / (size_t)ijdim);
        const size_t off3 = (size_t)l * s3;
        if (k <= P.kmin0 || k > P.kmax0) { W.kv[off3 + idx] = 0.0; return; }
        const double rw = S.rhogw[off3 + idx];
        if (rw == 0.0) { W.kv[off3 + idx] = 0.0; return; }
        const double gs2h   = S.gsgam2h[off3 + idx];
        const double rho_k  = S.rhog[off3 + idx] / S.gsgam2[off3 + idx];
        const double rho_m1 = S.rhog[off3 + idx - ijdim] / S.gsgam2[off3 + idx - ijdim];
        const double rhogh  = (S.afact[k] * rho_k + S.bfact[k] * rho_m1) * gs2h;
        W.kv[off3 + idx] = 0.5 * rw * rw / rhogh;
    });
}

static void kokkos_rceff_whole(const ExecSpace& exec, const NSW6Params* dP, int ijdim, int kdim, int lall,
                               const StatePtrs* dS, const WorkPtrs* dW) {
    const size_t s3 = (size_t)ijdim * kdim, n2 = (size_t)lall * ijdim;
    Kokkos::parallel_for("rceff", Kokkos::RangePolicy<ExecSpace>(exec, 0, n2), KOKKOS_LAMBDA(const size_t col) {
        const NSW6Params& P = *dP; const StatePtrs& S = *dS; const WorkPtrs& W = *dW;
        const int l  = (int)(col / (size_t)ijdim);
        const int ij = (int)(col - (size_t)l * ijdim);
        const size_t off3 = (size_t)l * s3, off4 = (size_t)l * (size_t)P.nqmax * s3, off2 = (size_t)l * (size_t)ijdim;
        rceff_column(P, ij, S.rho + off3, S.tem + off3, S.q + off4,
                     S.unccn + off3, S.factor_n_aie + off2,
                     W.rceff + off3, W.rctop + off2);
    });
}

static void kokkos_diag_w_and_bndcnd(const ExecSpace& exec, const NSW6Params* dP, int ijdim, int kdim, int lall,
                                     const StatePtrs* dS) {
    const size_t s3 = (size_t)ijdim * kdim, total = (size_t)lall * s3, n2 = (size_t)lall * ijdim;
    Kokkos::parallel_for("diag_w", Kokkos::RangePolicy<ExecSpace>(exec, 0, total), KOKKOS_LAMBDA(const size_t gid) {
        const NSW6Params& P = *dP; const StatePtrs& S = *dS;
        const int l = (int)(gid / s3);
        const size_t idx = gid - (size_t)l * s3;
        const int k = (int)(idx / (size_t)ijdim);
        const size_t off3 = (size_t)l * s3;
        if (k <= P.kmin0 || k > P.kmax0) { S.w[off3 + idx] = 0.0; return; }
        double rhogh = (S.afact[k] * S.rho[off3 + idx] + S.bfact[k] * S.rho[off3 + idx - ijdim]) * S.gsgam2h[off3 + idx];
        S.w[off3 + idx] = S.rhogw[off3 + idx] / rhogh;
    });

    Kokkos::parallel_for("bndcnd", Kokkos::RangePolicy<ExecSpace>(exec, 0, n2), KOKKOS_LAMBDA(const size_t col) {
        const NSW6Params& P = *dP; const StatePtrs& S = *dS;
        const int l = (int)(col / (size_t)ijdim), ij = (int)(col - (size_t)l * ijdim);
        const size_t off3 = (size_t)l * s3;
        const size_t idx_b  = (size_t)ij + (size_t)(P.kmin0 - 1) * ijdim;
        const size_t idx_b1 = (size_t)ij + (size_t)P.kmin0 * ijdim;
        const size_t idx_t  = (size_t)ij + (size_t)(P.kmax0 + 1) * ijdim;
        const size_t idx_t1 = (size_t)ij + (size_t)P.kmax0 * ijdim;

        double tb = S.tem[off3 + idx_b1];
        double pb = S.pre[off3 + idx_b1] + S.rho[off3 + idx_b1] * (S.phi[off3 + idx_b1] - S.phi[off3 + idx_b]);
        S.tem[off3 + idx_b] = tb; S.pre[off3 + idx_b] = pb; S.rho[off3 + idx_b] = pb / (P.Rdry * tb);
        S.vx[off3 + idx_b] = S.vx[off3 + idx_b1]; S.vy[off3 + idx_b] = S.vy[off3 + idx_b1]; S.vz[off3 + idx_b] = S.vz[off3 + idx_b1];
        S.w[off3 + idx_b] = 0.0;

        double tt = S.tem[off3 + idx_t1];
        double pt = S.pre[off3 + idx_t1] - S.rho[off3 + idx_t1] * (S.phi[off3 + idx_t] - S.phi[off3 + idx_t1]);
        if (pt <= 0.0) pt = S.pre[off3 + idx_t1];
        S.tem[off3 + idx_t] = tt; S.pre[off3 + idx_t] = pt; S.rho[off3 + idx_t] = pt / (P.Rdry * tt);
        S.vx[off3 + idx_t] = S.vx[off3 + idx_t1]; S.vy[off3 + idx_t] = S.vy[off3 + idx_t1]; S.vz[off3 + idx_t] = S.vz[off3 + idx_t1];
        S.w[off3 + idx_t] = 0.0;
    });
}

// =============================================================================
// MODE 0 (Kokkos-BB): Unfused per-region Kokkos, original NSW6 formulas
// =============================================================================
static void run_kokkos_mode0(const ExecSpace& exec, const NSW6Params* dP, int ijdim, int kdim, int lall, int do_rc,
                             const StatePtrs* dS, const WorkPtrs* dW) {
    const size_t s3 = (size_t)ijdim * kdim;
    const int kmin0 = 1, kmax0 = kdim - 2;

    for (int l = 0; l < lall; ++l) {
        const size_t off3 = (size_t)l * s3, off4 = (size_t)l * 6ULL * s3;
        const size_t off2 = (size_t)l * (size_t)ijdim, offv = (size_t)l * 4ULL * s3;
        const size_t offp = (size_t)l * 2ULL * (size_t)ijdim;

        Kokkos::parallel_for("m0_nf", Kokkos::RangePolicy<ExecSpace>(exec, 0, s3), KOKKOS_LAMBDA(const size_t idx) {
            const NSW6Params& P = *dP; const StatePtrs& S = *dS;
            int k = (int)(idx / (size_t)ijdim);
            if (k < P.kmin0 || k > P.kmax0) return;
            double rg = S.rhog[off3 + idx], temp = S.tem[off3 + idx], gs2 = S.gsgam2[off3 + idx];
            double rq[6], diffq = 0.0;
            for (int n = 0; n < 6; ++n) rq[n] = S.rhogq[off4 + idx + n * s3];
            for (int n = 1; n < 6; ++n) { double vc = (rq[n] > 0.0) ? rq[n] : 0.0; diffq += rq[n] - vc; rq[n] = vc; }
            rq[0] += diffq;
            double rc0 = (rq[0] > 0.0) ? rq[0] : 0.0; double diff0 = rq[0] - rc0; rq[0] = rc0;
            rg = rg * (1.0 - diff0);
            double rho_v = rg / gs2;
            double qv = rq[0]/rg, qc = rq[1]/rg, qr = rq[2]/rg, qi = rq[3]/rg, qs = rq[4]/rg, qg = rq[5]/rg;
            double qd = 1.0 - (qv + qc + qr + qi + qs + qg);
            S.rhog[off3 + idx] = rg;
            for (int n = 0; n < 6; ++n) { S.rhogq[off4 + idx + n * s3] = rq[n]; S.q[off4 + idx + n * s3] = rq[n] / rg; }
            S.rhoge[off3 + idx] = temp * rg * P.CVdry;
            S.rho[off3 + idx] = rho_v; S.pre[off3 + idx] = rho_v * (qd * P.Rdry + qv * P.Rvap) * temp;
        });

        Kokkos::parallel_for("m0_nsw6", Kokkos::RangePolicy<ExecSpace>(exec, 0, s3), KOKKOS_LAMBDA(const size_t idx) {
            const NSW6Params& P = *dP; const StatePtrs& S = *dS; const WorkPtrs& W = *dW;
            int k = (int)(idx / (size_t)ijdim), ij = (int)(idx - (size_t)k * ijdim);
            int d0=0, d1=0, d2=0, d3=0, d4=0;
            nsw6_pre_and_core_cell(P, ij, k, S.rhog+off3, S.rhoge+off3, S.rhogq+off4,
                S.rho+off3, S.tem+off3, S.pre+off3, S.q+off4, S.w+off3,
                S.unccn+off3, S.factor_n_aie+off2, S.gsgam2+off3, W.Vt+offv,
                &d0,&d1,&d2,&d3,&d4);
        });

        if (do_rc) {
            for (int k = kmin0; k <= kmax0; ++k) {
                Kokkos::parallel_for("m0_rceff_k", Kokkos::RangePolicy<ExecSpace>(exec, 0, ijdim), KOKKOS_LAMBDA(const int ij) {
                    const NSW6Params& P = *dP; const StatePtrs& S = *dS; const WorkPtrs& W = *dW;
                    size_t idx = (size_t)ij + (size_t)k * ijdim;
                    double qc = S.q[off4 + idx + P.i_qc * s3], r_v = 0.0;
                    if (qc >= 1.0e-5) {
                        double Nc_v = P.Nc_def * S.factor_n_aie[off2 + ij];
                        double xf = S.rho[off3 + idx] * qc / (Nc_v * 1.0e6);
                        r_v = P.rceff_coef_dgam * pow(P.rceff_coef_xf * xf + P.EPS, 1.0 / 3.0);
                    }
                    W.rceff[off3 + idx] = r_v;
                });
            }
        }

        Kokkos::parallel_for("m0_precip_h", Kokkos::RangePolicy<ExecSpace>(exec, 0, ijdim), KOKKOS_LAMBDA(const int ij) {
            const NSW6Params& P = *dP; const StatePtrs& S = *dS; const WorkPtrs& W = *dW;
            const int sp_nq[4] = { P.i_qr, P.i_qi, P.i_qs, P.i_qg };
            const double dt = P.dt;
            double *rhogw = S.rhogw + off3, *kv = W.kv + off3;
            const double *rhog = S.rhog + off3, *q = S.q + off4, *z = S.z + off3, *zh = S.zh + off3;
            const double *gs2 = S.gsgam2 + off3, *gs2h = S.gsgam2h + off3, *g2h = S.gam2h + off3, *Vt = W.Vt + offv;
            for (int k = 0; k < kdim; ++k) {
                size_t idx = (size_t)ij + (size_t)k * ijdim;
                if (k <= P.kmin0 || k > P.kmax0) { kv[idx] = 0.0; continue; }
                double rhogh = (S.afact[k]*(rhog[idx]/gs2[idx]) + S.bfact[k]*(rhog[idx-ijdim]/gs2[idx-ijdim])) * gs2h[idx];
                double rw = rhogw[idx];
                kv[idx] = 0.5 * rw * rw / rhogh;
            }
            const double dzh_kmin = z[ij + P.kmin0 * ijdim] - zh[ij + P.kmin0 * ijdim];
            double pfk = 0.0, pfw = 0.0;
            for (int k = P.kmin0 + 1; k <= P.kmax0; ++k) {
                size_t kdx = (size_t)ij + (size_t)k * ijdim;
                double zk = z[kdx], zkm1 = z[kdx - ijdim], dzh = zk - zkm1;
                double dzhm1 = (k - 1 == P.kmin0) ? dzh_kmin : (zkm1 - z[kdx - 2ULL * ijdim]);
                double fkv_k = 0.0, fw_k = 0.0;
                for (int s = 0; s < 4; ++s) {
                    double wk = Vt[kdx + s * s3]; if (wk == 0.0) continue;
                    double zdis = vadv1d_taylor_zdis(wk, Vt[kdx + ijdim + s * s3], Vt[kdx - ijdim + s * s3],
                                                     dzh, dzhm1, dzh + dzhm1, dt, k == P.kmin0 + 1);
                    if (zdis <= 0.0) continue;
                    double ztmp = zkm1 + zdis; int c = k; double zlo = zkm1;
                    for (int j = k; j <= P.kmax0; ++j) {
                        double zj = z[(size_t)ij + (size_t)j * ijdim], zjm1 = z[(size_t)ij + (size_t)(j - 1) * ijdim];
                        if (ztmp > zjm1 && ztmp <= zj) { c = j; zlo = zjm1; break; }
                    }
                    double rem = zdis - (zlo - zkm1);
                    for (int j = k; j <= c; ++j) {
                        size_t jdx = (size_t)ij + (size_t)j * ijdim;
                        double L = (j < c) ? (z[jdx] - z[jdx - ijdim]) : rem;
                        double qhj = (S.afact[j]*q[jdx + sp_nq[s]*s3] + S.bfact[j]*q[jdx - ijdim + sp_nq[s]*s3]) * (g2h[jdx]/gs2h[jdx]);
                        fkv_k -= L * qhj * kv[jdx];
                        fw_k  -= L * qhj * rhogw[jdx];
                    }
                }
                kv[kdx]    -= (fkv_k - pfk) * S.rdgzh[k];
                rhogw[kdx] -= (fw_k  - pfw) * S.rdgzh[k];
                pfk = fkv_k; pfw = fw_k;
            }
        });

        Kokkos::parallel_for("m0_precip_f", Kokkos::RangePolicy<ExecSpace>(exec, 0, ijdim), KOKKOS_LAMBDA(const int ij) {
            const NSW6Params& P = *dP; const StatePtrs& S = *dS; const WorkPtrs& W = *dW;
            const int sp_nq[4] = { P.i_qr, P.i_qi, P.i_qs, P.i_qg };
            const double dt = P.dt, inv_dt = 1.0 / dt;
            double *rhog = S.rhog+off3, *rhogvx = S.rhogvx+off3, *rhogvy = S.rhogvy+off3, *rhogvz = S.rhogvz+off3;
            double *rhoge = S.rhoge+off3, *rhogq = S.rhogq+off4;
            const double *rhogw = S.rhogw+off3, *tem = S.tem+off3, *vx = S.vx+off3, *vy = S.vy+off3, *vz = S.vz+off3;
            const double *z = S.z+off3, *zh = S.zh+off3, *gs2 = S.gsgam2+off3, *gs2h = S.gsgam2h+off3, *g2 = S.gam2+off3;
            const double *Vt = W.Vt+offv, *kv = W.kv+off3;
            const double ix_c = S.ix[off2+ij], iy_c = S.iy[off2+ij], iz_c = S.iz[off2+ij];
            const double jx_c = S.jx[off2+ij], jy_c = S.jy[off2+ij], jz_c = S.jz[off2+ij];

            auto ew = [&](int k, int s) {
                if (k == P.kmin0) return Vt[ij + P.kmin0*ijdim + s*s3];
                if (k >  P.kmax0) return Vt[ij + P.kmax0*ijdim + s*s3];
                size_t idx = (size_t)ij + (size_t)k*ijdim;
                return S.afact[k]*Vt[idx+s*s3] + S.bfact[k]*Vt[idx-ijdim+s*s3];
            };
            auto eif = [&](int ki, double fq[4], double& fe, double& fp, double& fk, double& fu, double& fv) {
                fq[0]=fq[1]=fq[2]=fq[3]=0; fe=fp=fk=fu=fv=0;
                if (ki > P.kmax0) return;
                size_t kdx = (size_t)ij + (size_t)ki*ijdim;
                double zki = zh[kdx], dzk = zh[kdx+ijdim] - zki, dzm1 = (ki>P.kmin0)?(zki-zh[kdx-ijdim]):dzk;
                for (int s = 0; s < 4; ++s) {
                    double wk = ew(ki, s); if (wk == 0.0) continue;
                    double zdis = vadv1d_taylor_zdis(wk, ew(ki+1,s), (ki>P.kmin0)?ew(ki-1,s):0.0, dzk, dzm1, dzk+dzm1, dt, ki==P.kmin0);
                    if (zdis <= 0.0) continue;
                    double ztmp = zki + zdis; int c = ki; double zhc = zki;
                    for (int j = ki; j <= P.kmax0; ++j) {
                        size_t jdx = (size_t)ij + (size_t)j*ijdim;
                        if (ztmp > zh[jdx] && ztmp <= zh[jdx+ijdim]) { c = j; zhc = zh[jdx]; break; }
                    }
                    double rem = zdis - (zhc - zki);
                    for (int j = ki; j <= c; ++j) {
                        size_t jdx = (size_t)ij + (size_t)j*ijdim;
                        double L = (j < c) ? (zh[jdx+ijdim] - zh[jdx]) : rem;
                        double w_sj = L * rhogq[jdx + sp_nq[s]*s3] * (g2[jdx]/gs2[jdx]);
                        fq[s] -= w_sj;
                        double rgj = rhog[jdx], rvx = rhogvx[jdx], rvy = rhogvy[jdx], rvz = rhogvz[jdx];
                        fe -= w_sj * P.CVdry * tem[jdx];
                        fp -= w_sj * P.GRAV * z[jdx];
                        fk -= w_sj * (0.5*(rvx*rvx+rvy*rvy+rvz*rvz)/(rgj*rgj));
                        fu -= w_sj * (vx[jdx]*ix_c + vy[jdx]*iy_c + vz[jdx]*iz_c);
                        fv -= w_sj * (vx[jdx]*jx_c + vy[jdx]*jy_c + vz[jdx]*jz_c);
                    }
                }
            };

            double fq_lo[4], fe_lo, fp_lo, fk_lo, fu_lo, fv_lo;
            eif(P.kmin0, fq_lo, fe_lo, fp_lo, fk_lo, fu_lo, fv_lo);
            S.precip[offp+ij] = -fq_lo[0]*inv_dt;
            S.precip[offp+ij+ijdim] = -(fq_lo[1]+fq_lo[2]+fq_lo[3])*inv_dt;
            S.precip_lh[off2+ij] = (fq_lo[1]+fq_lo[2]+fq_lo[3])*P.LHF*inv_dt;
            S.precip_rhoein[off2+ij] = -fe_lo*inv_dt;
            S.precip_rhophi[off2+ij] = -fp_lo*inv_dt;
            S.precip_rhokin[off2+ij] = -fk_lo*inv_dt;
            S.precip_energy[off2+ij] = S.precip_lh[off2+ij] + S.precip_rhoein[off2+ij] + S.precip_rhophi[off2+ij] + S.precip_rhokin[off2+ij];

            double d_rg=0, d_rge=0, d_rkht=0, d_rkhn=0, d_kvp=0;
            for (int k = P.kmin0; k <= P.kmax0; ++k) {
                double fq_hi[4], fe_hi, fp_hi, fk_hi, fu_hi, fv_hi;
                eif(k+1, fq_hi, fe_hi, fp_hi, fk_hi, fu_hi, fv_hi);
                size_t idx = (size_t)ij + (size_t)k*ijdim; double r = S.rdgz[k];
                double dq0=-(fq_hi[0]-fq_lo[0])*r, dq1=-(fq_hi[1]-fq_lo[1])*r, dq2=-(fq_hi[2]-fq_lo[2])*r, dq3=-(fq_hi[3]-fq_lo[3])*r;
                double sdq = dq0+dq1+dq2+dq3;
                double rg_o = rhog[idx], rvx_o = rhogvx[idx], rvy_o = rhogvy[idx], rvz_o = rhogvz[idx];
                double rkh_o = 0.5*(rvx_o*rvx_o+rvy_o*rvy_o+rvz_o*rvz_o)/rg_o, rg_n = rg_o + sdq;
                double du = -(fu_hi-fu_lo)*r, dv = -(fv_hi-fv_lo)*r;
                double rvx_n = rvx_o + du*ix_c + dv*jx_c, rvy_n = rvy_o + du*iy_c + dv*jy_c, rvz_n = rvz_o + du*iz_c + dv*jz_c;
                double rge_n = rhoge[idx] - (fe_hi-fe_lo)*r - (fp_hi-fp_lo)*r - sdq*P.GRAV*z[idx];
                double rkh_t = rkh_o - (fk_hi-fk_lo)*r, rkh_n = 0.5*(rvx_n*rvx_n+rvy_n*rvy_n+rvz_n*rvz_n)/rg_n;
                rhog[idx]=rg_n; rhogvx[idx]=rvx_n; rhogvy[idx]=rvy_n; rhogvz[idx]=rvz_n;
                rhogq[idx+sp_nq[0]*s3]+=dq0; rhogq[idx+sp_nq[1]*s3]+=dq1; rhogq[idx+sp_nq[2]*s3]+=dq2; rhogq[idx+sp_nq[3]*s3]+=dq3;
                double kvr = 0.0;
                if (k > P.kmin0) {
                    size_t im1 = idx - ijdim;
                    double rhogh_n = (S.afact[k]/gs2[idx]*gs2h[idx])*rg_n + (S.bfact[k]/gs2[im1]*gs2h[idx])*d_rg;
                    kvr = 0.5*rhogw[idx]*rhogw[idx]/rhogh_n;
                    double t2 = d_rkht + (S.cfact[k-1]*kv[idx] + S.dfact[k-1]*kv[im1]);
                    double t1 = d_rkhn + ((S.cfact[k-1]*gs2[im1]/gs2h[idx])*kvr + (S.dfact[k-1]*gs2[im1]/gs2h[im1])*d_kvp);
                    double rge_fin = d_rge + (t2 - t1);
                    rhoge[im1] = rge_fin; S.rho[off3+im1] = d_rg / gs2[im1]; S.tem[off3+im1] = rge_fin / (d_rg * P.CVdry);
                    for (int n = 0; n < 6; ++n) S.q[off4+im1+n*s3] = rhogq[im1+n*s3] / d_rg;
                }
                d_rg=rg_n; d_rge=rge_n; d_rkht=rkh_t; d_rkhn=rkh_n; d_kvp=kvr;
                for (int s=0;s<4;++s) fq_lo[s]=fq_hi[s]; fe_lo=fe_hi; fp_lo=fp_hi; fk_lo=fk_hi; fu_lo=fu_hi; fv_lo=fv_hi;
            }
            size_t ixk = (size_t)ij + (size_t)P.kmax0*ijdim;
            double rge_top = d_rge + (d_rkht + S.dfact[P.kmax0]*kv[ixk]) - (d_rkhn + (S.dfact[P.kmax0]*gs2[ixk]/gs2h[ixk])*d_kvp);
            rhoge[ixk] = rge_top; S.rho[off3+ixk] = d_rg / gs2[ixk]; S.tem[off3+ixk] = rge_top / (d_rg * P.CVdry);
            for (int n = 0; n < 6; ++n) S.q[off4+ixk+n*s3] = rhogq[ixk+n*s3] / d_rg;
        });

        Kokkos::parallel_for("m0_satadj", Kokkos::RangePolicy<ExecSpace>(exec, 0, s3), KOKKOS_LAMBDA(const size_t idx) {
            const NSW6Params& P = *dP; const StatePtrs& S = *dS;
            int k = (int)(idx / (size_t)ijdim), ij = (int)(idx - (size_t)k * ijdim);
            int d0=0, d1=0, d2=0, d3=0, d4=0;
            satadj_and_post_negfilter_cell(P, ij, k, S.rhog+off3, S.rhoge+off3, S.rhogq+off4,
                S.rho+off3, S.tem+off3, S.pre+off3, S.q+off4, S.gsgam2+off3, &d0,&d1,&d2,&d3,&d4);
        });

        Kokkos::parallel_for("m0_cnvvar", Kokkos::RangePolicy<ExecSpace>(exec, 0, ijdim), KOKKOS_LAMBDA(const int ij) {
            const NSW6Params& P = *dP; const StatePtrs& S = *dS;
            cnvvar_bndcnd_cldfrac_column(P, ij, S.rhog+off3, S.rhogvx+off3, S.rhogvy+off3, S.rhogvz+off3,
                S.rhogw+off3, S.rhoge+off3, S.rhogq+off4, S.rho+off3, S.pre+off3, S.tem+off3,
                S.vx+off3, S.vy+off3, S.vz+off3, S.w+off3, S.q+off4, S.gdcfrc+off3, S.gdclw+off3,
                S.phi+off3, S.gsgam2+off3, S.gsgam2h+off3, S.afact, S.bfact);
        });
    }
}

// =============================================================================
// MODE 1 (Kokkos-ENG): Idiomatic whole-domain Kokkos, ORIGINAL NSW6 math
// =============================================================================
static void run_kokkos_mode1(const ExecSpace& exec, const NSW6Params* dP, int ijdim, int kdim, int lall, int do_rc,
                             const StatePtrs* dS, const WorkPtrs* dW) {
    const size_t s3 = (size_t)ijdim * kdim, n3 = (size_t)lall * s3, n2 = (size_t)lall * ijdim;

    Kokkos::parallel_for("m1_nsw6", Kokkos::RangePolicy<ExecSpace>(exec, 0, n3), KOKKOS_LAMBDA(const size_t gid) {
        const NSW6Params& P = *dP; const StatePtrs& S = *dS; const WorkPtrs& W = *dW;
        const int l = (int)(gid / s3);
        const size_t idx = gid - (size_t)l * s3;
        const int k = (int)(idx / (size_t)ijdim), ij = (int)(idx - (size_t)k * ijdim);
        const size_t off3 = (size_t)l * s3, off4 = (size_t)l * 6ULL * s3, off2 = (size_t)l * ijdim, offv = (size_t)l * 4ULL * s3;
        int d0=0, d1=0, d2=0, d3=0, d4=0;
        nsw6_pre_and_core_cell(P, ij, k, S.rhog+off3, S.rhoge+off3, S.rhogq+off4,
            S.rho+off3, S.tem+off3, S.pre+off3, S.q+off4, S.w+off3,
            S.unccn+off3, S.factor_n_aie+off2, S.gsgam2+off3, W.Vt+offv,
            &d0,&d1,&d2,&d3,&d4);
    });

    if (do_rc) kokkos_rceff_whole(exec, dP, ijdim, kdim, lall, dS, dW);
    kokkos_init_kv(exec, dP, ijdim, kdim, lall, dS, dW);

    Kokkos::parallel_for("m1_precip_h", Kokkos::RangePolicy<ExecSpace>(exec, 0, n2), KOKKOS_LAMBDA(const size_t col) {
        const NSW6Params& P = *dP; const StatePtrs& S = *dS; const WorkPtrs& W = *dW;
        const int l = (int)(col / (size_t)ijdim), ij = (int)(col - (size_t)l * ijdim);
        const int sp_nq[4] = { P.i_qr, P.i_qi, P.i_qs, P.i_qg };
        const double dt = P.dt;
        const size_t off3 = (size_t)l * s3, off4 = (size_t)l * 6ULL * s3, offv = (size_t)l * 4ULL * s3;
        double *rhogw = S.rhogw + off3, *kv = W.kv + off3;
        const double *q = S.q + off4, *z = S.z + off3, *zh = S.zh + off3;
        const double *gs2h = S.gsgam2h + off3, *g2h = S.gam2h + off3, *Vt = W.Vt + offv;
        const double dzh_kmin = z[ij + P.kmin0 * ijdim] - zh[ij + P.kmin0 * ijdim];
        double pfk = 0.0, pfw = 0.0;
        for (int k = P.kmin0 + 1; k <= P.kmax0; ++k) {
            size_t kdx = (size_t)ij + (size_t)k * ijdim;
            double zk = z[kdx], zkm1 = z[kdx - ijdim], dzh = zk - zkm1;
            double dzhm1 = (k - 1 == P.kmin0) ? dzh_kmin : (zkm1 - z[kdx - 2ULL * ijdim]);
            double fkv_k = 0.0, fw_k = 0.0;
            for (int s = 0; s < 4; ++s) {
                double wk = Vt[kdx + s * s3]; if (wk == 0.0) continue;
                double zdis = vadv1d_taylor_zdis(wk, Vt[kdx + ijdim + s * s3], Vt[kdx - ijdim + s * s3],
                                                 dzh, dzhm1, dzh + dzhm1, dt, k == P.kmin0 + 1);
                if (zdis <= 0.0) continue;
                double ztmp = zkm1 + zdis; int c = k; double zlo = zkm1;
                for (int j = k; j <= P.kmax0; ++j) {
                    double zj = z[(size_t)ij + (size_t)j * ijdim], zjm1 = z[(size_t)ij + (size_t)(j - 1) * ijdim];
                    if (ztmp > zjm1 && ztmp <= zj) { c = j; zlo = zjm1; break; }
                }
                double rem = zdis - (zlo - zkm1);
                for (int j = k; j <= c; ++j) {
                    size_t jdx = (size_t)ij + (size_t)j * ijdim;
                    double L = (j < c) ? (z[jdx] - z[jdx - ijdim]) : rem;
                    double qhj = (S.afact[j]*q[jdx + sp_nq[s]*s3] + S.bfact[j]*q[jdx - ijdim + sp_nq[s]*s3]) * (g2h[jdx]/gs2h[jdx]);
                    fkv_k -= L * qhj * kv[jdx];
                    fw_k  -= L * qhj * rhogw[jdx];
                }
            }
            if (fkv_k != 0.0 || pfk != 0.0 || fw_k != 0.0 || pfw != 0.0) {
                kv[kdx]    -= (fkv_k - pfk) * S.rdgzh[k];
                rhogw[kdx] -= (fw_k  - pfw) * S.rdgzh[k];
            }
            pfk = fkv_k; pfw = fw_k;
        }
    });

    Kokkos::parallel_for("m1_precip_f", Kokkos::RangePolicy<ExecSpace>(exec, 0, n2), KOKKOS_LAMBDA(const size_t col) {
        const NSW6Params& P = *dP; const StatePtrs& S = *dS; const WorkPtrs& W = *dW;
        const int l = (int)(col / (size_t)ijdim), ij = (int)(col - (size_t)l * ijdim);
        const int sp_nq[4] = { P.i_qr, P.i_qi, P.i_qs, P.i_qg };
        const double dt = P.dt, inv_dt = 1.0 / dt;
        const size_t off3 = (size_t)l * s3, off4 = (size_t)l * 6ULL * s3;
        const size_t offv = (size_t)l * 4ULL * s3, offp = (size_t)l * 2ULL * (size_t)ijdim;
        double *rhog = S.rhog+off3, *rhogvx = S.rhogvx+off3, *rhogvy = S.rhogvy+off3, *rhogvz = S.rhogvz+off3;
        double *rhoge = S.rhoge+off3, *rhogq = S.rhogq+off4;
        const double *rhogw = S.rhogw+off3, *tem = S.tem+off3, *vx = S.vx+off3, *vy = S.vy+off3, *vz = S.vz+off3;
        const double *z = S.z+off3, *zh = S.zh+off3, *gs2 = S.gsgam2+off3, *gs2h = S.gsgam2h+off3, *g2 = S.gam2+off3;
        const double *Vt = W.Vt+offv, *kv = W.kv+off3;
        const double ix_c = S.ix[col], iy_c = S.iy[col], iz_c = S.iz[col];
        const double jx_c = S.jx[col], jy_c = S.jy[col], jz_c = S.jz[col];

        auto ew = [&](int k, int s) {
            if (k == P.kmin0) return Vt[ij + P.kmin0*ijdim + s*s3];
            if (k >  P.kmax0) return Vt[ij + P.kmax0*ijdim + s*s3];
            size_t idx = (size_t)ij + (size_t)k*ijdim;
            return S.afact[k]*Vt[idx+s*s3] + S.bfact[k]*Vt[idx-ijdim+s*s3];
        };
        auto eif = [&](int ki, double fq[4], double& fe, double& fp, double& fk, double& fu, double& fv) {
            fq[0]=fq[1]=fq[2]=fq[3]=0; fe=fp=fk=fu=fv=0;
            if (ki > P.kmax0) return;
            size_t kdx = (size_t)ij + (size_t)ki*ijdim;
            double zki = zh[kdx], dzk = zh[kdx+ijdim] - zki, dzm1 = (ki>P.kmin0)?(zki-zh[kdx-ijdim]):dzk;
            for (int s = 0; s < 4; ++s) {
                double wk = ew(ki, s); if (wk == 0.0) continue;
                double zdis = vadv1d_taylor_zdis(wk, ew(ki+1,s), (ki>P.kmin0)?ew(ki-1,s):0.0, dzk, dzm1, dzk+dzm1, dt, ki==P.kmin0);
                if (zdis <= 0.0) continue;
                double ztmp = zki + zdis; int c = ki; double zhc = zki;
                for (int j = ki; j <= P.kmax0; ++j) {
                    size_t jdx = (size_t)ij + (size_t)j*ijdim;
                    if (ztmp > zh[jdx] && ztmp <= zh[jdx+ijdim]) { c = j; zhc = zh[jdx]; break; }
                }
                double rem = zdis - (zhc - zki);
                for (int j = ki; j <= c; ++j) {
                    size_t jdx = (size_t)ij + (size_t)j*ijdim;
                    double L = (j < c) ? (zh[jdx+ijdim] - zh[jdx]) : rem;
                    double w_sj = L * rhogq[jdx + sp_nq[s]*s3] * (g2[jdx]/gs2[jdx]);
                    fq[s] -= w_sj;
                    double rgj = rhog[jdx], rvx = rhogvx[jdx], rvy = rhogvy[jdx], rvz = rhogvz[jdx];
                    fe -= w_sj * P.CVdry * tem[jdx];
                    fp -= w_sj * P.GRAV * z[jdx];
                    fk -= w_sj * (0.5*(rvx*rvx+rvy*rvy+rvz*rvz)/(rgj*rgj));
                    fu -= w_sj * (vx[jdx]*ix_c + vy[jdx]*iy_c + vz[jdx]*iz_c);
                    fv -= w_sj * (vx[jdx]*jx_c + vy[jdx]*jy_c + vz[jdx]*jz_c);
                }
            }
        };

        double fq_lo[4], fe_lo, fp_lo, fk_lo, fu_lo, fv_lo;
        eif(P.kmin0, fq_lo, fe_lo, fp_lo, fk_lo, fu_lo, fv_lo);
        S.precip[offp+ij] = -fq_lo[0]*inv_dt;
        S.precip[offp+ij+ijdim] = -(fq_lo[1]+fq_lo[2]+fq_lo[3])*inv_dt;
        S.precip_lh[col] = (fq_lo[1]+fq_lo[2]+fq_lo[3])*P.LHF*inv_dt;
        S.precip_rhoein[col] = -fe_lo*inv_dt;
        S.precip_rhophi[col] = -fp_lo*inv_dt;
        S.precip_rhokin[col] = -fk_lo*inv_dt;
        S.precip_energy[col] = S.precip_lh[col] + S.precip_rhoein[col] + S.precip_rhophi[col] + S.precip_rhokin[col];

        double d_rg=0, d_rge=0, d_rkht=0, d_rkhn=0, d_kvp=0;
        for (int k = P.kmin0; k <= P.kmax0; ++k) {
            double fq_hi[4], fe_hi, fp_hi, fk_hi, fu_hi, fv_hi;
            eif(k+1, fq_hi, fe_hi, fp_hi, fk_hi, fu_hi, fv_hi);
            size_t idx = (size_t)ij + (size_t)k*ijdim; double r = S.rdgz[k];
            double dq0=-(fq_hi[0]-fq_lo[0])*r, dq1=-(fq_hi[1]-fq_lo[1])*r, dq2=-(fq_hi[2]-fq_lo[2])*r, dq3=-(fq_hi[3]-fq_lo[3])*r;
            double sdq = dq0+dq1+dq2+dq3;
            double rg_o = rhog[idx], rvx_o = rhogvx[idx], rvy_o = rhogvy[idx], rvz_o = rhogvz[idx];
            double rkh_o = 0.5*(rvx_o*rvx_o+rvy_o*rvy_o+rvz_o*rvz_o)/rg_o, rg_n = rg_o + sdq;
            double du = -(fu_hi-fu_lo)*r, dv = -(fv_hi-fv_lo)*r;
            double rvx_n = rvx_o + du*ix_c + dv*jx_c, rvy_n = rvy_o + du*iy_c + dv*jy_c, rvz_n = rvz_o + du*iz_c + dv*jz_c;
            double rge_n = rhoge[idx] - (fe_hi-fe_lo)*r - (fp_hi-fp_lo)*r - sdq*P.GRAV*z[idx];
            double rkh_t = rkh_o - (fk_hi-fk_lo)*r, rkh_n = 0.5*(rvx_n*rvx_n+rvy_n*rvy_n+rvz_n*rvz_n)/rg_n;
            if (sdq != 0.0) rhog[idx] = rg_n;
            if (du != 0.0 || dv != 0.0) { rhogvx[idx]=rvx_n; rhogvy[idx]=rvy_n; rhogvz[idx]=rvz_n; }
            if (dq0 != 0.0) rhogq[idx+sp_nq[0]*s3]+=dq0;
            if (dq1 != 0.0) rhogq[idx+sp_nq[1]*s3]+=dq1;
            if (dq2 != 0.0) rhogq[idx+sp_nq[2]*s3]+=dq2;
            if (dq3 != 0.0) rhogq[idx+sp_nq[3]*s3]+=dq3;
            double kvr = 0.0;
            if (k > P.kmin0) {
                size_t im1 = idx - ijdim;
                double rhogh_n = (S.afact[k]/gs2[idx]*gs2h[idx])*rg_n + (S.bfact[k]/gs2[im1]*gs2h[idx])*d_rg;
                kvr = 0.5*rhogw[idx]*rhogw[idx]/rhogh_n;
                double t2 = d_rkht + (S.cfact[k-1]*kv[idx] + S.dfact[k-1]*kv[im1]);
                double t1 = d_rkhn + ((S.cfact[k-1]*gs2[im1]/gs2h[idx])*kvr + (S.dfact[k-1]*gs2[im1]/gs2h[im1])*d_kvp);
                rhoge[im1] = d_rge + (t2 - t1);
            }
            d_rg=rg_n; d_rge=rge_n; d_rkht=rkh_t; d_rkhn=rkh_n; d_kvp=kvr;
            for (int s=0;s<4;++s) fq_lo[s]=fq_hi[s]; fe_lo=fe_hi; fp_lo=fp_hi; fk_lo=fk_hi; fu_lo=fu_hi; fv_lo=fv_hi;
        }
        size_t ixk = (size_t)ij + (size_t)P.kmax0*ijdim;
        rhoge[ixk] = d_rge + (d_rkht + S.dfact[P.kmax0]*kv[ixk]) - (d_rkhn + (S.dfact[P.kmax0]*gs2[ixk]/gs2h[ixk])*d_kvp);
    });

    // Fused SatAdj (original Newton) + Post-NF + Diag (reads rhogq/rg directly!)
    Kokkos::parallel_for("m1_satadj_diag", Kokkos::RangePolicy<ExecSpace>(exec, 0, n3), KOKKOS_LAMBDA(const size_t gid) {
        const NSW6Params& P = *dP; const StatePtrs& S = *dS;
        const int l = (int)(gid / s3);
        const size_t idx = gid - (size_t)l * s3;
        const int k = (int)(idx / (size_t)ijdim);
        const size_t off3 = (size_t)l * s3, off4 = (size_t)l * 6ULL * s3;
        double *rhog = S.rhog+off3, *rhoge = S.rhoge+off3, *rhogq = S.rhogq+off4;

        double rg = rhog[idx], rge = rhoge[idx], gs2 = S.gsgam2[off3 + idx];
        double rgqv = rhogq[idx+P.i_qv*s3], rgqc = rhogq[idx+P.i_qc*s3], rgqr = rhogq[idx+P.i_qr*s3];
        double rgqi = rhogq[idx+P.i_qi*s3], rgqs = rhogq[idx+P.i_qs*s3], rgqg = rhogq[idx+P.i_qg*s3];
        double qv = rgqv/rg, qc = rgqc/rg, qr = rgqr/rg, qi = rgqi/rg, qs = rgqs/rg, qg = rgqg/rg;
        double ein_m = rge / rg + qv * P.LHV, qsum = qv + qc, qv0 = qsum, qc0 = 0.0;
        double rho_v = rg / gs2, tem0 = (ein_m - qv0 * P.LHV) / P.CVdry;

        if (k >= P.kmin0 && k <= P.kmax0) {
            double rtem00 = 1.0 / P.TEM00, inv_rho_Rv = 1.0 / (rho_v * P.Rvap);
            double tc0 = (tem0 > 10.0) ? tem0 : 10.0;
            double qsat = P.PSAT0 * exp(P.LovR_liq * (rtem00 - 1.0 / tc0)) * inv_rho_Rv / tem0;
            if (qsum - qsat > P.EPS) {
                for (int ite = 1; ite <= 100; ++ite) {
                    double tc = (tem0 > 10.0) ? tem0 : 10.0;
                    qsat = P.PSAT0 * exp(P.LovR_liq * (rtem00 - 1.0 / tc)) * inv_rho_Rv / tem0;
                    qv0 = qsat; qc0 = qsum - qsat;
                    double Em_new = tem0 * P.CVdry + qv0 * P.LHV;
                    double dq_dT  = (P.LovR_liq / (tem0 * tem0)) * qsat;
                    double dtemp  = (Em_new - ein_m) / (P.CVdry + dq_dT * P.LHV);
                    tem0 -= dtemp;
                    if (fabs(dtemp) < P.dtemp_criteria) break;
                }
            }
        }
        if (k < P.kmin0 || k > P.kmax0) {
            rhogq[idx+P.i_qv*s3]=rg*qv0; rhogq[idx+P.i_qc*s3]=rg*qc0;
            rhoge[idx]=rg*tem0*P.CVdry; S.rho[off3+idx]=rho_v;
            for (int n=0;n<6;++n) S.q[off4+idx+n*s3]=0.0;
            S.gdcfrc[off3+idx]=0.0; S.gdclw[off3+idx]=P.UNDEF;
            return;
        }
        rgqv = rg*qv0; rgqc = rg*qc0;
        double diffq = 0.0;
        auto cq = [&](double& v){ if (v < 0.0) { diffq += v; v = 0.0; } };
        cq(rgqc); cq(rgqr); cq(rgqi); cq(rgqs); cq(rgqg);
        rgqv += diffq; double diff0 = 0.0;
        if (rgqv < 0.0) { diff0 = rgqv; rgqv = 0.0; }
        rg *= (1.0 - diff0); rho_v = rg / gs2;
        rhog[idx]=rg; rhoge[idx]=tem0*rg*P.CVdry;
        rhogq[idx+P.i_qv*s3]=rgqv; rhogq[idx+P.i_qc*s3]=rgqc; rhogq[idx+P.i_qr*s3]=rgqr;
        rhogq[idx+P.i_qi*s3]=rgqi; rhogq[idx+P.i_qs*s3]=rgqs; rhogq[idx+P.i_qg*s3]=rgqg;
        qv=rgqv/rg; qc=rgqc/rg; qr=rgqr/rg; qi=rgqi/rg; qs=rgqs/rg; qg=rgqg/rg;
        double qd = 1.0 - (qv+qc+qr+qi+qs+qg);
        S.rho[off3+idx]=rho_v; S.tem[off3+idx]=tem0; S.pre[off3+idx]=rho_v*(qd*P.Rdry+qv*P.Rvap)*tem0;
        S.vx[off3+idx]=S.rhogvx[off3+idx]/rg; S.vy[off3+idx]=S.rhogvy[off3+idx]/rg; S.vz[off3+idx]=S.rhogvz[off3+idx]/rg;
        S.q[off4+idx+P.i_qv*s3]=qv; S.q[off4+idx+P.i_qc*s3]=qc; S.q[off4+idx+P.i_qr*s3]=qr;
        S.q[off4+idx+P.i_qi*s3]=qi; S.q[off4+idx+P.i_qs*s3]=qs; S.q[off4+idx+P.i_qg*s3]=qg;
        S.gdcfrc[off3+idx] = ((qc+qr+qi+qs+qg) >= 1.0e-10) ? 1.0 : 0.0;
        S.gdclw [off3+idx] = P.UNDEF;
    });

    kokkos_diag_w_and_bndcnd(exec, dP, ijdim, kdim, lall, dS);
}

// =============================================================================
// MODE 2 (Kokkos-MATH) & MODE 3 (Kokkos-MONO):
//   100% Pure Portable Kokkos + TeamPolicy Scratch Compaction + 3 Math Theorems
// =============================================================================
template <int OPT_MODE>
static void run_kokkos_math(const ExecSpace& exec, const NSW6Params* dP, int ijdim, int kdim, int lall, int do_rc,
                            const StatePtrs* dS, KokkosWorkspace& WK) {
    const WorkPtrs* dW = WK.dp;
    const int kmin0 = 1, kmax0 = kdim - 2, kint = kmax0 - kmin0 + 1;
    const size_t s3 = (size_t)ijdim * kdim, n3 = (size_t)lall * s3;
    const size_t n2 = (size_t)lall * ijdim, n_int = (size_t)lall * (size_t)kint * (size_t)ijdim;

    // 0. Async zero via Kokkos::deep_copy(exec, ...) — ZERO host-device fences!
    Kokkos::deep_copy(exec, WK.Vt, 0.0);
    Kokkos::deep_copy(exec, WK.col_spmask, 0u);
    Kokkos::deep_copy(exec, WK.col_ktop, 0);
    Kokkos::deep_copy(exec, WK.wet_count, 0u);
    Kokkos::deep_copy(exec, WK.wet_col_count, 0u);

    // 1A. Portable Kokkos::TeamPolicy Dry-Cell Screen + Team-Aggregated Compaction
    constexpr int TSIZE = 256;
    const int n_teams_int = (int)((n_int + TSIZE - 1) / TSIZE);
    const size_t shmem_b  = ScratchU32::shmem_size(2);
    Kokkos::parallel_for("km_screen",
        TeamPolicy(exec, n_teams_int, TSIZE).set_scratch_size(0, Kokkos::PerTeam(shmem_b)),
        KOKKOS_LAMBDA(const TeamMember& team) {
            const NSW6Params& P = *dP; const StatePtrs& S = *dS; const WorkPtrs& W = *dW;
            ScratchU32 s_meta(team.team_scratch(0), 2);
            Kokkos::single(Kokkos::PerTeam(team), [&]() { s_meta(0) = 0u; });
            team.team_barrier();

            const size_t gid = (size_t)team.league_rank() * TSIZE + (size_t)team.team_rank();
            bool is_wet = false;
            uint32_t packed = 0u;
            if (gid < n_int) {
                const size_t s_int = (size_t)kint * ijdim;
                const int l  = (int)(gid / s_int);
                const size_t rem = gid - (size_t)l * s_int;
                const int k  = (int)(rem / (size_t)ijdim) + kmin0;
                const int ij = (int)(rem - (size_t)(k - kmin0) * ijdim);
                const size_t off3 = (size_t)l * s3, off4 = (size_t)l * 6ULL * s3;
                const size_t idx  = (size_t)ij + (size_t)k * ijdim;

                const double* rhogq = S.rhogq + off4;
                bool no_hyd = (rhogq[idx + P.i_qc*s3] == 0.0 && rhogq[idx + P.i_qr*s3] == 0.0 &&
                               rhogq[idx + P.i_qi*s3] == 0.0 && rhogq[idx + P.i_qs*s3] == 0.0 &&
                               rhogq[idx + P.i_qg*s3] == 0.0);
                double rg = S.rhog[off3 + idx], temp = S.tem[off3 + idx], rgqv = rhogq[idx + P.i_qv * s3];
                bool clean_dry = no_hyd && (rgqv >= 0.0);
                if (clean_dry && temp < P.TEM00) {
                    double dens = rg / S.gsgam2[off3 + idx], qv = rgqv / rg;
                    double tc = (temp > 10.0) ? temp : 10.0;
                    double psati = P.PSAT0 * exp(P.LovR_ice * (1.0 / P.TEM00 - 1.0 / tc));
                    if (qv > psati / (dens * P.Rvap * temp)) clean_dry = false;
                }
                if (clean_dry) {
                    S.rhoge[off3 + idx] = temp * rg * P.CVdry;
                } else {
                    is_wet = true;
                    packed = pack_cell(ij, k, l);
                }
            }
            uint32_t local_pos = is_wet ? Kokkos::atomic_fetch_add(&s_meta(0), 1u) : 0u;
            team.team_barrier();
            Kokkos::single(Kokkos::PerTeam(team), [&]() {
                if (s_meta(0) > 0u) s_meta(1) = Kokkos::atomic_fetch_add(&W.wet_count[0], s_meta(0));
            });
            team.team_barrier();
            if (is_wet) W.wet_cells[s_meta(1) + local_pos] = packed;
        });

    // 1B. Pure-Kokkos Wet-Cell Microphysics Core (Dyadic Roh-Satoh via nsw6_math_core.cuh!)
    Kokkos::parallel_for("km_wet_core", Kokkos::RangePolicy<ExecSpace>(exec, 0, n_int), KOKKOS_LAMBDA(const size_t p) {
        const WorkPtrs& W = *dW;
        if (p >= W.wet_count[0]) return;
        const NSW6Params& P = *dP; const StatePtrs& S = *dS;
        const uint32_t packed = W.wet_cells[p];
        const int ij = (int)(packed & 0x3FFFFu);
        const int k  = (int)((packed >> 18) & 0xFFu);
        const int l  = (int)(packed >> 26);

        const size_t off3 = (size_t)l * s3, off4 = (size_t)l * 6ULL * s3;
        const size_t off2 = (size_t)l * ijdim, offv = (size_t)l * 4ULL * s3;
        const size_t idx  = (size_t)ij + (size_t)k * ijdim;

        int d0=0, d1=0, d2=0, d3=0, d4=0;
        nsw6_pre_and_core_cell(P, ij, k, S.rhog+off3, S.rhoge+off3, S.rhogq+off4,
            S.rho+off3, S.tem+off3, S.pre+off3, S.q+off4, S.w+off3,
            S.unccn+off3, S.factor_n_aie+off2, S.gsgam2+off3, W.Vt+offv,
            &d0,&d1,&d2,&d3,&d4);

        const double* rhogq = S.rhogq + off4;
        uint32_t spm = (rhogq[idx + P.i_qr*s3] > 0.0 ? 1u : 0u)
                     | ((P.precip_4water && rhogq[idx + P.i_qi*s3] > 0.0) ? 2u : 0u)
                     | (rhogq[idx + P.i_qs*s3] > 0.0 ? 4u : 0u)
                     | (rhogq[idx + P.i_qg*s3] > 0.0 ? 8u : 0u);
        if (spm) {
            size_t col = off2 + (size_t)ij;
            if ((W.col_spmask[col] & spm) != spm) Kokkos::atomic_or(&W.col_spmask[col], spm);
            if (k > W.col_ktop[col])              Kokkos::atomic_max(&W.col_ktop[col], k);
        }
    });

    if (do_rc) kokkos_rceff_whole(exec, dP, ijdim, kdim, lall, dS, dW);
    kokkos_init_kv(exec, dP, ijdim, kdim, lall, dS, dW);

    // 1C. Portable TeamPolicy Precipitating-Column Compaction
    const int n_teams_col = (int)((n2 + TSIZE - 1) / TSIZE);
    Kokkos::parallel_for("km_compact_cols",
        TeamPolicy(exec, n_teams_col, TSIZE).set_scratch_size(0, Kokkos::PerTeam(shmem_b)),
        KOKKOS_LAMBDA(const TeamMember& team) {
            const StatePtrs& S = *dS; const WorkPtrs& W = *dW;
            ScratchU32 s_meta(team.team_scratch(0), 2);
            Kokkos::single(Kokkos::PerTeam(team), [&]() { s_meta(0) = 0u; });
            team.team_barrier();

            const size_t col = (size_t)team.league_rank() * TSIZE + (size_t)team.team_rank();
            bool is_wet = false;
            uint32_t packed = 0u;
            if (col < n2) {
                const int l  = (int)(col / (size_t)ijdim);
                const int ij = (int)(col - (size_t)l * ijdim);
                if (W.col_spmask[col] == 0u) {
                    const size_t offp = (size_t)l * 2ULL * (size_t)ijdim;
                    S.precip[offp + ij] = 0.0; S.precip[offp + ij + ijdim] = 0.0;
                    S.precip_lh[col] = 0.0; S.precip_rhoein[col] = 0.0;
                    S.precip_rhophi[col] = 0.0; S.precip_rhokin[col] = 0.0; S.precip_energy[col] = 0.0;
                } else {
                    is_wet = true;
                    packed = pack_col(ij, l);
                }
            }
            uint32_t local_pos = is_wet ? Kokkos::atomic_fetch_add(&s_meta(0), 1u) : 0u;
            team.team_barrier();
            Kokkos::single(Kokkos::PerTeam(team), [&]() {
                if (s_meta(0) > 0u) s_meta(1) = Kokkos::atomic_fetch_add(&W.wet_col_count[0], s_meta(0));
            });
            team.team_barrier();
            if (is_wet) W.wet_cols[s_meta(1) + local_pos] = packed;
        });

    // 2A. Pure-Kokkos Fubini Half-Level Sedimentation (+ optional O(K) Monotone)
    Kokkos::parallel_for("km_precip_h", Kokkos::RangePolicy<ExecSpace>(exec, 0, n2), KOKKOS_LAMBDA(const size_t p) {
        const WorkPtrs& W = *dW;
        if (p >= W.wet_col_count[0]) return;
        const NSW6Params& P = *dP; const StatePtrs& S = *dS;
        const uint32_t pc = W.wet_cols[p];
        const int ij = (int)(pc & 0x3FFFFu), l = (int)(pc >> 18);
        const size_t col = (size_t)l * ijdim + ij;
        const uint32_t sp_mask = W.col_spmask[col];
        const int k_top = W.col_ktop[col], k_hi = (k_top + 1 <= P.kmax0) ? (k_top + 1) : P.kmax0;
        const int sp_nq[4] = { P.i_qr, P.i_qi, P.i_qs, P.i_qg };
        const double dt = P.dt;
        const size_t off3 = (size_t)l * s3, off4 = (size_t)l * 6ULL * s3, offv = (size_t)l * 4ULL * s3;
        double *rhogw = S.rhogw + off3, *kv = W.kv + off3;
        const double *q = S.q + off4, *z = S.z + off3, *zh = S.zh + off3;
        const double *gs2h = S.gsgam2h + off3, *g2h = S.gam2h + off3, *Vt = W.Vt + offv;
        const double dzh_kmin = z[ij + P.kmin0 * ijdim] - zh[ij + P.kmin0 * ijdim];
        double* zd_m = (OPT_MODE == 3) ? (W.zdis_m + offv) : nullptr;

        if constexpr (OPT_MODE == 3) {
            for (int s = 0; s < 4; ++s) {
                if (!(sp_mask & (1u << s))) continue;
                double xi_next = z[ij + (size_t)k_hi * ijdim];
                for (int k = k_hi; k >= P.kmin0 + 1; --k) {
                    size_t idx = (size_t)ij + (size_t)k * ijdim;
                    double wk = Vt[idx + s * s3], d_raw = 0.0, z_km1 = z[idx - ijdim];
                    if (wk != 0.0) {
                        double dzh_k = z[idx] - z_km1;
                        double dzh_km1 = (k - 1 == P.kmin0) ? dzh_kmin : (z_km1 - z[idx - 2ULL * ijdim]);
                        d_raw = vadv1d_taylor_zdis(wk, Vt[idx+ijdim+s*s3], Vt[idx-ijdim+s*s3], dzh_k, dzh_km1, dzh_k+dzh_km1, dt, k == P.kmin0 + 1);
                    }
                    double xi_k = z_km1 + d_raw;
                    double xi_p = (xi_k < xi_next) ? xi_k : xi_next;
                    double d_p  = xi_p - z_km1;
                    zd_m[idx + s * s3] = (d_p > 0.0) ? d_p : 0.0;
                    xi_next = xi_p;
                }
            }
        }

        double prev_fkv = 0.0, prev_fw = 0.0;
        for (int k = P.kmin0 + 1; k <= k_hi; ++k) {
            size_t kdx = (size_t)ij + (size_t)k * ijdim;
            double z_k = z[kdx], z_km1 = z[kdx - ijdim], dzh_k = z_k - z_km1;
            double dzh_km1 = (k - 1 == P.kmin0) ? dzh_kmin : (z_km1 - z[kdx - 2ULL * ijdim]);
            int kc[4] = { k-1, k-1, k-1, k-1 }, c_max = k - 1;
            double rem[4] = { 0, 0, 0, 0 };
            for (int s = 0; s < 4; ++s) {
                if (!(sp_mask & (1u << s))) continue;
                double zdis;
                if constexpr (OPT_MODE == 3) zdis = zd_m[kdx + s * s3];
                else {
                    double wk = Vt[kdx + s * s3]; if (wk == 0.0) continue;
                    zdis = vadv1d_taylor_zdis(wk, Vt[kdx+ijdim+s*s3], Vt[kdx-ijdim+s*s3], dzh_k, dzh_km1, dzh_k+dzh_km1, dt, k == P.kmin0 + 1);
                }
                if (zdis <= 0.0) continue;
                double ztmp = z_km1 + zdis;
                if (ztmp <= z_k) { kc[s] = k; rem[s] = zdis; if (k > c_max) c_max = k; }
                else {
                    int c = k; double z_lo_c = z_km1;
                    for (int j = k + 1; j <= P.kmax0; ++j) {
                        double zj = z[(size_t)ij + (size_t)j * ijdim];
                        if (ztmp <= zj) {
                            double zjm1 = z[(size_t)ij + (size_t)(j - 1) * ijdim];
                            if (ztmp > zjm1) { c = j; z_lo_c = zjm1; }
                            break;
                        }
                    }
                    kc[s] = c; rem[s] = zdis - (z_lo_c - z_km1);
                    if (c > c_max) c_max = c;
                }
            }
            double fkv_k = 0.0, fw_k = 0.0;
            if (c_max == k) {
                size_t jm1 = kdx - ijdim;
                double af = S.afact[k], bf = S.bfact[k], Wh = 0.0;
                for (int s = 0; s < 4; ++s) if (kc[s] == k) Wh += rem[s] * (af * q[kdx + sp_nq[s]*s3] + bf * q[jm1 + sp_nq[s]*s3]);
                if (Wh != 0.0) { double wr = Wh * (g2h[kdx] / gs2h[kdx]); fkv_k = -wr * kv[kdx]; fw_k = -wr * rhogw[kdx]; }
            } else if (c_max > k) {
                int j_end = (c_max <= k_hi) ? c_max : k_hi;
                for (int j = k; j <= j_end; ++j) {
                    size_t jdx = (size_t)ij + (size_t)j * ijdim, jm1 = jdx - ijdim;
                    double dj = z[jdx] - z[jm1], af = S.afact[j], bf = S.bfact[j], Wh = 0.0;
                    for (int s = 0; s < 4; ++s) {
                        int c = kc[s]; if (c < j) continue;
                        double L = (j < c) ? dj : rem[s];
                        Wh += L * (af * q[jdx + sp_nq[s]*s3] + bf * q[jm1 + sp_nq[s]*s3]);
                    }
                    if (Wh != 0.0) { double wr = Wh * (g2h[jdx] / gs2h[jdx]); fkv_k -= wr * kv[jdx]; fw_k -= wr * rhogw[jdx]; }
                }
            }
            if (fkv_k != 0.0 || prev_fkv != 0.0 || fw_k != 0.0 || prev_fw != 0.0) {
                kv[kdx]    -= (fkv_k - prev_fkv) * S.rdgzh[k];
                rhogw[kdx] -= (fw_k  - prev_fw)  * S.rdgzh[k];
            }
            prev_fkv = fkv_k; prev_fw = fw_k;
        }
    });

    // 2B. Pure-Kokkos Fubini Full-Level Sedimentation (+ optional O(K) Monotone)
    Kokkos::parallel_for("km_precip_f", Kokkos::RangePolicy<ExecSpace>(exec, 0, n2), KOKKOS_LAMBDA(const size_t p) {
        const WorkPtrs& W = *dW;
        if (p >= W.wet_col_count[0]) return;
        const NSW6Params& P = *dP; const StatePtrs& S = *dS;
        const uint32_t pc = W.wet_cols[p];
        const int ij = (int)(pc & 0x3FFFFu), l = (int)(pc >> 18);
        const size_t col = (size_t)l * ijdim + ij;
        const uint32_t sp_mask = W.col_spmask[col];
        const int k_hi = W.col_ktop[col];
        const int sp_nq[4] = { P.i_qr, P.i_qi, P.i_qs, P.i_qg };
        const double dt = P.dt, inv_dt = 1.0 / dt;
        const size_t off3 = (size_t)l * s3, off4 = (size_t)l * 6ULL * s3;
        const size_t offp = (size_t)l * 2ULL * (size_t)ijdim, offv = (size_t)l * 4ULL * s3;

        double *rhog = S.rhog+off3, *rhogvx = S.rhogvx+off3, *rhogvy = S.rhogvy+off3, *rhogvz = S.rhogvz+off3;
        double *rhoge = S.rhoge+off3, *rhogq = S.rhogq+off4;
        const double *rhogw = S.rhogw+off3, *tem = S.tem+off3, *vx = S.vx+off3, *vy = S.vy+off3, *vz = S.vz+off3;
        const double *z = S.z+off3, *zh = S.zh+off3, *gs2 = S.gsgam2+off3, *gs2h = S.gsgam2h+off3, *g2 = S.gam2+off3;
        const double *Vt = W.Vt+offv, *kv = W.kv+off3;
        const double ix_c = S.ix[col], iy_c = S.iy[col], iz_c = S.iz[col];
        const double jx_c = S.jx[col], jy_c = S.jy[col], jz_c = S.jz[col];

        auto ew = [&](int k, int s) {
            if (k == P.kmin0) return Vt[ij + P.kmin0*ijdim + s*s3];
            if (k >  P.kmax0) return Vt[ij + P.kmax0*ijdim + s*s3];
            size_t idx = (size_t)ij + (size_t)k*ijdim;
            return S.afact[k]*Vt[idx+s*s3] + S.bfact[k]*Vt[idx-ijdim+s*s3];
        };
        auto ez = [&](int k, int s) {
            double wk = ew(k, s); if (wk == 0.0) return 0.0;
            size_t idx = (size_t)ij + (size_t)k*ijdim;
            double zhk = zh[idx], dzk = zh[idx+ijdim] - zhk, dzm1 = (k>P.kmin0)?(zhk-zh[idx-ijdim]):dzk;
            return vadv1d_taylor_zdis(wk, ew(k+1,s), (k>P.kmin0)?ew(k-1,s):0.0, dzk, dzm1, dzk+dzm1, dt, k==P.kmin0);
        };
        double* zd_m = (OPT_MODE == 3) ? (W.zdis_m + offv) : nullptr;
        if constexpr (OPT_MODE == 3) {
            for (int s = 0; s < 4; ++s) {
                if (!(sp_mask & (1u << s))) continue;
                double xi_next = zh[ij + (size_t)(k_hi + 1) * ijdim];
                for (int k = k_hi; k >= P.kmin0; --k) {
                    size_t idx = (size_t)ij + (size_t)k * ijdim;
                    double zhk = zh[idx], xi_k = zhk + ez(k, s);
                    double xi_p = (xi_k < xi_next) ? xi_k : xi_next;
                    double d_p = xi_p - zhk;
                    zd_m[idx + s * s3] = (d_p > 0.0) ? d_p : 0.0;
                    xi_next = xi_p;
                }
            }
        }

        auto eif = [&](int ki, double fq[4], double& fe, double& fp, double& fk, double& fu, double& fv) {
            fq[0]=fq[1]=fq[2]=fq[3]=0; fe=fp=fk=fu=fv=0;
            if (ki > k_hi) return;
            size_t kdx = (size_t)ij + (size_t)ki*ijdim;
            double zki = zh[kdx], zkp1 = zh[kdx + ijdim];
            int kc[4] = { ki-1, ki-1, ki-1, ki-1 }, c_max = ki - 1;
            double rem[4] = { 0, 0, 0, 0 };
            for (int s = 0; s < 4; ++s) {
                if (!(sp_mask & (1u << s))) continue;
                double zdis = (OPT_MODE == 3) ? zd_m[kdx + s * s3] : ez(ki, s);
                if (zdis <= 0.0) continue;
                double ztmp = zki + zdis;
                if (ztmp <= zkp1) { kc[s] = ki; rem[s] = zdis; if (ki > c_max) c_max = ki; }
                else {
                    int c = ki; double zhc = zki;
                    for (int j = ki + 1; j <= P.kmax0; ++j) {
                        size_t jdx = (size_t)ij + (size_t)j * ijdim;
                        if (ztmp <= zh[jdx + ijdim]) {
                            double zj = zh[jdx]; if (ztmp > zj) { c = j; zhc = zj; } break;
                        }
                    }
                    kc[s] = c; rem[s] = zdis - (zhc - zki);
                    if (c > c_max) c_max = c;
                }
            }
            if (c_max == ki) {
                double rgs_j = g2[kdx] / gs2[kdx], W_tot = 0.0;
                for (int s = 0; s < 4; ++s) if (kc[s] == ki) { double w_sj = rem[s] * rhogq[kdx + sp_nq[s]*s3] * rgs_j; fq[s] = -w_sj; W_tot += w_sj; }
                if (W_tot != 0.0) {
                    double rgj = rhog[kdx], rvx = rhogvx[kdx], rvy = rhogvy[kdx], rvz = rhogvz[kdx];
                    double kh_r = 0.5*(rvx*rvx+rvy*rvy+rvz*rvz)/(rgj*rgj);
                    fe = -P.CVdry*(W_tot*tem[kdx]); fp = -P.GRAV*(W_tot*z[kdx]); fk = -W_tot*kh_r;
                    fu = -W_tot*(vx[kdx]*ix_c + vy[kdx]*iy_c + vz[kdx]*iz_c);
                    fv = -W_tot*(vx[kdx]*jx_c + vy[kdx]*jy_c + vz[kdx]*jz_c);
                }
            } else if (c_max > ki) {
                int j_end = (c_max <= k_hi) ? c_max : k_hi;
                double feu = 0.0, fpu = 0.0;
                for (int j = ki; j <= j_end; ++j) {
                    size_t jdx = (size_t)ij + (size_t)j * ijdim;
                    double dj = zh[jdx + ijdim] - zh[jdx], rgs_j = g2[jdx] / gs2[jdx], W_tot = 0.0;
                    for (int s = 0; s < 4; ++s) {
                        int c = kc[s]; if (c < j) continue;
                        double L = (j < c) ? dj : rem[s];
                        double w_sj = L * rhogq[jdx + sp_nq[s]*s3] * rgs_j;
                        fq[s] -= w_sj; W_tot += w_sj;
                    }
                    if (W_tot != 0.0) {
                        double rgj = rhog[jdx], rvx = rhogvx[jdx], rvy = rhogvy[jdx], rvz = rhogvz[jdx];
                        double kh_r = 0.5*(rvx*rvx+rvy*rvy+rvz*rvz)/(rgj*rgj);
                        feu -= W_tot * tem[jdx]; fpu -= W_tot * z[jdx]; fk -= W_tot * kh_r;
                        fu -= W_tot * (vx[jdx]*ix_c + vy[jdx]*iy_c + vz[jdx]*iz_c);
                        fv -= W_tot * (vx[jdx]*jx_c + vy[jdx]*jy_c + vz[jdx]*jz_c);
                    }
                }
                fe = P.CVdry * feu; fp = P.GRAV * fpu;
            }
        };

        double fq_lo[4], fe_lo, fp_lo, fk_lo, fu_lo, fv_lo;
        eif(P.kmin0, fq_lo, fe_lo, fp_lo, fk_lo, fu_lo, fv_lo);
        S.precip[offp+ij] = -fq_lo[0]*inv_dt;
        S.precip[offp+ij+ijdim] = -(fq_lo[1]+fq_lo[2]+fq_lo[3])*inv_dt;
        S.precip_lh[col] = (fq_lo[1]+fq_lo[2]+fq_lo[3])*P.LHF*inv_dt;
        S.precip_rhoein[col] = -fe_lo*inv_dt;
        S.precip_rhophi[col] = -fp_lo*inv_dt;
        S.precip_rhokin[col] = -fk_lo*inv_dt;
        S.precip_energy[col] = S.precip_lh[col] + S.precip_rhoein[col] + S.precip_rhophi[col] + S.precip_rhokin[col];

        double d_rg=0, d_rge=0, d_rkht=0, d_rkhn=0, d_kvp=0;
        for (int k = P.kmin0; k <= k_hi; ++k) {
            double fq_hi[4], fe_hi, fp_hi, fk_hi, fu_hi, fv_hi;
            eif(k+1, fq_hi, fe_hi, fp_hi, fk_hi, fu_hi, fv_hi);
            size_t idx = (size_t)ij + (size_t)k*ijdim; double r = S.rdgz[k];
            double dq0=-(fq_hi[0]-fq_lo[0])*r, dq1=-(fq_hi[1]-fq_lo[1])*r, dq2=-(fq_hi[2]-fq_lo[2])*r, dq3=-(fq_hi[3]-fq_lo[3])*r;
            double sdq = dq0+dq1+dq2+dq3;
            double rg_o = rhog[idx], rvx_o = rhogvx[idx], rvy_o = rhogvy[idx], rvz_o = rhogvz[idx];
            double rkh_o = 0.5*(rvx_o*rvx_o+rvy_o*rvy_o+rvz_o*rvz_o)/rg_o, rg_n = rg_o + sdq;
            double du = -(fu_hi-fu_lo)*r, dv = -(fv_hi-fv_lo)*r;
            double rvx_n = rvx_o + du*ix_c + dv*jx_c, rvy_n = rvy_o + du*iy_c + dv*jy_c, rvz_n = rvz_o + du*iz_c + dv*jz_c;
            double rge_n = rhoge[idx] - (fe_hi-fe_lo)*r - (fp_hi-fp_lo)*r - sdq*P.GRAV*z[idx];
            double rkh_t = rkh_o - (fk_hi-fk_lo)*r, rkh_n = 0.5*(rvx_n*rvx_n+rvy_n*rvy_n+rvz_n*rvz_n)/rg_n;
            if (sdq != 0.0) rhog[idx] = rg_n;
            if (du != 0.0 || dv != 0.0) { rhogvx[idx]=rvx_n; rhogvy[idx]=rvy_n; rhogvz[idx]=rvz_n; }
            if (dq0 != 0.0) rhogq[idx+sp_nq[0]*s3]+=dq0;
            if (dq1 != 0.0) rhogq[idx+sp_nq[1]*s3]+=dq1;
            if (dq2 != 0.0) rhogq[idx+sp_nq[2]*s3]+=dq2;
            if (dq3 != 0.0) rhogq[idx+sp_nq[3]*s3]+=dq3;
            double kvr = 0.0;
            if (k > P.kmin0) {
                size_t im1 = idx - ijdim;
                double rhogh_n = (S.afact[k]/gs2[idx]*gs2h[idx])*rg_n + (S.bfact[k]/gs2[im1]*gs2h[idx])*d_rg;
                kvr = 0.5*rhogw[idx]*rhogw[idx]/rhogh_n;
                double t2 = d_rkht + (S.cfact[k-1]*kv[idx] + S.dfact[k-1]*kv[im1]);
                double t1 = d_rkhn + ((S.cfact[k-1]*gs2[im1]/gs2h[idx])*kvr + (S.dfact[k-1]*gs2[im1]/gs2h[im1])*d_kvp);
                d_rge += (t2 - t1);
                rhoge[im1] = d_rge;
            }
            d_rg=rg_n; d_rge=rge_n; d_rkht=rkh_t; d_rkhn=rkh_n; d_kvp=kvr;
            for (int s=0;s<4;++s) fq_lo[s]=fq_hi[s]; fe_lo=fe_hi; fp_lo=fp_hi; fk_lo=fk_hi; fu_lo=fu_hi; fv_lo=fv_hi;
        }
        size_t idx = (size_t)ij + (size_t)k_hi * ijdim;
        double gs2_k = gs2[idx];
        if (k_hi < P.kmax0) {
            size_t ip1 = idx + ijdim;
            double gs2h_p1 = gs2h[ip1], kv_p1_t = kv[ip1];
            double rhogh_p1 = (S.afact[k_hi+1]/gs2[ip1]*gs2h_p1)*rhog[ip1] + (S.bfact[k_hi+1]/gs2_k*gs2h_p1)*d_rg;
            double rw_p1 = rhogw[ip1], kv_p1_n = 0.5*rw_p1*rw_p1/rhogh_p1;
            double t2 = d_rkht + (S.cfact[k_hi]*kv_p1_t + S.dfact[k_hi]*kv[idx]);
            double t1 = d_rkhn + ((S.cfact[k_hi]*gs2_k/gs2h_p1)*kv_p1_n + (S.dfact[k_hi]*gs2_k/gs2h[idx])*d_kvp);
            rhoge[idx] = d_rge + (t2 - t1);
            double gs2_p1 = gs2[ip1], kv_p2 = (k_hi + 1 < P.kmax0) ? kv[ip1 + ijdim] : 0.0;
            double d_p1 = (S.cfact[k_hi+1]*(1.0 - gs2_p1/gs2h[ip1+ijdim]))*kv_p2 + S.dfact[k_hi+1]*(kv_p1_t - (gs2_p1/gs2h_p1)*kv_p1_n);
            if (d_p1 != 0.0) rhoge[ip1] += d_p1;
        } else {
            rhoge[idx] = d_rge + (d_rkht + S.dfact[P.kmax0]*kv[idx]) - (d_rkhn + (S.dfact[P.kmax0]*gs2_k/gs2h[idx])*d_kvp);
        }
    });

    // 3. Pure-Kokkos Laguerre Householder-4 SatAdj + Post-NF + Conditional Diag
    Kokkos::parallel_for("km_satadj_diag", Kokkos::RangePolicy<ExecSpace>(exec, 0, n3), KOKKOS_LAMBDA(const size_t gid) {
        const NSW6Params& P = *dP; const StatePtrs& S = *dS; const WorkPtrs& W = *dW;
        const int l = (int)(gid / s3);
        const size_t idx = gid - (size_t)l * s3;
        const int k = (int)(idx / (size_t)ijdim), ij = (int)(idx - (size_t)k * ijdim);
        const size_t off3 = (size_t)l * s3, off4 = (size_t)l * 6ULL * s3;
        double *rhog = S.rhog+off3, *rhoge = S.rhoge+off3, *rhogq = S.rhogq+off4;

        double rg = rhog[idx], rge = rhoge[idx], gs2 = S.gsgam2[off3 + idx];
        bool rge_mod = false;
        int ktop = W.col_ktop[(size_t)l * ijdim + ij];
        if (k >= P.kmin0 && k <= P.kmax0 && (ktop == 0 || k >= ktop + 2)) {
            double gs2h_k = S.gsgam2h[off3+idx], gs2h_p1 = S.gsgam2h[off3+idx+ijdim];
            if (gs2 != gs2h_k || gs2 != gs2h_p1) {
                const double* kv = W.kv + off3;
                double de = (S.cfact[k]*(1.0 - gs2/gs2h_p1))*kv[idx+ijdim] + (S.dfact[k]*(1.0 - gs2/gs2h_k))*kv[idx];
                if (de != 0.0) { rge += de; rge_mod = true; }
            }
        }
        double rgqv = rhogq[idx+P.i_qv*s3], rgqc = rhogq[idx+P.i_qc*s3], rgqr = rhogq[idx+P.i_qr*s3];
        double rgqi = rhogq[idx+P.i_qi*s3], rgqs = rhogq[idx+P.i_qs*s3], rgqg = rhogq[idx+P.i_qg*s3];
        double qv = (rgqv!=0.0)?(rgqv/rg):0.0, qc = (rgqc!=0.0)?(rgqc/rg):0.0, qr = (rgqr!=0.0)?(rgqr/rg):0.0;
        double qi = (rgqi!=0.0)?(rgqi/rg):0.0, qs = (rgqs!=0.0)?(rgqs/rg):0.0, qg = (rgqg!=0.0)?(rgqg/rg):0.0;
        double ein_m = rge / rg + qv * P.LHV, qsum = qv + qc, qv0 = qsum, qc0 = 0.0;
        double rho_v = rg / gs2, tem0 = (ein_m - qv0 * P.LHV) / P.CVdry;
        bool q_mod = (qc != 0.0);

        if (k >= P.kmin0 && k <= P.kmax0 && qsum > 0.0) {
            double rtem00 = 1.0 / P.TEM00, inv_rho_Rv = 1.0 / (rho_v * P.Rvap);
            double tc0 = (tem0 > 10.0) ? tem0 : 10.0;
            double qsat = P.PSAT0 * exp(P.LovR_liq * (rtem00 - 1.0 / tc0)) * inv_rho_Rv / tem0;
            if (qsum - qsat > P.EPS) {
                q_mod = true;
                const double CV = P.CVdry, L = P.LHV, a = P.LovR_liq;
                for (int ite = 1; ite <= 20; ++ite) {
                    if (ite > 1) {
                        double tc = (tem0 > 10.0) ? tem0 : 10.0;
                        qsat = P.PSAT0 * exp(a * (rtem00 - 1.0 / tc)) * inv_rho_Rv / tem0;
                    }
                    double rT = 1.0 / tem0, b = a * rT, Lq = L * qsat;
                    double f  = CV * tem0 + Lq - ein_m;
                    double f1 = CV + Lq * rT * (b - 1.0);
                    double f2 = Lq * (rT * rT) * (b * b - 4.0 * b + 2.0);
                    double f3 = Lq * (rT * rT * rT) * (b * (b * (b - 9.0) + 18.0) - 6.0);
                    double dT_N = f / f1;
                    double num  = f * (f1 * f1 - 0.5 * f * f2);
                    double den  = f1 * (f1 * f1 - f * f2) + (1.0 / 6.0) * f * f * f3;
                    double dT   = (den != 0.0) ? (num / den) : dT_N;
                    if (!(dT * dT_N > 0.0 && fabs(dT) <= 2.0 * fabs(dT_N))) dT = dT_N;
                    tem0 -= dT; if (tem0 < 10.0) tem0 = 10.0;
                    if (fabs(dT) < 1.0e-12) {
                        double qn = qsat + qsat * rT * (b - 1.0) * (-dT);
                        qv0 = qn; qc0 = qsum - qn; break;
                    }
                }
            }
        }
        if (k < P.kmin0 || k > P.kmax0) {
            rhogq[idx+P.i_qv*s3] = rg*qv0; rhogq[idx+P.i_qc*s3] = rg*qc0;
            rhoge[idx] = rg * tem0 * P.CVdry; S.rho[off3+idx] = rho_v;
            for (int n=0;n<6;++n) S.q[off4+idx+n*s3] = 0.0;
            S.gdcfrc[off3+idx] = 0.0; S.gdclw[off3+idx] = P.UNDEF;
            return;
        }
        if (q_mod) { rgqv = rg*qv0; rgqc = rg*qc0; }
        double diffq = 0.0; bool neg = false;
        auto cq = [&](double& v){ if (v < 0.0) { diffq += v; v = 0.0; neg = true; } };
        cq(rgqc); cq(rgqr); cq(rgqi); cq(rgqs); cq(rgqg);
        rgqv += diffq; double diff0 = 0.0;
        if (rgqv < 0.0) { diff0 = rgqv; rgqv = 0.0; neg = true; }
        if (neg) {
            rg *= (1.0 - diff0); rho_v = rg / gs2; rhog[idx] = rg;
            rhogq[idx+P.i_qv*s3]=rgqv; rhogq[idx+P.i_qc*s3]=rgqc; rhogq[idx+P.i_qr*s3]=rgqr;
            rhogq[idx+P.i_qi*s3]=rgqi; rhogq[idx+P.i_qs*s3]=rgqs; rhogq[idx+P.i_qg*s3]=rgqg;
            qv=(rgqv!=0)?rgqv/rg:0; qc=(rgqc!=0)?rgqc/rg:0; qr=(rgqr!=0)?rgqr/rg:0;
            qi=(rgqi!=0)?rgqi/rg:0; qs=(rgqs!=0)?rgqs/rg:0; qg=(rgqg!=0)?rgqg/rg:0;
        } else if (q_mod) {
            rhogq[idx+P.i_qv*s3]=rgqv; rhogq[idx+P.i_qc*s3]=rgqc; qv=qv0; qc=qc0;
        }
        if (rge_mod || q_mod || neg) rhoge[idx] = tem0 * rg * P.CVdry;
        bool mass_dirty = (ktop > 0 && k <= ktop) || neg;
        bool rge_dirty  = (ktop > 0 && k <= ktop + 1) || rge_mod || q_mod || neg;
        if (mass_dirty) {
            S.rho[off3+idx] = rho_v;
            S.vx[off3+idx] = S.rhogvx[off3+idx]/rg; S.vy[off3+idx] = S.rhogvy[off3+idx]/rg; S.vz[off3+idx] = S.rhogvz[off3+idx]/rg;
            S.q[off4+idx+P.i_qv*s3]=qv; S.q[off4+idx+P.i_qc*s3]=qc; S.q[off4+idx+P.i_qr*s3]=qr;
            S.q[off4+idx+P.i_qi*s3]=qi; S.q[off4+idx+P.i_qs*s3]=qs; S.q[off4+idx+P.i_qg*s3]=qg;
        } else if (q_mod) { S.q[off4+idx+P.i_qv*s3]=qv; S.q[off4+idx+P.i_qc*s3]=qc; }
        if (rge_dirty) {
            double qd = 1.0 - (qv+qc+qr+qi+qs+qg);
            S.tem[off3+idx] = tem0; S.pre[off3+idx] = rho_v * (qd*P.Rdry + qv*P.Rvap) * tem0;
        }
        double qtot = qc+qr+qi+qs+qg;
        S.gdcfrc[off3+idx] = (qtot >= 1.0e-10) ? 1.0 : 0.0;
        S.gdclw[off3+idx]  = P.UNDEF;
    });

    kokkos_diag_w_and_bndcnd(exec, dP, ijdim, kdim, lall, dS);
}

int main(int argc, char** argv) {
    Kokkos::initialize(argc, argv);
    {
        ExecSpace exec;
        const int IJ = 16641, KD = 80, LALL = 10, BLOCKS = 8, REPEATS = 30;
        printf("=== PURE PORTABLE KOKKOS 4.7 FULL-SCALE NICAM NSW6 BENCHMARK ===\n");
        printf("Execution space: %s | Grid: %d x %d x %d = %.2fM cells | dt=30s | strict IEEE FP64\n",
               ExecSpace::name(), LALL, IJ, KD, (double)LALL * IJ * KD / 1e6);

        NSW6Params P1 = build_params(IJ, KD, 1, 1);
        NSW6Params P2 = build_params(IJ, KD, 2, 1);
        NSW6Params P3 = build_params(IJ, KD, 3, 1);
        auto make_dP = [&](const NSW6Params& P) {
            Kokkos::View<NSW6Params, MemSpace> d("dP");
            auto h = Kokkos::create_mirror_view(d);
            h() = P;
            Kokkos::deep_copy(d, h);
            return d;
        };
        auto dP1 = make_dP(P1), dP2 = make_dP(P2), dP3 = make_dP(P3);

        HostState H; H.init(IJ, KD, LALL, P1);
        KokkosState KS; KS.init(H);
        KokkosWorkspace WK; WK.init(IJ, KD, LALL);

        auto run_cfg = [&](int c) {
            if (c == 0)      run_kokkos_mode0(exec, dP1.data(), IJ, KD, LALL, 1, KS.dp, WK.dp);
            else if (c == 1) run_kokkos_mode1(exec, dP1.data(), IJ, KD, LALL, 1, KS.dp, WK.dp);
            else if (c == 2) run_kokkos_math<2>(exec, dP2.data(), IJ, KD, LALL, 1, KS.dp, WK);
            else if (c == 3) run_kokkos_math<2>(exec, dP2.data(), IJ, KD, LALL, 0, KS.dp, WK);
            else             run_kokkos_math<3>(exec, dP3.data(), IJ, KD, LALL, 1, KS.dp, WK);
        };

        for (int c = 0; c < 5; ++c) { KS.restore(exec); run_cfg(c); exec.fence(); }

        // 1. Single-step exact validation on region 0
        printf("\n=== SINGLE-STEP VALIDATION (Region 0: 16,641 x 80, vs Kokkos Mode 0) ===\n");
        size_t n3r = (size_t)IJ * KD, n2r = IJ;
        std::vector<double> ref_tem(n3r), ref_p(2 * n2r), ref_q(6 * n3r);
        auto h_tem = Kokkos::create_mirror_view(KS.tem);
        auto h_p   = Kokkos::create_mirror_view(KS.precip);
        auto h_q   = Kokkos::create_mirror_view(KS.q);
        for (int c : {0, 1, 2, 4}) {
            KS.restore(exec); run_cfg(c); exec.fence();
            Kokkos::deep_copy(h_tem, KS.tem);
            Kokkos::deep_copy(h_p,   KS.precip);
            Kokkos::deep_copy(h_q,   KS.q);
            if (c == 0) {
                std::copy(h_tem.data(), h_tem.data() + n3r, ref_tem.begin());
                std::copy(h_p.data(),   h_p.data()   + 2*n2r, ref_p.begin());
                std::copy(h_q.data(),   h_q.data()   + 6*n3r, ref_q.begin());
                continue;
            }
            double max_dT = 0, max_dpl = 0, max_dps = 0, max_dq = 0;
            long long n_fail = 0, n_neg = 0;
            for (size_t i = 0; i < n3r; ++i) {
                int k = (int)(i / IJ); if (k < 1 || k > KD - 2) continue;
                double dT = std::fabs(h_tem(i) - ref_tem[i]);
                if (dT > max_dT) max_dT = dT;
                if (dT > 1e-10 && dT > 1e-9 * std::fabs(ref_tem[i])) n_fail++;
                for (int s = 0; s < 6; ++s) {
                    double v = h_q(i + s * n3r), r = ref_q[i + s * n3r], dq = std::fabs(v - r);
                    if (dq > max_dq) max_dq = dq;
                    if (v < 0.0) n_neg++;
                    if (dq > 1e-10 && dq > 1e-9 * std::fabs(r)) n_fail++;
                }
            }
            for (size_t i = 0; i < n2r; ++i) {
                double dl = std::fabs(h_p(i) - ref_p[i]); if (dl > max_dpl) max_dpl = dl;
                double ds = std::fabs(h_p(i + n2r) - ref_p[i + n2r]); if (ds > max_dps) max_dps = ds;
            }
            const char* nm = (c == 1) ? "Kokkos-ENG " : (c == 2 ? "Kokkos-MATH" : "Kokkos-MONO");
            printf("  %s: max|dT|=%.3e K | max|dq|=%.3e | max|dP_liq|=%.3e | max|dP_sol|=%.3e | neg=%lld | N_FAIL=%lld\n",
                   nm, max_dT, max_dq, max_dpl, max_dps, n_neg, n_fail);
        }

        // 2. 8-Block Cyclically Rotated Paired Timing (95% Student-t CI)
        const char* cn[5] = { "Kokkos-BB", "Kokkos-ENG", "Kokkos-MATH(rc)", "Kokkos-MATH(nr)", "Kokkos-MONO(rc)" };
        double blk[5][8] = {};
        for (int b = 0; b < BLOCKS; ++b) {
            for (int oi = 0; oi < 5; ++oi) {
                int c = (b + oi) % 5;
                int nr = (c == 0) ? 10 : REPEATS;
                double acc = 0.0;
                for (int r = 0; r < nr; ++r) {
                    KS.restore(exec);
                    Kokkos::Timer t;
                    run_cfg(c);
                    exec.fence();
                    acc += t.seconds() * 1000.0;
                }
                blk[c][b] = acc / nr;
            }
        }

        printf("\n=== PURE KOKKOS PER-BLOCK MEANS (ms/step, 8 paired blocks on H200) ===\n");
        printf("Block | %11s | %11s | %15s | %15s | %15s\n", cn[0], cn[1], cn[2], cn[3], cn[4]);
        for (int b = 0; b < BLOCKS; ++b) {
            printf("  B%d  | %11.4f | %11.4f | %15.4f | %15.4f | %15.4f\n",
                   b + 1, blk[0][b], blk[1][b], blk[2][b], blk[3][b], blk[4][b]);
        }

        const double t975 = 2.364624251;
        printf("\n=== PURE KOKKOS MEANS & 95%% CI (df=7) ===\n");
        double mean[5] = {};
        for (int c = 0; c < 5; ++c) {
            double m = 0, v = 0;
            for (int b = 0; b < 8; ++b) m += blk[c][b];
            m /= 8.0; mean[c] = m;
            for (int b = 0; b < 8; ++b) v += (blk[c][b] - m) * (blk[c][b] - m);
            double hw = t975 * std::sqrt(v / 56.0);
            printf("  %-16s : %8.4f +/- %6.4f ms  (vs Kokkos-BB: %5.2fx)\n", cn[c], m, hw, mean[0] / m);
        }
        printf("\n=== PURE KOKKOS PAIRED CONTRASTS (Same Kokkos Framework!) ===\n");
        auto contrast = [&](int a, int b) {
            double m = 0, v = 0;
            for (int i = 0; i < 8; ++i) m += (blk[a][i] - blk[b][i]);
            m /= 8.0;
            for (int i = 0; i < 8; ++i) { double d = (blk[a][i] - blk[b][i]) - m; v += d * d; }
            double hw = t975 * std::sqrt(v / 56.0);
            printf("  %-16s -> %-16s : %+8.4f +/- %6.4f ms  (speedup = %.3fx)\n",
                   cn[a], cn[b], m, hw, mean[a] / mean[b]);
        };
        contrast(0, 1);
        contrast(1, 2);
        contrast(0, 2);
        contrast(0, 3);
        contrast(0, 4);
        contrast(4, 2);
    }
    Kokkos::finalize();
    return 0;
}
