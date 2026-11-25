// =============================================================================
// bench_v5_vs_v4.cu
//
// Strict Head-to-Head Verification & Timing of Hopper SM90a v5 vs v4.1:
//   1. Exact bit-level / strict-tolerance verification of v5 against v4.1
//      across ALL 23 output arrays (N_FAIL = 0!)
//   2. Paired CUDA Graph timing of v4.1 vs v5 (Mode 2 w/ rceff, Mode 2 no-rceff,
//      Mode 3 w/ rceff) + per-kernel timing breakdown
// =============================================================================
#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include "nsw6_hopper_v5.cuh"

static NSW6Params build_params(int ijdim, int kdim, int opt_mode, int compute_rceff, int collect_cnt) {
    NSW6Params P = {};
    P.ijdim = ijdim; P.kdim = kdim; P.kmin0 = 1; P.kmax0 = kdim - 2; P.nqmax = 6;
    P.i_qv = 0; P.i_qc = 1; P.i_qr = 2; P.i_qi = 3; P.i_qs = 4; P.i_qg = 5;
    P.opt_mode = opt_mode; P.ein_is_simple = 1; P.prcip_kin2ein = 1; P.precip_4water = 1;
    P.opt_indir = 0; P.opt_sclw = 0; P.opt_explicit_icegen = 1; P.compute_rceff = compute_rceff;
    P.is_top_tem = 1; P.is_top_epl = 0; P.is_btm_tem = 1; P.is_btm_epl = 0;
    P.mp_div_num = 1; P.collect_counters = collect_cnt;
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

struct DeviceState {
    size_t n3, n4, n2;
    HopperStatePtrs cur, init;
    void alloc_and_upload(const HostState& H) {
        n3 = (size_t)H.ijdim * H.kdim * H.lall; n4 = n3 * H.nqmax; n2 = (size_t)H.ijdim * H.lall;
        auto mk = [&](double*& dst, const std::vector<double>& src) {
            CUDA_CHECK(cudaMalloc(&dst, src.size() * sizeof(double)));
            CUDA_CHECK(cudaMemcpy(dst, src.data(), src.size() * sizeof(double), cudaMemcpyHostToDevice));
        };
        auto ro = [&](const double*& dst, const std::vector<double>& src) { double* t=nullptr; mk(t,src); dst=t; };
        ro(cur.z,H.z); ro(cur.zh,H.zh); ro(cur.phi,H.phi); ro(cur.gam2,H.gam2); ro(cur.gam2h,H.gam2h);
        ro(cur.gsgam2,H.gsgam2); ro(cur.gsgam2h,H.gsgam2h);
        ro(cur.ix,H.ix); ro(cur.iy,H.iy); ro(cur.iz,H.iz); ro(cur.jx,H.jx); ro(cur.jy,H.jy); ro(cur.jz,H.jz);
        ro(cur.unccn,H.unccn); ro(cur.factor_n_aie,H.factor_n_aie);
        mk(init.rhog,H.rhog); mk(init.rhogvx,H.rhogvx); mk(init.rhogvy,H.rhogvy); mk(init.rhogvz,H.rhogvz);
        mk(init.rhogw,H.rhogw); mk(init.rhoge,H.rhoge); mk(init.rhogq,H.rhogq);
        mk(init.rho,H.rho); mk(init.pre,H.pre); mk(init.tem,H.tem);
        mk(init.vx,H.vx); mk(init.vy,H.vy); mk(init.vz,H.vz); mk(init.w,H.w); mk(init.q,H.q);
        mk(cur.rhog,H.rhog); mk(cur.rhogvx,H.rhogvx); mk(cur.rhogvy,H.rhogvy); mk(cur.rhogvz,H.rhogvz);
        mk(cur.rhogw,H.rhogw); mk(cur.rhoge,H.rhoge); mk(cur.rhogq,H.rhogq);
        mk(cur.rho,H.rho); mk(cur.pre,H.pre); mk(cur.tem,H.tem);
        mk(cur.vx,H.vx); mk(cur.vy,H.vy); mk(cur.vz,H.vz); mk(cur.w,H.w); mk(cur.q,H.q);
        mk(cur.precip,H.precip); mk(cur.precip_lh,H.precip_lh); mk(cur.precip_rhoein,H.precip_rhoein);
        mk(cur.precip_rhophi,H.precip_rhophi); mk(cur.precip_rhokin,H.precip_rhokin); mk(cur.precip_energy,H.precip_energy);
        mk(cur.gdcfrc,H.gdcfrc); mk(cur.gdclw,H.gdclw);
    }
    void restore(cudaStream_t s) {
        auto cp = [&](double* d, const double* src, size_t n) { CUDA_CHECK(cudaMemcpyAsync(d,src,n*sizeof(double),cudaMemcpyDeviceToDevice,s)); };
        cp(cur.rhog,init.rhog,n3); cp(cur.rhogvx,init.rhogvx,n3); cp(cur.rhogvy,init.rhogvy,n3); cp(cur.rhogvz,init.rhogvz,n3);
        cp(cur.rhogw,init.rhogw,n3); cp(cur.rhoge,init.rhoge,n3); cp(cur.rhogq,init.rhogq,n4);
        cp(cur.rho,init.rho,n3); cp(cur.pre,init.pre,n3); cp(cur.tem,init.tem,n3);
        cp(cur.vx,init.vx,n3); cp(cur.vy,init.vy,n3); cp(cur.vz,init.vz,n3); cp(cur.w,init.w,n3); cp(cur.q,init.q,n4);
    }
};

int main() {
    const int IJ = 16641, KD = 80, LALL = 10;
    NSW6Params P2r = build_params(IJ, KD, 2, 1, 0);
    HostState H; H.init(IJ, KD, LALL, P2r);
    DeviceState D; D.alloc_and_upload(H);
    HopperWorkspaceV3 W4;
    HopperWorkspaceV5 W5;
    cudaStream_t strm; CUDA_CHECK(cudaStreamCreate(&strm));
    sm90_v4_upload_state_ptrs(D.cur, strm);

    // 1. Exact bit-for-bit verification of v5 against v4.1 across ALL 23 arrays!
    printf("=== EXACT VERIFICATION: v5 vs v4.1 (Full 13.31M cells, Mode 2 & Mode 3) ===\n");
    for (int mode : {2, 3}) {
        NSW6Params Pm = build_params(IJ, KD, mode, 1, 0);
        sm90_v2_upload_constants(Pm, H.rdgz.data(), H.rdgzh.data(), H.afact.data(), H.bfact.data(), H.cfact.data(), H.dfact.data());
        // Run v4.1
        D.restore(strm);
        sm90_v4_launch_step(Pm, LALL, D.cur, W4, strm);
        CUDA_CHECK(cudaStreamSynchronize(strm));

        auto dl = [&](const double* d, size_t n) { std::vector<double> v(n); CUDA_CHECK(cudaMemcpy(v.data(), d, n*8, cudaMemcpyDeviceToHost)); return v; };
        std::vector<std::pair<const char*, std::vector<double>>> ref = {
            {"rhog", dl(D.cur.rhog, D.n3)}, {"rhogvx", dl(D.cur.rhogvx, D.n3)}, {"rhogvy", dl(D.cur.rhogvy, D.n3)},
            {"rhogvz", dl(D.cur.rhogvz, D.n3)}, {"rhogw", dl(D.cur.rhogw, D.n3)}, {"rhoge", dl(D.cur.rhoge, D.n3)},
            {"rhogq", dl(D.cur.rhogq, D.n4)}, {"rho", dl(D.cur.rho, D.n3)}, {"pre", dl(D.cur.pre, D.n3)},
            {"tem", dl(D.cur.tem, D.n3)}, {"vx", dl(D.cur.vx, D.n3)}, {"vy", dl(D.cur.vy, D.n3)},
            {"vz", dl(D.cur.vz, D.n3)}, {"w", dl(D.cur.w, D.n3)}, {"q", dl(D.cur.q, D.n4)},
            {"precip", dl(D.cur.precip, 2*D.n2)}, {"precip_lh", dl(D.cur.precip_lh, D.n2)},
            {"precip_rhoein", dl(D.cur.precip_rhoein, D.n2)}, {"precip_rhophi", dl(D.cur.precip_rhophi, D.n2)},
            {"precip_rhokin", dl(D.cur.precip_rhokin, D.n2)}, {"precip_energy", dl(D.cur.precip_energy, D.n2)},
            {"gdcfrc", dl(D.cur.gdcfrc, D.n3)}, {"rceff", dl(W4.d_rceff, D.n3)}, {"rctop", dl(W4.d_rctop, D.n2)}
        };

        // Run v5 (twice with dirty Vt to prove no Vt-memset bug!)
        D.restore(strm);
        sm90_v5_launch_step(Pm, LALL, D.cur, W5, strm);
        D.restore(strm);
        sm90_v5_launch_step(Pm, LALL, D.cur, W5, strm);
        CUDA_CHECK(cudaStreamSynchronize(strm));

        std::vector<const double*> ptrs = {
            D.cur.rhog, D.cur.rhogvx, D.cur.rhogvy, D.cur.rhogvz, D.cur.rhogw, D.cur.rhoge,
            D.cur.rhogq, D.cur.rho, D.cur.pre, D.cur.tem, D.cur.vx, D.cur.vy, D.cur.vz, D.cur.w, D.cur.q,
            D.cur.precip, D.cur.precip_lh, D.cur.precip_rhoein, D.cur.precip_rhophi, D.cur.precip_rhokin,
            D.cur.precip_energy, D.cur.gdcfrc, W5.base.d_rceff, W5.base.d_rctop
        };
        double worst_abs = 0.0; long long n_fail = 0; const char* worst_arr = "none";
        for (size_t a = 0; a < ref.size(); ++a) {
            auto got = dl(ptrs[a], ref[a].second.size());
            for (size_t i = 0; i < got.size(); ++i) {
                double r = ref[a].second[i], g = got[i], diff = std::fabs(g - r);
                if (diff > worst_abs) { worst_abs = diff; worst_arr = ref[a].first; }
                if (std::isnan(g) || (diff > 1e-12 && diff > 1e-11 * std::fabs(r))) n_fail++;
            }
        }
        printf("  Mode %d (all 24 arrays x 13.31M cells): max|v5 - v4.1| = %.3e (%s) | N_FAIL = %lld\n",
               mode, worst_abs, worst_arr, n_fail);
    }

    // 2. Capture CUDA Graphs for v4.1 vs v5 and run 8 paired blocks!
    cudaGraph_t g[6]; cudaGraphExec_t ge[6];
    auto cap = [&](int idx, bool is_v5, int mode, int rc) {
        NSW6Params Pm = build_params(IJ, KD, mode, rc, 0);
        sm90_v2_upload_constants(Pm, H.rdgz.data(), H.rdgzh.data(), H.afact.data(), H.bfact.data(), H.cfact.data(), H.dfact.data());
        D.restore(strm);
        if (is_v5) sm90_v5_launch_step(Pm, LALL, D.cur, W5, strm);
        else       sm90_v4_launch_step(Pm, LALL, D.cur, W4, strm);
        CUDA_CHECK(cudaStreamSynchronize(strm));
        CUDA_CHECK(cudaStreamBeginCapture(strm, cudaStreamCaptureModeGlobal));
        if (is_v5) sm90_v5_launch_step(Pm, LALL, D.cur, W5, strm);
        else       sm90_v4_launch_step(Pm, LALL, D.cur, W4, strm);
        CUDA_CHECK(cudaStreamEndCapture(strm, &g[idx]));
        CUDA_CHECK(cudaGraphInstantiate(&ge[idx], g[idx], 0));
    };
    cap(0, false, 2, 1); cap(1, true, 2, 1);
    cap(2, false, 2, 0); cap(3, true, 2, 0);
    cap(4, false, 3, 1); cap(5, true, 3, 1);

    const char* nm[6] = {"v4.1-M2(rc)", "v5-M2(rc)", "v4.1-M2(nr)", "v5-M2(nr)", "v4.1-M3(rc)", "v5-M3(rc)"};
    const int md[6] = {2, 2, 2, 2, 3, 3};
    const int rc[6] = {1, 1, 0, 0, 1, 1};
    cudaEvent_t e0, e1; CUDA_CHECK(cudaEventCreate(&e0)); CUDA_CHECK(cudaEventCreate(&e1));

    const int BLOCKS = 8, REPS = 30;
    double blk[6][8] = {}, best[6];
    for (int i = 0; i < 6; ++i) best[i] = 1e9;

    for (int b = 0; b < BLOCKS; ++b) {
        for (int oi = 0; oi < 6; ++oi) {
            int c = (b + oi) % 6;
            NSW6Params Pm = build_params(IJ, KD, md[c], rc[c], 0);
            sm90_v2_upload_constants(Pm, H.rdgz.data(), H.rdgzh.data(), H.afact.data(), H.bfact.data(), H.cfact.data(), H.dfact.data());
            double acc = 0.0;
            for (int r = 0; r < REPS; ++r) {
                D.restore(strm);
                CUDA_CHECK(cudaEventRecord(e0, strm));
                CUDA_CHECK(cudaGraphLaunch(ge[c], strm));
                CUDA_CHECK(cudaEventRecord(e1, strm));
                CUDA_CHECK(cudaEventSynchronize(e1));
                float ms = 0; CUDA_CHECK(cudaEventElapsedTime(&ms, e0, e1));
                acc += ms;
                if (ms < best[c]) best[c] = ms;
            }
            blk[c][b] = acc / REPS;
        }
    }

    printf("\n=== PAIRED CUDA GRAPH TIMING: v5 vs v4.1 (8 blocks x 30 reps on H200) ===\n");
    for (int c = 0; c < 6; ++c) {
        double m = 0, v = 0;
        for (int b = 0; b < 8; ++b) m += blk[c][b];
        m /= 8.0;
        for (int b = 0; b < 8; ++b) v += (blk[c][b] - m) * (blk[c][b] - m);
        double hw = 2.364624251 * std::sqrt(v / 56.0);
        printf("  %-12s : mean = %7.4f +/- %6.4f ms | min = %7.4f ms\n", nm[c], m, hw, best[c]);
    }
    for (int p = 0; p < 6; p += 2) {
        double ma = 0, mb = 0;
        for (int b = 0; b < 8; ++b) { ma += blk[p][b]; mb += blk[p+1][b]; }
        ma /= 8.0; mb /= 8.0;
        printf("  %s -> %s : mean %.4f -> %.4f ms (%.3fx) | min %.4f -> %.4f ms (%.3fx, saved %.3f ms)\n",
               nm[p], nm[p+1], ma, mb, ma/mb, best[p], best[p+1], best[p]/best[p+1], best[p] - best[p+1]);
    }
    return 0;
}
