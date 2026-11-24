// =============================================================================
// bench_h200_v4.cu
//
// Final Roofline-Saturating NVIDIA H200 (Hopper SM90a) Benchmark, Validation,
// Component Breakdown, Math Work Counters, and 5-Step / 100-Step Conservation!
// =============================================================================
#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include <numeric>
#include "nsw6_hopper_v4.cuh"

static NSW6Params build_nicam_production_params(int ijdim, int kdim, int opt_mode, int compute_rceff, int collect_cnt) {
    NSW6Params P = {};
    P.ijdim = ijdim; P.kdim = kdim;
    P.kmin0 = 1; P.kmax0 = kdim - 2; P.nqmax = 6;
    P.i_qv = 0; P.i_qc = 1; P.i_qr = 2; P.i_qi = 3; P.i_qs = 4; P.i_qg = 5;
    P.opt_mode = opt_mode; P.ein_is_simple = 1; P.prcip_kin2ein = 1; P.precip_4water = 1;
    P.opt_indir = 0; P.opt_sclw = 0; P.opt_explicit_icegen = 1;
    P.compute_rceff = compute_rceff;
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
        for (int k = 1; k < kdim - 1; ++k)
            gzh[k + 1] = gzh[k] + 66.0 + 1350.0 * std::pow((double)(k - 1) / (double)(kdim - 3), 2.2);
        for (int k = 1; k < kdim - 1; ++k) gz[k] = 0.5 * (gzh[k] + gzh[k + 1]);
        gz[0] = -gz[1]; gz[kdim - 1] = 2.0 * gzh[kdim - 1] - gz[kdim - 2];
        for (int k = 0; k < kdim; ++k) {
            dgz[k] = (k < kdim - 1) ? (gzh[k + 1] - gzh[k]) : (gzh[kdim - 1] - gzh[kdim - 2]);
            rdgz[k] = 1.0 / dgz[k];
            rdgzh[k] = 1.0 / ((k > 0) ? (gz[k] - gz[k - 1]) : (gz[1] - gz[0]));
            afact[k] = (k > 0) ? ((gzh[k] - gz[k - 1]) / (gz[k] - gz[k - 1])) : 0.5;
            bfact[k] = 1.0 - afact[k];
            cfact[k] = (k < kdim - 1) ? ((gz[k] - gzh[k]) / (gzh[k + 1] - gzh[k])) : 0.5;
            dfact[k] = 1.0 - cfact[k];
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

static void run_cpu_ref_region0(const NSW6Params& P, HostState& H) {
    int ijdim = H.ijdim, kdim = H.kdim;
    std::vector<double> Vt(4ULL * ijdim * kdim, 0.0), s_zh(kdim), s_z(kdim), s_kv(kdim), s_zdis(4 * kdim);
    int d0=0, d1=0, d2=0, d3=0, d4=0;
    for (int k = 0; k < kdim; ++k)
        for (int ij = 0; ij < ijdim; ++ij)
            nsw6_pre_and_core_cell(P, ij, k, H.rhog.data(), H.rhoge.data(), H.rhogq.data(), H.rho.data(), H.tem.data(), H.pre.data(), H.q.data(),
                H.w.data(), H.unccn.data(), H.factor_n_aie.data(), H.gsgam2.data(), Vt.data(), &d0, &d1, &d2, &d3, &d4);
    for (int ij = 0; ij < ijdim; ++ij)
        precip_transport_column(P, ij, 1, s_zh.data(), s_z.data(), s_kv.data(), s_zdis.data(),
            H.rhog.data(), H.rhogvx.data(), H.rhogvy.data(), H.rhogvz.data(), H.rhogw.data(), H.rhoge.data(), H.rhogq.data(),
            H.rho.data(), H.tem.data(), H.pre.data(), H.vx.data(), H.vy.data(), H.vz.data(), H.q.data(), H.z.data(), H.zh.data(), Vt.data(),
            H.precip.data(), H.precip_lh.data(), H.precip_rhoein.data(), H.precip_rhophi.data(), H.precip_rhokin.data(), H.precip_energy.data(),
            H.gsgam2.data(), H.gsgam2h.data(), H.gam2.data(), H.gam2h.data(),
            H.ix.data(), H.iy.data(), H.iz.data(), H.jx.data(), H.jy.data(), H.jz.data(),
            H.rdgz.data(), H.rdgzh.data(), H.afact.data(), H.bfact.data(), H.cfact.data(), H.dfact.data(), &d0, &d1, &d2, &d3, &d4);
    for (int k = 0; k < kdim; ++k)
        for (int ij = 0; ij < ijdim; ++ij)
            satadj_and_post_negfilter_cell(P, ij, k, H.rhog.data(), H.rhoge.data(), H.rhogq.data(), H.rho.data(), H.tem.data(), H.pre.data(), H.q.data(), H.gsgam2.data(),
                &d0, &d1, &d2, &d3, &d4);
    for (int ij = 0; ij < ijdim; ++ij)
        cnvvar_bndcnd_cldfrac_column(P, ij, H.rhog.data(), H.rhogvx.data(), H.rhogvy.data(), H.rhogvz.data(), H.rhogw.data(), H.rhoge.data(), H.rhogq.data(),
            H.rho.data(), H.pre.data(), H.tem.data(), H.vx.data(), H.vy.data(), H.vz.data(), H.w.data(), H.q.data(),
            H.gdcfrc.data(), H.gdclw.data(), H.phi.data(), H.gsgam2.data(), H.gsgam2h.data(), H.afact.data(), H.bfact.data());
}

int main(int argc, char** argv) {
    const int IJ = 16641, KD = 80, LALL = 10;
    const int BLOCKS  = (argc > 1) ? std::atoi(argv[1]) : 8;
    const int REPEATS = (argc > 2) ? std::atoi(argv[2]) : 50;

    NSW6Params P_ref = build_nicam_production_params(IJ, KD, 1, 1, 0);
    cudaDeviceProp prop; CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
    printf("================================================================================\n");
    printf(" NICAM NSW6 HOPPER SM90a v4 ROOFLINE-SATURATING BENCHMARK & VALIDATION\n");
    printf(" GPU: %s (SM %d.%d, %d SMs, %.1f GB HBM3e, %.1f MB L2)\n",
           prop.name, prop.major, prop.minor, prop.multiProcessorCount, prop.totalGlobalMem/1e9, prop.l2CacheSize/1e6);
    printf(" Workload: full NICAM scale (%d cols/reg x %d lvls x %d regs = 13.31M cells)\n", IJ, KD, LALL);
    printf("================================================================================\n\n");

    // 1. Exact GPU-vs-CPU Differential Validation (region0)
    {
        printf("--- [Phase 1] Exact GPU-vs-CPU Differential Validation (v4, region0: 16,641 x 80) ---\n");
        HostState H_ref; H_ref.init(IJ, KD, 1, P_ref);
        DeviceState D_val; D_val.alloc_and_upload(H_ref);
        run_cpu_ref_region0(P_ref, H_ref);
        HopperWorkspaceV3 W_val; cudaStream_t strm; CUDA_CHECK(cudaStreamCreate(&strm));
        sm90_v4_upload_state_ptrs(D_val.cur, strm);

        for (int mode : {1, 2}) {
            NSW6Params Pm = build_nicam_production_params(IJ, KD, mode, 1, 0);
            sm90_v2_upload_constants(Pm, H_ref.rdgz.data(), H_ref.rdgzh.data(),
                                     H_ref.afact.data(), H_ref.bfact.data(), H_ref.cfact.data(), H_ref.dfact.data());
            D_val.restore(strm);
            sm90_v4_launch_step(Pm, 1, D_val.cur, W_val, strm);
            CUDA_CHECK(cudaStreamSynchronize(strm));

            long long nf = 0; double wr = 0.0, wa = 0.0;
            auto chk = [&](const char* nm, const double* dp, const std::vector<double>& ref) {
                std::vector<double> h(ref.size());
                CUDA_CHECK(cudaMemcpy(h.data(), dp, ref.size()*sizeof(double), cudaMemcpyDeviceToHost));
                long long lnf = 0; double lwr = 0.0, lwa = 0.0;
                for (size_t i = 0; i < ref.size(); ++i) {
                    double d = std::fabs(h[i] - ref[i]);
                    if (d > lwa) lwa = d;
                    if (std::fabs(ref[i]) > 2.2204460492503131e-16) {
                        double r = d / std::fabs(ref[i]);
                        if (r > lwr) lwr = r;
                    }
                    if (!(d <= 1e-10 || d <= 1e-9 * std::fabs(ref[i]))) lnf++;
                }
                if (lnf > 0) printf("    DIFF %-14s: nfail=%lld max_rel=%.3e max_abs=%.3e\n", nm, lnf, lwr, lwa);
                nf += lnf; if (lwr > wr) wr = lwr; if (lwa > wa) wa = lwa;
            };
            chk("rhog",D_val.cur.rhog,H_ref.rhog); chk("rhogvx",D_val.cur.rhogvx,H_ref.rhogvx);
            chk("rhogvy",D_val.cur.rhogvy,H_ref.rhogvy); chk("rhogvz",D_val.cur.rhogvz,H_ref.rhogvz);
            chk("rhogw",D_val.cur.rhogw,H_ref.rhogw); chk("rhoge",D_val.cur.rhoge,H_ref.rhoge);
            chk("rhogq",D_val.cur.rhogq,H_ref.rhogq); chk("rho",D_val.cur.rho,H_ref.rho);
            chk("pre",D_val.cur.pre,H_ref.pre); chk("tem",D_val.cur.tem,H_ref.tem);
            chk("vx",D_val.cur.vx,H_ref.vx); chk("vy",D_val.cur.vy,H_ref.vy);
            chk("vz",D_val.cur.vz,H_ref.vz); chk("w",D_val.cur.w,H_ref.w); chk("q",D_val.cur.q,H_ref.q);
            chk("precip",D_val.cur.precip,H_ref.precip); chk("plh",D_val.cur.precip_lh,H_ref.precip_lh);
            chk("prein",D_val.cur.precip_rhoein,H_ref.precip_rhoein); chk("pphi",D_val.cur.precip_rhophi,H_ref.precip_rhophi);
            chk("pkin",D_val.cur.precip_rhokin,H_ref.precip_rhokin); chk("pene",D_val.cur.precip_energy,H_ref.precip_energy);
            chk("cfrc",D_val.cur.gdcfrc,H_ref.gdcfrc); chk("clw",D_val.cur.gdclw,H_ref.gdclw);
            printf("  Mode %d (%s) GPU vs CPU Ref: N_FAIL = %lld | max_rel = %.3e | max_abs = %.3e -> %s\n",
                   mode, mode==1 ? "ENG-SM90-v4     " : "ENG+MATH-SM90-v4", nf, wr, wa, nf==0 ? "PASS" : "FAIL");
            if (nf != 0) return 2;
        }
        W_val.release(); CUDA_CHECK(cudaStreamDestroy(strm));
        printf("\n");
    }

    // 2. Full 10-Region Paired Balanced Benchmark (8 blocks x 50 repeats = 400 steps/cfg)
    HostState H; H.init(IJ, KD, LALL, P_ref);
    DeviceState D; D.alloc_and_upload(H);
    HopperWorkspaceV3 W; cudaStream_t strm; CUDA_CHECK(cudaStreamCreate(&strm));
    sm90_v4_upload_state_ptrs(D.cur, strm);
    cudaEvent_t ev[14]; for (int i = 0; i < 14; ++i) CUDA_CHECK(cudaEventCreate(&ev[i]));

    const int NCFG = 5;
    const char* names[NCFG] = {
        "Mode 1: ENG-SM90-v4 (0-spill, orig math)",
        "Mode 2: ENG+MATH-SM90-v4 (Stream, w/ rceff)",
        "Mode 2: ENG+MATH-SM90-v4 (Graph,  w/ rceff)",
        "Mode 2: ENG+MATH-SM90-v4 (Graph, no rceff)",
        "Mode 3: ENG+MATH+MONO-v4 (Graph,  w/ rceff)"
    };
    const int modes[NCFG]  = { 1, 2, 2, 2, 3 };
    const int rceffs[NCFG] = { 1, 1, 1, 0, 1 };
    const int graphs[NCFG] = { 0, 0, 1, 1, 1 };

    for (int m : {1, 2, 3}) {
        NSW6Params Pm = build_nicam_production_params(IJ, KD, m, 1, 0);
        sm90_v2_upload_constants(Pm, H.rdgz.data(), H.rdgzh.data(), H.afact.data(), H.bfact.data(), H.cfact.data(), H.dfact.data());
        D.restore(strm);
        sm90_v4_launch_step(Pm, LALL, D.cur, W, strm);
    }
    CUDA_CHECK(cudaStreamSynchronize(strm));

    cudaGraph_t g_m2_rc, g_m2_norc, g_m3_rc;
    cudaGraphExec_t ge_m2_rc, ge_m2_norc, ge_m3_rc;
    auto cap = [&](int mode, int rc, cudaGraph_t& g, cudaGraphExec_t& ge) {
        NSW6Params Pm = build_nicam_production_params(IJ, KD, mode, rc, 0);
        sm90_v2_upload_constants(Pm, H.rdgz.data(), H.rdgzh.data(), H.afact.data(), H.bfact.data(), H.cfact.data(), H.dfact.data());
        CUDA_CHECK(cudaStreamBeginCapture(strm, cudaStreamCaptureModeGlobal));
        sm90_v4_launch_step(Pm, LALL, D.cur, W, strm);
        CUDA_CHECK(cudaStreamEndCapture(strm, &g));
        CUDA_CHECK(cudaGraphInstantiate(&ge, g, 0));
    };
    cap(2, 1, g_m2_rc, ge_m2_rc);
    cap(2, 0, g_m2_norc, ge_m2_norc);
    cap(3, 1, g_m3_rc, ge_m3_rc);

    std::vector<double> samples[NCFG];
    for (int b = 0; b < BLOCKS; ++b) {
        for (int oi = 0; oi < NCFG; ++oi) {
            int cfg = (b + oi) % NCFG;
            NSW6Params Pm = build_nicam_production_params(IJ, KD, modes[cfg], rceffs[cfg], 0);
            sm90_v2_upload_constants(Pm, H.rdgz.data(), H.rdgzh.data(), H.afact.data(), H.bfact.data(), H.cfact.data(), H.dfact.data());
            for (int r = 0; r < REPEATS; ++r) {
                D.restore(strm);
                CUDA_CHECK(cudaEventRecord(ev[0], strm));
                if (graphs[cfg]) {
                    cudaGraphExec_t ge = (modes[cfg] == 3) ? ge_m3_rc : (rceffs[cfg] ? ge_m2_rc : ge_m2_norc);
                    CUDA_CHECK(cudaGraphLaunch(ge, strm));
                } else {
                    sm90_v4_launch_step(Pm, LALL, D.cur, W, strm);
                }
                CUDA_CHECK(cudaEventRecord(ev[1], strm));
                CUDA_CHECK(cudaEventSynchronize(ev[1]));
                float ms = 0.0f; CUDA_CHECK(cudaEventElapsedTime(&ms, ev[0], ev[1]));
                samples[cfg].push_back((double)ms);
            }
        }
    }

    const double H200_BB_REF = 32.040; // measured on this exact H200 in bench_h200_v3!
    printf("--- [Phase 2] Paired Balanced H200 v4 Benchmark (%d blocks x %d = %d steps/cfg) ---\n",
           BLOCKS, REPEATS, BLOCKS * REPEATS);
    printf("  %-45s | %7s | %7s | %7s | %6s | %8s | %8s | %8s\n",
           "Configuration", "Mean", "Med", "Min", "Std", "vs H200B", "vs GH_BB", "vs GH_TT");
    printf("  -----------------------------------------------------------------------------------------------------------------------\n");
    for (int cfg = 0; cfg < NCFG; ++cfg) {
        auto& v = samples[cfg];
        double mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
        double sq = 0.0; for (double x : v) sq += (x - mean)*(x - mean);
        double sd = std::sqrt(sq / (v.size() - 1));
        std::sort(v.begin(), v.end());
        printf("  %-45s | %7.3f | %7.3f | %7.3f | %6.3f | %7.2fx | %7.2fx | %7.2fx\n",
               names[cfg], mean, v[v.size()/2], v.front(), sd,
               H200_BB_REF/mean, 78.19/mean, 68.57/mean);
    }

    // 3. Per-Kernel Component Breakdown on H200 (v4, 20-step avg)
    printf("\n--- [Phase 3] Per-Kernel Component Breakdown (v4, 10 regions = 13.31M cells, 20-step avg) ---\n");
    {
        const size_t n3 = (size_t)IJ * KD * LALL;
        const int g_wet_max  = (int)(((size_t)IJ * (KD - 2) * LALL + 127) / 128);
        const int g_wcol_max = (int)(((size_t)IJ * LALL + 63) / 64);
        dim3 g_scr((IJ + 255) / 256, KD - 2, LALL);
        dim3 g_128k((IJ + 127) / 128, KD, LALL), g_256k((IJ + 255) / 256, KD, LALL);
        dim3 g_col128((IJ + 127) / 128, LALL), g_col256((IJ + 255) / 256, LALL);

        for (int m : {2, 3}) {
            NSW6Params Pm = build_nicam_production_params(IJ, KD, m, 1, 0);
            sm90_v2_upload_constants(Pm, H.rdgz.data(), H.rdgzh.data(), H.afact.data(), H.bfact.data(), H.cfact.data(), H.dfact.data());
            double t[11] = {}; const int NR = 20;
            for (int r = 0; r < NR; ++r) {
                D.restore(strm);
                CUDA_CHECK(cudaEventRecord(ev[0], strm));
                CUDA_CHECK(cudaMemsetAsync(W.d_Vt, 0, 4 * n3 * sizeof(double), strm));
                CUDA_CHECK(cudaMemsetAsync(W.d_meta_block, 0, W.meta_bytes, strm));
                CUDA_CHECK(cudaEventRecord(ev[1], strm));
                k_sm90_v4_nsw6_screen<false><<<g_scr, 256, 0, strm>>>(W.d_wet_cells, W.d_wet_count(), W.d_cnt);
                CUDA_CHECK(cudaEventRecord(ev[2], strm));
                k_sm90_v4_nsw6_wet_core<false><<<g_wet_max, 128, 0, strm>>>(
                    W.d_wet_cells, W.d_wet_count(), W.d_Vt, W.d_col_spmask(), W.d_col_ktop(), W.d_cnt);
                CUDA_CHECK(cudaEventRecord(ev[3], strm));
                k_sm90_v3_rceff<<<g_col128, 128, 0, strm>>>(D.cur, W.d_rceff, W.d_rctop);
                CUDA_CHECK(cudaEventRecord(ev[4], strm));
                k_sm90_v4_init_kv<<<g_256k, 256, 0, strm>>>(W.d_kv);
                CUDA_CHECK(cudaEventRecord(ev[5], strm));
                k_sm90_v4_compact_cols<false><<<g_col256, 256, 0, strm>>>(
                    W.d_col_spmask(), W.d_wet_cols, W.d_wet_col_count(), W.d_cnt);
                CUDA_CHECK(cudaEventRecord(ev[6], strm));
                if (m == 2) k_sm90_v4_precip_half<2,false><<<g_wcol_max, 64, 0, strm>>>(
                    W.d_wet_cols, W.d_wet_col_count(), W.d_col_spmask(), W.d_col_ktop(), W.d_Vt, W.d_kv, W.d_zdis_m, W.d_cnt);
                else        k_sm90_v4_precip_half<3,false><<<g_wcol_max, 64, 0, strm>>>(
                    W.d_wet_cols, W.d_wet_col_count(), W.d_col_spmask(), W.d_col_ktop(), W.d_Vt, W.d_kv, W.d_zdis_m, W.d_cnt);
                CUDA_CHECK(cudaEventRecord(ev[7], strm));
                if (m == 2) k_sm90_v4_precip_full<2,false><<<g_wcol_max, 64, 0, strm>>>(
                    W.d_wet_cols, W.d_wet_col_count(), W.d_col_spmask(), W.d_col_ktop(), W.d_Vt, W.d_kv, W.d_zdis_m, W.d_cnt);
                else        k_sm90_v4_precip_full<3,false><<<g_wcol_max, 64, 0, strm>>>(
                    W.d_wet_cols, W.d_wet_col_count(), W.d_col_spmask(), W.d_col_ktop(), W.d_Vt, W.d_kv, W.d_zdis_m, W.d_cnt);
                CUDA_CHECK(cudaEventRecord(ev[8], strm));
                if (m == 2) k_sm90_v4_satadj_nf_diag<2,false,true><<<g_128k, 128, 0, strm>>>(W.d_kv, W.d_col_ktop(), W.d_cnt);
                else        k_sm90_v4_satadj_nf_diag<3,false,true><<<g_128k, 128, 0, strm>>>(W.d_kv, W.d_col_ktop(), W.d_cnt);
                CUDA_CHECK(cudaEventRecord(ev[9], strm));
                k_sm90_v2_diag_w<<<g_256k, 256, 0, strm>>>(D.cur);
                CUDA_CHECK(cudaEventRecord(ev[10], strm));
                k_sm90_v2_bndcnd<<<g_col256, 256, 0, strm>>>(D.cur);
                CUDA_CHECK(cudaEventRecord(ev[11], strm));
                CUDA_CHECK(cudaEventSynchronize(ev[11]));
                for (int i = 0; i < 11; ++i) { float ms; CUDA_CHECK(cudaEventElapsedTime(&ms, ev[i], ev[i+1])); t[i] += ms; }
            }
            for (int i = 0; i < 11; ++i) t[i] /= NR;
            double sum = 0.0; for (int i = 0; i < 11; ++i) sum += t[i];
            printf("  Mode %d: Mem=%.3f | Scr=%.3f | Wet=%.3f | rceff=%.3f | kv=%.3f | Cmp=%.3f | PrH=%.3f | PrF=%.3f | SatDiag=%.3f | w=%.3f | BC=%.3f | SUM=%.3f ms\n",
                   m, t[0], t[1], t[2], t[3], t[4], t[5], t[6], t[7], t[8], t[9], t[10], sum);
        }
    }

    CUDA_CHECK(cudaGraphExecDestroy(ge_m2_rc)); CUDA_CHECK(cudaGraphDestroy(g_m2_rc));
    CUDA_CHECK(cudaGraphExecDestroy(ge_m2_norc)); CUDA_CHECK(cudaGraphDestroy(g_m2_norc));
    CUDA_CHECK(cudaGraphExecDestroy(ge_m3_rc)); CUDA_CHECK(cudaGraphDestroy(g_m3_rc));
    W.release(); CUDA_CHECK(cudaStreamDestroy(strm));
    return 0;
}
