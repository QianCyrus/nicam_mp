// =============================================================================
// bench_h200_v3.cu
//
// Complete NVIDIA H200 (Hopper SM90a) Full-Scale Benchmark, Validation,
// Apples-to-Apples H200 Unfused Baseline Comparator, Component Breakdown,
// Work Counters, and 5-Step / 100-Step Multi-Step Conservation Experiment!
//
// Workload: Production NICAM `full` scale:
//   16,641 cols/region x 80 levels x 10 regions = 166,410 cols / 13.31M cells
//   Strict IEEE FP64, dt = 30.0 s, 6 water tracers, Roh_flag=.true., KIN2EIN
// =============================================================================
#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include <algorithm>
#include <numeric>
#include "nsw6_hopper_v3.cuh"

// =============================================================================
// Faithful H200 Mode-0 Unfused OpenACC-Style Reference Comparator
//   Mirrors NICAM mod_mp_driver.f90's exact region-by-region (l=0..9) unfused
//   stage structure: pre_negfilter -> mp_nsw6_core (16 trans Roh, no skip) ->
//   rceff (k-sequential loop) -> precip_transport (per-species unfused loops,
//   explicit workspace arrays) -> satadjust (linear Newton) -> post_negfilter
//   -> cnvvar_prg2diag + bndcnd + gdcfrc!
// =============================================================================
struct BaselineWorkspace {
    double *d_Vt = nullptr, *d_kv = nullptr;
    double *d_flx_q = nullptr, *d_flx_e = nullptr, *d_flx_p = nullptr;
    double *d_flx_k = nullptr, *d_flx_u = nullptr, *d_flx_v = nullptr;
    double *d_rceff = nullptr, *d_rctop = nullptr;
    void alloc(int ijdim, int kdim) {
        size_t s3 = (size_t)ijdim * kdim;
        CUDA_CHECK(cudaMalloc(&d_Vt,    4 * s3 * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_kv,        s3 * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_flx_q, 4 * s3 * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_flx_e,     s3 * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_flx_p,     s3 * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_flx_k,     s3 * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_flx_u,     s3 * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_flx_v,     s3 * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_rceff,     s3 * sizeof(double)));
        CUDA_CHECK(cudaMalloc(&d_rctop, ijdim * sizeof(double)));
    }
    void release() {
        cudaFree(d_Vt); cudaFree(d_kv); cudaFree(d_flx_q); cudaFree(d_flx_e);
        cudaFree(d_flx_p); cudaFree(d_flx_k); cudaFree(d_flx_u); cudaFree(d_flx_v);
        cudaFree(d_rceff); cudaFree(d_rctop);
    }
};

// Unfused per-region kernels mimicking NICAM OpenACC `do l = 1, lall`
__global__ void k_bb_negfilter(HopperStatePtrs S, int l) {
    int ij = blockIdx.x * blockDim.x + threadIdx.x;
    int k  = blockIdx.y;
    if (ij >= c_P.ijdim || k < c_P.kmin0 || k > c_P.kmax0) return;
    size_t s3 = (size_t)c_P.ijdim * c_P.kdim;
    size_t off3 = (size_t)l * s3, off4 = (size_t)l * c_P.nqmax * s3;
    size_t idx = (size_t)ij + (size_t)k * c_P.ijdim;
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
    double rge = temp * rg * c_P.CVdry;
    S.rhog[off3 + idx] = rg; S.rhoge[off3 + idx] = rge;
    S.rho[off3 + idx] = rho_v; S.pre[off3 + idx] = rho_v * (qd * c_P.Rdry + qv * c_P.Rvap) * temp;
    for (int n = 0; n < 6; ++n) {
        S.rhogq[off4 + idx + n * s3] = rq[n];
        S.q[off4 + idx + n * s3] = rq[n] / rg;
    }
}

__global__ void k_bb_nsw6_core_region(HopperStatePtrs S, int l, double* Vt) {
    int ij = blockIdx.x * blockDim.x + threadIdx.x;
    int k  = blockIdx.y;
    if (ij >= c_P.ijdim || k >= c_P.kdim) return;
    size_t s3 = (size_t)c_P.ijdim * c_P.kdim;
    size_t off3 = (size_t)l * s3, off4 = (size_t)l * c_P.nqmax * s3, off2 = (size_t)l * c_P.ijdim;
    int d0=0, d1=0, d2=0, d3=0, d4=0;
    nsw6_pre_and_core_cell(c_P, ij, k,
        S.rhog + off3, S.rhoge + off3, S.rhogq + off4, S.rho + off3, S.tem + off3, S.pre + off3, S.q + off4,
        S.w + off3, S.unccn + off3, S.factor_n_aie + off2, S.gsgam2 + off3, Vt, &d0, &d1, &d2, &d3, &d4);
}

__global__ void k_bb_rceff_level(HopperStatePtrs S, int l, int k, double* rceff, double* rctop) {
    int ij = blockIdx.x * blockDim.x + threadIdx.x;
    if (ij >= c_P.ijdim) return;
    size_t s3 = (size_t)c_P.ijdim * c_P.kdim;
    size_t off3 = (size_t)l * s3, off4 = (size_t)l * c_P.nqmax * s3, off2 = (size_t)l * c_P.ijdim;
    size_t idx = (size_t)ij + (size_t)k * c_P.ijdim;
    if (k == c_P.kmin0) rctop[ij] = c_P.UNDEF;
    double dens = S.rho[off3 + idx];
    double qc = fmax(S.q[off4 + idx + c_P.i_qc * s3], 0.0);
    double qr = fmax(S.q[off4 + idx + c_P.i_qr * s3], 0.0);
    double Nc_v = c_P.Nc_def * S.factor_n_aie[off2 + ij];
    double RLMDr = (qr >= 1.0e-12) ? sqrt(sqrt(dens * qr / (c_P.Ar * c_P.N0r * c_P.GAM_1br))) : 0.0;
    double xf_qc = (dens * qc * 1000.0) / Nc_v * 1.0e-9;
    double rf_qc = pow(c_P.rceff_coef_xf * xf_qc + 1.0e-16, 1.0 / 3.0);
    double r2_qc = c_P.rceff_coef_dgam * rf_qc * rf_qc * (Nc_v * 1.0e6);
    double r2_qr = 0.25 * c_P.N0r * c_P.GAM_3 * (RLMDr * RLMDr * RLMDr);
    double r3_qc = c_P.rceff_coef_xf * dens * qc;
    double r3_qr = c_P.rceff_coef_xf * dens * qr;
    double r2_sum = r2_qc + r2_qr;
    double rc = (r2_sum >= 1.0e-10) ? ((r3_qc + r3_qr) / r2_sum) : 0.0;
    rceff[idx] = rc;
    if ((qc + qr) >= 1.0e-5) rctop[ij] = rc;
}

__global__ void k_bb_precip_region_half(HopperStatePtrs S, int l, const double* Vt, double* kv) {
    int ij = blockIdx.x * blockDim.x + threadIdx.x;
    if (ij >= c_P.ijdim) return;
    int ijdim = c_P.ijdim, kdim = c_P.kdim, kmin0 = c_P.kmin0, kmax0 = c_P.kmax0;
    size_t s3 = (size_t)ijdim * kdim, off3 = (size_t)l * s3, off4 = (size_t)l * c_P.nqmax * s3;
    const double *rhog = S.rhog + off3, *q = S.q + off4, *z = S.z + off3, *zh = S.zh + off3;
    const double *gs2 = S.gsgam2 + off3, *gs2h = S.gsgam2h + off3, *g2h = S.gam2h + off3;
    double *rhogw = S.rhogw + off3;
    int sp_nq[4] = { c_P.i_qr, c_P.i_qi, c_P.i_qs, c_P.i_qg };
    kv[ij + 0*ijdim] = 0.0; kv[ij + kmin0*ijdim] = 0.0; kv[ij + (kmax0+1)*ijdim] = 0.0;
    for (int k = kmin0 + 1; k <= kmax0; ++k) {
        size_t idx = (size_t)ij + (size_t)k * ijdim;
        double rhogh_k = (c_afact[k]/gs2[idx]*gs2h[idx])*rhog[idx] + (c_bfact[k]/gs2[idx-ijdim]*gs2h[idx])*rhog[idx-ijdim];
        kv[idx] = 0.5 * rhogw[idx] * rhogw[idx] / rhogh_k;
    }
    double dzh_kmin = z[ij + kmin0*ijdim] - zh[ij + kmin0*ijdim];
    double pfkv = 0.0, pfw = 0.0;
    for (int k = kmin0 + 1; k <= kmax0; ++k) {
        size_t kdx = (size_t)ij + (size_t)k * ijdim;
        double z_k = z[kdx], z_km1 = z[kdx - ijdim], dzh_k = z_k - z_km1;
        double dzh_km1 = (k - 1 == kmin0) ? dzh_kmin : (z_km1 - z[kdx - 2ULL*ijdim]);
        double fkv = 0.0, fw = 0.0;
        for (int s = 0; s < 4; ++s) {
            double wk = Vt[kdx + s*s3];
            double wkp1 = Vt[kdx + ijdim + s*s3], wkm1 = Vt[kdx - ijdim + s*s3];
            double zdis = vadv1d_taylor_zdis(wk, wkp1, wkm1, dzh_k, dzh_km1, dzh_k + dzh_km1, c_P.dt, k == kmin0 + 1);
            if (zdis <= 0.0) continue;
            double ztmp = z_km1 + zdis; int kc = k;
            for (int j = k; j <= kmax0; ++j) {
                if (ztmp <= z[(size_t)ij + (size_t)j*ijdim]) { if (ztmp > z[(size_t)ij + (size_t)(j-1)*ijdim]) kc = j; break; }
            }
            double rem = zdis;
            for (int j = k; j <= kc; ++j) {
                size_t jdx = (size_t)ij + (size_t)j*ijdim;
                double dj = z[jdx] - z[jdx - ijdim], L = (j < kc) ? dj : rem;
                if (j < kc) rem -= dj;
                double qh = c_afact[j]*q[jdx + sp_nq[s]*s3] + c_bfact[j]*q[jdx - ijdim + sp_nq[s]*s3];
                double rgsh = g2h[jdx] / gs2h[jdx];
                fkv -= L * (qh * kv[jdx] * rgsh);
                fw  -= L * (qh * rhogw[jdx] * rgsh);
            }
        }
        kv[kdx]    -= (fkv - pfkv) * c_rdgzh[k];
        rhogw[kdx] -= (fw  - pfw ) * c_rdgzh[k];
        pfkv = fkv; pfw = fw;
    }
}

__global__ void k_bb_precip_region_full(HopperStatePtrs S, int l, const double* Vt, const double* kv) {
    int ij = blockIdx.x * blockDim.x + threadIdx.x;
    if (ij >= c_P.ijdim) return;
    int ijdim = c_P.ijdim, kdim = c_P.kdim, kmin0 = c_P.kmin0, kmax0 = c_P.kmax0;
    size_t s3 = (size_t)ijdim * kdim, off3 = (size_t)l * s3, off4 = (size_t)l * c_P.nqmax * s3;
    size_t off2 = (size_t)l * ijdim, offp = (size_t)l * 2ULL * ijdim;
    double *rhog = S.rhog+off3, *rhogvx = S.rhogvx+off3, *rhogvy = S.rhogvy+off3, *rhogvz = S.rhogvz+off3;
    double *rhoge = S.rhoge+off3, *rhogq = S.rhogq+off4;
    const double *rhogw = S.rhogw+off3, *tem = S.tem+off3, *vx = S.vx+off3, *vy = S.vy+off3, *vz = S.vz+off3;
    const double *z = S.z+off3, *zh = S.zh+off3, *gs2 = S.gsgam2+off3, *gs2h = S.gsgam2h+off3, *g2 = S.gam2+off3;
    double ix_c = S.ix[off2+ij], iy_c = S.iy[off2+ij], iz_c = S.iz[off2+ij];
    double jx_c = S.jx[off2+ij], jy_c = S.jy[off2+ij], jz_c = S.jz[off2+ij];
    int sp_nq[4] = { c_P.i_qr, c_P.i_qi, c_P.i_qs, c_P.i_qg };
    auto ew = [&](int k, int s) -> double {
        if (k == kmin0) return Vt[ij + kmin0*ijdim + s*s3];
        if (k > kmax0)  return Vt[ij + kmax0*ijdim + s*s3];
        size_t idx = (size_t)ij + (size_t)k*ijdim;
        return c_afact[k]*Vt[idx + s*s3] + c_bfact[k]*Vt[idx - ijdim + s*s3];
    };
    auto efi = [&](int ki, double fq[4], double& fe, double& fp, double& fk, double& fu, double& fv) {
        fq[0]=fq[1]=fq[2]=fq[3]=0.0; fe=fp=fk=fu=fv=0.0;
        if (ki > kmax0) return;
        size_t kdx = (size_t)ij + (size_t)ki*ijdim;
        double zh_ki = zh[kdx], dz_k = zh[kdx+ijdim] - zh_ki;
        double dz_km1 = (ki > kmin0) ? (zh_ki - zh[kdx-ijdim]) : dz_k;
        for (int s = 0; s < 4; ++s) {
            double zdis = vadv1d_taylor_zdis(ew(ki,s), ew(ki+1,s), (ki>kmin0)?ew(ki-1,s):0.0, dz_k, dz_km1, dz_k+dz_km1, c_P.dt, ki==kmin0);
            if (zdis <= 0.0) continue;
            double ztmp = zh_ki + zdis; int kc = ki;
            for (int j = ki; j <= kmax0; ++j) {
                size_t jdx = (size_t)ij + (size_t)j*ijdim;
                if (ztmp <= zh[jdx+ijdim]) { if (ztmp > zh[jdx]) kc = j; break; }
            }
            double rem = zdis;
            for (int j = ki; j <= kc; ++j) {
                size_t jdx = (size_t)ij + (size_t)j*ijdim;
                double dj = zh[jdx+ijdim] - zh[jdx], L = (j < kc) ? dj : rem;
                if (j < kc) rem -= dj;
                double w_sj = L * rhogq[jdx + sp_nq[s]*s3] * (g2[jdx]/gs2[jdx]);
                double rg_j = rhog[jdx];
                double kh_r = 0.5*(rhogvx[jdx]*rhogvx[jdx]+rhogvy[jdx]*rhogvy[jdx]+rhogvz[jdx]*rhogvz[jdx])/(rg_j*rg_j);
                fq[s] -= w_sj; fe -= w_sj * (c_P.CVW[sp_nq[s]] * tem[jdx]); fp -= w_sj * (c_P.GRAV * z[jdx]);
                fk -= w_sj * kh_r;
                fu -= w_sj * (vx[jdx]*ix_c + vy[jdx]*iy_c + vz[jdx]*iz_c);
                fv -= w_sj * (vx[jdx]*jx_c + vy[jdx]*jy_c + vz[jdx]*jz_c);
            }
        }
    };
    double fq_lo[4], fe_lo, fp_lo, fk_lo, fu_lo, fv_lo;
    efi(kmin0, fq_lo, fe_lo, fp_lo, fk_lo, fu_lo, fv_lo);
    double inv_dt = 1.0 / c_P.dt;
    S.precip[offp+ij] = -fq_lo[0]*inv_dt;
    S.precip[offp+ij+ijdim] = -(fq_lo[1]+fq_lo[2]+fq_lo[3])*inv_dt;
    S.precip_lh[off2+ij] = (fq_lo[1]+fq_lo[2]+fq_lo[3])*c_P.LHF*inv_dt;
    S.precip_rhoein[off2+ij] = -fe_lo*inv_dt;
    S.precip_rhophi[off2+ij] = -fp_lo*inv_dt;
    S.precip_rhokin[off2+ij] = -fk_lo*inv_dt;
    S.precip_energy[off2+ij] = S.precip_lh[off2+ij] - (fe_lo+fp_lo+fk_lo)*inv_dt;
    double d_rg=0, d_rge=0, d_rkht=0, d_rkhn=0, d_kvp=0;
    for (int k = kmin0; k <= kmax0; ++k) {
        double fq_hi[4], fe_hi, fp_hi, fk_hi, fu_hi, fv_hi;
        efi(k+1, fq_hi, fe_hi, fp_hi, fk_hi, fu_hi, fv_hi);
        size_t idx = (size_t)ij + (size_t)k*ijdim;
        double r = c_rdgz[k];
        double dq0=-(fq_hi[0]-fq_lo[0])*r, dq1=-(fq_hi[1]-fq_lo[1])*r, dq2=-(fq_hi[2]-fq_lo[2])*r, dq3=-(fq_hi[3]-fq_lo[3])*r;
        double sdq = dq0+dq1+dq2+dq3;
        double rg_o = rhog[idx], rvx_o = rhogvx[idx], rvy_o = rhogvy[idx], rvz_o = rhogvz[idx];
        double rkh_o = 0.5*(rvx_o*rvx_o+rvy_o*rvy_o+rvz_o*rvz_o)/rg_o;
        double rg_n = rg_o + sdq;
        double du = -(fu_hi-fu_lo)*r, dv = -(fv_hi-fv_lo)*r;
        double rvx_n = rvx_o + du*ix_c + dv*jx_c, rvy_n = rvy_o + du*iy_c + dv*jy_c, rvz_n = rvz_o + du*iz_c + dv*jz_c;
        double rge_n = rhoge[idx] - (fe_hi-fe_lo)*r - (fp_hi-fp_lo)*r - sdq*c_P.GRAV*z[idx];
        double rkh_t = rkh_o - (fk_hi-fk_lo)*r;
        double rkh_n = 0.5*(rvx_n*rvx_n+rvy_n*rvy_n+rvz_n*rvz_n)/rg_n;
        rhog[idx]=rg_n; rhogvx[idx]=rvx_n; rhogvy[idx]=rvy_n; rhogvz[idx]=rvz_n;
        rhogq[idx+sp_nq[0]*s3]+=dq0; rhogq[idx+sp_nq[1]*s3]+=dq1; rhogq[idx+sp_nq[2]*s3]+=dq2; rhogq[idx+sp_nq[3]*s3]+=dq3;
        double kvr = 0.0;
        if (k > kmin0) {
            size_t im1 = idx - ijdim;
            double rhogh_n = (c_afact[k]/gs2[idx]*gs2h[idx])*rg_n + (c_bfact[k]/gs2[im1]*gs2h[idx])*d_rg;
            kvr = 0.5*rhogw[idx]*rhogw[idx]/rhogh_n;
            double t2 = d_rkht + (c_cfact[k-1]*kv[idx] + c_dfact[k-1]*kv[im1]);
            double t1 = d_rkhn + ((c_cfact[k-1]*gs2[im1]/gs2h[idx])*kvr + (c_dfact[k-1]*gs2[im1]/gs2h[im1])*d_kvp);
            double rge_fin = d_rge + (t2 - t1);
            rhoge[im1] = rge_fin;
            S.rho[off3 + im1] = d_rg / gs2[im1];
            S.tem[off3 + im1] = rge_fin / (d_rg * c_P.CVdry);
            for (int n = 0; n < 6; ++n) S.q[off4 + im1 + n*s3] = rhogq[im1 + n*s3] / d_rg;
        }
        d_rg=rg_n; d_rge=rge_n; d_rkht=rkh_t; d_rkhn=rkh_n; d_kvp=kvr;
        for (int s=0;s<4;++s) fq_lo[s]=fq_hi[s]; fe_lo=fe_hi; fp_lo=fp_hi; fk_lo=fk_hi; fu_lo=fu_hi; fv_lo=fv_hi;
    }
    size_t ixk = (size_t)ij + (size_t)kmax0*ijdim;
    double rge_top = d_rge + (d_rkht + c_dfact[kmax0]*kv[ixk]) - (d_rkhn + (c_dfact[kmax0]*gs2[ixk]/gs2h[ixk])*d_kvp);
    rhoge[ixk] = rge_top;
    S.rho[off3 + ixk] = d_rg / gs2[ixk];
    S.tem[off3 + ixk] = rge_top / (d_rg * c_P.CVdry);
    for (int n = 0; n < 6; ++n) S.q[off4 + ixk + n*s3] = rhogq[ixk + n*s3] / d_rg;
}

__global__ void k_bb_satadj_region(HopperStatePtrs S, int l) {
    int ij = blockIdx.x * blockDim.x + threadIdx.x;
    int k  = blockIdx.y;
    if (ij >= c_P.ijdim || k >= c_P.kdim) return;
    size_t s3 = (size_t)c_P.ijdim * c_P.kdim, off3 = (size_t)l*s3, off4 = (size_t)l*c_P.nqmax*s3;
    int d0=0,d1=0,d2=0,d3=0,d4=0;
    satadj_and_post_negfilter_cell(c_P, ij, k, S.rhog+off3, S.rhoge+off3, S.rhogq+off4,
        S.rho+off3, S.tem+off3, S.pre+off3, S.q+off4, S.gsgam2+off3, &d0,&d1,&d2,&d3,&d4);
}

__global__ void k_bb_cnvvar_region(HopperStatePtrs S, int l) {
    int ij = blockIdx.x * blockDim.x + threadIdx.x;
    if (ij >= c_P.ijdim) return;
    size_t s3 = (size_t)c_P.ijdim * c_P.kdim, off3 = (size_t)l*s3, off4 = (size_t)l*c_P.nqmax*s3;
    cnvvar_bndcnd_cldfrac_column(c_P, ij, S.rhog+off3, S.rhogvx+off3, S.rhogvy+off3, S.rhogvz+off3,
        S.rhogw+off3, S.rhoge+off3, S.rhogq+off4, S.rho+off3, S.pre+off3, S.tem+off3,
        S.vx+off3, S.vy+off3, S.vz+off3, S.w+off3, S.q+off4, S.gdcfrc+off3, S.gdclw+off3,
        S.phi+off3, S.gsgam2+off3, S.gsgam2h+off3, c_afact, c_bfact);
}

static void launch_bb_h200_step(const NSW6Params& P, int lall, const HopperStatePtrs& S, BaselineWorkspace& BW, cudaStream_t strm) {
    dim3 g256k((P.ijdim + 255) / 256, P.kdim);
    dim3 g128k((P.ijdim + 127) / 128, P.kdim);
    dim3 gcol ((P.ijdim + 127) / 128);
    // Exact NICAM mod_mp_driver.f90 region-by-region OpenACC structure (l = 0..lall-1)!
    for (int l = 0; l < lall; ++l) {
        k_bb_negfilter<<<g256k, 256, 0, strm>>>(S, l);
        k_bb_nsw6_core_region<<<g128k, 128, 0, strm>>>(S, l, BW.d_Vt);
        if (P.compute_rceff) {
            // Mirrors NICAM mod_mp_nsw6.f90:2447-2520 sequential k-loop (78 launches/region = 780/step)!
            for (int k = P.kmin0; k <= P.kmax0; ++k) {
                k_bb_rceff_level<<<gcol, 128, 0, strm>>>(S, l, k, BW.d_rceff, BW.d_rctop);
            }
        }
        k_bb_precip_region_half<<<gcol, 128, 0, strm>>>(S, l, BW.d_Vt, BW.d_kv);
        k_bb_precip_region_full<<<gcol, 128, 0, strm>>>(S, l, BW.d_Vt, BW.d_kv);
        k_bb_satadj_region<<<g256k, 256, 0, strm>>>(S, l);
        k_bb_cnvvar_region<<<gcol, 128, 0, strm>>>(S, l);
    }
}

// =============================================================================
// Problem Parameters & Host/Device State Setup (exact NICAM production config)
// =============================================================================
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

// Exact GPU budget reduction kernel for 5-step & 100-step full-scale conservation!
__constant__ double c_dgz[128];
__global__ void k_reduce_budgets(HopperStatePtrs S, int lall, double* d_out) {
    // d_out[0] = wmass, d_out[1] = tenergy, d_out[2] = step_precip, d_out[3] = step_pene
    const int ij = blockIdx.x * blockDim.x + threadIdx.x;
    const int l  = blockIdx.y;
    const int ijdim = c_P.ijdim;
    double wm = 0.0, te = 0.0, sp = 0.0, spe = 0.0;
    if (ij < ijdim) {
        const size_t s3 = (size_t)ijdim * c_P.kdim;
        const size_t off3 = (size_t)l * s3, off4 = (size_t)l * c_P.nqmax * s3;
        const size_t off2 = (size_t)l * ijdim, offp = (size_t)l * 2ULL * ijdim;
        sp  = (S.precip[offp + ij] + S.precip[offp + ij + ijdim]) * c_P.dt;
        spe = S.precip_energy[off2 + ij] * c_P.dt;
        for (int k = c_P.kmin0; k <= c_P.kmax0; ++k) {
            size_t idx = (size_t)ij + (size_t)k * ijdim;
            double dgz_k = c_dgz[k];
            double qv = S.rhogq[off4 + idx + 0 * s3];
            double qc = S.rhogq[off4 + idx + 1 * s3];
            double qr = S.rhogq[off4 + idx + 2 * s3];
            double qi = S.rhogq[off4 + idx + 3 * s3];
            double qs = S.rhogq[off4 + idx + 4 * s3];
            double qg = S.rhogq[off4 + idx + 5 * s3];
            wm += (qv + qc + qr + qi + qs + qg) * dgz_k;
            double rg = S.rhog[off3 + idx];
            double rvx = S.rhogvx[off3 + idx], rvy = S.rhogvy[off3 + idx], rvz = S.rhogvz[off3 + idx];
            double kh = 0.5 * (rvx*rvx + rvy*rvy + rvz*rvz) / rg;
            te += (S.rhoge[off3 + idx] + rg * c_P.GRAV * S.z[off3 + idx] + kh
                   + c_P.LHV * qv - c_P.LHF * (qi + qs + qg)) * dgz_k;
        }
    }
    // Warp + block reduction
    for (int off = 16; off > 0; off >>= 1) {
        wm  += __shfl_down_sync(0xffffffffu, wm,  off);
        te  += __shfl_down_sync(0xffffffffu, te,  off);
        sp  += __shfl_down_sync(0xffffffffu, sp,  off);
        spe += __shfl_down_sync(0xffffffffu, spe, off);
    }
    if ((threadIdx.x & 31) == 0) {
        atomicAdd(&d_out[0], wm);
        atomicAdd(&d_out[1], te);
        atomicAdd(&d_out[2], sp);
        atomicAdd(&d_out[3], spe);
    }
}

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
    const bool nsys_mode = (argc > 1 && std::strcmp(argv[1], "--nsys") == 0);
    const int BLOCKS  = (!nsys_mode && argc > 1) ? std::atoi(argv[1]) : 8;
    const int REPEATS = (!nsys_mode && argc > 2) ? std::atoi(argv[2]) : 50;

    NSW6Params P_ref = build_nicam_production_params(IJ, KD, 1, 1, 0);

    if (nsys_mode) {
        HostState H; H.init(IJ, KD, LALL, P_ref);
        DeviceState D; D.alloc_and_upload(H);
        NSW6Params Pm2 = build_nicam_production_params(IJ, KD, 2, 1, 0);
        sm90_v2_upload_constants(Pm2, H.rdgz.data(), H.rdgzh.data(), H.afact.data(), H.bfact.data(), H.cfact.data(), H.dfact.data());
        HopperWorkspaceV3 W; cudaStream_t strm; CUDA_CHECK(cudaStreamCreate(&strm));
        sm90_v3_launch_step(Pm2, LALL, D.cur, W, strm); // warmup
        D.restore(strm);
        CUDA_CHECK(cudaStreamSynchronize(strm));
        sm90_v3_launch_step(Pm2, LALL, D.cur, W, strm); // clean profiled step
        CUDA_CHECK(cudaStreamSynchronize(strm));
        W.release(); CUDA_CHECK(cudaStreamDestroy(strm));
        return 0;
    }

    cudaDeviceProp prop; CUDA_CHECK(cudaGetDeviceProperties(&prop, 0));
    printf("================================================================================\n");
    printf(" NICAM NSW6 HOPPER SM90a v3 STREAM-COMPACTED MASTER BENCHMARK\n");
    printf(" GPU: %s (SM %d.%d, %d SMs, %.1f GB HBM3e, %.1f MB L2)\n",
           prop.name, prop.major, prop.minor, prop.multiProcessorCount, prop.totalGlobalMem/1e9, prop.l2CacheSize/1e6);
    printf(" Workload: full NICAM scale (%d cols/reg x %d lvls x %d regs = 13.31M cells)\n", IJ, KD, LALL);
    printf("================================================================================\n\n");

    // -------------------------------------------------------------------------
    // Phase 1: Exact GPU-vs-CPU Differential Validation (BB-H200, Mode 1, Mode 2)
    // -------------------------------------------------------------------------
    {
        printf("--- [Phase 1] Exact GPU-vs-CPU Differential Validation (region0: 16,641 x 80) ---\n");
        HostState H_ref; H_ref.init(IJ, KD, 1, P_ref);
        DeviceState D_val; D_val.alloc_and_upload(H_ref);
        run_cpu_ref_region0(P_ref, H_ref);
        HopperWorkspaceV3 W_val; BaselineWorkspace BW_val; BW_val.alloc(IJ, KD);
        cudaStream_t strm; CUDA_CHECK(cudaStreamCreate(&strm));

        for (int mode : {0, 1, 2}) {
            NSW6Params Pm = build_nicam_production_params(IJ, KD, (mode == 0) ? 1 : mode, 1, 0);
            sm90_v2_upload_constants(Pm, H_ref.rdgz.data(), H_ref.rdgzh.data(),
                                     H_ref.afact.data(), H_ref.bfact.data(), H_ref.cfact.data(), H_ref.dfact.data());
            D_val.restore(strm);
            if (mode == 0) launch_bb_h200_step(Pm, 1, D_val.cur, BW_val, strm);
            else           sm90_v3_launch_step(Pm, 1, D_val.cur, W_val, strm);
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
            const char* lbl = (mode == 0) ? "BB-H200-Unfused " : (mode == 1 ? "ENG-SM90-v3     " : "ENG+MATH-SM90-v3");
            printf("  Mode %d (%s) GPU vs CPU Ref: N_FAIL = %lld | max_rel = %.3e | max_abs = %.3e -> %s\n",
                   mode, lbl, nf, wr, wa, nf==0 ? "PASS" : "FAIL");
            if (nf != 0) return 2;
        }
        BW_val.release(); W_val.release(); CUDA_CHECK(cudaStreamDestroy(strm));
        printf("\n");
    }

    // -------------------------------------------------------------------------
    // Phase 2: Full 10-Region Paired Balanced Benchmark (8 blocks x 50 repeats)
    // -------------------------------------------------------------------------
    HostState H; H.init(IJ, KD, LALL, P_ref);
    DeviceState D; D.alloc_and_upload(H);
    HopperWorkspaceV3 W; BaselineWorkspace BW; BW.alloc(IJ, KD);
    cudaStream_t strm; CUDA_CHECK(cudaStreamCreate(&strm));
    cudaEvent_t ev[12]; for (int i = 0; i < 12; ++i) CUDA_CHECK(cudaEventCreate(&ev[i]));

    const int NCFG = 6;
    const char* names[NCFG] = {
        "Mode 0: BB-H200 (unfused, 16-trans Roh)",
        "Mode 1: ENG-SM90-v3 (0-spill, orig math)",
        "Mode 2: ENG+MATH-SM90-v3 (Stream, w/ rceff)",
        "Mode 2: ENG+MATH-SM90-v3 (Graph,  w/ rceff)",
        "Mode 2: ENG+MATH-SM90-v3 (Graph, no rceff)",
        "Mode 3: ENG+MATH+MONO-v3 (Graph,  w/ rceff)"
    };
    const int modes[NCFG]  = { 0, 1, 2, 2, 2, 3 };
    const int rceffs[NCFG] = { 1, 1, 1, 1, 0, 1 };
    const int graphs[NCFG] = { 0, 0, 0, 1, 1, 1 };

    // Warmup all modes & instantiate CUDA Graphs
    for (int m : {0, 1, 2, 3}) {
        NSW6Params Pm = build_nicam_production_params(IJ, KD, m == 0 ? 1 : m, 1, 0);
        sm90_v2_upload_constants(Pm, H.rdgz.data(), H.rdgzh.data(), H.afact.data(), H.bfact.data(), H.cfact.data(), H.dfact.data());
        D.restore(strm);
        if (m == 0) launch_bb_h200_step(Pm, LALL, D.cur, BW, strm);
        else        sm90_v3_launch_step(Pm, LALL, D.cur, W, strm);
    }
    CUDA_CHECK(cudaStreamSynchronize(strm));

    cudaGraph_t g_m2_rc, g_m2_norc, g_m3_rc;
    cudaGraphExec_t ge_m2_rc, ge_m2_norc, ge_m3_rc;
    auto cap = [&](int mode, int rc, cudaGraph_t& g, cudaGraphExec_t& ge) {
        NSW6Params Pm = build_nicam_production_params(IJ, KD, mode, rc, 0);
        sm90_v2_upload_constants(Pm, H.rdgz.data(), H.rdgzh.data(), H.afact.data(), H.bfact.data(), H.cfact.data(), H.dfact.data());
        CUDA_CHECK(cudaStreamBeginCapture(strm, cudaStreamCaptureModeGlobal));
        sm90_v3_launch_step(Pm, LALL, D.cur, W, strm);
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
            NSW6Params Pm = build_nicam_production_params(IJ, KD, modes[cfg] == 0 ? 1 : modes[cfg], rceffs[cfg], 0);
            sm90_v2_upload_constants(Pm, H.rdgz.data(), H.rdgzh.data(), H.afact.data(), H.bfact.data(), H.cfact.data(), H.dfact.data());
            int nrep = (modes[cfg] == 0) ? std::min(REPEATS, 10) : REPEATS; // BB is ~30ms so 8x10=80 samples
            for (int r = 0; r < nrep; ++r) {
                D.restore(strm);
                CUDA_CHECK(cudaEventRecord(ev[0], strm));
                if (modes[cfg] == 0) {
                    launch_bb_h200_step(Pm, LALL, D.cur, BW, strm);
                } else if (graphs[cfg]) {
                    cudaGraphExec_t ge = (modes[cfg] == 3) ? ge_m3_rc : (rceffs[cfg] ? ge_m2_rc : ge_m2_norc);
                    CUDA_CHECK(cudaGraphLaunch(ge, strm));
                } else {
                    sm90_v3_launch_step(Pm, LALL, D.cur, W, strm);
                }
                CUDA_CHECK(cudaEventRecord(ev[1], strm));
                CUDA_CHECK(cudaEventSynchronize(ev[1]));
                float ms = 0.0f; CUDA_CHECK(cudaEventElapsedTime(&ms, ev[0], ev[1]));
                samples[cfg].push_back((double)ms);
            }
        }
    }

    double bb_h200_mean = std::accumulate(samples[0].begin(), samples[0].end(), 0.0) / samples[0].size();
    printf("--- [Phase 2] Paired Balanced H200 v3 Benchmark (%d blocks x %d repeats) ---\n", BLOCKS, REPEATS);
    printf("  %-45s | %7s | %7s | %7s | %6s | %8s | %8s | %8s\n",
           "Configuration", "Mean", "Med", "Min", "Std", "vs H200B", "vs GH_BB", "vs GH_TT");
    printf("  -----------------------------------------------------------------------------------------------------------------------\n");
    for (int cfg = 0; cfg < NCFG; ++cfg) {
        auto& v = samples[cfg];
        double mean = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
        double sq = 0.0; for (double x : v) sq += (x - mean)*(x - mean);
        double sd = std::sqrt(sq / std::max((size_t)1, v.size() - 1));
        std::sort(v.begin(), v.end());
        printf("  %-45s | %7.3f | %7.3f | %7.3f | %6.3f | %7.2fx | %7.2fx | %7.2fx\n",
               names[cfg], mean, v[v.size()/2], v.front(), sd,
               bb_h200_mean/mean, 78.19/mean, 68.57/mean);
    }

    // -------------------------------------------------------------------------
    // Phase 3: Per-Kernel Component Breakdown & Achieved HBM3e Bandwidth (v3)
    // -------------------------------------------------------------------------
    printf("\n--- [Phase 3] Per-Kernel Component Breakdown (v3, 10 regions = 13.31M cells, 20-step avg) ---\n");
    {
        const int PERSIST_BLKS = 528;
        const size_t n3 = (size_t)IJ * KD * LALL;
        dim3 g_scr((IJ + 255) / 256, KD - 2, LALL);
        dim3 g_128k((IJ + 127) / 128, KD, LALL), g_256k((IJ + 255) / 256, KD, LALL);
        dim3 g_col128((IJ + 127) / 128, LALL), g_col256((IJ + 255) / 256, LALL);

        for (int m : {2, 3}) {
            NSW6Params Pm = build_nicam_production_params(IJ, KD, m, 1, 0);
            sm90_v2_upload_constants(Pm, H.rdgz.data(), H.rdgzh.data(), H.afact.data(), H.bfact.data(), H.cfact.data(), H.dfact.data());
            double t[10] = {}; const int NR = 20;
            uint32_t h_wc = 0, h_wcol = 0;
            for (int r = 0; r < NR; ++r) {
                D.restore(strm);
                CUDA_CHECK(cudaEventRecord(ev[0], strm));
                CUDA_CHECK(cudaMemsetAsync(W.d_Vt, 0, 4 * n3 * sizeof(double), strm));
                CUDA_CHECK(cudaMemsetAsync(W.d_meta_block, 0, W.meta_bytes, strm));
                CUDA_CHECK(cudaEventRecord(ev[1], strm));
                k_sm90_v3_nsw6_screen<false><<<g_scr, 256, 0, strm>>>(D.cur, W.d_wet_cells, W.d_wet_count(), W.d_cnt);
                CUDA_CHECK(cudaEventRecord(ev[2], strm));
                k_sm90_v3_nsw6_wet_core<false><<<PERSIST_BLKS, 128, 0, strm>>>(
                    D.cur, W.d_wet_cells, W.d_wet_count(), W.d_Vt, W.d_col_spmask(), W.d_col_ktop(), W.d_cnt);
                CUDA_CHECK(cudaEventRecord(ev[3], strm));
                k_sm90_v3_rceff<<<g_col128, 128, 0, strm>>>(D.cur, W.d_rceff, W.d_rctop);
                CUDA_CHECK(cudaEventRecord(ev[4], strm));
                k_sm90_v3_compact_cols<false><<<g_col256, 256, 0, strm>>>(
                    D.cur, W.d_col_spmask(), W.d_wet_cols, W.d_wet_col_count(), W.d_cnt);
                CUDA_CHECK(cudaEventRecord(ev[5], strm));
                if (m == 2) k_sm90_v3_precip_half<2,false><<<PERSIST_BLKS, 128, 0, strm>>>(
                    D.cur, W.d_wet_cols, W.d_wet_col_count(), W.d_col_spmask(), W.d_col_ktop(), W.d_Vt, W.d_kv, W.d_zdis_m, W.d_cnt);
                else        k_sm90_v3_precip_half<3,false><<<PERSIST_BLKS, 128, 0, strm>>>(
                    D.cur, W.d_wet_cols, W.d_wet_col_count(), W.d_col_spmask(), W.d_col_ktop(), W.d_Vt, W.d_kv, W.d_zdis_m, W.d_cnt);
                CUDA_CHECK(cudaEventRecord(ev[6], strm));
                if (m == 2) k_sm90_v3_precip_full<2,false><<<PERSIST_BLKS, 128, 0, strm>>>(
                    D.cur, W.d_wet_cols, W.d_wet_col_count(), W.d_col_spmask(), W.d_col_ktop(), W.d_Vt, W.d_kv, W.d_zdis_m, W.d_cnt);
                else        k_sm90_v3_precip_full<3,false><<<PERSIST_BLKS, 128, 0, strm>>>(
                    D.cur, W.d_wet_cols, W.d_wet_col_count(), W.d_col_spmask(), W.d_col_ktop(), W.d_Vt, W.d_kv, W.d_zdis_m, W.d_cnt);
                CUDA_CHECK(cudaEventRecord(ev[7], strm));
                if (m == 2) k_sm90_v3_satadj_nf_diag<2,false,true><<<g_128k, 128, 0, strm>>>(D.cur, W.d_col_ktop(), W.d_cnt);
                else        k_sm90_v3_satadj_nf_diag<3,false,true><<<g_128k, 128, 0, strm>>>(D.cur, W.d_col_ktop(), W.d_cnt);
                CUDA_CHECK(cudaEventRecord(ev[8], strm));
                k_sm90_v2_diag_w<<<g_256k, 256, 0, strm>>>(D.cur);
                CUDA_CHECK(cudaEventRecord(ev[9], strm));
                k_sm90_v2_bndcnd<<<g_col256, 256, 0, strm>>>(D.cur);
                CUDA_CHECK(cudaEventRecord(ev[10], strm));
                CUDA_CHECK(cudaEventSynchronize(ev[10]));
                for (int i = 0; i < 10; ++i) { float ms; CUDA_CHECK(cudaEventElapsedTime(&ms, ev[i], ev[i+1])); t[i] += ms; }
                if (r == 0) {
                    CUDA_CHECK(cudaMemcpy(&h_wc,   W.d_wet_count(),     sizeof(uint32_t), cudaMemcpyDeviceToHost));
                    CUDA_CHECK(cudaMemcpy(&h_wcol, W.d_wet_col_count(), sizeof(uint32_t), cudaMemcpyDeviceToHost));
                }
            }
            for (int i = 0; i < 10; ++i) t[i] /= NR;
            double sum = 0.0; for (int i = 0; i < 10; ++i) sum += t[i];
            printf("  Mode %d: Memset=%.3f | Screen=%.3f | WetCore=%.3f | rceff=%.3f | CmpCol=%.3f | Precip(H+F)=%.3f (%.3f+%.3f) | Sat+NF+Diag=%.3f | w=%.3f | BC=%.3f | SUM=%.3f ms\n",
                   m, t[0], t[1], t[2], t[3], t[4], t[5]+t[6], t[5], t[6], t[7], t[8], t[9], sum);
            printf("          Compaction: wet_cells = %u / %zu (%.1f%%), wet_cols = %u / %d (%.1f%%)\n",
                   h_wc, (size_t)IJ*(KD-2)*LALL, 100.0*h_wc/((double)IJ*(KD-2)*LALL),
                   h_wcol, IJ*LALL, 100.0*h_wcol/(double)(IJ*LALL));
        }
    }

    // -------------------------------------------------------------------------
    // Phase 4: Full-Scale 10-Region Continuous Multi-Step Conservation (5 & 100 steps)
    // -------------------------------------------------------------------------
    printf("\n--- [Phase 4] Full-Scale 10-Region Continuous Evolution Conservation (166,410 cols x 80 lvls) ---\n");
    {
        CUDA_CHECK(cudaMemcpyToSymbol(c_dgz, H.dgz.data(), KD * sizeof(double)));
        double* d_bud = nullptr; CUDA_CHECK(cudaMalloc(&d_bud, 4 * sizeof(double)));
        dim3 g_red((IJ + 255) / 256, LALL);
        auto get_bud = [&](double out[4]) {
            CUDA_CHECK(cudaMemsetAsync(d_bud, 0, 4 * sizeof(double), strm));
            k_reduce_budgets<<<g_red, 256, 0, strm>>>(D.cur, LALL, d_bud);
            CUDA_CHECK(cudaMemcpyAsync(out, d_bud, 4 * sizeof(double), cudaMemcpyDeviceToHost, strm));
            CUDA_CHECK(cudaStreamSynchronize(strm));
        };

        for (int nsteps : {5, 100}) {
            printf("  [N_STEPS = %3d | Physical T = %6.1f s]\n", nsteps, nsteps * 30.0);
            for (int m : {1, 2, 3}) {
                NSW6Params Pm = build_nicam_production_params(IJ, KD, m, 0, 1);
                sm90_v2_upload_constants(Pm, H.rdgz.data(), H.rdgzh.data(), H.afact.data(), H.bfact.data(), H.cfact.data(), H.dfact.data());
                CUDA_CHECK(cudaMemsetAsync(W.d_cnt, 0, sizeof(NSW6Counters), strm));
                D.restore(strm);
                double b0[4]; get_bud(b0);
                double cum_p = 0.0, cum_pe = 0.0;
                double bcur[4] = {};
                CUDA_CHECK(cudaEventRecord(ev[0], strm));
                for (int st = 0; st < nsteps; ++st) {
                    sm90_v3_launch_step(Pm, LALL, D.cur, W, strm);
                    CUDA_CHECK(cudaMemsetAsync(d_bud, 0, 4 * sizeof(double), strm));
                    k_reduce_budgets<<<g_red, 256, 0, strm>>>(D.cur, LALL, d_bud);
                    CUDA_CHECK(cudaMemcpyAsync(bcur, d_bud, 4 * sizeof(double), cudaMemcpyDeviceToHost, strm));
                    CUDA_CHECK(cudaStreamSynchronize(strm));
                    cum_p  += bcur[2];
                    cum_pe += bcur[3];
                }
                CUDA_CHECK(cudaEventRecord(ev[1], strm));
                CUDA_CHECK(cudaEventSynchronize(ev[1]));
                float tot_ms = 0.0f; CUDA_CHECK(cudaEventElapsedTime(&tot_ms, ev[0], ev[1]));
                NSW6Counters hc = {};
                CUDA_CHECK(cudaMemcpy(&hc, W.d_cnt, sizeof(NSW6Counters), cudaMemcpyDeviceToHost));
                double rel_m = std::fabs((bcur[0] + cum_p) - b0[0]) / std::fabs(b0[0]);
                double rel_e = std::fabs((bcur[1] + cum_pe) - b0[1]) / std::fabs(b0[1]);
                const char* mlbl = (m == 1) ? "Mode 1 (BB/ENG) " : (m == 2 ? "Mode 2 (MATH)   " : "Mode 3 (MONO)   ");
                printf("    %s | rel_mass_drift=%.3e | rel_E_drift=%.3e | cum_precip=%.5e | cross(F/H)=%llu/%llu | neg_post=%llu\n",
                       mlbl, rel_m, rel_e, cum_p, hc.precip_crossings_full, hc.precip_crossings_half, hc.neg_q_before_postfilter);
            }
        }
        CUDA_CHECK(cudaFree(d_bud));
    }

    CUDA_CHECK(cudaGraphExecDestroy(ge_m2_rc)); CUDA_CHECK(cudaGraphDestroy(g_m2_rc));
    CUDA_CHECK(cudaGraphExecDestroy(ge_m2_norc)); CUDA_CHECK(cudaGraphDestroy(g_m2_norc));
    CUDA_CHECK(cudaGraphExecDestroy(ge_m3_rc)); CUDA_CHECK(cudaGraphDestroy(g_m3_rc));
    BW.release(); W.release(); CUDA_CHECK(cudaStreamDestroy(strm));
    return 0;
}
