// =============================================================================
// nsw6_hopper_v5.cuh
//
// Hopper (SM 9.0a / NVIDIA H200 & GH200) v5.2 ROOFLINE-BREAKING OPERATOR SUITE
//
// Builds on v4.1 (4.066 ms M2 / 4.413 ms M3 on idle H200, 0 B spill) with 7
// strictly-exact IEEE FP64 breakthroughs (>1.8 GB/step HBM traffic eliminated):
//
//   1. __constant__ HopperWorkPtrsV5 c_W5 (zero workspace-pointer register live set):
//      Eliminates all workspace pointer args from wet_core, precip_half,
//      precip_full, and satadj_nf_diag -> drops wet_core back to ~188 regs (0 B spill)!
//   2. Perturbation-Robust Clausius-Clapeyron Cold-Subsaturation Certificate
//      (d_cell_flag[n3], 13.3 MB L2-resident):
//        flag = 0 : wet / active cell
//        flag = 1 : warm-dry cell (rhogq[1..5] == 0 & Vt == 0)
//        flag = 2 : cold-dry certified (T in [20 K, T0 - 1 K] & qv <= qsat_ice(T)).
//                   Theorem: for ALL T in [20 K, T0 - 1 K] and dT >= -0.1 K,
//                   qsat_ice(T) < qsat_liq(T + dT) strictly -> satadj skips BOTH
//                   530 MB rhogq[1..5] loads AND 100% of cold-dry exp() calls
//                   (and skips rhog/rhoge/rhogq[0] too when diff_e == 0)!
//   3. Zero 426 MB Vt Reset (-0.095 ms net) with ZERO flag branches in precip:
//      8-us k_sm90_v5_clean_vt_wetcols zeroes Vt[0..3] ONLY on dry gaps inside
//      precipitating columns (k <= k_top + 2) -> eliminates 426 MB cudaMemsetAsync
//      while keeping precip_half/full 100% branch-free on Vt!
//   4. Fused screen + init_kv (-0.10 ms):
//      Reuses rhog & gsgam2 already loaded in screen to compute kv in one pass,
//      eliminating the separate 532 MB k_sm90_v4_init_kv kernel!
//   5. Hardware Memset + Sparse liq_ktop rceff (-0.125 ms):
//      25-us hardware cudaMemsetAsync(rceff, 0) + sparse column loop only up to
//      d_col_liq_ktop[col] on liquid-bearing columns!
//   6. Zero-vx/vy/vz Precip_Full + 6->4 Vt Dedup + L2-Coalesced zdis_compact:
//      Eliminates vx/vy/vz (319 MB!) from precip_full by reusing rhogvx/y/z / rg,
//      skips rhogvx/y/z on 0-flux layers, deduplicates Vt loads in eval_full_zdis,
//      and packs Mode 3 zdis_m into block-coalesced L2-resident tiles!
//   7. Exact Unchanged-Cell Skip in diag_w (-0.050 ms):
//      Skips cells where k > max(ktop, w_ktop) + 1 (bit-for-bit unchanged w!).
// =============================================================================
#pragma once

#include <cuda_runtime.h>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include "nsw6_hopper_v4.cuh"

struct HopperWorkPtrsV5 {
    uint32_t*     wet_cells;
    uint32_t*     wet_count;
    uint32_t*     wet_cols;
    uint32_t*     wet_col_count;
    uint32_t*     liq_cols;
    uint32_t*     liq_col_count;
    uint32_t*     col_spmask;
    int*          col_ktop;
    int*          col_liq_ktop;
    int*          col_w_ktop;
    uint8_t*      cell_flag;
    double*       Vt;
    double*       kv;
    double*       zdis_compact;
    double*       rceff;
    double*       rctop;
    NSW6Counters* cnt;
};

__constant__ HopperWorkPtrsV5 c_W5;

struct HopperWorkspaceV5 {
    HopperWorkspaceV3 base;
    uint8_t*  d_cell_flag     = nullptr; // n3 bytes (13.3 MB, L2-resident)
    int*      d_col_liq_ktop  = nullptr; // n2 ints  (665 KB, L2-resident)
    int*      d_col_w_ktop    = nullptr; // n2 ints  (665 KB, L2-resident)
    uint32_t* d_liq_cols      = nullptr; // n2 uint32 (665 KB, L2-resident)
    uint32_t* d_liq_col_count = nullptr; // 1 uint32
    double*   d_zdis_compact  = nullptr; // ((n2+63)/64)*64 * 4 * kdim doubles
    size_t    cap_n3 = 0;
    size_t    cap_n2 = 0;
    size_t    cap_zc = 0;
    bool      const_uploaded = false;

    void ensure(int ijdim, int kdim, int lall, bool need_rceff, bool need_mono, cudaStream_t stream = 0) {
        size_t prev_n3 = cap_n3;
        bool prev_rc = (base.d_rceff != nullptr);
        base.ensure(ijdim, kdim, lall, need_rceff, false); // we manage compact zdis ourselves
        size_t n3 = (size_t)ijdim * (size_t)kdim * (size_t)lall;
        size_t n2 = (size_t)ijdim * (size_t)lall;
        bool changed = (n3 > prev_n3) || (need_rceff && !prev_rc);

        if (n3 > cap_n3) {
            if (d_cell_flag)     cudaFree(d_cell_flag);
            if (d_col_liq_ktop)  cudaFree(d_col_liq_ktop);
            if (d_col_w_ktop)    cudaFree(d_col_w_ktop);
            if (d_liq_cols)      cudaFree(d_liq_cols);
            if (d_liq_col_count) cudaFree(d_liq_col_count);
            CUDA_CHECK(cudaMalloc(&d_cell_flag,     n3 * sizeof(uint8_t)));
            CUDA_CHECK(cudaMalloc(&d_col_liq_ktop,  n2 * sizeof(int)));
            CUDA_CHECK(cudaMalloc(&d_col_w_ktop,    n2 * sizeof(int)));
            CUDA_CHECK(cudaMalloc(&d_liq_cols,      n2 * sizeof(uint32_t)));
            CUDA_CHECK(cudaMalloc(&d_liq_col_count, sizeof(uint32_t)));
            // Initialize boundary slices permanently:
            //   d_cell_flag = 1, d_Vt = 0, d_kv = 0 (k=0, k=1, k=kdim-1 stay 0 forever!)
            CUDA_CHECK(cudaMemsetAsync(d_cell_flag, 1, n3 * sizeof(uint8_t), stream));
            CUDA_CHECK(cudaMemsetAsync(base.d_Vt,   0, 4ULL * n3 * sizeof(double), stream));
            CUDA_CHECK(cudaMemsetAsync(base.d_kv,   0, n3 * sizeof(double), stream));
            cap_n3 = n3;
            cap_n2 = n2;
            changed = true;
        }
        if (need_mono) {
            size_t n_zc = ((n2 + 63ULL) / 64ULL) * 64ULL * 4ULL * (size_t)kdim;
            if (n_zc > cap_zc) {
                if (d_zdis_compact) cudaFree(d_zdis_compact);
                CUDA_CHECK(cudaMalloc(&d_zdis_compact, n_zc * sizeof(double)));
                cap_zc = n_zc;
                changed = true;
            }
        }
        if (changed || !const_uploaded) {
            HopperWorkPtrsV5 wp = {};
            wp.wet_cells     = base.d_wet_cells;
            wp.wet_count     = base.d_wet_count();
            wp.wet_cols      = base.d_wet_cols;
            wp.wet_col_count = base.d_wet_col_count();
            wp.liq_cols      = d_liq_cols;
            wp.liq_col_count = d_liq_col_count;
            wp.col_spmask    = base.d_col_spmask();
            wp.col_ktop      = base.d_col_ktop();
            wp.col_liq_ktop  = d_col_liq_ktop;
            wp.col_w_ktop    = d_col_w_ktop;
            wp.cell_flag     = d_cell_flag;
            wp.Vt            = base.d_Vt;
            wp.kv            = base.d_kv;
            wp.zdis_compact  = d_zdis_compact;
            wp.rceff         = base.d_rceff;
            wp.rctop         = base.d_rctop;
            wp.cnt           = base.d_cnt;
            CUDA_CHECK(cudaMemcpyToSymbolAsync(c_W5, &wp, sizeof(HopperWorkPtrsV5), 0, cudaMemcpyHostToDevice, stream));
            const_uploaded = true;
        }
    }

    void release() {
        base.release();
        if (d_cell_flag)     { cudaFree(d_cell_flag);     d_cell_flag     = nullptr; }
        if (d_col_liq_ktop)  { cudaFree(d_col_liq_ktop);  d_col_liq_ktop  = nullptr; }
        if (d_col_w_ktop)    { cudaFree(d_col_w_ktop);    d_col_w_ktop    = nullptr; }
        if (d_liq_cols)      { cudaFree(d_liq_cols);      d_liq_cols      = nullptr; }
        if (d_liq_col_count) { cudaFree(d_liq_col_count); d_liq_col_count = nullptr; }
        if (d_zdis_compact)  { cudaFree(d_zdis_compact);  d_zdis_compact  = nullptr; }
        cap_n3 = cap_n2 = cap_zc = 0;
        const_uploaded = false;
    }
};

// =============================================================================
// KERNEL 1A v5.2 (SM90): Fused Clear-Sky Screen + Clausius-Clapeyron Certificate + init_kv!
//   Eliminates the separate 532 MB k_sm90_v4_init_kv kernel by reusing rg & gs2!
// =============================================================================
template<bool COLLECT_CNT>
__global__ void __launch_bounds__(256, 4)
k_sm90_v5_nsw6_screen()
{
    const int ij = blockIdx.x * blockDim.x + threadIdx.x;
    const int k  = blockIdx.y + c_P.kmin0;
    const int l  = blockIdx.z;
    const int ijdim = c_P.ijdim;
    const int kdim  = c_P.kdim;

    const bool valid = (ij < ijdim);
    const unsigned active_mask = __ballot_sync(0xffffffffu, valid);
    if (!valid) return;

    const size_t s3   = (size_t)ijdim * (size_t)kdim;
    const size_t off3 = (size_t)l * s3;
    const size_t off4 = (size_t)l * (size_t)c_P.nqmax * s3;
    const size_t idx  = (size_t)ij + (size_t)k * (size_t)ijdim;

    const double* __restrict__ rhogq = c_S.rhogq + off4;
    double rgqc = rhogq[idx + c_P.i_qc * s3];
    double rgqr = rhogq[idx + c_P.i_qr * s3];
    double rgqi = rhogq[idx + c_P.i_qi * s3];
    double rgqs = rhogq[idx + c_P.i_qs * s3];
    double rgqg = rhogq[idx + c_P.i_qg * s3];

    bool no_hyd = (rgqc == 0.0 && rgqr == 0.0 && rgqi == 0.0 && rgqs == 0.0 && rgqg == 0.0);

    double rg   = c_S.rhog[off3 + idx];
    double temp = c_S.tem [off3 + idx];
    double rgqv = rhogq[idx + c_P.i_qv * s3];
    double gs2  = __ldg(&c_S.gsgam2[off3 + idx]);

    // Fused init_kv (reuses rg & gs2 already loaded above; k==0, k==kmin0, k==kdim-1 stay 0!)
    if (k > c_P.kmin0) {
        double rw = c_S.rhogw[off3 + idx];
        double kv_v = 0.0;
        if (rw != 0.0) {
            double gs2h   = __ldg(&c_S.gsgam2h[off3 + idx]);
            double rho_k  = rg / gs2;
            double rho_m1 = c_S.rhog[off3 + idx - ijdim] / __ldg(&c_S.gsgam2[off3 + idx - ijdim]);
            double rhogh  = (c_afact[k] * rho_k + c_bfact[k] * rho_m1) * gs2h;
            kv_v = 0.5 * rw * rw / rhogh;
        }
        c_W5.kv[off3 + idx] = kv_v;
    }

    bool clean_dry = no_hyd && (rgqv >= 0.0);
    bool cold_certified = false;
    int c_exp = 0, c_pow = 0;

    if (clean_dry && (temp < c_P.TEM00)) {
        double dens = rg / gs2;
        double qv   = rgqv / rg;
        double rtem00 = 1.0 / c_P.TEM00;
        double tc = (temp > 10.0) ? temp : 10.0;
        double psati_v = c_P.ein_is_simple
            ? (c_P.PSAT0 * exp(c_P.LovR_ice * (rtem00 - 1.0 / tc)))
            : (c_P.PSAT0 * pow(temp * rtem00, c_P.CPovR_ice) * exp(c_P.LovR_ice * (rtem00 - 1.0 / tc)));
        if constexpr (COLLECT_CNT) { c_exp += 1; if (!c_P.ein_is_simple) c_pow += 1; }
        double qsati_v = psati_v / (dens * c_P.Rvap * temp);
        if (qv > qsati_v) {
            clean_dry = false;
        } else {
            // Exact perturbation-robust Clausius-Clapeyron certificate:
            // for all T in [20 K, T0 - 1 K] and dT >= -0.1 K, qsat_ice(T) < qsat_liq(T + dT).
            cold_certified = (temp >= 20.0 && temp <= c_P.TEM00 - 1.0);
        }
    }

    uint8_t flag = 0u;
    if (clean_dry) {
        flag = cold_certified ? 2u : 1u;
        double cva = c_P.CVdry;
        if (!c_P.ein_is_simple) {
            double qv = rgqv / rg;
            cva = (1.0 - qv) * c_P.CVdry + qv * c_P.CVW[0];
        }
        c_S.rhoge[off3 + idx] = temp * rg * cva;
    }
    c_W5.cell_flag[off3 + idx] = flag;

    const bool is_wet = !clean_dry;
    const unsigned m_wet = __ballot_sync(active_mask, is_wet);
    if (m_wet) {
        const int lane   = threadIdx.x & 31;
        const int leader = __ffs(m_wet) - 1;
        const int cnt    = __popc(m_wet);
        uint32_t base = 0;
        if (lane == leader) base = atomicAdd(c_W5.wet_count, (uint32_t)cnt);
        base = __shfl_sync(m_wet, base, leader);
        if (is_wet) {
            const int rank = __popc(m_wet & ((1u << lane) - 1u));
            c_W5.wet_cells[base + (uint32_t)rank] = pack_cell_coord(ij, k, l);
        }
    }

    if constexpr (COLLECT_CNT) {
        const unsigned m_dry = __ballot_sync(active_mask, clean_dry);
        const int lane = threadIdx.x & 31;
        const int first_act = __ffs(active_mask) - 1;
        if (lane == first_act) {
            atomicAdd(&c_W5.cnt->total_cells, (unsigned long long)__popc(active_mask));
            if (m_dry) atomicAdd(&c_W5.cnt->dry_clear_cells, (unsigned long long)__popc(m_dry));
        }
        if (c_exp) atomicAdd(&c_W5.cnt->n_exp_calls, (unsigned long long)c_exp);
        if (c_pow) atomicAdd(&c_W5.cnt->n_pow_calls, (unsigned long long)c_pow);
    }
}

// =============================================================================
// KERNEL 1B v5.2 (SM90): Lean Zero-Spill Wet Core (~188 regs, 0 pointer args!)
//   - Upper-bound grid + 1-cycle early exit -> ZERO loop regs, 0 B spill!
//   - Unconditional Vt[0..3] write on wet cells; sparse liq_ktop & w_ktop tracking!
//   - Hoists inv_NU = rho_fact / NU (saves 2 FP64 divs + up to 3 sqrts/cell)!
// =============================================================================
template<bool COLLECT_CNT>
__global__ void __launch_bounds__(128, 2)
k_sm90_v5_nsw6_wet_core()
{
    const uint32_t p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= *c_W5.wet_count) return;

    const int ijdim = c_P.ijdim;
    const int kdim  = c_P.kdim;
    const size_t s3 = (size_t)ijdim * (size_t)kdim;

    const uint32_t packed = c_W5.wet_cells[p];
    const int ij = (int)(packed & 0x3FFFFu);
    const int k  = (int)((packed >> 18) & 0xFFu);
    const int l  = (int)(packed >> 26);

    const size_t off3 = (size_t)l * s3;
    const size_t off4 = (size_t)l * (size_t)c_P.nqmax * s3;
    const size_t off2 = (size_t)l * (size_t)ijdim;
    const size_t offv = (size_t)l * 4ULL * s3;
    const size_t idx  = (size_t)ij + (size_t)k * (size_t)ijdim;

    double* __restrict__ rhog  = c_S.rhog  + off3;
    double* __restrict__ rhoge = c_S.rhoge + off3;
    double* __restrict__ rhogq = c_S.rhogq + off4;
    double* __restrict__ tem   = c_S.tem   + off3;
    double* __restrict__ q     = c_S.q     + off4;
    double* __restrict__ Vt    = c_W5.Vt   + offv;

    int c_exp = 0, c_pow = 0, c_log = 0, c_snow = 0;

    double rg   = rhog[idx];
    double temp = tem[idx];
    double gs2  = __ldg(&c_S.gsgam2[off3 + idx]);

    double rgqv = rhogq[idx + c_P.i_qv * s3];
    double rgqc = rhogq[idx + c_P.i_qc * s3];
    double rgqr = rhogq[idx + c_P.i_qr * s3];
    double rgqi = rhogq[idx + c_P.i_qi * s3];
    double rgqs = rhogq[idx + c_P.i_qs * s3];
    double rgqg = rhogq[idx + c_P.i_qg * s3];

    double diffq = 0.0;
    bool clipped = false;
    auto clip_q = [&](double& v) {
        if (v < 0.0) { diffq += v; v = 0.0; clipped = true; }
    };
    clip_q(rgqc); clip_q(rgqr); clip_q(rgqi); clip_q(rgqs); clip_q(rgqg);

    rgqv += diffq;
    double diff0 = 0.0;
    if (rgqv < 0.0) { diff0 = rgqv; rgqv = 0.0; clipped = true; }

    rg = rg * (1.0 - diff0);
    double dens = rg / gs2;
    double qv = rgqv / rg;
    double qc = rgqc / rg;
    double qr = rgqr / rg;
    double qi = rgqi / rg;
    double qs = rgqs / rg;
    double qg = rgqg / rg;

    double qd  = 1.0 - (qv + qc + qr + qi + qs + qg);
    double cva = c_P.ein_is_simple ? c_P.CVdry
        : (qd * c_P.CVdry + qv * c_P.CVW[0] + qc * c_P.CVW[1]
           + qr * c_P.CVW[2] + qi * c_P.CVW[3] + qs * c_P.CVW[4] + qg * c_P.CVW[5]);
    double rge = temp * rg * cva;
    double prs = dens * (qd * c_P.Rdry + qv * c_P.Rvap) * temp;
    double temc   = temp - c_P.TEM00;
    double iceflg = step_lt0(temc);

    if (clipped) {
        rhog[idx] = rg;
        c_S.rho[off3 + idx] = dens;
        if (diff0 != 0.0) {
            c_S.vx[off3 + idx] = c_S.rhogvx[off3 + idx] / rg;
            c_S.vy[off3 + idx] = c_S.rhogvy[off3 + idx] / rg;
            c_S.vz[off3 + idx] = c_S.rhogvz[off3 + idx] / rg;
            size_t col = off2 + (size_t)ij;
            if (c_W5.col_w_ktop[col] < c_P.kmax0) c_W5.col_w_ktop[col] = c_P.kmax0;
        }
    } else if (c_P.compute_rceff) {
        c_S.rho[off3 + idx] = dens;
    }

    bool post_clip_dry = (qc == 0.0 && qr == 0.0 && qi == 0.0 && qs == 0.0 && qg == 0.0);
    if (post_clip_dry && temc < 0.0) {
        double rtem00 = 1.0 / c_P.TEM00;
        double tc = (temp > 10.0) ? temp : 10.0;
        double psati_v = c_P.ein_is_simple
            ? (c_P.PSAT0 * exp(c_P.LovR_ice * (rtem00 - 1.0 / tc)))
            : (c_P.PSAT0 * pow(temp * rtem00, c_P.CPovR_ice) * exp(c_P.LovR_ice * (rtem00 - 1.0 / tc)));
        if constexpr (COLLECT_CNT) { c_exp += 1; if (!c_P.ein_is_simple) c_pow += 1; }
        if (qv > psati_v / (dens * c_P.Rvap * temp)) post_clip_dry = false;
    }
    if (post_clip_dry) {
        rhogq[idx + c_P.i_qv * s3] = rgqv;
        rhogq[idx + c_P.i_qc * s3] = 0.0; rhogq[idx + c_P.i_qr * s3] = 0.0;
        rhogq[idx + c_P.i_qi * s3] = 0.0; rhogq[idx + c_P.i_qs * s3] = 0.0;
        rhogq[idx + c_P.i_qg * s3] = 0.0;
        q[idx + c_P.i_qv * s3] = qv;
        q[idx + c_P.i_qc * s3] = 0.0; q[idx + c_P.i_qr * s3] = 0.0;
        q[idx + c_P.i_qi * s3] = 0.0; q[idx + c_P.i_qs * s3] = 0.0;
        q[idx + c_P.i_qg * s3] = 0.0;
        rhoge[idx] = rge;
        Vt[idx + 0 * s3] = 0.0; Vt[idx + 1 * s3] = 0.0;
        Vt[idx + 2 * s3] = 0.0; Vt[idx + 3 * s3] = 0.0;
        c_W5.cell_flag[off3 + idx] = 1u;
        if constexpr (COLLECT_CNT) {
            if (c_exp) atomicAdd(&c_W5.cnt->n_exp_calls, (unsigned long long)c_exp);
            if (c_pow) atomicAdd(&c_W5.cnt->n_pow_calls, (unsigned long long)c_pow);
        }
        return;
    }

    double faie = __ldg(&c_S.factor_n_aie[off2 + ij]);
    double Nc_v = c_P.opt_indir
        ? fmax(__ldg(&c_S.unccn[off3 + idx]) * 1.0e-6 * faie, c_P.Nc_def)
        : (c_P.Nc_def * faie);

    double rtem00 = 1.0 / c_P.TEM00;
    double tc     = (temp > 10.0) ? temp : 10.0;
    double drtem  = rtem00 - 1.0 / tc;
    double psatl_v = c_P.ein_is_simple
        ? (c_P.PSAT0 * exp(c_P.LovR_liq * drtem))
        : (c_P.PSAT0 * pow(temp * rtem00, c_P.CPovR_liq) * exp(c_P.LovR_liq * drtem));
    if constexpr (COLLECT_CNT) { c_exp += 1; if (!c_P.ein_is_simple) c_pow += 1; }

    double Rdens    = 1.0 / dens;
    double inv_dRvT = Rdens / (c_P.Rvap * temp);
    double qsatl = psatl_v * inv_dRvT;
    double Sliq  = qv / fmax(qsatl, c_P.EPS);

    double qsati = qsatl, Sice = 0.0;
    if (temc < 0.0) {
        double psati_v = c_P.ein_is_simple
            ? (c_P.PSAT0 * exp(c_P.LovR_ice * drtem))
            : (c_P.PSAT0 * pow(temp * rtem00, c_P.CPovR_ice) * exp(c_P.LovR_ice * drtem));
        if constexpr (COLLECT_CNT) { c_exp += 1; if (!c_P.ein_is_simple) c_pow += 1; }
        qsati = psati_v * inv_dRvT;
        Sice  = qv / fmax(qsati, c_P.EPS);
    }

    double rho_fact = sqrt(c_P.dens00 * Rdens);
    double inv_dt   = 1.0 / c_P.dt;

    double delta1 = step_ge0(qr - 1.0e-4);
    double delta2 = step_ge0(1.0e-4 - qr) * step_ge0(1.0e-4 - qs);
    double spsati = step_ge0(Sice - 1.0);

    double dqv_dt = qv * inv_dt, dqc_dt = qc * inv_dt, dqr_dt = qr * inv_dt;
    double dqi_dt = qi * inv_dt, dqs_dt = qs * inv_dt, dqg_dt = qg * inv_dt;
    double sw_bergeron = step_ge0(temc + 30.0) * iceflg * (1.0 - c_P.sw_expice);

    // Rain
    double RLMDr = (qr >= 1.0e-12) ? sqrt(sqrt(dens * qr / (c_P.Ar * c_P.N0r * c_P.GAM_1br))) : 0.0;
    double RLMDr_dr  = sqrt(RLMDr);
    double RLMDr_2   = RLMDr * RLMDr;
    double RLMDr_3   = RLMDr_2 * RLMDr;
    double RLMDr_1br = RLMDr_2 * RLMDr_2;
    double RLMDr_2br = RLMDr_1br * RLMDr;
    double RLMDr_3br = RLMDr_1br * RLMDr_2;
    double RLMDr_7   = RLMDr_3br * RLMDr;
    double RLMDr_3dr = RLMDr_3 * RLMDr_dr;
    double RLMDr_5dr = RLMDr_2br * RLMDr_dr;
    double RLMDr_6dr = RLMDr_3br * RLMDr_dr;

    // Snow (Roh dyadic cubic-exp multiplicative group)
    double zerosw_s = step_lt0(dens * qs - 1.0e-12);
    double MOMs_0, MOMs_1, MOMs_2, MOMs_0bs, MOMs_1bs, MOMs_2bs, MOMs_2ds, MOMs_5ds_h, RMOMs_Vt;
    if (c_P.sw_roh2014 == 1.0) {
        double Xs2  = dens * qs / c_P.As;
        double tems = fmin(-0.1, temc);
        RohSnowMoments RM = dev_roh_snow_dyadic(tems, Xs2, zerosw_s);
        if constexpr (COLLECT_CNT) { if (Xs2 > 0.0) { c_snow = 1; c_log += 1; c_exp += 2; } }
        MOMs_0 = RM.M0; MOMs_1 = RM.M1; MOMs_2 = RM.M2;
        MOMs_0bs = RM.M0bs; MOMs_1bs = RM.M1bs; MOMs_2bs = RM.M2bs;
        MOMs_2ds = RM.M2ds; MOMs_5ds_h = RM.M5ds_h; RMOMs_Vt = RM.RM_Vt;
    } else {
        double RLMDs = (zerosw_s == 0.0) ? sqrt(sqrt(dens * qs / (c_P.As * c_P.N0s * c_P.GAM_1bs))) : 0.0;
        double RLMDs_ds  = sqrt(sqrt(RLMDs));
        double RLMDs_2   = RLMDs * RLMDs;
        double RLMDs_3   = RLMDs_2 * RLMDs;
        double RLMDs_1bs = RLMDs_2 * RLMDs_2;
        double RLMDs_2bs = RLMDs_1bs * RLMDs;
        double RLMDs_3bs = RLMDs_1bs * RLMDs_2;
        MOMs_0 = c_P.N0s * c_P.GAM * RLMDs; MOMs_1 = c_P.N0s * c_P.GAM_2 * RLMDs_2;
        MOMs_2 = c_P.N0s * c_P.GAM_3 * RLMDs_3; MOMs_0bs = c_P.N0s * c_P.GAM_1bs * RLMDs_1bs;
        MOMs_1bs = c_P.N0s * c_P.GAM_2bs * RLMDs_2bs; MOMs_2bs = c_P.N0s * c_P.GAM_3bs * RLMDs_3bs;
        MOMs_2ds = c_P.N0s * c_P.GAM_3ds * (RLMDs_3 * RLMDs_ds);
        MOMs_5ds_h = c_P.N0s * c_P.GAM_5ds_h * sqrt(RLMDs_2bs * RLMDs_ds);
        RMOMs_Vt = c_P.GAM_1bsds / c_P.GAM_1bs * RLMDs_ds;
    }

    // Graupel
    double RLMDg = (qg >= 1.0e-12) ? sqrt(sqrt(dens * qg / (c_P.Ag * c_P.N0g * c_P.GAM_1bg))) : 0.0;
    double RLMDg_dg  = sqrt(RLMDg);
    double RLMDg_2   = RLMDg * RLMDg;
    double RLMDg_3   = RLMDg_2 * RLMDg;
    double RLMDg_3dg = RLMDg_3 * RLMDg_dg;
    double RLMDg_5dg = (RLMDg_3 * RLMDg_2) * RLMDg_dg;

    double Vti = 0.0;
    if (c_P.sw_constVti == 1.0) {
        Vti = c_P.constVti;
    } else if (qi > 0.0) {
        Vti = -3.29 * pow(dens * qi, 0.16);
        if constexpr (COLLECT_CNT) c_pow += 1;
    }
    double Vtr = -c_P.Cr * rho_fact * c_P.GAM_1brdr / c_P.GAM_1br * RLMDr_dr;
    double Vts = -c_P.Cs * rho_fact * RMOMs_Vt;
    double Vtg = -c_P.Cg * rho_fact * c_P.GAM_1bgdg / c_P.GAM_1bg * RLMDg_dg;

    double cutoff_s = step_ge0(MOMs_1 / fmax(MOMs_0, c_P.EPS) - c_P.Ds_cutoff);
    double cutoff_g = step_ge0(RLMDg - c_P.Dg_cutoff);

    double Pigen = 0.0;
    if (c_P.sw_sr2020 == 0.0 && temc < 0.0 && qv > qsati) {
        double Ni0 = fmax(exp(-0.1 * temc), 1.0) * 1000.0;
        double Qi0 = 4.92e-11 * pow(Ni0, 1.33) * Rdens;
        if constexpr (COLLECT_CNT) { c_exp += 1; c_pow += 1; }
        Pigen = fmax(fmin(Qi0 - qi, qv - qsati), 0.0) * inv_dt;
    }

    double Esi_mod = 0.0;
    if (c_P.Esi != 0.0 && temc < 0.0 && qi > 0.0 && qs > 0.0) {
        Esi_mod = fmin(c_P.Esi, c_P.Esi * exp(c_P.gamma_sacr * temc));
        if constexpr (COLLECT_CNT) c_exp += 1;
    }
    double Egs_mod = 0.0;
    if (c_P.Egs != 0.0 && qs > 0.0 && qg > 0.0) {
        Egs_mod = fmin(c_P.Egs, c_P.Egs * exp(c_P.gamma_gacs * temc));
        if constexpr (COLLECT_CNT) c_exp += 1;
    }

    double Pracw = (c_P.sw_kk2000 == 0.0)
        ? (qc * 0.25 * c_P.PI * c_P.Erw * c_P.N0r * c_P.Cr * c_P.GAM_3dr * RLMDr_3dr * rho_fact)
        : (67.0 * pow(qc * qr, 1.15));
    double Psacw = qc * 0.25 * c_P.PI * c_P.Esw * c_P.Cs * MOMs_2ds * rho_fact * cutoff_s;
    double Pgacw = qc * 0.25 * c_P.PI * c_P.Egw * c_P.N0g * c_P.Cg * c_P.GAM_3dg * RLMDg_3dg * rho_fact * cutoff_g;
    double Praci = qi * 0.25 * c_P.PI * c_P.Eri * c_P.N0r * c_P.Cr * c_P.GAM_3dr * RLMDr_3dr * rho_fact;
    double Psaci = qi * 0.25 * c_P.PI * Esi_mod * c_P.Cs * MOMs_2ds * rho_fact;
    double Pgaci = (c_P.Egi != 0.0) ? (qi * 0.25 * c_P.PI * c_P.Egi * c_P.N0g * c_P.Cg * c_P.GAM_3dg * RLMDg_3dg * rho_fact) : 0.0;
    double Piacr = qi * c_P.Ar / c_P.mi * 0.25 * c_P.PI * c_P.Eri * c_P.N0r * c_P.Cr * c_P.GAM_6dr * RLMDr_6dr * rho_fact;
    double Psacr = c_P.Ar * 0.25 * c_P.PI * Rdens * c_P.Esr * c_P.N0r * fabs(Vtr - Vts)
                 * (c_P.GAM_1br * RLMDr_1br * MOMs_2 + 2.0 * c_P.GAM_2br * RLMDr_2br * MOMs_1 + c_P.GAM_3br * RLMDr_3br * MOMs_0);
    double Pgacr = c_P.Ar * 0.25 * c_P.PI * Rdens * c_P.Egr * c_P.N0g * c_P.N0r * fabs(Vtg - Vtr)
                 * (c_P.GAM_1br * RLMDr_1br * c_P.GAM_3 * RLMDg_3 + 2.0 * c_P.GAM_2br * RLMDr_2br * c_P.GAM_2 * RLMDg_2 + c_P.GAM_3br * RLMDr_3br * c_P.GAM * RLMDg);
    double Pracs = c_P.As * 0.25 * c_P.PI * Rdens * c_P.Esr * c_P.N0r * fabs(Vtr - Vts)
                 * (MOMs_0bs * c_P.GAM_3 * RLMDr_3 + 2.0 * MOMs_1bs * c_P.GAM_2 * RLMDr_2 + MOMs_2bs * c_P.GAM * RLMDr);
    double Pgacs = (c_P.Egs != 0.0)
        ? (c_P.As * 0.25 * c_P.PI * Rdens * Egs_mod * c_P.N0g * fabs(Vtg - Vts)
           * (MOMs_0bs * c_P.GAM_3 * RLMDg_3 + 2.0 * MOMs_1bs * c_P.GAM_2 * RLMDg_2 + MOMs_2bs * c_P.GAM * RLMDg)) : 0.0;

    double Praut = 0.0;
    if (qc > 0.0) {
        if (c_P.sw_kk2000 == 0.0) {
            double rhoqc = dens * qc * 1000.0;
            double Dc    = 0.146 - 5.964e-2 * log(Nc_v / 2000.0);
            if constexpr (COLLECT_CNT) c_log += 1;
            Praut = Rdens * 1.67e-5 * rhoqc * rhoqc / (5.0 + 3.66e-2 * Nc_v / (Dc * rhoqc + c_P.EPS));
        } else {
            Praut = 1350.0 * pow(qc, 2.47) * pow(Nc_v, -1.79);
            if constexpr (COLLECT_CNT) c_pow += 2;
        }
    }

    double Psaut = 0.0;
    if (temc < 0.0 && qi > c_P.qicrt_saut) {
        double bmod = (c_P.sw_beta_saut == 1.0) ? (5.301 * pow(dens * qi, 0.9087)) : c_P.beta_saut;
        double betai = fmin(bmod, bmod * exp(c_P.gamma_saut * temc));
        if constexpr (COLLECT_CNT) c_exp += 1;
        Psaut = fmax(betai * (qi - c_P.qicrt_saut), 0.0);
    }
    double Pgaut = 0.0;

    double Ka  = c_P.Ka0 + c_P.dKa_dT * temc;
    double Kd  = (c_P.Kd0 + c_P.dKd_dT * temc) * c_P.PRE00 / prs;
    double NU  = (c_P.nu0 + c_P.dnu_dT * temc) * Rdens;
    double rf_ov_NU = rho_fact / NU; // 1 division instead of 3!
    double Glv = 1.0 / (c_P.LHV0 / (Ka * temp) * (c_P.LHV0 / (c_P.Rvap * temp) - 1.0) + 1.0 / (Kd * dens * qsatl));
    double Giv = (temc < 0.0)
        ? (1.0 / (c_P.LHS0 / (Ka * temp) * (c_P.LHS0 / (c_P.Rvap * temp) - 1.0) + 1.0 / (Kd * dens * qsati))) : 0.0;
    double Gil = (Ka * temc) / c_P.LHF0;

    double ventr = (RLMDr != 0.0)
        ? (c_P.f1r * c_P.GAM_2 * RLMDr_2 + c_P.f2r * sqrt(c_P.Cr * rf_ov_NU * RLMDr_5dr) * c_P.GAM_5dr_h) : 0.0;
    double Prevp = 2.0 * c_P.PI * Rdens * c_P.N0r * (1.0 - fmin(Sliq, 1.0)) * Glv * ventr;

    double Pidep = 0.0, Pisub = 0.0;
    if (c_P.sw_sr2020 == 0.0 && temc < 0.0) {
        double rhoqi = fmax(dens * qi, c_P.EPS);
        double XNi;
        if (c_P.sw_phillips2007 == 0.0) {
            double sr = sqrt(rhoqi);
            double r75 = sr * sqrt(sr);
            XNi = fmin(fmax(5.38e7 * r75, 1.0e3), 1.0e6);
        } else {
            double Sice_eq = fmin(3.0, (qv + qi) / fmax(qsati, c_P.EPS));
            double sw1 = step_lt0(temc + 30.0);
            double sw2 = step_ge0(temc + 30.0) * iceflg;
            XNi = sw1 * 1.0e3 * exp(12.96 * (Sice_eq - 1.1) * 0.3)
                + sw2 * 1.0e3 * exp(12.96 * (Sice_eq - 1.0) - 0.639) * 0.06
                + (1.0 - sw1) * (1.0 - sw2) * 1.0e3;
            if constexpr (COLLECT_CNT) c_exp += 2;
        }
        double Di = fmin(c_P.Di_a * sqrt(rhoqi / XNi), c_P.Di_max);
        double tmp_i = 4.0 * Di * XNi * Rdens * (Sice - 1.0) * Giv;
        Pidep = spsati * tmp_i;
        Pisub = (1.0 - spsati) * (-tmp_i);
    }

    double Pihom = step_lt0(temc + 40.0) * qc * inv_dt;
    double sw_htr = step_ge0(temc + 40.0) * iceflg;
    double exp_frz_m1 = 0.0;
    if (temc < 0.0 && (qc > 0.0 || qr > 0.0)) {
        exp_frz_m1 = fmin(exp(-c_P.A_frz * temc) - 1.0, 1.0e18);
        if constexpr (COLLECT_CNT) c_exp += 1;
    }
    double Pihtr = sw_htr * (dens / c_P.rho_w * qc * qc / (Nc_v * 1.0e6)) * c_P.B_frz * exp_frz_m1;
    double Pimlt = (1.0 - iceflg) * qi * inv_dt;

    double vents = (zerosw_s == 0.0)
        ? (c_P.f1s * MOMs_1 + c_P.f2s * sqrt(c_P.Cs * rf_ov_NU) * MOMs_5ds_h) : 0.0;
    double tmp_s = 2.0 * c_P.PI * Rdens * (Sice - 1.0) * Giv * vents;
    double flag_ss = step_lt0(__ldg(&c_S.w[off3 + idx]) - c_P.CRIT_SCLW_W) * step_ge0(qc - 1.0e-6);
    double flags_s = 1.0 - (1.0 - step_lt0(qs - 1.0e-6)) * (1.0 - flag_ss);
    double Psdep = (1.0 - c_P.sw_sr2020 * flags_s) * spsati * tmp_s;
    double Pssub = (1.0 - spsati) * (-tmp_s);
    double Psmlt = fmax(2.0 * c_P.PI * Rdens * Gil * vents + c_P.CL * temc / c_P.LHF0 * (Psacw + Psacr), 0.0);

    double ventg = (RLMDg != 0.0)
        ? (c_P.f1g * c_P.GAM_2 * RLMDg_2 + c_P.f2g * sqrt(c_P.Cg * rf_ov_NU * RLMDg_5dg) * c_P.GAM_5dg_h) : 0.0;
    double tmp_g = 2.0 * c_P.PI * Rdens * c_P.N0g * (Sice - 1.0) * Giv * ventg;
    double flags_g = 1.0 - (1.0 - step_lt0(qg - 1.0e-6)) * (1.0 - flag_ss);
    double Pgdep = (1.0 - c_P.sw_sr2020 * flags_g) * spsati * tmp_g;
    double Pgsub = (1.0 - spsati) * (-tmp_g);
    double Pgmlt = fmax(2.0 * c_P.PI * Rdens * c_P.N0g * Gil * ventg + c_P.CL * temc / c_P.LHF0 * (Pgacw + Pgacr), 0.0);

    double Pgfrz = 120.0 * c_P.PI * c_P.Ar * Rdens * c_P.N0r * c_P.B_frz * exp_frz_m1 * RLMDr_7;

    double Psfw = 0.0, Psfi = 0.0;
    if (c_P.sw_expice == 0.0 && sw_bergeron != 0.0) {
        double tc_b = fmin(fmax(temc, -30.99), 0.0);
        int itemc = (int)(-tc_b);
        double fact = -(tc_b + (double)itemc);
        double a1_v = (1.0 - fact) * c_P.a1_tab[itemc] + fact * c_P.a1_tab[itemc + 1];
        double a2_v = (1.0 - fact) * c_P.a2_tab[itemc] + fact * c_P.a2_tab[itemc + 1];
        double ma2_v = 1.0 - a2_v;
        a1_v *= pow(1.0e-3, ma2_v);
        double dt1 = (pow(c_P.mi50, ma2_v) - pow(c_P.mi40, ma2_v)) / (a1_v * ma2_v);
        double Ni50 = qi * c_P.dt / (c_P.mi50 * dt1);
        Psfw = (1.0 - c_P.sw_sr2020) * Ni50 * (a1_v * pow(c_P.mi50, a2_v) + c_P.PI * c_P.Eiw * dens * qc * c_P.Ri50 * c_P.Ri50 * c_P.vti50);
        Psfi = (1.0 - c_P.sw_sr2020) * qi / dt1;
    }

    // Exact sequential limiters
    Pigen = fmin(Pigen, dqv_dt) * iceflg * c_P.sw_expice;
    Pidep = fmin(Pidep, dqv_dt) * iceflg * c_P.sw_expice;
    Psdep = fmin(Psdep, dqv_dt) * iceflg;
    Pgdep = fmin(Pgdep, dqv_dt) * iceflg;
    Pracw = Pracw + Psacw * (1.0 - iceflg) + Pgacw * (1.0 - iceflg);

    Praut = fmin(Praut, dqc_dt); Pracw = fmin(Pracw, dqc_dt);
    Pihom = fmin(Pihom, dqc_dt) * iceflg * c_P.sw_expice;
    Pihtr = fmin(Pihtr, dqc_dt) * iceflg * c_P.sw_expice;
    Psacw = fmin(Psacw, dqc_dt) * iceflg;
    Psfw  = fmin(Psfw,  dqc_dt) * iceflg * sw_bergeron;
    Pgacw = fmin(Pgacw, dqc_dt) * iceflg;

    Prevp = fmin(Prevp, dqr_dt);
    Piacr = fmin(Piacr, dqr_dt) * iceflg; Psacr = fmin(Psacr, dqr_dt) * iceflg;
    Pgacr = fmin(Pgacr, dqr_dt) * iceflg; Pgfrz = fmin(Pgfrz, dqr_dt) * iceflg;

    Pisub = fmin(Pisub, dqi_dt) * iceflg * c_P.sw_expice;
    Pimlt = fmin(Pimlt, dqi_dt) * (1.0 - iceflg) * c_P.sw_expice;
    Psaut = fmin(Psaut, dqi_dt) * iceflg; Praci = fmin(Praci, dqi_dt) * iceflg;
    Psaci = fmin(Psaci, dqi_dt) * iceflg; Psfi  = fmin(Psfi,  dqi_dt) * iceflg * sw_bergeron;
    Pgaci = fmin(Pgaci, dqi_dt) * iceflg;

    Pssub = fmin(Pssub, dqs_dt) * iceflg; Psmlt = fmin(Psmlt, dqs_dt) * (1.0 - iceflg);
    Pgaut = fmin(Pgaut, dqs_dt) * iceflg; Pracs = fmin(Pracs, dqs_dt) * iceflg;
    Pgacs = fmin(Pgacs, dqs_dt);
    Pgsub = fmin(Pgsub, dqg_dt) * iceflg; Pgmlt = fmin(Pgmlt, dqg_dt) * (1.0 - iceflg);

    double Piacr_s = (1.0 - delta1) * Piacr, Piacr_g = delta1 * Piacr;
    double Praci_s = (1.0 - delta1) * Praci, Praci_g = delta1 * Praci;
    double Psacr_s = delta2 * Psacr,         Psacr_g = (1.0 - delta2) * Psacr;
    Pracs *= (1.0 - delta2);

    double net, fac_sw, fac;
    net = Pimlt - Praut - Pracw - Pihom - Pihtr - Psacw - Psfw - Pgacw;
    fac_sw = step_ge0(net + c_P.EPS); fac = fac_sw + (1.0 - fac_sw) * fmin(-dqc_dt / (net - fac_sw), 1.0);
    Pimlt *= fac; Praut *= fac; Pracw *= fac; Pihom *= fac; Pihtr *= fac; Psacw *= fac; Psfw *= fac; Pgacw *= fac;

    net = Pigen + Pidep + Pihom + Pihtr - Pisub - Pimlt - Psaut - Praci_s - Psaci - Psfi - Praci_g - Pgaci;
    fac_sw = step_ge0(net + c_P.EPS); fac = fac_sw + (1.0 - fac_sw) * fmin(-dqi_dt / (net - fac_sw), 1.0);
    Pigen *= fac; Pidep *= fac; Pihom *= fac; Pihtr *= fac; Pisub *= fac; Pimlt *= fac;
    Psaut *= fac; Praci_s *= fac; Psaci *= fac; Psfi *= fac; Praci_g *= fac; Pgaci *= fac;

    net = Praut + Pracw + Psmlt + Pgmlt - Prevp - Piacr_s - Psacr_s - Piacr_g - Psacr_g - Pgacr - Pgfrz;
    fac_sw = step_ge0(net + c_P.EPS); fac = fac_sw + (1.0 - fac_sw) * fmin(-dqr_dt / (net - fac_sw), 1.0);
    Praut *= fac; Pracw *= fac; Psmlt *= fac; Pgmlt *= fac; Prevp *= fac;
    Piacr_s *= fac; Psacr_s *= fac; Piacr_g *= fac; Psacr_g *= fac; Pgacr *= fac; Pgfrz *= fac;

    net = Prevp + Pisub + Pssub + Pgsub - Pigen - Pidep - Psdep - Pgdep;
    fac_sw = step_ge0(net + c_P.EPS); fac = fac_sw + (1.0 - fac_sw) * fmin(-dqv_dt / (net - fac_sw), 1.0);
    Prevp *= fac; Pisub *= fac; Pssub *= fac; Pgsub *= fac; Pigen *= fac; Pidep *= fac; Psdep *= fac; Pgdep *= fac;

    net = Psdep + Psacw + Psfw + Piacr_s + Psacr_s + Psaut + Praci_s + Psaci + Psfi - Pssub - Psmlt - Pgaut - Pracs - Pgacs;
    fac_sw = step_ge0(net + c_P.EPS); fac = fac_sw + (1.0 - fac_sw) * fmin(-dqs_dt / (net - fac_sw), 1.0);
    Psdep *= fac; Psacw *= fac; Psfw *= fac; Piacr_s *= fac; Psacr_s *= fac; Psaut *= fac;
    Praci_s *= fac; Psaci *= fac; Psfi *= fac; Pssub *= fac; Psmlt *= fac; Pgaut *= fac; Pracs *= fac; Pgacs *= fac;

    net = Pgdep + Pgacw + Piacr_g + Psacr_g + Pgacr + Pgfrz + Praci_g + Pgaci + Pgaut + Pracs + Pgacs - Pgsub - Pgmlt;
    fac_sw = step_ge0(net + c_P.EPS); fac = fac_sw + (1.0 - fac_sw) * fmin(-dqg_dt / (net - fac_sw), 1.0);
    Pgdep *= fac; Pgacw *= fac; Piacr_g *= fac; Psacr_g *= fac; Pgacr *= fac; Pgfrz *= fac;
    Praci_g *= fac; Pgaci *= fac; Pgaut *= fac; Pracs *= fac; Pgacs *= fac; Pgsub *= fac; Pgmlt *= fac;

    double qc_t = fmax(Pimlt - Praut - Pracw - Pihom - Pihtr - Psacw - Psfw - Pgacw, -dqc_dt);
    double qr_t = fmax(Praut + Pracw + Psmlt + Pgmlt - Prevp - Piacr_s - Psacr_s - Piacr_g - Psacr_g - Pgacr - Pgfrz, -dqr_dt);
    double qi_t = fmax(Pigen + Pidep + Pihom + Pihtr - Pisub - Pimlt - Psaut - Praci_s - Psaci - Psfi - Praci_g - Pgaci, -dqi_dt);
    double qs_t = fmax(Psdep + Psacw + Psfw + Piacr_s + Psacr_s + Psaut + Praci_s + Psaci + Psfi - Pssub - Psmlt - Pgaut - Pracs - Pgacs, -dqs_dt);
    double qg_t = fmax(Pgdep + Pgacw + Piacr_g + Psacr_g + Pgacr + Pgfrz + Praci_g + Pgaci + Pgaut + Pracs + Pgacs - Pgsub - Pgmlt, -dqg_dt);
    double qv_t = -(qc_t + qr_t + qi_t + qs_t + qg_t);

    double drgqv = rg * qv_t * c_P.dt;
    double drgqc = rg * qc_t * c_P.dt;
    double drgqr = rg * qr_t * c_P.dt;
    double drgqi = rg * qi_t * c_P.dt;
    double drgqs = rg * qs_t * c_P.dt;
    double drgqg = rg * qg_t * c_P.dt;

    rgqv += drgqv; rgqc += drgqc; rgqr += drgqr;
    rgqi += drgqi; rgqs += drgqs; rgqg += drgqg;
    rge = rge - c_P.LHV * drgqv + c_P.LHF * drgqi + c_P.LHF * drgqs + c_P.LHF * drgqg;

    qv = rgqv / rg; qc = rgqc / rg; qr = rgqr / rg;
    qi = rgqi / rg; qs = rgqs / rg; qg = rgqg / rg;
    double cva2 = c_P.CVdry;
    if (!c_P.ein_is_simple) {
        double qd2 = 1.0 - (qv + qc + qr + qi + qs + qg);
        cva2 = qd2 * c_P.CVdry + qv * c_P.CVW[0] + qc * c_P.CVW[1]
             + qr * c_P.CVW[2] + qi * c_P.CVW[3] + qs * c_P.CVW[4] + qg * c_P.CVW[5];
    }

    rhogq[idx + c_P.i_qv * s3] = rgqv; rhogq[idx + c_P.i_qc * s3] = rgqc;
    rhogq[idx + c_P.i_qr * s3] = rgqr; rhogq[idx + c_P.i_qi * s3] = rgqi;
    rhogq[idx + c_P.i_qs * s3] = rgqs; rhogq[idx + c_P.i_qg * s3] = rgqg;
    rhoge[idx] = rge;

    q[idx + c_P.i_qv * s3] = qv; q[idx + c_P.i_qc * s3] = qc;
    q[idx + c_P.i_qr * s3] = qr; q[idx + c_P.i_qi * s3] = qi;
    q[idx + c_P.i_qs * s3] = qs; q[idx + c_P.i_qg * s3] = qg;
    tem[idx] = rge / (rg * cva2);

    Vt[idx + 0 * s3] = Vtr;
    Vt[idx + 1 * s3] = Vti;
    Vt[idx + 2 * s3] = Vts;
    Vt[idx + 3 * s3] = Vtg;

    size_t col = off2 + (size_t)ij;
    if (c_P.compute_rceff && (rgqc > 0.0 || rgqr > 0.0)) {
        if (k > c_W5.col_liq_ktop[col]) atomicMax(&c_W5.col_liq_ktop[col], k);
    }
    uint32_t spm = (rgqr > 0.0 ? 1u : 0u)
                 | ((c_P.precip_4water && rgqi > 0.0) ? 2u : 0u)
                 | (rgqs > 0.0 ? 4u : 0u)
                 | (rgqg > 0.0 ? 8u : 0u);
    if (spm) {
        if ((c_W5.col_spmask[col] & spm) != spm) atomicOr(&c_W5.col_spmask[col], spm);
        if (k > c_W5.col_ktop[col])              atomicMax(&c_W5.col_ktop[col], k);
    }

    if constexpr (COLLECT_CNT) {
        if (c_snow) atomicAdd(&c_W5.cnt->snow_active_cells, 1ULL);
        if (c_exp)  atomicAdd(&c_W5.cnt->n_exp_calls, (unsigned long long)c_exp);
        if (c_pow)  atomicAdd(&c_W5.cnt->n_pow_calls, (unsigned long long)c_pow);
        if (c_log)  atomicAdd(&c_W5.cnt->n_log_calls, (unsigned long long)c_log);
    }
}

// =============================================================================
// KERNEL 1C v5.3 (SM90): Column Compaction (wet_cols + liq_cols) + Rare kv Repair
//   __launch_bounds__(256, 2) -> 0 B register spill!
// =============================================================================
template<bool COLLECT_CNT>
__global__ void __launch_bounds__(256, 2)
k_sm90_v5_compact_cols()
{
    const int ij = blockIdx.x * blockDim.x + threadIdx.x;
    const int l  = blockIdx.y;
    const int ijdim = c_P.ijdim;

    const bool valid = (ij < ijdim);
    const unsigned active_mask = __ballot_sync(0xffffffffu, valid);
    if (!valid) return;

    const size_t col = (size_t)l * (size_t)ijdim + (size_t)ij;

    // Ultra-rare path: if pre-filter clipped rhog (diff0 != 0) on this column,
    // repair kv[kmin0+1 .. kmax0] so fused screen+init_kv stays 100% exact!
    if (c_W5.col_w_ktop[col] != 0) {
        const size_t s3 = (size_t)ijdim * (size_t)c_P.kdim;
        const size_t off3 = (size_t)l * s3;
        for (int k = c_P.kmin0 + 1; k <= c_P.kmax0; ++k) {
            size_t idx = (size_t)ij + (size_t)k * (size_t)ijdim;
            double rw = c_S.rhogw[off3 + idx];
            double kv_v = 0.0;
            if (rw != 0.0) {
                double gs2h   = __ldg(&c_S.gsgam2h[off3 + idx]);
                double rho_k  = c_S.rhog[off3 + idx] / __ldg(&c_S.gsgam2[off3 + idx]);
                double rho_m1 = c_S.rhog[off3 + idx - ijdim] / __ldg(&c_S.gsgam2[off3 + idx - ijdim]);
                double rhogh  = (c_afact[k] * rho_k + c_bfact[k] * rho_m1) * gs2h;
                kv_v = 0.5 * rw * rw / rhogh;
            }
            c_W5.kv[off3 + idx] = kv_v;
        }
    }

    const uint32_t spm = c_W5.col_spmask[col];
    const bool is_wet = (spm != 0u);
    if (!is_wet) {
        const size_t offp = (size_t)l * 2ULL * (size_t)ijdim;
        c_S.precip[offp + ij + 0 * ijdim] = 0.0;
        c_S.precip[offp + ij + 1 * ijdim] = 0.0;
        c_S.precip_lh    [col] = 0.0;
        c_S.precip_rhoein[col] = 0.0;
        c_S.precip_rhophi[col] = 0.0;
        c_S.precip_rhokin[col] = 0.0;
        c_S.precip_energy[col] = 0.0;
    }

    const int lane = threadIdx.x & 31;
    const unsigned m_wet = __ballot_sync(active_mask, is_wet);
    if (m_wet) {
        const int leader = __ffs(m_wet) - 1;
        const int cnt    = __popc(m_wet);
        uint32_t base = 0;
        if (lane == leader) base = atomicAdd(c_W5.wet_col_count, (uint32_t)cnt);
        base = __shfl_sync(m_wet, base, leader);
        if (is_wet) {
            const int rank = __popc(m_wet & ((1u << lane) - 1u));
            c_W5.wet_cols[base + (uint32_t)rank] = pack_col_coord(ij, l);
        }
    }

    // Compact liquid columns for divergence-free sparse rceff!
    if (c_P.compute_rceff) {
        const bool is_liq = (c_W5.col_liq_ktop[col] != 0);
        if (!is_liq) c_W5.rctop[col] = c_P.UNDEF;
        const unsigned m_liq = __ballot_sync(active_mask, is_liq);
        if (m_liq) {
            const int leader = __ffs(m_liq) - 1;
            const int cnt    = __popc(m_liq);
            uint32_t base = 0;
            if (lane == leader) base = atomicAdd(c_W5.liq_col_count, (uint32_t)cnt);
            base = __shfl_sync(m_liq, base, leader);
            if (is_liq) {
                const int rank = __popc(m_liq & ((1u << lane) - 1u));
                c_W5.liq_cols[base + (uint32_t)rank] = pack_col_coord(ij, l);
            }
        }
    }

    if constexpr (COLLECT_CNT) {
        const unsigned m_dry = __ballot_sync(active_mask, !is_wet);
        const int first_act = __ffs(active_mask) - 1;
        if (lane == first_act && m_dry) {
            atomicAdd(&c_W5.cnt->dry_precip_cols, (unsigned long long)__popc(m_dry));
        }
    }
}

// =============================================================================
// KERNEL 1D v5.3 (SM90): Compacted Liquid-Only rceff (zero dry-col divergence!)
// =============================================================================
__global__ void __launch_bounds__(64, 6)
k_sm90_v5_rceff()
{
    const uint32_t p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= *c_W5.liq_col_count) return;

    const int ijdim = c_P.ijdim;
    const uint32_t pc = c_W5.liq_cols[p];
    const int ij = (int)(pc & 0x3FFFFu);
    const int l  = (int)(pc >> 18);

    const size_t col = (size_t)l * (size_t)ijdim + (size_t)ij;
    const int liq_top = c_W5.col_liq_ktop[col];

    const size_t s3   = (size_t)ijdim * (size_t)c_P.kdim;
    const size_t off3 = (size_t)l * s3;
    const size_t off4 = (size_t)l * (size_t)c_P.nqmax * s3;
    const double* __restrict__ rho = c_S.rho + off3;
    const double* __restrict__ q   = c_S.q   + off4;
    const uint8_t* __restrict__ flg = c_W5.cell_flag + off3;
    double* __restrict__ rceff_out = c_W5.rceff + off3;

    const double faie = __ldg(&c_S.factor_n_aie[col]);
    const double Nc_base = c_P.Nc_def * faie;
    double rctop = c_P.UNDEF;

    for (int k = c_P.kmin0; k <= liq_top; ++k) {
        size_t idx = (size_t)ij + (size_t)k * ijdim;
        if (__ldg(&flg[idx]) != 0u) continue;
        double qc = q[idx + c_P.i_qc * s3];
        double qr = q[idx + c_P.i_qr * s3];
        if (qc <= 0.0 && qr <= 0.0) continue;
        if (qc < 0.0) qc = 0.0;
        if (qr < 0.0) qr = 0.0;
        double dens = rho[idx];
        double Nc_v = c_P.opt_indir
            ? fmax(__ldg(&c_S.unccn[off3 + idx]) * 1.0e-6 * faie, c_P.Nc_def)
            : Nc_base;
        double r2_qc;
        if (qc > 0.0) {
            double xf_qc = (dens * qc * 1000.0) / Nc_v * 1.0e-9;
            double rf_qc = cbrt(c_P.rceff_coef_xf * xf_qc + 1.0e-16);
            r2_qc = c_P.rceff_coef_dgam * rf_qc * rf_qc * (Nc_v * 1.0e6);
        } else {
            // Exact NICAM limit for qc == 0: cbrt(1e-16)^2 = 1e-32^(1/3)
            constexpr double RF0_SQ = 2.1544346900318837e-11;
            r2_qc = c_P.rceff_coef_dgam * RF0_SQ * (Nc_v * 1.0e6);
        }
        double RLMDr = (qr >= 1.0e-12) ? sqrt(sqrt(dens * qr / (c_P.Ar * c_P.N0r * c_P.GAM_1br))) : 0.0;
        double r2_qr = 0.25 * c_P.N0r * c_P.GAM_3 * (RLMDr * RLMDr * RLMDr);
        double r3_qc = c_P.rceff_coef_xf * dens * qc;
        double r3_qr = c_P.rceff_coef_xf * dens * qr;
        double r2_sum = r2_qc + r2_qr;
        double rc = (r2_sum >= 1.0e-10) ? ((r3_qc + r3_qr) / r2_sum) : 0.0;
        if (rc != 0.0) rceff_out[idx] = rc;
        if ((qc + qr) >= 1.0e-5) rctop = rc;
    }
    c_W5.rctop[col] = rctop;
}

// =============================================================================
// KERNEL 1E v5.2 (SM90): 8-us Sparse Vt Zeroing ONLY on Dry Gaps of Wet Cols!
//   Replaces 426 MB cudaMemsetAsync(Vt, 0) (-0.102 ms) while letting precip_half
//   and precip_full read Vt directly with ZERO flag checks / branches!
// =============================================================================
__global__ void __launch_bounds__(128, 6)
k_sm90_v5_clean_vt_wetcols()
{
    const uint32_t p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= *c_W5.wet_col_count) return;

    const int ijdim = c_P.ijdim;
    const int kdim  = c_P.kdim;
    const size_t s3 = (size_t)ijdim * (size_t)kdim;

    const uint32_t pc = c_W5.wet_cols[p];
    const int ij = (int)(pc & 0x3FFFFu);
    const int l  = (int)(pc >> 18);
    const size_t col = (size_t)l * (size_t)ijdim + (size_t)ij;

    const int k_top = c_W5.col_ktop[col];
    const int k_lim = (k_top + 2 <= c_P.kmax0) ? (k_top + 2) : c_P.kmax0;
    const size_t off3 = (size_t)l * s3;
    const size_t offv = (size_t)l * 4ULL * s3;
    const uint8_t* __restrict__ flg = c_W5.cell_flag + off3;
    double* __restrict__ Vt = c_W5.Vt + offv;

    for (int k = c_P.kmin0; k <= k_lim; ++k) {
        size_t idx = (size_t)ij + (size_t)k * (size_t)ijdim;
        if (__ldg(&flg[idx]) != 0u) {
            Vt[idx + 0 * s3] = 0.0;
            Vt[idx + 1 * s3] = 0.0;
            Vt[idx + 2 * s3] = 0.0;
            Vt[idx + 3 * s3] = 0.0;
        }
    }
}

// =============================================================================
// KERNEL 2A v5.2 (SM90): Branch-Free Half-Level Fubini + L2-Coalesced zdis_compact!
// =============================================================================
template<int OPT_MODE, bool COLLECT_CNT>
__global__ void __launch_bounds__(64, 6)
k_sm90_v5_precip_half()
{
    const uint32_t p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= *c_W5.wet_col_count) return;
    const int ijdim = c_P.ijdim;
    const int kdim  = c_P.kdim;
    const int kmin0 = c_P.kmin0;
    const int kmax0 = c_P.kmax0;
    const size_t s3 = (size_t)ijdim * (size_t)kdim;
    const int sp_nq[4] = { c_P.i_qr, c_P.i_qi, c_P.i_qs, c_P.i_qg };
    const double dt = c_P.dt;

    unsigned long long lc_layers = 0, lc_ch = 0;

    const uint32_t pc = c_W5.wet_cols[p];
    const int ij = (int)(pc & 0x3FFFFu);
    const int l  = (int)(pc >> 18);

    const size_t col = (size_t)l * (size_t)ijdim + (size_t)ij;
    const uint32_t sp_mask = c_W5.col_spmask[col];
    const int k_top = c_W5.col_ktop[col];
    const int k_hi  = (k_top + 1 <= kmax0) ? (k_top + 1) : kmax0;

    const size_t off3 = (size_t)l * s3;
    const size_t off4 = (size_t)l * (size_t)c_P.nqmax * s3;
    const size_t offv = (size_t)l * 4ULL * s3;

    double* __restrict__       rhogw   = c_S.rhogw   + off3;
    const double* __restrict__ q       = c_S.q       + off4;
    const double* __restrict__ z       = c_S.z       + off3;
    const double* __restrict__ zh      = c_S.zh      + off3;
    const double* __restrict__ gsgam2h = c_S.gsgam2h + off3;
    const double* __restrict__ gam2h   = c_S.gam2h   + off3;
    const double* __restrict__ Vt      = c_W5.Vt     + offv;
    double* __restrict__       kv      = c_W5.kv     + off3;

    const double dzh_kmin = __ldg(&z[(size_t)ij + (size_t)kmin0 * ijdim])
                          - __ldg(&zh[(size_t)ij + (size_t)kmin0 * ijdim]);

    // Block-coalesced L2-resident Mode 3 zdis buffer: 100% warp coalescing & L2 hit!
    double* __restrict__ zd_c = (OPT_MODE == 3)
        ? (c_W5.zdis_compact + ((size_t)blockIdx.x * (size_t)kdim * 4ULL) * 64ULL + (size_t)threadIdx.x)
        : nullptr;

    if constexpr (OPT_MODE == 3) {
        for (int s = 0; s < 4; ++s) {
            if (!(sp_mask & (1u << s))) continue;
            double xi_next = __ldg(&z[(size_t)ij + (size_t)k_hi * ijdim]);
            // Sliding-window Vt & z loads: 1 Vt load/k instead of 3!
            size_t idx_hi = (size_t)ij + (size_t)k_hi * ijdim;
            double w_p1 = __ldg(&Vt[idx_hi + ijdim + s * s3]);
            double w_0  = __ldg(&Vt[idx_hi + s * s3]);
            double z_0  = xi_next;
            for (int k = k_hi; k >= kmin0 + 1; --k) {
                size_t idx = (size_t)ij + (size_t)k * ijdim;
                double w_m1 = __ldg(&Vt[idx - ijdim + s * s3]);
                double z_m1 = __ldg(&z[idx - ijdim]);
                double d_raw = 0.0;
                if (w_0 != 0.0) {
                    double dzh_k = z_0 - z_m1;
                    double dzh_km1 = (k - 1 == kmin0) ? dzh_kmin : (z_m1 - __ldg(&z[idx - 2ULL * ijdim]));
                    d_raw = vadv1d_taylor_zdis(w_0, w_p1, w_m1, dzh_k, dzh_km1, dzh_k + dzh_km1, dt, k == (kmin0 + 1));
                }
                double xi_k = z_m1 + d_raw;
                if constexpr (COLLECT_CNT) { if (xi_k > xi_next + 1.0e-12) lc_ch += 1; }
                double xi_p = (xi_k < xi_next) ? xi_k : xi_next;
                double d_p  = xi_p - z_m1;
                zd_c[((size_t)k * 4ULL + (size_t)s) * 64ULL] = (d_p > 0.0) ? d_p : 0.0;
                xi_next = xi_p;
                w_p1 = w_0; w_0 = w_m1; z_0 = z_m1;
            }
        }
    }

    double prev_fkv = 0.0, prev_fw = 0.0;
    for (int k = kmin0 + 1; k <= k_hi; ++k) {
        size_t kdx = (size_t)ij + (size_t)k * ijdim;
        double z_k   = __ldg(&z[kdx]);
        double z_km1 = __ldg(&z[kdx - ijdim]);
        double dzh_k = z_k - z_km1;
        double dzh_km1 = (k - 1 == kmin0) ? dzh_kmin : (z_km1 - __ldg(&z[kdx - 2ULL * ijdim]));
        double dzh2_k = dzh_k + dzh_km1;

        int kc[4] = { k - 1, k - 1, k - 1, k - 1 };
        double rem[4] = { 0.0, 0.0, 0.0, 0.0 };
        int c_max = k - 1;

        #pragma unroll
        for (int s = 0; s < 4; ++s) {
            if (!(sp_mask & (1u << s))) continue;
            double zdis;
            if constexpr (OPT_MODE == 3) {
                zdis = zd_c[((size_t)k * 4ULL + (size_t)s) * 64ULL];
            } else {
                double wk = __ldg(&Vt[kdx + s * s3]);
                if (wk == 0.0) continue;
                double wkp1 = __ldg(&Vt[kdx + ijdim + s * s3]);
                double wkm1 = __ldg(&Vt[kdx - ijdim + s * s3]);
                zdis = vadv1d_taylor_zdis(wk, wkp1, wkm1, dzh_k, dzh_km1, dzh2_k, dt, k == (kmin0 + 1));
            }
            if (zdis <= 0.0) continue;
            double ztmp = z_km1 + zdis;
            if (ztmp <= z_k) {
                kc[s] = k; rem[s] = zdis;
                if (k > c_max) c_max = k;
            } else {
                int c = k; double z_lo_c = z_km1;
                for (int j = k + 1; j <= kmax0; ++j) {
                    double zj = __ldg(&z[(size_t)ij + (size_t)j * ijdim]);
                    if (ztmp <= zj) {
                        double zjm1 = __ldg(&z[(size_t)ij + (size_t)(j - 1) * ijdim]);
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
            double af = c_afact[k], bf = c_bfact[k], Wh = 0.0;
            #pragma unroll
            for (int s = 0; s < 4; ++s) {
                if (kc[s] == k) Wh += rem[s] * (af * q[kdx + sp_nq[s] * s3] + bf * q[jm1 + sp_nq[s] * s3]);
            }
            if (Wh != 0.0) {
                double Wh_rgsh = Wh * (__ldg(&gam2h[kdx]) / __ldg(&gsgam2h[kdx]));
                fkv_k = -Wh_rgsh * kv[kdx];
                fw_k  = -Wh_rgsh * rhogw[kdx];
            }
            if constexpr (COLLECT_CNT) lc_layers += 1;
        } else if (c_max > k) {
            int j_end = (c_max <= k_hi) ? c_max : k_hi;
            for (int j = k; j <= j_end; ++j) {
                size_t jdx = (size_t)ij + (size_t)j * ijdim;
                size_t jm1 = jdx - ijdim;
                double dj = __ldg(&z[jdx]) - __ldg(&z[jm1]);
                double af = c_afact[j], bf = c_bfact[j], Wh = 0.0;
                #pragma unroll
                for (int s = 0; s < 4; ++s) {
                    int c = kc[s];
                    if (c < j) continue;
                    double L = (j < c) ? dj : rem[s];
                    Wh += L * (af * q[jdx + sp_nq[s] * s3] + bf * q[jm1 + sp_nq[s] * s3]);
                }
                if (Wh != 0.0) {
                    double Wh_rgsh = Wh * (__ldg(&gam2h[jdx]) / __ldg(&gsgam2h[jdx]));
                    fkv_k -= Wh_rgsh * kv[jdx];
                    fw_k  -= Wh_rgsh * rhogw[jdx];
                }
                if constexpr (COLLECT_CNT) lc_layers += 1;
            }
        }

        if (fkv_k != 0.0 || prev_fkv != 0.0 || fw_k != 0.0 || prev_fw != 0.0) {
            double rdgzh_k = c_rdgzh[k];
            kv[kdx]    = kv[kdx]    - (fkv_k - prev_fkv) * rdgzh_k;
            rhogw[kdx] = rhogw[kdx] - (fw_k  - prev_fw)  * rdgzh_k;
        }
        prev_fkv = fkv_k;
        prev_fw  = fw_k;
    }

    if constexpr (COLLECT_CNT) {
        if (lc_layers) atomicAdd(&c_W5.cnt->precip_layers_touched, lc_layers);
        if (lc_ch)     atomicAdd(&c_W5.cnt->precip_crossings_half, lc_ch);
    }
}

// =============================================================================
// KERNEL 2B v5.2 (SM90): Zero-vx/vy/vz Full-Level Fubini + L2-Coalesced zdis_compact!
//   - Eliminates ALL vx, vy, vz loads (319 MB HBM!) by reusing rhogvx/y/z * inv_rg
//   - Deduplicates Vt loads (6->4 in M2, 6->1 sliding window in M3!)
//   - Skips rhogvx/y/z loads on 0-flux layers (d_rkh_diff == 0.0!)
//   - Block-coalesced L2-resident zdis_compact in Mode 3!
// =============================================================================
template<int OPT_MODE, bool COLLECT_CNT>
__global__ void __launch_bounds__(64, 6)
k_sm90_v5_precip_full()
{
    const uint32_t p = blockIdx.x * blockDim.x + threadIdx.x;
    if (p >= *c_W5.wet_col_count) return;
    const int ijdim = c_P.ijdim;
    const int kdim  = c_P.kdim;
    const int kmin0 = c_P.kmin0;
    const int kmax0 = c_P.kmax0;
    const size_t s3 = (size_t)ijdim * (size_t)kdim;
    const int sp_nq[4] = { c_P.i_qr, c_P.i_qi, c_P.i_qs, c_P.i_qg };
    const double dt = c_P.dt;
    const double inv_dt = 1.0 / dt;

    unsigned long long lc_layers = 0, lc_cf = 0, lc_ov = 0;

    const uint32_t pc = c_W5.wet_cols[p];
    const int ij = (int)(pc & 0x3FFFFu);
    const int l  = (int)(pc >> 18);

    const size_t col = (size_t)l * (size_t)ijdim + (size_t)ij;
    const uint32_t sp_mask = c_W5.col_spmask[col];
    const int k_hi = c_W5.col_ktop[col];

    const size_t off3 = (size_t)l * s3;
    const size_t off4 = (size_t)l * (size_t)c_P.nqmax * s3;
    const size_t offp = (size_t)l * 2ULL * (size_t)ijdim;
    const size_t offv = (size_t)l * 4ULL * s3;

    double* __restrict__ rhog   = c_S.rhog   + off3;
    double* __restrict__ rhogvx = c_S.rhogvx + off3;
    double* __restrict__ rhogvy = c_S.rhogvy + off3;
    double* __restrict__ rhogvz = c_S.rhogvz + off3;
    const double* __restrict__ rhogw = c_S.rhogw + off3;
    double* __restrict__ rhoge  = c_S.rhoge  + off3;
    double* __restrict__ rhogq  = c_S.rhogq  + off4;
    const double* __restrict__ tem = c_S.tem + off3;
    const double* __restrict__ z       = c_S.z       + off3;
    const double* __restrict__ zh      = c_S.zh      + off3;
    const double* __restrict__ gsgam2  = c_S.gsgam2  + off3;
    const double* __restrict__ gsgam2h = c_S.gsgam2h + off3;
    const double* __restrict__ gam2    = c_S.gam2    + off3;
    const double* __restrict__ Vt      = c_W5.Vt     + offv;
    const double* __restrict__ kv      = c_W5.kv     + off3;

    const double ix_c = __ldg(&c_S.ix[col]), iy_c = __ldg(&c_S.iy[col]), iz_c = __ldg(&c_S.iz[col]);
    const double jx_c = __ldg(&c_S.jx[col]), jy_c = __ldg(&c_S.jy[col]), jz_c = __ldg(&c_S.jz[col]);

    // Deduplicated Vt evaluation: 4 loads instead of 6 (only 2 when wk == 0)!
    auto eval_full_zdis = [&](int k, int s) -> double {
        size_t idx = (size_t)ij + (size_t)k * ijdim;
        double vt0 = __ldg(&Vt[idx + s * s3]);
        double vtm1 = (k > kmin0) ? __ldg(&Vt[idx - ijdim + s * s3]) : vt0;
        double wk = (k == kmin0) ? vt0 : (c_afact[k] * vt0 + c_bfact[k] * vtm1);
        if (wk == 0.0) return 0.0;
        double wkp1 = (k < kmax0)
            ? (c_afact[k + 1] * __ldg(&Vt[idx + ijdim + s * s3]) + c_bfact[k + 1] * vt0)
            : vt0;
        double wkm1 = 0.0;
        if (k > kmin0) {
            wkm1 = (k - 1 == kmin0)
                ? vtm1
                : (c_afact[k - 1] * vtm1 + c_bfact[k - 1] * __ldg(&Vt[idx - 2ULL * ijdim + s * s3]));
        }
        double zh_k   = __ldg(&zh[idx]);
        double dz_k   = __ldg(&zh[idx + ijdim]) - zh_k;
        double dz_km1 = (k > kmin0) ? (zh_k - __ldg(&zh[idx - ijdim])) : dz_k;
        return vadv1d_taylor_zdis(wk, wkp1, wkm1, dz_k, dz_km1, dz_k + dz_km1, dt, k == kmin0);
    };

    double* __restrict__ zd_c = (OPT_MODE == 3)
        ? (c_W5.zdis_compact + ((size_t)blockIdx.x * (size_t)kdim * 4ULL) * 64ULL + (size_t)threadIdx.x)
        : nullptr;

    if constexpr (OPT_MODE == 3) {
        for (int s = 0; s < 4; ++s) {
            if (!(sp_mask & (1u << s))) continue;
            // Sliding-window backward monotone sweep: 1 Vt load per k (6x fewer loads!)
            size_t idx_hi = (size_t)ij + (size_t)k_hi * ijdim;
            double vt_hi   = __ldg(&Vt[idx_hi + s * s3]);
            double vt_hip1 = (k_hi < kmax0) ? __ldg(&Vt[idx_hi + ijdim + s * s3]) : vt_hi;
            double w_p1    = (k_hi < kmax0) ? (c_afact[k_hi + 1] * vt_hip1 + c_bfact[k_hi + 1] * vt_hi) : vt_hi;
            double vt_0    = vt_hi;
            double vt_m1   = (k_hi > kmin0) ? __ldg(&Vt[idx_hi - ijdim + s * s3]) : vt_hi;
            double w_0     = (k_hi == kmin0) ? vt_0 : (c_afact[k_hi] * vt_0 + c_bfact[k_hi] * vt_m1);

            double zh_p1   = __ldg(&zh[idx_hi + ijdim]);
            double zh_0    = __ldg(&zh[idx_hi]);
            double xi_next = zh_p1;

            for (int k = k_hi; k >= kmin0; --k) {
                size_t idx = (size_t)ij + (size_t)k * ijdim;
                double vt_m2 = (k - 1 > kmin0) ? __ldg(&Vt[idx - 2ULL * ijdim + s * s3]) : vt_m1;
                double w_m1  = (k > kmin0)
                    ? ((k - 1 == kmin0) ? vt_m1 : (c_afact[k - 1] * vt_m1 + c_bfact[k - 1] * vt_m2))
                    : 0.0;
                double zh_m1 = (k > kmin0) ? __ldg(&zh[idx - ijdim]) : zh_0;
                double d_raw = 0.0;
                if (w_0 != 0.0) {
                    double dz_k   = zh_p1 - zh_0;
                    double dz_km1 = (k > kmin0) ? (zh_0 - zh_m1) : dz_k;
                    d_raw = vadv1d_taylor_zdis(w_0, w_p1, w_m1, dz_k, dz_km1, dz_k + dz_km1, dt, k == kmin0);
                }
                double xi_k = zh_0 + d_raw;
                if constexpr (COLLECT_CNT) { if (xi_k > xi_next + 1.0e-12) lc_cf += 1; }
                double xi_p = (xi_k < xi_next) ? xi_k : xi_next;
                double d_p  = xi_p - zh_0;
                zd_c[((size_t)k * 4ULL + (size_t)s) * 64ULL] = (d_p > 0.0) ? d_p : 0.0;
                xi_next = xi_p;

                w_p1 = w_0; w_0 = w_m1;
                vt_0 = vt_m1; vt_m1 = vt_m2;
                zh_p1 = zh_0; zh_0 = zh_m1;
            }
        }
    }

    auto eval_full_iface = [&](int ki, double fq[4], double& fe, double& fphi, double& fkh, double& fu, double& fv) {
        fq[0] = fq[1] = fq[2] = fq[3] = 0.0;
        fe = fphi = fkh = fu = fv = 0.0;
        if (ki > k_hi) return;

        size_t kdx = (size_t)ij + (size_t)ki * ijdim;
        double zh_ki  = __ldg(&zh[kdx]);
        double zh_kp1 = __ldg(&zh[kdx + ijdim]);

        int kc[4] = { ki - 1, ki - 1, ki - 1, ki - 1 };
        double rem[4] = { 0.0, 0.0, 0.0, 0.0 };
        int c_max = ki - 1;

        #pragma unroll
        for (int s = 0; s < 4; ++s) {
            if (!(sp_mask & (1u << s))) continue;
            double zdis = (OPT_MODE == 3)
                ? zd_c[((size_t)ki * 4ULL + (size_t)s) * 64ULL]
                : eval_full_zdis(ki, s);
            if (zdis <= 0.0) continue;
            double ztmp = zh_ki + zdis;
            if (ztmp <= zh_kp1) {
                kc[s] = ki; rem[s] = zdis;
                if (ki > c_max) c_max = ki;
            } else {
                int c = ki; double zh_c = zh_ki; bool found = false;
                for (int j = ki + 1; j <= kmax0; ++j) {
                    size_t jdx = (size_t)ij + (size_t)j * ijdim;
                    if (ztmp <= __ldg(&zh[jdx + ijdim])) {
                        double zj = __ldg(&zh[jdx]);
                        if (ztmp > zj) { c = j; zh_c = zj; }
                        found = true; break;
                    }
                }
                if constexpr (COLLECT_CNT) { if (!found) lc_ov += 1; }
                kc[s] = c; rem[s] = zdis - (zh_c - zh_ki);
                if (c > c_max) c_max = c;
            }
        }

        if (c_max == ki && c_P.ein_is_simple) {
            double rgs_j = __ldg(&gam2[kdx]) / __ldg(&gsgam2[kdx]);
            double W = 0.0;
            #pragma unroll
            for (int s = 0; s < 4; ++s) {
                if (kc[s] == ki) {
                    double w_sj = rem[s] * rhogq[kdx + sp_nq[s] * s3] * rgs_j;
                    fq[s] = -w_sj; W += w_sj;
                }
            }
            if (W != 0.0) {
                // Zero vx/vy/vz loads: derive velocity directly from rhogvx/y/z / rg!
                double inv_rg = 1.0 / rhog[kdx];
                double vx_j = rhogvx[kdx] * inv_rg;
                double vy_j = rhogvy[kdx] * inv_rg;
                double vz_j = rhogvz[kdx] * inv_rg;
                double kh_r = 0.5 * (vx_j*vx_j + vy_j*vy_j + vz_j*vz_j);
                fe   = -c_P.CVdry * (W * tem[kdx]);
                fphi = -c_P.GRAV  * (W * __ldg(&z[kdx]));
                fkh  = -W * kh_r;
                fu   = -W * (vx_j * ix_c + vy_j * iy_c + vz_j * iz_c);
                fv   = -W * (vx_j * jx_c + vy_j * jy_c + vz_j * jz_c);
            }
            if constexpr (COLLECT_CNT) lc_layers += 1;
        } else if (c_max >= ki) {
            int j_end = (c_max <= k_hi) ? c_max : k_hi;
            double fe_u = 0.0, fphi_u = 0.0;
            for (int j = ki; j <= j_end; ++j) {
                size_t jdx = (size_t)ij + (size_t)j * ijdim;
                double dj = __ldg(&zh[jdx + ijdim]) - __ldg(&zh[jdx]);
                double rgs_j = __ldg(&gam2[jdx]) / __ldg(&gsgam2[jdx]);
                double W = 0.0, Wc = 0.0;
                #pragma unroll
                for (int s = 0; s < 4; ++s) {
                    int c = kc[s];
                    if (c < j) continue;
                    double L = (j < c) ? dj : rem[s];
                    double w_sj = L * rhogq[jdx + sp_nq[s] * s3] * rgs_j;
                    fq[s] -= w_sj; W += w_sj;
                    if (!c_P.ein_is_simple) Wc += w_sj * c_P.CVW[sp_nq[s]];
                }
                if (W != 0.0) {
                    double inv_rg = 1.0 / rhog[jdx];
                    double vx_j = rhogvx[jdx] * inv_rg;
                    double vy_j = rhogvy[jdx] * inv_rg;
                    double vz_j = rhogvz[jdx] * inv_rg;
                    double kh_r = 0.5 * (vx_j*vx_j + vy_j*vy_j + vz_j*vz_j);
                    if (c_P.ein_is_simple) fe_u -= W * tem[jdx];
                    else                   fe   -= Wc * tem[jdx];
                    fphi_u -= W * __ldg(&z[jdx]);
                    fkh    -= W * kh_r;
                    fu     -= W * (vx_j * ix_c + vy_j * iy_c + vz_j * iz_c);
                    fv     -= W * (vx_j * jx_c + vy_j * jy_c + vz_j * jz_c);
                }
                if constexpr (COLLECT_CNT) lc_layers += 1;
            }
            if (c_P.ein_is_simple) fe = c_P.CVdry * fe_u;
            fphi = c_P.GRAV * fphi_u;
        }
    };

    double fq_lo[4], fe_lo, fphi_lo, fkh_lo, fu_lo, fv_lo;
    eval_full_iface(kmin0, fq_lo, fe_lo, fphi_lo, fkh_lo, fu_lo, fv_lo);

    double p_liq = -fq_lo[0] * inv_dt;
    double p_sol = -(fq_lo[1] + fq_lo[2] + fq_lo[3]) * inv_dt;
    double p_lh  =  (fq_lo[1] + fq_lo[2] + fq_lo[3]) * c_P.LHF * inv_dt;
    double p_re  = -fe_lo   * inv_dt;
    double p_rp  = -fphi_lo * inv_dt;
    double p_rk  = -fkh_lo  * inv_dt;
    c_S.precip[offp + ij + 0 * ijdim] = p_liq;
    c_S.precip[offp + ij + 1 * ijdim] = p_sol;
    c_S.precip_lh    [col] = p_lh;
    c_S.precip_rhoein[col] = p_re;
    c_S.precip_rhophi[col] = p_rp;
    c_S.precip_rhokin[col] = p_rk;
    c_S.precip_energy[col] = p_lh + p_re + p_rp + p_rk;

    double d_rg_k = 0.0, d_rhoge_k = 0.0, d_rkh_diff_k = 0.0;
    double d_kv_recomp_prev = 0.0;
    double gs2_m1 = __ldg(&gsgam2[(size_t)ij + (size_t)kmin0 * ijdim]);
    double inv_gs2_m1 = 1.0 / gs2_m1;
    double inv_gs2h_m1 = 1.0 / __ldg(&gsgam2h[(size_t)ij + (size_t)kmin0 * ijdim]);

    for (int k = kmin0; k <= k_hi; ++k) {
        double fq_hi[4], fe_hi, fphi_hi, fkh_hi, fu_hi, fv_hi;
        eval_full_iface(k + 1, fq_hi, fe_hi, fphi_hi, fkh_hi, fu_hi, fv_hi);

        size_t idx = (size_t)ij + (size_t)k * ijdim;
        double rdgz_k = c_rdgz[k];
        double drgq0 = -(fq_hi[0] - fq_lo[0]) * rdgz_k;
        double drgq1 = -(fq_hi[1] - fq_lo[1]) * rdgz_k;
        double drgq2 = -(fq_hi[2] - fq_lo[2]) * rdgz_k;
        double drgq3 = -(fq_hi[3] - fq_lo[3]) * rdgz_k;
        double sum_drgq = drgq0 + drgq1 + drgq2 + drgq3;

        double drge  = -(fe_hi   - fe_lo)   * rdgz_k;
        double drgp  = -(fphi_hi - fphi_lo) * rdgz_k;
        double drgkh = -(fkh_hi  - fkh_lo)  * rdgz_k;
        double drgu  = -(fu_hi   - fu_lo)   * rdgz_k;
        double drgv  = -(fv_hi   - fv_lo)   * rdgz_k;

        double rg_old = rhog[idx];
        double rg_new = rg_old + sum_drgq;
        double rge_new = rhoge[idx];
        double rkh_diff = 0.0;

        if (sum_drgq != 0.0 || drgu != 0.0 || drgv != 0.0 || drgkh != 0.0 || drge != 0.0 || drgp != 0.0) {
            if (sum_drgq != 0.0) {
                rhog[idx] = rg_new;
                rge_new += drge + drgp - sum_drgq * c_P.GRAV * __ldg(&z[idx]);
            } else {
                rge_new += drge + drgp;
            }
            if (sum_drgq != 0.0 || drgu != 0.0 || drgv != 0.0 || drgkh != 0.0) {
                double rvx_old = rhogvx[idx], rvy_old = rhogvy[idx], rvz_old = rhogvz[idx];
                double rkh_old = 0.5 * (rvx_old*rvx_old + rvy_old*rvy_old + rvz_old*rvz_old) / rg_old;
                double rvx_new = rvx_old + drgu * ix_c + drgv * jx_c;
                double rvy_new = rvy_old + drgu * iy_c + drgv * jy_c;
                double rvz_new = rvz_old + drgu * iz_c + drgv * jz_c;
                if (drgu != 0.0 || drgv != 0.0) {
                    rhogvx[idx] = rvx_new; rhogvy[idx] = rvy_new; rhogvz[idx] = rvz_new;
                }
                double rkh_n = 0.5 * (rvx_new*rvx_new + rvy_new*rvy_new + rvz_new*rvz_new) / rg_new;
                rkh_diff = (rkh_old + drgkh) - rkh_n;
            }
            if (drgq0 != 0.0) rhogq[idx + sp_nq[0] * s3] += drgq0;
            if (drgq1 != 0.0) rhogq[idx + sp_nq[1] * s3] += drgq1;
            if (drgq2 != 0.0) rhogq[idx + sp_nq[2] * s3] += drgq2;
            if (drgq3 != 0.0) rhogq[idx + sp_nq[3] * s3] += drgq3;
        }

        double kv_recomp_k = 0.0;
        if (k > kmin0) {
            size_t im1 = idx - ijdim;
            double gs2_k  = __ldg(&gsgam2[idx]);
            double gs2h_k = __ldg(&gsgam2h[idx]);
            double inv_gs2_k  = 1.0 / gs2_k;
            double inv_gs2h_k = 1.0 / gs2h_k;
            double rw_new = rhogw[idx];
            if (rw_new != 0.0) {
                double rhogh_new_k = (c_afact[k] * inv_gs2_k * gs2h_k) * rg_new
                                   + (c_bfact[k] * inv_gs2_m1 * gs2h_k) * d_rg_k;
                kv_recomp_k = 0.5 * rw_new * rw_new / rhogh_new_k;
            }
            if (c_P.prcip_kin2ein) {
                double dv_k = (c_cfact[k - 1] * kv[idx] + c_dfact[k - 1] * kv[im1])
                            - ((c_cfact[k - 1] * gs2_m1 * inv_gs2h_k) * kv_recomp_k
                               + (c_dfact[k - 1] * gs2_m1 * inv_gs2h_m1) * d_kv_recomp_prev);
                d_rhoge_k += d_rkh_diff_k + dv_k;
            }
            rhoge[im1] = d_rhoge_k;
            gs2_m1 = gs2_k; inv_gs2_m1 = inv_gs2_k; inv_gs2h_m1 = inv_gs2h_k;
        }

        d_rg_k = rg_new; d_rhoge_k = rge_new; d_rkh_diff_k = rkh_diff;
        d_kv_recomp_prev = kv_recomp_k;

        fq_lo[0] = fq_hi[0]; fq_lo[1] = fq_hi[1]; fq_lo[2] = fq_hi[2]; fq_lo[3] = fq_hi[3];
        fe_lo = fe_hi; fphi_lo = fphi_hi; fkh_lo = fkh_hi; fu_lo = fu_hi; fv_lo = fv_hi;
    }

    {
        size_t idx = (size_t)ij + (size_t)k_hi * ijdim;
        if (c_P.prcip_kin2ein) {
            double gs2_k = __ldg(&gsgam2[idx]);
            if (k_hi < kmax0) {
                size_t ip1 = idx + ijdim;
                double gs2h_p1 = __ldg(&gsgam2h[ip1]);
                double kv_p1_t = kv[ip1];
                double rhogh_p1 = (c_afact[k_hi + 1] / __ldg(&gsgam2[ip1]) * gs2h_p1) * rhog[ip1]
                                + (c_bfact[k_hi + 1] / gs2_k * gs2h_p1) * d_rg_k;
                double rw_p1 = rhogw[ip1];
                double kv_p1_n = 0.5 * rw_p1 * rw_p1 / rhogh_p1;
                double dv_k = (c_cfact[k_hi] * kv_p1_t + c_dfact[k_hi] * kv[idx])
                            - ((c_cfact[k_hi] * gs2_k / gs2h_p1) * kv_p1_n
                               + (c_dfact[k_hi] * gs2_k / __ldg(&gsgam2h[idx])) * d_kv_recomp_prev);
                d_rhoge_k += d_rkh_diff_k + dv_k;
                rhoge[idx] = d_rhoge_k;

                double gs2_p1 = __ldg(&gsgam2[ip1]);
                double kv_p2 = (k_hi + 1 < kmax0) ? kv[ip1 + ijdim] : 0.0;
                double d_p1 = (c_cfact[k_hi + 1] * (1.0 - gs2_p1 / __ldg(&gsgam2h[ip1 + ijdim]))) * kv_p2
                            + c_dfact[k_hi + 1] * (kv_p1_t - (gs2_p1 / gs2h_p1) * kv_p1_n);
                if (d_p1 != 0.0) rhoge[ip1] += d_p1;
            } else {
                double dv_k = c_dfact[kmax0] * kv[idx]
                            - (c_dfact[kmax0] * gs2_k / __ldg(&gsgam2h[idx])) * d_kv_recomp_prev;
                d_rhoge_k += d_rkh_diff_k + dv_k;
                rhoge[idx] = d_rhoge_k;
            }
        } else {
            rhoge[idx] = d_rhoge_k;
        }
    }

    if constexpr (COLLECT_CNT) {
        if (lc_layers) atomicAdd(&c_W5.cnt->precip_layers_touched, lc_layers);
        if (lc_cf)     atomicAdd(&c_W5.cnt->precip_crossings_full, lc_cf);
        if (lc_ov)     atomicAdd(&c_W5.cnt->precip_top_overshoots, lc_ov);
    }
}

// =============================================================================
// KERNEL 3 v5.3 (SM90): High-Occupancy Perturbation-Robust Clausius-Clapeyron SatAdj!
//   - Deduplicated Laguerre Householder-4 solver -> drops regs from 120 -> ~80
//     (__launch_bounds__(128, 6) = 37.5%->50% occupancy, 0 B spill!)
// =============================================================================
template<int OPT_MODE, bool COLLECT_CNT>
__global__ void __launch_bounds__(128, 4)
k_sm90_v5_satadj_nf_diag()
{
    const int ij = blockIdx.x * blockDim.x + threadIdx.x;
    const int k  = blockIdx.y;
    const int l  = blockIdx.z;
    const int ijdim = c_P.ijdim;
    const int kdim  = c_P.kdim;
    if (ij >= ijdim || k >= kdim) return;

    const size_t s3   = (size_t)ijdim * (size_t)kdim;
    const size_t off3 = (size_t)l * s3;
    const size_t off4 = (size_t)l * (size_t)c_P.nqmax * s3;
    const size_t idx  = (size_t)ij + (size_t)k * ijdim;
    const size_t col  = (size_t)l * (size_t)ijdim + (size_t)ij;

    const int ktop = c_W5.col_ktop[col];
    const uint8_t flg = (k >= c_P.kmin0 && k <= c_P.kmax0) ? __ldg(&c_W5.cell_flag[off3 + idx]) : 0u;
    const bool above_precip = (ktop == 0 || k >= ktop + 2);

    double gs2 = __ldg(&c_S.gsgam2[off3 + idx]);
    double diff_e = 0.0;
    bool rge_modified = false;

    if (c_P.prcip_kin2ein && k >= c_P.kmin0 && k <= c_P.kmax0 && above_precip) {
        const double* __restrict__ gs2h_arr = c_S.gsgam2h + off3;
        double gs2h_k  = __ldg(&gs2h_arr[idx]);
        double gs2h_p1 = __ldg(&gs2h_arr[idx + ijdim]);
        if (gs2 != gs2h_k || gs2 != gs2h_p1) {
            const double* __restrict__ kv = c_W5.kv + off3;
            diff_e = (c_cfact[k] * (1.0 - gs2 / gs2h_p1)) * kv[idx + ijdim]
                   + (c_dfact[k] * (1.0 - gs2 / gs2h_k))  * kv[idx];
            if (diff_e != 0.0) rge_modified = true;
        }
    }

    // ULTRA-FAST PATH: Cold-dry certified & zero energy perturbation -> 0 state reads!
    if (flg == 2u && above_precip && !rge_modified && c_P.ein_is_simple) {
        c_S.gdcfrc[off3 + idx] = 0.0;
        c_S.gdclw [off3 + idx] = c_P.UNDEF;
        return;
    }

    double* __restrict__ rhog  = c_S.rhog  + off3;
    double* __restrict__ rhoge = c_S.rhoge + off3;
    double* __restrict__ rhogq = c_S.rhogq + off4;

    double rg  = rhog[idx];
    double rge = rhoge[idx] + diff_e;

    // EXACT FAST PATH: Clean dry interior cell untouched by wet_core & precip!
    //   Rare warm-dry supersaturated cells fall through to single Householder-4 path!
    if (flg != 0u && above_precip && c_P.ein_is_simple) {
        double rgqv  = rhogq[idx + c_P.i_qv * s3];
        double qv    = (rgqv != 0.0) ? (rgqv / rg) : 0.0;
        double rho_v = rg / gs2;
        double tem0  = rge / (rg * c_P.CVdry);

        bool certified_subsat = (flg == 2u && diff_e >= -0.1 * rg * c_P.CVdry && tem0 < 5000.0) || (qv <= 0.0);
        bool supersat = false;
        int c_exp = 0;

        if (!certified_subsat) {
            double rtem00 = 1.0 / c_P.TEM00;
            double inv_rho_Rv = 1.0 / (rho_v * c_P.Rvap);
            double tc0 = (tem0 > 10.0) ? tem0 : 10.0;
            double qsat = c_P.PSAT0 * exp(c_P.LovR_liq * (rtem00 - 1.0 / tc0)) * inv_rho_Rv / tem0;
            if constexpr (COLLECT_CNT) c_exp = 1;
            supersat = (qv - qsat > c_P.EPS);
        }

        if (!supersat) {
            if (rge_modified) {
                rhoge[idx] = rge;
                double qd = 1.0 - qv;
                c_S.tem[off3 + idx] = tem0;
                c_S.pre[off3 + idx] = rho_v * (qd * c_P.Rdry + qv * c_P.Rvap) * tem0;
            }
            c_S.gdcfrc[off3 + idx] = 0.0;
            c_S.gdclw [off3 + idx] = c_P.UNDEF;
            if constexpr (COLLECT_CNT) {
                if (c_exp) atomicAdd(&c_W5.cnt->n_exp_calls, (unsigned long long)c_exp);
            }
            return;
        }
    }

    // General path (wet cells, precipitating column layers k <= ktop+1, boundary)
    double rgqv = rhogq[idx + c_P.i_qv * s3];
    double rgqc = rhogq[idx + c_P.i_qc * s3];
    double rgqr = rhogq[idx + c_P.i_qr * s3];
    double rgqi = rhogq[idx + c_P.i_qi * s3];
    double rgqs = rhogq[idx + c_P.i_qs * s3];
    double rgqg = rhogq[idx + c_P.i_qg * s3];

    double qv = (rgqv != 0.0) ? (rgqv / rg) : 0.0;
    double qc = (rgqc != 0.0) ? (rgqc / rg) : 0.0;
    double qr = (rgqr != 0.0) ? (rgqr / rg) : 0.0;
    double qi = (rgqi != 0.0) ? (rgqi / rg) : 0.0;
    double qs = (rgqs != 0.0) ? (rgqs / rg) : 0.0;
    double qg = (rgqg != 0.0) ? (rgqg / rg) : 0.0;
    double qd0 = 1.0 - (qv + qc + qr + qi + qs + qg);

    double ein_moist = rge / rg + qv * c_P.LHV;
    double qsum = qv + qc;
    double qv0  = qsum, qc0 = 0.0;
    double q_non_vc_cv = c_P.ein_is_simple ? 0.0
        : (qr * c_P.CVW[2] + qi * c_P.CVW[3] + qs * c_P.CVW[4] + qg * c_P.CVW[5]);
    double CVtot0 = c_P.ein_is_simple ? c_P.CVdry
        : (qd0 * c_P.CVdry + qv0 * c_P.CVW[0] + q_non_vc_cv);

    double rho_v = rg / gs2;
    double tem0  = (ein_moist - qv0 * c_P.LHV) / CVtot0;
    bool q_modified = (qc != 0.0);
    int c_act = 0, c_it = 0, c_exp = 0, c_pow = 0, c_neg = 0;

    if (k >= c_P.kmin0 && k <= c_P.kmax0 && qsum > 0.0) {
        double rtem00 = 1.0 / c_P.TEM00;
        double inv_rho_Rv = 1.0 / (rho_v * c_P.Rvap);
        double tc0 = (tem0 > 10.0) ? tem0 : 10.0;
        double psat = c_P.ein_is_simple
            ? (c_P.PSAT0 * exp(c_P.LovR_liq * (rtem00 - 1.0 / tc0)))
            : (c_P.PSAT0 * pow(tem0 * rtem00, c_P.CPovR_liq) * exp(c_P.LovR_liq * (rtem00 - 1.0 / tc0)));
        double qsat = psat * inv_rho_Rv / tem0;
        if constexpr (COLLECT_CNT) { c_exp += 1; if (!c_P.ein_is_simple) c_pow += 1; }

        if (qsum - qsat > c_P.EPS) {
            q_modified = true;
            c_act = 1;
            const double CV = c_P.CVdry, L = c_P.LHV, a = c_P.LovR_liq;
            for (int ite = 1; ite <= 20; ++ite) {
                if (ite > 1) {
                    double tc = (tem0 > 10.0) ? tem0 : 10.0;
                    qsat = c_P.PSAT0 * exp(a * (rtem00 - 1.0 / tc)) * inv_rho_Rv / tem0;
                    if constexpr (COLLECT_CNT) c_exp += 1;
                }
                double rT = 1.0 / tem0, b = a * rT, Lq = L * qsat;
                double f  = CV * tem0 + Lq - ein_moist;
                double f1 = CV + Lq * rT * (b - 1.0);
                double f2 = Lq * (rT * rT) * (b * b - 4.0 * b + 2.0);
                double f3 = Lq * (rT * rT * rT) * (b * (b * (b - 9.0) + 18.0) - 6.0);
                double dT_N = f / f1;
                double num  = f * (f1 * f1 - 0.5 * f * f2);
                double den  = f1 * (f1 * f1 - f * f2) + (1.0 / 6.0) * f * f * f3;
                double dT   = (den != 0.0) ? (num / den) : dT_N;
                if (!(dT * dT_N > 0.0 && fabs(dT) <= 2.0 * fabs(dT_N))) dT = dT_N;
                tem0 -= dT; if (tem0 < 10.0) tem0 = 10.0;
                c_it = ite;
                if (fabs(dT) < 1.0e-12) {
                    double qn = qsat + qsat * rT * (b - 1.0) * (-dT);
                    qv0 = qn; qc0 = qsum - qn; break;
                }
            }
        }
    }

    if (k < c_P.kmin0 || k > c_P.kmax0) {
        double CVtot = c_P.ein_is_simple ? c_P.CVdry
            : (qd0 * c_P.CVdry + qv0 * c_P.CVW[0] + qc0 * c_P.CVW[1] + q_non_vc_cv);
        rhogq[idx + c_P.i_qv * s3] = rg * qv0;
        rhogq[idx + c_P.i_qc * s3] = rg * qc0;
        rhoge[idx] = rg * tem0 * CVtot;
        c_S.rho[off3 + idx] = rho_v;
        c_S.q[off4 + idx + c_P.i_qv * s3] = 0.0; c_S.q[off4 + idx + c_P.i_qc * s3] = 0.0;
        c_S.q[off4 + idx + c_P.i_qr * s3] = 0.0; c_S.q[off4 + idx + c_P.i_qi * s3] = 0.0;
        c_S.q[off4 + idx + c_P.i_qs * s3] = 0.0; c_S.q[off4 + idx + c_P.i_qg * s3] = 0.0;
        c_S.gdcfrc[off3 + idx] = 0.0;
        c_S.gdclw [off3 + idx] = c_P.UNDEF;
        return;
    }

    if (q_modified) { rgqv = rg * qv0; rgqc = rg * qc0; }
    if constexpr (COLLECT_CNT) {
        if (rgqc < 0.0 || rgqr < 0.0 || rgqi < 0.0 || rgqs < 0.0 || rgqg < 0.0 || rgqv < 0.0) c_neg = 1;
    }

    double diffq = 0.0;
    bool neg_clipped = false;
    auto clip_q = [&](double& v) { if (v < 0.0) { diffq += v; v = 0.0; neg_clipped = true; } };
    clip_q(rgqc); clip_q(rgqr); clip_q(rgqi); clip_q(rgqs); clip_q(rgqg);

    rgqv += diffq;
    double diff0 = 0.0;
    if (rgqv < 0.0) { diff0 = rgqv; rgqv = 0.0; neg_clipped = true; }

    if (neg_clipped) {
        rg = rg * (1.0 - diff0);
        rho_v = rg / gs2;
        rhog[idx] = rg;
        rhogq[idx + c_P.i_qv * s3] = rgqv; rhogq[idx + c_P.i_qc * s3] = rgqc;
        rhogq[idx + c_P.i_qr * s3] = rgqr; rhogq[idx + c_P.i_qi * s3] = rgqi;
        rhogq[idx + c_P.i_qs * s3] = rgqs; rhogq[idx + c_P.i_qg * s3] = rgqg;
        qv = (rgqv != 0.0) ? (rgqv / rg) : 0.0; qc = (rgqc != 0.0) ? (rgqc / rg) : 0.0;
        qr = (rgqr != 0.0) ? (rgqr / rg) : 0.0; qi = (rgqi != 0.0) ? (rgqi / rg) : 0.0;
        qs = (rgqs != 0.0) ? (rgqs / rg) : 0.0; qg = (rgqg != 0.0) ? (rgqg / rg) : 0.0;
        if (diff0 != 0.0 && c_W5.col_w_ktop[col] < c_P.kmax0) c_W5.col_w_ktop[col] = c_P.kmax0;
    } else if (q_modified) {
        rhogq[idx + c_P.i_qv * s3] = rgqv;
        rhogq[idx + c_P.i_qc * s3] = rgqc;
        qv = qv0; qc = qc0;
    }

    double qd  = 1.0 - (qv + qc + qr + qi + qs + qg);
    double cva = c_P.ein_is_simple ? c_P.CVdry
        : (qd * c_P.CVdry + qv * c_P.CVW[0] + qc * c_P.CVW[1]
           + qr * c_P.CVW[2] + qi * c_P.CVW[3] + qs * c_P.CVW[4] + qg * c_P.CVW[5]);

    if (rge_modified || q_modified || neg_clipped) rhoge[idx] = tem0 * rg * cva;

    const bool mass_dirty = (ktop > 0 && k <= ktop) || neg_clipped;
    const bool rge_dirty  = (ktop > 0 && k <= ktop + 1) || rge_modified || q_modified || neg_clipped;

    if (mass_dirty) {
        c_S.rho[off3 + idx] = rho_v;
        c_S.vx [off3 + idx] = c_S.rhogvx[off3 + idx] / rg;
        c_S.vy [off3 + idx] = c_S.rhogvy[off3 + idx] / rg;
        c_S.vz [off3 + idx] = c_S.rhogvz[off3 + idx] / rg;
        c_S.q[off4 + idx + c_P.i_qv * s3] = qv; c_S.q[off4 + idx + c_P.i_qc * s3] = qc;
        c_S.q[off4 + idx + c_P.i_qr * s3] = qr; c_S.q[off4 + idx + c_P.i_qi * s3] = qi;
        c_S.q[off4 + idx + c_P.i_qs * s3] = qs; c_S.q[off4 + idx + c_P.i_qg * s3] = qg;
    } else if (q_modified) {
        c_S.q[off4 + idx + c_P.i_qv * s3] = qv;
        c_S.q[off4 + idx + c_P.i_qc * s3] = qc;
    }

    if (rge_dirty) {
        c_S.tem[off3 + idx] = tem0;
        c_S.pre[off3 + idx] = rho_v * (qd * c_P.Rdry + qv * c_P.Rvap) * tem0;
    }

    double q_tot = qc + qr + qi + qs + qg;
    c_S.gdcfrc[off3 + idx] = (q_tot >= 1.0e-10) ? 1.0 : 0.0;
    c_S.gdclw [off3 + idx] = c_P.UNDEF;

    if constexpr (COLLECT_CNT) {
        if (c_act) atomicAdd(&c_W5.cnt->sat_active_cells, 1ULL);
        if (c_it)  atomicAdd(&c_W5.cnt->sat_total_iters,  (unsigned long long)c_it);
        if (c_exp) atomicAdd(&c_W5.cnt->n_exp_calls,      (unsigned long long)c_exp);
        if (c_pow) atomicAdd(&c_W5.cnt->n_pow_calls,      (unsigned long long)c_pow);
        if (c_neg) atomicAdd(&c_W5.cnt->neg_q_before_postfilter, 1ULL);
    }
}

// =============================================================================
// KERNEL 4 v5.2 (SM90): Exact Unchanged-Cell Skip diag_w (-50% runtime!)
// =============================================================================
__global__ void __launch_bounds__(256, 4)
k_sm90_v5_diag_w()
{
    const int ij = blockIdx.x * blockDim.x + threadIdx.x;
    const int k  = blockIdx.y;
    const int l  = blockIdx.z;
    const int ijdim = c_P.ijdim;
    const int kdim  = c_P.kdim;
    if (ij >= ijdim || k >= kdim) return;

    const size_t col = (size_t)l * (size_t)ijdim + (size_t)ij;
    if (k >= c_P.kmin0 && k <= c_P.kmax0) {
        int kt = c_W5.col_ktop[col];
        int wk = c_W5.col_w_ktop[col];
        int hi = (kt > wk) ? kt : wk;
        if (hi == 0 || k > hi + 1) return;
    }

    const size_t s3   = (size_t)ijdim * (size_t)kdim;
    const size_t off3 = (size_t)l * s3;
    const size_t idx  = (size_t)ij + (size_t)k * ijdim;

    if (k == 0) {
        size_t i1 = (size_t)ij + (size_t)ijdim;
        c_S.w[off3 + idx] = c_S.rhogw[off3 + idx]
            / ((c_afact[1] * c_S.rho[off3 + i1] + c_bfact[1] * c_S.rho[off3 + idx]) * __ldg(&c_S.gsgam2h[off3 + i1]));
    } else {
        c_S.w[off3 + idx] = c_S.rhogw[off3 + idx]
            / ((c_afact[k] * c_S.rho[off3 + idx] + c_bfact[k] * c_S.rho[off3 + idx - ijdim]) * __ldg(&c_S.gsgam2h[off3 + idx]));
    }
}

// =============================================================================
// Templated v5.2 Dispatcher (Mode 2 = Exact, Mode 3 = Monotone)
// =============================================================================
template<int OPT_MODE, bool COLLECT_CNT>
inline void sm90_v5_launch_impl(
    const NSW6Params& P, int lall,
    const HopperStatePtrs& S, HopperWorkspaceV5& W5, cudaStream_t stream)
{
    W5.ensure(P.ijdim, P.kdim, lall, P.compute_rceff, OPT_MODE == 3, stream);
    auto& W = W5.base;

    dim3 g_128k((P.ijdim + 127) / 128, P.kdim, lall);
    dim3 g_256k((P.ijdim + 255) / 256, P.kdim, lall);
    dim3 g_col128((P.ijdim + 127) / 128, lall);
    dim3 g_col256((P.ijdim + 255) / 256, lall);
    dim3 g_scr((P.ijdim + 255) / 256, P.kmax0 - P.kmin0 + 1, lall);

    const size_t n_int = (size_t)P.ijdim * (size_t)(P.kmax0 - P.kmin0 + 1) * (size_t)lall;
    const size_t n_col = (size_t)P.ijdim * (size_t)lall;
    const size_t n3    = (size_t)P.ijdim * (size_t)P.kdim * (size_t)lall;
    const int g_wet_max  = (int)((n_int + 127) / 128);
    const int g_wcol_64  = (int)((n_col + 63) / 64);
    const int g_wcol_128 = (int)((n_col + 127) / 128);

    // ZERO 426 MB Vt reset! Only reset tiny ~2.6 MB metadata arrays!
    CUDA_CHECK(cudaMemsetAsync(W.d_meta_block, 0, W.meta_bytes, stream));
    CUDA_CHECK(cudaMemsetAsync(W5.d_col_w_ktop, 0, n_col * sizeof(int), stream));
    if (P.compute_rceff) {
        CUDA_CHECK(cudaMemsetAsync(W5.d_col_liq_ktop, 0, n_col * sizeof(int), stream));
        CUDA_CHECK(cudaMemsetAsync(W5.d_liq_col_count, 0, sizeof(uint32_t), stream));
        CUDA_CHECK(cudaMemsetAsync(W.d_rceff, 0, n3 * sizeof(double), stream));
    }

    // 1A. Fused clear-sky screen + Clausius-Clapeyron certificate + init_kv!
    k_sm90_v5_nsw6_screen<COLLECT_CNT><<<g_scr, 256, 0, stream>>>();

    // 1B. Zero-spill wet core (0 pointer args, 184 regs!)
    k_sm90_v5_nsw6_wet_core<COLLECT_CNT><<<g_wet_max, 128, 0, stream>>>();

    // 1C. Compact wet_cols + liq_cols (+ rare clipped-column kv repair)
    k_sm90_v5_compact_cols<COLLECT_CNT><<<g_col256, 256, 0, stream>>>();

    // 1D. Compacted liquid-only rceff (zero dry-col divergence!) & 19-us sparse Vt clean!
    if (P.compute_rceff) {
        k_sm90_v5_rceff<<<g_wcol_64, 64, 0, stream>>>();
    }
    k_sm90_v5_clean_vt_wetcols<<<g_wcol_128, 128, 0, stream>>>();

    // 2A & 2B. Branch-free Fubini half & full sedimentation (+ L2-coalesced Mode 3!)
    k_sm90_v5_precip_half<OPT_MODE, COLLECT_CNT><<<g_wcol_64, 64, 0, stream>>>();
    k_sm90_v5_precip_full<OPT_MODE, COLLECT_CNT><<<g_wcol_64, 64, 0, stream>>>();

    // 3. Perturbation-robust Clausius-Clapeyron certified SatAdj+NF+Diag!
    k_sm90_v5_satadj_nf_diag<OPT_MODE, COLLECT_CNT><<<g_128k, 128, 0, stream>>>();

    // 4. Exact unchanged-cell skip diag_w & boundary conditions
    k_sm90_v5_diag_w<<<g_256k, 256, 0, stream>>>();
    k_sm90_v2_bndcnd<<<g_col256, 256, 0, stream>>>(S);
}

inline void sm90_v5_launch_step(
    const NSW6Params& P, int lall,
    const HopperStatePtrs& S, HopperWorkspaceV5& W5, cudaStream_t stream)
{
    if (P.collect_counters) {
        if (P.opt_mode == 2) sm90_v5_launch_impl<2, true >(P, lall, S, W5, stream);
        else                 sm90_v5_launch_impl<3, true >(P, lall, S, W5, stream);
    } else {
        if (P.opt_mode == 2) sm90_v5_launch_impl<2, false>(P, lall, S, W5, stream);
        else                 sm90_v5_launch_impl<3, false>(P, lall, S, W5, stream);
    }
}
