// =============================================================================
// nsw6_cpu_bridge.cpp
//
// Exact CPU host driver invoking the IDENTICAL nsw6_math_core.cuh functions
// used on NVIDIA GH200! Enables bit-exact differential verification against
// pristine Fortran baseline/nicam_mp right on macOS (gfortran + clang++).
// =============================================================================

#include <vector>
#include <cstring>
#include "nsw6_math_core.cuh"

namespace {
static NSW6Counters g_cpu_cnt = {};
static std::vector<double> g_cpu_Vt;
static std::vector<double> g_cpu_rceff;
static std::vector<double> g_cpu_rctop;
}

extern "C" {

void nsw6_opt_reset_counters(void* /*stream_ptr*/) {
    std::memset(&g_cpu_cnt, 0, sizeof(g_cpu_cnt));
}

void nsw6_opt_get_counters(NSW6Counters* out_cnt, void* /*stream_ptr*/) {
    if (out_cnt) *out_cnt = g_cpu_cnt;
}

void nsw6_opt_driver_launch(
    const NSW6Params* P_host,
    double* rhog, double* rhogvx, double* rhogvy, double* rhogvz,
    double* rhogw, double* rhoge, double* rhogq,
    double* rho, double* pre, double* tem,
    double* vx, double* vy, double* vz, double* w, double* q,
    double* precip, double* precip_lh, double* precip_rhoein,
    double* precip_rhophi, double* precip_rhokin, double* precip_energy,
    double* gdcfrc, double* gdclw,
    const double* unccn, const double* factor_n_aie,
    const double* z, const double* zh,
    const double* ix, const double* iy, const double* iz,
    const double* jx, const double* jy, const double* jz,
    const double* phi, const double* gam2, const double* gam2h,
    const double* gsgam2, const double* gsgam2h,
    const double* GRD_rdgz, const double* GRD_rdgzh,
    const double* GRD_afact, const double* GRD_bfact,
    const double* GRD_cfact, const double* GRD_dfact,
    void* /*stream_ptr*/)
{
    const NSW6Params& P = *P_host;
    size_t n3 = (size_t)P.ijdim * (size_t)P.kdim;
    if (g_cpu_Vt.size() < 4 * n3) g_cpu_Vt.resize(4 * n3);
    if (P.compute_rceff && g_cpu_rceff.size() < n3) {
        g_cpu_rceff.resize(n3);
        g_cpu_rctop.resize((size_t)P.ijdim);
    }

    std::vector<double> s_zh(P.kdim), s_z(P.kdim), s_kv(P.kdim), s_zdis(4 * P.kdim);

    for (int m = 0; m < P.mp_div_num; ++m) {
        // 1. Pre-negative_filter + NSW6 Core + state update
        for (int k = 0; k < P.kdim; ++k) {
            for (int ij = 0; ij < P.ijdim; ++ij) {
                int c_dry = 0, c_snow = 0, c_exp = 0, c_pow = 0, c_log = 0;
                nsw6_pre_and_core_cell(
                    P, ij, k,
                    rhog, rhoge, rhogq, rho, tem, pre, q,
                    w, unccn, factor_n_aie, gsgam2, g_cpu_Vt.data(),
                    &c_dry, &c_snow, &c_exp, &c_pow, &c_log);
                if (P.collect_counters && k >= P.kmin0 && k <= P.kmax0) {
                    g_cpu_cnt.total_cells       += 1;
                    g_cpu_cnt.dry_clear_cells   += (unsigned long long)c_dry;
                    g_cpu_cnt.snow_active_cells += (unsigned long long)c_snow;
                    g_cpu_cnt.n_exp_calls       += (unsigned long long)c_exp;
                    g_cpu_cnt.n_pow_calls       += (unsigned long long)c_pow;
                    g_cpu_cnt.n_log_calls       += (unsigned long long)c_log;
                }
            }
        }

        // Optional: effective radius column kernel
        if (P.compute_rceff && m == P.mp_div_num - 1) {
            for (int ij = 0; ij < P.ijdim; ++ij) {
                rceff_column(P, ij, rho, tem, q, unccn, factor_n_aie,
                             g_cpu_rceff.data(), g_cpu_rctop.data());
            }
        }

        // 2. Multispecies Precipitation Transport + KIN2EIN
        for (int ij = 0; ij < P.ijdim; ++ij) {
            int c_dry_col = 0, c_layers = 0, c_cf = 0, c_ch = 0, c_ov = 0;
            precip_transport_column(
                P, ij, 1,
                s_zh.data(), s_z.data(), s_kv.data(), s_zdis.data(),
                rhog, rhogvx, rhogvy, rhogvz, rhogw, rhoge, rhogq,
                rho, tem, pre, vx, vy, vz, q, z, zh, g_cpu_Vt.data(),
                precip, precip_lh, precip_rhoein, precip_rhophi, precip_rhokin, precip_energy,
                gsgam2, gsgam2h, gam2, gam2h,
                ix, iy, iz, jx, jy, jz,
                GRD_rdgz, GRD_rdgzh, GRD_afact, GRD_bfact, GRD_cfact, GRD_dfact,
                &c_dry_col, &c_layers, &c_cf, &c_ch, &c_ov);
            if (P.collect_counters) {
                g_cpu_cnt.dry_precip_cols       += (unsigned long long)c_dry_col;
                g_cpu_cnt.precip_layers_touched += (unsigned long long)c_layers;
                g_cpu_cnt.precip_crossings_full += (unsigned long long)c_cf;
                g_cpu_cnt.precip_crossings_half += (unsigned long long)c_ch;
                g_cpu_cnt.precip_top_overshoots += (unsigned long long)c_ov;
            }
        }

        // 3. Fused Laguerre-Householder-4 Saturation Adjustment + Post-negative_filter
        for (int k = 0; k < P.kdim; ++k) {
            for (int ij = 0; ij < P.ijdim; ++ij) {
                int c_act = 0, c_it = 0, c_exp = 0, c_pow = 0, c_neg = 0;
                satadj_and_post_negfilter_cell(
                    P, ij, k,
                    rhog, rhoge, rhogq, rho, tem, pre, q, gsgam2,
                    &c_act, &c_it, &c_exp, &c_pow, &c_neg);
                if (P.collect_counters && k >= P.kmin0 && k <= P.kmax0) {
                    g_cpu_cnt.sat_active_cells        += (unsigned long long)c_act;
                    g_cpu_cnt.sat_total_iters         += (unsigned long long)c_it;
                    g_cpu_cnt.n_exp_calls             += (unsigned long long)c_exp;
                    g_cpu_cnt.n_pow_calls             += (unsigned long long)c_pow;
                    g_cpu_cnt.neg_q_before_postfilter += (unsigned long long)c_neg;
                }
            }
        }
    }

    // 4. Fused CNVVAR_prg2diag_in + BNDCND_thermo + GDCFRC + GDCLW
    for (int ij = 0; ij < P.ijdim; ++ij) {
        cnvvar_bndcnd_cldfrac_column(
            P, ij,
            rhog, rhogvx, rhogvy, rhogvz, rhogw, rhoge, rhogq,
            rho, pre, tem, vx, vy, vz, w, q, gdcfrc, gdclw,
            phi, gsgam2, gsgam2h, GRD_afact, GRD_bfact);
    }
}

} // extern "C"
