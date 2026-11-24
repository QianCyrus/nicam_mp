// =============================================================================
// nsw6_cuda_kernels.cu
//
// Production Hopper SM90a (NVIDIA GH200 / H200) Strict-IEEE FP64 CUDA C++
// Operator Suite for NICAM NSW6 Cloud Microphysics (ANZIAM 2027).
//
// Backported directly from the validated, zero-spill Hopper SM90a v4.1 suite
// (4.069 ms/step with rceff, 3.897 ms/step miniapp on H200, N_FAIL = 0).
// Called directly from Fortran mod_mp_opt_bridge.f90 via iso_c_binding on the
// active OpenACC CUDA stream (acc_get_cuda_stream(1)) with ZERO host-device sync!
// =============================================================================

#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "nsw6_hopper_v4.cuh"

namespace {
static HopperWorkspaceV3 g_ws;
static bool g_grid_uploaded = false;
static int  g_grid_kdim     = 0;
static const double* g_last_rdgz = nullptr;
} // namespace

extern "C" {

void nsw6_opt_reset_counters(void* stream_ptr) {
    cudaStream_t stream = reinterpret_cast<cudaStream_t>(stream_ptr);
    if (g_ws.d_cnt) {
        cudaMemsetAsync(g_ws.d_cnt, 0, sizeof(NSW6Counters), stream);
    }
}

void nsw6_opt_get_counters(NSW6Counters* out_cnt, void* stream_ptr) {
    cudaStream_t stream = reinterpret_cast<cudaStream_t>(stream_ptr);
    if (g_ws.d_cnt && out_cnt) {
        cudaMemcpyAsync(out_cnt, g_ws.d_cnt, sizeof(NSW6Counters), cudaMemcpyDeviceToHost, stream);
        cudaStreamSynchronize(stream);
    }
}

// Batched multi-region or single-region (lall >= 1) Hopper SM90a v4.1 entry point
void nsw6_opt_driver_launch_batched(
    const NSW6Params* P_host,
    int lall,
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
    void* stream_ptr)
{
    NSW6Params P = *P_host;
    cudaStream_t stream = reinterpret_cast<cudaStream_t>(stream_ptr);

    // 1. Stream-ordered constant-memory parameter upload (0 host sync!)
    CUDA_CHECK(cudaMemcpyToSymbolAsync(c_P, &P, sizeof(NSW6Params), 0, cudaMemcpyHostToDevice, stream));

    // Upload vertical grid arrays from OpenACC device pointers ONCE per process
    if (!g_grid_uploaded || P.kdim != g_grid_kdim || GRD_rdgz != g_last_rdgz) {
        size_t kb = (size_t)P.kdim * sizeof(double);
        CUDA_CHECK(cudaMemcpyToSymbolAsync(c_rdgz,  GRD_rdgz,  kb, 0, cudaMemcpyDeviceToDevice, stream));
        CUDA_CHECK(cudaMemcpyToSymbolAsync(c_rdgzh, GRD_rdgzh, kb, 0, cudaMemcpyDeviceToDevice, stream));
        CUDA_CHECK(cudaMemcpyToSymbolAsync(c_afact, GRD_afact, kb, 0, cudaMemcpyDeviceToDevice, stream));
        CUDA_CHECK(cudaMemcpyToSymbolAsync(c_bfact, GRD_bfact, kb, 0, cudaMemcpyDeviceToDevice, stream));
        CUDA_CHECK(cudaMemcpyToSymbolAsync(c_cfact, GRD_cfact, kb, 0, cudaMemcpyDeviceToDevice, stream));
        CUDA_CHECK(cudaMemcpyToSymbolAsync(c_dfact, GRD_dfact, kb, 0, cudaMemcpyDeviceToDevice, stream));
        g_grid_uploaded = true;
        g_grid_kdim     = P.kdim;
        g_last_rdgz     = GRD_rdgz;
    }

    // 2. Pack state device pointers into constant memory c_S (frees 64 regs!)
    HopperStatePtrs S = {};
    S.rhog = rhog; S.rhogvx = rhogvx; S.rhogvy = rhogvy; S.rhogvz = rhogvz;
    S.rhogw = rhogw; S.rhoge = rhoge; S.rhogq = rhogq;
    S.rho = rho; S.pre = pre; S.tem = tem;
    S.vx = vx; S.vy = vy; S.vz = vz; S.w = w; S.q = q;
    S.precip = precip; S.precip_lh = precip_lh; S.precip_rhoein = precip_rhoein;
    S.precip_rhophi = precip_rhophi; S.precip_rhokin = precip_rhokin; S.precip_energy = precip_energy;
    S.gdcfrc = gdcfrc; S.gdclw = gdclw;
    S.unccn = unccn; S.factor_n_aie = factor_n_aie; S.z = z; S.zh = zh;
    S.ix = ix; S.iy = iy; S.iz = iz; S.jx = jx; S.jy = jy; S.jz = jz;
    S.phi = phi; S.gam2 = gam2; S.gam2h = gam2h; S.gsgam2 = gsgam2; S.gsgam2h = gsgam2h;
    sm90_v4_upload_state_ptrs(S, stream);

    // 3. Dispatch zero-spill Hopper SM90a v4.1 pipeline (Modes 1, 2, 3)
    for (int m = 0; m < P.mp_div_num; ++m) {
        sm90_v4_launch_step(P, lall, S, g_ws, stream);
    }
}

// Drop-in single-region ABI called by mod_mp_opt_bridge.f90
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
    void* stream_ptr)
{
    nsw6_opt_driver_launch_batched(
        P_host, 1,
        rhog, rhogvx, rhogvy, rhogvz, rhogw, rhoge, rhogq,
        rho, pre, tem, vx, vy, vz, w, q,
        precip, precip_lh, precip_rhoein, precip_rhophi, precip_rhokin, precip_energy,
        gdcfrc, gdclw,
        unccn, factor_n_aie, z, zh,
        ix, iy, iz, jx, jy, jz,
        phi, gam2, gam2h, gsgam2, gsgam2h,
        GRD_rdgz, GRD_rdgzh, GRD_afact, GRD_bfact, GRD_cfact, GRD_dfact,
        stream_ptr);
}

} // extern "C"
