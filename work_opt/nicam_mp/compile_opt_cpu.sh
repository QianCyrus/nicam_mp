#!/bin/bash -e
# Local CPU build (Apple Silicon macOS gfortran+clang++ OR Linux gfortran+g++)
# Links the EXACT SAME nsw6_math_core.cuh used on NVIDIA GH200 into the real
# NICAM Fortran mp_miniapp for bit-exact differential verification!

mkdir -p workspace lib bin

CXX=${CXX:-clang++}
FC=${FC:-gfortran}

echo "=== [1/2] Compiling nsw6_cpu_bridge.cpp (exact same nsw6_math_core.cuh as CUDA!) ==="
${CXX} -O3 -std=c++17 -fPIC -c cuda_kernels/nsw6_cpu_bridge.cpp -o workspace/nsw6_cpu_bridge.o
ar rcs lib/libnsw6_cpu_opt.a workspace/nsw6_cpu_bridge.o

echo "=== [2/2] Building NICAM-MP Fortran library with ${FC} ==="
FFLAGS="-cpp -O3 -fconvert=big-endian -frecord-marker=4 -fno-range-check -DMINIAPP -Jworkspace -Iworkspace"

SRCS=(
  nhm/share/mod_precision.f90
  nhm/share/mod_stdio.f90
  nhm/share/mod_adm.f90
  nhm/share/mod_prof.f90
  nhm/share/mod_cnst.f90
  nhm/share/mod_misc.f90
  nhm/share/mod_grd.f90
  nhm/share/mod_ch_includes.f90
  nhm/share/mod_runconf.f90
  nhm/share/mod_thrmdyn.f90
  nhm/share/mod_bndcnd.f90
  nhm/share/mod_cnvvar.f90
  nhm/share/mod_isotope.f90
  nhm/physics/mod_satadjust.f90
  nhm/physics/mod_vadv1d.f90
  nhm/physics/mod_precip_transport.f90
  nhm/physics/mod_mp_nsw6.f90
  nhm/physics/mod_mp_opt_bridge.f90
  nhm/physics/mod_mp_driver.f90
)

OBJS=()
for src in "${SRCS[@]}"; do
  base=$(basename "${src%.f90}")
  obj="workspace/${base}.o"
  ${FC} ${FFLAGS} -c "${src}" -o "${obj}"
  OBJS+=("${obj}")
done

ar rcs lib/libnicammp_cpu.a "${OBJS[@]}"
echo "=== Built lib/libnicammp_cpu.a and lib/libnsw6_cpu_opt.a successfully! ==="
