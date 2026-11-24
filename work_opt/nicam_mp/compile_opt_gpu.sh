#!/bin/bash -e
# Build optimized NICAM NSW6 GPU binary (NVHPC 25.9 OpenACC + CUDA 13.0 sm_90 fused operators)
# Runs on Miyabi-G login/compute node. Preserves pristine baseline/nicam_mp untouched!

VERSION=25.9
module purge || true
module load nvidia/${VERSION}
module load nv-hpcx/${VERSION} || module load nv-hpcx
module load netcdf || true
module load netcdf-fortran || true
module load hdf5 || true

CUDA_ROOT=${CUDA_HOME:-/work/opt/local/aarch64/cores/nvidia/${VERSION}/Linux_aarch64/${VERSION}/cuda/13.0}
NVCC=${CUDA_ROOT}/bin/nvcc

echo "=== [1/2] Compiling native strict-IEEE FP64 CUDA C++ operators (sm_90) ==="
mkdir -p workspace lib bin
${NVCC} -O3 -std=c++17 -arch=sm_90a --prec-div=true --prec-sqrt=true -Xcompiler -fPIC \
    -c cuda_kernels/nsw6_cuda_kernels.cu -o workspace/nsw6_cuda_kernels.o
ar rcs lib/libnsw6_cuda_opt.a workspace/nsw6_cuda_kernels.o

echo "=== [2/2] Building NICAM-MP Fortran/OpenACC + CUDA fused binary ==="
export ENABLE_OPENMP=F
export ENABLE_OPENACC=T
export ENABLE_NETCDF=T
export ENABLE_MORETIMER=${ENABLE_MORETIMER:-F}
export NICAM_SYS=Miyabi-G
export NEW_BINDIR=./bin_gpu_opt
export MAKELOG="make_gpu_opt.log"

MY_LIBPATH=/work/opt/local/aarch64/apps/nvidia/${VERSION}/netcdf-fortran/4.6.1/lib
MY_INCLUDEPATH=/work/opt/local/aarch64/apps/nvidia/${VERSION}/netcdf-fortran/4.6.1/include
export NETCDF_LIBS="-L${MY_LIBPATH} -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lm -lz"
export NETCDF_INCLUDE="-I${MY_INCLUDEPATH}"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:${MY_LIBPATH}:${CUDA_ROOT}/lib64"

export OPT_KERNEL_LIBS="$(pwd)/lib/libnsw6_cuda_opt.a -L${CUDA_ROOT}/lib64 -lcudart -lstdc++"

make clean
# Rebuild libnsw6_cuda_opt.a if make clean touched lib/
${NVCC} -O3 -std=c++17 -arch=sm_90a --prec-div=true --prec-sqrt=true -Xcompiler -fPIC \
    -c cuda_kernels/nsw6_cuda_kernels.cu -o workspace/nsw6_cuda_kernels.o
ar rcs lib/libnsw6_cuda_opt.a workspace/nsw6_cuda_kernels.o

make -j 8 OPT_KERNEL_LIBS="${OPT_KERNEL_LIBS}" &> ${MAKELOG}
rm -rf ${NEW_BINDIR}
cp -rf ./bin ${NEW_BINDIR}
echo "=== Build complete: ${NEW_BINDIR}/mp_miniapp ==="
