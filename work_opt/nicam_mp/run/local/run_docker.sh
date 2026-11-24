#!/bin/bash -e
# Host-side entry point to run the NSW6 mini-app inside fortran-nvhpc:26.3.
# (Lives in run/local/ — mounts the mini-app root two levels up.)
#
#   BIN=bin_cpu ./run_docker.sh                  # CPU build, region0
#   GPU=T BIN=bin_gpu CNF=mp_miniapp_full.cnf ./run_docker.sh   # GPU build, full size
#
# Env:
#   BIN     directory holding the mp_miniapp binary (bin_cpu | bin_gpu | bin)
#   CNF     config file in this directory (default mp_miniapp.cnf)
#   GPU     'T' to pass --gpus all to docker
#   OMP     OMP_NUM_THREADS (default 1)

IMAGE=${IMAGE:-fortran-nvhpc:26.3}
RUN_DIR=$(cd "$(dirname "$0")" && pwd)
TOP_DIR=$(cd "${RUN_DIR}/../.." && pwd)

BIN=${BIN:-bin_cpu}
CNF=${CNF:-mp_miniapp.cnf}
OMP=${OMP:-1}

GPUOPT=()
if [ "${GPU:-F}" = "T" ]; then GPUOPT=(--gpus all); fi

docker run --rm "${GPUOPT[@]}" \
  --user "$(id -u):$(id -g)" \
  -v "${TOP_DIR}:/work" -w /work/run/local \
  -e OMP_NUM_THREADS="${OMP}" \
  "${IMAGE}" bash -lc "
    source /etc/profile.d/lmod.sh 2>/dev/null || true
    module load nvhpc-hpcx/26.3
    mkdir -p output
    mpirun -bind-to none -n 1 /work/${BIN}/mp_miniapp ${CNF}
  "
