#!/bin/bash -e
# Host-side entry point: build the GPU (OpenACC) mini-app inside fortran-nvhpc:26.3.
# Twin of nicam/compile_docker.sh.

IMAGE=${IMAGE:-fortran-nvhpc:26.3}
PROJECT_DIR=$(cd "$(dirname "$0")" && pwd)

docker run --rm \
  --user "$(id -u):$(id -g)" \
  -v "${PROJECT_DIR}:/work" -w /work \
  -e ENABLE_OPENMP=F \
  -e ENABLE_OPENACC=T \
  -e NICAM_SYS=Linux64-nvidia-ompi-time \
  -e NEW_BINDIR=./bin_gpu \
  -e MAKELOG=make_gpu.log \
  "${IMAGE}" bash /work/compile_docker_inner.sh
