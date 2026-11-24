#!/bin/bash -e
# Runs inside fortran-nvhpc:26.3.  Invoked by compile_docker.sh / compile_docker_cpu.sh.
# Twin of nicam/compile_docker_inner.sh, minus NetCDF (the mini-app links nothing).

source /etc/profile.d/lmod.sh 2>/dev/null || true
module load nvhpc-hpcx/26.3

: "${NEW_BINDIR:=./bin_docker}"
: "${MAKELOG:=make_docker.log}"
: "${MAKEJOBS:=8}"

cd /work

make allclean >/dev/null 2>&1 || true
make -j${MAKEJOBS} &> "${MAKELOG}" || { tail -40 "${MAKELOG}"; exit 1; }

rm -rf "${NEW_BINDIR}"
cp -rf ./bin "${NEW_BINDIR}"

echo "Build finished. Binaries in ${NEW_BINDIR}"
