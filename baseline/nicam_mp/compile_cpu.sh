#!/bin/bash -e
# Fortran CPU/OpenMP-only build of the MP mini-app for Miyabi-G.
# CPU/OpenMP build (ENABLE_OPENMP=T, ENABLE_OPENACC=F), own BINDIR + log.

VERSION=25.9

module load nvidia/${VERSION}
module load nv-hpcx
module load netcdf
module load netcdf-fortran
module load hdf5

export ENABLE_OPENMP=T          # <- was F
export ENABLE_OPENACC=F         # <- was T
export ENABLE_NETCDF=T

export NEW_BINDIR=./bin_cpu    # <- own dir, don't clobber bin_gpu
export NICAM_SYS=Miyabi-G            # same sysdep; the Makedef honors the toggles

export MAKELOG="make_cpu.log"        # <- own log

MY_LIBPATH=/work/opt/local/aarch64/apps/nvidia/${VERSION}/netcdf-fortran/4.6.1/lib
MY_INCLUDEPATH=/work/opt/local/aarch64/apps/nvidia/${VERSION}/netcdf-fortran/4.6.1/include
export NETCDF_LIBS="-L${MY_LIBPATH} -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lm -lz"
export NETCDF_INCLUDE="-I${MY_INCLUDEPATH}"
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:${MY_LIBPATH}"

make allclean                   # REQUIRED between GPU<->CPU builds (flags differ per object)
make -j 8 &> ${MAKELOG}   # capped: unbounded -j gets killed on Miyabi login nodes
rm -rf ${NEW_BINDIR}
cp -rf ./bin ${NEW_BINDIR}
