#!/bin/bash
################################################################################
# Miyabi-G (GH200) Continuous Multi-Step Conservation Experiment (5 & 100 steps)
# Quantifies water mass & Ooyama total moist energy conservation of:
#   Mode 0 (BB) / Mode 1 (ENG) / Mode 2 (ENG+MATH) vs Mode 3 (ENG+MATH+MONO)
# on the real NICAM `full` (166,410 cols x 80 lvls) atmospheric state!
################################################################################
#PBS -W group_list=jh260050
#PBS -q debug-g
#PBS -l select=1:mpiprocs=1
#PBS -l walltime=00:15:00
#PBS -j oe

set -euo pipefail
module purge
module load nvidia/25.9
module load nv-hpcx/25.9
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:/work/opt/local/aarch64/apps/nvidia/25.9/netcdf-fortran/4.6.1/lib:/work/opt/local/aarch64/apps/nvidia/25.9/netcdf/4.9.2/lib:/work/opt/local/aarch64/apps/nvidia/25.9/hdf5/1.14.6/lib"

cd "$PBS_O_WORKDIR"
mkdir -p output logs

BIN=../../bin_gpu_opt

for NSTEPS in 5 100; do
  sed -e "s|num_repeat *= *[0-9]*|num_repeat    = ${NSTEPS}|" \
      -e "s|validate *= *\.[a-zA-Z]*\.|validate      = .false.|" \
      mp_miniapp_full.cnf > mp_miniapp_full_cont${NSTEPS}.cnf

  for MODE in 0 1 2 3; do
    echo "============================================================"
    echo " [MULTISTEP] NSTEPS=${NSTEPS}, NICAM_OPT_MODE=${MODE}"
    echo "============================================================"
    NICAM_CONTINUOUS_STEPS=1 NICAM_OPT_MODE=${MODE} NICAM_COMPUTE_RCEFF=0 NICAM_COLLECT_COUNTERS=1 \
      mpirun -n 1 ${BIN}/mp_miniapp mp_miniapp_full_cont${NSTEPS}.cnf | tee logs/multistep_${NSTEPS}_m${MODE}.log
  done
done
