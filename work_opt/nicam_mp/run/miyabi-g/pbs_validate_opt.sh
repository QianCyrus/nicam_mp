#!/bin/bash
################################################################################
# Miyabi-G (GH200) Validation & Math Work-Counter Gate for ANZIAM 2027
# Validates Modes 0 (BB), 1 (ENG), 2 (ENG+MATH) on real NICAM `full` dataset
# Gate: N_FAIL == 0 across all 23 arrays (34 rows, 360.6M elements, rel 1e-9 / abs 1e-10)
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

# Token budget guard
command -v show_token >/dev/null 2>&1 && show_token || true

BIN=../../bin_gpu_opt
# Build full-scale validation config with official GPU tolerance envelope (rel 1e-9, abs 1e-10)
sed -e "s|num_repeat *= *[0-9]*|num_repeat    = 1|" \
    -e "s|tolerance *= *[0-9.D+-]*|tolerance     = 1.0D-9|" \
    -e "s|tolerance_abs *= *[0-9.D+-]*|tolerance_abs = 1.0D-10|" \
    -e "s|validate *= *\.[a-zA-Z]*\.|validate      = .true.|" \
    mp_miniapp_full.cnf > mp_miniapp_full_val.cnf

for MODE in 0 1 2; do
  echo "============================================================"
  echo " [VALIDATE] NICAM_OPT_MODE=${MODE} (full 10 regions, rceff=1)"
  echo "============================================================"
  NICAM_OPT_MODE=${MODE} NICAM_COMPUTE_RCEFF=1 NICAM_COLLECT_COUNTERS=1 \
    mpirun -n 1 ${BIN}/mp_miniapp mp_miniapp_full_val.cnf | tee logs/validate_mode${MODE}.log
  grep -q "VALIDATION PASS" logs/validate_mode${MODE}.log
done

echo "ALL EXACT VALIDATION GATES PASSED (Modes 0, 1, 2 -> N_FAIL = 0)!"
