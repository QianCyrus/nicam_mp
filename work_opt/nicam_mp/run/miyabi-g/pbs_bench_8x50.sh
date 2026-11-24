#!/bin/bash
################################################################################
# Miyabi-G (GH200) Paired Balanced 8-Block x 50-Repeat Ablation Benchmark
# Clean 4-way attribution on real NICAM `full` (166,410 cols x 80 lvls):
#   Mode 0: BB            (pristine Fortran/OpenACC baseline)
#   Mode 1: ENG           (engineering-only fused SM90, original math)
#   Mode 2: ENG+MATH      (+ Roh dyadic + Laguerre-HH4 + Fubini + stream compaction)
#   Mode 3: ENG+MATH+MONO (+ monotone backward-envelope departure projection)
################################################################################
#PBS -W group_list=jh260050
#PBS -q debug-g
#PBS -l select=1:mpiprocs=1
#PBS -l walltime=00:25:00
#PBS -j oe

set -euo pipefail
module purge
module load nvidia/25.9
module load nv-hpcx/25.9
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}:/work/opt/local/aarch64/apps/nvidia/25.9/netcdf-fortran/4.6.1/lib:/work/opt/local/aarch64/apps/nvidia/25.9/netcdf/4.9.2/lib:/work/opt/local/aarch64/apps/nvidia/25.9/hdf5/1.14.6/lib"

cd "$PBS_O_WORKDIR"
mkdir -p output logs

command -v show_token >/dev/null 2>&1 && show_token || true

BIN=../../bin_gpu_opt
sed -e "s|num_repeat *= *[0-9]*|num_repeat    = 50|" \
    -e "s|validate *= *\.[a-zA-Z]*\.|validate      = .false.|" \
    mp_miniapp_full.cnf > mp_miniapp_full_bench50.cnf

MODES=(0 1 2 3)
echo "block,mode,rceff,mp_driver_ms" > logs/bench_8x50_gh200.csv

for BLK in $(seq 0 7); do
  for IDX in 0 1 2 3; do
    MODE=${MODES[$(( (BLK + IDX) % 4 ))]}
    LOG="logs/bench_b${BLK}_m${MODE}.log"
    NICAM_OPT_MODE=${MODE} NICAM_COMPUTE_RCEFF=1 NICAM_COLLECT_COUNTERS=0 \
      mpirun -n 1 ${BIN}/mp_miniapp mp_miniapp_full_bench50.cnf > "${LOG}" 2>&1
    MS=$(grep -i "mp_driver" "${LOG}" | tail -n 1 | awk '{print $(NF-1)}' || echo "NA")
    echo "Block ${BLK} Mode ${MODE}: ${MS}"
    echo "${BLK},${MODE},1,${MS}" >> logs/bench_8x50_gh200.csv
  done
done

echo "Benchmark finished! Results in logs/bench_8x50_gh200.csv"
command -v show_token >/dev/null 2>&1 && show_token || true
