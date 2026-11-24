#! /bin/bash
################################################################################
#
# NSW6 cloud-microphysics mini-app on Miyabi-G (1 GPU)
#
################################################################################
#PBS -W group_list=<group name>
#PBS -q debug-g
#PBS -l select=1:mpiprocs=1
#PBS -l walltime=00:10:00
#PBS -j oe

module purge
module load nvidia/25.9
module load nv-hpcx/25.9
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/work/opt/local/aarch64/apps/nvidia/25.9/netcdf-fortran/4.6.1/lib:/work/opt/local/aarch64/apps/nvidia/25.9/netcdf/4.9.2/lib:/work/opt/local/aarch64/apps/nvidia/25.9/hdf5/1.14.6/lib

cd $PBS_O_WORKDIR
mkdir -p output

# Binary from ../../compile.sh (GPU) or ../../compile_cpu.sh (CPU).
# Config: mp_miniapp_gpu.cnf = region0 + the GPU tolerance envelope (rel 1e-9 /
# abs 1e-10, see README); use mp_miniapp_full.cnf + bin_gpu for the full-size
# performance run; a CPU run with mp_miniapp.cnf must report exactly 0.0.
BIN=${BIN:-../../bin_gpu}
CNF=${CNF:-mp_miniapp_gpu.cnf}

mpirun -n 1 ${BIN}/mp_miniapp ${CNF} || exit 1
