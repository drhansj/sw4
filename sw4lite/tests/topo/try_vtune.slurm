#!/bin/bash
#SBATCH -N 1
#SBATCH -p regular
#SBATCH -A m1759
#SBATCH -t 01:30:00 
#SBATCH -J rahuljob
#SBATCH -o rahuljob.o%j
#SBATCH -C knl,quad,flat
#SBATCH --perf=vtune

module load vtune
module load numactl

#echo "Running with 8 MPI ranks and 32 logical cores per rank, ${8*32} threads per node."

# OpenMP runtime parameters
export OMP_NUM_THREADS=64
export OMP_PLACES=cores
export OMP_PROC_BIND=spread


srun -n 1 -c 256 --cpu_bind=cores -- numactl -p 1 amplxe-cl -collect hpc-performance -data-limit=0 -r topo_vtune -- ./sw4lite_c skinny-rev.in
