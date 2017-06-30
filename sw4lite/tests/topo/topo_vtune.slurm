#!/bin/bash
#SBATCH -N 1
#SBATCH -p regular
#SBATCH -A m1759
#SBATCH -t 01:30:00 
#SBATCH -J rahuljob
#SBATCH -o rahuljob.o%j
#SBATCH -C knl,quad,cache
#SBATCH --perf=vtune

module load vtune
module load numactl

#echo "Running with 8 MPI ranks and 32 logical cores per rank, ${8*32} threads per node."

# OpenMP runtime parameters
export OMP_NUM_THREADS=64
export OMP_PLACES=threads
export OMP_PROC_BIND=spread


srun -n 1 -c 256 --cpu_bind=cores amplxe-cl -collect hpc-performance -data-limit=10240 -r topo_vtune_hpcPerf -- ./sw4lite_c skinny-rev.in
srun -n 1 -c 256 --cpu_bind=cores amplxe-cl -collect general-exploration -data-limit=10240 -r topo_vtune_genExplore -- ./sw4lite_c skinny-rev.in

