#!/bin/bash
#SBATCH -N 1
#SBATCH -p regular
#SBATCH -A m1759
#SBATCH -t 01:30:00 
#SBATCH -J rahuljob
#SBATCH -o rahuljob.o%j
#SBATCH -C knl,quad,cache

#echo "Running with 8 MPI ranks and 32 logical cores per rank, ${8*32} threads per node."

# OpenMP runtime parameters
export OMP_NUM_THREADS=68
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

srun -n 1 -c 272 --cpu_bind=cores advixe-cl --collect survey -data-limit=10240 --project-dir topo_advisor --search-dir src:r=../../optimize_mp_c_cori -- ./sw4lite_c skinny-rev.in
