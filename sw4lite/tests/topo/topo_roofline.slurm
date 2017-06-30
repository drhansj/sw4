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
export OMP_NUM_THREADS=64
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

srun -n 1 -c 256 --cpu_bind=cores advixe-cl --collect survey --trace-mpi -data-limit=0 --project-dir topo_roofline --search-dir src:r=../../optimize_mp_c_cori -- ./sw4lite_c skinny-rev.in
srun -n 1 -c 256 --cpu_bind=cores advixe-cl --collect tripcounts -flops-and-masks -no-auto-finalize -data-limit=0 --project-dir topo_roofline -- ./sw4lite_c skinny-rev.in
