#!/bin/bash
#SBATCH -N 1
#SBATCH -p debug
#SBATCH -A m1759
#SBATCH -S 4
#SBATCH -t 00:20:00 
#SBATCH -C haswell

export OMP_PLACES=cores
export OMP_PROC_BIND=spread

for mpi_rank in 1 2 4 8 16 32; do
    threads_per_rank=$((32/${mpi_rank}))
    export OMP_NUM_THREADS=${threads_per_rank}
    logical_cores=$(( ${threads_per_rank}*2 ))
    echo "mpi_rank = ${mpi_rank} threads_per_rank = ${threads_per_rank} logical_cores = ${logical_cores}"
    echo "export OMP_NUM_THREADS=${threads_per_rank}"
    echo "srun -n ${mpi_rank} -c ${logical_cores} --cpu_bind=cores ./sw4lite_c skinny-rev.in"

    srun -n ${mpi_rank} -c ${logical_cores} --cpu_bind=cores ./sw4lite_c skinny-rev.in
done
