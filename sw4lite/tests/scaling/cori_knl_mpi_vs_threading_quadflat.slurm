#!/bin/bash
#SBATCH -N 1
#SBATCH -p debug
#SBATCH -A m2545
# #SBATCH -S 2 (FOR MULTI-NODE)
#SBATCH -t 00:10:00 
#SBATCH -C knl,quad,flat

# Sam William's recommended settings
module load craype-hugepages2M
# export KMP_AFFINITY=compact,granularity=thread,1
export KMP_BLOCKTIME=infinite

# Num used hw threads per core
numht=1

# For KNL
numcores=64
# Num of cores x 4
numtotalht=256

# do the loop over mpi ranks
for numranks in 1 2 4 8 16 32 64; do

  # Total number of ranks across all the nodes we requested with -N above
  numtotal=$(( ${numranks}*${SLURM_JOB_NUM_NODES} ))
  # Cores per mpi rank requested from srun: numranks / numtotalht
  numcorespermpi=$(echo ${numranks} ${numtotalht} | awk '{print $2/$1}')
  # OMP threads is ht * numranks / numcores
  numomp=$(( ${numht}*$(echo ${numranks} ${numcores} | awk '{print $2/$1}') ))

  echo "Running with ${numranks} MPI ranks and ${numomp} threads, using ${numht} threads per core."

  # Intel OpenMP runtime parameters
  echo "OMP_NUM_THREADS=${numomp}"
  export OMP_NUM_THREADS=${numomp}
  export OMP_PLACES=cores"(64)"

  # Run the job with this MPI + OpenMP configuration
  MPI_COMMAND="srun -n ${numranks} -c ${numcorespermpi} --cpu_bind=cores" 

  # Use this command to check OMP affinity 
  #    RUN_COMMAND="numactl -m 1 check-hybrid.intel.cori"
  # Run sw4lite with the settings in uni.in
  RUN_COMMAND="numactl -m 1 ../../optimize_mp/sw4g uni.in"

  # Echo and run the command
  COMMAND="${MPI_COMMAND} ${RUN_COMMAND}" 
  echo ${COMMAND}
  # Run the command, redirect output for each loop iteration
  ${COMMAND} > mpivsomp_mpi${numranks}_omp${numomp}_ht${numht}_qf.txt

done
