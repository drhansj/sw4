#!/bin/bash -l

export KMP_BLOCKTIME=infinite

# Num used hw threads per core
numht=1

# For KNL
numcores=64
# Num of cores x 4
numtotalht=256

# do the loop over mpi ranks
for numranks in 1 2 4 8 16 32 64; do

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
  MPI_COMMAND="mpirun -np ${numranks}" 

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


# # do the loop
# for numprocs in 1 2 4 8 16 32 64; do
# # for numprocs in 1 2 4; do
# 
# # Number of ranks x threads is constant
# numomp=$(echo ${numprocs} | awk '{print 64/$1}')
# numht=2
# 
#     echo "Running with ${numprocs} MPI ranks and ${numomp} cores, ${numht} threads per core."
# 
#     # Intel OpenMP runtime parameters
#     export OMP_NUM_THREADS=$(( ${numomp}*${numht} ))
#     export KMP_PLACE_THREADS=1s${numomp}c${numht}t
# 
#     # Run the job with this MPI + OpenMP configuration
#     MPI_COMMAND="mpirun -np ${numprocs}"
#     RUN_COMMAND="numactl -m 1 ./sw4lite_c uni.in"
#     COMMAND="${MPI_COMMAND} ${RUN_COMMAND}" 
#     echo ${COMMAND}
# #    echo ${COMMAND} > timing_mpi${numprocs}_omp${numomp}_ht${numht}.txt
#     ${COMMAND} > mpivsomp_mpi${numprocs}_omp${numomp}_ht${numht}.txt
# 
# done
