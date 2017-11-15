#!/bin/bash
#SBATCH -N 64
#SBATCH -p debug
#SBATCH -A m2868
#SBATCH -S 4     # Special cores per nodes (4 are idle)
#SBATCH -t 00:30:00 
#SBATCH -L SCRATCH
#SBATCH -C knl,quad,cache
#DW jobdw capacity=10GB access_mode=striped type=scratch
#DW stage_in source=/global/cscratch1/sd/u6336/sw4-ecp/performance/rfile/USGSBayAreaVM-08.3.0-corder.rfile destination=$DW_JOB_STRIPED/USGSBayAreaVM-08.3.0-corder.rfile type=file
#DW stage_in source=/global/cscratch1/sd/u6336/sw4-ecp/performance/rfile/m6.5-20.0x13.0.s500.v5.1.srf destination=$DW_JOB_STRIPED/m6.5-20.0x13.0.s500.v5.1.srf type=file
## Don't need this, it gets generated
## #DW stage_in source=/global/cscratch1/sd/u6336/sw4-ecp/performance/rfile/hayward-h17-cori.in destination=$DW_JOB_STRIPED/hayward-h17-cori.in type=file

numnodes=64
numprocspernode=64
echo "numprocspernode=" $numprocspernode

# 1st case
CASE=hayward-h17-cori-small-bb
# Replace the paths with BB paths
echo "replacing \${INPUT_DIR} with ${DW_JOB_STRIPED} in ${CASE}.pre"
export INPUT_DIR=${DW_JOB_STRIPED}
envsubst < ${CASE}.pre > ${INPUT_DIR}/${CASE}.in
cp ${INPUT_DIR}/${CASE}.in .

# for trial in `seq 1 10`;
for trial in `seq 1 3`;
do

# for size in 32 48 64 96 128; 
for size in 32 64; 
do
numprocs=$(( ${size}*${numprocspernode} ))

# hyper-threads
numlc=1

echo "Running on ${numnodes} nodes with ${numprocs} MPI ranks."

# Run the job with this MPI + OpenMP configuration
MPI_COMMAND="srun -N ${numnodes} -n ${numprocs} -c ${numlc} --cpu_bind=cores" 

RUN_COMMAND="../sw4 ${DW_JOB_STRIPED}/${CASE}.in" 
OUT_FILE="${CASE}.out.n${numprocs}.t${trial}.knl"
echo ${MPI_COMMAND} ${RUN_COMMAND} ">&" "${OUT_FILE}"
${MPI_COMMAND} ${RUN_COMMAND} >& ${OUT_FILE}

done
done

