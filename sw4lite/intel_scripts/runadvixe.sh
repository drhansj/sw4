#!/bin/bash

source ~/sourceme.sh

# export LD_LIBRARY_PATH=/panfs/users/student13/vpic_project/build/HACK_MASTER_CLE_PTH_INT_OPT_V16_AVX512_KNL_VTUNE/vpic/build/lib:${LD_LIBRARY_PATH}

MYAPP="numactl --preferred=1 ./sw4lite_c small-advixe-rev.in"
PROJDIR="/home/student15/advisor/results_roofline"
BINSRCH="/home/student15/sw4-ecp/sw4lite/optimize_mp_c"
SRCSRCH="/home/student15/sw4-ecp/sw4lite/src"
SYMSRCH="/home/student15/sw4-ecp/sw4lite/optimize_mp_c"

export ADVIXE_EXPERIMENTAL=roofline

#COLL="survey"
#COLLARGS="-trace-mpi -no-auto-finalize"
COLL="tripcounts"
COLLARGS="-trace-mpi -flops-and-masks -no-auto-finalize"

#export OMP_NUM_THREADS=68
export OMP_NUM_THREADS=1
mpirun -n 1 -gtool "advixe-cl --collect ${COLL} ${COLLARGS} -search-dir bin:r=${BINSRCH} -search-dir src:r=${SRCSRCH} -search-dir sym:r=${SYMSRCH} -project-dir=${PROJDIR}:0" ${MYAPP}
#mpirun -n 68 -gtool "advixe-cl --collect ${COLL} ${COLLARGS} -search-dir bin:r=${BINSRCH} -search-dir src:r=${SRCSRCH} -search-dir sym:r=${SYMSRCH} -project-dir=${PROJDIR}:0" ./mpi_hello
