#!/bin/bash

source ~/sourceme.sh

export LD_LIBRARY_PATH=/panfs/users/student13/vpic_project/build/HACK_MASTER_CLE_PTH_INT_OPT_V16_AVX512_KNL_VTUNE/vpic/build/lib:${LD_LIBRARY_PATH}

MYAPP="numactl --preferred=1 ./test_lpi_np_00068.Linux --tpp 4"
PROJDIR="/panfs/users/student13/runs/ddr/Results_advixe_1"
BINSRCH="/panfs/users/student13/runs/ddr"
SRCSRCH="/panfs/users/student13/vpic_project/src"
SYMSRCH="/panfs/users/student13/vpic_project/build/HACK_MASTER_CLE_PTH_INT_OPT_V16_AVX512_KNL_VTUNE/vpic/build/lib"

#COLL="survey"
#COLLARGS="-trace-mpi -no-auto-finalize"
COLL="tripcounts"
COLLARGS="-trace-mpi -flops-and-masks -no-auto-finalize"

mpirun -n 68 -gtool "advixe-cl --collect ${COLL} ${COLLARGS} -search-dir bin:r=${BINSRCH} -search-dir src:r=${SRCSRCH} -search-dir sym:r=${SYMSRCH} -project-dir=${PROJDIR}:0" ${MYAPP}
#mpirun -n 68 -gtool "advixe-cl --collect ${COLL} ${COLLARGS} -search-dir bin:r=${BINSRCH} -search-dir src:r=${SRCSRCH} -search-dir sym:r=${SYMSRCH} -project-dir=${PROJDIR}:0" ./mpi_hello
