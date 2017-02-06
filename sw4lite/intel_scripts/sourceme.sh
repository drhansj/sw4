#!/bin/bash
source /opt/intel/compiler/2017/bin/compilervars.sh intel64
source /opt/intel/vtune/2017/vtune_amplifier_xe/amplxe-vars.sh 
source /opt/intel/advisor/2017u1/advisor_2017/advixe-vars.sh 

source /opt/intel/impi/5.1.3.181/bin/mpivars.sh

export ADVIXE_EXPERIMENTAL=roofline
export PATH=/opt/crtdc/gcc/gcc-4.9.2/bin:${PATH}
export LD_LIBRARY_PATH=/opt/crtdc/gcc/gcc-4.9.2/lib64:${LD_LIBRARY_PATH}
export INCLUDE=/opt/crtdc/gcc/gcc-4.9.2/include/c++/4.9.2:${INCLUDE}
export MANPATH=/opt/crtdc/gcc/gcc-4.9.2/man:${MANPATH}

# cmake 3.5.2
export PATH=/opt/crtdc/cmake/3.0.2/bin:${PATH}

export SCRATCH=/panfs/users/student15

export OMP_NUM_THREADS=4
export KMP_BLOCKTIME=0
export I_MPI_PIN_DOMAIN="4:compact"
export I_MPI_DEBUG=5
export KMP_AFFINITY="verbose"
