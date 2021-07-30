#!/bin/bash
#============ QSUB Options ============
#QSUB -q gr20001b
#QSUB -ug gr20001
#QSUB -W 300:00
#QSUB -A p=72:t=2:c=1
#--- p : # of processes
#--- t : # of threads
#--- c : # of cores (= t)

#============ Shell Script ============
set -x
cd $QSUB_WORKDIR

# automatically
# export OMP_NUM_THREADS=$QSUB_THREADS

date
mpiexec.hydra ./do.exe
date
