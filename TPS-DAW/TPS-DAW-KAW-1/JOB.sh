#!/bin/bash
#============ QSUB Options ============
#QSUB -q gr20001a
#QSUB -ug gr20001
#QSUB -W 300:00
#QSUB -A p=1:t=8:c=2
#--- p : # of processes
#--- t : # of threads
#--- c : # of cores (= t)

#============ Shell Script ============
set -x

# automatically
# export OMP_NUM_THREADS=$QSUB_THREADS

date
aprun -n $QSUB_PROCS -d $QSUB_THREADS -N $QSUB_PPN ./do.exe
date
