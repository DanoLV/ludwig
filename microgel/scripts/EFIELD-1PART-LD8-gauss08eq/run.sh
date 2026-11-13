#!/bin/bash
export LD_LIBRARY_PATH=/usr/local/petsc-cuda-hypre/lib:/usr/local/ompi/lib:$LD_LIBRARY_PATH
B=$(dirname "$(readlink -f "$0")")
run() { cd $B/corr_$1 && ./Ludwig.exe > run.log 2>&1; rc=$?; cp run.log logs/output.txt
        echo "corr_$1 exit=$rc terminado=$(grep -ac 'finished normally' run.log)"; }
run no & run yes & wait
echo LISTO
