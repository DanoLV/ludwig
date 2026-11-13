#!/bin/bash
export LD_LIBRARY_PATH=/usr/local/petsc-cuda-hypre/lib:/usr/local/ompi/lib:$LD_LIBRARY_PATH
cd "$(dirname "$(readlink -f "$0")")"
./Ludwig.exe > run.log 2>&1; rc=$?; cp run.log logs/output.txt
echo "exit=$rc terminado=$(grep -ac 'finished normally' run.log)"
