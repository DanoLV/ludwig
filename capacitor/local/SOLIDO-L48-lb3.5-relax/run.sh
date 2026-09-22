#!/bin/bash
# Lanza la relajacion en segundo plano:  ./run.sh
# Seguimiento:  tail -f run.log   |   python3 convergencia.py
cd "$(dirname "$0")" || exit 1

export LD_LIBRARY_PATH=/usr/local/petsc-cuda-hypre/lib:/usr/local/ompi/lib:$LD_LIBRARY_PATH
export OMP_NUM_THREADS=1
export CUDA_VISIBLE_DEVICES=0

# ludwig.c exige ./proceced_data/ y escribe alli, EN CADA PASO y por particula,
# particle_Esub.csv / particle_force.csv (GB para 2544 particulas). Sin carga
# no tienen informacion: se descartan enlazandolos a /dev/null.
mkdir -p logs proceced_data
ln -sf /dev/null proceced_data/particle_Esub.csv
ln -sf /dev/null proceced_data/particle_force.csv

nohup mpirun -np 1 ./Ludwig.exe > run.log 2>&1 &
echo $! > logs/pid
echo "Ludwig lanzado, PID $(cat logs/pid). Log: run.log"
