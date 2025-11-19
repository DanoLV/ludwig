##############################################################################
#
#  unix-mpicc-default.mk
#
#  Compile for parallel execution assuming an mpi compiler
#  wrapper "mpicc" is available.
#
##############################################################################

BUILD   = parallel
MODEL   = -D_D3Q19_
TARGET  =
# TARGET  = nvcc

HAVE_PETSC = true
PETSC_INC  = -I/usr/local/petsc-cuda/include #-I/home/bater/Sim/petsc/arch-linux-c-opt/include
PETSC_LIB  = -L/usr/local/petsc/lib -lpetsc #-I/home/bater/Sim/petsc/arch-linux-c-opt/lib

CC      = mpicc -fopenmp
# CFLAGS  = -O0 -g -Wall #-O2
CFLAGS = -O3 -march=native -Wall -DNSIMDVL=4 -DNDEBUG \
         -fno-fast-math -fno-associative-math

LAUNCH_MPIRUN_CMD = mpirun -np 1 -mca pml ucx


# BUILD   = serial                  # here "serial"
# MODEL   = -D_D3Q19_               # preprocessor macro for model
# TARGET  =

# CC      = mpicc -fopenmp
# CFLAGS  = -O2 -g -Wall