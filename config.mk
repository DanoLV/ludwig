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
TARGET  = nvcc

HAVE_PETSC = true
PETSC_INC  = -I/usr/local/petsc-cuda/include #-I/home/bater/Sim/petsc/arch-linux-c-opt/include
PETSC_LIB  = -L/usr/local/petsc/lib -lpetsc #-I/home/bater/Sim/petsc/arch-linux-c-opt/lib

# Original CPU-only configuration
# CC      = mpicc -fopenmp
# CFLAGS  = -O3 -march=native -Wall -DNSIMDVL=4 -DNDEBUG \
#          -fno-fast-math -fno-associative-math

# CUDA configuration for RTX 4070 (sm_89)
CC      = nvcc
CFLAGS  = -ccbin=/usr/local/ompi/bin/mpicc -O3 -DADDR_SOA -DNSIMDVL=4 -DNDEBUG -arch=sm_89 -x cu -dc

AR      = ar
ARFLAGS = -cr
LDFLAGS = -arch=sm_89

MPI_INC_PATH = -I/usr/local/ompi/include
MPI_LIB_PATH = -L/usr/local/ompi/lib -Xlinker -rpath -Xlinker /usr/local/ompi/lib -lmpi

LAUNCH_MPIRUN_CMD = mpirun -np 1 -mca pml ucx


# BUILD   = serial                  # here "serial"
# MODEL   = -D_D3Q19_               # preprocessor macro for model
# TARGET  =

# CC      = mpicc -fopenmp
# CFLAGS  = -O2 -g -Wall