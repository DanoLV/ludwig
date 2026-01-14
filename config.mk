##############################################################################
#
#  unix-mpicc-default.mk
#
#  Compile for parallel execution assuming an mpi compiler
#  wrapper "mpicc" is available.
#
##############################################################################

# BUILD   = serial                  # here "serial"
# MODEL   = -D_D3Q19_               # preprocessor macro for model
# TARGET  =
# CC      = mpicc -fopenmp
# CFLAGS  = -O0 -g -Wall

BUILD   = parallel
MODEL   = -D_D3Q19_

# PETSc configuration 
HAVE_PETSC = true
PETSC_INC  = -I/usr/local/petsc-cuda/include 
PETSC_LIB  = -L/usr/local/petsc-cuda/lib -Xlinker -rpath -Xlinker /usr/local/petsc-cuda/lib -lpetsc 
# PETSC_LIB  = -L/usr/local/petsc-cuda/lib -lpetsc 

# # Original CPU-only configuration
# TARGET  =
# CC      = mpicc -fopenmp
# # CFLAGS  = -O3 -march=native -Wall -DNSIMDVL=4 -DNDEBUG \
# #         -fno-fast-math -fno-associative-math
# CFLAGS  = -O0 -g -Wall # Debug flags

# CUDA configuration for RTX 4070 (sm_89)
TARGET  = nvcc
CC      = nvcc
CFLAGS  = -ccbin=/usr/local/ompi/bin/mpicc \
			-O3 -DADDR_SOA -DNSIMDVL=4 -DNDEBUG \
			-arch=sm_89 -x cu -dc
# # Para Debug
# CFLAGS = -ccbin=/usr/local/ompi/bin/mpicc \
#          -O0 -g -G -DADDR_SOA -DNSIMDVL=4 \
#          -arch=sm_89 -x cu -dc \
#          -Xcompiler -Wall  # Debug flags
AR      = ar
ARFLAGS = -cr
LDFLAGS = -arch=sm_89

# MPI configuration
MPI_INC_PATH = -I/usr/local/ompi/include
MPI_LIB_PATH = -L/usr/local/ompi/lib -Xlinker -rpath -Xlinker /usr/local/ompi/lib -lmpi

LAUNCH_MPIRUN_CMD = mpirun -np 1 -mca pml ucx
