###############################################################################
#
#  Makefile
#
#  Various targets for utility programs to perform pre-processing
#  and post-processing tasks.
#
#  Note all the builds are parallel, but the utilities are
#  intended to be run in serial.
#
#  Edinburgh Soft Matter and Statistical Physics Group and
#  Edinburgh Parallel Computing Centre
#
#  (c) 2016-2023 The University of Edinburgh
#
#  Contributing authors:
#  Kevin Stratford (kevin@epcc.ed.ac.uk)
#
###############################################################################

include ../Makefile.mk

# Compilation etc

SRC  = $(ROOT_DIR)./src
INCL = -I$(SRC) $(TARGET_INC_PATH) $(MPI_INC_PATH) $(PETSC_INC)

BLIBS = $(MPI_LIB_PATH) $(MPI_LIB) $(TARGET_LIB_PATH) $(TARGET_LIB) $(PETSC_LIB)
# libludwig.a now references cuFFT (psi_fft.o, psi_refined.o), so every
# utility that links it needs -lcufft.
LIBS  = $(SRC)/libludwig.a ${BLIBS} -lm -lcufft


default:
	$(MAKE) build

build:
	$(MAKE) extract
	$(MAKE) capillary
	$(MAKE) colloid_init
	$(MAKE) extract_colloids
	$(MAKE) coll_squ_subgrid_init
	$(MAKE) multi_poly_init
	$(MAKE) microgel_poly_cross_density
	$(MAKE) polarizer

colloid_init: colloid_init.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

extract_colloids: extract_colloids.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

extract_squirmer_polymer: extract_squirmer_polymer.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

capillary: capillary.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS) -lm

extract: extract.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

coll_squ_subgrid_init: coll_squ_subgrid_init.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

multi_poly_init: multi_poly_init.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

microgel_poly: microgel_poly.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

microgel_poly_filled: microgel_poly_filled.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

microgel_poly_cross_density: microgel_poly_cross_density.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)	
	
polarizer: polarizer.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

PM_self_force: PM_self_force.o
	$(CC) $(LDFLAGS) -o $@ $^ $(LIBS)

# Default rules

.PHONY : clean
clean:
	rm -f *.o colloid_init extract_colloids capillary extract \
	coll_squ_subgrid_init multi_poly_init polarizer
	rm -f *gcda *gcno

.SUFFIXES:
.SUFFIXES: .c .o

.c.o:
	$(CC) $(CFLAGS) $(MODEL) $(INCL) -c $*.c
