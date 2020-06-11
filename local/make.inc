################################################################
#
# Section I: Preprocessor options
#
################################################################

#0/1 = turn pardiso off/on
PARDISO=0

#0/1 = turn MUMPS off/on
ifdef MUMPS
  USEMUMPS=1
else
  USEMUMPS=0
endif

#0/1 = turn parallel MG solver off/on
USEPARMG=1

#0/1 = turn OpenMP compilation off/on
USEOPENMP=0

#0/1 = turn verification off/on
USEVERIFICATION=0

#Please set precompiler flags here
FPPFLAGS=-DPARDISO=$(PARDISO) -DMUMPS=$(USEMUMPS) -DPARMG=$(USEPARMG) -DVERIFICATION=$(USEVERIFICATION)

################################################################
#
# Section II: Compiler options
#
################################################################

#Default optimization level (O=optimized, g=debug)
OPTLEVEL = O

F90FLAGS =
CFLAGS =

ifeq ($(OPTLEVEL), O) #optimized
     ifeq ($(COMPTYPE), i) #intel
     	  F90FLAGS += -O3 -xHOST
     endif
     ifeq ($(COMPTYPE), g) #gnu
     	  F90FLAGS += -ffree-line-length-0 -O3 
     endif
     ifeq ($(COMPTYPE), c) #cray
     	  F90FLAGS +=
     endif
endif


ifeq ($(OPTLEVEL), g) #debug
     ifeq ($(COMPTYPE), i) #intel
     	  F90FLAGS += -g -traceback -CB
     endif
     ifeq ($(COMPTYPE), g) #gnu
     	  F90FLAGS += -ffree-line-length-0 -g -fbacktrace -fcheck=all -pedantic -Wall 
     endif
     ifeq ($(COMPTYPE), c) #cray
     	  F90FLAGS += -g -O0
     endif
endif


ifeq ($(USEOPENMP), 1)
     ifeq ($(COMPTYPE), i) #intel
     	  F90FLAGS += -qopenmp
	  CFLAGS   += -qopenmp
     endif
     ifeq ($(COMPTYPE), g) #gnu
     	  F90FLAGS += -fopenmp
	  CFLAGS   += -fopenmp
     endif
endif


################################################################
#
# Section III: Libraries and where to find them
#
################################################################

IDIRS := -I$(SPC_LOCAL)/include/$(OPTLEVEL)

#LIBS := -lbsplines -lfutils -lpppack -lpputils2 \
#        -lhdf5_fortran -lhdf5 -lz -ldl -lpthread
LIBS  := -lfutils -lhdf5_fortran -lhdf5 -lz -ldl -lpthread     
LDIRS := -L$(SPC_LOCAL)/lib/$(OPTLEVEL) -L$(HDF5)/lib

#
# setup of MKL library
# working well with intel compiler and intelmpi
# MKL_LIB has to be set in .bashrc 
#
ifdef MKL_LIB

     # common libraries
     LIBS += -lmkl_scalapack_lp64 -lmkl_core -lmkl_blacs_intelmpi_lp64 -lm

     # compiler dependent libraries
     ifeq ($(COMPTYPE), i) #intel
     	  LIBS += -lmkl_intel_lp64
     endif
     ifeq ($(COMPTYPE), g) #gnu
     	  LIBS += -lmkl_gf_lp64
     endif

     # different library depending on threaded or sequential
     ifeq ($(USEOPENMP), 1)
     	LIBS += -lmkl_intel_thread
     else
	LIBS += -lmkl_sequential
     endif

     # MKL needs intel omp library when compiled with gnu
     ifeq ($(USEOPENMP), 1)
	ifeq ($(COMPTYPE), g) #gnu
     	    LIBS += -liomp5
	endif
     endif

     LDIRS += -L$(MKL_LIB) 
     IDIRS += -I$(MKL_LIB)/../../include

endif

#
# Add mumps libraries if preprocessor set
#
ifeq ($(USEMUMPS),1)
     LIBS  += -lzmumps -ldmumps -lmumps_common -lpord
     IDIRS += -I$(MUMPS)/include
     LDIRS += -L$(MUMPS)/lib
     ifdef PARMETIS
     	   LIBS  += -lparmetis -lmetis 
	   LDIRS += -L$(PARMETIS)/lib -L$(METIS)/lib
     endif
     ifdef SCOTCH
     	   LIBS  += -lesmumps -lscotch -lscotcherr 
           LDIRS += -L$(SCOTCH)/lib
     endif
endif

#
# Add parallel multigrid libraries if preprocessor set
#
#ifeq ($(USEPARMG),1)
#  LIBS += -lparmgrid
#endif

################################################################
#
# Section IV: Set up compiler and compiler flags
#
################################################################

CC  = $(SPC_MPICC)
FC  = $(SPC_MPIF90)
F90 = $(SPC_MPIF90)
GBS_EXEC = $(BINDIR)/gbs_$(PLAT)_$(COMP)
LDFLAGS = $(F90FLAGS)

#
# Intel + Dora + OpenMP needs dynamic linking
#
ifeq ($(PLAT), dora)
    ifeq ($(COMPTYPE), i) #intel
        ifeq ($(USEOPENMP),1)
            LDFLAGS += -dynamic
        endif
    endif
endif

################################################################
#
# Section V: Set up inclusion of modules and libraries during
#            compiling / linking phase
#
################################################################

#Flag for finding external modules in MODDIR
ifeq ($(COMPTYPE), i) #intel
     EXTMOD    = -module $(MODDIR)
endif
ifeq ($(COMPTYPE), g) #gnu
     EXTMOD    = -J $(MODDIR)
endif
ifeq ($(COMPTYPE), c) #cray
     EXTMOD    = -em -J $(MODDIR)
endif

#Flag for finding external libraries in LDIR
EXTLIBS   = $(LDIRS) -Wl,--start-group $(LIBS) -Wl,--end-group

#Flag for finding external include files 
EXTINC    = $(IDIRS)
