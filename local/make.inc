################################################################
#
# Section I: Compiler options
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
# Section II: Libraries and where to find them
#
################################################################

IDIRS := -I$(INTEL_LOCAL)/include/$(OPTLEVEL)

LIBS  := -lfutils -lhdf5_fortran -lhdf5 -lz -ldl -lpthread
LDIRS := -L$(INTEL_LOCAL)/lib/$(OPTLEVEL) -L$(HDF5)/lib

# Add Multiple-Precision Library
LIBS += -lfm


ifdef FFTW3DIR
      LIBS  += -lfftw3 -lfftw3_mpi
      LDIRS += -L$(FFTW3DIR)/lib64
      IDIRS += -I$(FFTW3DIR)/include
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
