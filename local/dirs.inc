#Local code, binaries, pputils library
PREFIX   = $(HOME)/gyacomo
SRCDIR   = $(PREFIX)/src
BINDIR   = $(PREFIX)/bin
OBJDIR   = $(PREFIX)/obj
LIBDIR   = $(HOME)/lib
#LIBDIR   = $(HOME)/lib_gnu 
MODDIR   = $(PREFIX)/mod

FMDIR      = $(LIBDIR)/FM
FUTILS_DIR = $(LIBDIR)/futils/src
# For local install :
HDF5_LIB   = $(HDF5_ROOT)/lib
#HDF5_LIB   = /usr/local/hdf5-1.8.22/lib
FFTW3DIR   = $(LIBDIR)/fftw-3.3.10
#LAPACKDIR  = $(LIBDIR)/lapack-3.10.0
# For marconi (use of modules):
#  -comment HDF5_LIB, FFTW3DIR and LAPACKDIR
#  -uncomment the following lines if the module are loaded
#FFTW3DIR   = $(FFTW_HOME)#mandatory
#LAPACKDIR  = $(LAPACK_HOME)#optional
