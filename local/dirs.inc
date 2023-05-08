#Local code, binaries, pputils library
PREFIX   = $(HOME)/gyacomo
SRCDIR   = $(PREFIX)/src
BINDIR   = $(PREFIX)/bin
OBJDIR   = $(PREFIX)/obj
LIBDIR   = $(PREFIX)/lib
MODDIR   = $(PREFIX)/mod

FMDIR      = $(HOME)/lib/FM
FUTILS_DIR = $(HOME)/lib/futils/src
# For local install :
HDF5_LIB   = /usr/local/hdf5-1.8.22/lib
FFTW3DIR   = $(HOME)/lib/fftw-3.3.10
# For marconi:
#  -comment HDF5_LIB and FFTW3DIR
#  -uncomment the following lines
#FFTW3DIR   = $(FFTW_HOME)#for Marconi
#LAPACKDIR  = $(HOME)/lib/lapack-3.10.0