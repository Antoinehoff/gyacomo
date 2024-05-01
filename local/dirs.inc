# Variables for main path to /gyacomo
PREFIX   = $(shell pwd)

# Auxiliary variables (for SPCPC)
LIBDIR   = $(PREFIX)/lib
SRCDIR   = $(PREFIX)/src
BINDIR   = $(PREFIX)/bin
OBJDIR   = $(PREFIX)/obj
MODDIR   = $(PREFIX)/mod
FMDIR    = $(LIBDIR)/FM
FUTILS_DIR = $(LIBDIR)/futils/src
FFTW3DIR   = $(LIBDIR)/fftw-3.3.10
#LAPACKDIR  = $(LIBDIR)/lapack-3.10.0
HDF5_LIB   = $(HDF5_ROOT)/lib

# Overwrite path according to the OS environment (Linux, MacOS or Marconi)
ENVTYPE = Linux

ifeq ($(ENVTYPE), MacOS) #MacOS
	LAPACKDIR = /opt/homebrew/opt/lapack
	HDF5_LIB   = /usr/local/hdf5-1.8.22/lib
endif
ifeq ($(ENVTYPE), Marconi) #Marconi env.
	FFTW3DIR   = $(FFTW_HOME)
endif