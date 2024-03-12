#Local code, binaries, pputils library
PREFIX   = $(HOME)/gyacomo
SRCDIR   = $(PREFIX)/src
BINDIR   = $(PREFIX)/bin
OBJDIR   = $(PREFIX)/obj
LIBDIR   = $(HOME)/lib
MODDIR   = $(PREFIX)/mod

#OS environment: Linux, Marconi or MacOS
ENVTYPE = Linux

ifeq ($(ENVTYPE), Linux) #Linux env.
	FMDIR      = $(LIBDIR)/FM
	FUTILS_DIR = $(LIBDIR)/futils/src
	FFTW3DIR   = $(LIBDIR)/fftw-3.3.10
	LAPACKDIR  = $(HOME)/lib/lapack-3.10.0
	HDF5_LIB   = $(HDF5_ROOT)/lib
endif
ifeq ($(ENVTYPE), Marconi) #Marconi env.
	FMDIR      = $(LIBDIR)/FM
	FUTILS_DIR = $(LIBDIR)/futils/src
	FFTW3DIR   = $(FFTW_HOME)
endif
ifeq ($(ENVTYPE), MacOS) #MacOS
	FMDIR      = $(LIBDIR)/FM
	FUTILS_DIR = $(LIBDIR)/futils/src
	FFTW3DIR   = $(LIBDIR)/fftw-3.3.10
	LAPACKDIR = /opt/homebrew/opt/lapack
	HDF5_LIB   = /usr/local/hdf5-1.8.22/lib
endif
