# Variables for main path and libraries
PREFIX   = $(HOME)/gyacomo # this must be the path to your /gyaomo folder
LIBDIR   = $(HOME)/lib     # this is the path to the library folder

# OS environment: Linux, Marconi or MacOS (use Linux per default)
ENVTYPE = Linux

# Auxiliary variables
SRCDIR   = $(PREFIX)/src   # this is the path to the /gyacomo/src
BINDIR   = $(PREFIX)/bin   # this is the path to the /gyacomo/bin
OBJDIR   = $(PREFIX)/obj   # this is the path to the /gyacomo/obj
MODDIR   = $(PREFIX)/mod   # this is the path to the /gyacomo/src

ifeq ($(ENVTYPE), Linux) #Linux env.
	FMDIR      = $(LIBDIR)/FM
	FUTILS_DIR = $(LIBDIR)/futils/src
	FFTW3DIR   = $(LIBDIR)/fftw-3.3.10
	#LAPACKDIR  = $(LIBDIR)/lapack-3.10.0
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
