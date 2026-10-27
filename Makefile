include local/dirs.inc
include local/make.inc

#------------------------------------------------------------------------------
# Compiler family detection for the generic "mpif90" targets
#
# mpif90 wraps gfortran on some machines and ifort on others, so these targets
# cannot assume either: -xHOST, -fast and -fpp are Intel spellings that gfortran
# rejects outright. Ask the compiler what it is instead of guessing.
#
# These use "=" (recursive), not ":=", so they are expanded inside the recipes,
# where the target-specific F90 is in scope.
#
# The named-compiler targets below (marconi*, daint, gfsp/gfdp/gdebug) state
# their compiler on purpose and keep their own explicit flags.
#------------------------------------------------------------------------------
F90_IS_GNU = $(findstring GNU,$(shell $(F90) --version 2>/dev/null | head -1))
# The sources are named .f90, so preprocessing must be requested explicitly.
FPP_AUTO   = $(if $(F90_IS_GNU),-cpp,-fpp)
OPT_AUTO   = $(if $(F90_IS_GNU),-O3 -ffree-line-length-0,-O3 -xHOST)
FAST_AUTO  = $(if $(F90_IS_GNU),-Ofast -ffree-line-length-0,-fast)
DBG_AUTO   = $(if $(F90_IS_GNU),-O0 -g -Wall -fcheck=all -fbacktrace -ffree-line-length-0,-C -g -traceback -warn all -check all -debug all -init=zero,snan -fpe0)
# Where the compiler writes (and looks for) the .mod files.
MODOUT_AUTO= $(if $(F90_IS_GNU),-J $(MODDIR),-module $(MODDIR))

# local/make.inc only sets EXTMOD when COMPTYPE is given. Without it the .mod
# files would land in whatever directory make runs from, i.e. the repository
# root, so fall back to MODDIR. An EXTMOD from COMPTYPE (e.g. the Cray one) or
# from a target below still takes precedence.
EXTMOD ?= $(MODOUT_AUTO)

# Standard version with optimized compilation and double precision
all: dp

# Double precision version
dp: F90 = mpif90
dp: F90FLAGS = $(OPT_AUTO)
dp: FPPFLAGS = $(FPP_AUTO)
dp: EXEC = $(BINDIR)/gyacomo23_dp
dp: dirs src/srcinfo.h
dp: compile

# Single precision version
sp: F90 = mpif90
sp: F90FLAGS = -DSINGLE_PRECISION $(OPT_AUTO)
sp: FPPFLAGS = $(FPP_AUTO)
sp: EXEC = $(BINDIR)/gyacomo23_sp
sp: dirs src/srcinfo.h
sp: compile

# Fast compilation
fast: F90 = mpif90
fast: F90FLAGS = $(FAST_AUTO)
fast: FPPFLAGS = $(FPP_AUTO)
fast: EXEC = $(BINDIR)/gyacomo23_fast
fast: dirs src/srcinfo.h
fast: compile

# test exec name
test: F90 = mpif90
test: F90FLAGS = -O1
test: FPPFLAGS = $(FPP_AUTO)
test: EXEC = $(BINDIR)/gyacomo23_test
test: dirs src/srcinfo.h
test: compile

# Debug version with all flags
debug: F90 = mpif90
debug: F90FLAGS = $(DBG_AUTO)
debug: FPPFLAGS = $(FPP_AUTO)
debug: EXEC = $(BINDIR)/gyacomo23_debug
debug: dirs src/srcinfo.h
debug: compile

# For compiling on marconi
marconi: F90 = mpiifort
marconi: F90FLAGS = -O3 -xHOST
marconi: FPPFLAGS = -fpp
marconi: EXEC = $(BINDIR)/gyacomo23_dp
marconi: dirs src/srcinfo.h
marconi: compile

# For compiling on marconi in single prec.
marconi_sp: F90 = mpiifort
marconi_sp: F90FLAGS = -DSINGLE_PRECISION -O3 -xHOST
marconi_sp: FPPFLAGS = -fpp
marconi_sp: EXEC = $(BINDIR)/gyacomo23_sp
marconi_sp: dirs src/srcinfo.h
marconi_sp: compile

# For compiling on marconi in single prec.
marconi_dbg: F90 = mpiifort
marconi_dbg: F90FLAGS = -DSINGLE_PRECISION -g -traceback -ftrapuv -warn all -debug all
marconi_dbg: FPPFLAGS = -fpp
marconi_dbg: EXEC = $(BINDIR)/gyacomo23_dbg
marconi_dbg: dirs src/srcinfo.h
marconi_dbg: compile

# For compiling on daint
daint: F90 = ftn
daint: F90FLAGS = -O3
daint: FPPFLAGS = -eZ
daint: EXEC = $(BINDIR)/gyacomo23_dp
daint: dirs src/srcinfo.h
daint: compile

# gfortran sp version, for compilation with gfortran
gfsp: F90 = mpif90
gfsp: F90FLAGS = -DSINGLE_PRECISION -O3 -std=legacy -ffree-line-length-0
gfsp: FPPFLAGS = -cpp
gfsp: EXTMOD   = -J $(MODDIR)
gfsp: EXEC = $(BINDIR)/gyacomo23_sp
gfsp: dirs src/srcinfo.h
gfsp: compile

# gfortran dp version, for compilation with gfortran
gfdp: F90 = mpif90
gfdp: F90FLAGS = -O3 -std=legacy -ffree-line-length-0
gfdp: FPPFLAGS = -cpp
gfdp: EXTMOD   = -J $(MODDIR)
gfdp: EXEC = $(BINDIR)/gyacomo23_dp
gfdp: dirs src/srcinfo.h
gfdp: compile

# gfortran version, compile with gfortran
gdebug: F90 = mpif90
gdebug: F90FLAGS = -O0 -g -Wall -fcheck=all -fbacktrace -std=legacy -ffree-line-length-0
gdebug: FPPFLAGS = -cpp
gdebug: EXTMOD   = -J $(MODDIR)
gdebug: EXEC = $(BINDIR)/gyacomo23_debug
gdebug: dirs src/srcinfo.h
gdebug: compile

# test SVD
test_svd: F90 = mpif90
test_svd: F90FLAGS = -DLAPACK -DTEST_SVD -O3
test_svd: FPPFLAGS = -fpp
test_svd: EXEC = $(BINDIR)/gyacomo23_test_svd
test_svd: dirs src/srcinfo.h
test_svd: compile

# to verify the path and other stuff
verif:
	@echo this is your PREFIX variable: $(PREFIX)
	@echo this is your HDF5_LIB variable: $(HDF5_LIB)
	@echo this is your FUTILS_DIR variable: $(FUTILS_DIR)
	@echo this is your FFTW3DIR variable: $(FFTW3DIR)
	@echo this is your FMDIR variable: $(FMDIR)

new_prob:
	@sh new_prob.sh

# subroutines
dirs:
	mkdir -p $(BINDIR)
	mkdir -p $(OBJDIR)
	mkdir -p $(MODDIR)
src/srcinfo.h:
	( cd src/srcinfo; $(MAKE);)
clean: cleanobj cleanmod cleanunit
	@rm -f src/srcinfo.h
	@rm -f src/srcinfo/srcinfo.h
cleanobj:
	@rm -f $(OBJDIR)/*o
cleanmod:
	@rm -f $(MODDIR)/*mod
	@rm -f *.mod
mvmod:
	mv *.mod mod/.

# attach git info
$(OBJDIR)/diagnose.o : src/srcinfo.h

# Main source dependencies
FOBJ=$(OBJDIR)/adaptive_dt_mod.o $(OBJDIR)/advance_field_mod.o $(OBJDIR)/array_mod.o $(OBJDIR)/auxval.o \
$(OBJDIR)/basic_mod.o $(OBJDIR)/coeff_mod.o $(OBJDIR)/closure_mod.o $(OBJDIR)/circular_mod.o \
$(OBJDIR)/collision_mod.o $(OBJDIR)/nonlinear_mod.o $(OBJDIR)/control.o \
$(OBJDIR)/diagnostics_mod.o $(OBJDIR)/endrun.o \
$(OBJDIR)/ExB_shear_flow_mod.o \
$(OBJDIR)/fields_mod.o $(OBJDIR)/fourier_mod.o $(OBJDIR)/geometry_mod.o \
$(OBJDIR)/ghosts_mod.o $(OBJDIR)/grid_mod.o \
$(OBJDIR)/initial_mod.o $(OBJDIR)/lag_interp_mod.o $(OBJDIR)/linear_growth_mod.o \
$(OBJDIR)/main.o \
$(OBJDIR)/memory.o $(OBJDIR)/miller_mod.o $(OBJDIR)/model_mod.o \
$(OBJDIR)/moments_eq_rhs_mod.o $(OBJDIR)/numerics_mod.o $(OBJDIR)/parallel_mod.o \
$(OBJDIR)/ppexit.o $(OBJDIR)/prec_const_mod.o \
$(OBJDIR)/processing_mod.o $(OBJDIR)/readinputs.o $(OBJDIR)/restarts_mod.o \
$(OBJDIR)/solve_EM_fields.o $(OBJDIR)/species_mod.o $(OBJDIR)/stepon.o $(OBJDIR)/tesend.o \
$(OBJDIR)/time_integration_mod.o $(OBJDIR)/utility_mod.o $(OBJDIR)/units_mod.o \
$(OBJDIR)/CLA_mod.o $(OBJDIR)/vmec_mod.o

# Add finally the definition of the Gyacomo directory
# (helps if the code needs to load a file as for the init Ricci)
 compile: F90FLAGS += -D__GYACDIR__=\"$(shell pwd)\"
# To compile the executable
 compile: $(FOBJ)
	$(F90) $(LDFLAGS) $(OBJDIR)/*.o $(EXTMOD) $(EXTINC) $(EXTLIBS) -o $(EXEC)

# Modules compilation
 $(OBJDIR)/adaptive_dt_mod.o : src/adaptive_dt_mod.f90 \
   $(OBJDIR)/prec_const_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/grid_mod.o \
	 $(OBJDIR)/fields_mod.o $(OBJDIR)/array_mod.o $(OBJDIR)/model_mod.o \
	 $(OBJDIR)/time_integration_mod.o $(OBJDIR)/closure_mod.o $(OBJDIR)/fourier_mod.o \
	 $(OBJDIR)/parallel_mod.o $(OBJDIR)/species_mod.o $(OBJDIR)/geometry_mod.o \
	 $(OBJDIR)/nonlinear_mod.o $(OBJDIR)/stepon.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/adaptive_dt_mod.f90 -o $@

 $(OBJDIR)/advance_field_mod.o : src/advance_field_mod.f90 \
   $(OBJDIR)/grid_mod.o $(OBJDIR)/array_mod.o $(OBJDIR)/initial_mod.o \
	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/time_integration_mod.o $(OBJDIR)/basic_mod.o \
	 $(OBJDIR)/fields_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/closure_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/advance_field_mod.f90 -o $@

 $(OBJDIR)/array_mod.o : src/array_mod.f90 \
	 $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/array_mod.f90 -o $@

 $(OBJDIR)/auxval.o : src/auxval.f90 \
 	 $(OBJDIR)/fourier_mod.o $(OBJDIR)/memory.o $(OBJDIR)/model_mod.o \
	 $(OBJDIR)/geometry_mod.o  $(OBJDIR)/grid_mod.o $(OBJDIR)/numerics_mod.o \
	 $(OBJDIR)/parallel_mod.o $(OBJDIR)/processing_mod.o $(OBJDIR)/CLA_mod.o \
	 $(OBJDIR)/ExB_shear_flow_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/auxval.f90 -o $@

 $(OBJDIR)/basic_mod.o : src/basic_mod.f90 \
 	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/parallel_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/basic_mod.f90 -o $@

 $(OBJDIR)/calculus_mod.o : src/calculus_mod.f90  \
 	 $(OBJDIR)/basic_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/grid_mod.o \
	 $(OBJDIR)/parallel_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/calculus_mod.f90 -o $@

 $(OBJDIR)/circular_mod.o : src/circular_mod.f90 \
 	 $(OBJDIR)/basic_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/lag_interp_mod.o\
	 $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/circular_mod.f90 -o $@

 $(OBJDIR)/coeff_mod.o : src/coeff_mod.f90 \
 	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/model_mod.o \
	 $(OBJDIR)/basic_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/coeff_mod.f90 -o $@

 $(OBJDIR)/closure_mod.o : src/closure_mod.f90 \
 	 $(OBJDIR)/model_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/grid_mod.o \
	 $(OBJDIR)/array_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/CLA_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/closure_mod.f90 -o $@

 $(OBJDIR)/collision_mod.o : src/collision_mod.f90 \
 	 $(OBJDIR)/array_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/cosolver_interface_mod.o\
	 $(OBJDIR)/fields_mod.o $(OBJDIR)/grid_mod.o $(OBJDIR)/model_mod.o \
	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/species_mod.o $(OBJDIR)/time_integration_mod.o \
	 $(OBJDIR)/utility_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/collision_mod.f90 -o $@

 $(OBJDIR)/control.o : src/control.f90 \
	 $(OBJDIR)/auxval.o $(OBJDIR)/geometry_mod.o $(OBJDIR)/prec_const_mod.o \
     $(OBJDIR)/basic_mod.o $(OBJDIR)/ppexit.o $(OBJDIR)/readinputs.o $(OBJDIR)/tesend.o \
     $(OBJDIR)/diagnostics_mod.o $(OBJDIR)/linear_growth_mod.o $(OBJDIR)/stepon.o \
     $(OBJDIR)/adaptive_dt_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/control.f90 -o $@

 $(OBJDIR)/cosolver_interface_mod.o : src/cosolver_interface_mod.f90 \
   $(OBJDIR)/grid_mod.o $(OBJDIR)/array_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/species_mod.o\
	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/closure_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/cosolver_interface_mod.f90 -o $@

 $(OBJDIR)/diagnostics_mod.o : src/diagnostics_mod.f90 \
 	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/processing_mod.o $(OBJDIR)/array_mod.o \
     $(OBJDIR)/basic_mod.o $(OBJDIR)/fields_mod.o \
	 $(OBJDIR)/grid_mod.o $(OBJDIR)/initial_mod.o $(OBJDIR)/model_mod.o \
	 $(OBJDIR)/time_integration_mod.o $(OBJDIR)/parallel_mod.o $(OBJDIR)/units_mod.o \
	 $(OBJDIR)/linear_growth_mod.o $(OBJDIR)/adaptive_dt_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/diagnostics_mod.f90 -o $@

 $(OBJDIR)/endrun.o : src/endrun.f90 \
 	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/basic_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/endrun.f90 -o $@

 $(OBJDIR)/ExB_shear_flow_mod.o : src/ExB_shear_flow_mod.f90 \
	 $(OBJDIR)/basic_mod.o $(OBJDIR)/grid_mod.o $(OBJDIR)/prec_const_mod.o\
	 $(OBJDIR)/geometry_mod.o $(OBJDIR)/numerics_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/ExB_shear_flow_mod.f90 -o $@

 $(OBJDIR)/fields_mod.o : src/fields_mod.f90 \
   $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/fields_mod.f90 -o $@

 $(OBJDIR)/fourier_mod.o : src/fourier_mod.f90 \
 	 $(OBJDIR)/basic_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/utility_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/fourier_mod.f90 -o $@

 $(OBJDIR)/geometry_mod.o : src/geometry_mod.f90 \
 	 $(OBJDIR)/array_mod.o $(OBJDIR)/calculus_mod.o $(OBJDIR)/circular_mod.o \
	 $(OBJDIR)/miller_mod.o $(OBJDIR)/grid_mod.o $(OBJDIR)/model_mod.o \
	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/utility_mod.o $(OBJDIR)/vmec_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/geometry_mod.f90 -o $@

 $(OBJDIR)/vmec_mod.o : src/vmec_mod.f90 \
 	 $(OBJDIR)/basic_mod.o $(OBJDIR)/parallel_mod.o $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/vmec_mod.f90 -o $@

 $(OBJDIR)/ghosts_mod.o : src/ghosts_mod.f90 \
 	 $(OBJDIR)/basic_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/grid_mod.o\
	 $(OBJDIR)/geometry_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/time_integration_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/ghosts_mod.f90 -o $@

 $(OBJDIR)/grid_mod.o : src/grid_mod.f90 \
 	 $(OBJDIR)/basic_mod.o $(OBJDIR)/fourier_mod.o $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/grid_mod.f90 -o $@

 $(OBJDIR)/initial_mod.o : src/initial_mod.f90 \
 	 $(OBJDIR)/array_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/fields_mod.o \
     $(OBJDIR)/model_mod.o $(OBJDIR)/numerics_mod.o $(OBJDIR)/ExB_shear_flow_mod.o \
	 $(OBJDIR)/solve_EM_fields.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/ghosts_mod.o \
	 $(OBJDIR)/grid_mod.o $(OBJDIR)/collision_mod.o $(OBJDIR)/closure_mod.o \
	 $(OBJDIR)/nonlinear_mod.o \
	 $(OBJDIR)/restarts_mod.o $(OBJDIR)/time_integration_mod.o $(OBJDIR)/utility_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/initial_mod.f90 -o $@

 $(OBJDIR)/lag_interp_mod.o : src/lag_interp_mod.f90 \
 	 $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/lag_interp_mod.f90 -o $@

 $(OBJDIR)/linear_growth_mod.o : src/linear_growth_mod.f90 \
	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/grid_mod.o \
	 $(OBJDIR)/fields_mod.o $(OBJDIR)/parallel_mod.o $(OBJDIR)/model_mod.o \
	 $(OBJDIR)/time_integration_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/linear_growth_mod.f90 -o $@

 $(OBJDIR)/main.o : src/main.f90 $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/main.f90 -o $@

 $(OBJDIR)/memory.o : src/memory.f90 $ \
 	 $(OBJDIR)/array_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/collision_mod.o\
     $(OBJDIR)/fields_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/time_integration_mod.o \
	 $(OBJDIR)/grid_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/memory.f90 -o $@

 $(OBJDIR)/miller_mod.o : src/miller_mod.f90 \
 	 $(OBJDIR)/basic_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/lag_interp_mod.o\
	 $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/miller_mod.f90 -o $@

 $(OBJDIR)/model_mod.o : src/model_mod.f90 \
 	 $(OBJDIR)/grid_mod.o $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/model_mod.f90 -o $@

 $(OBJDIR)/moments_eq_rhs_mod.o : src/moments_eq_rhs_mod.f90 \
 	 $(OBJDIR)/array_mod.o $(OBJDIR)/calculus_mod.o $(OBJDIR)/fields_mod.o \
	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/grid_mod.o $(OBJDIR)/model_mod.o \
	 $(OBJDIR)/time_integration_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/moments_eq_rhs_mod.f90 -o $@

 $(OBJDIR)/nonlinear_mod.o : src/nonlinear_mod.f90 \
 	 $(OBJDIR)/array_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/ExB_shear_flow_mod.o \
	 $(OBJDIR)/fourier_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/grid_mod.o $(OBJDIR)/model_mod.o\
	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/time_integration_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/nonlinear_mod.f90 -o $@

 $(OBJDIR)/numerics_mod.o : src/numerics_mod.f90 \
 	 $(OBJDIR)/prec_const_mod.o  $(OBJDIR)/basic_mod.o $(OBJDIR)/coeff_mod.o \
	 $(OBJDIR)/utility_mod.o $(OBJDIR)/species_mod.o $(OBJDIR)/geometry_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/numerics_mod.f90 -o $@

 $(OBJDIR)/parallel_mod.o : src/parallel_mod.f90 \
 	 $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/parallel_mod.f90 -o $@

 $(OBJDIR)/ppexit.o : src/ppexit.f90 \
 	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/coeff_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/ppexit.f90 -o $@

 $(OBJDIR)/prec_const_mod.o : src/prec_const_mod.f90
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/prec_const_mod.f90 -o $@

 $(OBJDIR)/processing_mod.o : src/processing_mod.f90 \
 	 $(OBJDIR)/array_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/grid_mod.o \
   $(OBJDIR)/fields_mod.o $(OBJDIR)/basic_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/processing_mod.f90 -o $@

 $(OBJDIR)/readinputs.o : src/readinputs.f90  \
   $(OBJDIR)/diagnostics_mod.o $(OBJDIR)/initial_mod.o $(OBJDIR)/model_mod.o  \
	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/grid_mod.o $(OBJDIR)/time_integration_mod.o \
  	 $(OBJDIR)/units_mod.o $(OBJDIR)/linear_growth_mod.o $(OBJDIR)/adaptive_dt_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/readinputs.f90 -o $@

 $(OBJDIR)/restarts_mod.o : src/restarts_mod.f90   \
	 $(OBJDIR)/grid_mod.o $(OBJDIR)/time_integration_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/restarts_mod.f90 -o $@

 $(OBJDIR)/solve_EM_fields.o : src/solve_EM_fields.f90 \
   $(OBJDIR)/array_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/grid_mod.o \
	 $(OBJDIR)/ghosts_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/array_mod.o \
	 $(OBJDIR)/time_integration_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/parallel_mod.o \
	 $(OBJDIR)/processing_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/solve_EM_fields.f90 -o $@

 $(OBJDIR)/species_mod.o : src/species_mod.f90 \
   $(OBJDIR)/basic_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/species_mod.f90 -o $@

 $(OBJDIR)/stepon.o : src/stepon.f90 \
 	 $(OBJDIR)/initial_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/advance_field_mod.o \
     $(OBJDIR)/basic_mod.o $(OBJDIR)/nonlinear_mod.o $(OBJDIR)/grid_mod.o \
	 $(OBJDIR)/array_mod.o $(OBJDIR)/numerics_mod.o $(OBJDIR)/fields_mod.o \
	 $(OBJDIR)/ghosts_mod.o $(OBJDIR)/moments_eq_rhs_mod.o $(OBJDIR)/solve_EM_fields.o\
	 $(OBJDIR)/utility_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/time_integration_mod.o \
	 $(OBJDIR)/CLA_mod.o $(OBJDIR)/ExB_shear_flow_mod.o $(OBJDIR)/closure_mod.o \
	 $(OBJDIR)/collision_mod.o $(OBJDIR)/processing_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/stepon.f90 -o $@

 $(OBJDIR)/tesend.o : src/tesend.f90 \
 	 $(OBJDIR)/basic_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/linear_growth_mod.o \
	 $(OBJDIR)/adaptive_dt_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/tesend.f90 -o $@

 $(OBJDIR)/time_integration_mod.o : src/time_integration_mod.f90 \
   $(OBJDIR)/basic_mod.o  $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/time_integration_mod.f90 -o $@

 $(OBJDIR)/utility_mod.o : src/utility_mod.f90  \
   $(OBJDIR)/basic_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/time_integration_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/utility_mod.f90 -o $@

 $(OBJDIR)/units_mod.o : src/units_mod.f90  \
	$(OBJDIR)/basic_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/geometry_mod.o 
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/units_mod.f90 -o $@

 $(OBJDIR)/CLA_mod.o : src/CLA_mod.f90  \
   $(OBJDIR)/basic_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/grid_mod.o \
   $(OBJDIR)/time_integration_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/CLA_mod.f90 -o $@

#------------------------------------------------------------------------------
# Unit tests -- "make check"
#
# The tests are built into their own obj/ and mod/ subtrees through a recursive
# make that overrides OBJDIR and MODDIR. That isolation is deliberate: the
# normal targets write their .mod files to mod/ with the flags of the last
# build, and a test build must not depend on, or clobber, that.
#
# Each test program links against every object except main.o, so a test can USE
# any module of the code without this Makefile needing a per-test dependency
# list. The price is one link per test, which is cheap next to the compile.
#
# Runtime checks are on: an out of bounds access in the code under test should
# fail the test rather than go unnoticed.
#
# Set TESTRUNNER if your MPI needs a launcher even for a single rank, e.g.
#   make check TESTRUNNER="mpirun -np 1"
#------------------------------------------------------------------------------
UNITSRCDIR = $(PREFIX)/unit
UNITBINDIR = $(BINDIR)/unit
UNITOBJDIR = $(OBJDIR)/tests
TESTRUNNER =

CHK_AUTO   = $(if $(F90_IS_GNU),-O1 -g -fcheck=all -fbacktrace -ffree-line-length-0,-O1 -g -traceback -check all)

# One executable per test program in unit/ (test_harness_mod is the shared
# assertion helper, not a program).
UNIT_PROGS = test_prec_const test_utility test_lag_interp test_calculus \
             test_numerics test_coeff test_grid test_time_integration \
             test_ExB_shear_flow test_linear_growth test_bpar test_staggered_z \
             test_adaptive_dt test_vmec test_adiab_e test_nonlinear_closure
UNIT_EXECS = $(addprefix $(UNITBINDIR)/,$(UNIT_PROGS))
UNIT_OBJS  = $(addprefix $(UNITOBJDIR)/,$(addsuffix .o,$(UNIT_PROGS)))

# allocate_array does not deallocate first, so time_integration_readinputs can
# only run once per process. The scheme is a command line argument instead and
# the binary is run once per scheme.
TI_SCHEMES = EE RK2 RK3 SSP_RK3 ODE23 RK4 DOPRI5

# Same story for the ExB suite: Setup_ExB_shear_flow allocates without
# deallocating, so the linearity is a command line argument and the binary is
# run once per mode.
EXB_MODES = linear nonlinear off

# And for the delta B_par suite: readinputs runs once per process, so the
# flag/beta combination is a command line argument and the binary is run once
# per combination.
BPAR_MODES = on off nobeta

# And for the adaptive time step suite: one deck (linear hyperdiffusion or a
# nonlinear FFT setup) per process, so the mode is a command line argument.
ADT_MODES = hyperdiff cfl

# And for the nonlinear Laguerre closure suite: one closure per process, so the
# closure is a command line argument.
NLC_MODES = full_sum anti_laguerre_aliasing anti_laguerre_aliasing_sym truncation truncation_nmax1

check:
	@$(MAKE) --no-print-directory run-unit-tests \
	    OBJDIR=$(PREFIX)/obj/unit MODDIR=$(PREFIX)/mod/unit

run-unit-tests: F90      = mpif90
run-unit-tests: F90FLAGS = $(CHK_AUTO) -D__GYACDIR__=\"$(PREFIX)\"
run-unit-tests: FPPFLAGS = $(FPP_AUTO)
run-unit-tests: EXTMOD   = $(MODOUT_AUTO)
run-unit-tests: dirs unitdirs src/srcinfo.h $(UNIT_EXECS)
	@echo ""
	@echo "== running unit tests =="
	@echo ""
	@fail=0; \
	for t in $(filter-out $(UNITBINDIR)/test_time_integration $(UNITBINDIR)/test_ExB_shear_flow $(UNITBINDIR)/test_bpar $(UNITBINDIR)/test_adaptive_dt $(UNITBINDIR)/test_nonlinear_closure,$(UNIT_EXECS)); do \
	  $(TESTRUNNER) $$t || fail=`expr $$fail + 1`; \
	done; \
	for s in $(TI_SCHEMES); do \
	  $(TESTRUNNER) $(UNITBINDIR)/test_time_integration $$s || fail=`expr $$fail + 1`; \
	done; \
	for m in $(EXB_MODES); do \
	  $(TESTRUNNER) $(UNITBINDIR)/test_ExB_shear_flow $$m || fail=`expr $$fail + 1`; \
	done; \
	for m in $(BPAR_MODES); do \
	  $(TESTRUNNER) $(UNITBINDIR)/test_bpar $$m || fail=`expr $$fail + 1`; \
	done; \
	for m in $(ADT_MODES); do \
	  $(TESTRUNNER) $(UNITBINDIR)/test_adaptive_dt $$m || fail=`expr $$fail + 1`; \
	done; \
	for m in $(NLC_MODES); do \
	  $(TESTRUNNER) $(UNITBINDIR)/test_nonlinear_closure $$m || fail=`expr $$fail + 1`; \
	done; \
	if [ $$fail -eq 0 ]; then \
	  echo "== all unit tests passed =="; \
	else \
	  echo "== $$fail unit test suite(s) FAILED =="; exit 1; \
	fi

unitdirs:
	mkdir -p $(UNITBINDIR)
	mkdir -p $(UNITOBJDIR)

cleanunit:
	@rm -rf $(PREFIX)/obj/unit $(PREFIX)/mod/unit $(BINDIR)/unit

# Shared assertion helpers
$(UNITOBJDIR)/test_harness_mod.o : $(UNITSRCDIR)/test_harness_mod.f90 \
   $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) $< -o $@

# Every test program may reach any module, so all of them are prerequisites.
$(UNITOBJDIR)/%.o : $(UNITSRCDIR)/%.f90 $(UNITOBJDIR)/test_harness_mod.o $(FOBJ)
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) $< -o $@

# The object list is globbed at recipe time so that objects built as a side
# effect of $(FOBJ) (calculus_mod, cosolver_interface_mod) are picked up too.
# main.o is dropped because it carries the PROGRAM unit of the solver.
$(UNITBINDIR)/% : $(UNITOBJDIR)/%.o
	$(F90) $(LDFLAGS) $< $(UNITOBJDIR)/test_harness_mod.o \
	  $(filter-out $(OBJDIR)/main.o,$(wildcard $(OBJDIR)/*.o)) \
	  $(EXTMOD) $(EXTINC) $(EXTLIBS) -o $@

.PHONY: check run-unit-tests unitdirs cleanunit
# Without this, make treats the test objects as intermediates of the chained
# pattern rules above and deletes them, recompiling every test on every run.
.SECONDARY: $(UNIT_OBJS)
