include local/dirs.inc
include local/make.inc

# Standard version with optimized compilation and double precision
all: dp

# Double precision version
dp: F90 = mpif90
dp: F90FLAGS = -O3 -xHOST
dp: EXEC = $(BINDIR)/gyacomo23_dp
dp: dirs src/srcinfo.h
dp: compile

# Single precision version
sp: F90 = mpif90
sp: F90FLAGS = -DSINGLE_PRECISION -O3 -xHOST
sp: EXEC = $(BINDIR)/gyacomo23_sp
sp: dirs src/srcinfo.h
sp: compile

# Fast compilation
fast: F90 = mpif90
fast: F90FLAGS = -fast
fast: EXEC = $(BINDIR)/gyacomo23_fast
fast: dirs src/srcinfo.h
fast: compile

# Debug version with all flags
debug: F90 = mpif90
debug: F90FLAGS = -C -g -traceback -warn all -check all -debug all -init=zero,snan -fpe0
debug: EXEC = $(BINDIR)/gyacomo23_debug
debug: dirs src/srcinfo.h
debug: compile

# For compiling on marconi
marconi: F90 = mpiifort
marconi: F90FLAGS = -O3 -xHOST
marconi: EXEC = $(BINDIR)/gyacomo23_dp
marconi: dirs src/srcinfo.h
marconi: compile

# For compiling on marconi in single prec.
marconi_sp: F90 = mpiifort
marconi_sp: F90FLAGS = -O3 -xHOST
marconi_sp: F90FLAGS = -DSINGLE_PRECISION -O3 -xHOST
marconi_sp: EXEC = $(BINDIR)/gyacomo23_sp
marconi_sp: dirs src/srcinfo.h
marconi_sp: compile

# For compiling on marconi in single prec.
marconi_dbg: F90 = mpiifort
marconi_dbg: F90FLAGS = -DSINGLE_PRECISION -g -traceback -ftrapuv -warn all -debug all
marconi_dbg: EXEC = $(BINDIR)/gyacomo23_dbg
marconi_dbg: dirs src/srcinfo.h
marconi_dbg: compile

# For compiling on daint
daint: F90 = ftn
daint: F90FLAGS = -O3
daint: EXEC = $(BINDIR)/gyacomo23_dp
daint: dirs src/srcinfo.h
daint: compile

# gfortran sp version, for compilation with gfortran
gfsp: F90 = mpif90
gfsp: F90FLAGS = -DSINGLE_PRECISION -O3 -std=legacy -ffree-line-length-0
gfsp: EXTMOD   = -J $(MODDIR)
gfsp: EXEC = $(BINDIR)/gyacomo23_sp
gfsp: dirs src/srcinfo.h
gfsp: compile

# gfortran dp version, for compilation with gfortran
gfdp: F90 = mpif90
gfdp: F90FLAGS = -O3 -std=legacy -ffree-line-length-0
gfdp: EXTMOD   = -J $(MODDIR)
gfdp: EXEC = $(BINDIR)/gyacomo23_dp
gfdp: dirs src/srcinfo.h
gfdp: compile

# gfortran version, compile with gfortran
gdebug: F90 = mpif90
gdebug: F90FLAGS = -C -g -std=legacy -ffree-line-length-0
gdebug: EXTMOD   = -J $(MODDIR)
gdebug: EXEC = $(BINDIR)/gyacomo23_debug
gdebug: dirs src/srcinfo.h
gdebug: compile

# test SVD
test_svd: F90 = mpif90
# test_svd: F90FLAGS = -DTEST_SVD -g -traceback -ftrapuv -warn all -debug all
test_svd: F90FLAGS = -DTEST_SVD -O3
test_svd: EXEC = $(BINDIR)/gyacomo23_test_svd
test_svd: dirs src/srcinfo.h
test_svd: compile
# subroutines
dirs:
	mkdir -p $(BINDIR)
	mkdir -p $(OBJDIR)
	mkdir -p $(MODDIR)
src/srcinfo.h:
	( cd src/srcinfo; $(MAKE);)
clean: cleanobj cleanmod
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
FOBJ=$(OBJDIR)/advance_field_mod.o $(OBJDIR)/array_mod.o $(OBJDIR)/auxval.o \
$(OBJDIR)/basic_mod.o $(OBJDIR)/coeff_mod.o $(OBJDIR)/closure_mod.o $(OBJDIR)/circular_mod.o \
$(OBJDIR)/collision_mod.o $(OBJDIR)/nonlinear_mod.o $(OBJDIR)/control.o \
$(OBJDIR)/diagnostics_mod.o $(OBJDIR)/endrun.o \
$(OBJDIR)/ExB_shear_flow_mod.o \
$(OBJDIR)/fields_mod.o $(OBJDIR)/fourier_mod.o $(OBJDIR)/geometry_mod.o \
$(OBJDIR)/ghosts_mod.o $(OBJDIR)/grid_mod.o \
$(OBJDIR)/initial_mod.o $(OBJDIR)/lag_interp_mod.o $(OBJDIR)/main.o \
$(OBJDIR)/memory.o $(OBJDIR)/miller_mod.o $(OBJDIR)/model_mod.o \
$(OBJDIR)/moments_eq_rhs_mod.o $(OBJDIR)/numerics_mod.o $(OBJDIR)/parallel_mod.o \
$(OBJDIR)/ppexit.o $(OBJDIR)/prec_const_mod.o \
$(OBJDIR)/processing_mod.o $(OBJDIR)/readinputs.o $(OBJDIR)/restarts_mod.o \
$(OBJDIR)/solve_EM_fields.o $(OBJDIR)/species_mod.o $(OBJDIR)/stepon.o $(OBJDIR)/tesend.o \
$(OBJDIR)/time_integration_mod.o $(OBJDIR)/utility_mod.o $(OBJDIR)/CLA_mod.o

# Add finally the definition of the Gyacomo directory
# (helps if the code needs to load a file as for the init Ricci)
 compile: F90FLAGS += -D__GYACDIR__=\"$(shell pwd)\"
# To compile the executable
 compile: $(FOBJ)
	$(F90) $(LDFLAGS) $(OBJDIR)/*.o $(EXTMOD) $(EXTINC) $(EXTLIBS) -o $(EXEC)

# Modules compilation
 $(OBJDIR)/advance_field_mod.o : src/advance_field_mod.F90 \
   $(OBJDIR)/grid_mod.o $(OBJDIR)/array_mod.o $(OBJDIR)/initial_mod.o \
	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/time_integration_mod.o $(OBJDIR)/basic_mod.o \
	 $(OBJDIR)/fields_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/closure_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/advance_field_mod.F90 -o $@

 $(OBJDIR)/array_mod.o : src/array_mod.F90 \
	 $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/array_mod.F90 -o $@

 $(OBJDIR)/auxval.o : src/auxval.F90 \
 	 $(OBJDIR)/fourier_mod.o $(OBJDIR)/memory.o $(OBJDIR)/model_mod.o \
	 $(OBJDIR)/geometry_mod.o  $(OBJDIR)/grid_mod.o $(OBJDIR)/numerics_mod.o \
	 $(OBJDIR)/parallel_mod.o $(OBJDIR)/processing_mod.o $(OBJDIR)/CLA_mod.o \
	 $(OBJDIR)/ExB_shear_flow_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/auxval.F90 -o $@

 $(OBJDIR)/basic_mod.o : src/basic_mod.F90 \
 	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/parallel_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/basic_mod.F90 -o $@

 $(OBJDIR)/calculus_mod.o : src/calculus_mod.F90  \
 	 $(OBJDIR)/basic_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/grid_mod.o \
	 $(OBJDIR)/parallel_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/calculus_mod.F90 -o $@

 $(OBJDIR)/circular_mod.o : src/circular_mod.F90 \
 	 $(OBJDIR)/basic_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/lag_interp_mod.o\
	 $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/circular_mod.F90 -o $@

 $(OBJDIR)/coeff_mod.o : src/coeff_mod.F90 \
 	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/model_mod.o \
	 $(OBJDIR)/basic_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/coeff_mod.F90 -o $@

 $(OBJDIR)/closure_mod.o : src/closure_mod.F90 \
 	 $(OBJDIR)/model_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/grid_mod.o \
	 $(OBJDIR)/array_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/CLA_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/closure_mod.F90 -o $@

 $(OBJDIR)/collision_mod.o : src/collision_mod.F90 \
 	 $(OBJDIR)/array_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/cosolver_interface_mod.o\
	 $(OBJDIR)/fields_mod.o $(OBJDIR)/grid_mod.o $(OBJDIR)/model_mod.o \
	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/species_mod.o $(OBJDIR)/time_integration_mod.o \
	 $(OBJDIR)/utility_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/collision_mod.F90 -o $@

 $(OBJDIR)/control.o : src/control.F90 \
	 $(OBJDIR)/auxval.o $(OBJDIR)/geometry_mod.o $(OBJDIR)/prec_const_mod.o \
     $(OBJDIR)/basic_mod.o $(OBJDIR)/ppexit.o $(OBJDIR)/readinputs.o $(OBJDIR)/tesend.o \
     $(OBJDIR)/diagnostics_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/control.F90 -o $@

 $(OBJDIR)/cosolver_interface_mod.o : src/cosolver_interface_mod.F90 \
   $(OBJDIR)/grid_mod.o $(OBJDIR)/array_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/species_mod.o\
	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/closure_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/cosolver_interface_mod.F90 -o $@

 $(OBJDIR)/diagnostics_mod.o : src/diagnostics_mod.F90 \
 	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/processing_mod.o $(OBJDIR)/array_mod.o \
     $(OBJDIR)/basic_mod.o $(OBJDIR)/fields_mod.o \
	 $(OBJDIR)/grid_mod.o $(OBJDIR)/initial_mod.o $(OBJDIR)/model_mod.o \
	 $(OBJDIR)/time_integration_mod.o $(OBJDIR)/parallel_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/diagnostics_mod.F90 -o $@

 $(OBJDIR)/endrun.o : src/endrun.F90 \
 	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/basic_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/endrun.F90 -o $@

 $(OBJDIR)/ExB_shear_flow_mod.o : src/ExB_shear_flow_mod.F90 \
	 $(OBJDIR)/basic_mod.o $(OBJDIR)/grid_mod.o $(OBJDIR)/prec_const_mod.o\
	 $(OBJDIR)/geometry_mod.o $(OBJDIR)/numerics_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/ExB_shear_flow_mod.F90 -o $@

 $(OBJDIR)/fields_mod.o : src/fields_mod.F90 \
   $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/fields_mod.F90 -o $@

 $(OBJDIR)/fourier_mod.o : src/fourier_mod.F90 \
 	 $(OBJDIR)/basic_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/utility_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/fourier_mod.F90 -o $@

 $(OBJDIR)/geometry_mod.o : src/geometry_mod.F90 \
 	 $(OBJDIR)/array_mod.o $(OBJDIR)/calculus_mod.o $(OBJDIR)/circular_mod.o \
	 $(OBJDIR)/miller_mod.o $(OBJDIR)/grid_mod.o $(OBJDIR)/model_mod.o \
	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/utility_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/geometry_mod.F90 -o $@

 $(OBJDIR)/ghosts_mod.o : src/ghosts_mod.F90 \
 	 $(OBJDIR)/basic_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/grid_mod.o\
	 $(OBJDIR)/geometry_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/time_integration_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/ghosts_mod.F90 -o $@

 $(OBJDIR)/grid_mod.o : src/grid_mod.F90 \
 	 $(OBJDIR)/basic_mod.o $(OBJDIR)/fourier_mod.o $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/grid_mod.F90 -o $@

 $(OBJDIR)/initial_mod.o : src/initial_mod.F90 \
 	 $(OBJDIR)/array_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/fields_mod.o \
     $(OBJDIR)/model_mod.o $(OBJDIR)/numerics_mod.o \
	 $(OBJDIR)/solve_EM_fields.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/ghosts_mod.o \
	 $(OBJDIR)/grid_mod.o $(OBJDIR)/collision_mod.o $(OBJDIR)/closure_mod.o \
	 $(OBJDIR)/nonlinear_mod.o \
	 $(OBJDIR)/restarts_mod.o $(OBJDIR)/time_integration_mod.o $(OBJDIR)/utility_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/initial_mod.F90 -o $@

 $(OBJDIR)/lag_interp_mod.o : src/lag_interp_mod.F90 \
 	 $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/lag_interp_mod.F90 -o $@

 $(OBJDIR)/main.o : src/main.F90 $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/main.F90 -o $@

 $(OBJDIR)/memory.o : src/memory.F90 $ \
 	 $(OBJDIR)/array_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/collision_mod.o\
     $(OBJDIR)/fields_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/time_integration_mod.o \
	 $(OBJDIR)/grid_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/memory.F90 -o $@

 $(OBJDIR)/miller_mod.o : src/miller_mod.F90 \
 	 $(OBJDIR)/basic_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/lag_interp_mod.o\
	 $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/miller_mod.F90 -o $@

 $(OBJDIR)/model_mod.o : src/model_mod.F90 \
 	 $(OBJDIR)/grid_mod.o $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/model_mod.F90 -o $@

 $(OBJDIR)/moments_eq_rhs_mod.o : src/moments_eq_rhs_mod.F90 \
 	 $(OBJDIR)/array_mod.o $(OBJDIR)/calculus_mod.o $(OBJDIR)/fields_mod.o \
	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/grid_mod.o $(OBJDIR)/model_mod.o \
	 $(OBJDIR)/time_integration_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/moments_eq_rhs_mod.F90 -o $@

 $(OBJDIR)/nonlinear_mod.o : src/nonlinear_mod.F90 \
 	 $(OBJDIR)/array_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/ExB_shear_flow_mod.o \
	 $(OBJDIR)/fourier_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/grid_mod.o $(OBJDIR)/model_mod.o\
	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/time_integration_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/nonlinear_mod.F90 -o $@

 $(OBJDIR)/numerics_mod.o : src/numerics_mod.F90 \
 	 $(OBJDIR)/prec_const_mod.o  $(OBJDIR)/basic_mod.o $(OBJDIR)/coeff_mod.o \
	 $(OBJDIR)/utility_mod.o $(OBJDIR)/species_mod.o $(OBJDIR)/geometry_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/numerics_mod.F90 -o $@

 $(OBJDIR)/parallel_mod.o : src/parallel_mod.F90 \
 	 $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/parallel_mod.F90 -o $@

 $(OBJDIR)/ppexit.o : src/ppexit.F90 \
 	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/coeff_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/ppexit.F90 -o $@

 $(OBJDIR)/prec_const_mod.o : src/prec_const_mod.F90
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/prec_const_mod.F90 -o $@

 $(OBJDIR)/processing_mod.o : src/processing_mod.F90 \
 	 $(OBJDIR)/array_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/grid_mod.o \
   $(OBJDIR)/fields_mod.o $(OBJDIR)/basic_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/processing_mod.F90 -o $@

 $(OBJDIR)/readinputs.o : src/readinputs.F90  \
   $(OBJDIR)/diagnostics_mod.o $(OBJDIR)/initial_mod.o $(OBJDIR)/model_mod.o  \
  $(OBJDIR)/prec_const_mod.o $(OBJDIR)/grid_mod.o $(OBJDIR)/time_integration_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/readinputs.F90 -o $@

 $(OBJDIR)/restarts_mod.o : src/restarts_mod.F90   \
	 $(OBJDIR)/grid_mod.o $(OBJDIR)/time_integration_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/restarts_mod.F90 -o $@

 $(OBJDIR)/solve_EM_fields.o : src/solve_EM_fields.F90 \
   $(OBJDIR)/array_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/grid_mod.o \
	 $(OBJDIR)/ghosts_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/array_mod.o \
	 $(OBJDIR)/time_integration_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/parallel_mod.o \
	 $(OBJDIR)/processing_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/solve_EM_fields.F90 -o $@

 $(OBJDIR)/species_mod.o : src/species_mod.F90 \
   $(OBJDIR)/basic_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/species_mod.F90 -o $@

 $(OBJDIR)/stepon.o : src/stepon.F90 \
 	 $(OBJDIR)/initial_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/advance_field_mod.o \
   $(OBJDIR)/basic_mod.o $(OBJDIR)/nonlinear_mod.o $(OBJDIR)/grid_mod.o \
	 $(OBJDIR)/array_mod.o $(OBJDIR)/numerics_mod.o $(OBJDIR)/fields_mod.o \
	 $(OBJDIR)/ghosts_mod.o $(OBJDIR)/moments_eq_rhs_mod.o $(OBJDIR)/solve_EM_fields.o\
	 $(OBJDIR)/utility_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/time_integration_mod.o \
	 $(OBJDIR)/CLA_mod.o $(OBJDIR)/ExB_shear_flow_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/stepon.F90 -o $@

 $(OBJDIR)/tesend.o : src/tesend.F90 \
 	 $(OBJDIR)/basic_mod.o $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/tesend.F90 -o $@

 $(OBJDIR)/time_integration_mod.o : src/time_integration_mod.F90 \
   $(OBJDIR)/basic_mod.o  $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/time_integration_mod.F90 -o $@

 $(OBJDIR)/utility_mod.o : src/utility_mod.F90  \
   $(OBJDIR)/basic_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/time_integration_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/utility_mod.F90 -o $@

 $(OBJDIR)/CLA_mod.o : src/CLA_mod.F90  \
   $(OBJDIR)/basic_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/grid_mod.o \
   $(OBJDIR)/time_integration_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/CLA_mod.F90 -o $@
