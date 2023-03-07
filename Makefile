include local/dirs.inc
include local/make.inc
#Different namings depending on the make input
EXEC = $(BINDIR)/gyacomo	#all
EFST = $(BINDIR)/gyacomo_fast   #fast
EDBG = $(BINDIR)/gyacomo_debug	#debug
EALP = $(BINDIR)/gyacomo_alpha  #alpha
EGFT = $(BINDIR)/gyacomo_gfort  #gfort version
# F90 = mpiifort
F90 = mpif90
# #F90 = ftn #for piz-daint cluster
# # Add Multiple-Precision Library
# EXTLIBS += -L$(FMDIR)/lib
# EXTINC += -I$(FMDIR)/mod
# # Add local fftw dir
# EXTLIBS += -L$(FFTWDIR)/lib
# EXTINC += -I$(FFTWDIR)/include
# # Add lapack
# EXTLIBS += -L$(LAPACKDIR)/lib
# EXTINC += -I$(LAPACKDIR)/mod
# Standard version with optimized compilation
all: dirs src/srcinfo.h
all: F90FLAGS = -O3 -xHOST
all: $(EXEC)
# Fast compilation
fast: dirs src/srcinfo.h
fast: F90FLAGS = -fast
fast: $(EFST)
# Debug version with all flags
debug: dirs src/srcinfo.h
debug: F90FLAGS = -C -g -traceback -ftrapuv -warn all -debug all
# debug: F90FLAGS = -g -traceback -check all -ftrapuv -warn all -debug all
debug: $(EDBG)
# Alpha version, optimized as all but creates another binary
alpha: dirs src/srcinfo.h
alpha: F90FLAGS = -O3 -xHOST
alpha: $(EALP)
# gfortran version, compile with gfortran
gfort: dirs src/srcinfo.h
gfort: F90FLAGS = -g -std=legacy -ffree-line-length-0
gfort: EXTMOD   = -J $(MODDIR)
gfort: $(EGFT)
install: dirs src/srcinfo.h $(EXEC) mvmod

run: all
	(cd wk; $(EXEC);)

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

cleanbin:
	@rm -f $(EXEC)

mvmod:
	mv *.mod mod/.

$(OBJDIR)/diagnose.o : src/srcinfo.h

FOBJ=$(OBJDIR)/advance_field_mod.o $(OBJDIR)/array_mod.o $(OBJDIR)/auxval.o \
$(OBJDIR)/basic_mod.o $(OBJDIR)/coeff_mod.o $(OBJDIR)/closure_mod.o \
$(OBJDIR)/collision_mod.o $(OBJDIR)/nonlinear_mod.o $(OBJDIR)/control.o \
$(OBJDIR)/diagnose.o $(OBJDIR)/diagnostics_par_mod.o $(OBJDIR)/endrun.o \
$(OBJDIR)/fields_mod.o $(OBJDIR)/fourier_mod.o $(OBJDIR)/geometry_mod.o \
$(OBJDIR)/ghosts_mod.o $(OBJDIR)/grid_mod.o $(OBJDIR)/inital.o \
$(OBJDIR)/initial_par_mod.o $(OBJDIR)/lag_interp_mod.o $(OBJDIR)/main.o \
$(OBJDIR)/memory.o $(OBJDIR)/miller_mod.o $(OBJDIR)/model_mod.o \
$(OBJDIR)/moments_eq_rhs_mod.o $(OBJDIR)/numerics_mod.o $(OBJDIR)/parallel_mod.o \
$(OBJDIR)/ppexit.o $(OBJDIR)/prec_const_mod.o \
$(OBJDIR)/processing_mod.o $(OBJDIR)/readinputs.o $(OBJDIR)/restarts_mod.o \
$(OBJDIR)/solve_EM_fields.o $(OBJDIR)/species_mod.o $(OBJDIR)/stepon.o $(OBJDIR)/tesend.o \
$(OBJDIR)/time_integration_mod.o $(OBJDIR)/utility_mod.o

 $(EXEC): $(FOBJ)
	$(F90) $(LDFLAGS) $(OBJDIR)/*.o $(EXTMOD) $(EXTINC) $(EXTLIBS) -o $@

 $(EFST): $(FOBJ)
	$(F90) $(LDFLAGS) $(OBJDIR)/*.o $(EXTMOD) $(EXTINC) $(EXTLIBS) -o $@

 $(EDBG): $(FOBJ)
	$(F90) $(LDFLAGS) $(OBJDIR)/*.o $(EXTMOD) $(EXTINC) $(EXTLIBS) -o $@

 $(EALP): $(FOBJ)
	$(F90) $(LDFLAGS) $(OBJDIR)/*.o $(EXTMOD) $(EXTINC) $(EXTLIBS) -o $@

 $(EGFT): $(FOBJ)
	$(F90) $(LDFLAGS) $(OBJDIR)/*.o $(EXTMOD) $(EXTINC) $(EXTLIBS) -o $@

 $(OBJDIR)/advance_field_mod.o : src/advance_field_mod.F90 \
   $(OBJDIR)/grid_mod.o $(OBJDIR)/array_mod.o $(OBJDIR)/initial_par_mod.o \
	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/time_integration_mod.o $(OBJDIR)/basic_mod.o \
	 $(OBJDIR)/fields_mod.o $(OBJDIR)/model_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/advance_field_mod.F90 -o $@

 $(OBJDIR)/array_mod.o : src/array_mod.F90 \
	 $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/array_mod.F90 -o $@

 $(OBJDIR)/auxval.o : src/auxval.F90 \
 	 $(OBJDIR)/fourier_mod.o $(OBJDIR)/memory.o $(OBJDIR)/model_mod.o \
	 $(OBJDIR)/geometry_mod.o  $(OBJDIR)/grid_mod.o $(OBJDIR)/numerics_mod.o \
	 $(OBJDIR)/parallel_mod.o $(OBJDIR)/processing_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/auxval.F90 -o $@

 $(OBJDIR)/basic_mod.o : src/basic_mod.F90 \
 	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/parallel_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/basic_mod.F90 -o $@

 $(OBJDIR)/calculus_mod.o : src/calculus_mod.F90  \
 	 $(OBJDIR)/basic_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/grid_mod.o \
	 $(OBJDIR)/parallel_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/calculus_mod.F90 -o $@

 $(OBJDIR)/coeff_mod.o : src/coeff_mod.F90 \
 	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/model_mod.o \
	 $(OBJDIR)/basic_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/coeff_mod.F90 -o $@

 $(OBJDIR)/closure_mod.o : src/closure_mod.F90 \
 	 $(OBJDIR)/model_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/grid_mod.o \
	 $(OBJDIR)/array_mod.o $(OBJDIR)/fields_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/closure_mod.F90 -o $@

 $(OBJDIR)/collision_mod.o : src/collision_mod.F90 \
 	 $(OBJDIR)/array_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/cosolver_interface_mod.o\
	 $(OBJDIR)/fields_mod.o $(OBJDIR)/grid_mod.o $(OBJDIR)/model_mod.o \
	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/species_mod.o $(OBJDIR)/time_integration_mod.o \
	 $(OBJDIR)/utility_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/collision_mod.F90 -o $@

 $(OBJDIR)/control.o : src/control.F90 \
	 $(OBJDIR)/auxval.o $(OBJDIR)/geometry_mod.o $(OBJDIR)/prec_const_mod.o \
   $(OBJDIR)/basic_mod.o $(OBJDIR)/ppexit.o $(OBJDIR)/readinputs.o $(OBJDIR)/tesend.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/control.F90 -o $@

 $(OBJDIR)/cosolver_interface_mod.o : src/cosolver_interface_mod.F90 \
   $(OBJDIR)/grid_mod.o $(OBJDIR)/array_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/species_mod.o\
	 $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/cosolver_interface_mod.F90 -o $@

 $(OBJDIR)/diagnose.o : src/diagnose.F90 \
 	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/processing_mod.o $(OBJDIR)/array_mod.o \
   $(OBJDIR)/basic_mod.o $(OBJDIR)/diagnostics_par_mod.o $(OBJDIR)/fields_mod.o \
	 $(OBJDIR)/grid_mod.o $(OBJDIR)/initial_par_mod.o $(OBJDIR)/model_mod.o \
	 $(OBJDIR)/time_integration_mod.o\
	 $(OBJDIR)/parallel_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/diagnose.F90 -o $@

 $(OBJDIR)/diagnostics_par_mod.o : src/diagnostics_par_mod.F90 \
 	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/basic_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/diagnostics_par_mod.F90 -o $@

 $(OBJDIR)/endrun.o : src/endrun.F90 \
 	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/basic_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/endrun.F90 -o $@

 $(OBJDIR)/fields_mod.o : src/fields_mod.F90 \
   $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/fields_mod.F90 -o $@

 $(OBJDIR)/fourier_mod.o : src/fourier_mod.F90 \
 	 $(OBJDIR)/basic_mod.o $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/fourier_mod.F90 -o $@

 $(OBJDIR)/geometry_mod.o : src/geometry_mod.F90 \
 	 $(OBJDIR)/array_mod.o $(OBJDIR)/calculus_mod.o $(OBJDIR)/miller_mod.o \
   $(OBJDIR)/grid_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/prec_const_mod.o \
	 $(OBJDIR)/utility_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/geometry_mod.F90 -o $@

 $(OBJDIR)/ghosts_mod.o : src/ghosts_mod.F90 \
 	 $(OBJDIR)/basic_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/grid_mod.o\
   $(OBJDIR)/geometry_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/time_integration_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/ghosts_mod.F90 -o $@

 $(OBJDIR)/grid_mod.o : src/grid_mod.F90 \
 	 $(OBJDIR)/basic_mod.o $(OBJDIR)/fourier_mod.o $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/grid_mod.F90 -o $@

 $(OBJDIR)/inital.o : src/inital.F90 \
 	 $(OBJDIR)/array_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/fields_mod.o \
	 $(OBJDIR)/initial_par_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/numerics_mod.o \
	 $(OBJDIR)/solve_EM_fields.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/ghosts_mod.o \
	 $(OBJDIR)/grid_mod.o \
	 $(OBJDIR)/restarts_mod.o $(OBJDIR)/time_integration_mod.o $(OBJDIR)/utility_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/inital.F90 -o $@

 $(OBJDIR)/initial_par_mod.o : src/initial_par_mod.F90 \
 	 $(OBJDIR)/basic_mod.o $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/initial_par_mod.F90 -o $@

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
 	 $(OBJDIR)/basic_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/lag_interp_mod.o \
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
 	 $(OBJDIR)/array_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/fourier_mod.o \
	 $(OBJDIR)/fields_mod.o $(OBJDIR)/grid_mod.o $(OBJDIR)/model_mod.o\
	 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/time_integration_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/nonlinear_mod.F90 -o $@

 $(OBJDIR)/numerics_mod.o : src/numerics_mod.F90 \
 	 $(OBJDIR)/prec_const_mod.o  $(OBJDIR)/basic_mod.o $(OBJDIR)/coeff_mod.o \
	 $(OBJDIR)/utility_mod.o
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
   $(OBJDIR)/diagnostics_par_mod.o $(OBJDIR)/initial_par_mod.o $(OBJDIR)/model_mod.o  \
  $(OBJDIR)/prec_const_mod.o $(OBJDIR)/grid_mod.o $(OBJDIR)/time_integration_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/readinputs.F90 -o $@

 $(OBJDIR)/restarts_mod.o : src/restarts_mod.F90   \
	 $(OBJDIR)/diagnostics_par_mod.o $(OBJDIR)/grid_mod.o $(OBJDIR)/time_integration_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/restarts_mod.F90 -o $@

 $(OBJDIR)/solve_EM_fields.o : src/solve_EM_fields.F90 \
   $(OBJDIR)/array_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/grid_mod.o \
	 $(OBJDIR)/ghosts_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/array_mod.o \
	 $(OBJDIR)/time_integration_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/parallel_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/solve_EM_fields.F90 -o $@

 $(OBJDIR)/species_mod.o : src/species_mod.F90 \
   $(OBJDIR)/basic_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/species_mod.F90 -o $@

 $(OBJDIR)/stepon.o : src/stepon.F90 \
 	 $(OBJDIR)/initial_par_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/advance_field_mod.o \
   $(OBJDIR)/basic_mod.o $(OBJDIR)/nonlinear_mod.o $(OBJDIR)/grid_mod.o \
	 $(OBJDIR)/array_mod.o $(OBJDIR)/numerics_mod.o $(OBJDIR)/fields_mod.o \
	 $(OBJDIR)/ghosts_mod.o $(OBJDIR)/moments_eq_rhs_mod.o $(OBJDIR)/solve_EM_fields.o\
	 $(OBJDIR)/utility_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/time_integration_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/stepon.F90 -o $@

 $(OBJDIR)/tesend.o : src/tesend.F90 \
 	 $(OBJDIR)/basic_mod.o $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/tesend.F90 -o $@

 $(OBJDIR)/time_integration_mod.o : src/time_integration_mod.F90 \
   $(OBJDIR)/basic_mod.o  $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/time_integration_mod.F90 -o $@

 $(OBJDIR)/utility_mod.o : src/utility_mod.F90  \
   $(OBJDIR)/grid_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/utility_mod.F90 -o $@
