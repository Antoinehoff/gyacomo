include local/dirs.inc
include local/make.inc


EXEC = $(BINDIR)/monli1D

F90 = mpif90


all: dirs src/srcinfo.h $(EXEC)

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

$(OBJDIR)/diagnose.o : src/srcinfo.h

FOBJ=$(OBJDIR)/advance_field.o $(OBJDIR)/array_mod.o $(OBJDIR)/auxval.o $(OBJDIR)/basic_mod.o \
$(OBJDIR)/control.o $(OBJDIR)/diagnose.o $(OBJDIR)/diagnostics_par_mod.o $(OBJDIR)/endrun.o \
$(OBJDIR)/fft_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/gradients.o $(OBJDIR)/inital.o \
$(OBJDIR)/initial_par_mod.o $(OBJDIR)/main.o $(OBJDIR)/memory.o $(OBJDIR)/model_mod.o \
$(OBJDIR)/moments_eq_rhs.o $(OBJDIR)/numerics.o $(OBJDIR)/poisson.o $(OBJDIR)/ppexit.o \
$(OBJDIR)/ppinit.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/readinputs.o $(OBJDIR)/space_grid_mod.o \
$(OBJDIR)/stepon.o $(OBJDIR)/temp_eq_rhs.o $(OBJDIR)/tesend.o $(OBJDIR)/theta_eq_rhs.o \
$(OBJDIR)/time_integration_mod.o $(OBJDIR)/utility_mod.o $(OBJDIR)/vpar_eq_rhs.o

$(EXEC): $(FOBJ)
	$(F90) $(LDFLAGS) $(OBJDIR)/*.o $(EXTMOD) $(EXTINC) $(EXTLIBS) -o $@

$(OBJDIR)/advance_field.o : src/advance_field.F90 $(OBJDIR)/space_grid_mod.o $(OBJDIR)/array_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/time_integration_mod.o $(OBJDIR)/basic_mod.o 
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/advance_field.F90 -o $@
$(OBJDIR)/array_mod.o : src/array_mod.F90 $(OBJDIR)/prec_const_mod.o 
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/array_mod.F90 -o $@
$(OBJDIR)/auxval.o : src/auxval.F90 $(OBJDIR)/fft_mod.o $(OBJDIR)/gradients.o $(OBJDIR)/memory.o $(OBJDIR)/model_mod.o  $(OBJDIR)/space_grid_mod.o 
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/auxval.F90 -o $@
$(OBJDIR)/basic_mod.o : src/basic_mod.F90 $(OBJDIR)/prec_const_mod.o 
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/basic_mod.F90 -o $@
$(OBJDIR)/control.o : src/control.F90 $(OBJDIR)/auxval.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/ppexit.o $(OBJDIR)/ppinit.o $(OBJDIR)/readinputs.o $(OBJDIR)/tesend.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/control.F90 -o $@
$(OBJDIR)/diagnose.o : src/diagnose.F90 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/diagnostics_par_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/initial_par_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/time_integration_mod.o 
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/diagnose.F90 -o $@
$(OBJDIR)/diagnostics_par_mod.o : src/diagnostics_par_mod.F90 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/basic_mod.o 
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/diagnostics_par_mod.F90 -o $@
$(OBJDIR)/endrun.o : src/endrun.F90 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/basic_mod.o 
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/endrun.F90 -o $@
$(OBJDIR)/fft_mod.o : src/fft_mod.F90 $(OBJDIR)/array_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/prec_const_mod.o 
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/fft_mod.F90 -o $@
$(OBJDIR)/fields_mod.o : src/fields_mod.F90 $(OBJDIR)/prec_const_mod.o 
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/fields_mod.F90 -o $@
$(OBJDIR)/gradients.o : src/gradients.F90 $(OBJDIR)/basic_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/space_grid_mod.o 
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/gradients.F90 -o $@
$(OBJDIR)/inital.o : src/inital.F90 $(OBJDIR)/array_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/initial_par_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/poisson.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/space_grid_mod.o $(OBJDIR)/time_integration_mod.o 
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/inital.F90 -o $@
$(OBJDIR)/initial_par_mod.o : src/initial_par_mod.F90 $(OBJDIR)/basic_mod.o $(OBJDIR)/prec_const_mod.o 
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/initial_par_mod.F90 -o $@
$(OBJDIR)/main.o : src/main.F90 $(OBJDIR)/prec_const_mod.o 
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/main.F90 -o $@
$(OBJDIR)/memory.o : src/memory.F90 $ $(OBJDIR)/array_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/time_integration_mod.o $(OBJDIR)/space_grid_mod.o 
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/memory.F90 -o $@
$(OBJDIR)/model_mod.o : src/model_mod.F90 $(OBJDIR)/prec_const_mod.o 
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/model_mod.F90 -o $@
$(OBJDIR)/moments_eq_rhs.o : src/moments_eq_rhs.F90 $(OBJDIR)/array_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/space_grid_mod.o $(OBJDIR)/time_integration_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/moments_eq_rhs.F90 -o $@
$(OBJDIR)/numerics.o : src/numerics.F90 $(OBJDIR)/array_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/space_grid_mod.o $(OBJDIR)/time_integration_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/numerics.F90 -o $@
$(OBJDIR)/poisson.o : src/poisson.F90 $(OBJDIR)/array_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/space_grid_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/array_mod.o $(OBJDIR)/time_integration_mod.o $(OBJDIR)/basic_mod.o 
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/poisson.F90 -o $@
$(OBJDIR)/ppexit.o : src/ppexit.F90 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/basic_mod.o 
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/ppexit.F90 -o $@
$(OBJDIR)/ppinit.o : src/ppinit.F90 $(OBJDIR)/array_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/space_grid_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/array_mod.o $(OBJDIR)/time_integration_mod.o $(OBJDIR)/basic_mod.o 
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/ppinit.F90 -o $@
$(OBJDIR)/prec_const_mod.o : src/prec_const_mod.F90 
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/prec_const_mod.F90 -o $@
$(OBJDIR)/readinputs.o : src/readinputs.F90  $(OBJDIR)/diagnostics_par_mod.o $(OBJDIR)/initial_par_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/space_grid_mod.o $(OBJDIR)/time_integration_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/readinputs.F90 -o $@
$(OBJDIR)/space_grid_mod.o : src/space_grid_mod.F90 $(OBJDIR)/basic_mod.o $(OBJDIR)/prec_const_mod.o 
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/space_grid_mod.F90 -o $@
$(OBJDIR)/stepon.o : src/stepon.F90 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/advance_field.o $(OBJDIR)/basic_mod.o $(OBJDIR)/numerics.o $(OBJDIR)/space_grid_mod.o $(OBJDIR)/array_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/moments_eq_rhs.o $(OBJDIR)/poisson.o $(OBJDIR)/temp_eq_rhs.o $(OBJDIR)/theta_eq_rhs.o $(OBJDIR)/time_integration_mod.o $(OBJDIR)/utility_mod.o $(OBJDIR)/vpar_eq_rhs.o $(OBJDIR)/model_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/stepon.F90 -o $@
$(OBJDIR)/temp_eq_rhs.o : src/temp_eq_rhs.F90 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/space_grid_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/array_mod.o $(OBJDIR)/time_integration_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/model_mod.o 
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/temp_eq_rhs.F90 -o $@
$(OBJDIR)/tesend.o : src/tesend.F90 $(OBJDIR)/basic_mod.o $(OBJDIR)/prec_const_mod.o 
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/tesend.F90 -o $@
$(OBJDIR)/theta_eq_rhs.o : src/theta_eq_rhs.F90 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/space_grid_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/array_mod.o $(OBJDIR)/time_integration_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/model_mod.o 
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/theta_eq_rhs.F90 -o $@
$(OBJDIR)/time_integration_mod.o : src/time_integration_mod.F90 $(OBJDIR)/basic_mod.o $(OBJDIR)/prec_const_mod.o 
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/time_integration_mod.F90 -o $@
$(OBJDIR)/utility_mod.o : src/utility_mod.F90  $(OBJDIR)/basic_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/space_grid_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/utility_mod.F90 -o $@
$(OBJDIR)/vpar_eq_rhs.o : src/vpar_eq_rhs.F90 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/space_grid_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/array_mod.o $(OBJDIR)/time_integration_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/model_mod.o 
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/vpar_eq_rhs.F90 -o $@

