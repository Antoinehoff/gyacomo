include local/dirs.inc
include local/make.inc

MKLROOT = /usr/local/intel/composerxe/composer_xe_2015/mkl

EXEC = $(BINDIR)/helaz

F90 = mpif90

# Add Multiple-Precision Library
EXTLIBS += -L$(FMDIR)/lib
EXTINC += -I$(FMDIR)/mod

EXTLIBS += -L$(FFTWDIR)/lib64
EXTINC += -I$(FFTWDIR)/include

all: dirs src/srcinfo.h $(EXEC)

run: all
	(cd wk; $(EXEC);)

dirs:
	mkdir -p $(BINDIR)
	mkdir -p $(OBJDIR)
	mkdir -p $(MODDIR)
	mkdir -p $(CHCKPTDIR)

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
$(OBJDIR)/coeff_mod.o $(OBJDIR)/compute_Sapj.o $(OBJDIR)/control.o $(OBJDIR)/fourier_mod.o \
$(OBJDIR)/diagnose.o $(OBJDIR)/diagnostics_par_mod.o $(OBJDIR)/endrun.o $(OBJDIR)/fields_mod.o \
$(OBJDIR)/inital.o $(OBJDIR)/initial_par_mod.o $(OBJDIR)/main.o $(OBJDIR)/memory.o \
$(OBJDIR)/model_mod.o $(OBJDIR)/mkl_dfti.o $(OBJDIR)/moments_eq_rhs.o $(OBJDIR)/poisson.o \
$(OBJDIR)/ppexit.o $(OBJDIR)/ppinit.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/readinputs.o \
$(OBJDIR)/grid_mod.o $(OBJDIR)/stepon.o $(OBJDIR)/tesend.o $(OBJDIR)/time_integration_mod.o \
$(OBJDIR)/utility_mod.o

 $(EXEC): $(FOBJ)
	$(F90) $(LDFLAGS) $(OBJDIR)/*.o $(EXTMOD) $(EXTINC) $(EXTLIBS) -o $@

 $(OBJDIR)/advance_field.o : src/advance_field.F90 $(OBJDIR)/grid_mod.o $(OBJDIR)/array_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/time_integration_mod.o $(OBJDIR)/basic_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/advance_field.F90 -o $@

 $(OBJDIR)/array_mod.o : src/array_mod.F90 $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/array_mod.F90 -o $@

 $(OBJDIR)/auxval.o : src/auxval.F90 $(OBJDIR)/memory.o $(OBJDIR)/model_mod.o  $(OBJDIR)/grid_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/auxval.F90 -o $@

 $(OBJDIR)/basic_mod.o : src/basic_mod.F90 $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/basic_mod.F90 -o $@

 $(OBJDIR)/coeff_mod.o : src/coeff_mod.F90 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/basic_mod.o
		$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/coeff_mod.F90 -o $@

 $(OBJDIR)/compute_Sapj.o : src/compute_Sapj.F90 $(OBJDIR)/array_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/fourier_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/grid_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/time_integration_mod.o
		$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/compute_Sapj.F90 -o $@

 $(OBJDIR)/control.o : src/control.F90 $(OBJDIR)/auxval.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/ppexit.o $(OBJDIR)/ppinit.o $(OBJDIR)/readinputs.o $(OBJDIR)/tesend.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/control.F90 -o $@

 $(OBJDIR)/fourier_mod.o : src/fourier_mod.F90 $(OBJDIR)/basic_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/mkl_dfti.o $(OBJDIR)/grid_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/fourier_mod.F90 -o $@

 $(OBJDIR)/diagnose.o : src/diagnose.F90 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/array_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/diagnostics_par_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/grid_mod.o $(OBJDIR)/initial_par_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/time_integration_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/diagnose.F90 -o $@

 $(OBJDIR)/diagnostics_par_mod.o : src/diagnostics_par_mod.F90 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/basic_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/diagnostics_par_mod.F90 -o $@

 $(OBJDIR)/endrun.o : src/endrun.F90 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/basic_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/endrun.F90 -o $@

 $(OBJDIR)/fields_mod.o : src/fields_mod.F90 $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/fields_mod.F90 -o $@

 $(OBJDIR)/grid_mod.o : src/grid_mod.F90 $(OBJDIR)/basic_mod.o $(OBJDIR)/prec_const_mod.o
		$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/grid_mod.F90 -o $@

 $(OBJDIR)/inital.o : src/inital.F90 $(OBJDIR)/array_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/initial_par_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/poisson.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/grid_mod.o $(OBJDIR)/time_integration_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/inital.F90 -o $@

 $(OBJDIR)/initial_par_mod.o : src/initial_par_mod.F90 $(OBJDIR)/basic_mod.o $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/initial_par_mod.F90 -o $@

 $(OBJDIR)/main.o : src/main.F90 $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/main.F90 -o $@

 $(OBJDIR)/memory.o : src/memory.F90 $ $(OBJDIR)/array_mod.o $(OBJDIR)/basic_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/time_integration_mod.o $(OBJDIR)/grid_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/memory.F90 -o $@

 $(OBJDIR)/mkl_dfti.o : $(MKLROOT)/include/mkl_dfti.f90
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) $(MKLROOT)/include/mkl_dfti.f90 -o $@

 $(OBJDIR)/model_mod.o : src/model_mod.F90 $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/model_mod.F90 -o $@

 $(OBJDIR)/moments_eq_rhs.o : src/moments_eq_rhs.F90 $(OBJDIR)/array_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/grid_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/time_integration_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/moments_eq_rhs.F90 -o $@

 $(OBJDIR)/poisson.o : src/poisson.F90 $(OBJDIR)/array_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/grid_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/array_mod.o $(OBJDIR)/time_integration_mod.o $(OBJDIR)/basic_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/poisson.F90 -o $@

 $(OBJDIR)/ppexit.o : src/ppexit.F90 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/basic_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/ppexit.F90 -o $@

 $(OBJDIR)/ppinit.o : src/ppinit.F90 $(OBJDIR)/array_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/grid_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/array_mod.o $(OBJDIR)/time_integration_mod.o $(OBJDIR)/basic_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/ppinit.F90 -o $@

 $(OBJDIR)/prec_const_mod.o : src/prec_const_mod.F90
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/prec_const_mod.F90 -o $@

 $(OBJDIR)/readinputs.o : src/readinputs.F90  $(OBJDIR)/diagnostics_par_mod.o $(OBJDIR)/initial_par_mod.o $(OBJDIR)/model_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/grid_mod.o $(OBJDIR)/time_integration_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/readinputs.F90 -o $@

 $(OBJDIR)/stepon.o : src/stepon.F90 $(OBJDIR)/prec_const_mod.o $(OBJDIR)/advance_field.o $(OBJDIR)/basic_mod.o $(OBJDIR)/grid_mod.o $(OBJDIR)/array_mod.o $(OBJDIR)/fields_mod.o $(OBJDIR)/moments_eq_rhs.o $(OBJDIR)/poisson.o $(OBJDIR)/time_integration_mod.o $(OBJDIR)/utility_mod.o $(OBJDIR)/model_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/stepon.F90 -o $@

 $(OBJDIR)/tesend.o : src/tesend.F90 $(OBJDIR)/basic_mod.o $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/tesend.F90 -o $@

 $(OBJDIR)/time_integration_mod.o : src/time_integration_mod.F90 $(OBJDIR)/basic_mod.o $(OBJDIR)/prec_const_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/time_integration_mod.F90 -o $@

 $(OBJDIR)/utility_mod.o : src/utility_mod.F90  $(OBJDIR)/basic_mod.o $(OBJDIR)/prec_const_mod.o $(OBJDIR)/grid_mod.o
	$(F90) -c $(F90FLAGS) $(FPPFLAGS) $(EXTMOD) $(EXTINC) src/utility_mod.F90 -o $@
