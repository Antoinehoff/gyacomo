SUBROUTINE diagnose_profiler(kstep)
  USE basic
  USE diagnostics_par
  USE futils, ONLY: creatd,creatg,append,closef
  USE time_integration
  USE utility
  USE prec_const
  IMPLICIT NONE
  INTEGER, INTENT(in) :: kstep
  INTEGER :: dims(1) = (/0/)
  !____________________________________________________________________________!
  IF ((kstep .EQ. 0)) THEN
    CALL init_outfile(comm0,   prffile0,prffile,fidprf)

    ! data time measurement
    CALL creatd(fidprf, 0, dims, "/data/Tc_rhs",        "cumulative rhs computation time")
    CALL creatd(fidprf, 0, dims, "/data/Tc_adv_field",  "cumulative adv. fields computation time")
    CALL creatd(fidprf, 0, dims, "/data/Tc_clos",       "cumulative closure computation time")
    CALL creatd(fidprf, 0, dims, "/data/Tc_ghost",       "cumulative communication time")
    CALL creatd(fidprf, 0, dims, "/data/Tc_coll",       "cumulative collision computation time")
    CALL creatd(fidprf, 0, dims, "/data/Tc_poisson",    "cumulative poisson computation time")
    CALL creatd(fidprf, 0, dims, "/data/Tc_Sapj",       "cumulative Sapj computation time")
    CALL creatd(fidprf, 0, dims, "/data/Tc_checkfield", "cumulative checkfield computation time")
    CALL creatd(fidprf, 0, dims, "/data/Tc_diag",       "cumulative sym computation time")
    CALL creatd(fidprf, 0, dims, "/data/Tc_process",    "cumulative process computation time")
    CALL creatd(fidprf, 0, dims, "/data/Tc_step",       "cumulative total step computation time")
    CALL creatd(fidprf, 0, dims, "/data/time",          "current simulation time")
  ENDIF
  !_____________________________________________________________________________
  !                   2.   Periodic diagnostics
  !
  IF (kstep .GE. 0) THEN
  ! Time measurement data
  CALL append(fidprf, "/data/Tc_rhs",              tc_rhs,ionode=0)
  CALL append(fidprf, "/data/Tc_adv_field",  tc_adv_field,ionode=0)
  CALL append(fidprf, "/data/Tc_clos",            tc_clos,ionode=0)
  CALL append(fidprf, "/data/Tc_ghost",          tc_ghost,ionode=0)
  CALL append(fidprf, "/data/Tc_coll",            tc_coll,ionode=0)
  CALL append(fidprf, "/data/Tc_poisson",      tc_poisson,ionode=0)
  CALL append(fidprf, "/data/Tc_Sapj",            tc_Sapj,ionode=0)
  CALL append(fidprf, "/data/Tc_checkfield",tc_checkfield,ionode=0)
  CALL append(fidprf, "/data/Tc_diag",            tc_diag,ionode=0)
  CALL append(fidprf, "/data/Tc_process",      tc_process,ionode=0)
  CALL append(fidprf, "/data/Tc_step",            tc_step,ionode=0)
  CALL append(fidprf, "/data/time",                  time,ionode=0)
  !_____________________________________________________________________________
  !                   3.   Final diagnostics
  ELSEIF (kstep .EQ. -1) THEN
    !   Close diagnostic files
    CALL mpi_barrier(MPI_COMM_WORLD, ierr)
    CALL closef(fidprf)

  END IF
END SUBROUTINE diagnose_profiler
