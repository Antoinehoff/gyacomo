SUBROUTINE diagnose_timetraces(kstep)
  USE basic
  USE grid
  USE diagnostics_par
  USE futils
  USE array
  USE model
  USE initial_par
  USE fields
  USE time_integration
  USE utility
  USE prec_const
  USE collision, ONLY: coll_outputinputs
  USE geometry
  USE processing
  IMPLICIT NONE
  INTEGER, INTENT(in) :: kstep
  INTEGER, parameter  :: BUFSIZE = 2
  INTEGER :: rank = 0
  INTEGER :: dims(1) = (/0/)
  !____________________________________________________________________________!
  IF ((kstep .EQ. 0)) THEN
    CALL init_outfile(comm0,   ttrfile0,ttrfile,fidttr)

    !  var0d group (gyro transport)
    CALL creatd(fidttr, rank, dims,  "/data/time",     "Time t*c_s/R")
    CALL creatd(fidttr, rank, dims, "/data/cstep", "iteration number")

    IF (write_gamma) THEN
     CALL creatd(fidttr, rank, dims, "/data/gflux_ri", "Radial gyro ion transport")
     CALL creatd(fidttr, rank, dims, "/data/pflux_ri", "Radial part ion transport")
    ENDIF
    IF (write_hf) THEN
     CALL creatd(fidttr, rank, dims, "/data/hflux_x", "Radial part ion heat flux")
    ENDIF
    IF (cstep==0) THEN
     iframe0d=0
    ENDIF
    CALL attach(fidttr,"/data/" , "frames", iframe0d)
  ENDIF

  !_____________________________________________________________________________
  !                   2.   Periodic diagnostics
  !
  IF (kstep .GE. 0) THEN

  !                       0d time traces arrays
  IF ( MOD(cstep, nsave_0d) == 0 ) THEN
   ! Processing data
   CALL append(fidttr,  "/data/time",           time,ionode=0)
   CALL append(fidttr, "/data/cstep", real(cstep,dp),ionode=0)
   CALL getatt(fidttr,      "/data/",       "frames",iframe0d)
   iframe0d=iframe0d+1
   CALL attach(fidttr,"/data/" , "frames", iframe0d)
   ! Ion transport data
   IF (write_gamma) THEN
     CALL compute_radial_ion_transport
     CALL append(fidttr, "/data/gflux_ri",gflux_ri,ionode=0)
     CALL append(fidttr, "/data/pflux_ri",pflux_ri,ionode=0)
   ENDIF
   IF (write_hf) THEN
     CALL compute_radial_heatflux
     CALL append(fidttr, "/data/hflux_x",hflux_x,ionode=0)
   ENDIF
  END IF
  !_____________________________________________________________________________
  !                   3.   Final diagnostics
  ELSEIF (kstep .EQ. -1) THEN
    !   Close diagnostic files
    CALL mpi_barrier(MPI_COMM_WORLD, ierr)
    CALL closef(fidttr)

  END IF
END SUBROUTINE diagnose_timetraces
!____________________________________________________________________________!
!                       AUXILIARY ROUTINES
