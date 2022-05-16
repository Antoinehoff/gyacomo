SUBROUTINE diagnose(kstep)
  !   Diagnostics, writing simulation state to disk

  USE basic
  USE diagnostics_par
  IMPLICIT NONE

  INTEGER, INTENT(in) :: kstep

  CALL cpu_time(t0_diag) ! Measuring time

  !! Basic diagnose loop for reading input file, displaying advancement and ending
  IF ((kstep .EQ. 0)) THEN
    INQUIRE(unit=lu_in, name=input_fname)
    CLOSE(lu_in)
  ENDIF
  IF (kstep .GE. 0) THEN
    ! Terminal info
    IF (MOD(cstep, INT(1.0/dt)) == 0 .AND. (my_id .EQ. 0)) THEN
     WRITE(*,"(F6.0,A,F6.0)") time,"/",tmax
    ENDIF
  ELSEIF (kstep .EQ. -1) THEN
    CALL cpu_time(finish)
     ! Display computational time cost
     IF (my_id .EQ. 0) CALL display_h_min_s(finish-start)
  END IF
  !! Specific diagnostic calls
  CALL diagnose_full(kstep)
  ! IF(nsave_5d .GT. 0) CALL diagnose_moments(kstep)
  ! IF(nsave_3d .GT. 0) CALL diagnose_momspectrum(kstep)
  ! IF(nsave_3d .GT. 0) CALL diagnose_fields(kstep)
  ! IF(nsave_0d .GT. 0) CALL diagnose_profiler(kstep)
  ! IF(nsave_0d .GT. 0) CALL diagnose_gridgeom(kstep)
  ! IF(nsave_0d .GT. 0) CALL diagnose_timetraces(kstep)

  CALL cpu_time(t1_diag); tc_diag = tc_diag + (t1_diag - t0_diag)


END SUBROUTINE diagnose

SUBROUTINE init_outfile(comm,file0,file,fid)
  USE diagnostics_par, ONLY : write_doubleprecision, diag_par_outputinputs, input_fname
  USE basic,           ONLY : my_id, jobnum, lu_in, basic_outputinputs
  USE grid,            ONLY : grid_outputinputs
  USE geometry,        ONLY : geometry_outputinputs
  USE model,           ONLY : model_outputinputs
  USE collision,       ONLY : coll_outputinputs
  USE initial_par,     ONLY : initial_outputinputs
  USE time_integration,ONLY : time_integration_outputinputs
  USE futils,          ONLY : creatf, creatg, creatd, attach, putfile
  IMPLICIT NONE
  !input
  INTEGER,            INTENT(IN)    :: comm
  CHARACTER(len=256), INTENT(IN)    :: file0
  CHARACTER(len=256), INTENT(OUT)   :: file
  INTEGER,            INTENT(OUT)   :: fid
  CHARACTER(len=256)                :: str,fname
  INCLUDE 'srcinfo.h'

  ! Writing output filename
  WRITE(file,'(a,a1,i2.2,a3)') TRIM(file0)   ,'_',jobnum,'.h5'
  !                      1.1   Initial run
  ! Main output file creation
  IF (write_doubleprecision) THEN
    CALL creatf(file, fid, real_prec='d', mpicomm=comm)
  ELSE
    CALL creatf(file, fid, mpicomm=comm)
  END IF
  IF (my_id .EQ. 0) WRITE(*,'(3x,a,a)')  TRIM(file), ' created'
  !  basic data group
  CALL creatg(fid, "/data", "data")
  !  File group
  CALL creatg(fid, "/files", "files")
  CALL attach(fid, "/files",  "jobnum", jobnum)

  !  Add input namelist variables as attributes of /data/input, defined in srcinfo.h
  ! IF (my_id .EQ. 0) WRITE(*,*) 'VERSION=', VERSION
  ! IF (my_id .EQ. 0) WRITE(*,*)  'BRANCH=', BRANCH
  ! IF (my_id .EQ. 0) WRITE(*,*)  'AUTHOR=', AUTHOR
  ! IF (my_id .EQ. 0) WRITE(*,*)    'HOST=', HOST

  ! Add the code info and parameters to the file
  WRITE(str,'(a,i2.2)') "/data/input"
  CALL creatd(fid, 0,(/0/),TRIM(str),'Input parameters')
  CALL attach(fid, TRIM(str),     "version",  VERSION) !defined in srcinfo.h
  CALL attach(fid, TRIM(str),      "branch",   BRANCH) !defined in srcinfo.h
  CALL attach(fid, TRIM(str),      "author",   AUTHOR) !defined in srcinfo.h
  CALL attach(fid, TRIM(str),    "execdate", EXECDATE) !defined in srcinfo.h
  CALL attach(fid, TRIM(str),        "host",     HOST) !defined in srcinfo.h

  CALL basic_outputinputs(fid,str)
  CALL grid_outputinputs(fid, str)
  CALL geometry_outputinputs(fid, str)
  CALL diag_par_outputinputs(fid, str)
  CALL model_outputinputs(fid, str)
  CALL coll_outputinputs(fid, str)
  CALL initial_outputinputs(fid, str)
  CALL time_integration_outputinputs(fid, str)

  !  Save STDIN (input file) of this run
  IF(jobnum .LE. 99) THEN
     WRITE(str,'(a,i2.2)') "/files/STDIN.",jobnum
  ELSE
     WRITE(str,'(a,i3.2)') "/files/STDIN.",jobnum
  END IF
  CALL putfile(fid, TRIM(str), TRIM(input_fname),ionode=0)
END SUBROUTINE init_outfile

!! Auxiliary routines hidden in headers
INCLUDE 'diag_headers/diagnose_full.h'
! INCLUDE 'diag_headers/diagnose_moments.h'
! INCLUDE 'diag_headers/diagnose_momspectrum.h'
! INCLUDE 'diag_headers/diagnose_fields.h'
! INCLUDE 'diag_headers/diagnose_profiler.h'
! INCLUDE 'diag_headers/diagnose_gridgeom.h'
! INCLUDE 'diag_headers/diagnose_timetraces.h'
