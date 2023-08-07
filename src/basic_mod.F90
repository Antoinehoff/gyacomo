MODULE basic
  !   Basic module for time dependent problems
  use, intrinsic :: iso_c_binding
  use prec_const, ONLY : xp
  IMPLICIT none
  PRIVATE
  ! INPUT PARAMETERS
  INTEGER,  PUBLIC, PROTECTED :: nrun       = 1        ! Number of time steps to run
  real(xp), PUBLIC, PROTECTED :: tmax       = 100000.0 ! Maximum simulation time
  real(xp), PUBLIC, PROTECTED :: dt         = 1.0      ! Time step
  real(xp), PUBLIC, PROTECTED :: maxruntime = 1e9      ! Maximum simulation CPU time
  INTEGER,  PUBLIC, PROTECTED :: job2load   = 99       ! jobnum of the checkpoint to load
  ! Auxiliary variables
  real(xp), PUBLIC, PROTECTED :: time   = 0            ! Current simulation time (Init from restart file)

  INTEGER, PUBLIC, PROTECTED  :: jobnum  = 0           ! Job number
  INTEGER, PUBLIC, PROTECTED  :: step    = 0           ! Calculation step of this run
  INTEGER, PUBLIC, PROTECTED  :: cstep   = 0           ! Current step number (Init from restart file)
  LOGICAL, PUBLIC             :: nlend   = .FALSE.     ! Signal end of run
  LOGICAL, PUBLIC             :: crashed = .FALSE.     ! Signal end of crashed run
  INTEGER, PUBLIC :: iframe0d ! counting the number of times 0d datasets are outputed (for diagnose)
  INTEGER, PUBLIC :: iframe1d ! counting the number of times 1d datasets are outputed (for diagnose)
  INTEGER, PUBLIC :: iframe2d ! counting the number of times 2d datasets are outputed (for diagnose)
  INTEGER, PUBLIC :: iframe3d ! counting the number of times 3d datasets are outputed (for diagnose)
  INTEGER, PUBLIC :: iframe5d ! counting the number of times 5d datasets are outputed (for diagnose)
  !  List of logical file units
  INTEGER, PUBLIC, PROTECTED  :: lu_in   = 90              ! File duplicated from STDIN
  INTEGER, PUBLIC, PROTECTED  :: lu_stop = 91              ! stop file, see subroutine TESEND
  ! To measure computation time
  type :: chrono
    real(xp) :: tstart !start of the chrono
    real(xp) :: tstop  !stop 
    real(xp) :: ttot   !cumulative time
  end type chrono
  ! Define the chronos for each relevant routines
  type(chrono), PUBLIC, PROTECTED :: chrono_runt, chrono_mrhs, chrono_advf, chrono_pois, chrono_sapj,&
   chrono_diag, chrono_chck, chrono_step, chrono_clos, chrono_ghst, chrono_coll, chrono_napj, &
   chrono_grad, chrono_ExBs
#ifdef TEST_SVD
  ! A chrono for SVD tests
  type(chrono), PUBLIC, PROTECTED :: chrono_CLA
#endif
  ! This sets if the outputs is done through a large gather or using parallelization from futils
  !  it is recommended to set it to .true.
  LOGICAL, PUBLIC, PROTECTED :: GATHERV_OUTPUT = .true.
  ! Routines interfaces
  PUBLIC :: allocate_array, basic_outputinputs,basic_data,&
            speak, str, increase_step, increase_cstep, increase_time, display_h_min_s,&
            set_basic_cp, daytim, start_chrono, stop_chrono, change_dt
  ! Interface for allocating arrays, these routines allocate and initialize directly to zero
  INTERFACE allocate_array
    MODULE PROCEDURE allocate_array_xp1,allocate_array_xp2,allocate_array_xp3, &
                     allocate_array_xp4, allocate_array_xp5, allocate_array_xp6, allocate_array_xp7
    MODULE PROCEDURE allocate_array_dc1,allocate_array_dc2,allocate_array_dc3, &
                     allocate_array_dc4, allocate_array_dc5, allocate_array_dc6, allocate_array_dc7
    MODULE PROCEDURE allocate_array_i1,allocate_array_i2,allocate_array_i3,allocate_array_i4
    MODULE PROCEDURE allocate_array_l1,allocate_array_l2,allocate_array_l3,allocate_array_l4
  END INTERFACE allocate_array

  INTERFACE str
    MODULE PROCEDURE str_xp, str_int
  END INTERFACE

CONTAINS
  !================================================================================
  SUBROUTINE basic_data
    !   Read basic data for input file

    use prec_const
    IMPLICIT NONE

    NAMELIST /BASIC/  nrun, dt, tmax, maxruntime, job2load

    CALL find_input_file

    READ(lu_in,basic)

    !Init chronometers
    chrono_mrhs%ttot = 0
    chrono_pois%ttot = 0
    chrono_sapj%ttot = 0
    chrono_napj%ttot = 0
    chrono_grad%ttot = 0
    chrono_advf%ttot = 0
    chrono_ghst%ttot = 0
    chrono_clos%ttot = 0
    chrono_chck%ttot = 0
    chrono_diag%ttot = 0
    chrono_step%ttot = 0
    chrono_ExBs%ttot = 0
#ifdef TEST_SVD
    chrono_CLA%ttot = 0
#endif
  END SUBROUTINE basic_data


  SUBROUTINE basic_outputinputs(fid)
    !
    !    Write the input parameters to the results_xx.h5 file
    !
    USE prec_const
    USE futils, ONLY: attach, creatd
    IMPLICIT NONE
    INTEGER, INTENT(in) :: fid
    CHARACTER(len=256)  :: str
    WRITE(str,'(a)') '/data/input/basic'
    CALL creatd(fid, 0,(/0/),TRIM(str),'Basic Input')
    CALL attach(fid, TRIM(str), "start_iframe0d", iframe0d)
    CALL attach(fid, TRIM(str), "start_iframe2d", iframe2d)
    CALL attach(fid, TRIM(str), "start_iframe3d", iframe3d)
    CALL attach(fid, TRIM(str), "start_iframe5d", iframe5d)
    CALL attach(fid, TRIM(str),  "start_time",     time)
    CALL attach(fid, TRIM(str), "start_cstep",    cstep-1)
    CALL attach(fid, TRIM(str),          "dt",       dt)
    CALL attach(fid, TRIM(str),        "tmax",     tmax)
    CALL attach(fid, TRIM(str),        "nrun",     nrun)
    CALL attach(fid, TRIM(str),    "cpu_time",       -1)
  END SUBROUTINE basic_outputinputs
  !! Increments private attributes
  SUBROUTINE increase_step
    IMPLICIT NONE
    step  = step  + 1
  END SUBROUTINE
  SUBROUTINE increase_cstep
    IMPLICIT NONE
    cstep  = cstep  + 1
  END SUBROUTINE
  SUBROUTINE change_dt(new_dt)
    IMPLICIT NONE
    REAL(xp), INTENT(IN) :: new_dt   
    dt = new_dt
  END SUBROUTINE
  SUBROUTINE increase_time
    IMPLICIT NONE
    time  = time  + dt
  END SUBROUTINE
  SUBROUTINE set_basic_cp(cstep_cp,time_cp,jobnum_cp)
    IMPLICIT NONE
    REAL(xp), INTENT(IN) :: time_cp
    INTEGER,  INTENT(IN) :: cstep_cp, jobnum_cp
    cstep  = cstep_cp
    time   = time_cp
    jobnum = jobnum_cp+1
  END SUBROUTINE
  !! Chrono handling
  SUBROUTINE start_chrono(timer)
    IMPLICIT NONE
    type(chrono) :: timer
    CALL cpu_time(timer%tstart)
  END SUBROUTINE
  SUBROUTINE stop_chrono(timer)
    IMPLICIT NONE
    type(chrono) :: timer
    CALL cpu_time(timer%tstop)
    timer%ttot = timer%ttot + (timer%tstop-timer%tstart)
  END SUBROUTINE
  !================================================================================
  ! routine to speak in the terminal
  SUBROUTINE speak(message)
    USE parallel, ONLY: my_id
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(in) :: message
    IF(my_id .EQ. 0) write(*,*) message
  END SUBROUTINE
  !================================================================================
  SUBROUTINE find_input_file
    USE parallel, ONLY: my_id
    IMPLICIT NONE
    CHARACTER(len=32) :: str_, input_file
    INTEGER :: nargs, fileid, l, ierr
    LOGICAL :: mlexist
    nargs = COMMAND_ARGUMENT_COUNT()
    IF((nargs .EQ. 1) .OR. (nargs .EQ. 4)) THEN
      CALL GET_COMMAND_ARGUMENT(nargs, str_, l, ierr)
      READ(str_(1:l),'(i3)')  fileid
      WRITE(input_file,'(a,a1,i2.2,a3)') 'fort','_',fileid,'.90'

      INQUIRE(file=input_file, exist=mlexist)
      IF( mlexist ) THEN
        IF(my_id.EQ.0) write(*,*) 'Reading input ', input_file,'...'
        OPEN(lu_in, file=input_file)
      ELSE
        IF(my_id.EQ.0) write(*,*) 'Reading input fort.90...'
        OPEN(lu_in, file='fort.90')
      ENDIF
    ENDIF
  END SUBROUTINE find_input_file
  !================================================================================
  SUBROUTINE daytim(str)
    !   Print date and time
    USE parallel, ONLY: my_id
    use prec_const
    IMPLICIT NONE

    CHARACTER(len=*) , INTENT(in) :: str
    CHARACTER(len=16) :: d, t, dat, time
    !________________________________________________________________________________
    !
    CALL DATE_AND_TIME(d,t)
    dat=d(7:8) // '/' // d(5:6) // '/' // d(1:4)
    time=t(1:2) // ':' // t(3:4) // ':' // t(5:10)
    IF (my_id .EQ. 0) &
      WRITE(*,'(a,1x,a,1x,a)') str, dat(1:10), time(1:12)
    !
  END SUBROUTINE daytim
  !================================================================================
  SUBROUTINE display_h_min_s(time)
    USE parallel, ONLY: my_id
    IMPLICIT NONE
    real(xp) :: time
    integer  :: days, hours, mins, secs
    days = FLOOR(time/24./3600.);
    hours= FLOOR(time/3600.);
    mins = FLOOR(time/60.);
    secs = FLOOR(time);

    IF ( days .GT. 0 ) THEN !display day h min s
      hours = (time/3600./24. - days) * 24
      mins  = (time/3600. - days*24. - hours) * 60
      secs  = (time/60. - days*24.*60 - hours*60 - mins) * 60
      IF (my_id .EQ. 0) WRITE(*,*) 'CPU Time = ', days, '[day]', hours, '[h]', mins, '[min]', secs, '[s]'
      IF (my_id .EQ. 0) WRITE(*,*) '(',time,'[s])'

    ELSEIF ( hours .GT. 0 ) THEN !display h min s
      mins  = (time/3600. - hours) * 60
      secs  = (time/60. - hours*60 - mins) * 60
      IF (my_id .EQ. 0) WRITE(*,*) 'CPU Time = ', hours, '[h]', mins, '[min]', secs, '[s]'
      IF (my_id .EQ. 0) WRITE(*,*) '(',time,'[s])'

    ELSEIF ( mins .GT. 0 ) THEN !display min s
      secs  = (time/60. - mins) * 60
      IF (my_id .EQ. 0) WRITE(*,*) 'CPU Time = ', mins, '[min]', secs, '[s]'
      IF (my_id .EQ. 0) WRITE(*,*) '(',time,'[s])'

    ELSE ! display s
      IF (my_id .EQ. 0) WRITE(*,*) 'CPU Time = ', FLOOR(time), '[s]'

    ENDIF
  END SUBROUTINE display_h_min_s
!================================================================================

  function str_xp(k) result( str_ )
  !   "Convert an integer to string."
      REAL(xp), intent(in) :: k
      character(len=10):: str_
      write (str_, "(G10.2)") k
      str_ = adjustl(str_)
  end function str_xp

  function str_int(k) result( str_ )
  !   "Convert an integer to string."
      integer, intent(in) :: k
      character(len=10)   :: str_
      write (str_, "(I8)") k
      str_ = adjustl(str_)
  end function str_int

! To allocate arrays of doubles, integers, etc. at run time
  SUBROUTINE allocate_array_xp1(a,is1,ie1)
    IMPLICIT NONE
    real(xp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1
    ALLOCATE(a(is1:ie1))
    a=0.0_xp
  END SUBROUTINE allocate_array_xp1

  SUBROUTINE allocate_array_xp2(a,is1,ie1,is2,ie2)
    IMPLICIT NONE
    real(xp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2
    ALLOCATE(a(is1:ie1,is2:ie2))
    a=0.0_xp
  END SUBROUTINE allocate_array_xp2

  SUBROUTINE allocate_array_xp3(a,is1,ie1,is2,ie2,is3,ie3)
    IMPLICIT NONE
    real(xp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2,is3,ie3
    ALLOCATE(a(is1:ie1,is2:ie2,is3:ie3))
    a=0.0_xp
  END SUBROUTINE allocate_array_xp3

  SUBROUTINE allocate_array_xp4(a,is1,ie1,is2,ie2,is3,ie3,is4,ie4)
    IMPLICIT NONE
    real(xp), DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2,is3,ie3,is4,ie4
    ALLOCATE(a(is1:ie1,is2:ie2,is3:ie3,is4:ie4))
    a=0.0_xp
  END SUBROUTINE allocate_array_xp4

  SUBROUTINE allocate_array_xp5(a,is1,ie1,is2,ie2,is3,ie3,is4,ie4,is5,ie5)
    IMPLICIT NONE
    real(xp), DIMENSION(:,:,:,:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2,is3,ie3,is4,ie4,is5,ie5
    ALLOCATE(a(is1:ie1,is2:ie2,is3:ie3,is4:ie4,is5:ie5))
    a=0.0_xp
  END SUBROUTINE allocate_array_xp5

  SUBROUTINE allocate_array_xp6(a,is1,ie1,is2,ie2,is3,ie3,is4,ie4,is5,ie5,is6,ie6)
    IMPLICIT NONE
    real(xp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2,is3,ie3,is4,ie4,is5,ie5,is6,ie6
    ALLOCATE(a(is1:ie1,is2:ie2,is3:ie3,is4:ie4,is5:ie5,is6:ie6))
    a=0.0_xp
  END SUBROUTINE allocate_array_xp6

  SUBROUTINE allocate_array_xp7(a,is1,ie1,is2,ie2,is3,ie3,is4,ie4,is5,ie5,is6,ie6,is7,ie7)
    IMPLICIT NONE
    REAL(xp), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2,is3,ie3,is4,ie4,is5,ie5,is6,ie6,is7,ie7
    ALLOCATE(a(is1:ie1,is2:ie2,is3:ie3,is4:ie4,is5:ie5,is6:ie6,is7:ie7))
    a=0.0_xp
  END SUBROUTINE allocate_array_xp7
  !========================================

  SUBROUTINE allocate_array_dc1(a,is1,ie1)
    IMPLICIT NONE
    COMPLEX(xp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1
    ALLOCATE(a(is1:ie1))
    a=CMPLX(0.0_xp,0.0_xp)
  END SUBROUTINE allocate_array_dc1

  SUBROUTINE allocate_array_dc2(a,is1,ie1,is2,ie2)
    IMPLICIT NONE
    COMPLEX(xp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2
    ALLOCATE(a(is1:ie1,is2:ie2))
    a=CMPLX(0.0_xp,0.0_xp)
  END SUBROUTINE allocate_array_dc2

  SUBROUTINE allocate_array_dc3(a,is1,ie1,is2,ie2,is3,ie3)
    IMPLICIT NONE
    COMPLEX(xp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2,is3,ie3
    ALLOCATE(a(is1:ie1,is2:ie2,is3:ie3))
    a=CMPLX(0.0_xp,0.0_xp)
  END SUBROUTINE allocate_array_dc3

  SUBROUTINE allocate_array_dc4(a,is1,ie1,is2,ie2,is3,ie3,is4,ie4)
    IMPLICIT NONE
    COMPLEX(xp), DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2,is3,ie3,is4,ie4
    ALLOCATE(a(is1:ie1,is2:ie2,is3:ie3,is4:ie4))
    a=CMPLX(0.0_xp,0.0_xp)
  END SUBROUTINE allocate_array_dc4

  SUBROUTINE allocate_array_dc5(a,is1,ie1,is2,ie2,is3,ie3,is4,ie4,is5,ie5)
    IMPLICIT NONE
    COMPLEX(xp), DIMENSION(:,:,:,:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2,is3,ie3,is4,ie4,is5,ie5
    ALLOCATE(a(is1:ie1,is2:ie2,is3:ie3,is4:ie4,is5:ie5))
    a=CMPLX(0.0_xp,0.0_xp)
  END SUBROUTINE allocate_array_dc5

  SUBROUTINE allocate_array_dc6(a,is1,ie1,is2,ie2,is3,ie3,is4,ie4,is5,ie5,is6,ie6)
    IMPLICIT NONE
    COMPLEX(xp), DIMENSION(:,:,:,:,:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2,is3,ie3,is4,ie4,is5,ie5,is6,ie6
    ALLOCATE(a(is1:ie1,is2:ie2,is3:ie3,is4:ie4,is5:ie5,is6:ie6))
    a=CMPLX(0.0_xp,0.0_xp)
  END SUBROUTINE allocate_array_dc6

  SUBROUTINE allocate_array_dc7(a,is1,ie1,is2,ie2,is3,ie3,is4,ie4,is5,ie5,is6,ie6,is7,ie7)
    IMPLICIT NONE
    COMPLEX(xp), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2,is3,ie3,is4,ie4,is5,ie5,is6,ie6,is7,ie7
    ALLOCATE(a(is1:ie1,is2:ie2,is3:ie3,is4:ie4,is5:ie5,is6:ie6,is7:ie7))
    a=CMPLX(0.0_xp,0.0_xp)
  END SUBROUTINE allocate_array_dc7
  !========================================

  SUBROUTINE allocate_array_i1(a,is1,ie1)
    IMPLICIT NONE
    INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1
    ALLOCATE(a(is1:ie1))
    a=0
  END SUBROUTINE allocate_array_i1

  SUBROUTINE allocate_array_i2(a,is1,ie1,is2,ie2)
    IMPLICIT NONE
    INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2
    ALLOCATE(a(is1:ie1,is2:ie2))
    a=0
  END SUBROUTINE allocate_array_i2

  SUBROUTINE allocate_array_i3(a,is1,ie1,is2,ie2,is3,ie3)
    IMPLICIT NONE
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2,is3,ie3
    ALLOCATE(a(is1:ie1,is2:ie2,is3:ie3))
    a=0
  END SUBROUTINE allocate_array_i3

  SUBROUTINE allocate_array_i4(a,is1,ie1,is2,ie2,is3,ie3,is4,ie4)
    IMPLICIT NONE
    INTEGER, DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2,is3,ie3,is4,ie4
    ALLOCATE(a(is1:ie1,is2:ie2,is3:ie3,is4:ie4))
    a=0
  END SUBROUTINE allocate_array_i4

  SUBROUTINE allocate_array_i5(a,is1,ie1,is2,ie2,is3,ie3,is4,ie4,is5,ie5)
    IMPLICIT NONE
    INTEGER, DIMENSION(:,:,:,:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2,is3,ie3,is4,ie4,is5,ie5
    ALLOCATE(a(is1:ie1,is2:ie2,is3:ie3,is4:ie4,is5:ie5))
    a=0
  END SUBROUTINE allocate_array_i5

  !========================================

  SUBROUTINE allocate_array_l1(a,is1,ie1)
    IMPLICIT NONE
    LOGICAL, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1
    ALLOCATE(a(is1:ie1))
    a=.false.
  END SUBROUTINE allocate_array_l1

  SUBROUTINE allocate_array_l2(a,is1,ie1,is2,ie2)
    IMPLICIT NONE
    LOGICAL, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2
    ALLOCATE(a(is1:ie1,is2:ie2))
    a=.false.
  END SUBROUTINE allocate_array_l2

  SUBROUTINE allocate_array_l3(a,is1,ie1,is2,ie2,is3,ie3)
    IMPLICIT NONE
    LOGICAL, DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2,is3,ie3
    ALLOCATE(a(is1:ie1,is2:ie2,is3:ie3))
    a=.false.
  END SUBROUTINE allocate_array_l3

  SUBROUTINE allocate_array_l4(a,is1,ie1,is2,ie2,is3,ie3,is4,ie4)
    IMPLICIT NONE
    LOGICAL, DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2,is3,ie3,is4,ie4
    ALLOCATE(a(is1:ie1,is2:ie2,is3:ie3,is4:ie4))
    a=.false.
  END SUBROUTINE allocate_array_l4

  SUBROUTINE allocate_array_l5(a,is1,ie1,is2,ie2,is3,ie3,is4,ie4,is5,ie5)
    IMPLICIT NONE
    LOGICAL, DIMENSION(:,:,:,:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2,is3,ie3,is4,ie4,is5,ie5
    ALLOCATE(a(is1:ie1,is2:ie2,is3:ie3,is4:ie4,is5:ie5))
    a=.false.
  END SUBROUTINE allocate_array_l5

END MODULE basic
