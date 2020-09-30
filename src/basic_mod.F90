MODULE basic
  !   Basic module for time dependent problems

  use prec_const
  IMPLICIT none

  INTEGER  :: nrun   = 1           ! Number of time steps to run
  real(dp) :: tmax   = 100000.0    ! Maximum simulation time
  real(dp) :: dt     = 1.0         ! Time step
  real(dp) :: time   = 0           ! Current simulation time (Init from restart file)

  INTEGER :: jobnum  = 0           ! Job number
  INTEGER :: step    = 0           ! Calculation step of this run
  INTEGER :: cstep   = 0           ! Current step number (Init from restart file)
  LOGICAL :: RESTART = .FALSE.     ! Signal end of run
  LOGICAL :: nlend   = .FALSE.     ! Signal end of run

  INTEGER :: ierr                  ! flag for MPI error
  INTEGER :: iframe1d              ! counting the number of times 1d datasets are outputed (for diagnose)
  INTEGER :: iframe2d              ! counting the number of times 2d datasets are outputed (for diagnose)
  INTEGER :: iframe3d              ! counting the number of times 3d datasets are outputed (for diagnose)
  INTEGER :: iframe5d              ! counting the number of times 5d datasets are outputed (for diagnose)

  !  List of logical file units
  INTEGER :: lu_in   = 90              ! File duplicated from STDIN
  INTEGER :: lu_job  = 91              ! myjob file

  real :: start, finish ! To measure computation time

  INTERFACE allocate_array
    MODULE PROCEDURE allocate_array_dp1,allocate_array_dp2,allocate_array_dp3,allocate_array_dp4
    MODULE PROCEDURE allocate_array_dc1,allocate_array_dc2,allocate_array_dc3,allocate_array_dc4, allocate_array_dc5
    MODULE PROCEDURE allocate_array_i1,allocate_array_i2,allocate_array_i3,allocate_array_i4
    MODULE PROCEDURE allocate_array_l1,allocate_array_l2,allocate_array_l3,allocate_array_l4
  END INTERFACE allocate_array

CONTAINS

  !================================================================================
  SUBROUTINE basic_data
    !   Read basic data for input file

    use prec_const
    IMPLICIT NONE

    NAMELIST /BASIC/  nrun, dt, tmax, RESTART

    READ(lu_in,basic)

  END SUBROUTINE basic_data
  !================================================================================
  SUBROUTINE daytim(str)
    !   Print date and time

    use prec_const
    IMPLICIT NONE

    CHARACTER(len=*) , INTENT(in) :: str
    CHARACTER(len=16) :: d, t, dat, time
    !________________________________________________________________________________
    !
    CALL DATE_AND_TIME(d,t)
    dat=d(7:8) // '/' // d(5:6) // '/' // d(1:4)
    time=t(1:2) // ':' // t(3:4) // ':' // t(5:10)
    WRITE(*,'(a,1x,a,1x,a)') str, dat(1:10), time(1:12)
    !
  END SUBROUTINE daytim
  !================================================================================

! To allocate arrays of doubles, integers, etc. at run time

  SUBROUTINE allocate_array_dp1(a,is1,ie1)
    IMPLICIT NONE
    real(dp), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1
    ALLOCATE(a(is1:ie1))
    a=0.0_dp
  END SUBROUTINE allocate_array_dp1

  SUBROUTINE allocate_array_dp2(a,is1,ie1,is2,ie2)
    IMPLICIT NONE
    real(dp), DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2
    ALLOCATE(a(is1:ie1,is2:ie2))
    a=0.0_dp
  END SUBROUTINE allocate_array_dp2

  SUBROUTINE allocate_array_dp3(a,is1,ie1,is2,ie2,is3,ie3)
    IMPLICIT NONE
    real(dp), DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2,is3,ie3
    ALLOCATE(a(is1:ie1,is2:ie2,is3:ie3))
    a=0.0_dp
  END SUBROUTINE allocate_array_dp3

  SUBROUTINE allocate_array_dp4(a,is1,ie1,is2,ie2,is3,ie3,is4,ie4)
    IMPLICIT NONE
    real(dp), DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2,is3,ie3,is4,ie4
    ALLOCATE(a(is1:ie1,is2:ie2,is3:ie3,is4:ie4))
    a=0.0_dp
  END SUBROUTINE allocate_array_dp4

  SUBROUTINE allocate_array_dp5(a,is1,ie1,is2,ie2,is3,ie3,is4,ie4,is5,ie5)
    IMPLICIT NONE
    real(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2,is3,ie3,is4,ie4,is5,ie5
    ALLOCATE(a(is1:ie1,is2:ie2,is3:ie3,is4:ie4,is5:ie5))
    a=0.0_dp
  END SUBROUTINE allocate_array_dp5

  !========================================

  SUBROUTINE allocate_array_dc1(a,is1,ie1)
    IMPLICIT NONE
    DOUBLE COMPLEX, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1
    ALLOCATE(a(is1:ie1))
    a=CMPLX(0.0_dp,0.0_dp)
  END SUBROUTINE allocate_array_dc1

  SUBROUTINE allocate_array_dc2(a,is1,ie1,is2,ie2)
    IMPLICIT NONE
    DOUBLE COMPLEX, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2
    ALLOCATE(a(is1:ie1,is2:ie2))
    a=CMPLX(0.0_dp,0.0_dp)
  END SUBROUTINE allocate_array_dc2

  SUBROUTINE allocate_array_dc3(a,is1,ie1,is2,ie2,is3,ie3)
    IMPLICIT NONE
    DOUBLE COMPLEX, DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2,is3,ie3
    ALLOCATE(a(is1:ie1,is2:ie2,is3:ie3))
    a=CMPLX(0.0_dp,0.0_dp)
  END SUBROUTINE allocate_array_dc3

  SUBROUTINE allocate_array_dc4(a,is1,ie1,is2,ie2,is3,ie3,is4,ie4)
    IMPLICIT NONE
    DOUBLE COMPLEX, DIMENSION(:,:,:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2,is3,ie3,is4,ie4
    ALLOCATE(a(is1:ie1,is2:ie2,is3:ie3,is4:ie4))
    a=CMPLX(0.0_dp,0.0_dp)
  END SUBROUTINE allocate_array_dc4

  SUBROUTINE allocate_array_dc5(a,is1,ie1,is2,ie2,is3,ie3,is4,ie4,is5,ie5)
    IMPLICIT NONE
    DOUBLE COMPLEX, DIMENSION(:,:,:,:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2,is3,ie3,is4,ie4,is5,ie5
    ALLOCATE(a(is1:ie1,is2:ie2,is3:ie3,is4:ie4,is5:ie5))
    a=CMPLX(0.0_dp,0.0_dp)
  END SUBROUTINE allocate_array_dc5

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
    real(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE, INTENT(INOUT) :: a
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
    real(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE, INTENT(INOUT) :: a
    INTEGER, INTENT(IN) :: is1,ie1,is2,ie2,is3,ie3,is4,ie4,is5,ie5
    ALLOCATE(a(is1:ie1,is2:ie2,is3:ie3,is4:ie4,is5:ie5))
    a=.false.
  END SUBROUTINE allocate_array_l5

END MODULE basic
