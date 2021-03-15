MODULE prec_const

  use mpi

  use, intrinsic :: iso_fortran_env, only: REAL32, REAL64, &
                                           stdin=>input_unit, &
                                           stdout=>output_unit, &
                                           stderr=>error_unit

  ! Precision for real and complex
  INTEGER, PARAMETER :: sp = REAL32 !Single precision, should not be used
  INTEGER, PARAMETER :: dp = REAL64 !real(dp), enforced through the code

  INTEGER, private :: dp_r, dp_p !Range and Aprecision of doubles
  INTEGER, private :: sp_r, sp_p !Range and precision of singles


  INTEGER, private :: MPI_SP !Single precision for MPI
  INTEGER, private :: MPI_DP !Double precision in MPI
  INTEGER, private :: MPI_SUM_DP !Sum reduction operation for DP datatype
  INTEGER, private :: MPI_MAX_DP !Max reduction operation for DP datatype
  INTEGER, private :: MPI_MIN_DP !Min reduction operation for DP datatype


  ! Some useful constants, to avoid recomputing then too often
  REAL(dp),    PARAMETER :: PI=3.141592653589793238462643383279502884197_dp
  REAL(dp),    PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_dp
  REAL(dp),    PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_dp
  REAL(dp),    PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_dp
  REAL(dp),    PARAMETER :: INVSQRT2=0.7071067811865475244008443621048490392848359377_dp
  REAL(dp),    PARAMETER :: SQRT3=1.73205080756887729352744634150587236694281_dp
  REAL(dp),    PARAMETER :: onetwelfth=0.08333333333333333333333333333333333333333333333_dp
  REAL(dp),    PARAMETER :: onetwentyfourth=0.04166666666666666666666666666666666666666666667_dp
  REAL(dp),    PARAMETER :: onethird=0.33333333333333333333333333333333333333333333333_dp
  REAL(dp),    PARAMETER :: onesixth=0.1666666666666666666666666666666666666666666667_dp
  REAL(dp),    PARAMETER :: fivesixths=0.8333333333333333333333333333333333333333333333_dp
  REAL(dp),    PARAMETER :: sevensixths=1.1666666666666666666666666666666666666666666667_dp
  REAL(dp),    PARAMETER :: elevensixths=1.833333333333333333333333333333333333333333333_dp
  REAL(dp),    PARAMETER :: nineeighths=1.125_dp
  REAL(dp),    PARAMETER :: onesixteenth=0.0625_dp
  REAL(dp),    PARAMETER :: ninesixteenths=0.5625_dp
  REAL(dp),    PARAMETER :: thirteentwelfths = 1.083333333333333333333333333333333333333333333_dp
  COMPLEX(dp), PARAMETER :: imagu = (0._dp,1._dp)
  CONTAINS

    SUBROUTINE INIT_PREC_CONST

      IMPLICIT NONE
      integer :: ierr,me

      REAL(sp) :: a = 1_sp
      REAL(dp) :: b = 1_dp

      !Get range and precision of ISO FORTRAN sizes
      sp_r = range(a)
      sp_p = precision(a)

      dp_r = range(b)
      dp_p = precision(b)

      CALL mpi_comm_rank(MPI_COMM_WORLD,me,ierr)

      !Create MPI datatypes that support the specific size
      CALL MPI_Type_create_f90_real(sp_p,sp_r,MPI_sp,ierr)
      CALL MPI_Type_create_f90_real(dp_p,dp_r,MPI_dp,ierr)

    END SUBROUTINE INIT_PREC_CONST

END MODULE prec_const
