MODULE prec_const

  use mpi

  use, intrinsic :: iso_fortran_env, only: REAL32, REAL64, &
                                           stdin=>input_unit, &
                                           stdout=>output_unit, &
                                           stderr=>error_unit
  use, intrinsic :: iso_c_binding

  ! Define single and double precision
  INTEGER, PARAMETER :: sp = REAL32 !Single precision
  INTEGER, PARAMETER :: dp = REAL64 !Double precision  
  INTEGER, private :: dp_r, dp_p !Range and Aprecision of doubles
  INTEGER, private :: sp_r, sp_p !Range and precision of singles
  INTEGER, private :: MPI_SP !Single precision for MPI
  INTEGER, private :: MPI_DP !Double precision in MPI
  INTEGER, private :: MPI_SUM_DP !Sum reduction operation for DP datatype
  INTEGER, private :: MPI_MAX_DP !Max reduction operation for DP datatype
  INTEGER, private :: MPI_MIN_DP !Min reduction operation for DP datatype

  ! Define a generic precision parameter for the entire program
#ifdef SINGLE_PRECISION
    INTEGER, PARAMETER :: xp       = REAL32
    INTEGER, PARAMETER :: c_xp_c   = C_FLOAT_COMPLEX
    INTEGER, PARAMETER :: c_xp_r   = C_FLOAT
    INTEGER, PARAMETER :: mpi_xp_r = MPI_FLOAT
    INTEGER, PARAMETER :: mpi_xp_c = MPI_COMPLEX
#else
    INTEGER, PARAMETER :: xp       = REAL64
    INTEGER, PARAMETER :: c_xp_c   = C_DOUBLE_COMPLEX
    INTEGER, PARAMETER :: c_xp_r   = C_DOUBLE
    INTEGER, PARAMETER :: mpi_xp_r = MPI_DOUBLE
    INTEGER, PARAMETER :: mpi_xp_c = MPI_DOUBLE_COMPLEX
#endif
  ! Auxiliary variables (unused)
  INTEGER, private   :: xp_r, xp_p !Range and precision of single
  INTEGER, private   :: MPI_XP     !Double precision in MPI
  INTEGER, private   :: MPI_SUM_XP !Sum reduction operation for xp datatype
  INTEGER, private   :: MPI_MAX_XP !Max reduction operation for xp datatype
  INTEGER, private   :: MPI_MIN_XP !Min reduction operation for xp datatype


  ! Some useful constants, to avoid recomputing them too often
  REAL(xp),    PARAMETER :: PI=3.141592653589793238462643383279502884197_xp
  REAL(xp),    PARAMETER :: PIO2=1.57079632679489661923132169163975144209858_xp
  REAL(xp),    PARAMETER :: TWOPI=6.283185307179586476925286766559005768394_xp
  REAL(xp),    PARAMETER :: ONEOPI=0.3183098861837906912164442019275156781077_xp
  REAL(xp),    PARAMETER :: SQRT2=1.41421356237309504880168872420969807856967_xp
  REAL(xp),    PARAMETER :: SQRT6=SQRT(6._xp)
  REAL(xp),    PARAMETER :: INVSQRT2=0.7071067811865475244008443621048490392848359377_xp
  REAL(xp),    PARAMETER :: SQRT3=1.73205080756887729352744634150587236694281_xp
  REAL(xp),    PARAMETER :: onetwelfth=0.08333333333333333333333333333333333333333333333_xp
  REAL(xp),    PARAMETER :: onetwentyfourth=0.04166666666666666666666666666666666666666666667_xp
  REAL(xp),    PARAMETER :: onethird=0.33333333333333333333333333333333333333333333333_xp
  REAL(xp),    PARAMETER :: twothird=0.66666666666666666666666666666666666666666666666_xp
  REAL(xp),    PARAMETER :: onesixth=0.1666666666666666666666666666666666666666666667_xp
  REAL(xp),    PARAMETER :: fivesixths=0.8333333333333333333333333333333333333333333333_xp
  REAL(xp),    PARAMETER :: sevensixths=1.1666666666666666666666666666666666666666666667_xp
  REAL(xp),    PARAMETER :: elevensixths=1.833333333333333333333333333333333333333333333_xp
  REAL(xp),    PARAMETER :: nineeighths=1.125_xp
  REAL(xp),    PARAMETER :: onesixteenth=0.0625_xp
  REAL(xp),    PARAMETER :: ninesixteenths=0.5625_xp
  REAL(xp),    PARAMETER :: thirteentwelfths = 1.083333333333333333333333333333333333333333333_xp
  COMPLEX(xp), PARAMETER :: imagu = (0._xp,1._xp)
  CONTAINS

    SUBROUTINE INIT_PREC_CONST

      IMPLICIT NONE
      integer :: ierr,me

      ! REAL(sp) :: a = 1_sp
      ! REAL(dp) :: b = 1_dp
      !Get range and precision of ISO FORTRAN sizes
      ! sp_r = range(a)
      ! sp_p = precision(a)
      ! dp_r = range(b)
      ! dp_p = precision(b)
      
      REAL(xp) :: c = 1_xp
      xp_r = range(c)
      xp_p = precision(c)

      CALL mpi_comm_rank(MPI_COMM_WORLD,me,ierr)

      !Create MPI datatypes that support the specific size
      ! CALL MPI_Type_create_f90_real(sp_p,sp_r,MPI_sp,ierr)
      ! CALL MPI_Type_create_f90_real(dp_p,dp_r,MPI_xp,ierr)
      CALL MPI_Type_create_f90_real(xp_p,xp_r,MPI_xp,ierr)

    END SUBROUTINE INIT_PREC_CONST

END MODULE prec_const
