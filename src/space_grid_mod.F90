MODULE space_grid
  ! Grid module for spatial discretization

  USE prec_const
  IMPLICIT NONE
  PRIVATE

  !   GRID Namelist
  INTEGER, PUBLIC, PROTECTED :: pmax=3         ! The maximal Hermite-moment computed
  INTEGER, PUBLIC, PROTECTED :: nz=1           ! Number of total internal grid points in z
  REAL(dp), PUBLIC, PROTECTED :: zmin=0._dp    ! z coordinate for left boundary
  REAL(dp), PUBLIC, PROTECTED :: zmax=1._dp    ! z coordinate for right boundary


  ! Indices of s -> start , e-> end
  INTEGER, PUBLIC, PROTECTED :: ips, ipe
  INTEGER, PUBLIC, PROTECTED :: izs, ize

  ! Toroidal direction
  REAL(dp), PUBLIC, PROTECTED :: deltaz
  real(dp), PUBLIC, PROTECTED :: deltazi
  real(dp), PUBLIC, PROTECTED :: deltazih
  real(dp), PUBLIC, PROTECTED :: deltazisq

  ! Grids containing position
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC, PROTECTED :: parray
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC, PROTECTED :: zarray
  
  INTEGER, PUBLIC :: neq ! Number of equations in poisson solver (number of rows/cols of laplace matrix)
  
  ! Public Functions
  PUBLIC :: set_pgrid, set_zgrid
  PUBLIC :: space_grid_readinputs, space_grid_outputinputs

contains

  subroutine set_pgrid
    ! Initialize p grid (Hermite moment) parameters

    use prec_const
    implicit none
    integer :: ip

    ips=3 ! Only start at 3 since N_e^0=1 and N_e^1 = N_e^2 = 0 
    ipe=pmax

    allocate(parray(ips:ipe))
    DO ip = ips,ipe
       parray(ip) = ip ! simply contains the moment number
    END DO

  end subroutine set_pgrid


  subroutine set_zgrid
    ! Initialize z  grid

    use prec_const
    implicit none

    integer :: iz
    
    ! Start and end indices of grid
    izs=1
    ize=nz

    ! Grid spacings, precompute some inverses
    deltaz    = (zmax-zmin)/real(nz,dp) !zmax*2._dp*pi/real(nz,dp) 
    deltazi   = 1._dp/deltaz ! inverse
    deltazih   = 0.5_dp/deltaz ! half of inverse 
    deltazisq = 1._dp/deltaz/deltaz ! inverse squared

    ! Discretized z positions
    allocate(zarray(izs:ize))
    DO iz = izs,ize
       zarray(iz) = zmin + real(iz-1,dp)*deltaz
    END DO

  end subroutine set_zgrid


  SUBROUTINE space_grid_readinputs
    ! Read the input parameters

    USE basic, ONLY : lu_in

    USE prec_const
    IMPLICIT NONE

    NAMELIST /GRID/ pmax, nz, zmin, zmax

    READ(lu_in,grid)
    WRITE(*,grid)

  END SUBROUTINE space_grid_readinputs


  SUBROUTINE space_grid_outputinputs(fidres, str)
    ! Write the input parameters to the results_xx.h5 file

    USE futils, ONLY: attach

    USE prec_const
    IMPLICIT NONE

    INTEGER, INTENT(in) :: fidres
    CHARACTER(len=256), INTENT(in) :: str

     CALL attach(fidres, TRIM(str), "pmax", pmax)
     CALL attach(fidres, TRIM(str), "nz", nz)
     CALL attach(fidres, TRIM(str), "zmin", zmin)
     CALL attach(fidres, TRIM(str), "zmax", zmax)

   END SUBROUTINE space_grid_outputinputs


END MODULE space_grid
