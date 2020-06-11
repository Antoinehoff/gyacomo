MODULE fourier_grid
  ! Grid module for spatial discretization

  USE prec_const
  IMPLICIT NONE
  PRIVATE

  !   GRID Namelist
  INTEGER,  PUBLIC, PROTECTED ::  pmax=3         ! The maximal Hermite-moment computed
  INTEGER,  PUBLIC, PROTECTED ::  jmax=2         ! The maximal Laguerre-moment computed
  INTEGER,  PUBLIC, PROTECTED ::  nkr=1          ! Number of total internal grid points in kr
  INTEGER,  PUBLIC, PROTECTED ::  nkz=1          ! Number of total internal grid points in kz
  REAL(dp), PUBLIC, PROTECTED :: krmin=0._dp    ! kr coordinate for left boundary
  REAL(dp), PUBLIC, PROTECTED :: krmax=1._dp    ! kr coordinate for right boundary  
  REAL(dp), PUBLIC, PROTECTED :: kzmin=0._dp    ! kz coordinate for left boundary
  REAL(dp), PUBLIC, PROTECTED :: kzmax=1._dp    ! kz coordinate for right boundary


  ! Indices of s -> start , e-> end
  INTEGER, PUBLIC, PROTECTED ::  ips,  ipe,  ijs,  ije
  INTEGER, PUBLIC, PROTECTED :: ikrs, ikre, ikzs, ikze

  ! Toroidal direction
  REAL(dp), PUBLIC, PROTECTED ::    deltakr,    deltakz
  real(dp), PUBLIC, PROTECTED ::   deltakri,   deltakzi
  real(dp), PUBLIC, PROTECTED ::  deltakrih,  deltakzih
  real(dp), PUBLIC, PROTECTED :: deltakrisq, deltakzisq

  ! Grids containing position
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC, PROTECTED :: parray
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC, PROTECTED :: zarray
  
  INTEGER, PUBLIC :: neq ! Number of equations in poisson solver (number of rows/cols of laplace matrix)
  
  ! Public Functions
  PUBLIC :: set_pgrid, set_zgrid
  PUBLIC :: space_grid_readinputs, space_grid_outputinputs

contains

  subroutine set_pjgrid
    ! Initialize pj grid (Hermite moment) parameters

    use prec_const
    implicit none
    integer :: ip, ij
    !since N_e^00=1 and N_e^10 = N_e^20 = = N_e^01 = 0
    ips=0
    ijs=0
    ipe=pmax
    jpe=jmax

    allocate(pjarray(ips:ipe,jps:jpe))
    DO ip = ips,ipe
      DO ij = ijs,ije
        pjarray(ip,ij) = ip ! simply contains the moment number
      END DO
    END DO

  end subroutine set_pjgrid


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


  SUBROUTINE fourier_grid_readinputs
    ! Read the input parameters

    USE basic, ONLY : lu_in

    USE prec_const
    IMPLICIT NONE

    NAMELIST /GRID/ pmax, nz, zmin, zmax

    READ(lu_in,grid)
    WRITE(*,grid)

  END SUBROUTINE fourier_grid_readinputs


  SUBROUTINE fourier_grid_outputinputs(fidres, str)
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

   END SUBROUTINE fourier_grid_outputinputs


END MODULE fourier_grid
