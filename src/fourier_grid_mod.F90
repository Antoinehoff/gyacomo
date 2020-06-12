MODULE fourier_grid
  ! Grid module for spatial discretization

  USE prec_const
  IMPLICIT NONE
  PRIVATE

  !   GRID Namelist
  INTEGER,  PUBLIC, PROTECTED :: pmaxe = 2      ! The maximal electron Hermite-moment computed
  INTEGER,  PUBLIC, PROTECTED :: jmaxe = 2      ! The maximal electron Laguerre-moment computed  
  INTEGER,  PUBLIC, PROTECTED :: pmaxi = 2      ! The maximal ion Hermite-moment computed
  INTEGER,  PUBLIC, PROTECTED :: jmaxi = 2      ! The maximal ion Laguerre-moment computed
  INTEGER,  PUBLIC, PROTECTED :: nkr   = 10     ! Number of total internal grid points in kr
  REAL(dp), PUBLIC, PROTECTED :: krmin = 0._dp  ! kr coordinate for left boundary
  REAL(dp), PUBLIC, PROTECTED :: krmax = 1._dp  ! kr coordinate for right boundary
  INTEGER,  PUBLIC, PROTECTED :: nkz   = 10     ! Number of total internal grid points in kz
  REAL(dp), PUBLIC, PROTECTED :: kzmin = 0._dp  ! kz coordinate for left boundary
  REAL(dp), PUBLIC, PROTECTED :: kzmax = 1._dp  ! kz coordinate for right boundary

  ! Indices of s -> start , e-> end
  INTEGER, PUBLIC, PROTECTED ::  ipjs, ipje
  INTEGER, PUBLIC, PROTECTED ::  Nmome, Nmomi, Nmomtot
  INTEGER, PUBLIC, PROTECTED ::  ikrs, ikre, ikzs, ikze

  ! Toroidal direction
  REAL(dp), PUBLIC, PROTECTED ::  deltakr
  REAL(dp), PUBLIC, PROTECTED ::  deltakz

  ! Grids containing position in fourier space
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC, PROTECTED :: kzarray
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC, PROTECTED :: krarray
  
  ! Grid containing the polynomials degrees
  INTEGER,  DIMENSION(:), ALLOCATABLE, PUBLIC, PROTECTED :: pjarray

  ! Public Functions
  PUBLIC :: set_krgrid, set_kzgrid, set_pj
  PUBLIC :: fourier_grid_readinputs, fourier_grid_outputinputs
  PUBLIC :: bare, bari
  PUBLIC :: rabe, rabi

contains

  subroutine set_krgrid
    use prec_const
    implicit none
    integer :: ikr  
    ! Start and end indices of grid
    ikrs = 1
    ikre = nkr    
    ! Grid spacings, precompute some inverses
    deltakr = (krmax - krmin) / real(nkr,dp)
    ! Discretized kr positions
    allocate(krarray(ikrs:ikre))
    DO ikr = ikrs,ikre
      krarray(ikr) = krmin + real(ikr-1,dp) * deltakr
    END DO
  end subroutine set_krgrid

  subroutine set_kzgrid
    use prec_const
    implicit none
    integer :: ikz
    ! Start and end indices of grid
    ikzs = 1
    ikze = nkz    
    ! Grid spacings, precompute some inverses
    deltakz = (kzmax - kzmin) / real(nkz,dp)
    ! Discretized kz positions
    allocate(kzarray(ikzs:ikze))
    DO ikz = ikzs,ikze
       kzarray(ikz) = kzmin + real(ikz-1,dp) * deltakz
    END DO
  end subroutine set_kzgrid

  subroutine set_pj
    use prec_const
    implicit none
    integer :: ipj
    ! number of electrons moments
    Nmome   = (Pmaxe + 1)*(Jmaxe + 1)
    ! number of ions moments
    Nmomi   = (Pmaxi + 1)*(Jmaxi + 1)
    ! total number of moments
    Nmomtot = Nmome + Nmomi
    ipjs = 1
    ipje = Nmomtot
    ! Polynomials degrees pj = (Jmaxs + 1)*p + j + 1
    allocate(pjarray(ipjs:ipje))
    DO ipj = ipjs,ipje
      pjarray(ipj) = ipj
    END DO
  end subroutine set_pj

  SUBROUTINE fourier_grid_readinputs
    ! Read the input parameters

    USE basic, ONLY : lu_in

    USE prec_const
    IMPLICIT NONE

    NAMELIST /GRID/ pmaxe, jmaxe, pmaxi, jmaxi, &
                    nkr, krmin, krmax, nkz, kzmin, kzmax

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
    CALL attach(fidres, TRIM(str), "pmaxe", pmaxe)
    CALL attach(fidres, TRIM(str), "jmaxe", jmaxe)
    CALL attach(fidres, TRIM(str), "pmaxi", pmaxi)
    CALL attach(fidres, TRIM(str), "jmaxi", jmaxi)
    CALL attach(fidres, TRIM(str),   "nkr",   nkr)
    CALL attach(fidres, TRIM(str), "krmin", krmin)
    CALL attach(fidres, TRIM(str), "krmax", krmax)
    CALL attach(fidres, TRIM(str),   "nkz",   nkz)
    CALL attach(fidres, TRIM(str), "kzmin", kzmin)
    CALL attach(fidres, TRIM(str), "kzmax", kzmax)
  END SUBROUTINE fourier_grid_outputinputs

  SUBROUTINE bare(p,j,idx)
    USE prec_const
    IMPLICIT NONE
    INTEGER, INTENT(in) :: p,j
    INTEGER, INTENT(out):: idx

    idx = (jmaxe + 1)*p + j + 1

  END SUBROUTINE bare

  SUBROUTINE bari(p,j,idx)
    INTEGER, INTENT(in) :: p,j
    INTEGER, INTENT(out):: idx

    idx = (jmaxi + 1)*p + j + 1

  END SUBROUTINE bari

  SUBROUTINE rabe(idx, p, j)
    INTEGER, INTENT(in) :: idx
    INTEGER, INTENT(out):: p,j
    p = FLOOR(real((idx-1) / (jmaxe + 1)))
    j = idx - p * (jmaxe+1)
  END SUBROUTINE rabe

  SUBROUTINE rabi(idx, p, j)
    INTEGER, INTENT(in):: idx
    INTEGER, INTENT(out) :: p,j
    p = FLOOR(real((idx-Nmome - 1) / (jmaxi + 1)))
    j = (idx-Nmome) - p * (jmaxi+1)
  END SUBROUTINE rabi

END MODULE fourier_grid
