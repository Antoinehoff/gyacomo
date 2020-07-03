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

  ! Indices of s -> start , e-> END
  INTEGER, PUBLIC, PROTECTED ::  ips_e,ipe_e, ijs_e,ije_e
  INTEGER, PUBLIC, PROTECTED ::  ips_i,ipe_i, ijs_i,ije_i
  INTEGER, PUBLIC, PROTECTED ::  ns_e,ne_e, ns_i,ne_i     !start/end indices for coll. mat.
  INTEGER, PUBLIC, PROTECTED ::  ikrs, ikre, ikzs, ikze

  ! Toroidal direction
  REAL(dp), PUBLIC, PROTECTED ::  deltakr
  REAL(dp), PUBLIC, PROTECTED ::  deltakz

  ! Grids containing position in fourier space
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC, PROTECTED :: kzarray
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC, PROTECTED :: krarray
  
  ! Grid containing the polynomials degrees
  INTEGER,  DIMENSION(:), ALLOCATABLE, PUBLIC, PROTECTED :: parray_e
  INTEGER,  DIMENSION(:), ALLOCATABLE, PUBLIC, PROTECTED :: parray_i
  INTEGER,  DIMENSION(:), ALLOCATABLE, PUBLIC, PROTECTED :: jarray_e
  INTEGER,  DIMENSION(:), ALLOCATABLE, PUBLIC, PROTECTED :: jarray_i

  ! Public Functions
  PUBLIC :: set_krgrid, set_kzgrid, set_pj
  PUBLIC :: fourier_grid_readinputs, fourier_grid_outputinputs
  PUBLIC :: bare, bari

CONTAINS

  SUBROUTINE set_krgrid
    USE prec_const
    IMPLICIT NONE
    INTEGER :: ikr  
    ! Start and END indices of grid
    ikrs = 1
    ikre = nkr    
    ! Grid spacings, precompute some inverses
    deltakr = (krmax - krmin) / REAL(nkr,dp)
    ! Discretized kr positions
    ALLOCATE(krarray(ikrs:ikre))
    DO ikr = ikrs,ikre
      krarray(ikr) = krmin + REAL(ikr-1,dp) * deltakr
    END DO
  END SUBROUTINE set_krgrid

  SUBROUTINE set_kzgrid
    USE prec_const
    IMPLICIT NONE
    INTEGER :: ikz
    ! Start and END indices of grid
    ikzs = 1
    ikze = nkz    
    ! Grid spacings, precompute some inverses
    deltakz = (kzmax - kzmin) / REAL(nkz,dp)
    ! Discretized kz positions
    ALLOCATE(kzarray(ikzs:ikze))
    DO ikz = ikzs,ikze
       kzarray(ikz) = kzmin + REAL(ikz-1,dp) * deltakz
    END DO
  END SUBROUTINE set_kzgrid

  SUBROUTINE set_pj
    USE prec_const
    IMPLICIT NONE
    INTEGER :: ip, ij
    ips_e = 1; ipe_e = pmaxe + 1
    ips_i = 1; ipe_i = pmaxi + 1    
    ALLOCATE(parray_e(ips_e:ipe_e))
    ALLOCATE(parray_i(ips_i:ipe_i))
    DO ip = ips_e,ipe_e; parray_e(ip) = ip-1; END DO
    DO ip = ips_i,ipe_i; parray_i(ip) = ip-1; END DO

    ijs_e = 1; ije_e = jmaxe + 1
    ijs_i = 1; ije_i = jmaxi + 1
    ALLOCATE(jarray_e(ijs_e:ije_e))
    ALLOCATE(jarray_i(ijs_i:ije_i))
    DO ij = ijs_e,ije_e; jarray_e(ij) = ij-1; END DO
    DO ij = ijs_i,ije_i; jarray_i(ij) = ij-1; END DO

    ns_e = bare(ips_e,ijs_e); ne_e = bare(ipe_e,ije_e)
    ns_i = bari(ips_i,ijs_i); ne_i = bari(ipe_i,ije_i)

  END SUBROUTINE set_pj

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

  
  FUNCTION bare(ip,ij)
    IMPLICIT NONE
    INTEGER :: bare, ip, ij
    bare = (jmaxe+1)*ip + ij + 1
  END FUNCTION

  FUNCTION bari(ip,ij)
    IMPLICIT NONE
    INTEGER :: bari, ip, ij
    bari = bare(pmaxe,jmaxe) + (jmaxi+1)*ip + ij + 1
  END FUNCTION


END MODULE fourier_grid
