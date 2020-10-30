MODULE grid
  ! Grid module for spatial discretization
  USE prec_const
  IMPLICIT NONE
  PRIVATE

  !   GRID Namelist
  INTEGER,  PUBLIC, PROTECTED :: pmaxe = 1      ! The maximal electron Hermite-moment computed
  INTEGER,  PUBLIC, PROTECTED :: jmaxe = 1      ! The maximal electron Laguerre-moment computed
  INTEGER,  PUBLIC, PROTECTED :: pmaxi = 1      ! The maximal ion Hermite-moment computed
  INTEGER,  PUBLIC, PROTECTED :: jmaxi = 1      ! The maximal ion Laguerre-moment computed
  INTEGER,  PUBLIC, PROTECTED :: maxj  = 1      ! The maximal ion Laguerre-moment computed
  INTEGER,  PUBLIC, PROTECTED :: Nr    = 16     ! Number of total internal grid points in r
  REAL(dp), PUBLIC, PROTECTED :: Lr    = 1._dp  ! horizontal length of the spatial box
  INTEGER,  PUBLIC, PROTECTED :: Nz    = 16     ! Number of total internal grid points in z
  REAL(dp), PUBLIC, PROTECTED :: Lz    = 1._dp  ! vertical length of the spatial box
  INTEGER,  PUBLIC, PROTECTED :: Nkr   = 8      ! Number of total internal grid points in kr
  REAL(dp), PUBLIC, PROTECTED :: Lkr   = 1._dp  ! horizontal length of the fourier box
  INTEGER,  PUBLIC, PROTECTED :: Nkz   = 16     ! Number of total internal grid points in kz
  REAL(dp), PUBLIC, PROTECTED :: Lkz   = 1._dp  ! vertical length of the fourier box
  REAL(dp), PUBLIC, PROTECTED :: kpar  = 0_dp   ! parallel wave vector component

  ! For Orszag filter
  REAL(dp), PUBLIC, PROTECTED :: two_third_krmax
  REAL(dp), PUBLIC, PROTECTED :: two_third_kzmax

  ! 1D Antialiasing arrays (2/3 rule)
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: AA_r
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: AA_z

  ! Grids containing position in physical space
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: rarray
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: zarray
  REAL(dp), PUBLIC, PROTECTED ::  deltar,  deltaz
  INTEGER,  PUBLIC, PROTECTED  ::  irs,  ire,  izs,  ize
  INTEGER,  PUBLIC :: ir,iz ! counters

  ! Grids containing position in fourier space
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: krarray
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: kzarray
  REAL(dp), PUBLIC, PROTECTED ::  deltakr, deltakz
  INTEGER,  PUBLIC, PROTECTED ::  ikrs, ikre, ikzs, ikze, a
  INTEGER,  PUBLIC, PROTECTED :: ikr_0, ikz_0 ! Indices of k-grid origin
  INTEGER,  PUBLIC :: ikr, ikz, ip, ij ! counters

  ! Grid containing the polynomials degrees
  INTEGER,  DIMENSION(:), ALLOCATABLE, PUBLIC :: parray_e
  INTEGER,  DIMENSION(:), ALLOCATABLE, PUBLIC :: parray_i
  INTEGER,  DIMENSION(:), ALLOCATABLE, PUBLIC :: jarray_e
  INTEGER,  DIMENSION(:), ALLOCATABLE, PUBLIC :: jarray_i
  INTEGER, PUBLIC, PROTECTED ::  ips_e,ipe_e, ijs_e,ije_e
  INTEGER, PUBLIC, PROTECTED ::  ips_i,ipe_i, ijs_i,ije_i

  ! Public Functions
  PUBLIC :: set_pj
  PUBLIC :: set_rgrid,  set_zgrid
  PUBLIC :: set_krgrid, set_kzgrid
  PUBLIC :: grid_readinputs, grid_outputinputs
  PUBLIC :: bare, bari
CONTAINS

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
    maxj  = MAX(jmaxi, jmaxe)
  END SUBROUTINE set_pj

  SUBROUTINE set_rgrid
    USE prec_const
    IMPLICIT NONE
    INTEGER :: ir
    ! Start and END indices of grid
    irs = 1
    ire = Nr
    ! Grid spacing
    IF ( Nr .GT. 1 ) THEN ! To avoid case with 0 intervals
      deltar = Lr / REAL(Nr-1,dp)
    ELSE
      deltar = 1;
    ENDIF
    ! Discretized r positions ordered as dx*(0 1 2 -3 -2 -1)
    ALLOCATE(rarray(irs:ire))
    IF (NR .GT. 1) THEN
      DO ir = irs,ire
        rarray(ir) = deltar*(MODULO(ir-1,Nr/2)-Nr/2*FLOOR(2.*real(ir-1)/real(Nr)))
      END DO
      rarray(Nr/2+1) = -rarray(Nr/2+1)
    ELSE
      rarray(1) = 0._dp
    ENDIF

  END SUBROUTINE set_rgrid

  SUBROUTINE set_zgrid
    USE prec_const
    IMPLICIT NONE
    INTEGER :: iz
    ! Start and END indices of grid
    izs = 1
    ize = Nr
    ! Grid spacing
    IF ( Nz .GT. 1 ) THEN ! To avoid case with 0 intervals
      deltaz = Lz / REAL(Nz-1,dp)
    ELSE
      deltaz = 1;
    ENDIF
    ! Discretized r positions ordered as dx*(0 1 2 -3 -2 -1)
    ALLOCATE(zarray(irs:ire))
    DO iz = izs,ize
      zarray(iz) = deltaz*(MODULO(iz-1,Nz/2)-Nz/2*FLOOR(2.*real(iz-1)/real(Nz)))
    END DO
    zarray(Nz/2+1) = -zarray(Nz/2+1)

  END SUBROUTINE set_zgrid

  SUBROUTINE set_krgrid
    USE prec_const
    IMPLICIT NONE
    INTEGER :: ikr

    Nkr = Nr/2+1 ! Defined only on positive kr since fields are real
    ! Start and END indices of grid
    ikrs = 1
    ikre = Nkr
    ! Grid spacings
    deltakr = 2._dp*PI/Nr/deltar

    ! Discretized kr positions ordered as dk*(0 1 2)
    ALLOCATE(krarray(ikrs:ikre))
    DO ikr = ikrs,ikre
      krarray(ikr) = REAL(ikr-1,dp) * deltakr
      if (krarray(ikr) .EQ. 0) THEN
        ikr_0 = ikr
      ENDIF
    END DO

    ! Orszag 2/3 filter
    two_third_krmax = 2._dp/3._dp*deltakr*Nkr
    ALLOCATE(AA_r(ikrs:ikre))
    DO ikr = ikrs,ikre
      IF ( (krarray(ikr) .GT. -two_third_krmax) .AND. (krarray(ikr) .LT. two_third_krmax) ) THEN
        AA_r(ikr) = 1._dp;
      ELSE
        AA_r(ikr) = 0._dp;
      ENDIF
    END DO
  END SUBROUTINE set_krgrid

  SUBROUTINE set_kzgrid
    USE prec_const
    USE model, ONLY : kr0KH, ikz0KH
    IMPLICIT NONE

    Nkz = Nz;
    ! Start and END indices of grid
    ikzs = 1
    ikze = Nkz
    ! Grid spacings
    deltakz = 2._dp*PI/Nkz/deltaz

    ! Discretized kz positions ordered as dk*(0 1 2 -3 -2 -1)
    ALLOCATE(kzarray(ikzs:ikze))
    DO ikz = ikzs,ikze
      kzarray(ikz) = deltakz*(MODULO(ikz-1,Nkz/2)-Nkz/2*FLOOR(2.*real(ikz-1)/real(Nkz)))
      if (kzarray(ikz) .EQ. 0) THEN
        ikz_0 = ikz
      ENDIF
    END DO
    kzarray(Nz/2+1) = -kzarray(Nz/2+1)

    ! Orszag 2/3 filter
    two_third_kzmax = 2._dp/3._dp*deltakz*(Nkz/2);
    ALLOCATE(AA_z(ikzs:ikze))
    DO ikz = ikzs,ikze
      IF ( (kzarray(ikz) .GT. -two_third_kzmax) .AND. (kzarray(ikz) .LT. two_third_kzmax) ) THEN
        AA_z(ikz) = 1._dp;
      ELSE
        AA_z(ikz) = 0._dp;
      ENDIF
    END DO

    ! Put kz0RT to the nearest grid point on kz
    ikz0KH = NINT(kr0KH/deltakr)+1
    kr0KH  = kzarray(ikz0KH)
    WRITE(*,*) 'ikz0KH = ', ikz0KH
    WRITE(*,*) 'kr0KH = ', kr0KH

  END SUBROUTINE set_kzgrid

  SUBROUTINE grid_readinputs
    ! Read the input parameters
    USE prec_const
    IMPLICIT NONE
    INTEGER :: lu_in   = 90              ! File duplicated from STDIN

    NAMELIST /GRID/ pmaxe, jmaxe, pmaxi, jmaxi, &
                    Nr,  Lr,  Nz,  Lz, kpar
    READ(lu_in,grid)

  END SUBROUTINE grid_readinputs


  SUBROUTINE grid_outputinputs(fidres, str)
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
    CALL attach(fidres, TRIM(str),   "nr",   nr)
    CALL attach(fidres, TRIM(str),   "Lr",   Lr)
    CALL attach(fidres, TRIM(str),   "nz",   nz)
    CALL attach(fidres, TRIM(str),   "Lz",   Lz)
    CALL attach(fidres, TRIM(str),   "nkr",   nkr)
    CALL attach(fidres, TRIM(str),   "Lkr",   Lkr)
    CALL attach(fidres, TRIM(str),   "nkz",   nkz)
    CALL attach(fidres, TRIM(str),   "Lkz",   Lkz)
  END SUBROUTINE grid_outputinputs

  FUNCTION bare(p_,j_)
    IMPLICIT NONE
    INTEGER :: bare, p_, j_
    bare = (jmaxe+1)*p_ + j_ + 1
  END FUNCTION

  FUNCTION bari(p_,j_)
    IMPLICIT NONE
    INTEGER :: bari, p_, j_
    bari = (jmaxi+1)*p_ + j_ + 1
  END FUNCTION

END MODULE grid
