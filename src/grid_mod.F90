MODULE grid
  ! Grid module for spatial discretization
  USE prec_const
  USE basic

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
  integer(C_INTPTR_T), PUBLIC :: local_nkr, local_nz
  integer(C_INTPTR_T), PUBLIC :: local_nkr_offset, local_nz_offset
  INTEGER,             PUBLIC :: local_nkp
  INTEGER,             PUBLIC :: local_np_e, local_np_i
  integer(C_INTPTR_T), PUBLIC :: local_np_e_offset, local_np_i_offset
  INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: counts_np_e, counts_np_i
  INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: displs_np_e, displs_np_i

  ! Grids containing position in fourier space
  REAL(dp), DIMENSION(:),   ALLOCATABLE, PUBLIC :: krarray
  REAL(dp), DIMENSION(:),   ALLOCATABLE, PUBLIC :: kzarray
  REAL(dp), DIMENSION(:),   ALLOCATABLE, PUBLIC :: kparray     ! kperp array
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: kparray_2D  ! kperp array in 2D
  INTEGER,  DIMENSION(:),   ALLOCATABLE, PUBLIC :: ikparray    ! kperp indices array
  INTEGER,  DIMENSION(:,:), ALLOCATABLE, PUBLIC :: ikparray_2D ! to convert (ikr,ikz) -> ikp
  REAL(dp), PUBLIC, PROTECTED ::  deltakr, deltakz
  INTEGER,  PUBLIC, PROTECTED ::  ikrs, ikre, ikzs, ikze, ikps, ikpe
  INTEGER,  PUBLIC, PROTECTED :: ikr_0, ikz_0 ! Indices of k-grid origin
  INTEGER,  PUBLIC :: ikr, ikz, ip, ij, ikp ! counters

  ! Grid containing the polynomials degrees
  INTEGER,  DIMENSION(:), ALLOCATABLE, PUBLIC :: parray_e
  INTEGER,  DIMENSION(:), ALLOCATABLE, PUBLIC :: parray_i
  INTEGER,  DIMENSION(:), ALLOCATABLE, PUBLIC :: jarray_e
  INTEGER,  DIMENSION(:), ALLOCATABLE, PUBLIC :: jarray_i
  INTEGER, PUBLIC, PROTECTED ::  ips_e,ipe_e, ijs_e,ije_e ! Start and end indices for pol. deg.
  INTEGER, PUBLIC, PROTECTED ::  ips_i,ipe_i, ijs_i,ije_i
  INTEGER, PUBLIC, PROTECTED ::  ipsg_e,ipeg_e, ijsg_e,ijeg_e ! Ghosts start and end indices
  INTEGER, PUBLIC, PROTECTED ::  ipsg_i,ipeg_i, ijsg_i,ijeg_i

  ! Public Functions
  PUBLIC :: init_1Dgrid_distr
  PUBLIC :: set_pgrid, set_jgrid
  PUBLIC :: set_krgrid, set_kzgrid, set_kpgrid
  PUBLIC :: grid_readinputs, grid_outputinputs
  PUBLIC :: bare, bari

  ! Precomputations
  real(dp), PUBLIC, PROTECTED    :: pmaxe_dp, pmaxi_dp, jmaxe_dp,jmaxi_dp

CONTAINS

  SUBROUTINE init_1Dgrid_distr

    ! write(*,*) Nr
    local_nkr        = (Nr/2+1)/num_procs_kr
    ! write(*,*) local_nkr
    local_nkr_offset = rank_r*local_nkr

    if (rank_r .EQ. num_procs_kr-1) local_nkr = (Nr/2+1)-local_nkr_offset

  END SUBROUTINE init_1Dgrid_distr

  SUBROUTINE set_pgrid
    USE prec_const
    IMPLICIT NONE
    INTEGER :: ip, istart, iend, in

    ! Local data distribution
    CALL decomp1D(pmaxe+1, num_procs_p, rank_p, ips_e, ipe_e)
    CALL decomp1D(pmaxi+1, num_procs_p, rank_p, ips_i, ipe_i)
    local_np_e = ipe_e - ips_e + 1
    local_np_i = ipe_i - ips_i + 1
    ! List of shift and local numbers between the different processes (used in scatterv and gatherv)
    ALLOCATE(counts_np_e (1:num_procs_p))
    ALLOCATE(counts_np_i (1:num_procs_p))
    ALLOCATE(displs_np_e (1:num_procs_p))
    ALLOCATE(displs_np_i (1:num_procs_p))
    DO in = 0,num_procs_p-1
      CALL decomp1D(pmaxe+1, num_procs_p, in, istart, iend)
      counts_np_e(in+1) = iend-istart+1
      displs_np_e(in+1) = istart-1
      CALL decomp1D(pmaxi+1, num_procs_p, in, istart, iend)
      counts_np_i(in+1) = iend-istart+1
      displs_np_i(in+1) = istart-1
    ENDDO
    write(*,*) rank_p, ': counts = ', counts_np_e
    write(*,*) rank_p, ': disps = ',  displs_np_e

    ! local grid computation
    ALLOCATE(parray_e(ips_e:ipe_e))
    ALLOCATE(parray_i(ips_i:ipe_i))
    DO ip = ips_e,ipe_e; parray_e(ip) = (ip-1); END DO
    DO ip = ips_i,ipe_i; parray_i(ip) = (ip-1); END DO

    ! Ghosts boundaries
    ipsg_e = ips_e - 2; ipeg_e = ipe_e + 2;
    ipsg_i = ips_i - 2; ipeg_i = ipe_i + 2;
    ! Precomputations
    pmaxe_dp   = real(pmaxe,dp)
    pmaxi_dp   = real(pmaxi,dp)
  END SUBROUTINE set_pgrid

  SUBROUTINE set_jgrid
    USE prec_const
    IMPLICIT NONE
    INTEGER :: ij

    ijs_e = 1; ije_e = jmaxe + 1
    ijs_i = 1; ije_i = jmaxi + 1
    ALLOCATE(jarray_e(ijs_e:ije_e))
    ALLOCATE(jarray_i(ijs_i:ije_i))
    DO ij = ijs_e,ije_e; jarray_e(ij) = ij-1; END DO
    DO ij = ijs_i,ije_i; jarray_i(ij) = ij-1; END DO

    maxj  = MAX(jmaxi, jmaxe)

    ! Ghosts boundaries
    ijsg_e = ijs_e - 1; ijeg_e = ije_e + 1;
    ijsg_i = ijs_i - 1; ijeg_i = ije_i + 1;
    ! Precomputations
    jmaxe_dp   = real(jmaxe,dp)
    jmaxi_dp   = real(jmaxi,dp)
  END SUBROUTINE set_jgrid


  SUBROUTINE set_krgrid
    USE prec_const
    IMPLICIT NONE
    INTEGER :: i_

    Nkr = Nr/2+1 ! Defined only on positive kr since fields are real
    ! Start and END indices of grid
    ikrs = local_nkr_offset + 1
    ikre = ikrs + local_nkr - 1
    ALLOCATE(krarray(ikrs:ikre))

    ! Grid spacings
    IF (Lr .EQ. 0) THEN
      deltakr = 1._dp
    ELSE
      deltakr = 2._dp*PI/Lr
    ENDIF

    ! Discretized kr positions ordered as dk*(0 1 2 3)
    DO ikr = ikrs,ikre
      krarray(ikr) = REAL(ikr-1,dp) * deltakr
      IF (krarray(ikr) .EQ. 0) ikr_0 = ikr
    END DO

    ! Orszag 2/3 filter
    two_third_krmax = 2._dp/3._dp*deltakr*Nkr
    ALLOCATE(AA_r(ikrs:ikre))
    DO ikr = ikrs,ikre
      IF ( (krarray(ikr) .LT. two_third_krmax) ) THEN
        AA_r(ikr) = 1._dp;
      ELSE
        AA_r(ikr) = 0._dp;
      ENDIF
    END DO
  END SUBROUTINE set_krgrid

  SUBROUTINE set_kzgrid
    USE prec_const
    IMPLICIT NONE
    INTEGER :: i_, counter

    Nkz = Nz;
    ! Start and END indices of grid
    ikzs = 1
    ikze = Nkz
    ALLOCATE(kzarray(ikzs:ikze))

    ! Grid spacings and discretized kz positions ordered as dk*(0 1 2 3 -2 -1)
    IF (Lz .EQ. 0) THEN ! 1D linear case
      deltakz    = 1._dp
      kzarray(1) = 0
      ikz_0      = 1
    ELSE
      deltakz = 2._dp*PI/Lz
      DO ikz = ikzs,ikze
        kzarray(ikz) = deltakz*(MODULO(ikz-1,Nkz/2)-Nkz/2*FLOOR(2.*real(ikz-1)/real(Nkz)))
        if (ikz .EQ. Nz/2+1)     kzarray(ikz) = -kzarray(ikz)
        IF (kzarray(ikz) .EQ. 0) ikz_0 = ikz
      END DO
    ENDIF

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
  END SUBROUTINE set_kzgrid

  SUBROUTINE set_kpgrid !Precompute the grid of kperp
    IMPLICIT NONE
    INTEGER :: ikz_sym, tmp_, counter
    ! 2D to 1D indices array convertor
    ALLOCATE(ikparray_2D(ikrs:ikre,ikzs:ikze))
    ALLOCATE( kparray_2D(ikrs:ikre,ikzs:ikze))
    ! local number of different kperp
    local_nkp = local_nkr * (local_nkr-1)/2 + 1
    ! Allocate 1D array of kperp values and indices
    ALLOCATE(ikparray(1:local_nkr))
    ALLOCATE( kparray(1:local_nkr))

    ! Fill the arrays
    tmp_ = 0; counter = 1
    DO ikz = ikzs,ikze
      DO ikr = ikrs,ikre
        ! Set a symmetry on kz
        IF (ikz .LE. Nkz/2+1) THEN
          ikz_sym = ikz
        ELSE
          ikz_sym = Nkz+2 - ikz
        ENDIF
        ! Formula to find the 2D to 1D kperp equivalences ordered as
        !      10
        !    6 9
        !  3 5 8
        !1 2 4 7  etc...
        ikp = MAX(ikr-1,ikz_sym-1)*MIN(ikr-1,ikz_sym-1)/2 + min(ikr-1,ikz_sym-1)
        ikparray_2D(ikr,ikz) = ikp
        kparray_2D(ikr,ikz)  = SQRT(krarray(ikr)**2 + kzarray(ikz)**2)
      ENDDO
    ENDDO

  END SUBROUTINE

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

  SUBROUTINE decomp1D( n, numprocs, myid, s, e )
      INTEGER :: n, numprocs, myid, s, e
      INTEGER :: nlocal
      INTEGER :: deficit

      nlocal   = n / numprocs
      s        = myid * nlocal + 1
      deficit  = MOD(n,numprocs)
      s        = s + MIN(myid,deficit)
      IF (myid .LT. deficit) nlocal = nlocal + 1
      e = s + nlocal - 1
      IF (e .GT. n .OR. myid .EQ. numprocs-1) e = n
  END SUBROUTINE decomp1D

END MODULE grid
