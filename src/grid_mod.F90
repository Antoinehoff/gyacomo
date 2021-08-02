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
  INTEGER,  PUBLIC, PROTECTED :: Nx    = 16     ! Number of total internal grid points in x
  REAL(dp), PUBLIC, PROTECTED :: Lx    = 1._dp  ! horizontal length of the spatial box
  INTEGER,  PUBLIC, PROTECTED :: Ny    = 16     ! Number of total internal grid points in y
  REAL(dp), PUBLIC, PROTECTED :: Ly    = 1._dp  ! vertical length of the spatial box
  INTEGER,  PUBLIC, PROTECTED :: Nz    = 1      ! Number of total perpendicular planes
  REAL(dp), PUBLIC, PROTECTED :: q0    = 1._dp  ! q factor
  INTEGER,  PUBLIC, PROTECTED :: Nkx   = 8      ! Number of total internal grid points in kx
  REAL(dp), PUBLIC, PROTECTED :: Lkx   = 1._dp  ! horizontal length of the fourier box
  INTEGER,  PUBLIC, PROTECTED :: Nky   = 16     ! Number of total internal grid points in ky
  REAL(dp), PUBLIC, PROTECTED :: Lky   = 1._dp  ! vertical length of the fourier box
  REAL(dp), PUBLIC, PROTECTED :: kpar  = 0_dp   ! parallel wave vector component
  ! For Orszag filter
  REAL(dp), PUBLIC, PROTECTED :: two_third_kxmax
  REAL(dp), PUBLIC, PROTECTED :: two_third_kymax
  REAL(dp), PUBLIC, PROTECTED :: two_third_kpmax

  ! 1D Antialiasing arrays (2/3 rule)
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: AA_x
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: AA_y

  ! Grids containing position in physical space
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: xarray
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: yarray
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: zarray
  REAL(dp), PUBLIC, PROTECTED ::  deltax,  deltay, deltaz
  INTEGER,  PUBLIC, PROTECTED  ::  ixs,  ixe,  iys,  iye, izs, ize
  INTEGER,  PUBLIC :: ir,iz ! counters
  integer(C_INTPTR_T), PUBLIC :: local_nkx, local_nky
  integer(C_INTPTR_T), PUBLIC :: local_nkx_offset, local_nky_offset
  INTEGER,             PUBLIC :: local_nkp
  INTEGER,             PUBLIC :: local_np_e, local_np_i
  integer(C_INTPTR_T), PUBLIC :: local_np_e_offset, local_np_i_offset
  INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: counts_np_e, counts_np_i
  INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: displs_np_e, displs_np_i

  ! Grids containing position in fourier space
  REAL(dp), DIMENSION(:),   ALLOCATABLE, PUBLIC :: kxarray
  REAL(dp), DIMENSION(:),   ALLOCATABLE, PUBLIC :: kyarray
  REAL(dp), DIMENSION(:),   ALLOCATABLE, PUBLIC :: kparray     ! kperp array
  REAL(dp), PUBLIC, PROTECTED ::  deltakx, deltaky, kx_max, ky_max, kp_max
  INTEGER,  PUBLIC, PROTECTED ::  ikxs, ikxe, ikys, ikye, ikps, ikpe
  INTEGER,  PUBLIC, PROTECTED :: ikx_0, iky_0, ikx_max, iky_max ! Indices of k-grid origin and max
  INTEGER,  PUBLIC :: ikx, iky, ip, ij, ikp ! counters
  LOGICAL,  PUBLIC, PROTECTED :: contains_kx0   = .false. ! rank of the proc containing kx=0 indices
  LOGICAL,  PUBLIC, PROTECTED :: contains_kxmax = .false. ! rank of the proc containing kx=max indices

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
  PUBLIC :: set_kxgrid, set_kygrid, set_kpgrid, set_zgrid
  PUBLIC :: grid_readinputs, grid_outputinputs
  PUBLIC :: bare, bari

  ! Precomputations
  real(dp), PUBLIC, PROTECTED    :: pmaxe_dp, pmaxi_dp, jmaxe_dp,jmaxi_dp

CONTAINS

  SUBROUTINE init_1Dgrid_distr

    ! write(*,*) Nx
    local_nkx        = (Nx/2+1)/num_procs_kx
    ! write(*,*) local_nkx
    local_nkx_offset = rank_kx*local_nkx

    if (rank_kx .EQ. num_procs_kx-1) local_nkx = (Nx/2+1)-local_nkx_offset

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
    !DGGK operator uses moments at index p=2 (ip=3) for the p=0 term so the
    ! process that contains ip=1 MUST contain ip=3 as well for both e and i.
    IF(((ips_e .EQ. 1) .OR. (ips_i .EQ. 1)) .AND. ((ipe_e .LT. 3) .OR. (ipe_i .LT. 3)))&
     WRITE(*,*) "Warning : distribution along p may not work with DGGK"
    ENDDO

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


  SUBROUTINE set_kxgrid
    USE prec_const
    IMPLICIT NONE
    INTEGER :: i_

    Nkx = Nx/2+1 ! Defined only on positive kx since fields are real
    ! Start and END indices of grid
    ikxs = local_nkx_offset + 1
    ikxe = ikxs + local_nkx - 1
    ALLOCATE(kxarray(ikxs:ikxe))

    ! Grid spacings
    IF (Lx .EQ. 0) THEN
      deltakx = 1._dp
      kx_max  = 0._dp
    ELSE
      deltakx = 2._dp*PI/Lx
      kx_max  = (Nx/2+1)*deltakx
    ENDIF

    ! Creating a grid ordered as dk*(0 1 2 3)
    DO ikx = ikxs,ikxe
      kxarray(ikx) = REAL(ikx-1,dp) * deltakx
      ! Finding kx=0
      IF (kxarray(ikx) .EQ. 0) THEN
        ikx_0 = ikx
        contains_kx0 = .true.
      ENDIF
      ! Finding kxmax idx
      IF (kxarray(ikx) .EQ. kx_max) THEN
        ikx_max = ikx
        contains_kxmax = .true.
      ENDIF
    END DO

    ! Orszag 2/3 filter
    two_third_kxmax = 2._dp/3._dp*deltakx*Nkx
    ALLOCATE(AA_x(ikxs:ikxe))
    DO ikx = ikxs,ikxe
      IF ( (kxarray(ikx) .LT. two_third_kxmax) ) THEN
        AA_x(ikx) = 1._dp;
      ELSE
        AA_x(ikx) = 0._dp;
      ENDIF
    END DO
  END SUBROUTINE set_kxgrid

  SUBROUTINE set_kygrid
    USE prec_const
    IMPLICIT NONE
    INTEGER :: i_, counter

    Nky = Ny;
    ! Start and END indices of grid
    ikys = 1
    ikye = Nky
    ALLOCATE(kyarray(ikys:ikye))

    IF (Ly .EQ. 0) THEN ! 1D linear case
      deltaky    = 1._dp
      kyarray(1) = 0
      iky_0      = 1
      iky_max    = 1
    ELSE
      deltaky = 2._dp*PI/Ly
      ky_max  = (Ny/2)*deltakx
      ! Creating a grid ordered as dk*(0 1 2 3 -2 -1)
      DO iky = ikys,ikye
        kyarray(iky) = deltaky*(MODULO(iky-1,Nky/2)-Nky/2*FLOOR(2.*real(iky-1)/real(Nky)))
        if (iky .EQ. Ny/2+1)     kyarray(iky) = -kyarray(iky)
        ! Finding ky=0
        IF (kyarray(iky) .EQ. 0) iky_0 = iky
        ! Finding kymax
        IF (kyarray(ikx) .EQ. ky_max) ikx_max = ikx
      END DO
    ENDIF

    ! Orszag 2/3 filter
    two_third_kymax = 2._dp/3._dp*deltaky*(Nky/2);
    ALLOCATE(AA_y(ikys:ikye))
    DO iky = ikys,ikye
      IF ( (kyarray(iky) .GT. -two_third_kymax) .AND. (kyarray(iky) .LT. two_third_kymax) ) THEN
        AA_y(iky) = 1._dp;
      ELSE
        AA_y(iky) = 0._dp;
      ENDIF
    END DO
  END SUBROUTINE set_kygrid

  SUBROUTINE set_kpgrid !Precompute the grid of kperp
    IMPLICIT NONE
    INTEGER :: iky_sym, tmp_, counter
    REAL(dp):: local_kp_min, local_kp_max
    ! Find the min and max kperp to load subsequent GK matrices
    local_kp_min = kxarray(ikxs) !smallest local kperp is on the kx axis
    local_kp_max = SQRT(kxarray(ikxe)**2 + kyarray(Nky/2+1)**2)
    ikps = ikxs
    ikpe = INT(CEILING(local_kp_max/deltakx))+2
    ! local number of different kperp
    local_nkp = ikpe - ikps + 1
    ! Allocate 1D array of kperp values and indices
    ALLOCATE(kparray(ikps:ikpe))
    DO ikp = ikps,ikpe
      kparray(ikp) = REAL(ikp-1,dp) * deltakx
    ENDDO
    write(*,*) rank_kx, ': ikps = ', ikps, 'ikpe = ',ikpe
    two_third_kpmax = SQRT(two_third_kxmax**2+two_third_kymax**2)
    kp_max = 3._dp/2._dp * two_third_kpmax
  END SUBROUTINE

  SUBROUTINE set_zgrid
    USE prec_const
    IMPLICIT NONE
    INTEGER :: i_
    ! Start and END indices of grid
    izs = 1
    ize = Nz
    ALLOCATE(zarray(izs:ize))
    IF (Nz .EQ. 1) THEN ! full perp case
      deltaz    = 1._dp
      zarray(1) = 0
    ELSE
      deltaz = q0*2._dp*PI/REAL(Nz+1,dp)
      DO iz = izs,ize
        zarray(iz) = REAL((iz-1),dp)*deltaz
      ENDDO
    ENDIF
    if(my_id.EQ.0) write(*,*) '#parallel planes = ', Nz
  END SUBROUTINE set_zgrid

  SUBROUTINE grid_readinputs
    ! Read the input parameters
    USE prec_const
    IMPLICIT NONE
    INTEGER :: lu_in   = 90              ! File duplicated from STDIN

    NAMELIST /GRID/ pmaxe, jmaxe, pmaxi, jmaxi, &
                    Nx,  Lx,  Ny,  Ly, Nz, q0
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
    CALL attach(fidres, TRIM(str),   "Nx",   Nx)
    CALL attach(fidres, TRIM(str),   "Lx",   Lx)
    CALL attach(fidres, TRIM(str),   "Ny",   Ny)
    CALL attach(fidres, TRIM(str),   "Ly",   Ly)
    CALL attach(fidres, TRIM(str),   "Nz",   Nz)
    CALL attach(fidres, TRIM(str),   "q0",   q0)
    CALL attach(fidres, TRIM(str),   "Nkx",   Nkx)
    CALL attach(fidres, TRIM(str),   "Lkx",   Lkx)
    CALL attach(fidres, TRIM(str),   "Nky",   Nky)
    CALL attach(fidres, TRIM(str),   "Lky",   Lky)
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
