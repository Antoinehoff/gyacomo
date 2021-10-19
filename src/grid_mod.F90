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
  INTEGER,  PUBLIC, PROTECTED :: maxj  = 1      ! The maximal Laguerre-moment
  INTEGER,  PUBLIC, PROTECTED :: dmaxe = 1      ! The maximal full GF set of e-moments v^dmax
  INTEGER,  PUBLIC, PROTECTED :: dmaxi = 1      ! The maximal full GF set of i-moments v^dmax
  INTEGER,  PUBLIC, PROTECTED :: Nx    = 16     ! Number of total internal grid points in x
  REAL(dp), PUBLIC, PROTECTED :: Lx    = 1._dp  ! horizontal length of the spatial box
  INTEGER,  PUBLIC, PROTECTED :: Ny    = 16     ! Number of total internal grid points in y
  REAL(dp), PUBLIC, PROTECTED :: Ly    = 1._dp  ! vertical length of the spatial box
  INTEGER,  PUBLIC, PROTECTED :: Nz    = 1      ! Number of total perpendicular planes
  REAL(dp), PUBLIC, PROTECTED :: q0    = 1._dp  ! safety factor
  REAL(dp), PUBLIC, PROTECTED :: shear = 0._dp  ! magnetic field shear
  REAL(dp), PUBLIC, PROTECTED :: eps   = 0._dp ! inverse aspect ratio
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
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: zarray, zarray_full
  REAL(dp), PUBLIC, PROTECTED ::  deltax,  deltay, deltaz
  INTEGER,  PUBLIC, PROTECTED  ::  ixs,  ixe,  iys,  iye, izs, ize
  INTEGER,  PUBLIC :: ir,iz ! counters
  integer(C_INTPTR_T), PUBLIC :: local_nkx, local_nky
  integer(C_INTPTR_T), PUBLIC :: local_nkx_offset, local_nky_offset
  INTEGER,             PUBLIC :: local_nkp
  INTEGER,             PUBLIC :: local_np_e, local_np_i
  INTEGER,             PUBLIC :: total_np_e, total_np_i
  integer(C_INTPTR_T), PUBLIC :: local_np_e_offset, local_np_i_offset
  INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: counts_np_e, counts_np_i
  INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: displs_np_e, displs_np_i

  ! Grids containing position in fourier space
  REAL(dp), DIMENSION(:),   ALLOCATABLE, PUBLIC :: kxarray, kxarray_full
  REAL(dp), DIMENSION(:),   ALLOCATABLE, PUBLIC :: kyarray, kyarray_full
  REAL(dp), DIMENSION(:),   ALLOCATABLE, PUBLIC :: kparray     ! kperp array
  REAL(dp), PUBLIC, PROTECTED ::  deltakx, deltaky, kx_max, ky_max, kp_max
  INTEGER,  PUBLIC, PROTECTED ::  ikxs, ikxe, ikys, ikye, ikps, ikpe
  INTEGER,  PUBLIC, PROTECTED :: ikx_0, iky_0, ikx_max, iky_max ! Indices of k-grid origin and max
  INTEGER,  PUBLIC :: ikx, iky, ip, ij, ikp ! counters
  LOGICAL,  PUBLIC, PROTECTED :: contains_kx0   = .false. ! rank of the proc containing kx=0 indices
  LOGICAL,  PUBLIC, PROTECTED :: contains_kxmax = .false. ! rank of the proc containing kx=max indices

  ! Grid containing the polynomials degrees
  INTEGER,  DIMENSION(:), ALLOCATABLE, PUBLIC :: parray_e, parray_e_full
  INTEGER,  DIMENSION(:), ALLOCATABLE, PUBLIC :: parray_i, parray_i_full
  INTEGER,  DIMENSION(:), ALLOCATABLE, PUBLIC :: jarray_e, jarray_e_full
  INTEGER,  DIMENSION(:), ALLOCATABLE, PUBLIC :: jarray_i, jarray_i_full
  INTEGER, PUBLIC, PROTECTED ::  ips_e,ipe_e, ijs_e,ije_e ! Start and end indices for pol. deg.
  INTEGER, PUBLIC, PROTECTED ::  ips_i,ipe_i, ijs_i,ije_i
  INTEGER, PUBLIC, PROTECTED ::  ipsg_e,ipeg_e, ijsg_e,ijeg_e ! Ghosts start and end indices
  INTEGER, PUBLIC, PROTECTED ::  ipsg_i,ipeg_i, ijsg_i,ijeg_i
  INTEGER, PUBLIC, PROTECTED ::  deltape, ip0_e, ip1_e, ip2_e ! Pgrid spacing and moment 0,1,2 index
  INTEGER, PUBLIC, PROTECTED ::  deltapi, ip0_i, ip1_i, ip2_i

  ! Public Functions
  PUBLIC :: init_1Dgrid_distr
  PUBLIC :: set_pgrid, set_jgrid
  PUBLIC :: set_kxgrid, set_kygrid, set_kpgrid, set_zgrid
  PUBLIC :: grid_readinputs, grid_outputinputs
  PUBLIC :: bare, bari

  ! Precomputations
  real(dp), PUBLIC, PROTECTED    :: pmaxe_dp, pmaxi_dp, jmaxe_dp,jmaxi_dp

CONTAINS


  SUBROUTINE grid_readinputs
    ! Read the input parameters
    USE prec_const
    IMPLICIT NONE
    INTEGER :: lu_in   = 90              ! File duplicated from STDIN

    NAMELIST /GRID/ pmaxe, jmaxe, pmaxi, jmaxi, &
                    Nx,  Lx,  Ny,  Ly, Nz, q0, shear, eps
    READ(lu_in,grid)

    !! Compute the maximal degree of full GF moments set
    !   i.e. : all moments N_a^pj s.t. p+2j<=d are simulated (see GF closure)
    dmaxe = min(pmaxe,2*jmaxe+1)
    dmaxi = min(pmaxi,2*jmaxi+1)

    ! If no parallel dim (Nz=1), the moment hierarchy is separable between odds and even P
    !! and since the energy is injected in P=0 and P=2 for density/temperature gradients
    !! there is no need of simulating the odd p which will only be damped.
    !! We define in this case a grid Parray = 0,2,4,...,Pmax i.e. deltap = 2 instead of 1
    !! to spare computation
    IF(Nz .EQ. 1) THEN
      deltape = 2; deltapi = 2;
    ELSE
      deltape = 1; deltapi = 1;
    ENDIF

  END SUBROUTINE grid_readinputs

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

    ! Total number of Hermite polynomials we will evolve
    total_np_e = (Pmaxe/deltape) + 1
    total_np_i = (Pmaxi/deltapi) + 1
    ! Build the full grids on process 0 to diagnose it without comm
    ALLOCATE(parray_e_full(1:total_np_e))
    ALLOCATE(parray_i_full(1:total_np_i))
    ! P
    DO ip = 1,total_np_e; parray_e_full(ip) = (ip-1)*deltape; END DO
    DO ip = 1,total_np_i; parray_i_full(ip) = (ip-1)*deltapi; END DO
    !! Parallel data distribution
    ! Local data distribution
    CALL decomp1D(total_np_e, num_procs_p, rank_p, ips_e, ipe_e)
    CALL decomp1D(total_np_i, num_procs_p, rank_p, ips_i, ipe_i)
    local_np_e = ipe_e - ips_e + 1
    local_np_i = ipe_i - ips_i + 1
    ! List of shift and local numbers between the different processes (used in scatterv and gatherv)
    ALLOCATE(counts_np_e (1:num_procs_p))
    ALLOCATE(counts_np_i (1:num_procs_p))
    ALLOCATE(displs_np_e (1:num_procs_p))
    ALLOCATE(displs_np_i (1:num_procs_p))
    DO in = 0,num_procs_p-1
      CALL decomp1D(total_np_e, num_procs_p, in, istart, iend)
      counts_np_e(in+1) = iend-istart+1
      displs_np_e(in+1) = istart-1
      CALL decomp1D(total_np_i, num_procs_p, in, istart, iend)
      counts_np_i(in+1) = iend-istart+1
      displs_np_i(in+1) = istart-1
    ENDDO

    ! local grid computation
    ALLOCATE(parray_e(ips_e:ipe_e))
    ALLOCATE(parray_i(ips_i:ipe_i))
    DO ip = ips_e,ipe_e
      parray_e(ip) = (ip-1)*deltape
      ! Storing indices of particular degrees for DG and fluid moments computations
      IF(parray_e(ip) .EQ. 0) ip0_e = ip
      IF(parray_e(ip) .EQ. 1) ip1_e = ip
      IF(parray_e(ip) .EQ. 2) ip2_e = ip
    END DO
    DO ip = ips_i,ipe_i
      parray_i(ip) = (ip-1)*deltapi
      IF(parray_i(ip) .EQ. 0) ip0_i = ip
      IF(parray_i(ip) .EQ. 1) ip1_i = ip
      IF(parray_i(ip) .EQ. 2) ip2_i = ip
    END DO
    !DGGK operator uses moments at index p=2 (ip=3) for the p=0 term so the
    ! process that contains ip=1 MUST contain ip=3 as well for both e and i.
    IF(((ips_e .EQ. ip0_e) .OR. (ips_i .EQ. ip0_e)) .AND. ((ipe_e .LT. ip2_e) .OR. (ipe_i .LT. ip2_i)))&
     WRITE(*,*) "Warning : distribution along p may not work with DGGK"

    ! Ghosts boundaries
    ipsg_e = ips_e - 2/deltape; ipeg_e = ipe_e + 2/deltape;
    ipsg_i = ips_i - 2/deltapi; ipeg_i = ipe_i + 2/deltapi;
    ! Precomputations
    pmaxe_dp   = real(pmaxe,dp)
    pmaxi_dp   = real(pmaxi,dp)
  END SUBROUTINE set_pgrid

  SUBROUTINE set_jgrid
    USE prec_const
    IMPLICIT NONE
    INTEGER :: ij

    ! Build the full grids on process 0 to diagnose it without comm
    ALLOCATE(jarray_e_full(1:jmaxe+1))
    ALLOCATE(jarray_i_full(1:jmaxi+1))
    ! J
    DO ij = 1,jmaxe+1; jarray_e_full(ij) = (ij-1); END DO
    DO ij = 1,jmaxi+1; jarray_i_full(ij) = (ij-1); END DO
    ! Local data
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
    USE model, ONLY: NON_LIN
    IMPLICIT NONE
    INTEGER :: i_

    Nkx = Nx/2+1 ! Defined only on positive kx since fields are real
    ! Grid spacings
    IF (Lx .EQ. 0) THEN
      deltakx = 1._dp
      kx_max  = 0._dp
    ELSE
      deltakx = 2._dp*PI/Lx
      kx_max  = Nkx*deltakx
    ENDIF

    ! Build the full grids on process 0 to diagnose it without comm
    ALLOCATE(kxarray_full(1:Nkx))
    DO ikx = 1,Nkx
     kxarray_full(ikx) = REAL(ikx-1,dp) * deltakx
    END DO

    !! Parallel distribution
    ! Start and END indices of grid
    ikxs = local_nkx_offset + 1
    ikxe = ikxs + local_nkx - 1
    ALLOCATE(kxarray(ikxs:ikxe))

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
    two_third_kxmax = 2._dp/3._dp*deltakx*(Nkx-1)
    ALLOCATE(AA_x(ikxs:ikxe))
    DO ikx = ikxs,ikxe
      IF ( (kxarray(ikx) .LT. two_third_kxmax) .OR. (.NOT. NON_LIN)) THEN
        AA_x(ikx) = 1._dp;
      ELSE
        AA_x(ikx) = 0._dp;
      ENDIF
    END DO
  END SUBROUTINE set_kxgrid

  SUBROUTINE set_kygrid
    USE prec_const
    USE model, ONLY: NON_LIN
    IMPLICIT NONE
    INTEGER :: i_, counter

    Nky = Ny;
    ALLOCATE(kyarray_full(1:Nky))
    ! Local data
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
    ! Build the full grids on process 0 to diagnose it without comm
    ! ky
    IF (Nky .GT. 1) THEN
     DO iky = 1,Nky
       kyarray_full(iky) = deltaky*(MODULO(iky-1,Nky/2)-Nky/2*FLOOR(2.*real(iky-1)/real(Nky)))
       IF (iky .EQ. Ny/2+1) kyarray_full(iky) = -kyarray_full(iky)
     END DO
    ELSE
     kyarray_full(1) =  0
    ENDIF
    ! Orszag 2/3 filter
    two_third_kymax = 2._dp/3._dp*deltaky*(Nky/2-1);
    ALLOCATE(AA_y(ikys:ikye))
    DO iky = ikys,ikye
      IF ( (kyarray(iky) .GT. -two_third_kymax) .AND. &
           (kyarray(iky) .LT. two_third_kymax) .OR. (.NOT. NON_LIN)) THEN
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
    two_third_kpmax = SQRT(two_third_kxmax**2+two_third_kymax**2)
    kp_max = 3._dp/2._dp * two_third_kpmax
  END SUBROUTINE

  SUBROUTINE set_zgrid
    USE prec_const
    IMPLICIT NONE
    INTEGER :: i_
    ! Build the full grids on process 0 to diagnose it without comm
    ALLOCATE(zarray_full(1:Nz))
    ! z
    IF (Nz .GT. 1) THEN
      DO iz = 1,Nz
        zarray_full(iz) = deltaz*(iz-1)
      END DO
    ELSE
      zarray_full(1) =  0
    ENDIF
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
    CALL attach(fidres, TRIM(str),"shear",shear)
    CALL attach(fidres, TRIM(str),  "eps",  eps)
    CALL attach(fidres, TRIM(str),  "Nkx",  Nkx)
    CALL attach(fidres, TRIM(str),  "Lkx",  Lkx)
    CALL attach(fidres, TRIM(str),  "Nky",  Nky)
    CALL attach(fidres, TRIM(str),  "Lky",  Lky)
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
