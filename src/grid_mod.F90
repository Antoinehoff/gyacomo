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
  INTEGER,  PUBLIC, PROTECTED :: Nexc  = 1      ! factor to increase Lx when shear>0 (Lx = Nexc/kymin/shear)
  INTEGER,  PUBLIC, PROTECTED :: Ny    = 16     ! Number of total internal grid points in y
  REAL(dp), PUBLIC, PROTECTED :: Ly    = 1._dp  ! vertical length of the spatial box
  INTEGER,  PUBLIC, PROTECTED :: Nz    = 1      ! Number of total perpendicular planes
  INTEGER,  PUBLIC, PROTECTED :: Npol  = 1      ! number of poloidal turns
  INTEGER,  PUBLIC, PROTECTED :: Odz   = 4      ! order of z interp and derivative schemes
  INTEGER,  PUBLIC, PROTECTED :: Nkx   = 8      ! Number of total internal grid points in kx
  REAL(dp), PUBLIC, PROTECTED :: Lkx   = 1._dp  ! horizontal length of the fourier box
  INTEGER,  PUBLIC, PROTECTED :: Nky   = 16     ! Number of total internal grid points in ky
  REAL(dp), PUBLIC, PROTECTED :: Lky   = 1._dp  ! vertical length of the fourier box
  REAL(dp), PUBLIC, PROTECTED :: kpar  = 0_dp   ! parallel wave vector component
  ! For Orszag filter
  REAL(dp), PUBLIC, PROTECTED :: two_third_kxmax
  REAL(dp), PUBLIC, PROTECTED :: two_third_kymax
  REAL(dp), PUBLIC :: two_third_kpmax

  ! 1D Antialiasing arrays (2/3 rule)
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: AA_x
  REAL(dp), DIMENSION(:), ALLOCATABLE, PUBLIC :: AA_y

  ! Grids containing position in physical space
  REAL(dp), DIMENSION(:),   ALLOCATABLE, PUBLIC :: xarray
  REAL(dp), DIMENSION(:),   ALLOCATABLE, PUBLIC :: yarray
  ! Local and global z grids, 2D since it has to store odd and even grids
  REAL(dp), DIMENSION(:,:), ALLOCATABLE, PUBLIC :: zarray
  REAL(dp), DIMENSION(:),   ALLOCATABLE, PUBLIC :: zarray_full
  ! local z weights for computing simpson rule
  INTEGER,  DIMENSION(:),   ALLOCATABLE, PUBLIC :: zweights_SR
  REAL(dp), PUBLIC, PROTECTED  ::  deltax,  deltay, deltaz, inv_deltaz
  REAL(dp), PUBLIC, PROTECTED  ::  diff_kx_coeff, diff_ky_coeff, diff_dz_coeff
  INTEGER,  PUBLIC, PROTECTED  ::  ixs,  ixe,  iys,  iye,  izs,  ize
  INTEGER,  PUBLIC, PROTECTED  ::  izgs, izge ! ghosts
  LOGICAL,  PUBLIC, PROTECTED  ::  SG = .true.! shifted grid flag
  INTEGER,  PUBLIC :: ir,iz ! counters
  ! Data about parallel distribution for ky.kx
  integer(C_INTPTR_T), PUBLIC :: local_nkx, local_nky
  integer(C_INTPTR_T), PUBLIC :: local_nkx_offset, local_nky_offset
  INTEGER,             PUBLIC :: local_nkp
  INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: counts_nkx, counts_nky
  INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: displs_nkx, displs_nky
  ! "" for p
  INTEGER,             PUBLIC :: local_np_e, local_np_i
  INTEGER,             PUBLIC :: total_np_e, total_np_i, Np_e, Np_i
  integer(C_INTPTR_T), PUBLIC :: local_np_e_offset, local_np_i_offset
  INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: rcv_p_e, rcv_p_i
  INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: dsp_p_e, dsp_p_i
  ! "" for z
  INTEGER,             PUBLIC :: local_nz
  INTEGER,             PUBLIC :: total_nz
  integer(C_INTPTR_T), PUBLIC :: local_nz_offset
  INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: counts_nz
  INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC :: displs_nz
  ! "" for j (not parallelized)
  INTEGER,             PUBLIC :: local_nj_e, local_nj_i, Nj_e, Nj_i
  ! Grids containing position in fourier space
  REAL(dp), DIMENSION(:),     ALLOCATABLE, PUBLIC :: kxarray, kxarray_full
  REAL(dp), DIMENSION(:),     ALLOCATABLE, PUBLIC :: kyarray, kyarray_full
  ! Kperp array depends on kx, ky, z (geometry), eo (even or odd zgrid)
  REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, PUBLIC :: kparray
  REAL(dp), PUBLIC, PROTECTED ::  deltakx, deltaky, kx_max, ky_max, kx_min, ky_min!, kp_max
  REAL(dp), PUBLIC, PROTECTED ::  local_kxmax, local_kymax
  INTEGER,  PUBLIC, PROTECTED ::  ikxs, ikxe, ikys, ikye!, ikps, ikpe
  INTEGER,  PUBLIC, PROTECTED ::  ikx_0, iky_0, ikx_max, iky_max ! Indices of k-grid origin and max
  INTEGER,  PUBLIC            ::  ikx, iky, ip, ij, ikp, pp2, eo ! counters
  LOGICAL,  PUBLIC, PROTECTED ::  contains_kx0   = .false. ! flag if the proc contains kx=0 index
  LOGICAL,  PUBLIC, PROTECTED ::  contains_ky0   = .false. ! flag if the proc contains ky=0 index
  LOGICAL,  PUBLIC, PROTECTED ::  contains_kymax = .false. ! flag if the proc contains kx=kxmax index
  LOGICAL,  PUBLIC, PROTECTED ::  contains_zmax  = .false. ! flag if the proc contains z=pi-dz index
  LOGICAL,  PUBLIC, PROTECTED ::  contains_zmin  = .false. ! flag if the proc contains z=-pi index
  LOGICAL,  PUBLIC, PROTECTED ::       SINGLE_KY = .false. ! to check if it is a single non 0 ky simulation
  ! Grid containing the polynomials degrees
  INTEGER,  DIMENSION(:), ALLOCATABLE, PUBLIC :: parray_e, parray_e_full
  INTEGER,  DIMENSION(:), ALLOCATABLE, PUBLIC :: parray_i, parray_i_full
  INTEGER,  DIMENSION(:), ALLOCATABLE, PUBLIC :: jarray_e, jarray_e_full
  INTEGER,  DIMENSION(:), ALLOCATABLE, PUBLIC :: jarray_i, jarray_i_full
  INTEGER,  PUBLIC, PROTECTED ::  ips_e,ipe_e, ijs_e,ije_e ! Start and end indices for pol. deg.
  INTEGER,  PUBLIC, PROTECTED ::  ips_i,ipe_i, ijs_i,ije_i
  INTEGER,  PUBLIC, PROTECTED ::  ipgs_e,ipge_e, ijgs_e,ijge_e ! Ghosts start and end indices
  INTEGER,  PUBLIC, PROTECTED ::  ipgs_i,ipge_i, ijgs_i,ijge_i
  INTEGER,  PUBLIC, PROTECTED ::  deltape, ip0_e, ip1_e, ip2_e, ip3_e ! Pgrid spacing and moment 0,1,2 index
  INTEGER,  PUBLIC, PROTECTED ::  deltapi, ip0_i, ip1_i, ip2_i, ip3_i
  LOGICAL,  PUBLIC, PROTECTED ::  CONTAINS_ip0_e, CONTAINS_ip0_i
  LOGICAL,  PUBLIC, PROTECTED ::  CONTAINS_ip1_e, CONTAINS_ip1_i
  LOGICAL,  PUBLIC, PROTECTED ::  CONTAINS_ip2_e, CONTAINS_ip2_i
  LOGICAL,  PUBLIC, PROTECTED ::  CONTAINS_ip3_e, CONTAINS_ip3_i
  LOGICAL,  PUBLIC, PROTECTED ::  SOLVE_POISSON, SOLVE_AMPERE
  INTEGER,  PUBLIC, PROTECTED ::  ij0_i, ij0_e
! Usefull inverse numbers
  REAL(dp), PUBLIC, PROTECTED :: inv_Nx, inv_Ny, inv_Nz

  ! Public Functions
  PUBLIC :: init_1Dgrid_distr
  PUBLIC :: set_pgrid, set_jgrid
  PUBLIC :: set_kxgrid, set_kygrid, set_zgrid
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
                    Nx, Lx, Nexc, Ny, Ly, Nz, Npol, SG
    READ(lu_in,grid)

    IF(Nz .EQ. 1) & ! overwrite SG option if Nz = 1 for safety of use
      SG      = .FALSE.

    !! Compute the maximal degree of full GF moments set
    !   i.e. : all moments N_a^pj s.t. p+2j<=d are simulated (see GF closure)
    dmaxe = min(pmaxe,2*jmaxe+1)
    dmaxi = min(pmaxi,2*jmaxi+1)

    ! Usefull precomputations
    inv_Nx = 1._dp/REAL(Nx,dp)
    inv_Ny = 1._dp/REAL(Ny,dp)

  END SUBROUTINE grid_readinputs

  SUBROUTINE init_1Dgrid_distr
    ! write(*,*) Nx
    local_nky        = (Ny/2+1)/num_procs_ky
    ! write(*,*) local_nkx
    local_nky_offset = rank_ky*local_nky
    if (rank_ky .EQ. num_procs_ky-1) local_nky = (Ny/2+1)-local_nky_offset
  END SUBROUTINE init_1Dgrid_distr

  SUBROUTINE set_pgrid
    USE prec_const
    USE model, ONLY: beta ! To know if we solve ampere or not and put odd  p moments
    IMPLICIT NONE
    INTEGER :: ip, istart, iend, in

    ! If no parallel dim (Nz=1) and no EM effects (beta=0), the moment hierarchy
    !! is separable between odds and even P and since the energy is injected in
    !! P=0 and P=2 for density/temperature gradients there is no need of
    !! simulating the odd p which will only be damped.
    !! We define in this case a grid Parray = 0,2,4,...,Pmax i.e. deltap = 2
    !! instead of 1 to spare computation
    IF((Nz .EQ. 1) .AND. (beta .EQ. 0._dp)) THEN
      deltape = 2; deltapi = 2;
      pp2     = 1; ! index p+2 is ip+1
    ELSE
      deltape = 1; deltapi = 1;
      pp2     = 2; ! index p+2 is ip+2
    ENDIF

    ! Total number of Hermite polynomials we will evolve
    total_np_e = (Pmaxe/deltape) + 1
    total_np_i = (Pmaxi/deltapi) + 1
    Np_e = total_np_e ! Reduced names (redundant)
    Np_i = total_np_i
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
    ! Ghosts boundaries
    ipgs_e = ips_e - 2/deltape; ipge_e = ipe_e + 2/deltape;
    ipgs_i = ips_i - 2/deltapi; ipge_i = ipe_i + 2/deltapi;
    ! List of shift and local numbers between the different processes (used in scatterv and gatherv)
    ALLOCATE(rcv_p_e (1:num_procs_p))
    ALLOCATE(rcv_p_i (1:num_procs_p))
    ALLOCATE(dsp_p_e (1:num_procs_p))
    ALLOCATE(dsp_p_i (1:num_procs_p))
    DO in = 0,num_procs_p-1
      CALL decomp1D(total_np_e, num_procs_p, in, istart, iend)
      rcv_p_e(in+1) = iend-istart+1
      dsp_p_e(in+1) = istart-1
      CALL decomp1D(total_np_i, num_procs_p, in, istart, iend)
      rcv_p_i(in+1) = iend-istart+1
      dsp_p_i(in+1) = istart-1
    ENDDO

    !! local grid computations
    ! Flag to avoid unnecessary logical operations
    CONTAINS_ip0_e = .FALSE.; CONTAINS_ip1_e = .FALSE.
    CONTAINS_ip2_e = .FALSE.; CONTAINS_ip3_e = .FALSE.
    CONTAINS_ip0_i = .FALSE.; CONTAINS_ip1_i = .FALSE.
    CONTAINS_ip2_i = .FALSE.; CONTAINS_ip3_i = .FALSE.
    SOLVE_POISSON  = .FALSE.; SOLVE_AMPERE   = .FALSE.
    ALLOCATE(parray_e(ipgs_e:ipge_e))
    ALLOCATE(parray_i(ipgs_i:ipge_i))
    DO ip = ipgs_e,ipge_e
      parray_e(ip) = (ip-1)*deltape
      ! Storing indices of particular degrees for fluid moments computations
      SELECT CASE (parray_e(ip))
      CASE(0); ip0_e = ip; CONTAINS_ip0_e = .TRUE.
      CASE(1); ip1_e = ip; CONTAINS_ip1_e = .TRUE.
      CASE(2); ip2_e = ip; CONTAINS_ip2_e = .TRUE.
      CASE(3); ip3_e = ip; CONTAINS_ip3_e = .TRUE.
      END SELECT
    END DO
    DO ip = ipgs_i,ipge_i
      parray_i(ip) = (ip-1)*deltapi
      ! Storing indices of particular degrees for fluid moments computations
      SELECT CASE (parray_i(ip))
      CASE(0); ip0_i = ip; CONTAINS_ip0_i = .TRUE.
      CASE(1); ip1_i = ip; CONTAINS_ip1_i = .TRUE.
      CASE(2); ip2_i = ip; CONTAINS_ip2_i = .TRUE.
      CASE(3); ip3_i = ip; CONTAINS_ip3_i = .TRUE.
      END SELECT
    END DO
    IF(CONTAINS_ip0_e .AND. CONTAINS_ip0_i) SOLVE_POISSON = .TRUE.
    IF(CONTAINS_ip1_e .AND. CONTAINS_ip1_i) SOLVE_AMPERE  = .TRUE.
    !DGGK operator uses moments at index p=2 (ip=3) for the p=0 term so the
    ! process that contains ip=1 MUST contain ip=3 as well for both e and i.
    IF(((ips_e .EQ. ip0_e) .OR. (ips_i .EQ. ip0_e)) .AND. ((ipe_e .LT. ip2_e) .OR. (ipe_i .LT. ip2_i)))&
     WRITE(*,*) "Warning : distribution along p may not work with DGGK"
    ! Precomputations
    pmaxe_dp   = real(pmaxe,dp)
    pmaxi_dp   = real(pmaxi,dp)

    ! Overwrite SOLVE_AMPERE flag if beta is zero
    IF(beta .EQ. 0._dp) THEN
      SOLVE_AMPERE = .FALSE.
    ENDIF
  END SUBROUTINE set_pgrid

  SUBROUTINE set_jgrid
    USE prec_const
    IMPLICIT NONE
    INTEGER :: ij
    ! Total number of J degrees
    Nj_e = jmaxe+1
    Nj_i = jmaxi+1
    ! Build the full grids on process 0 to diagnose it without comm
    ALLOCATE(jarray_e_full(1:jmaxe+1))
    ALLOCATE(jarray_i_full(1:jmaxi+1))
    ! J
    DO ij = 1,jmaxe+1; jarray_e_full(ij) = (ij-1); END DO
    DO ij = 1,jmaxi+1; jarray_i_full(ij) = (ij-1); END DO
    ! Local data
    ijs_e = 1; ije_e = jmaxe + 1
    ijs_i = 1; ije_i = jmaxi + 1
    ! Ghosts boundaries
    ijgs_e = ijs_e - 1; ijge_e = ije_e + 1;
    ijgs_i = ijs_i - 1; ijge_i = ije_i + 1;
    ! Local number of J
    local_nj_e = ijge_e - ijgs_e + 1
    local_nj_i = ijge_i - ijgs_i + 1
    ALLOCATE(jarray_e(ijgs_e:ijge_e))
    ALLOCATE(jarray_i(ijgs_i:ijge_i))
    DO ij = ijgs_e,ijge_e; jarray_e(ij) = ij-1; END DO
    DO ij = ijgs_i,ijge_i; jarray_i(ij) = ij-1; END DO
    ! Precomputations
    maxj       = MAX(jmaxi, jmaxe)
    jmaxe_dp   = real(jmaxe,dp)
    jmaxi_dp   = real(jmaxi,dp)
    ! j=0 indices
    DO ij = ijs_e,ije_e; IF(jarray_e(ij) .EQ. 0) ij0_e = ij; END DO
    DO ij = ijs_i,ije_i; IF(jarray_i(ij) .EQ. 0) ij0_i = ij; END DO
  END SUBROUTINE set_jgrid


  SUBROUTINE set_kygrid
    USE prec_const
    USE model, ONLY: LINEARITY, N_HD
    IMPLICIT NONE
    INTEGER :: in, istart, iend
    Nky = Ny/2+1 ! Defined only on positive kx since fields are real
    ! Grid spacings
    IF (Ny .EQ. 1) THEN
      deltaky = 2._dp*PI/Ly
      ky_max  = deltaky
      ky_min  = deltaky
    ELSE
      deltaky = 2._dp*PI/Ly
      ky_max  = (Nky-1)*deltaky
      ky_min  = deltaky
    ENDIF
    ! Build the full grids on process 0 to diagnose it without comm
    ALLOCATE(kyarray_full(1:Nky))
    DO iky = 1,Nky
     kyarray_full(iky) = REAL(iky-1,dp) * deltaky
    END DO
    !! Parallel distribution
    ikys = local_nky_offset + 1
    ikye = ikys + local_nky - 1
    ALLOCATE(kyarray(ikys:ikye))
    local_kymax = 0._dp
    ! List of shift and local numbers between the different processes (used in scatterv and gatherv)
    ALLOCATE(counts_nky (1:num_procs_ky))
    ALLOCATE(displs_nky (1:num_procs_ky))
    DO in = 0,num_procs_ky-1
      CALL decomp1D(Nky, num_procs_ky, in, istart, iend)
      counts_nky(in+1) = iend-istart+1
      displs_nky(in+1) = istart-1
    ENDDO
    ! Creating a grid ordered as dk*(0 1 2 3)
    DO iky = ikys,ikye
      IF(Ny .EQ. 1) THEN
        kyarray(iky)      = deltaky
        kyarray_full(iky) = deltaky
        SINGLE_KY         = .TRUE.
      ELSE
        kyarray(iky) = REAL(iky-1,dp) * deltaky
      ENDIF
      ! Finding kx=0
      IF (kyarray(iky) .EQ. 0) THEN
        iky_0 = iky
        contains_ky0 = .true.
      ENDIF
      ! Finding local kxmax value
      IF (ABS(kyarray(iky)) .GT. local_kymax) THEN
        local_kymax = ABS(kyarray(iky))
      ENDIF
      ! Finding kxmax idx
      IF (kyarray(iky) .EQ. ky_max) THEN
        iky_max = iky
        contains_kymax = .true.
      ENDIF
    END DO
    ! Orszag 2/3 filter
    two_third_kymax = 2._dp/3._dp*deltaky*(Nky-1)
    ! For hyperdiffusion
    IF(LINEARITY.EQ.'linear') THEN
      diff_ky_coeff= (1._dp/ky_max)**N_HD
    ELSE
      diff_ky_coeff= (1._dp/two_third_kymax)**N_HD
    ENDIF

    ALLOCATE(AA_y(ikys:ikye))
    DO iky = ikys,ikye
      IF ( (kyarray(iky) .LT. two_third_kymax) .OR. (LINEARITY .EQ. 'linear')) THEN
        AA_y(iky) = 1._dp;
      ELSE
        AA_y(iky) = 0._dp;
      ENDIF
    END DO
  END SUBROUTINE set_kygrid

  SUBROUTINE set_kxgrid(shear)
    USE prec_const
    USE model, ONLY: LINEARITY, N_HD
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: shear
    REAL    :: Lx_adapted
    IF(shear .GT. 0) THEN
      IF(my_id.EQ.0) write(*,*) 'Magnetic shear detected: set up sheared kx grid..'
      ! mininal size of box in x to respect dkx = 2pi shear dky
      Lx_adapted = Ly/(2._dp*pi*shear*Npol)
      ! Put Nexc to 0 so that it is computed from a target value Lx
      IF(Nexc .EQ. 0) THEN
        Nexc = CEILING(0.9 * Lx/Lx_adapted)
        IF(my_id.EQ.0) write(*,*) 'Adapted Nexc =', Nexc
      ENDIF
      ! x length is adapted
      Lx = Lx_adapted*Nexc
    ENDIF
    Nkx = Nx;
    ! Local data
    ! Start and END indices of grid
    ikxs = 1
    ikxe = Nkx
    local_nkx = ikxe - ikxs + 1
    ALLOCATE(kxarray(ikxs:ikxe))
    ALLOCATE(kxarray_full(1:Nkx))
    IF (Nx .EQ. 1) THEN
      deltakx         = 1._dp
      kxarray(1)      = 0._dp
      ikx_0           = 1
      contains_kx0    = .true.
      kx_max          = 0._dp
      ikx_max         = 1
      kx_min          = 0._dp
      kxarray_full(1) = 0._dp
      local_kxmax     = 0._dp
    ELSE ! Build apprpopriate grid
      deltakx      = 2._dp*PI/Lx
      IF(MODULO(Nkx,2) .EQ. 0) THEN ! Even number of Nkx (-2 -1 0 1 2 3)
        kx_max = (Nkx/2)*deltakx
        kx_min = -kx_max+deltakx
        ! Creating a grid ordered as dk*(0 1 2 3 -2 -1)
        local_kxmax = 0._dp
        DO ikx = ikxs,ikxe
          kxarray(ikx) = deltakx*(MODULO(ikx-1,Nkx/2)-Nkx/2*FLOOR(2.*real(ikx-1)/real(Nkx)))
          if (ikx .EQ. Nx/2+1)     kxarray(ikx) = -kxarray(ikx)
          ! Finding kx=0
          IF (kxarray(ikx) .EQ. 0) THEN
            ikx_0 = ikx
            contains_kx0 = .true.
          ENDIF
          ! Finding local kxmax
          IF (ABS(kxarray(ikx)) .GT. local_kxmax) THEN
            local_kxmax = ABS(kxarray(ikx))
          ENDIF
          ! Finding kxmax
          IF (kxarray(ikx) .EQ. kx_max) ikx_max = ikx
        END DO
        ! Build the full grids on process 0 to diagnose it without comm
        ! kx
        DO ikx = 1,Nkx
            kxarray_full(ikx) = deltakx*(MODULO(ikx-1,Nkx/2)-Nkx/2*FLOOR(2.*real(ikx-1)/real(Nkx)))
            IF (ikx .EQ. Nx/2+1) kxarray_full(ikx) = -kxarray_full(ikx)
        END DO
      ELSE ! Odd number of kx (-2 -1 0 1 2)
        kx_max = (Nkx-1)/2*deltakx
        kx_min = -kx_max
        ! Creating a grid ordered as dk*(0 1 2 -2 -1)
        local_kxmax = 0._dp
        DO ikx = ikxs,ikxe
          IF(ikx .LE. (Nkx-1)/2+1) THEN
            kxarray(ikx) = deltakx*(ikx-1)
          ELSE
            kxarray(ikx) = deltakx*(ikx-Nkx-1)
          ENDIF
          ! Finding kx=0
          IF (kxarray(ikx) .EQ. 0) THEN
            ikx_0 = ikx
            contains_kx0 = .true.
          ENDIF
          ! Finding local kxmax
          IF (ABS(kxarray(ikx)) .GT. local_kxmax) THEN
            local_kxmax = ABS(kxarray(ikx))
          ENDIF
          ! Finding kxmax
          IF (kxarray(ikx) .EQ. kx_max) ikx_max = ikx
        END DO
        ! Build the full grids on process 0 to diagnose it without comm
        ! kx
        DO ikx = 1,Nkx
          IF(ikx .LE. (Nkx-1)/2+1) THEN
            kxarray_full(ikx) = deltakx*(ikx-1)
          ELSE
            kxarray_full(ikx) = deltakx*(ikx-Nkx-1)
          ENDIF
        END DO
      ENDIF
    ENDIF
    ! Orszag 2/3 filter
    two_third_kxmax = 2._dp/3._dp*kx_max;

    ! For hyperdiffusion
    IF(LINEARITY.EQ.'linear') THEN
      diff_kx_coeff= (1._dp/kx_max)**N_HD
    ELSE
      diff_kx_coeff= (1._dp/two_third_kxmax)**N_HD
    ENDIF

    ! Antialiasing filter
    ALLOCATE(AA_x(ikxs:ikxe))
    DO ikx = ikxs,ikxe
      IF ( ((kxarray(ikx) .GT. -two_third_kxmax) .AND. &
           (kxarray(ikx) .LT. two_third_kxmax))   .OR. (LINEARITY .EQ. 'linear')) THEN
        AA_x(ikx) = 1._dp;
      ELSE
        AA_x(ikx) = 0._dp;
      ENDIF
    END DO
  END SUBROUTINE set_kxgrid


  SUBROUTINE set_zgrid
    USE prec_const
    USE model, ONLY: mu_z
    IMPLICIT NONE
    REAL    :: grid_shift, Lz, zmax, zmin
    INTEGER :: istart, iend, in
    total_nz = Nz
    ! Length of the flux tube (in ballooning angle)
    Lz         = 2_dp*pi*Npol
    ! Z stepping (#interval = #points since periodic)
    deltaz        = Lz/REAL(Nz,dp)
    inv_deltaz    = 1._dp/deltaz
    ! Parallel hyperdiffusion coefficient
    IF(mu_z .GT. 0) THEN
      diff_dz_coeff = REAL((imagu*deltaz/2._dp)**4) ! adaptive fourth derivative (~GENE)
    ELSE
      diff_dz_coeff = 1._dp    ! non adaptive (positive sign to compensate mu_z neg)
    ENDIF
    IF (SG) THEN
      grid_shift = deltaz/2._dp
    ELSE
      grid_shift = 0._dp
    ENDIF
    ! Build the full grids on process 0 to diagnose it without comm
    ALLOCATE(zarray_full(1:Nz))
    IF (Nz .EQ. 1) Npol = 0
    zmax = 0; zmin = 0;
    DO iz = 1,total_nz ! z in [-pi pi-dz] x Npol
      zarray_full(iz) = REAL(iz-1,dp)*deltaz - Lz/2._dp
      IF(zarray_full(iz) .GT. zmax) zmax = zarray_full(iz)
      IF(zarray_full(iz) .LT. zmin) zmin = zarray_full(iz)
    END DO
    !! Parallel data distribution
    ! Local data distribution
    CALL decomp1D(total_nz, num_procs_z, rank_z, izs, ize)
    local_nz = ize - izs + 1
    ! Ghosts boundaries (depend on the order of z operators)
    IF(Nz .EQ. 1) THEN
      izgs = izs;     izge = ize;
      zarray_full(izs) = 0;
    ELSEIF(Nz .GE. 4) THEN
      izgs = izs - 2; izge = ize + 2;
    ELSE
      ERROR STOP 'Error stop: Nz is not appropriate!!'
    ENDIF
    ! List of shift and local numbers between the different processes (used in scatterv and gatherv)
    ALLOCATE(counts_nz (1:num_procs_z))
    ALLOCATE(displs_nz (1:num_procs_z))
    DO in = 0,num_procs_z-1
      CALL decomp1D(total_nz, num_procs_z, in, istart, iend)
      counts_nz(in+1) = iend-istart+1
      displs_nz(in+1) = istart-1
    ENDDO
    ! Local z array
    ALLOCATE(zarray(izgs:izge,0:1))
    DO iz = izgs,izge
      IF(iz .EQ. 0) THEN
        zarray(iz,0) = zarray_full(total_nz)
        zarray(iz,1) = zarray_full(total_nz) + grid_shift
      ELSEIF(iz .EQ. -1) THEN
        zarray(iz,0) = zarray_full(total_nz-1)
        zarray(iz,1) = zarray_full(total_nz-1) + grid_shift
      ELSEIF(iz .EQ. total_nz + 1) THEN
        zarray(iz,0) = zarray_full(1)
        zarray(iz,1) = zarray_full(1) + grid_shift
      ELSEIF(iz .EQ. total_nz + 2) THEN
        zarray(iz,0) = zarray_full(2)
        zarray(iz,1) = zarray_full(2) + grid_shift
      ELSE
        zarray(iz,0) = zarray_full(iz)
        zarray(iz,1) = zarray_full(iz) + grid_shift
      ENDIF
    ENDDO
    IF(abs(zarray(izs,0) - zmin) .LT. EPSILON(zmin)) &
      contains_zmin = .TRUE.
    IF(abs(zarray(ize,0) - zmax) .LT. EPSILON(zmax)) &
      contains_zmax = .TRUE.
    ! Weitghs for Simpson rule
    ALLOCATE(zweights_SR(izs:ize))
    DO iz = izs,ize
      IF(MODULO(iz,2) .EQ. 1) THEN ! odd iz
        zweights_SR(iz) = 2._dp
      ELSE ! even iz
        zweights_SR(iz) = 4._dp
      ENDIF
    ENDDO
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
    CALL attach(fidres, TRIM(str), "Nexc", Nexc)
    CALL attach(fidres, TRIM(str),   "Ny",   Ny)
    CALL attach(fidres, TRIM(str),   "Ly",   Ly)
    CALL attach(fidres, TRIM(str),   "Nz",   Nz)
    CALL attach(fidres, TRIM(str),  "Nkx",  Nkx)
    CALL attach(fidres, TRIM(str),  "Lkx",  Lkx)
    CALL attach(fidres, TRIM(str),  "Nky",  Nky)
    CALL attach(fidres, TRIM(str),  "Lky",  Lky)
    CALL attach(fidres, TRIM(str),   "SG",   SG)
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
