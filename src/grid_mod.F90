MODULE grid
  ! Grid module for spatial discretization
  USE prec_const
  USE basic
  USE parallel, ONLY: my_id, comm_ky
  USE iso_c_binding
  IMPLICIT NONE
  PRIVATE

  !   GRID Input
  INTEGER,  PUBLIC, PROTECTED :: pmax  = 1      ! The maximal Hermite-moment computed
  INTEGER,  PUBLIC, PROTECTED :: jmax  = 1      ! The maximal Laguerre-moment computed
  INTEGER,  PUBLIC, PROTECTED :: maxj  = 1     ! The maximal Laguerre-moment
  INTEGER,  PUBLIC, PROTECTED :: dmax  = 1      ! The maximal full GF set of i-moments v^dmax
  INTEGER,  PUBLIC, PROTECTED :: Nx    = 4      ! Number of total internal grid points in x
  REAL(xp), PUBLIC, PROTECTED :: Lx    = 120_xp ! horizontal length of the spatial box
  INTEGER,  PUBLIC, PROTECTED :: Nexc  = 1      ! factor to increase Lx when shear>0 (Lx = Nexc/kymin/shear)
  INTEGER,  PUBLIC, PROTECTED :: Ny    = 4      ! Number of total internal grid points in y
  REAL(xp), PUBLIC, PROTECTED :: Ly    = 120_xp ! vertical length of the spatial box
  INTEGER,  PUBLIC, PROTECTED :: Nz    = 4      ! Number of total perpendicular planes
  INTEGER,  PUBLIC, PROTECTED :: Odz   = 4      ! order of z interp and derivative schemes
  INTEGER,  PUBLIC, PROTECTED :: Nkx            ! Number of total internal grid points in kx
  INTEGER,  PUBLIC, PROTECTED :: Nky            ! Number of total internal grid points in ky
  REAL(xp), PUBLIC, PROTECTED :: kpar  = 0_xp   ! parallel wave vector component
  ! Grid arrays
  INTEGER,  DIMENSION(:),   ALLOCATABLE, PUBLIC,PROTECTED :: parray,  parray_full
  INTEGER,  DIMENSION(:),   ALLOCATABLE, PUBLIC,PROTECTED :: jarray,  jarray_full
  REAL(xp), DIMENSION(:),   ALLOCATABLE, PUBLIC,PROTECTED :: kxarray, kxarray_full
  REAL(xp), DIMENSION(:),   ALLOCATABLE, PUBLIC,PROTECTED :: kyarray, kyarray_full
  REAL(xp), DIMENSION(:,:), ALLOCATABLE, PUBLIC,PROTECTED :: zarray
  REAL(xp), DIMENSION(:),   ALLOCATABLE, PUBLIC,PROTECTED :: zarray_full
  REAL(xp), DIMENSION(:,:,:,:), ALLOCATABLE, PUBLIC,PROTECTED :: kparray !kperp
  ! Indexation variables
  INTEGER,  PUBLIC, PROTECTED ::  ias ,iae  ! species index
  INTEGER,  PUBLIC, PROTECTED ::  ips ,ipe  ! Hermite
  INTEGER,  PUBLIC, PROTECTED ::  ijs ,ije  ! Laguerre
  INTEGER,  PUBLIC, PROTECTED ::  ikys,ikye ! Fourier y mode
  INTEGER,  PUBLIC, PROTECTED ::  ikxs,ikxe ! Fourier x mode
  INTEGER,  PUBLIC, PROTECTED ::  izs ,ize  ! z-grid
  INTEGER,  PUBLIC, PROTECTED ::  ieven, iodd ! indices for the staggered grids
  INTEGER,  PUBLIC, PROTECTED ::  ip0, ip1, ip2, ip3, ij0, ij1, pp2
  INTEGER,  PUBLIC, PROTECTED ::  ikx0, iky0, ikx_max, iky_max ! Indices of k-grid origin and max
  ! Total numbers of points for Hermite and Laguerre
  INTEGER, PUBLIC, PROTECTED :: total_na
  INTEGER, PUBLIC, PROTECTED :: total_np
  INTEGER, PUBLIC, PROTECTED :: total_nj
  INTEGER, PUBLIC, PROTECTED :: total_nky
  INTEGER, PUBLIC, PROTECTED :: total_nkx
  INTEGER, PUBLIC, PROTECTED :: total_nz
  ! Local numbers of points (without ghosts)
  INTEGER, PUBLIC, PROTECTED :: local_na
  INTEGER, PUBLIC, PROTECTED :: local_np
  INTEGER, PUBLIC, PROTECTED :: local_nj
  INTEGER, PUBLIC, PROTECTED :: local_nky
  INTEGER, PUBLIC, PROTECTED :: local_nkx
  INTEGER, PUBLIC, PROTECTED :: local_nz
  INTEGER, PUBLIC, PROTECTED :: local_nkp
  INTEGER, PUBLIC, PROTECTED :: ngp, ngj, ngx, ngy, ngz ! number of ghosts points
  INTEGER, PUBLIC, PROTECTED :: nzgrid  ! one or two depending on the staggered grid option
  ! Local offsets
  INTEGER, PUBLIC, PROTECTED :: local_na_offset
  INTEGER, PUBLIC, PROTECTED :: local_np_offset
  INTEGER, PUBLIC, PROTECTED :: local_nj_offset
  INTEGER, PUBLIC, PROTECTED :: local_nky_offset
  INTEGER, PUBLIC, PROTECTED :: local_nkx_offset
  INTEGER, PUBLIC, PROTECTED :: local_nz_offset
  ! C-pointer type for FFTW3
  integer(C_INTPTR_T), PUBLIC,PROTECTED :: local_nkx_ptr, local_nky_ptr
  integer(C_INTPTR_T), PUBLIC,PROTECTED :: local_nkx_ptr_offset, local_nky_ptr_offset
  ! Grid spacing and limits
  REAL(xp), PUBLIC, PROTECTED ::  deltap, deltaz, inv_deltaz
  REAL(xp), PUBLIC, PROTECTED ::  deltakx, deltaky, kx_max, ky_max, kx_min, ky_min!, kp_max
  INTEGER , PUBLIC, PROTECTED ::  local_pmin,  local_pmax
  INTEGER , PUBLIC, PROTECTED ::  local_jmin,  local_jmax
  REAL(xp), PUBLIC, PROTECTED ::  local_kymin, local_kymax
  REAL(xp), PUBLIC, PROTECTED ::  local_kxmin, local_kxmax
  REAL(xp), DIMENSION(:), ALLOCATABLE, PUBLIC, PROTECTED ::  local_zmin,  local_zmax
  ! local z weights for computing simpson rule
  REAL(xp),  DIMENSION(:),   ALLOCATABLE, PUBLIC,PROTECTED :: zweights_SR
  ! Numerical diffusion scaling
  REAL(xp), PUBLIC, PROTECTED  ::  diff_p_coeff, diff_j_coeff
  REAL(xp), PUBLIC, PROTECTED  ::  diff_kx_coeff, diff_ky_coeff, diff_dz_coeff
  LOGICAL,  PUBLIC, PROTECTED  ::  SG = .true.! shifted grid flag
  ! Array to know the distribution of data among all processes (for MPI comm)
  INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC,PROTECTED :: counts_total_nkx, counts_nky, counts_nz
  INTEGER, DIMENSION(:), ALLOCATABLE, PUBLIC,PROTECTED :: displs_total_nkx, displs_nky, displs_nz
  ! Kperp array depends on kx, ky, z (geometry), eo (even or odd zgrid)
  LOGICAL,  PUBLIC, PROTECTED ::  contains_kx0   = .false. ! flag if the proc contains kx=0 index
  LOGICAL,  PUBLIC, PROTECTED ::  contains_ky0   = .false. ! flag if the proc contains ky=0 index
  LOGICAL,  PUBLIC, PROTECTED ::  contains_kymax = .false. ! flag if the proc contains kx=kxmax index
  LOGICAL,  PUBLIC, PROTECTED ::  contains_zmax  = .false. ! flag if the proc contains z=pi-dz index
  LOGICAL,  PUBLIC, PROTECTED ::  contains_zmin  = .false. ! flag if the proc contains z=-pi index
  LOGICAL,  PUBLIC, PROTECTED ::       SINGLE_KY = .false. ! to check if it is a single non 0 ky simulation
  ! Grid containing the polynomials degrees
  LOGICAL,  PUBLIC, PROTECTED ::  CONTAINSp0, CONTAINSp1, CONTAINSp2, CONTAINSp3
  LOGICAL,  PUBLIC, PROTECTED ::  SOLVE_POISSON, SOLVE_AMPERE
  ! Usefull inverse numbers
  REAL(xp), PUBLIC, PROTECTED :: inv_Nx, inv_Ny, inv_Nz
  ! For Orszag filter
  REAL(xp), PUBLIC, PROTECTED :: two_third_kxmax
  REAL(xp), PUBLIC, PROTECTED :: two_third_kymax
  REAL(xp), PUBLIC, PROTECTED :: two_third_kpmax
  ! 1D Antialiasing arrays (2/3 rule)
  REAL(xp), DIMENSION(:), ALLOCATABLE, PUBLIC,PROTECTED :: AA_x
  REAL(xp), DIMENSION(:), ALLOCATABLE, PUBLIC,PROTECTED :: AA_y

  ! Public Functions
  PUBLIC :: init_1Dgrid_distr
  PUBLIC :: set_grids, set_kparray
  PUBLIC :: grid_readinputs, grid_outputinputs
  PUBLIC :: bar

  ! Precomputations
  real(xp), PUBLIC, PROTECTED    :: pmax_xp, jmax_xp

CONTAINS


  SUBROUTINE grid_readinputs
    ! Read the input parameters
    USE prec_const
    IMPLICIT NONE
    INTEGER :: lun   = 90              ! File duplicated from STDIN
    NAMELIST /GRID/ pmax, jmax, Nx, Lx, Ny, Ly, Nz, SG, Nexc
    READ(lun,grid)
    IF(Nz .EQ. 1) & ! overwrite SG option if Nz = 1 for safety of use
      SG      = .FALSE.
    !! Compute the maximal degree of full GF moments set
    !   i.e. : all moments N_a^pj s.t. p+2j<=d are simulated (see GF closure)
    dmax = min(pmax,2*jmax+1)
    ! Usefull precomputations
    inv_Nx = 1._xp/REAL(Nx,xp)
    inv_Ny = 1._xp/REAL(Ny,xp)
  END SUBROUTINE grid_readinputs
  !!!! GRID REPRESENTATION
  ! We define the grids that contain ghosts (p,j,z) with indexing from 1 to Nlocal + nghost
  ! e.g. nghost = 4, nlocal = 4
  ! |x|x|_|_|_|_|x|x|
  !  1 2 3 4 5 6 7 8    (x are ghosts)
  ! the interior points are contained between 1+Ng/2 and Nlocal+Ng/2
  ! The other grids are simply
  ! |_|_|_|_|
  !  1 2 3 4
  SUBROUTINE set_grids(shear,Npol,LINEARITY,N_HD,EM,Na)
    USE fourier, ONLY: init_grid_distr_and_plans
    REAL(xp), INTENT(IN) :: shear
    INTEGER,  INTENT(IN) :: Npol
    CHARACTER(len=*), INTENT(IN) :: LINEARITY
    INTEGER, INTENT(IN)  :: N_HD
    LOGICAL, INTENT(IN)  :: EM
    INTEGER, INTENT(IN)  :: Na
    CALL set_agrid(Na)
    CALL set_pgrid(EM)
    CALL set_jgrid
    !! Parallel distribution of kx ky grid
    IF (LINEARITY .NE. 'linear') THEN
      IF (my_id .EQ. 0) write(*,*) 'FFTW3 y-grid distribution'
      CALL init_grid_distr_and_plans(Nx,Ny,comm_ky,local_nkx_ptr,local_nkx_ptr_offset,local_nky_ptr,local_nky_ptr_offset)
    ELSE
      CALL init_1Dgrid_distr
      IF (my_id .EQ. 0) write(*,*) 'Manual y-grid distribution'
    ENDIF
    CALL set_kygrid(LINEARITY,N_HD)
    CALL set_kxgrid(shear,Npol,LINEARITY,N_HD)
    CALL set_zgrid (Npol)
    ! print*, 'p:',parray
    ! print*, 'j:',jarray
    ! print*, 'ky:',kyarray
    ! print*, 'kx:',kxarray
    ! print*, 'z:',zarray
    ! print*, parray(ip0)
    ! print*, jarray(ij0)
    ! print*, kyarray(iky0)
    ! print*, kxarray(ikx0)
  END SUBROUTINE set_grids

  SUBROUTINE init_1Dgrid_distr
    USE parallel, ONLY: num_procs_ky, rank_ky
    ! write(*,*) Nx
    local_nky_ptr        = (Ny/2+1)/num_procs_ky
    ! write(*,*) local_nkx_ptr
    local_nky_ptr_offset = rank_ky*local_nky_ptr
    if (rank_ky .EQ. num_procs_ky-1) local_nky_ptr = (Ny/2+1)-local_nky_ptr_offset
  END SUBROUTINE init_1Dgrid_distr

  SUBROUTINE set_agrid(Na) ! you're a sorcerer harry
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Na
    ias = 1
    iae = Na
    total_Na = Na
    local_Na = Na
    local_Na_offset = ias - 1
  END SUBROUTINE

  SUBROUTINE set_pgrid(EM)
    USE prec_const
    USE parallel, ONLY: num_procs_p, rank_p
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: EM
    INTEGER :: ip
    ! If no parallel dim (Nz=1) and no EM effects (beta=0), the moment hierarchy
    !! is separable between odds and even P and since the energy is injected in
    !! P=0 and P=2 for density/temperature gradients there is no need of
    !! simulating the odd p which will only be damped.
    !! We define in this case a grid Parray = 0,2,4,...,Pmax i.e. deltap = 2
    !! instead of 1 to spare computation
    IF((Nz .EQ. 1) .AND. .NOT. EM) THEN
      deltap  = 2
      Ngp     = 2  ! two ghosts cells for p+/-2 only
      pp2     = 1  ! index p+2 is ip+1
    ELSE
      deltap  = 1
      Ngp     = 4  ! four ghosts cells for p+/-1 and p+/-2 terms
      pp2     = 2  ! index p+2 is ip+2
    ENDIF
    ! Total number of Hermite polynomials we will evolve
    total_np = (Pmax/deltap) + 1
    ! Build the full grids on process 0 to diagnose it without comm
    ALLOCATE(parray_full(total_np))
    ! P
    DO ip = 1,total_np; parray_full(ip) = (ip-1)*deltap; END DO
    !! Parallel data distribution
    ! Local data distribution
    CALL decomp1D(total_np, num_procs_p, rank_p, ips, ipe)
    local_np       = ipe - ips + 1
    local_np_offset = ips - 1
    !! local grid computations
    ! Flag to avoid unnecessary logical operations
    CONTAINSp0 = .FALSE.; CONTAINSp1 = .FALSE.
    CONTAINSp2 = .FALSE.; CONTAINSp3 = .FALSE.
    SOLVE_POISSON  = .FALSE.; SOLVE_AMPERE   = .FALSE.
    ALLOCATE(parray(local_np+ngp))
    ! Fill the interior (no ghosts)
    DO ip = 1,local_np+ngp
      parray(ip) = (ip-1-ngp/2+local_np_offset)*deltap
      ! Storing indices of particular degrees for fluid moments computations
      SELECT CASE (parray(ip))
      CASE(0); ip0 = ip; CONTAINSp0 = .TRUE.
      CASE(1); ip1 = ip; CONTAINSp1 = .TRUE.
      CASE(2); ip2 = ip; CONTAINSp2 = .TRUE.
      CASE(3); ip3 = ip; CONTAINSp3 = .TRUE.
      END SELECT
    END DO
    local_pmax = parray(local_np+ngp/2)
    local_pmin = parray(1+ngp/2)
    IF(CONTAINSp0) SOLVE_POISSON = .TRUE.
    IF(CONTAINSp1) SOLVE_AMPERE  = .TRUE.
    !DGGK operator uses moments at index p=2 (ip=3) for the p=0 term so the
    ! process that contains ip=1 MUST contain ip=3 as well for both e and i.
    IF(CONTAINSp0 .AND. .NOT. (CONTAINSp2))&
     WRITE(*,*) "Warning : distribution along p may not work with DGGK"
    ! Precomputations
    pmax_xp       = real(pmax,xp)
    diff_p_coeff  = pmax_xp*(1._xp/pmax_xp)**6
    ! Overwrite SOLVE_AMPERE flag if beta is zero
    IF(.NOT. EM) THEN
      SOLVE_AMPERE = .FALSE.
    ENDIF
  END SUBROUTINE set_pgrid

  SUBROUTINE set_jgrid
    USE prec_const
    IMPLICIT NONE
    INTEGER :: ij
    ! Total number of J degrees
    total_nj   = jmax+1
    local_jmax = jmax
    Ngj= 2      ! 2-points ghosts for j+\-1 terms
    ! Build the full grids on process 0 to diagnose it without comm
    ALLOCATE(jarray_full(total_nj))
    ! J
    DO ij = 1,total_nj; jarray_full(ij) = (ij-1); END DO
    ! Indices of local data
    ijs = 1; ije = jmax + 1
    ! Local number of J
    local_nj        = ije - ijs + 1
    local_nj_offset = ijs - 1
    ALLOCATE(jarray(local_nj+ngj))
    DO ij = 1,local_nj+ngj
      jarray(ij) = ij-1-ngj/2+local_nj_offset
    END DO
    local_jmax = jarray(local_nj+ngj/2)
    local_jmin = jarray(1+ngj/2)
    ! Precomputations
    jmax_xp      = real(jmax,xp)
    diff_j_coeff = jmax_xp*(1._xp/jmax_xp)**6
    ! j=0 and j=1 indices
    DO ij = 1,local_nj+ngj
      IF(jarray(ij) .EQ. 0) ij0 = ij
      IF(jarray(ij) .EQ. 1) ij1 = ij
    END DO
  END SUBROUTINE set_jgrid

  SUBROUTINE set_kygrid(LINEARITY,N_HD)
    USE prec_const
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) ::LINEARITY
    INTEGER, INTENT(IN) :: N_HD
    INTEGER :: iky
    Nky       = Ny/2+1 ! Defined only on positive kx since fields are real
    total_nky = Nky
    ! Grid spacings
    IF (Ny .EQ. 1) THEN
      deltaky = 2._xp*PI/Ly
      ky_max  = deltaky
      ky_min  = deltaky
    ELSE
      deltaky = 2._xp*PI/Ly
      ky_max  = (Nky-1)*deltaky
      ky_min  = deltaky
    ENDIF
    Ngy = 0 ! no ghosts cells in ky
    ! Build the full grids on process 0 to diagnose it without comm
    ALLOCATE(kyarray_full(Nky))
    DO iky = 1,Nky
     kyarray_full(iky) = REAL(iky-1,xp) * deltaky
    END DO
    ikys = local_nky_ptr_offset + 1
    ikye = ikys + local_nky_ptr - 1
    local_nky = ikye - ikys + 1
    local_nky_offset = local_nky_ptr_offset
    ALLOCATE(kyarray(local_nky))
    local_kymax = 0._xp
    ! Creating a grid ordered as dk*(0 1 2 3)
    ! We loop over the natural iky numbers (|1 2 3||4 5 6||... Nky|)
    DO iky = 1,local_nky
      ! We shift the natural iky index by the offset to obtain the mpi dependent
      ! indexation (|1 2 3||1 2 3|... local_nky|)
      IF(Ny .EQ. 1) THEN
        kyarray(iky)      = deltaky
        kyarray_full(iky) = deltaky
        SINGLE_KY         = .TRUE.
      ELSE
        kyarray(iky) = kyarray_full(iky+local_nky_offset)
      ENDIF
      ! Finding kx=0
      IF (kyarray(iky) .EQ. 0) THEN
        iky0 = iky
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
    two_third_kymax = 2._xp/3._xp*deltaky*(Nky-1)
    ALLOCATE(AA_y(local_nky))
    DO iky = 1,local_nky
      IF ( (kyarray(iky) .LT. two_third_kymax) .OR. (LINEARITY .EQ. 'linear')) THEN
        AA_y(iky) = 1._xp;
      ELSE
        AA_y(iky) = 0._xp;
      ENDIF
    END DO
    ! For hyperdiffusion
    IF(LINEARITY.EQ.'linear') THEN
      diff_ky_coeff= (1._xp/ky_max)**N_HD
    ELSE
      diff_ky_coeff= (1._xp/two_third_kymax)**N_HD
    ENDIF
  END SUBROUTINE set_kygrid

  SUBROUTINE set_kxgrid(shear,Npol,LINEARITY,N_HD)
    USE prec_const
    IMPLICIT NONE
    REAL(xp), INTENT(IN) :: shear
    INTEGER,  INTENT(IN) :: Npol
    CHARACTER(len=*), INTENT(IN) ::LINEARITY
    INTEGER, INTENT(IN)  :: N_HD
    INTEGER :: ikx
    REAL(xp):: Lx_adapted
    IF(shear .GT. 0) THEN
      IF(my_id.EQ.0) write(*,*) 'Magnetic shear detected: set up sheared kx grid..'
      ! mininal size of box in x to respect dkx = 2pi shear dky
      Lx_adapted = Ly/(2._xp*pi*shear*Npol)
      ! Put Nexc to 0 so that it is computed from a target value Lx
      IF(Nexc .EQ. 0) THEN
        Nexc = CEILING(0.9 * Lx/Lx_adapted)
        IF(my_id.EQ.0) write(*,*) 'Adapted Nexc =', Nexc
      ENDIF
      ! x length is adapted
      Lx = Lx_adapted*Nexc
    ENDIF
    Nkx       = Nx
    total_nkx = Nx
    ! Local data
    ! Start and END indices of grid
    ikxs = 1
    ikxe = total_nkx
    local_nkx_ptr = ikxe - ikxs + 1
    local_nkx     = ikxe - ikxs + 1
    local_nkx_offset = ikxs - 1
    ALLOCATE(kxarray(local_nkx))
    ALLOCATE(kxarray_full(total_nkx))
    IF (Nx .EQ. 1) THEN
      deltakx         = 1._xp
      kxarray(1)      = 0._xp
      ikx0           = 1
      contains_kx0    = .true.
      kx_max          = 0._xp
      ikx_max         = 1
      kx_min          = 0._xp
      kxarray_full(1) = 0._xp
      local_kxmax     = 0._xp
    ELSE ! Build apprpopriate grid
      deltakx      = 2._xp*PI/Lx
      IF(MODULO(total_nkx,2) .EQ. 0) THEN ! Even number of kx (-2 -1 0 1 2 3)
        kx_max = (total_nkx/2)*deltakx
        kx_min = -kx_max+deltakx
        ! Creating a grid ordered as dk*(0 1 2 3 -2 -1)
        DO ikx = 1,total_nkx
          kxarray_full(ikx) = deltakx*REAL(MODULO(ikx-1,total_nkx/2)-(total_nkx/2)*FLOOR(2.*real(ikx-1)/real(total_nkx)),xp)
          IF (ikx .EQ. total_nkx/2+1) kxarray_full(ikx) = -kxarray_full(ikx)
        END DO
        ! Set local grid (not parallelized so same as full one)
        local_kxmax = 0._xp
        DO ikx = 1,local_nkx
          kxarray(ikx) = kxarray_full(ikx+local_nkx_offset)
          ! Finding kx=0
          IF (kxarray(ikx) .EQ. 0) THEN
            ikx0 = ikx
            contains_kx0 = .true.
          ENDIF
          ! Finding local kxmax
          IF (ABS(kxarray(ikx)) .GT. local_kxmax) THEN
            local_kxmax = ABS(kxarray(ikx))
            ikx_max = ikx
          ENDIF
        END DO
      ELSE ! Odd number of kx (-2 -1 0 1 2)
        kx_max = (total_nkx-1)/2*deltakx
      ENDIF
    ENDIF
    ! Orszag 2/3 filter
    two_third_kxmax = 2._xp/3._xp*kx_max;
    ! Antialiasing filter
    ALLOCATE(AA_x(local_nkx))
    DO ikx = 1,local_nkx
      IF ( ((kxarray(ikx) .GT. -two_third_kxmax) .AND. &
           (kxarray(ikx) .LT. two_third_kxmax))   .OR. (LINEARITY .EQ. 'linear')) THEN
        AA_x(ikx) = 1._xp;
      ELSE
        AA_x(ikx) = 0._xp;
      ENDIF
    END DO
    ! For hyperdiffusion
    IF(LINEARITY.EQ.'linear') THEN
      diff_kx_coeff= (1._xp/kx_max)**N_HD
    ELSE
      diff_kx_coeff= (1._xp/two_third_kxmax)**N_HD
    ENDIF
  END SUBROUTINE set_kxgrid

  SUBROUTINE set_zgrid(Npol)
    USE prec_const
    USE parallel, ONLY: num_procs_z, rank_z
    IMPLICIT NONE
    REAL(xp):: grid_shift, Lz, zmax, zmin
    INTEGER :: istart, iend, in, Npol, iz, ig, eo, iglob
    total_nz = Nz
    ! Length of the flux tube (in ballooning angle)
    Lz         = 2._xp*pi*REAL(Npol,xp)
    ! Z stepping (#interval = #points since periodic)
    deltaz        = Lz/REAL(Nz,xp)
    inv_deltaz    = 1._xp/deltaz
    ! Parallel hyperdiffusion coefficient
    diff_dz_coeff = (deltaz/2._xp)**4 ! adaptive fourth derivative (~GENE)
    IF (SG) THEN
      CALL speak('--2 staggered z grids--')
      grid_shift = deltaz/2._xp
      ! indices for even p and odd p grids (used in kernel, jacobian, gij etc.)
      ieven  = 1
      iodd   = 2
      nzgrid = 2
    ELSE
      grid_shift = 0._xp
      ieven  = 1
      iodd   = 1
      nzgrid = 1
    ENDIF
    ! Build the full grids on process 0 to diagnose it without comm
    ALLOCATE(zarray_full(total_nz))
    IF (Nz .EQ. 1) Npol = 0
    zmax = 0; zmin = 0;
    DO iz = 1,total_nz ! z in [-pi pi-dz] x Npol
      zarray_full(iz) = REAL(iz-1,xp)*deltaz - Lz/2._xp
      IF(zarray_full(iz) .GT. zmax) zmax = zarray_full(iz)
      IF(zarray_full(iz) .LT. zmin) zmin = zarray_full(iz)
    END DO
    !! Parallel data distribution
    IF( (Nz .EQ. 1) .AND. (num_procs_z .GT. 1) ) &
    ERROR STOP '>> ERROR << Cannot have multiple core in z-direction (Nz = 1)'
    ! Local data distribution
    CALL decomp1D(total_nz, num_procs_z, rank_z, izs, ize)
    local_nz        = ize - izs + 1
    local_nz_offset = izs - 1
    ! Ghosts boundaries (depend on the order of z operators)
    IF(Nz .EQ. 1) THEN
      ngz              = 0
      zarray_full(izs) = 0
    ELSEIF(Nz .GE. 4) THEN
      ngz =4
      IF(mod(Nz,2) .NE. 0 ) THEN
        ERROR STOP '>> ERROR << Nz must be an even number for Simpson integration rule !!!!'
     ENDIF
    ELSE
      ERROR STOP '>> ERROR << Nz is not appropriate!!'
    ENDIF
    ! List of shift and local numbers between the different processes (used in scatterv and gatherv)
    ALLOCATE(counts_nz (num_procs_z))
    ALLOCATE(displs_nz (num_procs_z))
    DO in = 0,num_procs_z-1
      CALL decomp1D(total_nz, num_procs_z, in, istart, iend)
      counts_nz(in+1) = iend-istart+1
      displs_nz(in+1) = istart-1
    ENDDO
    ! Local z array
    ALLOCATE(zarray(local_nz+ngz,nzgrid))
    !! interior point loop
    DO iz = 1,local_nz
      DO eo = 1,nzgrid
        zarray(iz+ngz/2,eo) = zarray_full(iz+local_nz_offset) + REAL(eo-1,xp)*grid_shift
      ENDDO
    ENDDO
    CALL allocate_array(local_zmax,1,nzgrid)
    CALL allocate_array(local_zmin,1,nzgrid)
    DO eo = 1,nzgrid
      ! Find local extrema
      local_zmax(eo) = zarray(local_nz+ngz/2,eo)
      local_zmin(eo) = zarray(1+ngz/2,eo)
      ! Fill the ghosts
      ! Continue angles
      ! DO ig = 1,ngz/2
      !   zarray(ig,eo)          = local_zmin(eo)-REAL(ngz/2-(ig-1),xp)*deltaz
      !   zarray(local_nz+ngz/2+ig,eo) = local_zmax(eo)+REAL(ig,xp)*deltaz
      ! ENDDO
      ! Periodic z \in (-pi pi-dz)
      DO ig = 1,ngz/2 ! first ghost cells
        iglob = ig+local_nz_offset-ngz/2
        IF (iglob .LE. 0) &
          iglob = iglob + total_nz
        zarray(ig,eo) = zarray_full(iglob)
      ENDDO
      DO ig = local_nz+ngz/2,local_nz+ngz ! last ghost cells
        iglob = ig+local_nz_offset-ngz/2
        IF (iglob .GT. total_nz) &
          iglob = iglob - total_nz
        zarray(ig,eo) = zarray_full(iglob)
      ENDDO
      ! Set up the flags to know if the process contains the tip and/or the tail
      ! of the z domain (important for z-boundary condition)
      IF(abs(local_zmin(eo) - (zmin+REAL(eo-1,xp)*grid_shift)) .LT. EPSILON(zmin)) &
        contains_zmin = .TRUE.
      IF(abs(local_zmax(eo) - (zmax+REAL(eo-1,xp)*grid_shift)) .LT. EPSILON(zmax)) &
        contains_zmax = .TRUE.
    ENDDO
    ! local weights for Simpson rule
    ALLOCATE(zweights_SR(local_nz))
    IF(total_nz .EQ. 1) THEN
      zweights_SR = 1._xp
    ELSE
      DO iz = 1,local_nz
        IF(MODULO(iz+local_nz_offset,2) .EQ. 1) THEN ! odd iz
          zweights_SR(iz) = onethird*deltaz*2._xp
        ELSE ! even iz
          zweights_SR(iz) = onethird*deltaz*4._xp
        ENDIF
      ENDDO
    ENDIF
  END SUBROUTINE set_zgrid

  SUBROUTINE set_kparray(gxx, gxy, gyy,hatB)
    REAL(xp), DIMENSION(local_nz+ngz,nzgrid), INTENT(IN) :: gxx,gxy,gyy,hatB
    INTEGER     :: eo,iz,iky,ikx
    REAL(xp)    :: kx, ky
    CALL allocate_array( kparray, 1,local_nky, 1,local_nkx, 1,local_nz+ngz, 1,nzgrid)
    DO eo = 1,nzgrid
      DO iz = 1,local_nz+ngz
        DO iky = 1,local_nky
          ky = kyarray(iky)
          DO ikx = 1,local_nkx
            kx = kxarray(ikx)
            ! there is a factor 1/B from the normalization; important to match GENE
            ! this factor comes from $b_a$ argument in the Bessel. Kperp is not used otherwise.
            kparray(iky, ikx, iz, eo) = &
            SQRT( gxx(iz,eo)*kx**2 + 2._xp*gxy(iz,eo)*kx*ky + gyy(iz,eo)*ky**2)/ hatB(iz,eo)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    two_third_kpmax = 2._xp/3._xp * MAXVAL(kparray)
  END SUBROUTINE

  SUBROUTINE grid_outputinputs(fid)
    ! Write the input parameters to the results_xx.h5 file
    USE futils, ONLY: attach, creatd
    IMPLICIT NONE
    INTEGER, INTENT(in) :: fid
    CHARACTER(len=256)  :: str
    WRITE(str,'(a)') '/data/input/grid'
    CALL creatd(fid, 0,(/0/),TRIM(str),'Grid Input')
    CALL attach(fid, TRIM(str),  "pmax", pmax)
    CALL attach(fid, TRIM(str),"deltap", deltap)
    CALL attach(fid, TRIM(str),    "Np",  total_np)
    CALL attach(fid, TRIM(str),  "jmax", jmax)
    CALL attach(fid, TRIM(str),    "Nj",  total_nj)
    CALL attach(fid, TRIM(str),   "Nkx",  Nkx)
    CALL attach(fid, TRIM(str),    "Nx",   Nx)
    CALL attach(fid, TRIM(str),    "Lx",   Lx)
    CALL attach(fid, TRIM(str),  "Nexc", Nexc)
    CALL attach(fid, TRIM(str),    "Ny",   Ny)
    CALL attach(fid, TRIM(str),   "Nky",  Nky)
    CALL attach(fid, TRIM(str),    "Ly",   Ly)
    CALL attach(fid, TRIM(str),    "Nz",   Nz)
    CALL attach(fid, TRIM(str),   "total_nkx",  total_nkx)
    CALL attach(fid, TRIM(str),   "Nky",  Nky)
    CALL attach(fid, TRIM(str),    "SG",   SG)
  END SUBROUTINE grid_outputinputs

  FUNCTION bar(p_,j_)
    IMPLICIT NONE
    INTEGER :: bar, p_, j_
    bar = (jmax+1)*p_ + j_ + 1
  END FUNCTION

  SUBROUTINE decomp1D( n, numprocs, myid, is, ie )
      INTEGER :: n, numprocs, myid, is, ie
      INTEGER :: nlocal
      INTEGER :: deficit

      nlocal   = n / numprocs
      is        = myid * nlocal + 1
      deficit  = MOD(n,numprocs)
      is        = is + MIN(myid,deficit)
      IF (myid .LT. deficit) nlocal = nlocal + 1
      ie = is + nlocal - 1
      IF ((ie .GT. n) .OR. (myid .EQ. numprocs-1)) ie = n
  END SUBROUTINE decomp1D

END MODULE grid
