MODULE grid
  ! Grid module for spatial discretization
  USE prec_const, ONLY: xp
  USE basic,      ONLY: allocate_array, speak
  USE parallel,   ONLY: my_id, comm_ky
  USE iso_c_binding
  IMPLICIT NONE
  PRIVATE

  !   GRID parameters
  INTEGER,  PUBLIC, PROTECTED :: pmax  = 2      ! The maximal Hermite-moment computed
  INTEGER,  PUBLIC, PROTECTED :: jmax  = 1      ! The maximal Laguerre-moment computed
  INTEGER,  PUBLIC, PROTECTED :: dmax  = 2      ! The maximal full GF set of i-moments v^dmax
  INTEGER,  PUBLIC, PROTECTED :: Nx    = 64     ! Number of total internal grid points in x
  REAL(xp), PUBLIC, PROTECTED :: Lx    = 120_xp ! horizontal length of the spatial box
  INTEGER,  PUBLIC, PROTECTED :: Nexc  = 1      ! factor to increase Lx when shear>0 (Lx = Nexc/kymin/shear)
  INTEGER,  PUBLIC, PROTECTED :: Ny    = 64     ! Number of total internal grid points in y
  REAL(xp), PUBLIC, PROTECTED :: Ly    = 120_xp ! vertical length of the spatial box
  INTEGER,  PUBLIC, PROTECTED :: Nz    = 1      ! Number of total perpendicular planes
  !   Modes number (determined from Nx and Ny)
  INTEGER,  PUBLIC, PROTECTED :: Nkx            ! Number of total internal grid points in kx
  INTEGER,  PUBLIC, PROTECTED :: Nky            ! Number of total internal grid points in ky
  ! Grid arrays
  INTEGER,  DIMENSION(:),   ALLOCATABLE, PUBLIC,PROTECTED :: parray,  parray_full
  INTEGER,  DIMENSION(:),   ALLOCATABLE, PUBLIC,PROTECTED :: jarray,  jarray_full
  REAL(xp), DIMENSION(:,:), ALLOCATABLE, PUBLIC,PROTECTED :: kxarray  ! moving kx grid (ExB)
  REAL(xp), DIMENSION(:),   ALLOCATABLE, PUBLIC,PROTECTED :: kxarray0 ! inirial kyarray
  REAL(xp), DIMENSION(:),   ALLOCATABLE, PUBLIC,PROTECTED :: kxarray_full
  REAL(xp), DIMENSION(:),   ALLOCATABLE, PUBLIC,PROTECTED :: xarray
  REAL(xp), DIMENSION(:),   ALLOCATABLE, PUBLIC,PROTECTED :: kyarray, kyarray_full
  REAL(xp), DIMENSION(:),   ALLOCATABLE, PUBLIC,PROTECTED :: ikyarray, inv_ikyarray !mode indices arrays
  REAL(xp), DIMENSION(:,:), ALLOCATABLE, PUBLIC,PROTECTED :: zarray
  REAL(xp), DIMENSION(:),   ALLOCATABLE, PUBLIC,PROTECTED :: zarray_full
  REAL(xp), DIMENSION(:,:,:,:), ALLOCATABLE, PUBLIC,PROTECTED :: kparray  !moving kperp grid 
  REAL(xp), DIMENSION(:,:,:,:), ALLOCATABLE, PUBLIC,PROTECTED :: kparray0 !initial kperp grid
  REAL(xp), DIMENSION(:,:,:,:), ALLOCATABLE, PUBLIC,PROTECTED :: kp2array !moving kperp^2 grid
  ! Kronecker delta for p=0, p=1, p=2, j=0, j=1
  REAL(xp), DIMENSION(:),   ALLOCATABLE, PUBLIC,PROTECTED :: kroneck_p0, kroneck_p1, kroneck_p2, kroneck_p3
  REAL(xp), DIMENSION(:),   ALLOCATABLE, PUBLIC,PROTECTED :: kroneck_j0, kroneck_j1
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
  INTEGER, PUBLIC, PROTECTED :: local_nkx  ! = total_Nkx = Nx, not parallel
  INTEGER, PUBLIC, PROTECTED :: local_nx   ! = local_nx_ptr from FFTW (used only for Fourier)
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
  INTEGER, PUBLIC, PROTECTED :: local_nx_offset
  INTEGER, PUBLIC, PROTECTED :: local_nz_offset
  ! C-pointer type for FFTW3
  integer(C_INTPTR_T), PUBLIC,PROTECTED :: local_nx_ptr, local_nky_ptr
  integer(C_INTPTR_T), PUBLIC,PROTECTED :: local_nx_ptr_offset, local_nky_ptr_offset
  ! Grid spacing and limits
  REAL(xp), PUBLIC, PROTECTED ::  deltap, deltaz, inv_deltaz, inv_dkx
  REAL(xp), PUBLIC, PROTECTED ::  deltakx, deltaky, deltax, kx_max, ky_max, kx_min, ky_min!, kp_max
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
  PUBLIC :: init_grids_data, set_grids, update_grids
  PUBLIC :: set_kxgrid, set_kparray ! Can be called from geometry to include shear effects
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

  !! Init the local and global number of points in all directions
  SUBROUTINE init_grids_data(Na,EM,LINEARITY) 
    USE fourier,  ONLY: init_grid_distr_and_plans
    USE parallel, ONLY: num_procs_p, rank_p, num_procs_ky, rank_ky, num_procs_z, rank_z
    USE utility,  ONLY: decomp1D
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: Na        ! number of species coming from the model module
    LOGICAL, INTENT(IN) :: EM        ! electromagnetic effects (to skip odd Hermite or not)
    CHARACTER(len=*), INTENT(IN) :: LINEARITY    ! Linear or nonlinear run
    INTEGER :: in, istart, iend
    !!----------------- SPECIES INDICES (not parallelized)
    ias = 1
    iae = Na
    total_Na = Na
    local_Na = Na
    local_Na_offset = ias - 1
    !!----------------- HERMITE INDICES (parallelized)
    ! If no parallel dim (Nz=1) and no EM effects (beta=0), the moment hierarchy
    !! is separable between odds and even P
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
    ! Local data distribution
    CALL decomp1D(total_np, num_procs_p, rank_p, ips, ipe)
    local_np        = ipe - ips + 1
    local_np_offset = ips - 1
    ! Allocate the grid arrays
    ALLOCATE(parray_full(total_np))
    ALLOCATE(parray(local_np+ngp))
    !!----------------- LAGUERRE INDICES (not parallelized)
    ! Total number of J degrees
    total_nj   = jmax+1
    local_jmax = jmax
    Ngj= 2      ! 2-points ghosts for j+\-1 terms
    ! Indices of local data
    ijs = 1; ije = jmax + 1
    ! Local number of J
    local_nj        = ije - ijs + 1
    local_nj_offset = ijs - 1
    ! allocate global and local
    ALLOCATE(jarray_full(total_nj))
    ALLOCATE(jarray(local_nj+ngj))
    ALLOCATE(kroneck_j0(local_nj+ngj));
    ALLOCATE(kroneck_j1(local_nj+ngj));
    !!----------------- SPATIAL GRIDS (parallelized)
    !! Parallel distribution of kx ky grid
    IF (LINEARITY .NE. 'linear') THEN ! we let FFTW distribute if we use it
      IF (my_id .EQ. 0) write(*,*) 'FFTW3 y-grid distribution'
      CALL init_grid_distr_and_plans(Nx,Ny,comm_ky,local_nx_ptr,local_nx_ptr_offset,local_nky_ptr,local_nky_ptr_offset)
    ELSE ! otherwise we distribute equally
      IF (my_id .EQ. 0) write(*,*) 'Manual y-grid distribution'
      ! balanced distribution among the processes
      CALL  decomp1D( Ny/2+1, num_procs_ky, rank_ky, ikys, ikye )
      local_nky_ptr        = ikye - ikys + 1
      local_nky_ptr_offset = ikys - 1
    ENDIF
    !!----------------- BINORMAL KY INDICES (parallelized)
    Nky              = Ny/2+1 ! Defined only on positive kx since fields are real
    total_nky        = Nky
    Ngy              = 0 ! no ghosts cells in ky
    ikys             = local_nky_ptr_offset + 1
    ikye             = ikys + local_nky_ptr - 1
    local_nky        = ikye - ikys + 1
    local_nky_offset = local_nky_ptr_offset
    ALLOCATE(kyarray_full(Nky))
    ALLOCATE(kyarray(local_nky))
    ALLOCATE(ikyarray(Nky))
    ALLOCATE(inv_ikyarray(Nky))
    ALLOCATE(AA_y(local_nky))
    !!---------------- RADIAL KX INDICES (not parallelized)
    Nkx              = Nx
    total_nkx        = Nx
    ikxs             = 1
    ikxe             = total_nkx
    local_nkx        = ikxe - ikxs + 1
    local_nkx_offset = ikxs - 1
    ALLOCATE(kxarray_full(total_nkx))
    ALLOCATE(kxarray0(local_Nkx))
    ALLOCATE(kxarray(local_nky,local_Nkx))
    ALLOCATE(AA_x(local_nkx))
    !!---------------- RADIAL X GRID (only for Fourier routines)
    local_nx        = local_nx_ptr
    local_nx_offset = local_nx_ptr_offset
    ALLOCATE(xarray(Nx))
    !!---------------- PARALLEL Z GRID (parallelized)
    total_nz = Nz
    IF (SG) THEN
      CALL speak('--2 staggered z grids--')
      ! indices for even p and odd p grids (used in kernel, jacobian, gij etc.)
      ieven  = 1
      iodd   = 2
      nzgrid = 2
    ELSE
      ieven  = 1
      iodd   = 1
      nzgrid = 1
    ENDIF
    !! Parallel data distribution
    IF( (Nz .EQ. 1) .AND. (num_procs_z .GT. 1) ) &
    ERROR STOP '>> ERROR << Cannot have multiple core in z-direction (Nz = 1)'
    ! Local data distribution
    CALL decomp1D(total_nz, num_procs_z, rank_z, izs, ize)
    local_nz        = ize - izs + 1
    local_nz_offset = izs - 1
    ! Ghosts boundaries (depend on the order of z operators)
    IF(Nz .EQ. 1) THEN
      ngz = 0
    ELSEIF(Nz .GE. 4) THEN
      ngz =4
      IF(mod(Nz,2) .NE. 0 ) THEN
        ERROR STOP '>> ERROR << Nz must be an even number for Simpson integration rule !!!!'
     ENDIF
    ELSE
      ERROR STOP '>> ERROR << Nz is not appropriate!!'
    ENDIF
    ALLOCATE(zarray(local_nz+ngz,nzgrid))
    ALLOCATE(zarray_full(total_nz))
    ALLOCATE(counts_nz (num_procs_z))
    ALLOCATE(displs_nz (num_procs_z))
    CALL allocate_array(local_zmax,1,nzgrid)
    CALL allocate_array(local_zmin,1,nzgrid)
    ALLOCATE(zweights_SR(local_nz))
    ! List of shift and local numbers between the different processes (used in scatterv and gatherv)
    DO in = 0,num_procs_z-1
      CALL decomp1D(total_nz, num_procs_z, in, istart, iend)
      counts_nz(in+1) = iend-istart+1
      displs_nz(in+1) = istart-1
    ENDDO
    !!---------------- Kperp grid
    CALL allocate_array( kparray, 1,local_nky, 1,local_nkx, 1,local_nz+ngz, 1,nzgrid)
    CALL allocate_array(kparray0, 1,local_nky, 1,local_nkx, 1,local_nz+ngz, 1,nzgrid)
    CALL allocate_array(kp2array, 1,local_nky, 1,local_nkx, 1,local_nz+ngz, 1,nzgrid)
  END SUBROUTINE
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
    REAL(xp), INTENT(IN) :: shear, Npol
    CHARACTER(len=*), INTENT(IN) :: LINEARITY
    INTEGER, INTENT(IN)  :: N_HD
    LOGICAL, INTENT(IN)  :: EM
    INTEGER, INTENT(IN)  :: Na
    CALL set_agrid(Na)
    CALL set_pgrid(EM)
    CALL set_jgrid
    CALL set_kygrid(LINEARITY,N_HD)
    CALL set_kxgrid(shear,Npol,1._xp,LINEARITY,N_HD) ! this will be redone after geometry if Cyq0_x0 .NE. 1
    CALL set_xgrid
    CALL set_zgrid (Npol)
  END SUBROUTINE set_grids

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
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: EM
    INTEGER :: ip
    ! Build the full grids on process 0 to diagnose it without comm
    ! P
    DO ip = 1,total_np; parray_full(ip) = (ip-1)*deltap; END DO
    !! local grid computations
    ! Fill pgrid array
    DO ip = 1,local_np+ngp
      parray(ip) = (ip-1-ngp/2+local_np_offset)*deltap
    ENDDO
    local_pmax = parray(local_np+ngp/2)
    local_pmin = parray(1+ngp/2)
    ! Allocate kronecker arrays for p=0,1,2,3
    ALLOCATE(kroneck_p0(local_np+ngp)); kroneck_p0 = 0._xp
    ALLOCATE(kroneck_p1(local_np+ngp)); kroneck_p1 = 0._xp
    ALLOCATE(kroneck_p2(local_np+ngp)); kroneck_p2 = 0._xp
    ALLOCATE(kroneck_p3(local_np+ngp)); kroneck_p3 = 0._xp
    DO ip = 1,local_np+ngp
      SELECT CASE (parray(ip))
      CASE(0)
        ip0            = ip
        kroneck_p0(ip) = 1._xp
      CASE(1)
        ip1            = ip
        kroneck_p1(ip) = 1._xp
      CASE(2)
        ip2            = ip
        kroneck_p2(ip) = 1._xp
      CASE(3)
        ip3            = ip
        kroneck_p3(ip) = 1._xp
      END SELECT
    END DO
    ! Set local flags to avoid unnecessary logical operations
    CONTAINSp0 = .FALSE.; CONTAINSp1 = .FALSE.
    CONTAINSp2 = .FALSE.; CONTAINSp3 = .FALSE.
    DO ip = 1,local_np+ngp
      SELECT CASE (parray(ip))
      CASE(0); CONTAINSp0 = .TRUE.
      CASE(1); CONTAINSp1 = .TRUE.
      CASE(2); CONTAINSp2 = .TRUE.
      CASE(3); CONTAINSp3 = .TRUE.
      END SELECT
    END DO
    ! Flags that sets if Poisson and Ampere have to be solved locally
    SOLVE_POISSON  = .FALSE.; SOLVE_AMPERE   = .FALSE.
    DO ip = 1+ngp/2,local_np+ngp/2
      SELECT CASE (parray(ip))
      CASE(0); SOLVE_POISSON = .TRUE.
      CASE(1); SOLVE_AMPERE  = .TRUE.
      END SELECT
    END DO
    ! Precomputations
    pmax_xp       = real(pmax,xp)
    diff_p_coeff  = pmax_xp*(1._xp/(pmax_xp+1._xp))**6
    ! Overwrite SOLVE_AMPERE flag if beta is zero
    IF(.NOT. EM) THEN
      SOLVE_AMPERE = .FALSE.
    ENDIF
  END SUBROUTINE set_pgrid

  SUBROUTINE set_jgrid
    USE prec_const
    IMPLICIT NONE
    INTEGER :: ij
    ! Build the full grids on process 0 to diagnose it without comm
    ! J
    DO ij = 1,total_nj
      jarray_full(ij) = (ij-1)
    END DO
    DO ij = 1,local_nj+ngj
      jarray(ij) = ij-1-ngj/2+local_nj_offset
    END DO
    local_jmax = jarray(local_nj+ngj/2)
    local_jmin = jarray(1+ngj/2)
    ! Precomputations
    jmax_xp      = real(jmax,xp)
    diff_j_coeff = jmax_xp*(1._xp/(jmax_xp+1._xp))**6
    ! j=0 and j=1 indices
    DO ij = 1,local_nj+ngj
      IF(jarray(ij) .EQ. 0) ij0 = ij
      IF(jarray(ij) .EQ. 1) ij1 = ij
    END DO
    ! Kronecker arrays for j
    kroneck_j0 = 0._xp
    kroneck_j1 = 0._xp
    DO ij = 1,local_nj+ngj
      SELECT CASE(jarray(ij))
      CASE(0)
        kroneck_j0(ij) = 1._xp
      CASE(1)
        kroneck_j1(ij) = 1._xp
      END SELECT  
    END DO
  END SUBROUTINE set_jgrid

  SUBROUTINE set_kygrid(LINEARITY,N_HD)
    USE prec_const
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) ::LINEARITY
    INTEGER, INTENT(IN) :: N_HD
    INTEGER :: iky
    ! Grid spacings
    IF (Ny .EQ. 1) THEN
      ERROR STOP "Gyacomo cannot run with only one ky"
    ELSE
      deltaky = 2._xp*PI/Ly
      ky_max  = (Nky-1)*deltaky
      ky_min  = deltaky
    ENDIF
    ! Build the full grids on process 0 to diagnose it without comm
    DO iky = 1,Nky
     kyarray_full(iky) = REAL(iky-1,xp) * deltaky
     ! full indices grid (for ExB NL factor)
     ikyarray(iky)     = REAL((iky)-1,xp)
     IF(ikyarray(iky) .GT. 0) THEN
       inv_ikyarray(iky) = 1._xp/ikyarray(iky)
     ELSE
       inv_ikyarray(iky) = 0._xp
     ENDIF
    END DO
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
        kyarray(iky)      = kyarray_full(iky+local_nky_offset)
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

  SUBROUTINE set_kxgrid(shear,Npol,Cyq0_x0,LINEARITY,N_HD)
    USE prec_const
    IMPLICIT NONE
    REAL(xp), INTENT(IN) :: shear, Npol, Cyq0_x0
    CHARACTER(len=*), INTENT(IN) ::LINEARITY
    INTEGER, INTENT(IN)  :: N_HD
    INTEGER :: ikx, iky
    REAL(xp):: Lx_adapted
    IF(shear .GT. 0) THEN
      IF(my_id.EQ.0) write(*,*) 'Magnetic shear detected: set up sheared kx grid..'
      ! mininal size of box in x to respect dkx = 2pi shear dky
      Lx_adapted = Ly/(2._xp*pi*shear*Npol*Cyq0_x0)
      ! Put Nexc to 0 so that it is computed from a target value Lx
      IF(Nexc .EQ. 0) THEN
        Nexc = CEILING(0.9 * Lx/Lx_adapted)
        IF(my_id.EQ.0) write(*,*) 'Adapted Nexc =', Nexc
      ENDIF
      ! x length is adapted
      Lx = Lx_adapted*Nexc
    ENDIF
    deltakx = 2._xp*PI/Lx
    inv_dkx = 1._xp/deltakx 
    IF(MODULO(total_nkx,2) .EQ. 0) THEN ! Even number of kx (-2 -1 0 1 2 3)
      ! Creating a grid ordered as dk*(0 1 2 3 -2 -1)
      DO ikx = 1,total_nkx
        kxarray_full(ikx) = deltakx*REAL(MODULO(ikx-1,total_nkx/2)-(total_nkx/2)*FLOOR(2.*real(ikx-1)/real(total_nkx)),xp)
        IF (ikx .EQ. total_nkx/2+1) kxarray_full(ikx) = -kxarray_full(ikx)
      END DO
      kx_max = MAXVAL(kxarray_full)!(total_nkx/2)*deltakx
      kx_min = MINVAL(kxarray_full)!-kx_max+deltakx
      ! Set local grid (not parallelized so same as full one)
      local_kxmax = 0._xp
      DO iky = 1,local_nky
        DO ikx = 1,local_nkx
          kxarray0(ikx)    = kxarray_full(ikx+local_nkx_offset)
          kxarray(iky,ikx) = kxarray_full(ikx+local_nkx_offset)
          ! Finding kx=0
          IF (kxarray(1,ikx) .EQ. 0) THEN
            ikx0 = ikx
            contains_kx0 = .true.
          ENDIF
          ! Finding local kxmax
          IF (ABS(kxarray(iky,ikx)) .GT. local_kxmax) THEN
            local_kxmax = ABS(kxarray(iky,ikx))
            ikx_max = ikx
          ENDIF
        END DO
      END DO
    ELSE ! Odd number of kx (-2 -1 0 1 2)
      ERROR STOP "Gyacomo is safer with an even Kx number"
    ENDIF
    ! Orszag 2/3 filter
    two_third_kxmax = 2._xp/3._xp*(kx_max-deltakx);
    ! Antialiasing filter
    DO iky = 1,local_nky
      DO ikx = 1,local_nkx
        IF ( ((kxarray(iky,ikx) .GT. -two_third_kxmax) .AND. &
            (kxarray(iky,ikx) .LT. two_third_kxmax))   .OR. (LINEARITY .EQ. 'linear')) THEN
          AA_x(ikx) = 1._xp;
        ELSE
          AA_x(ikx) = 0._xp;
        ENDIF
      END DO
    END DO
    ! For hyperdiffusion
    IF(LINEARITY.EQ.'linear') THEN
      diff_kx_coeff= (1._xp/kx_max)**N_HD
    ELSE
      diff_kx_coeff= (1._xp/two_third_kxmax)**N_HD
    ENDIF
  END SUBROUTINE set_kxgrid

  !----------- Radial x grid
  ! used only for the computation  of the ExB shear NL factor
  SUBROUTINE set_xgrid
    USE prec_const, ONLY: xp, pi
    IMPLICIT NONE
    INTEGER :: ix
    REAL    :: L_
    L_     = 2._xp*pi/deltakx
    deltax = L_/REAL(Nx,xp)
    ! full xgrid
    DO ix = 1,Nx
      xarray(ix) = REAL(ix-1,xp)*deltax
    ENDDO
  END SUBROUTINE set_xgrid

  !----------- Parallel z grid
  SUBROUTINE set_zgrid(Npol)
    USE prec_const
    IMPLICIT NONE
    REAL(xp):: grid_shift, Lz, zmax, zmin, Npol
    INTEGER :: iz, ig, eo
    ! Length of the flux tube (in ballooning angle)
    Lz         = 2._xp*pi*Npol
    ! Z stepping (#interval = #points since periodic)
    deltaz        = Lz/REAL(Nz,xp)
    inv_deltaz    = 1._xp/deltaz
    ! Parallel hyperdiffusion coefficient
    diff_dz_coeff = (deltaz/2._xp)**4 ! adaptive fourth derivative (~GENE)
    IF (SG) THEN
      grid_shift = deltaz/2._xp
    ELSE
      grid_shift = 0._xp
    ENDIF
    ! Build the full grids on process 0 to diagnose it without comm   
    IF (Nz .EQ. 1) &
      Npol = 0._xp
    zmax = 0; zmin = 0;
    DO iz = 1,total_nz ! z in [-pi pi-dz] x Npol
      zarray_full(iz) = REAL(iz-1,xp)*deltaz - Lz/2._xp
      IF(zarray_full(iz) .GT. zmax) zmax = zarray_full(iz)
      IF(zarray_full(iz) .LT. zmin) zmin = zarray_full(iz)
    END DO
    IF(Nz .EQ. 1) &
      zarray_full(izs) = 0
    ! Local z array
    DO iz = 1,local_nz
      DO eo = 1,nzgrid
        zarray(iz+ngz/2,eo) = zarray_full(iz+local_nz_offset) + REAL(eo-1,xp)*grid_shift
      ENDDO
    ENDDO
    DO eo = 1,nzgrid
      ! Find local extrema
      local_zmax(eo) = zarray(local_nz+ngz/2,eo)
      local_zmin(eo) = zarray(1+ngz/2,eo)
      ! continue in z<pi and z>pi
      DO ig = 1,ngz/2
        zarray(local_nz+ngz/2+ig,eo) = zarray(local_nz+ngz/2,eo) + ig*deltaz
        zarray(ig,eo) = zarray(1+ngz/2,eo) - (3-ig)*deltaz
      ENDDO
      ! Set up the flags to know if the process contains the tip and/or the tail
      ! of the z domain (important for z-boundary condition)
      IF(abs(local_zmin(eo) - (zmin+REAL(eo-1,xp)*grid_shift)) .LT. EPSILON(zmin)) &
        contains_zmin = .TRUE.
      IF(abs(local_zmax(eo) - (zmax+REAL(eo-1,xp)*grid_shift)) .LT. EPSILON(zmax)) &
        contains_zmax = .TRUE.
    ENDDO
    ! local weights for Simpson rule
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

  SUBROUTINE set_kparray(gxx, gxy, gyy, inv_hatB2)
    IMPLICIT NONE
    REAL(xp), DIMENSION(local_nz+ngz,nzgrid), INTENT(IN) :: gxx,gxy,gyy,inv_hatB2
    INTEGER     :: eo,iz,iky,ikx
    REAL(xp)    :: kx, ky
    DO eo = 1,nzgrid
      DO iz = 1,local_nz+ngz
        DO iky = 1,local_nky
          ky = kyarray(iky)
          DO ikx = 1,local_nkx
            kx = kxarray(iky,ikx)
            ! there is a factor 1/B from the normalization; important to match GENE
            ! this factor comes from $b_a$ argument in the Bessel. Kperp is not used otherwise.
            kp2array(iky, ikx, iz, eo) = &
              (gxx(iz,eo)*kx**2 + 2._xp*gxy(iz,eo)*kx*ky + gyy(iz,eo)*ky**2)*inv_hatB2(iz,eo)

            kparray (iky, ikx, iz, eo)  = SQRT(kp2array(iky, ikx, iz, eo))
            kparray0(iky, ikx, iz, eo)  = SQRT(kp2array(iky, ikx, iz, eo))
            ENDDO
        ENDDO
      ENDDO
    ENDDO
    two_third_kpmax = 2._xp/3._xp * MAXVAL(kparray)
  END SUBROUTINE

  SUBROUTINE update_grids (dkx_ExB,gxx,gxy,gyy,inv_hatB2)
    IMPLICIT NONE
    REAL(xp), DIMENSION(local_nky),INTENT(IN) :: dkx_ExB ! ExB correction dkx = gamma_E ky dtshift
    REAL(xp), DIMENSION(local_nz+ngz,nzgrid), INTENT(IN) :: gxx,gxy,gyy,inv_hatB2
    INTEGER     :: eo,iz,iky,ikx
    REAL(xp)    :: kx, ky
    ! Update the kx grid
    DO ikx = 1,total_Nkx
      DO iky = 1,local_nky
        kxarray(iky,ikx) = kxarray0(ikx) - dkx_ExB(iky)
      ENDDO
    ENDDO
    ! Update the kperp grid
    DO eo = 1,nzgrid
      DO iz = 1,local_nz+ngz
        DO ikx = 1,local_nkx
          DO iky = 1,local_nky
            ky = kyarray(iky)
            kx = kxarray(iky,ikx)
            kp2array(iky,ikx,iz,eo) = (gxx(iz,eo)*kx**2 + 2._xp*gxy(iz,eo)*kx*ky + gyy(iz,eo)*ky**2)*inv_hatB2(iz,eo)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
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

END MODULE grid
