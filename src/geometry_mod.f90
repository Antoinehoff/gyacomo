module geometry
! computes geometrical quantities
! Adapted from B.J.Frei MOLIX code (2021)
  use prec_const, ONLY: xp
implicit none
  PRIVATE
  ! Geometry input parameters
  CHARACTER(len=16), &
               PUBLIC, PROTECTED :: geom      = 's-alpha'
  REAL(xp),    PUBLIC, PROTECTED :: q0        = 1.4_xp  ! safety factor
  REAL(xp),    PUBLIC, PROTECTED :: shear     = 0._xp   ! magnetic field shear
  REAL(xp),    PUBLIC, PROTECTED :: eps       = 0.18_xp ! inverse aspect ratio
  REAL(xp),    PUBLIC, PROTECTED :: alpha_MHD = 0._xp   ! shafranov shift effect alpha = -q2 R dbeta/dr
  ! parameters for Miller geometry
  REAL(xp),    PUBLIC, PROTECTED :: kappa     = 1._xp ! elongation (1 for circular)
  REAL(xp),    PUBLIC, PROTECTED :: s_kappa   = 0._xp ! r normalized derivative skappa = r/kappa dkappa/dr
  REAL(xp),    PUBLIC, PROTECTED :: delta     = 0._xp ! triangularity
  REAL(xp),    PUBLIC, PROTECTED :: s_delta   = 0._xp ! '' sdelta = r/sqrt(1-delta2) ddelta/dr
  REAL(xp),    PUBLIC, PROTECTED :: zeta      = 0._xp ! squareness
  REAL(xp),    PUBLIC, PROTECTED :: s_zeta    = 0._xp ! '' szeta = r dzeta/dr
  REAL(xp),    PUBLIC, PROTECTED :: theta1    = 0._xp ! for Miller global
  REAL(xp),    PUBLIC, PROTECTED :: theta2    = 0._xp !
  ! to apply shift in the parallel z-BC if shearless
  REAL(xp),    PUBLIC, PROTECTED :: shift_y   = 0._xp ! for Arno <3
  REAL(xp),    PUBLIC, PROTECTED :: Npol      = 1._xp ! number of poloidal turns (real for 3D zpinch studies)
  LOGICAL ,    PUBLIC, PROTECTED :: PB_PHASE  = .false. ! To activate the parallel boundary phase factor
  ! Chooses the type of parallel BC we use for the unconnected kx modes (active for non-zero shear only)
  !  'periodic'     : Connect a disconnected kx to a mode on the other cadran
  !  'dirichlet'    : Connect a disconnected kx to 0
  !  'disconnected' : Connect all kx to 0
  !  'shearless'    : Connect all kx to itself
  CHARACTER(len=256), &
               PUBLIC, PROTECTED :: parallel_bc = 'dirichlet'
  ! VMEC (stellarator) flux-tube geometry, geom = 'vmec' (see vmec_mod.f90)
  CHARACTER(len=256), &
               PUBLIC, PROTECTED :: vmec_file      = ''      ! HDF5 file from scripts/python_utilities/wout2h5.py
  REAL(xp),    PUBLIC, PROTECTED :: torflux        = 0.5_xp  ! s = psi_t/psi_t,edge of the flux surface
  REAL(xp),    PUBLIC, PROTECTED :: alpha0         = 0._xp   ! field line label alpha = theta_pest - iota zeta
  REAL(xp),    PUBLIC, PROTECTED :: zeta_center    = 0._xp   ! toroidal angle at the centre of the tube
  REAL(xp),    PUBLIC, PROTECTED :: nfield_periods = -1._xp  ! toroidal extent in field periods (<= 0: one full turn)
  INTEGER,     PUBLIC, PROTECTED :: surface_option = 0       ! 0: exact s, 1: nearest half-grid, 2: nearest full-grid surface
  LOGICAL,     PUBLIC, PROTECTED :: z_equal_arc    = .false. ! equal-arc-length z (not implemented yet)
  CHARACTER(len=16), &
               PUBLIC, PROTECTED :: Lref_opt       = 'R'     ! reference length: 'R' (major radius) or 'a' (minor radius)
  CHARACTER(len=16), &
               PUBLIC, PROTECTED :: twist_shift_opt= 'local' ! 'local': metric at the tube ends (Martin et al. 2018), 'global': 2 pi shat iota P/nfp
  ! Twist-and-shift factor Gamma: a mode (kx,ky) at z = +pi connects to (kx + Gamma ky, ky) at z = -pi
  !  Gamma = 2 pi shear Npol Cyq0_x0 for the tokamak models, from the metric for VMEC
  REAL(xp),    PUBLIC, PROTECTED :: twist_shift_fac = 0._xp

  ! GENE unused additional parameters for miller_mod
  REAL(xp), PUBLIC, PROTECTED :: edge_opt      = 0.0_xp ! meant to redistribute the points in z
  REAL(xp), PUBLIC, PROTECTED :: major_R       = 1.0_xp ! major radius
  REAL(xp), PUBLIC, PROTECTED :: major_Z       = 0.0_xp ! vertical elevation
  REAL(xp), PUBLIC, PROTECTED :: xpdx_pm_geom  = 1._xp ! amplitude mag. eq. pressure grad.
  REAL(xp), PUBLIC, PROTECTED ::          C_y  = 0._xp ! defines y coordinate : Cy (q theta - phi)
  REAL(xp), PUBLIC, PROTECTED ::         C_xy  = 1._xp ! defines x coordinate : B = Cxy Vx x Vy
  REAL(xp), PUBLIC            ::       Cyq0_x0 = 1._xp ! factor that is not 1 when non circular geom (Miller)
  ! Geometrical auxiliary variables
  LOGICAL,     PUBLIC, PROTECTED :: SHEARED  = .false. ! flag for shear magn. geom or not
  ! Curvature
  REAL(xp),    PUBLIC, DIMENSION(:,:,:,:), ALLOCATABLE :: Ckxky  ! dimensions: kx, ky, z, odd/even p
  ! Jacobian
  REAL(xp),    PUBLIC, DIMENSION(:,:), ALLOCATABLE :: Jacobian ! dimensions: z, odd/even p
  COMPLEX(xp), PUBLIC, PROTECTED        :: iInt_Jacobian ! Inverse integrated Jacobian
  ! Metric (local arrays)
  REAL(xp),    PUBLIC, DIMENSION(:,:), ALLOCATABLE :: gxx, gxy, gxz, gyy, gyz, gzz ! dimensions: z, odd/even p
  REAL(xp),    PUBLIC, DIMENSION(:,:), ALLOCATABLE :: dxdr, dxdZ, Rc, phic, Zc
  ! derivatives of magnetic field strength
  REAL(xp),    PUBLIC, DIMENSION(:,:), ALLOCATABLE :: dBdx, dBdy, dBdz, dlnBdz
  ! Relative magnetic field strength
  REAL(xp),    PUBLIC, DIMENSION(:,:), ALLOCATABLE :: hatB, inv_hatB2 ! normalized bckg magnetic gradient amplitude and its inv squared
  ! Relative strength of major radius
  REAL(xp),    PUBLIC, DIMENSION(:,:), ALLOCATABLE :: hatR, hatZ
  ! Toroidal and poloidal (VMEC) angles along the tube (= z for the tokamak models)
  REAL(xp),    PUBLIC, DIMENSION(:,:), ALLOCATABLE :: zeta_arr, theta_arr
  ! Some geometrical coefficients
  REAL(xp),    PUBLIC, DIMENSION(:,:) , ALLOCATABLE :: gradz_coeff  ! 1 / [ J_{xyz} \hat{B} ]
  ! full arrays for output
  REAL(xp),    PUBLIC, DIMENSION(:), ALLOCATABLE :: gxx_full, gxy_full, gxz_full, gyy_full, gyz_full, gzz_full
  REAL(xp),    PUBLIC, DIMENSION(:), ALLOCATABLE :: dxdr_full, dxdZ_full, Rc_full, phic_full, Zc_full
  REAL(xp),    PUBLIC, DIMENSION(:), ALLOCATABLE :: dBdx_full, dBdy_full, dBdz_full, dlnBdz_full, hatB_full, hatR_full, hatZ_full, gradz_coeff_full
  ! Array to map the index of mode (kx,ky,-pi) to (kx+2pi*s*ky,ky,pi) for sheared periodic boundary condition
  INTEGER,     PUBLIC, DIMENSION(:,:), ALLOCATABLE :: ikx_zBC_L, ikx_zBC_R
  ! Geometric factor in front of the parallel phi derivative (not implemented)
  ! REAL(xp),    PUBLIC, DIMENSION(:,:), ALLOCATABLE :: Gamma_phipar
  ! pb_phase, for parallel boundary phase, contains the factor that occurs when taking into account
  !   that q0 is defined in the middle of the fluxtube whereas the radial position spans in [0,Lx)
  !   This shift introduces a (-1)^(Nexc*iky) phase change that is included in GENE
  COMPLEX(xp), PUBLIC, DIMENSION(:),   ALLOCATABLE :: pb_phase_L, pb_phase_R

  ! Functions
  PUBLIC :: geometry_readinputs, geometry_outputinputs,&
            eval_magnetic_geometry, set_ikx_zBC_map, evaluate_magn_curv
CONTAINS


  SUBROUTINE geometry_readinputs
    USE basic, ONLY: lu_in, speak
    USE vmec,  ONLY: vmec_read_h5, vmec_select_surface, vmec_setup_tube,&
                     vmec_q0, vmec_shat, vmec_a_minor, vmec_R_major, vmec_rho0
    ! Read the input parameters
    IMPLICIT NONE
    NAMELIST /GEOMETRY/ geom, q0, shear, eps,&
      kappa, s_kappa,delta, s_delta, zeta, s_zeta,& ! For Miller
      theta1, theta2,& ! for Miller global (not ready)
      parallel_bc, shift_y, Npol, PB_PHASE,&
      vmec_file, torflux, alpha0, zeta_center, nfield_periods, surface_option,& ! For VMEC
      z_equal_arc, Lref_opt, twist_shift_opt
    READ(lu_in,geometry)
    PB_PHASE = .false.
    SELECT CASE(geom)
    CASE('vmec','VMEC')
      ! The equilibrium is read here (not in eval_magnetic_geometry) because the
      ! units and grid modules need eps, q0 and shear before the geometry is built
      IF(ABS(Npol-1._xp) .GT. EPSILON(Npol)) &
        ERROR STOP ">> ERROR << Npol must be 1 for VMEC geometry (z in [-pi,pi] spans nfield_periods field periods)"
      IF(z_equal_arc) ERROR STOP ">> ERROR << z_equal_arc is not implemented yet for VMEC geometry"
      CALL vmec_read_h5(vmec_file)
      CALL vmec_select_surface(torflux, surface_option)
      CALL vmec_setup_tube(alpha0, zeta_center, nfield_periods, Lref_opt)
      q0    = REAL(vmec_q0,xp)
      shear = REAL(vmec_shat,xp)
      eps   = REAL(vmec_a_minor*vmec_rho0/vmec_R_major,xp)
    END SELECT
    ! Provisional (the final value follows from the twist-and-shift factor in eval_magnetic_geometry)
    IF(shear .GT. 0._xp) SHEARED = .TRUE.
    SELECT CASE(geom)
    CASE('z-pinch','Z-pinch','zpinch','Zpinch')
      parallel_bc = 'shearless'
    END SELECT
    SELECT CASE(parallel_bc)
      CASE ('dirichlet')
      CASE ('periodic')
      CASE ('cyclic')
      CASE ('shearless')
      CASE ('disconnected')
      CASE DEFAULT
        ERROR STOP '>> ERROR << Parallel BC not recognized'
    END SELECT
    CALL speak('Parallel BC : '//parallel_bc,2)
  END SUBROUTINE geometry_readinputs

  subroutine eval_magnetic_geometry
    USE grid,     ONLY: total_nky, total_nz, local_nkx, local_nky, local_nz, ngz,&
                        set_kparray, nzgrid, zweights_SR, ieven, set_kxgrid, zarray
    USE basic,    ONLY: speak, str
    USE miller,   ONLY: set_miller_parameters, get_miller
    USE circular, ONLY: get_circ
    USE vmec,     ONLY: get_vmec_geometry, vmec_sigma_y, vmec_rho0, vmec_Lref_over_a
    USE calculus, ONLY: simpson_rule_z
    USE model,    ONLY: beta, LINEARITY, N_HD, ExBrate
    USE species,  ONLY: Ptot
    USE prec_const, ONLY: PI
    ! evalute metrix, elementwo_third_kpmaxts, jacobian and gradient
    implicit none
    COMPLEX(xp), DIMENSION(local_nz) :: integrant

    ! Allocate arrays
    CALL geometry_allocate_mem(local_nky,local_nkx,local_nz,ngz,nzgrid)

    ! Set MHD pressure coefficient, flagged by model module's MHD_PD
    alpha_MHD = q0**2*beta*Ptot

    !
    IF( (total_nky .EQ. 1) .AND. (total_nz .EQ. 1)) THEN !1D perp linear run
      CALL speak('1D perpendicular geometry',2)
      call eval_1D_geometry
    ELSE
      SELECT CASE(geom)
        CASE('s-alpha','salpha')
          CALL speak('s-alpha geometry',2)
          IF(FLOOR(Npol) .NE. CEILING(Npol)) ERROR STOP "ERROR STOP: Npol must be integer for s-alpha geometry"
          IF(MODULO(FLOOR(Npol),2) .EQ. 0)   ERROR STOP "Npol must be odd for s-alpha"
          call eval_salpha_geometry
          C_y     = 1._xp
          Cyq0_x0 = 1._xp
        CASE('z-pinch','Z-pinch','zpinch','Zpinch')
          CALL speak('Z-pinch geometry',2)
          call eval_zpinch_geometry
          SHEARED = .FALSE.
          shear   = 0._xp
          q0      = 0._xp
          eps     = 0._xp
          kappa   = 1._xp
          C_y     = 1._xp
          Cyq0_x0 = 1._xp
        CASE('miller','Miller','Miller_global','miller_global')
          CALL speak('Miller geometry',2)
          IF(FLOOR(Npol) .NE. CEILING(Npol)) ERROR STOP "ERROR STOP: Npol must be integer for Miller geometry"
          IF(MODULO(FLOOR(Npol),2) .EQ. 0)   THEN
            write(*,*) "Npol must be odd for Miller, (Npol = ",Npol,")"
            ERROR STOP
          ENDIF
          call set_miller_parameters(kappa,s_kappa,delta,s_delta,zeta,s_zeta,theta1,theta2)
          call get_miller(geom,eps,major_R,major_Z,q0,shear,FLOOR(Npol),alpha_MHD,edge_opt,&
                          C_y,C_xy,Cyq0_x0,xpdx_pm_geom,gxx,gxy,gxz,gyy,gyz,gzz,&
                          dBdx,dBdy,dBdz,hatB,jacobian,hatR,hatZ,dxdR,dxdZ)
        CASE('circular','circ')
          CALL speak('Circular geometry',2)
          IF(FLOOR(Npol) .NE. CEILING(Npol)) ERROR STOP "ERROR STOP: Npol must be integer for circular geometry"
          IF(MODULO(FLOOR(Npol),2) .EQ. 0)   THEN
            write(*,*) "Npol must be odd for circular, (Npol = ",Npol,")"
            ERROR STOP
          ENDIF
          call get_circ(eps,major_R,major_Z,q0,shear,&
                          C_y,C_xy,Cyq0_x0,gxx,gxy,gxz,gyy,gyz,gzz,&
                          dBdx,dBdy,dBdz,hatB,jacobian,hatR,hatZ,dxdR,dxdZ)  
        CASE('vmec','VMEC')
          CALL speak('VMEC flux-tube geometry',2)
          IF(ABS(ExBrate) .GT. 0._xp) ERROR STOP '>> ERROR << ExB shear flow is not supported with VMEC geometry yet'
          CALL get_vmec_geometry(local_nz,ngz,nzgrid,zarray,twist_shift_opt,&
                                 gxx,gxy,gxz,gyy,gyz,gzz,dBdx,dBdy,dBdz,hatB,Jacobian,&
                                 hatR,hatZ,Rc,phic,Zc,dxdR,dxdZ,zeta_arr,theta_arr,twist_shift_fac)
          ! B = C_xy grad x X grad y with C_xy = 1 by construction (see vmec_mod)
          C_xy    = 1._xp
          C_y     = REAL(vmec_sigma_y,xp)*REAL(vmec_rho0/vmec_Lref_over_a,xp) ! dy/dalpha in L_ref units
          Cyq0_x0 = 1._xp
        CASE DEFAULT
          ERROR STOP '>> ERROR << geometry not recognized!!'
        END SELECT
    ENDIF
    ! inv squared of the magnetic gradient bckg amplitude (for fast kperp update)
    inv_hatB2 = 1/hatB/hatB
    ! Twist-and-shift factor Gamma (Delta kx = Gamma ky at the parallel boundary)
    SELECT CASE(geom)
    CASE('vmec','VMEC') ! twist_shift_fac, zeta_arr, theta_arr set by get_vmec_geometry
    CASE DEFAULT        ! tokamak models: 2 pi shear Npol, times the Miller factor
      twist_shift_fac = 2._xp*PI*shear*Npol*Cyq0_x0
      zeta_arr  = zarray
      theta_arr = zarray
    END SELECT
    SHEARED = ABS(twist_shift_fac) .GT. 1e-5_xp
    CALL speak('Twist-and-shift factor Gamma = '//str(twist_shift_fac),1)
    ! Reset kx grid (to account for the twist-and-shift factor of the geometry)
    CALL speak('Reset kx grid according to the twist-and-shift factor',1)
    CALL set_kxgrid(twist_shift_fac,LINEARITY,N_HD) 
    !
    ! Evaluate perpendicular wavenumber
    !  k_\perp^2 = g^{xx} k_x^2 + 2 g^{xy}k_x k_y + k_y^2 g^{yy}
    !  normalized to rhos_
    CALL set_kparray(gxx,gxy,gyy,inv_hatB2)
    ! Evaluate magnetic curvature operator Ckxky
    CALL evaluate_magn_curv
    ! coefficient in the front of parallel derivative
    gradz_coeff = 1._xp /(jacobian*hatB)
    ! d/dz(ln B) to correspond to the formulation in Hoffmann et al. 2023
    dlnBdz      = dBdz/hatB
    !
    ! set the mapping for parallel boundary conditions
    CALL set_ikx_zBC_map
    !
    ! Compute the inverse z integrated Jacobian (useful for flux averaging)
    integrant = Jacobian((1+ngz/2):(local_nz+ngz/2),ieven) ! Convert into complex array
    CALL simpson_rule_z(local_nz,zweights_SR,integrant,iInt_Jacobian)
    iInt_Jacobian = 1._xp/iInt_Jacobian ! reverse it
  END SUBROUTINE eval_magnetic_geometry

  SUBROUTINE evaluate_magn_curv
    USE grid,     ONLY: local_nkx, local_nky, local_nz, ngz,&
                        kxarray, kyarray, nzgrid
    IMPLICIT NONE
    REAL(xp) :: kx,ky
    real(xp) :: Cx, Cy, g_xz, g_yz, g_zz
    INTEGER  :: eo,iz,iky,ikx
    DO eo = 1,nzgrid
      DO iz = 1,local_nz+ngz
        ! !covariant metric
        g_xz = jacobian(iz,eo)**2*(gxy(iz,eo)*gyz(iz,eo)-gyy(iz,eo)*gxz(iz,eo))
        g_yz = jacobian(iz,eo)**2*(gxy(iz,eo)*gxz(iz,eo)-gxx(iz,eo)*gyz(iz,eo))
        g_zz = jacobian(iz,eo)**2*(gxx(iz,eo)*gyy(iz,eo)-gxy(iz,eo)**2)
        ! Formulation of curv. coeff with the covariant metric
        ! note: the minus sign for Cx is also present in GENE geometry module in the
        !       set_curvature routine and is done afterwards with an if statement..
        Cx =-(dBdy(iz,eo) - g_yz/g_zz*dBdz(iz,eo))/hatB(iz,eo)
        Cy = (dBdx(iz,eo) - g_xz/g_zz*dBdz(iz,eo))/hatB(iz,eo)
        DO iky = 1,local_nky
          ky = kyarray(iky)
           DO ikx= 1,local_nkx
             kx = kxarray(iky,ikx)
             Ckxky(iky, ikx, iz,eo) = (Cx*kx + Cy*ky)/C_xy
           ENDDO
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE
  !
  !--------------------------------------------------------------------------------
  !

  SUBROUTINE eval_salpha_geometry
    USE grid, ONLY : local_nz,ngz,zarray,nzgrid
  ! evaluate s-alpha geometry model
  implicit none
  REAL(xp) :: z
  INTEGER  :: iz, eo

  DO eo = 1,nzgrid
   DO iz = 1,local_nz+ngz
    z = zarray(iz,eo)

    ! metric
      gxx(iz,eo) = 1._xp
      gxy(iz,eo) = shear*z - alpha_MHD*SIN(z)
      gxz(iz,eo) = 0._xp
      gyy(iz,eo) = 1._xp + (shear*z - alpha_MHD*SIN(z))**2
      gyz(iz,eo) = 1._xp/eps
      gzz(iz,eo) = 1._xp/eps**2
      dxdR(iz,eo)= COS(z)
      dxdZ(iz,eo)= SIN(z)

    ! Poloidal plane coordinates
      hatR(iz,eo) = 1._xp + eps*COS(z)
      hatZ(iz,eo) = eps*SIN(z)

    ! toroidal coordinates
      Rc  (iz,eo) = hatR(iz,eo)
      phic(iz,eo) = z
      Zc  (iz,eo) = hatZ(iz,eo)

    ! Relative strengh of modulus of B
      hatB(iz,eo) = 1._xp/(1._xp + eps*COS(z))

    ! Jacobian
      Jacobian(iz,eo) = q0/hatB(iz,eo)

    ! Derivative of the magnetic field strenght
      dBdx(iz,eo) = -COS(z)*hatB(iz,eo)**2 ! LB = 1
      dBdy(iz,eo) =  0._xp
      dBdz(iz,eo) =  eps*SIN(z)*hatB(iz,eo)**2

    ! Curvature factor
      C_xy = 1._xp

   ENDDO
  ENDDO
  END SUBROUTINE eval_salpha_geometry
  !
  !--------------------------------------------------------------------------------
  !

  SUBROUTINE eval_zpinch_geometry
  USE grid, ONLY : local_nz,ngz,zarray,nzgrid
  implicit none
  REAL(xp) :: z
  INTEGER  :: iz, eo

  DO eo = 1,nzgrid
   DO iz = 1,local_nz+ngz
    z = zarray(iz,eo)

    ! metric
      gxx(iz,eo) = 1._xp
      gxy(iz,eo) = 0._xp
      gxz(iz,eo) = 0._xp
      gyy(iz,eo) = 1._xp ! 1/R but R is the normalization length
      gyz(iz,eo) = 0._xp
      gzz(iz,eo) = 1._xp
      dxdR(iz,eo)= COS(z)
      dxdZ(iz,eo)= SIN(z)

    ! Relative strengh of radius
      hatR(iz,eo) = 1._xp ! R but R is the normalization length
      hatZ(iz,eo) = 1._xp

    ! toroidal coordinates
      Rc  (iz,eo) = hatR(iz,eo)
      phic(iz,eo) = z
      Zc  (iz,eo) = hatZ(iz,eo)

    ! Jacobian
      Jacobian(iz,eo) = 1._xp ! R but R is the normalization length

    ! Relative strengh of modulus of B
      hatB   (iz,eo) = 1._xp

    ! Derivative of the magnetic field strenght
      dBdx(iz,eo) = -hatB(iz,eo) ! LB = 1
      dBdy(iz,eo) = 0._xp
      dBdz(iz,eo) = 0._xp ! Gene put a factor hatB or 1/hatR in this
  ENDDO
  ENDDO

  ! Curvature factor
    C_xy = 1._xp
  END SUBROUTINE eval_zpinch_geometry
    !
    !--------------------------------------------------------------------------------
    ! NOT TESTED
  subroutine eval_1D_geometry
    USE grid, ONLY : local_nz,ngz,zarray, nzgrid
    ! evaluate 1D perp geometry model
    implicit none
    REAL(xp) :: z
    INTEGER  :: iz, eo
    DO eo = 1,nzgrid
      DO iz = 1,local_nz+ngz
      z = zarray(iz,eo)

      ! metric
      gxx(iz,eo) = 1._xp
      gxy(iz,eo) = 0._xp
      gyy(iz,eo) = 1._xp

      ! Relative strengh of radius
      hatR(iz,eo) = 1._xp

      ! Jacobian
      Jacobian(iz,eo) = 1._xp

      ! Relative strengh of modulus of B
      hatB(iz,eo) = 1._xp

      ENDDO
    ENDDO
   END SUBROUTINE eval_1D_geometry

   !
   !--------------------------------------------------------------------------------
   !

 SUBROUTINE set_ikx_zBC_map
   USE grid,       ONLY: local_nky,total_nkx,contains_zmin,contains_zmax, Nexc,&
                         local_nky_offset, kx_max, kx_min, kyarray, kxarray, deltakx
   USE prec_const, ONLY: imagu, pi
   IMPLICIT NONE
   REAL(xp) :: shift
   INTEGER :: ikx,iky, mn_y, isg, m_src, m_tgt, m_min, m_max
   ALLOCATE(ikx_zBC_L(local_nky,total_nkx))
   ALLOCATE(ikx_zBC_R(local_nky,total_nkx))
   ALLOCATE(pb_phase_L(local_nky))
   ALLOCATE(pb_phase_R(local_nky))
   !! No shear case (simple id mapping) or not at the end of the z domain
   !3            | 1    2    3    4    5    6 |  ky = 3 dky
   !2   ky       | 1    2    3    4    5    6 |  ky = 2 dky
   !1   A        | 1    2    3    4    5    6 |  ky = 1 dky
   !0   | -> kx  | 1____2____3____4____5____6 |  ky = 0 dky
   !(e.g.) kx =    0   0.1  0.2  0.3 -0.2 -0.1  (dkx=free)
   DO iky = 1,local_nky
     DO ikx = 1,total_nkx
       ikx_zBC_L(iky,ikx) = ikx ! connect to itself per default
       ikx_zBC_R(iky,ikx) = ikx
     ENDDO
     pb_phase_L(iky) = 1._xp ! no phase change per default
     pb_phase_R(iky) = 1._xp
   ENDDO
   ! Parallel boundary are not trivial for sheared case and if
   !  the user does not ask explicitly for shearless bc
   IF(SHEARED .AND. (parallel_bc .NE. 'shearless')) THEN
     ! The kx grid satisfies dkx = |Gamma| dky/Nexc (set_kxgrid), so the shift
     ! Gamma*ky is exactly sign(Gamma)*mn_y*Nexc grid steps: the connection and the
     ! out-of-domain test are done on the integer mode numbers m = kx/dkx
     isg   = NINT(SIGN(1._xp, twist_shift_fac))
     m_min = NINT(kx_min/deltakx)
     m_max = NINT(kx_max/deltakx)
     !!!!!!!!!! LEFT PARALLEL BOUNDARY
     ! Modify connection map only at border of z (matters for MPI z-parallelization)
     IF(contains_zmin) THEN ! Check if the process is at the start of the fluxtube
       DO iky = 1,local_nky
        ! get the real mode number (iky starts at 1 and is shifted from paral)
        mn_y = iky-1+local_nky_offset
        ! Shift of kx across the parallel boundary (2 pi shear Npol ky for a tokamak)
        shift = twist_shift_fac*kyarray(iky)
          DO ikx = 1,total_nkx
            ! Usual formula for shifting indices using that dkx = |Gamma|*dky/Nexc
            ! Wrap to the other side of the kx grid with MODULO so that shifts
            ! larger than one grid width still land on a valid index
            ikx_zBC_L(iky,ikx) = MODULO(ikx-isg*mn_y*Nexc-1, total_nkx) + 1
            ! Check if it points out of the kx domain
            m_src = NINT(kxarray(iky,ikx)/deltakx)
            m_tgt = m_src - isg*mn_y*Nexc
            IF( (m_tgt .LT. m_min) .OR. (m_tgt .GT. m_max) ) THEN ! outside of the frequ domain
              SELECT CASE(parallel_bc)
              CASE ('dirichlet','disconnected') ! connected to 0
                ikx_zBC_L(iky,ikx) = -99
              CASE ('periodic')
                ikx_zBC_L(iky,ikx) = ikx
              CASE ('cyclic')! reroute it by cycling through modes
                ikx_zBC_L(iky,ikx) = MODULO(ikx_zBC_L(iky,ikx)-1,total_nkx)+1
              CASE DEFAULT
                ERROR STOP "The parallel BC is not recognized (dirichlet, periodic, cyclic or disconnected)"
              END SELECT
            ENDIF
          ENDDO
          ! phase present in GENE from a shift of the x origin by Lx/2 (useless?)
          ! We also put the user defined shift in the y direction (see Volcokas et al. 2022)
          IF (PB_PHASE) &
            pb_phase_L(iky) = (-1._xp)**(Nexc*mn_y-1)*EXP(imagu*REAL(mn_y,xp)*2._xp*pi*shift_y)
       ENDDO
     ENDIF
     ! Option for disconnecting every modes, viz. connecting all boundary to 0
     IF(parallel_bc .EQ. 'disconnected') ikx_zBC_L = -99
     !!!!!!!!!! RIGHT PARALLEL BOUNDARY
     IF(contains_zmax) THEN ! Check if the process is at the end of the flux-tube
       DO iky = 1,local_nky
        ! get the real mode number (iky starts at 1 and is shifted from paral)
        mn_y = iky-1+local_nky_offset
        ! Shift of kx across the parallel boundary (2 pi shear Npol ky for a tokamak)
        shift = twist_shift_fac*kyarray(iky)
          DO ikx = 1,total_nkx
            ! Usual formula for shifting indices
            ! Wrap to the other side of the kx grid with MODULO so that shifts
            ! larger than one grid width still land on a valid index
            ikx_zBC_R(iky,ikx) = MODULO(ikx+isg*mn_y*Nexc-1, total_nkx) + 1
            ! Check if it points out of the kx domain
            m_src = NINT(kxarray(iky,ikx)/deltakx)
            m_tgt = m_src + isg*mn_y*Nexc
            IF( (m_tgt .LT. m_min) .OR. (m_tgt .GT. m_max) ) THEN ! outside of the frequ domain
              SELECT CASE(parallel_bc)
              CASE ('dirichlet','disconnected') ! connected to 0
                ikx_zBC_R(iky,ikx) = -99
              CASE ('periodic') ! connected to itself as for shearless
                ikx_zBC_R(iky,ikx) = ikx
              CASE ('cyclic')
                ! write(*,*) 'check',ikx,iky, kxarray(iky,ikx) + shift, '>', kx_max
                ikx_zBC_R(iky,ikx) = MODULO(ikx_zBC_R(iky,ikx)-1,total_nkx)+1
              CASE DEFAULT
                ERROR STOP "The parallel BC is not recognized (dirichlet, periodic, cyclic or disconnected)"
              END SELECT
            ENDIF
          ENDDO
          ! phase present in GENE from a shift ofthe x origin by Lx/2 (useless?)
          ! We also put the user defined shift in the y direction (see Volcokas et al. 2022)
          IF (PB_PHASE) &
            pb_phase_R(iky) = (-1._xp)**(Nexc*mn_y-1)*EXP(-imagu*REAL(mn_y,xp)*2._xp*pi*shift_y)
       ENDDO
     ENDIF
     ! Option for disconnecting every modes, viz. connecting all boundary to 0
     IF(parallel_bc .EQ. 'disconnected') ikx_zBC_R = -99
    ENDIF
    ! write(*,*) kxarray
    ! write(*,*) kyarray
    ! write(*,*) 'ikx_zBC_L :-----------'
    ! DO iky = 1,local_nky
    !   print*, 'iky=', iky
    !   print*, ikx_zBC_L(iky,:)
    ! enddo
    ! print*, pb_phase_L
    ! write(*,*) 'ikx_zBC_R :-----------'
    ! DO iky = 1,local_nky
    !   print*, 'iky=', iky
    !   print*, ikx_zBC_R(iky,:)
    ! enddo
    ! print*, pb_phase_R
    ! print*, shift_y
    ! print*, kx_min, kx_max
    ! stop
    !!!!!!! Example of maps ('x' means connected to 0 value, in the table it is -99)
    ! dirichlet connection map BC of the RIGHT boundary (z=pi*Npol-dz)
    !3            | 4    x    x    x    2    3 |  ky = 3 dky
    !2   ky       | 3    4    x    x    1    2 |  ky = 2 dky
    !1   A        | 2    3    4    x    6    1 |  ky = 1 dky
    !0   | -> kx  | 1____2____3____4____5____6 |  ky = 0 dky
    !kx =           0   0.1  0.2  0.3 -0.2 -0.1  (dkx=2pi*shear*npol*dky)

    ! periodic connection map BC of the RIGHT boundary (z=pi*Npol-dz)
    !3            | 4    2    3    4    2    3 |  ky = 3 dky
    !2   ky       | 3    4    3    4    1    2 |  ky = 2 dky
    !1   A        | 2    3    4    4    6    1 |  ky = 1 dky
    !0   | -> kx  | 1____2____3____4____5____6 |  ky = 0 dky
    !kx =           0   0.1  0.2  0.3 -0.2 -0.1  (dkx=2pi*shear*npol*dky)

    ! cyclic connection map BC of the LEFT boundary (z=-pi*Npol)
    !3            | 4    5    6    1    2    3 |  ky = 3 dky
    !2   ky       | 5    6    1    2    3    4 |  ky = 2 dky
    !1   A        | 6    1    2    3    4    5 |  ky = 1 dky
    !0   | -> kx  | 1____2____3____4____5____6 |  ky = 0 dky
    !(e.g.) kx =    0   0.1  0.2  0.3 -0.2 -0.1  (dkx=2pi*shear*npol*dky)

    ! shearless connection map BC of the LEFT/RIGHT boundary (z=+/-pi*Npol)
    !3            | 1    2    3    4    5    6 |  ky = 3 dky
    !2   ky       | 1    2    3    4    5    6 |  ky = 2 dky
    !1   A        | 1    2    3    4    5    6 |  ky = 1 dky
    !0   | -> kx  | 1____2____3____4____5____6 |  ky = 0 dky
    !(e.g.) kx =    0   0.1  0.2  0.3 -0.2 -0.1  (dkx=2pi*shear*npol*dky)

    ! disconnected connection map BC of the LEFT/RIGHT boundary (z=+/-pi*Npol)
    !3            | x    x    x    x    x    x |  ky = 3 dky
    !2   ky       | x    x    x    x    x    x |  ky = 2 dky
    !1   A        | x    x    x    x    x    x |  ky = 1 dky
    !0   | -> kx  | x____x____x____x____x____x |  ky = 0 dky
    !(e.g.) kx =    0   0.1  0.2  0.3 -0.2 -0.1  (dkx=2pi*shear*npol*dky)
END SUBROUTINE set_ikx_zBC_map

!
!--------------------------------------------------------------------------------
!

   SUBROUTINE geometry_allocate_mem(local_nky,local_nkx,local_nz,ngz,nzgrid)
     INTEGER, INTENT(IN) :: local_nky,local_nkx,local_nz,ngz,nzgrid
       ! Curvature and geometry
       ALLOCATE( Ckxky(local_nky,local_nkx,local_nz+ngz,nzgrid))
       ALLOCATE(   Jacobian(local_nz+ngz,nzgrid))
       ALLOCATE(        gxx(local_nz+ngz,nzgrid))
       ALLOCATE(        gxy(local_nz+ngz,nzgrid))
       ALLOCATE(        gxz(local_nz+ngz,nzgrid))
       ALLOCATE(        gyy(local_nz+ngz,nzgrid))
       ALLOCATE(        gyz(local_nz+ngz,nzgrid))
       ALLOCATE(        gzz(local_nz+ngz,nzgrid))
       ALLOCATE(       dBdx(local_nz+ngz,nzgrid))
       ALLOCATE(       dBdy(local_nz+ngz,nzgrid))
       ALLOCATE(       dBdz(local_nz+ngz,nzgrid))
       ALLOCATE(     dlnBdz(local_nz+ngz,nzgrid))
       ALLOCATE(  inv_hatB2(local_nz+ngz,nzgrid))
       ALLOCATE(       hatB(local_nz+ngz,nzgrid))
       ALLOCATE(       hatR(local_nz+ngz,nzgrid))
       ALLOCATE(       hatZ(local_nz+ngz,nzgrid))
       ALLOCATE(         Rc(local_nz+ngz,nzgrid))
       ALLOCATE(       phic(local_nz+ngz,nzgrid))
       ALLOCATE(         Zc(local_nz+ngz,nzgrid))
       ALLOCATE(       dxdR(local_nz+ngz,nzgrid))
       ALLOCATE(       dxdZ(local_nz+ngz,nzgrid))
       ALLOCATE(gradz_coeff(local_nz+ngz,nzgrid))
       ALLOCATE(   zeta_arr(local_nz+ngz,nzgrid))
       ALLOCATE(  theta_arr(local_nz+ngz,nzgrid))

   END SUBROUTINE geometry_allocate_mem

   SUBROUTINE geometry_outputinputs(fid)
     ! Write the input parameters to the results_xx.h5 file
     USE futils, ONLY: attach, creatd
     USE vmec,   ONLY: vmec_nfield_periods, vmec_a_minor, vmec_R_major, vmec_B_ref, vmec_L_ref,&
                       vmec_rho0, vmec_iota, vmec_nfp, vmec_sign_torflux, vmec_sigma_y, vmec_sigma_z
     IMPLICIT NONE
     INTEGER, INTENT(in) :: fid
     CHARACTER(len=256)  :: str
     WRITE(str,'(a)') '/data/input/geometry'
     CALL creatd(fid, 0,(/0/),TRIM(str),'Geometry Input')
     CALL attach(fid, TRIM(str),"geometry",  geom)
     CALL attach(fid, TRIM(str),      "q0",    q0)
     CALL attach(fid, TRIM(str),   "shear", shear)
     CALL attach(fid, TRIM(str),     "eps",   eps)
     CALL attach(fid, TRIM(str),    "Npol",  Npol)
     CALL attach(fid, TRIM(str),"parallel_bc", parallel_bc)
     CALL attach(fid, TRIM(str),"twist_shift_fac", twist_shift_fac)
     SELECT CASE(geom)
     CASE('vmec','VMEC')
       CALL attach(fid, TRIM(str),       "vmec_file", vmec_file)
       CALL attach(fid, TRIM(str),         "torflux", torflux)
       CALL attach(fid, TRIM(str),          "alpha0", alpha0)
       CALL attach(fid, TRIM(str),     "zeta_center", zeta_center)
       CALL attach(fid, TRIM(str),  "nfield_periods", REAL(vmec_nfield_periods,xp))
       CALL attach(fid, TRIM(str),  "surface_option", surface_option)
       CALL attach(fid, TRIM(str),        "Lref_opt", Lref_opt)
       CALL attach(fid, TRIM(str), "twist_shift_opt", twist_shift_opt)
       CALL attach(fid, TRIM(str),         "a_minor", REAL(vmec_a_minor,xp))
       CALL attach(fid, TRIM(str),         "R_major", REAL(vmec_R_major,xp))
       CALL attach(fid, TRIM(str),           "B_ref", REAL(vmec_B_ref,xp))
       CALL attach(fid, TRIM(str),           "L_ref", REAL(vmec_L_ref,xp))
       CALL attach(fid, TRIM(str),            "rho0", REAL(vmec_rho0,xp))
       CALL attach(fid, TRIM(str),            "iota", REAL(vmec_iota,xp))
       CALL attach(fid, TRIM(str),             "nfp", vmec_nfp)
       CALL attach(fid, TRIM(str),    "sign_torflux", vmec_sign_torflux)
       CALL attach(fid, TRIM(str),         "sigma_y", vmec_sigma_y)
       CALL attach(fid, TRIM(str),         "sigma_z", vmec_sigma_z)
     END SELECT
   END SUBROUTINE geometry_outputinputs
end module geometry
