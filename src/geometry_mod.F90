module geometry
! computes geometrical quantities
! Adapted from B.J.Frei MOLIX code (2021)

  use prec_const
  use model
  use grid
  use array
  use fields
  use basic
  use calculus, ONLY: simpson_rule_z
  use miller, ONLY: set_miller_parameters, get_miller
implicit none
  PRIVATE
  ! Geometry input parameters
  CHARACTER(len=16), &
               PUBLIC, PROTECTED :: geom
  REAL(dp),    PUBLIC, PROTECTED :: q0        = 1.4_dp  ! safety factor
  REAL(dp),    PUBLIC, PROTECTED :: shear     = 0._dp   ! magnetic field shear
  REAL(dp),    PUBLIC, PROTECTED :: eps       = 0.18_dp ! inverse aspect ratio
  REAL(dp),    PUBLIC, PROTECTED :: alpha_MHD = 0 ! shafranov shift effect alpha = -q2 R dbeta/dr
  ! parameters for Miller geometry
  REAL(dp),    PUBLIC, PROTECTED :: kappa     = 1._dp ! elongation
  REAL(dp),    PUBLIC, PROTECTED :: s_kappa   = 0._dp ! r normalized derivative skappa = r/kappa dkappa/dr
  REAL(dp),    PUBLIC, PROTECTED :: delta     = 0._dp ! triangularity
  REAL(dp),    PUBLIC, PROTECTED :: s_delta   = 0._dp ! '' sdelta = r/sqrt(1-delta2) ddelta/dr
  REAL(dp),    PUBLIC, PROTECTED :: zeta      = 0._dp ! squareness
  REAL(dp),    PUBLIC, PROTECTED :: s_zeta    = 0._dp ! '' szeta = r dzeta/dr

  ! GENE unused additional parameters for miller_mod
  REAL(dp), PUBLIC, PROTECTED :: edge_opt      = 0 ! meant to redistribute the points in z
  REAL(dp), PUBLIC, PROTECTED :: major_R       = 1 ! major radius
  REAL(dp), PUBLIC, PROTECTED :: major_Z       = 0 ! vertical elevation
  REAL(dp), PUBLIC, PROTECTED :: dpdx_pm_geom  = 0 ! amplitude mag. eq. pressure grad.
  REAL(dp), PUBLIC, PROTECTED ::          C_y  = 0 ! defines y coordinate : Cy (q theta - phi)
  REAL(dp), PUBLIC, PROTECTED ::         C_xy  = 0 ! defines x coordinate : B = Cxy Vx x Vy

  ! Geometrical auxiliary variables
  LOGICAL,     PUBLIC, PROTECTED :: SHEARED  = .false. ! flag for shear magn. geom or not
  ! Curvature
  REAL(dp),    PUBLIC, DIMENSION(:,:,:,:), ALLOCATABLE :: Ckxky  ! dimensions: kx, ky, z, odd/even p
  ! Jacobian
  REAL(dp),    PUBLIC, DIMENSION(:,:), ALLOCATABLE :: Jacobian ! dimensions: z, odd/even p
  COMPLEX(dp), PUBLIC, PROTECTED        :: iInt_Jacobian ! Inverse integrated Jacobian
  ! Metric
  REAL(dp),    PUBLIC, DIMENSION(:,:), ALLOCATABLE :: gxx, gxy, gxz, gyy, gyz, gzz ! dimensions: z, odd/even p
  REAL(dp),    PUBLIC, DIMENSION(:,:), ALLOCATABLE :: dxdr, dxdZ, Rc, phic, Zc
  ! derivatives of magnetic field strength
  REAL(dp),    PUBLIC, DIMENSION(:,:), ALLOCATABLE :: gradxB, gradyB, gradzB
  ! Relative magnetic field strength
  REAL(dp),    PUBLIC, DIMENSION(:,:), ALLOCATABLE :: hatB, hatB_NL
  ! Relative strength of major radius
  REAL(dp),    PUBLIC, DIMENSION(:,:), ALLOCATABLE :: hatR, hatZ
  ! Some geometrical coefficients
  REAL(dp),    PUBLIC, DIMENSION(:,:) , ALLOCATABLE :: gradz_coeff  ! 1 / [ J_{xyz} \hat{B} ]
  ! Array to map the index of mode (kx,ky,-pi) to (kx+2pi*s*ky,ky,pi) for sheared periodic boundary condition
  INTEGER,     PUBLIC, DIMENSION(:,:), ALLOCATABLE :: ikx_zBC_L, ikx_zBC_R

  ! Functions
  PUBLIC :: geometry_readinputs, geometry_outputinputs,&
            eval_magnetic_geometry, set_ikx_zBC_map
CONTAINS


  SUBROUTINE geometry_readinputs
    ! Read the input parameters
    IMPLICIT NONE
    NAMELIST /GEOMETRY/ geom, q0, shear, eps,&
      kappa, s_kappa,delta, s_delta, zeta, s_zeta ! For miller
    READ(lu_in,geometry)
    IF(shear .NE. 0._dp) SHEARED = .true.

  END SUBROUTINE geometry_readinputs

  subroutine eval_magnetic_geometry
    ! evalute metrix, elementwo_third_kpmaxts, jacobian and gradient
    implicit none
    REAL(dp) :: kx,ky
    COMPLEX(dp), DIMENSION(izs:ize) :: integrant
    real(dp) :: G1,G2,G3,Cx,Cy

    ! Allocate arrays
    CALL geometry_allocate_mem
    !
    IF( (Ny .EQ. 1) .AND. (Nz .EQ. 1)) THEN !1D perp linear run
      IF( my_id .eq. 0 ) WRITE(*,*) '1D perpendicular geometry'
      call eval_1D_geometry
    ELSE
      SELECT CASE(geom)
        CASE('s-alpha')
          IF( my_id .eq. 0 ) WRITE(*,*) 's-alpha-B geometry'
          call eval_salphaB_geometry
        CASE('Z-pinch')
          IF( my_id .eq. 0 ) WRITE(*,*) 'Z-pinch geometry'
          call eval_zpinch_geometry
          SHEARED = .FALSE.
          shear   = 0._dp
        CASE('miller')
          IF( my_id .eq. 0 ) WRITE(*,*) 'Miller geometry'
          call set_miller_parameters(kappa,s_kappa,delta,s_delta,zeta,s_zeta)
          call get_miller(eps,major_R,major_Z,q0,shear,alpha_MHD,edge_opt,&
               C_y,C_xy,dpdx_pm_geom,gxx,gyy,gzz,gxy,gxz,gyz,&
               gradxB,gradyB,hatB,jacobian,gradzB,hatR,hatZ,dxdR,dxdZ,&
               Ckxky,gradz_coeff)
        CASE DEFAULT
          STOP 'geometry not recognized!!'
        END SELECT
    ENDIF
    !
    ! Evaluate perpendicular wavenumber
    !  k_\perp^2 = g^{xx} k_x^2 + 2 g^{xy}k_x k_y + k_y^2 g^{yy}
    !  normalized to rhos_
    DO eo = 0,1
      DO iz = izgs,izge
        DO iky = ikys, ikye
          ky = kyarray(iky)
          DO ikx = ikxs, ikxe
            kx = kxarray(ikx)
            kparray(iky, ikx, iz, eo) = &
            SQRT( gxx(iz,eo)*kx**2 + 2._dp*gxy(iz,eo)*kx*ky + gyy(iz,eo)*ky**2)/hatB(iz,eo)
            ! there is a factor 1/B from the normalization; important to match GENE
          ENDDO
        ENDDO
      ENDDO
      ! Curvature operator (Frei et al. 2022 eq 2.15)
      DO iz = izgs,izge
        G1 = gxy(iz,eo)*gxy(iz,eo)-gxx(iz,eo)*gyy(iz,eo)
        G2 = gxy(iz,eo)*gxz(iz,eo)-gxx(iz,eo)*gyz(iz,eo)
        G3 = gyy(iz,eo)*gxz(iz,eo)-gxy(iz,eo)*gyz(iz,eo)
        Cx = (G1*gradyB(iz,eo) + G2*gradzB(iz,eo))/hatB(iz,eo)
        Cy = (G3*gradzB(iz,eo) - G1*gradxB(iz,eo))/hatB(iz,eo)

        DO iky = ikys, ikye
          ky = kyarray(iky)
           DO ikx= ikxs, ikxe
             kx = kxarray(ikx)
             Ckxky(iky, ikx, iz,eo) = (Cx*kx + Cy*ky)
           ENDDO
        ENDDO
        ! coefficient in the front of parallel derivative
        gradz_coeff(iz,eo) = 1._dp / jacobian(iz,eo) / hatB(iz,eo)
        ! Factor in front of the nonlinear term
        hatB_NL(iz,eo) = Jacobian(iz,eo)&
            *(gxx(iz,eo)*gyy(iz,eo) - gxy(iz,eo)**2)/hatB(iz,eo)
      ENDDO
    ENDDO


    ! set the mapping for parallel boundary conditions
    CALL set_ikx_zBC_map

    two_third_kpmax = 2._dp/3._dp * MAXVAL(kparray)
    !
    ! Compute the inverse z integrated Jacobian (useful for flux averaging)
    integrant = Jacobian(izs:ize,0) ! Convert into complex array
    CALL simpson_rule_z(integrant,iInt_Jacobian)
    iInt_Jacobian = 1._dp/iInt_Jacobian ! reverse it
  END SUBROUTINE eval_magnetic_geometry
  !
  !--------------------------------------------------------------------------------
  !

  SUBROUTINE eval_salphaB_geometry
  ! evaluate s-alpha geometry model
  implicit none
  REAL(dp) :: z, kx, ky, Gx, Gy
  alpha_MHD = 0._dp

  parity: DO eo = 0,1
  zloop: DO iz = izgs,izge
    z = zarray(iz,eo)

    ! metric
      gxx(iz,eo) = 1._dp
      gxy(iz,eo) = shear*z - alpha_MHD*SIN(z)
      gxz(iz,eo) = 0._dp
      gyy(iz,eo) = 1._dp + (shear*z - alpha_MHD*SIN(z))**2
      gyz(iz,eo) = 1._dp/eps
      gzz(iz,eo) = 0._dp
      dxdR(iz,eo)= COS(z)
      dxdZ(iz,eo)= SIN(z)

    ! Relative strengh of radius
      hatR(iz,eo) = 1._dp + eps*COS(z)
      hatZ(iz,eo) = 1._dp + eps*SIN(z)

    ! toroidal coordinates
      Rc  (iz,eo) = hatR(iz,eo)
      phic(iz,eo) = z
      Zc  (iz,eo) = hatZ(iz,eo)

    ! Jacobian
      Jacobian(iz,eo) = q0*hatR(iz,eo)

    ! Relative strengh of modulus of B
      hatB(iz,eo) = 1._dp / hatR(iz,eo)

    ! Derivative of the magnetic field strenght
      gradxB(iz,eo) = -COS(z) ! Gene put a factor hatB^2 or 1/hatR^2 in this
      gradyB(iz,eo) =  0._dp
      gradzB(iz,eo) =  eps*SIN(z)/hatR(iz,eo)  ! Gene put a factor hatB or 1/hatR in this

    ! Curvature operator
    Gx =           (gxz(iz,eo) * gxy(iz,eo) - gxx(iz,eo) * gyz(iz,eo)) *eps*SIN(Z) ! Kx
    Gy = -COS(z) + (gxz(iz,eo) * gyy(iz,eo) - gxy(iz,eo) * gyz(iz,eo)) *eps*SIN(Z) ! Ky
    DO iky = ikys, ikye
        ky = kyarray(iky)
         DO ikx= ikxs, ikxe
           kx = kxarray(ikx)
           Ckxky(iky, ikx, iz,eo) = (-SIN(z)*kx - COS(z)*ky -(shear*z-alpha_MHD*SIN(z))*SIN(z)*ky)/ hatB(iz,eo)
           ! Ckxky(iky, ikx, iz,eo) = (Gx*kx + Gy*ky) * hatB(iz,eo) ! .. multiply by hatB to cancel the 1/ hatB factor in moments_eqs_rhs.f90 routine
         ENDDO
      ENDDO
    ! coefficient in the front of parallel derivative
      gradz_coeff(iz,eo) = 1._dp / Jacobian(iz,eo) / hatB(iz,eo)

  ENDDO zloop
  ENDDO parity

  END SUBROUTINE eval_salphaB_geometry
  !
  !--------------------------------------------------------------------------------
  !

  SUBROUTINE eval_zpinch_geometry
  implicit none
  REAL(dp) :: z, kx, ky, alpha_MHD
  alpha_MHD = 0._dp

  parity: DO eo = 0,1
  zloop: DO iz = izgs,izge
    z = zarray(iz,eo)

    ! metric
      gxx(iz,eo) = 1._dp
      gxy(iz,eo) = 0._dp
      gxz(iz,eo) = 0._dp
      gyy(iz,eo) = 1._dp
      gyz(iz,eo) = 0._dp
      gzz(iz,eo) = 1._dp
      dxdR(iz,eo)= COS(z)
      dxdZ(iz,eo)= SIN(z)

    ! Relative strengh of radius
      hatR(iz,eo) = 1._dp
      hatZ(iz,eo) = 1._dp

    ! toroidal coordinates
      Rc  (iz,eo) = hatR(iz,eo)
      phic(iz,eo) = z
      Zc  (iz,eo) = hatZ(iz,eo)

    ! Jacobian
      Jacobian(iz,eo) = 1._dp

    ! Relative strengh of modulus of B
      hatB   (iz,eo) = 1._dp
      hatB_NL(iz,eo) = 1._dp

    ! Derivative of the magnetic field strenght
      gradxB(iz,eo) = -1._dp ! Gene put a factor hatB^2 or 1/hatR^2 in this
      gradyB(iz,eo) = 0._dp
      gradzB(iz,eo) = 0._dp ! Gene put a factor hatB or 1/hatR in this

    ! Curvature operator
      DO iky = ikys, ikye
        ky = kyarray(iky)
         DO ikx= ikxs, ikxe
           kx = kxarray(ikx)
           Ckxky(iky, ikx, iz,eo) = -ky
         ENDDO
      ENDDO
    ! coefficient in the front of parallel derivative
      gradz_coeff(iz,eo) = 1._dp / Jacobian(iz,eo) / hatB(iz,eo)

  ENDDO zloop
  ENDDO parity

  END SUBROUTINE eval_zpinch_geometry
    !
    !--------------------------------------------------------------------------------
    !
  subroutine eval_1D_geometry
    ! evaluate 1D perp geometry model
    implicit none
    REAL(dp) :: z, kx, ky

    parity: DO eo = 0,1
    zloop: DO iz = izs,ize
     z = zarray(iz,eo)

    ! metric
       gxx(iz,eo) = 1._dp
       gxy(iz,eo) = 0._dp
       gyy(iz,eo) = 1._dp

    ! Relative strengh of radius
       hatR(iz,eo) = 1._dp

    ! Jacobian
       Jacobian(iz,eo) = 1._dp

    ! Relative strengh of modulus of B
       hatB(iz,eo) = 1._dp

    ! Curvature operator
       DO iky = ikys, ikye
         ky = kyarray(iky)
          DO ikx= ikxs, ikxe
            kx = kxarray(ikx)
            Ckxky(ikx, iky, iz,eo) = -kx ! .. multiply by hatB to cancel the 1/ hatB factor in moments_eqs_rhs.f90 routine
          ENDDO
       ENDDO
    ! coefficient in the front of parallel derivative
       gradz_coeff(iz,eo) = 1._dp
  ENDDO zloop
  ENDDO parity

   END SUBROUTINE eval_1D_geometry

   !
   !--------------------------------------------------------------------------------
   !

 SUBROUTINE set_ikx_zBC_map
 IMPLICIT NONE
 REAL :: shift, kx_shift
 ! For periodic CHI BC or 0 dirichlet
 LOGICAL :: PERIODIC_CHI_BC = .TRUE.
 ALLOCATE(ikx_zBC_R(ikys:ikye,ikxs:ikxe))
 ALLOCATE(ikx_zBC_L(ikys:ikye,ikxs:ikxe))

 ! No periodic connection for extension of the domain
 IF(Nexc .GT. 1) PERIODIC_CHI_BC = .TRUE.

 !! No shear case (simple id mapping)
 !3            | 1    2    3    4    5    6 |  ky = 3 dky
 !2   ky       | 1    2    3    4    5    6 |  ky = 2 dky
 !1   A        | 1    2    3    4    5    6 |  ky = 1 dky
 !0   | -> kx  | 1____2____3____4____5____6 |  ky = 0 dky
 !(e.g.) kx =    0   0.1  0.2  0.3 -0.2 -0.1  (dkx=free)
 DO iky = ikys,ikye
   DO ikx = ikxs,ikxe
     ikx_zBC_L(iky,ikx) = ikx
     ikx_zBC_R(iky,ikx) = ikx
   ENDDO
 ENDDO
 ! Modify connection map only at border of z
 IF(SHEARED) THEN
   ! connection map BC of the RIGHT boundary (z=pi*Npol-dz) (even NZ)
   !3            | 4    x    x    x    2    3 |  ky = 3 dky
   !2   ky       | 3    4    x    x    1    2 |  ky = 2 dky
   !1   A        | 2    3    4    x    6    1 |  ky = 1 dky
   !0   | -> kx  | 1____2____3____4____5____6 |  ky = 0 dky
   !kx =           0   0.1  0.2  0.3 -0.2 -0.1  (dkx=2pi*shear*npol*dky)

   ! connection map BC of the RIGHT boundary (z=pi*Npol-dz) (ODD NZ)
   !3            | x    x    x    2    3 |  ky = 3 dky
   !2   ky       | 3    x    x    1    2 |  ky = 2 dky
   !1   A        | 2    3    x    5    1 |  ky = 1 dky
   !0   | -> kx  | 1____2____3____4____5 |  ky = 0 dky
   !kx =           0   0.1  0.2 -0.2 -0.1  (dkx=2pi*shear*npol*dky)
   IF(contains_zmax) THEN ! Check if the process is at the end of the FT
     DO iky = ikys,ikye
       shift = 2._dp*PI*shear*kyarray(iky)*Npol
         DO ikx = ikxs,ikxe
           kx_shift = kxarray(ikx) + shift
           ! We use EPSILON() to treat perfect equality case
           IF( ((kx_shift-EPSILON(kx_shift)) .GT. kx_max) .AND. (.NOT. PERIODIC_CHI_BC) )THEN ! outside of the frequ domain
             ikx_zBC_R(iky,ikx) = -99
           ELSE
             ikx_zBC_R(iky,ikx) = ikx+(iky-1)*Nexc
             IF( ikx_zBC_R(iky,ikx) .GT. Nkx ) &
              ikx_zBC_R(iky,ikx) = ikx_zBC_R(iky,ikx) - Nkx
           ENDIF
         ENDDO
     ENDDO
   ENDIF
   ! connection map BC of the LEFT boundary (z=-pi*Npol)
   !3            | x    5    6    1    x    x |  ky = 3 dky
   !2   ky       | 5    6    1    2    x    x |  ky = 2 dky
   !1   A        | 6    1    2    3    x    5 |  ky = 1 dky
   !0   | -> kx  | 1____2____3____4____5____6 |  ky = 0 dky
   !(e.g.) kx =    0   0.1  0.2  0.3 -0.2 -0.1  (dkx=2pi*shear*npol*dky)
   IF(contains_zmin) THEN ! Check if the process is at the start of the FT
     DO iky = ikys,ikye
       shift = 2._dp*PI*shear*kyarray(iky)*Npol
         DO ikx = ikxs,ikxe
           kx_shift = kxarray(ikx) - shift
           ! We use EPSILON() to treat perfect equality case
           IF( ((kx_shift+EPSILON(kx_shift)) .LT. kx_min) .AND. (.NOT. PERIODIC_CHI_BC) ) THEN ! outside of the frequ domain
             ikx_zBC_L(iky,ikx) = -99
           ELSE
             ikx_zBC_L(iky,ikx) = ikx-(iky-1)*Nexc
             IF( ikx_zBC_L(iky,ikx) .LT. 1 ) &
              ikx_zBC_L(iky,ikx) = ikx_zBC_L(iky,ikx) + Nkx
           ENDIF
         ENDDO
     ENDDO
   ENDIF
 ELSE
ENDIF
END SUBROUTINE set_ikx_zBC_map

!
!--------------------------------------------------------------------------------
!

   SUBROUTINE geometry_allocate_mem

       ! Curvature and geometry
       CALL allocate_array( Ckxky,   ikys,ikye, ikxs,ikxe,izgs,izge,0,1)
       CALL allocate_array(   Jacobian,izgs,izge, 0,1)
       CALL allocate_array(        gxx,izgs,izge, 0,1)
       CALL allocate_array(        gxy,izgs,izge, 0,1)
       CALL allocate_array(        gxz,izgs,izge, 0,1)
       CALL allocate_array(        gyy,izgs,izge, 0,1)
       CALL allocate_array(        gyz,izgs,izge, 0,1)
       CALL allocate_array(        gzz,izgs,izge, 0,1)
       CALL allocate_array(     gradxB,izgs,izge, 0,1)
       CALL allocate_array(     gradyB,izgs,izge, 0,1)
       CALL allocate_array(     gradzB,izgs,izge, 0,1)
       CALL allocate_array(       hatB,izgs,izge, 0,1)
       CALL allocate_array(    hatB_NL,izgs,izge, 0,1)
       CALL allocate_array(       hatR,izgs,izge, 0,1)
       CALL allocate_array(       hatZ,izgs,izge, 0,1)
       CALL allocate_array(         Rc,izgs,izge, 0,1)
       CALL allocate_array(       phic,izgs,izge, 0,1)
       CALL allocate_array(         Zc,izgs,izge, 0,1)
       CALL allocate_array(       dxdR,izgs,izge, 0,1)
       CALL allocate_array(       dxdZ,izgs,izge, 0,1)
       call allocate_array(gradz_coeff,izgs,izge, 0,1)
       CALL allocate_array( kparray, ikys,ikye, ikxs,ikxe,izgs,izge,0,1)

   END SUBROUTINE geometry_allocate_mem

   SUBROUTINE geometry_outputinputs(fidres, str)
     ! Write the input parameters to the results_xx.h5 file

     USE futils, ONLY: attach

     USE prec_const
     IMPLICIT NONE

     INTEGER, INTENT(in) :: fidres
     CHARACTER(len=256), INTENT(in) :: str
     CALL attach(fidres, TRIM(str),"geometry",  geom)
     CALL attach(fidres, TRIM(str),      "q0",    q0)
     CALL attach(fidres, TRIM(str),   "shear", shear)
     CALL attach(fidres, TRIM(str),     "eps",   eps)
   END SUBROUTINE geometry_outputinputs
end module geometry
