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

implicit none
  PRIVATE
  ! Geometry parameters
  CHARACTER(len=16), &
               PUBLIC, PROTECTED :: geom
  REAL(dp),    PUBLIC, PROTECTED :: q0       = 1.4_dp  ! safety factor
  REAL(dp),    PUBLIC, PROTECTED :: shear    = 0._dp   ! magnetic field shear
  REAL(dp),    PUBLIC, PROTECTED :: eps      = 0.18_dp ! inverse aspect ratio

  ! Geometrical operators
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
  REAL(dp),    PUBLIC, DIMENSION(:,:), ALLOCATABLE :: hatB
  ! Relative strength of major radius
  REAL(dp),    PUBLIC, DIMENSION(:,:), ALLOCATABLE :: hatR, hatZ
  ! Some geometrical coefficients
  REAL(dp),    PUBLIC, DIMENSION(:,:) , ALLOCATABLE :: gradz_coeff  ! 1 / [ J_{xyz} \hat{B} ]

  ! Functions
  PUBLIC :: geometry_readinputs, geometry_outputinputs, eval_magnetic_geometry

CONTAINS


  SUBROUTINE geometry_readinputs
    ! Read the input parameters
    IMPLICIT NONE
    NAMELIST /GEOMETRY/ geom, q0, shear, eps
    READ(lu_in,geometry)

  END SUBROUTINE geometry_readinputs

  subroutine eval_magnetic_geometry
    ! evalute metrix, elementwo_third_kpmaxts, jacobian and gradient
    implicit none
    REAL(dp) :: kx,ky
    COMPLEX(dp), DIMENSION(izs:ize) :: integrant
    INTEGER :: fid

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
        CASE DEFAULT
          ERROR STOP 'Error stop: geometry not recognized!!'
        END SELECT
    ENDIF
    !
    ! Evaluate perpendicular wavenumber
    !  k_\perp^2 = g^{xx} k_x^2 + 2 g^{xy}k_x k_y + k_y^2 g^{yy}
    !  normalized to rhos_
    DO eo = 0,1
     DO iky = ikys, ikye
       ky = kyarray(iky)
        DO ikx = ikxs, ikxe
          kx = kxarray(ikx)
          DO iz = izgs,izge
             kparray(ikx, iky, iz, eo) = &
              SQRT( gxx(iz,eo)*kx**2 + 2._dp*gxy(iz,eo)*kx*ky + gyy(iz,eo)*ky**2)/hatB(iz,eo)
              ! there is a factor 1/B from the normalization; important to match GENE
          ENDDO
       ENDDO
    ENDDO
    ENDDO

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
  REAL(dp) :: z, kx, ky, alpha_MHD
  alpha_MHD = 0._dp

  parity: DO eo = 0,1
  zloop: DO iz = izgs,izge
    z = zarray(iz,eo)

    ! metric
      gxx(iz,eo) = 1._dp
      gxy(iz,eo) = shear*z - alpha_MHD*SIN(z)
      gxz(iz,eo) = 0._dp
      gyy(iz,eo) = 1._dp + (shear*z - alpha_MHD*SIN(z))**2
      gyz(iz,eo) = 1._dp/(eps + EPSILON(eps)) !avoid 1/0 in Zpinch config
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
      gradyB(iz,eo) = 0._dp
      gradzB(iz,eo) = eps * SIN(z) / hatR(iz,eo) ! Gene put a factor hatB or 1/hatR in this

    ! Curvature operator
      DO iky = ikys, ikye
        ky = kyarray(iky)
         DO ikx= ikxs, ikxe
           kx = kxarray(ikx)
           Ckxky(ikx, iky, iz,eo) = (-SIN(z)*kx - (COS(z) + (shear*z - alpha_MHD*SIN(z))* SIN(z))*ky) * hatB(iz,eo) ! .. multiply by hatB to cancel the 1/ hatB factor in moments_eqs_rhs.f90 routine
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
    ! evaluate s-alpha geometry model
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
        hatB(iz,eo) = 1._dp

      ! Derivative of the magnetic field strenght
        gradxB(iz,eo) = 0._dp ! Gene put a factor hatB^2 or 1/hatR^2 in this
        gradyB(iz,eo) = 0._dp
        gradzB(iz,eo) = 0._dp ! Gene put a factor hatB or 1/hatR in this

      ! Curvature operator
        DO iky = ikys, ikye
          ky = kyarray(iky)
           DO ikx= ikxs, ikxe
             kx = kxarray(ikx)
             Ckxky(ikx, iky, iz,eo) = - ky * hatB(iz,eo) ! .. multiply by hatB to cancel the 1/ hatB factor in moments_eqs_rhs.f90 routine
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

   SUBROUTINE geometry_allocate_mem

       ! Curvature and geometry
       CALL allocate_array( Ckxky,   ikxs,ikxe, ikys,ikye,izgs,izge,0,1)
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
       CALL allocate_array(       hatR,izgs,izge, 0,1)
       CALL allocate_array(       hatZ,izgs,izge, 0,1)
       CALL allocate_array(         Rc,izgs,izge, 0,1)
       CALL allocate_array(       phic,izgs,izge, 0,1)
       CALL allocate_array(         Zc,izgs,izge, 0,1)
       CALL allocate_array(       dxdR,izgs,izge, 0,1)
       CALL allocate_array(       dxdZ,izgs,izge, 0,1)
       call allocate_array(gradz_coeff,izgs,izge, 0,1)
       CALL allocate_array( kparray, ikxs,ikxe, ikys,ikye,izgs,izge,0,1)

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
