module geometry
! computes geometrical quantities
! Adapted from B.J.Frei MOLIX code (2021)

  use prec_const
  use model
  use grid
  use array
  use fields
  use basic

implicit none
  public

contains

  subroutine eval_magnetic_geometry
    ! evalute metrix, elements, jacobian and gradient
    implicit none
    !
     IF( my_id .eq. 0 ) WRITE(*,*) 'circular geometry'
     ! circular model
     call eval_circular_geometry
    !
  END SUBROUTINE eval_magnetic_geometry
  !
  !--------------------------------------------------------------------------------
  !
  subroutine eval_circular_geometry
    ! evaluate circular geometry model

    ! Ref: Lapilonne et al., PoP, 2009
    implicit none
    REAL(dp) :: z, kx, ky

    zloop: DO iz = izs,ize
      z = zarray(iz)

      ! Metric elements
      gxx(iz) = 1._dp
      gxy(iz) = shear*z -eps*SIN(z)
      gyy(iz) = 1._dp +(shear*z)**2 -2._dp*eps*COS(z) -2._dp*shear*eps*z*SIN(z)

      ! Relative strengh of radius
      hatR(iz) = 1._dp + eps*cos(z)

      hatB(iz) = SQRT(gxx(iz)*gyy(iz) - gxy(iz)**2)

      ! Jacobian
      Jacobian(iz) = q0*hatR(iz)**2

      ! Derivative of the magnetic field strenght
      gradxB(iz) = -(COS(z) + eps*SIN(z)**2)/hatB(iz)/hatB(iz)
      gradzB(iz) =  eps*SIN(z)
      ! coefficient in the front of parallel derivative
      gradz_coeff(iz) = 1._dp / Jacobian(iz) / hatB(iz)

      ! Curvature operator
      DO iky = ikys, ikye
        ky = kyarray(iky)
        DO ikx= ikxs, ikxe
          kx = kxarray(ikx)
          Ckxky(ikx,iky,iz) = -SIN(z)*kx -(COS(z)+shear*z*SIN(z))*ky
        ENDDO
      ENDDO
      IF (Nky .EQ. 1) THEN ! linear 1D run we switch kx and ky for parallel opt
        DO iky = ikys, ikye
          ky = kyarray(iky)
          DO ikx= ikxs, ikxe
            kx = kxarray(ikx)
            Ckxky(ikx,iky,iz) = -SIN(z)*ky -(COS(z)+shear*z*SIN(z))*kx
          ENDDO
        ENDDO
      ENDIF
    ENDDO zloop
  END SUBROUTINE eval_circular_geometry
  !
  !--------------------------------------------------------------------------------
  ! UNUSED
  ! subroutine eval_salpha_geometry
  !  ! evaluate s-alpha geometry model
  !  implicit none
  !  REAL(dp) :: z, kx, ky
  !
  !  zloop: DO iz = izs,ize
  !   z = zarray(iz)
  !   gxx(iz) = 1._dp
  !   gxy(iz) = shear*z
  !   gyy(iz) = 1._dp + (shear*z)**2
  !   gyz(iz) = 1._dp/eps
  !   gxz(iz) = 0._dp
  !
  !  ! Relative strengh of radius
  !     hatR(iz) = 1._dp + eps*cos(z)
  !
  !  ! Jacobian
  !     Jacobian(iz) = q0*hatR(iz)
  !
  !  ! Relative strengh of modulus of B
  !     hatB(iz) = 1._dp / hatR( iz)
  !
  !  ! Derivative of the magnetic field strenght
  !     gradxB(iz) = - cos( z) / hatR(iz)
  !     gradzB(iz) = eps * sin(z)
  !
  !  ! Gemoetrical coefficients for the curvature operator
  !  ! Note: Gamma2 and Gamma3 are obtained directly form Gamma1 in the expression of the curvature operator implemented here
  !  !
  !     Gamma1(iz) = gxy(iz) * gxy(iz) - gxx(iz) * gyy(iz)
  !     Gamma2(iz) = gxz(iz) * gxy(iz) - gxx(iz) * gyz(iz)
  !     Gamma3(iz) = gxz(iz) * gyy(iz) - gxy(iz) * gyz(iz)
  !
  !  ! Curvature operator
  !     DO iky = ikys, ikye
  !       ky = kyarray(iky)
  !        DO ikx= ikxl, ikxr
  !          kx = kxarray(ikx,iky)
  !          Cxy(ikx, iky, iz) = (-sin(z)*kx - (cos(z) + shear* z* sin(z))*ky) * hatB(iz) ! .. multiply by hatB to cancel the 1/ hatB factor in moments_eqs_rhs.f90 routine
  !        ENDDO
  !     ENDDO
  !
  !  ! coefficient in the front of parallel derivative
  !     gradz_coeff(iz) = 1._dp / Jacobian(iz) / hatB(iz)
  !  ENDDO
  !
  !  ! Evaluate perpendicular wavenumber
  !  !  k_\perp^2 = g^{xx} k_x^2 + 2 g^{xy}k_x k_y + k_y^2 g^{yy}
  !  !  normalized to rhos_
  !     DO iky = ikys, ikye
  !       ky = kyarray(iky)
  !        DO ikx = ikxl, ikxr
  !          kx = kxarray(ikx,iky)
  !           kperp_array(ikx, iky, iz) = sqrt( gxx(iz)*kx**2 + 2._dp* gxy(iz) * kx*ky + gyy(iz)* ky**2) !! / hatB( iz) ! there is a factor 1/B from the normalization; important to match GENE
  !        ENDDO
  !     ENDDO
  !  ENDDO zloop
  !  !
  !  IF( me .eq. 0 ) PRINT*, 'max kperp = ', maxval( kperp_array)
  !
  ! END SUBROUTINE eval_salpha_geometry
  !
  !--------------------------------------------------------------------------------
  !
end module geometry
