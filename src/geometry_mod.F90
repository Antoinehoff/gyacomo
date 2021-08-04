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
     IF( my_id .eq. 0 ) WRITE(*,*) 's-alpha geometry'
     ! circular model
     call eval_s_alpha_geometry
    !
  END SUBROUTINE eval_magnetic_geometry
  !
  !--------------------------------------------------------------------------------
  !
  subroutine eval_s_alpha_geometry
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
  END SUBROUTINE eval_s_alpha_geometry

end module geometry
