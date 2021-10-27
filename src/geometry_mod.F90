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
    REAL(dp) :: kx,ky
    !
    IF( (Ny .EQ. 1) .AND. (Nz .EQ. 1)) THEN !1D perp linear run
      IF( my_id .eq. 0 ) WRITE(*,*) '1D perpendicular geometry'
      call eval_1D_geometry
    ELSE
     IF( my_id .eq. 0 ) WRITE(*,*) 's-alpha-B geometry'
     call eval_salphaB_geometry
    ENDIF
    !
    ! Evaluate perpendicular wavenumber
    !  k_\perp^2 = g^{xx} k_x^2 + 2 g^{xy}k_x k_y + k_y^2 g^{yy}
    !  normalized to rhos_
    DO iz = izs,ize
       DO iky = ikys, ikye
         ky = kyarray(iky)
          DO ikx = ikxs, ikxe
            kx = kxarray(ikx)
             kparray(ikx, iky, iz) = &
              SQRT( gxx(iz)*kx**2 + 2._dp*gxy(iz)*kx*ky + gyy(iz)*ky**2)/hatB(iz)
              ! there is a factor 1/B from the normalization; important to match GENE
          ENDDO
       ENDDO
    ENDDO
  END SUBROUTINE eval_magnetic_geometry
  !
  !--------------------------------------------------------------------------------
  !

  subroutine eval_salphaB_geometry
  ! evaluate s-alpha geometry model
  implicit none
  REAL(dp) :: z, kx, ky

  zloop: DO iz = izs,ize
    z = zarray(iz)

    ! metric
      gxx(iz) = 1._dp
      gxy(iz) = shear*z
      gyy(iz) = 1._dp + (shear*z)**2

    ! Relative strengh of radius
      hatR(iz) = 1._dp + eps*COS(z)

    ! Jacobian
      Jacobian(iz) = q0*hatR(iz)

    ! Relative strengh of modulus of B
      hatB(iz) = 1._dp / hatR(iz)

    ! Derivative of the magnetic field strenght
      gradxB(iz) = -COS(z)
      gradzB(iz) = eps * SIN(z) / hatR(iz)

    ! Curvature operator
      DO iky = ikys, ikye
        ky = kyarray(iky)
         DO ikx= ikxs, ikxe
           kx = kxarray(ikx)
           Ckxky(ikx, iky, iz) = (-SIN(z)*kx - (COS(z) + shear* z* SIN(z))*ky) * hatB(iz) ! .. multiply by hatB to cancel the 1/ hatB factor in moments_eqs_rhs.f90 routine
         ENDDO
      ENDDO
    ! coefficient in the front of parallel derivative
      gradz_coeff(iz) = 1._dp / Jacobian(iz) / hatB(iz)

  ENDDO zloop

  END SUBROUTINE eval_salphaB_geometry
  !
  !--------------------------------------------------------------------------------
  !

  subroutine eval_1D_geometry
    ! evaluate 1D perp geometry model
    implicit none
    REAL(dp) :: z, kx, ky

    zloop: DO iz = izs,ize
     z = zarray(iz)

    ! metric
       gxx(iz) = 1._dp
       gxy(iz) = 0._dp
       gyy(iz) = 1._dp

    ! Relative strengh of radius
       hatR(iz) = 1._dp

    ! Jacobian
       Jacobian(iz) = 1._dp

    ! Relative strengh of modulus of B
       hatB(iz) = 1._dp

    ! Curvature operator
       DO iky = ikys, ikye
         ky = kyarray(iky)
          DO ikx= ikxs, ikxe
            kx = kxarray(ikx)
            Ckxky(ikx, iky, iz) = -kx ! .. multiply by hatB to cancel the 1/ hatB factor in moments_eqs_rhs.f90 routine
          ENDDO
       ENDDO
    ! coefficient in the front of parallel derivative
       gradz_coeff(iz) = 1._dp
    ENDDO zloop

   END SUBROUTINE eval_1D_geometry
   !
   !--------------------------------------------------------------------------------
   !
end module geometry
