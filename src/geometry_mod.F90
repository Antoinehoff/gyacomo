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
  public
  COMPLEX(dp), PROTECTED :: iInt_Jacobian

contains

  subroutine eval_magnetic_geometry
    ! evalute metrix, elementwo_third_kpmaxts, jacobian and gradient
    implicit none
    REAL(dp) :: kx,ky
    COMPLEX(dp), DIMENSION(izs:ize) :: integrant
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
    DO eo = 0,1
    DO iz = izs,ize
       DO iky = ikys, ikye
         ky = kyarray(iky)
          DO ikx = ikxs, ikxe
            kx = kxarray(ikx)
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
    integrant = Jacobian(:,0) ! Convert into complex array
    CALL simpson_rule_z(integrant,iInt_Jacobian)
    iInt_Jacobian = 1._dp/iInt_Jacobian ! reverse it
  END SUBROUTINE eval_magnetic_geometry
  !
  !--------------------------------------------------------------------------------
  !

  subroutine eval_salphaB_geometry
  ! evaluate s-alpha geometry model
  implicit none
  REAL(dp) :: z, kx, ky

  parity: DO eo = 0,1
  zloop: DO iz = izs,ize
    z = zarray(iz,eo)

    ! metric
      gxx(iz,eo) = 1._dp
      gxy(iz,eo) = shear*z
      gyy(iz,eo) = 1._dp + (shear*z)**2

    ! Relative strengh of radius
      hatR(iz,eo) = 1._dp + eps*COS(z)

    ! Jacobian
      Jacobian(iz,eo) = q0*hatR(iz,eo)

    ! Relative strengh of modulus of B
      hatB(iz,eo) = 1._dp / hatR(iz,eo)

    ! Derivative of the magnetic field strenght
      gradxB(iz,eo) = -COS(z)
      gradzB(iz,eo) = eps * SIN(z) / hatR(iz,eo)

    ! Curvature operator
      DO iky = ikys, ikye
        ky = kyarray(iky)
         DO ikx= ikxs, ikxe
           kx = kxarray(ikx)
           Ckxky(ikx, iky, iz,eo) = (-SIN(z)*kx - (COS(z) + shear* z* SIN(z))*ky) * hatB(iz,eo) ! .. multiply by hatB to cancel the 1/ hatB factor in moments_eqs_rhs.f90 routine
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
end module geometry
