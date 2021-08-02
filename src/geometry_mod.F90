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
   REAL(dp) :: shear   = 0._dp
   REAL(dp) :: epsilon = 0.18_dp

   ! Metric elements

   DO iz = izs,ize
      gxx(iz) = 1._dp
      gxy(iz) = shear*zarray(iz) - epsilon*sin(zarray(iz))
      gyy(iz) = 1._dp + (shear*zarray(iz))**2 &
           - 2._dp * epsilon *COS(zarray(iz)) - 2._dp*shear * epsilon * zarray(iz)*SIN(zarray(iz))
      gxz( iz) = - SIN(zarray(iz))
      gyz( iz) = ( 1._dp - 2._dp * epsilon *COS(zarray( iz)) - epsilon*shear * zarray( iz) * SIN(zarray(iz)) ) /epsilon
   ENDDO

   ! Relative strengh of radius
   DO iz = izs,ize
       hatR(iz) = 1._dp + epsilon*cos(zarray(iz))
    ENDDO

    DO iz = izs, ize
       hatB( iz) = sqrt(gxx(iz) * gyy(iz) - gxy(iz)*  gxy(iz))
   ENDDO


   ! Jacobian
    DO iz = izs,ize
       Jacobian(iz) = q0*hatR(iz)*hatR(iz)
    ENDDO

    ! Derivative of the magnetic field strenght
    DO iz = izs, ize
       gradxB(iz) = - ( COS( zarray(iz))   + epsilon* SIN( zarray(iz)) * SIN(zarray(iz) )) / hatB(iz)/hatB(iz)
       gradzB( iz) = epsilon * SIN(zarray(iz)) *( 1._dp - epsilon*COS(zarray(iz)) ) / hatB(iz)/hatB(iz)
    ENDDO

    ! Gemoetrical coefficients for the curvature operator
    ! Note: Gamma2 and Gamma3 are obtained directly form Gamma1 in the expression of the curvature operator implemented here
   DO iz = izs, ize
      Gamma1( iz) = gxy(iz) * gxy(iz) - gxx(iz) * gyy(iz)
      Gamma2( iz) = gxz( iz) * gxy( iz) - gxx( iz) * gyz( iz)
      Gamma3( iz) = gxz( iz) * gyy( iz) - gxy(iz) * gyz(  iz)
   ENDDO

    ! Curvature operator
    DO iz = izs, ize
       DO iky = ikys, ikye
          DO ikx= ikxs, ikxe
             Ckxky( ikx, iky,  iz) = Gamma1(iz)/hatB(iz)*((kxarray(ikx) &
                  + shear*zarray(iz)*kyarray(iky)) * gradzB(iz)/epsilon &
                  - gradxB(iz)* kyarray(iky))
          ENDDO
       ENDDO
    ENDDO

    ! coefficient in the front of parallel derivative
    DO iz = izs, ize
       gradz_coeff( iz) = 1._dp / Jacobian( iz) / hatB( iz)
    ENDDO

  END SUBROUTINE eval_circular_geometry

end module geometry
