SUBROUTINE theta_eq_rhs  

  ! This can be used to evaluate the RHS of the Hermite Hierarchy Equation .... 

 ! theta is ln(density)
  USE time_integration
  USE array
  USE fields
  USE space_grid
  USE model


  use prec_const
  IMPLICIT NONE

  INTEGER :: iz
  real(dp) :: tmp


  do iz = izs,ize ! Compute coefficients
    theta_Cpl(iz) = -sqrt_exp_temp(iz)*INVSQRT2*vpar_n(iz)
    theta_Hpl(iz) = -sqrt_exp_temp(iz)*SQRT2
  end do


  select case (gradient_scheme)
  ! case('up2') ! Upwind scheme, check the sign of coefficients
  !   do iz = izs,ize
  !     if (theta_Cpl(iz)>0) then
  !       tmp = theta_Cpl(iz)*thetaz(iz)
  !     else
  !       tmp = theta_Cpl(iz)*thetazm(iz)
  !     endif

  !     theta_rhs(iz, updatetlevel) = tmp &
  !                                   +theta_Hpl(iz)*vparzm(iz) ! theta_Hpl is always negative
      
  !   end do

  case('fa2','fd4','we4')
    do iz = izs,ize
      theta_rhs(iz, updatetlevel) = diff_theta*thetazz(iz) & ! numerical diffusion added for stability
                                    +theta_Cpl(iz)*thetaz(iz) &
                                    +theta_Hpl(iz)*vparz_n(iz)
!theta_rhs(iz, updatetlevel) = 0._dp ! debug !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end do

  end select

END SUBROUTINE theta_eq_rhs
