SUBROUTINE temp_eq_rhs  

  USE basic
  USE time_integration
  USE array
  USE fields
  USE space_grid
  USE model
  USE gradients

  use prec_const
  IMPLICIT NONE

  INTEGER :: iz
  real(dp) :: tmp


  do iz = izs,ize ! Compute coefficients
    temp_Cpl(iz) = -sqrt_exp_temp(iz)*SQRT3*INVSQRT2*moments(3,iz,updatetlevel)
    temp_Fpl(iz) = -2._dp*SQRT2*(SQRT3*moments(3,iz,updatetlevel)+2._dp)
!    temp_Fpl(iz) = -sqrt_exp_temp(iz)*SQRT2*(SQRT3*moments(3,iz,updatetlevel)+2._dp)
!    temp_Fpl(iz) = -SQRT2*(SQRT3*moments(3,iz,updatetlevel)+2._dp)
    temp_Hpl(iz) = -sqrt_exp_temp(iz)*2._dp*SQRT2
    temp_Ipl(iz) = -sqrt_exp_temp(iz)*SQRT2*SQRT3
  end do

  

  select case (gradient_scheme)
  ! case('up2') ! Upwind scheme, check the sign of coefficients
  !   do iz = izs,ize
  !     if (temp_Cpl(iz)>0) then
  !       tmp = temp_Cpl(iz)*thetaz(iz)
  !     else
  !       tmp = temp_Cpl(iz)*thetazm(iz)
  !     endif

  !     if (temp_Fpl(iz)>0) then
  !       tmp = tmp + temp_Fpl(iz)*tempz(iz)
  !     else
  !       tmp = tmp + temp_Fpl(iz)*tempzm(iz)
  !     endif

  !     temp_rhs(iz,updatetlevel) = tmp &
  !                                 +temp_Hpl(iz)*vparzm(iz) &    ! temp_Hpl is always negative
  !                                 +temp_Ipl(iz)*momentszm(3,iz) ! temp_Ipl is always negative
  !   end do

  case('fa2','fd4','we4')
    do iz = izs,ize
      temp_rhs(iz, updatetlevel) = diff_temp*tempzz(iz) & ! numerical diffusion added for stability
                                   +temp_Cpl(iz)*thetaz(iz) &
                                   +temp_Fpl(iz)*sqrt_exp_tempz(iz) &
                                   +temp_Hpl(iz)*vparz_n(iz) &
                                   +temp_Ipl(iz)*momentsz(3,iz) 
!                                   - SQRT2*2._dp*nu*temp(iz,updatetlevel) ! Lenard-Bernstein Operator conserves energy
    end do

  end select

  ! Compute interpolation on v-grid, used by vpar_eq_rhs
  CALL interp_n2v(temp_rhs(:, updatetlevel), temp_rhs_v)

END SUBROUTINE temp_eq_rhs
