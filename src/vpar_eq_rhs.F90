SUBROUTINE vpar_eq_rhs  

  USE time_integration
  USE array
  USE fields
  USE space_grid
  USE gradients
  USE model, ONLY: diff_vpar, gradient_scheme,nu

  use prec_const
  IMPLICIT NONE

  INTEGER :: iz
  real(dp) :: tmp


  do iz = izs,ize ! Compute coefficients
    vpar_Cpl(iz) = -sqrt_exp_temp_v(iz)*0.5_dp*INVSQRT2
    vpar_Dpl(iz) = INVSQRT2/sqrt_exp_temp_v(iz)   ! ... should be +, otherwise instability
    vpar_Epl(iz) = -vpar(iz,updatetlevel)*0.5_dp
    vpar_Fpl(iz) = -(1._dp+2._dp*vpar(iz,updatetlevel)*vpar(iz,updatetlevel))*INVSQRT2
!    vpar_Fpl(iz) = -sqrt_exp_temp_v(iz)*(1._dp+2._dp*vpar(iz,updatetlevel)*vpar(iz,updatetlevel))*0.5_dp*INVSQRT2
    vpar_Hpl(iz) = -sqrt_exp_temp_v(iz)*SQRT2*vpar(iz,updatetlevel)
  end do



  select case (gradient_scheme)
  ! case('up2') ! Upwind scheme, check the sign of coefficients
  !   do iz = izs,ize
  !     if (vpar_Hpl(iz)>0) then
  !       tmp = vpar_Hpl(iz)*vparz(iz)
  !     else
  !       tmp = vpar_Hpl(iz)*vparzm(iz)
  !     endif

  !     vpar_rhs(iz,updatetlevel) = tmp &
  !                                 +vpar_Cpl(iz)*thetazm(iz) &               ! vpar_Cpl is always negative
  !                                 +vpar_Dpl(iz)*phizm(iz) &                 ! vpar_Dpl is always negative
  !                                 +vpar_Epl(iz)*temp_rhs(iz,updatetlevel) &
  !                                 +vpar_Fpl(iz)*tempzm(iz)                  ! vpar_Fpl is always negative
  !   end do

  case('fa2','fd4','we4')
    do iz = izs,ize
      vpar_rhs(iz, updatetlevel) = diff_vpar*vparzz(iz) & ! numerical diffusion added for stability
                                   +vpar_Cpl(iz)*thetaz_v(iz) &
                                   +vpar_Dpl(iz)*phiz_v(iz) &
                                   +vpar_Epl(iz)*temp_rhs_v(iz) &
                                   +vpar_Fpl(iz)*sqrt_exp_tempz_v(iz) &
                                   +vpar_Hpl(iz)*vparz(iz) &
                                   -SQRT2*nu*vpar(iz,updatetlevel) ! add collisional friction from Lenard Bernstein Collision Operator
    end do

  end select

  ! Compute interpolation on n-grid, used by moments_eq_rhs
  CALL interp_v2n(vpar_rhs(:, updatetlevel), vpar_rhs_n)


END SUBROUTINE vpar_eq_rhs
