SUBROUTINE evaluation_auxfield_total

  USE space_grid
  USE fields
  USE array
  USE time_integration
  USE gradients
  USE prec_const

  IMPLICIT NONE
  integer :: iz
  integer :: ip

  DO iz = izs, ize
    ! Actual temperature = exp ( ln temp )
    ! Actual density = exp ( ln theta )
    sqrt_exp_temp(iz) = exp(temp(iz,updatetlevel)*0.5_dp)
  END DO

  CALL parallel_gradients

  CALL interpolation

END SUBROUTINE evaluation_auxfield_total


!Compute parallel gradients for RHS of equations
SUBROUTINE parallel_gradients
  USE array
  USE fields
  USE gradients
  USE model, only: gradient_scheme
  USE space_grid
  USE time_integration, only: updatetlevel

  use prec_const
  IMPLICIT NONE

  integer :: ip

  CALL gradz(theta(:,updatetlevel), thetaz)
  CALL gradz(sqrt_exp_temp, sqrt_exp_tempz)
  CALL gradz(vpar(:,updatetlevel), vparz)
  CALL gradz(phi(:), phiz)
  do ip=ips,ipe
    CALL gradz(moments(ip,:,updatetlevel), momentsz(ip,:))
  enddo

  ! Derivatives that are on the opposite grid
  CALL gradz_n2v(theta(:,updatetlevel), thetaz_v)
  CALL gradz_n2v(sqrt_exp_temp, sqrt_exp_tempz_v)
  CALL gradz_v2n(vpar(:,updatetlevel), vparz_n)
  CALL gradz_n2v(phi, phiz_v)

  select case (gradient_scheme)
!  case('up2') ! If upwind scheme : need spatial derivaties for the "minus side"
!    CALL gradzm(theta(:,updatetlevel), thetazm)
!    CALL gradzm(temp(:,updatetlevel), tempzm)
!    CALL gradzm(vpar(:,updatetlevel), vparzm)
!    CALL gradzm(phi(:), phizm)
!    do ip=ips,ipe
!      CALL gradzm(moments(ip,:,updatetlevel), momentszm(ip,:))
!    enddo

  case('fa2','fd4','we4') ! If no upwind scheme : numerical diffusion needs 2nd order spatial derivatives
    CALL gradzz(theta(:,updatetlevel), thetazz)
    CALL gradzz(temp(:,updatetlevel), tempzz)
    CALL gradzz(vpar(:,updatetlevel), vparzz)    
    do ip=ips,ipe
      CALL gradzz(moments(ip,:,updatetlevel), momentszz(ip,:))
    enddo
  end select


END SUBROUTINE parallel_gradients



SUBROUTINE interpolation
  USE array
  USE fields
  USE gradients
  USE time_integration, only: updatetlevel

  use prec_const
  IMPLICIT NONE


  CALL interp_n2v(sqrt_exp_temp, sqrt_exp_temp_v)
  
  CALL interp_v2n(vpar(:,updatetlevel), vpar_n)



END SUBROUTINE interpolation

