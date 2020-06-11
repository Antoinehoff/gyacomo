SUBROUTINE memory
  ! Allocate arrays (done dynamically otherwise size is unknown)

  USE array
  USE basic
  USE fields
  USE model, ONLY: gradient_scheme
  USE space_grid
  USE time_integration  

  use prec_const
  IMPLICIT NONE



  ! Principal fields
  CALL allocate_array(theta,izs,ize,1,ntimelevel)
  CALL allocate_array(temp,izs,ize,1,ntimelevel)
  CALL allocate_array(vpar,izs,ize,1,ntimelevel)
  CALL allocate_array(moments,ips,ipe,izs,ize,1,ntimelevel)
  CALL allocate_array(phi,izs,ize)
  ! On opposite grid
  CALL allocate_array(vpar_n,izs,ize)

  ! Right hand sides of the time evolution equations
  CALL allocate_array(theta_rhs,izs,ize,1,ntimelevel)
  CALL allocate_array(temp_rhs,izs,ize,1,ntimelevel)
  CALL allocate_array(vpar_rhs,izs,ize,1,ntimelevel)
  CALL allocate_array(moments_rhs,ips,ipe,izs,ize,1,ntimelevel)
  ! On opposite grid
  CALL allocate_array(temp_rhs_v,izs,ize)
  CALL allocate_array(vpar_rhs_n,izs,ize)
  
  ! Auxiliary quantities
  CALL allocate_array(sqrt_exp_temp,izs,ize)
  CALL allocate_array(sqrt_exp_temp_v,izs,ize)

  ! Spatial 1rst derivatives on respective grids
  CALL allocate_array(thetaz,izs,ize)
  CALL allocate_array(sqrt_exp_tempz,izs,ize)
  CALL allocate_array(vparz,izs,ize)
  CALL allocate_array(phiz,izs,ize)
  CALL allocate_array(momentsz,ips,ipe,izs,ize)
  ! On opposite grid
  CALL allocate_array(thetaz_v,izs,ize)
  CALL allocate_array(sqrt_exp_tempz_v,izs,ize)
  CALL allocate_array(vparz_n,izs,ize)
  CALL allocate_array(phiz_v,izs,ize)
  

  ! select case (gradient_scheme)
  ! case('up2') ! If upwind scheme : need spatial derivaties for the "minus side"
  !   CALL allocate_array(thetazm,izs,ize)
  !   CALL allocate_array(tempzm,izs,ize)
  !   CALL allocate_array(vparzm,izs,ize)
  !   CALL allocate_array(phizm,izs,ize)
  !   CALL allocate_array(momentszm,ips,ipe,izs,ize)

  !   CALL allocate_array(thetazz,0,0)
  !   CALL allocate_array(tempzz,0,0)
  !   CALL allocate_array(vparzz,0,0)
  !   CALL allocate_array(momentszz,0,0,0,0)
  
  ! case('fa2','fd4') ! If no upwind scheme : numerical diffusion needs 2nd order spatial derivatives
  !   CALL allocate_array(thetazm,0,0)
  !   CALL allocate_array(tempzm,0,0)
  !   CALL allocate_array(vparzm,0,0)
  !   CALL allocate_array(phizm,0,0)
  !   CALL allocate_array(momentszm,0,0,0,0)
  
    CALL allocate_array(thetazz,izs,ize)
    CALL allocate_array(tempzz,izs,ize)
    CALL allocate_array(vparzz,izs,ize)
    CALL allocate_array(momentszz,ips,ipe,izs,ize)
  ! end select


  ! Intermediate steps in rhs of equations
    ! theta
  CALL allocate_array(theta_Cpl,izs,ize)
  CALL allocate_array(theta_Hpl,izs,ize)
    ! temp
  CALL allocate_array(temp_Cpl,izs,ize)
  CALL allocate_array(temp_Fpl,izs,ize)
  CALL allocate_array(temp_Hpl,izs,ize)
  CALL allocate_array(temp_Ipl,izs,ize)
    ! vpar
  CALL allocate_array(vpar_Cpl,izs,ize)
  CALL allocate_array(vpar_Dpl,izs,ize)
  CALL allocate_array(vpar_Epl,izs,ize)
  CALL allocate_array(vpar_Fpl,izs,ize)
  CALL allocate_array(vpar_Hpl,izs,ize)
    ! moments
  CALL allocate_array(moments_Apl,izs,ize)
  CALL allocate_array(moments_Bpl,izs,ize)
  CALL allocate_array(moments_Cpl,izs,ize)
  CALL allocate_array(moments_Dpl,izs,ize)
  CALL allocate_array(moments_Epl,izs,ize)
  CALL allocate_array(moments_Fpl,izs,ize)
  CALL allocate_array(moments_Gpl,izs,ize)
  CALL allocate_array(moments_Hpl,izs,ize)


  ! Poisson solver
  ! CALL init(neq, 1, laplace_smumps)
  CALL allocate_array(laplace_rhs,1,neq)  
  CALL allocate_array(laplace_sol,1,neq)  


  ! FFT
  CALL allocate_array(fft_filter,izs,ize)
  CALL allocate_array(fft_input,izs,ize)


  ! WENO gradients
  CALL allocate_array(tmp1,izs,ize)
  CALL allocate_array(tmp2,izs,ize)
  CALL allocate_array(omegak_times_q3k,izs,ize,0,2)


END SUBROUTINE memory
