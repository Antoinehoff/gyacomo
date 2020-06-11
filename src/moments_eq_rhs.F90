SUBROUTINE moments_eq_rhs  

  USE basic
  USE time_integration
  USE array
  USE fields
  USE space_grid
  USE model

  use prec_const
  IMPLICIT NONE

  INTEGER:: ip, iz
  REAL(dp) :: ip_dp
  REAL(dp) :: tmp, moml, sqrtT
  
  do ip = ips, ipe
    ip_dp = real(ip,dp) ! From int to double (compute once)


    do iz = izs,ize ! Compute coefficients

      sqrtT = sqrt_exp_temp(iz)

      ! N_e^{p} term
      moml = moments(ip,iz,updatetlevel)      
      moments_Apl(iz) = -nu*ip_dp*moml   ! Collisional damping from Lenard-Bernstein Operator
      moments_Bpl(iz) = -moml
      moments_Cpl(iz) = -sqrtT*vpar_n(iz)*INVSQRT2*moml
      moments_Dpl(iz) = 0._dp
      moments_Epl(iz) = -ip_dp*0.5_dp*moml
      moments_Fpl(iz) = -vpar_n(iz)*2._dp*SQRT2*ip_dp*moml
!      moments_Fpl(iz) = -sqrtT*vpar_n(iz)*SQRT2*ip_dp*moml
      moments_Gpl(iz) = 0._dp
      moments_Hpl(iz) = -sqrtT*SQRT2*(ip_dp+1._dp)*moml

      ! N_e^{p+1} term
      if (ip+1 .le. ipe) then
        moml = moments(ip+1,iz,updatetlevel)
        moments_Cpl(iz) = moments_Cpl(iz) -sqrtT*sqrt(ip_dp+1._dp)*0.5_dp*moml
        moments_Fpl(iz) = moments_Fpl(iz) -ip_dp*sqrt(ip_dp+1._dp)*moml
!        moments_Fpl(iz) = moments_Fpl(iz) -sqrtT*ip_dp*sqrt(ip_dp+1._dp)*0.5_dp*moml
      endif


      ! N_e^{p-1} term
      if (ip-1 .ge. ips) then
         moml = moments(ip-1,iz,updatetlevel)
         if (ip .ne. 3) then 
            moments_Apl(iz) = moments_Apl(iz) -nu*SQRT2*sqrt(ip_dp)*moml*vpar_n(iz) ! Collisions
         endif
         moments_Cpl(iz) = moments_Cpl(iz) -sqrtT*sqrt(ip_dp)*0.5_dp*moml
         moments_Dpl(iz) = moments_Dpl(iz) +sqrt(ip_dp)/sqrtT*moml
         moments_Epl(iz) = moments_Epl(iz) -vpar_n(iz)*sqrt(ip_dp*0.5_dp)*moml
         moments_Fpl(iz) = moments_Fpl(iz) -sqrt(ip_dp)*(2._dp*ip_dp-1._dp+2._dp*vpar_n(iz)*vpar_n(iz))*moml
!         moments_Fpl(iz) = moments_Fpl(iz) -sqrtT*sqrt(ip_dp)*0.5_dp*(2._dp*ip_dp-1._dp+2._dp*vpar_n(iz)*vpar_n(iz))*moml
         moments_Gpl(iz) = moments_Gpl(iz) -SQRT2*sqrt(ip_dp)*moml
         moments_Hpl(iz) = moments_Hpl(iz) -sqrtT*vpar_n(iz)*2._dp*sqrt(ip_dp)*moml
      endif


      ! N_e^{p-2} term
      if (ip-2 .ge. ips) then
        moml = moments(ip-2,iz,updatetlevel)
        moments_Epl(iz) = moments_Epl(iz) -sqrt(ip_dp*(ip_dp-1._dp))*0.5_dp*moml
        moments_Fpl(iz) = moments_Fpl(iz) -2._dp*sqrt(2._dp*ip_dp*(ip_dp-1._dp))*moml
!        moments_Fpl(iz) = moments_Fpl(iz) -sqrtT*sqrt(2._dp*ip_dp*(ip_dp-1._dp))*moml
        moments_Hpl(iz) = moments_Hpl(iz) -sqrtT*sqrt(2._dp*ip_dp*(ip_dp-1._dp))*moml
      endif


      ! N_e^{p-3} term
      if (ip-3 .ge. ips) then
        moments_Fpl(iz) = moments_Fpl(iz) -sqrt(ip_dp*(ip_dp-1._dp)*(ip_dp-2._dp))*moments(ip-3,iz,updatetlevel)
!        moments_Fpl(iz) = moments_Fpl(iz) -sqrtT*sqrt(ip_dp*(ip_dp-1._dp)*(ip_dp-2._dp))*0.5_dp*moments(ip-3,iz,updatetlevel)
      else if (ip .eq. ips) then ! ip-3=0 so N_e^0 is equal to 1
        moments_Fpl(iz) = moments_Fpl(iz) -SQRT3*SQRT2
!        moments_Fpl(iz) = moments_Fpl(iz) -sqrtT*SQRT3*SQRT2*0.5_dp
      endif

    end do

    select case (gradient_scheme)
    ! case('up2') ! Upwind scheme, check the sign of coefficients
    !   do iz = izs,ize
    !     if (moments_Cpl(iz)>0) then
    !       tmp = moments_Cpl(iz)*thetaz(iz)
    !     else
    !       tmp = moments_Cpl(iz)*thetazm(iz)
    !     endif

    !     if (moments_Dpl(iz)>0) then
    !       tmp = tmp + moments_Dpl(iz)*phiz(iz)
    !     else
    !       tmp = tmp + moments_Dpl(iz)*phizm(iz)
    !     endif

    !     if (moments_Fpl(iz)>0) then
    !       tmp = tmp + moments_Fpl(iz)*tempz(iz)
    !     else
    !       tmp = tmp + moments_Fpl(iz)*tempzm(iz)
    !     endif
        
    !     if (moments_Hpl(iz)>0) then
    !       tmp = tmp + moments_Hpl(iz)*vparz(iz)
    !     else
    !       tmp = tmp + moments_Hpl(iz)*vparzm(iz)
    !     endif

    !     if (vpar(iz,updatetlevel)<0) then ! I_p^l
    !       moments_Ipl = -sqrtT*SQRT2*vpar(iz,updatetlevel)*momentsz(ip,iz)
    !     else
    !       moments_Ipl = -sqrtT*SQRT2*vpar(iz,updatetlevel)*momentszm(ip,iz)
    !     endif
    !     if (ip+1 .le. ipe) moments_Ipl = moments_Ipl -sqrtT*sqrt(ip_dp+1._dp)*momentszm(ip+1,iz) ! coefficient always negative
    !     if (ip-1 .ge. ips) moments_Ipl = moments_Ipl -sqrtT*sqrt(ip_dp)*momentszm(ip-1,iz) ! coefficient always negative


    !     moments_rhs(ip,iz,updatetlevel) = tmp &
    !                                      +moments_Apl(iz) &
    !                                      +moments_Bpl(iz)*theta_rhs(iz,updatetlevel) &
    !                                      +moments_Epl(iz)*temp_rhs(iz,updatetlevel) &
    !                                      +moments_Gpl(iz)*vpar_rhs(iz,updatetlevel) &
    !                                      +moments_Ipl

    !   end do
    case('fa2','fd4','we4')
      do iz = izs,ize

        moments_Ipl = -sqrtT*SQRT2*vpar_n(iz)*momentsz(ip,iz)
        if (ip+1 .le. ipe) moments_Ipl = moments_Ipl -sqrtT*sqrt(ip_dp+1._dp)*momentsz(ip+1,iz) ! coefficient always negative
        if (ip-1 .ge. ips) moments_Ipl = moments_Ipl -sqrtT*sqrt(ip_dp)*momentsz(ip-1,iz) ! coefficient always negative

        moments_rhs(ip,iz,updatetlevel) = diff_moments*momentszz(ip,iz) & ! numerical diffusion added for stability
                                         +moments_Apl(iz) &
                                         +moments_Bpl(iz)*theta_rhs(iz,updatetlevel) &
                                         +moments_Cpl(iz)*thetaz(iz) &
                                         +moments_Dpl(iz)*phiz(iz) &
                                         +moments_Epl(iz)*temp_rhs(iz,updatetlevel) &
                                         +moments_Fpl(iz)*sqrt_exp_tempz(iz) &
                                         +moments_Gpl(iz)*vpar_rhs_n(iz) &
                                         +moments_Hpl(iz)*vparz_n(iz) &
                                         +moments_Ipl
!        moments_rhs(ip,iz,updatetlevel) = 0._dp ! debug
      end do
    end select

  end do

END SUBROUTINE moments_eq_rhs
