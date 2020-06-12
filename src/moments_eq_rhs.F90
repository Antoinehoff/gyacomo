SUBROUTINE moments_eq_rhs  

  USE basic
  USE time_integration
  USE array
  USE fields
  USE fourier_grid
  USE model

  use prec_const
  IMPLICIT NONE

  INTEGER:: ip,ij, ipj, ikr,ikz
  REAL(dp) :: ip_dp, ij_dp
  REAL(dp) :: tmp, moml, sqrtT
  
  do ipj = ipjs, ipje

    if (ipj .le. Nmome + 1) then ! electrons moments
      CALL rabe(ipj, ip, ij) ! Compute p,j electrons moments degree
    else ! ions moments
      CALL rabi(ipj, ip, ij) ! Compute p,j ions moments degree
    endif

    ip_dp = real(ip,dp) 
    ij_dp = real(ij,dp) 

    do ikr = ikrs,ikre 
      do ikz = ikzs,ikze

      ! N_e^{p} term

      ! N_e^{p+1} term


      ! N_e^{p-1} term


      ! N_e^{p-2} term

    end do
  end do
  
    do ikr = ikrs,ikre
      do ikz = ikzs,ikze
        moments_rhs(ipj,ikr,ikz,updatetlevel) = 0._dp ! debug
      end do
    end do
  end do

END SUBROUTINE moments_eq_rhs
