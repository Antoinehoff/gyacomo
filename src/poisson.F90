SUBROUTINE poisson
  ! Solve poisson equation to get phi

  USE time_integration, ONLY: updatetlevel
  USE array
  USE fields
  USE fourier_grid
  use model, ONLY : tau_e, tau_i, sigma_e, sigma_i, q_e, q_i, lambdaD

  USE prec_const
  IMPLICIT NONE

  INTEGER     :: ikr,ikz, ini,ine, i0j
  REAL(dp)    :: ini_dp, ine_dp
  REAL(dp)    :: kr, kz, kperp2
  REAL(dp)    :: Kn, Kn2 ! sub kernel factor for recursive build
  REAL(dp)    :: b_e2, b_i2, alphaD
  REAL(dp)    :: sigmae2_taue_o2, sigmai2_taui_o2, qe2_taue, qi2_taui ! To avoid redondant computation
  COMPLEX(dp) :: gammaD, gammaphi
  COMPLEX(dp) :: sum_kernel2_e,    sum_kernel2_i    ! Store to compute sum Kn^2
  COMPLEX(dp) :: sum_kernel_mom_e, sum_kernel_mom_i ! Store to compute sum Kn*Napn
  
  !Precompute species dependant factors
  sigmae2_taue_o2 = sigma_e**2 * tau_e/2. 
  sigmai2_taui_o2 = sigma_i**2 * tau_i/2.
  qe2_taue        = (q_e**2)/tau_e
  qi2_taui        = (q_i**2)/tau_i

  DO ikr=ikrs,ikre
    DO ikz=ikzs,ikze
      kr    = krarray(ikr)
      kz    = kzarray(ikz)
      kperp2 = kr**2 + kz**2

      !! Electrons sum(Kernel * Ne0n)  (skm) and sum(Kernel**2) (sk2)
      b_e2  =  kperp2 * sigmae2_taue_o2 ! non dim kernel argument (kperp2 sigma_a sqrt(2 tau_a))

      ! Initialization for n = 0 (ine = 1)
      Kn  = 1. 
      Kn2 = 1.
      sum_kernel_mom_e = moments_e(1,1,ikr,ikz,updatetlevel)
      sum_kernel2_e    = Kn2

      ! loop over n only if the max polynomial degree is 1 or more
      if (jmaxe .GT. 0) then
        DO ine=1,jmaxe+1 ! ine = n+1
          ine_dp = REAL(ine-1,dp)

          Kn  = Kn  * b_e2/(ine_dp+1)       ! We update iteratively the kernel functions (to spare factorial computations)
          Kn2 = Kn2 *(b_e2/(ine_dp+1))**2   ! the exp term is still missing but does not depend on ine ...

          sum_kernel_mom_e  = sum_kernel_mom_e  + Kn * moments_e(1,ine,ikr,ikz,updatetlevel)
          sum_kernel2_e     = sum_kernel2_e     + Kn2 ! ... sum recursively ...
        END DO
      endif
      sum_kernel_mom_e = sum_kernel_mom_e * exp(-b_e2) ! ... multiply the final sum with the missing exponential term
      sum_kernel2_e    = sum_kernel2_e    * exp(-2.*b_e2)    ! and its squared using exp(x)^2 = exp(2x).

      !! Ions sum(Kernel * Ni0n)  (skm) and sum(Kernel**2) (sk2)
      b_i2  = kperp2 * sigmai2_taui_o2

      ! Initialization for n = 0 (ini = 1)
      Kn  = 1.
      Kn2 = 1.
      sum_kernel_mom_i = moments_i(1, 1, ikr, ikz, updatetlevel)
      sum_kernel2_i    = Kn2

      ! loop over n only if the max polynomial degree is 1 or more
      if (jmaxi .GT. 0) then
        DO ini=1,jmaxi+1
          ini_dp = REAL(ini-1,dp) ! Real index (0 to jmax)

          Kn  = Kn * b_i2/(ini_dp+1)       ! We update iteratively the kernel functions this time for ions ...
          Kn2 = Kn2 *(b_i2/(ini_dp+1.))**2

          sum_kernel_mom_i  = sum_kernel_mom_i  + Kn * moments_i(1,ini,ikr,ikz,updatetlevel)
          sum_kernel2_i     = sum_kernel2_i     + Kn2 ! ... sum recursively ...
        END DO
      endif
      sum_kernel_mom_i = sum_kernel_mom_i * exp(-b_i2) ! ... multiply the final sum with the missing exponential term.
      sum_kernel2_i    = sum_kernel2_i    * exp(-2.*b_i2) 

      !!! Assembling the poisson equation
      alphaD   = kperp2 * lambdaD**2.
      gammaD   = alphaD + qe2_taue * (1. - sum_kernel2_e) &
                        + qi2_taui * (1. - sum_kernel2_i)

      gammaphi = q_e * sum_kernel_mom_e + q_i * sum_kernel_mom_i ! gamma_D * phi term
      
      phi(ikr, ikz) =  gammaphi/gammaD
    
    END DO
  END DO

END SUBROUTINE poisson