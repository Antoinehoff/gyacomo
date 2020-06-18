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
  COMPLEX(dp) :: kr, kz, kperp2
  COMPLEX(dp) :: xkm_, xk2_ ! sub kernel factor for recursive build
  COMPLEX(dp) :: b_e2, b_i2, alphaD
  COMPLEX(dp) :: sigmaetaue2, sigmaitaui2 ! To avoid redondant computation
  COMPLEX(dp) :: gammaD, gammaphi
  COMPLEX(dp) :: sum_kernel2_e,    sum_kernel2_i
  COMPLEX(dp) :: sum_kernel_mom_e, sum_kernel_mom_i

  sigmaetaue2 = sigma_e**2 * tau_e/2.
  sigmaitaui2 = sigma_i**2 * tau_i/2.

  DO ikr=ikrs,ikre
    DO ikz=ikzs,ikze
      kr    = krarray(ikr)
      kz    = kzarray(ikz)
      kperp2 = kr**2 + kz**2

      ! Compute electrons sum(Kernel * Ne0n)  (skm) and sum(Kernel**2) (sk2)
      b_e2  =  kperp2 * sigmaetaue2 ! non dim kernel argument (kperp2 sigma_a sqrt(2 tau_a)/2)

      xkm_ = 1. ! Initialization for n = 0
      xk2_ = 1.
      
      sum_kernel_mom_e = moments_e(1,1,ikr,ikz,updatetlevel)
      sum_kernel2_e    = xk2_

      if (jmaxe .GT. 0) then
        DO ine=2,jmaxe+1
          ine_dp = REAL(ine-1,dp)

          xkm_ = xkm_ * b_e2/ine_dp
          xk2_ = xk2_ *(b_e2/ine_dp)**2

          sum_kernel_mom_e  = sum_kernel_mom_e  + xkm_ * moments_e(1,ine,ikr,ikz,updatetlevel)
          sum_kernel2_e     = sum_kernel2_e     + xk2_
        END DO
      endif
      sum_kernel2_e    = sum_kernel2_e    * exp(-b_e2)
      sum_kernel_mom_e = sum_kernel_mom_e * exp(-2.*b_e2)

      ! Compute ions sum(Kernel * Ni0n)  (skm) and sum(Kernel**2) (sk2)
      b_i2  = kperp2 * sigmaitaui2

      xkm_ = 1.
      xk2_ = 1.

      sum_kernel_mom_i = moments_i(1, 1, ikr, ikz, updatetlevel)
      sum_kernel2_i    = xk2_
      if (jmaxi .GT. 0) then
        DO ini=2,jmaxi + 1
          ini_dp = REAL(ini-1,dp) ! Real index (0 to jmax)

          xkm_ = xkm_ * b_i2/ini_dp
          xk2_ = xk2_ *(b_i2/ini_dp)**2

          sum_kernel_mom_i  = sum_kernel_mom_i  + xkm_ * moments_i(1,ini,ikr,ikz,updatetlevel)
          sum_kernel2_i     = sum_kernel2_i     + xk2_
        END DO
      endif
      sum_kernel2_i    = sum_kernel2_i    * exp(-b_i2)
      sum_kernel_mom_i = sum_kernel_mom_i * exp(-2.*b_i2)

      ! Assembling the poisson equation
      alphaD   = kperp2 * lambdaD**2.
      gammaD   = alphaD + (q_e**2)/tau_e * sum_kernel2_e &
                        + (q_i**2)/tau_i * sum_kernel2_i

      gammaphi = q_e * sum_kernel_mom_e + q_i * sum_kernel_mom_i
      
      phi(ikr, ikz) =  gammaphi/gammaD
    
    END DO
  END DO

END SUBROUTINE poisson