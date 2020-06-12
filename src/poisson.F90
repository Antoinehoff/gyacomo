
subroutine poisson
  ! Solve poisson equation to get phi

  USE time_integration, ONLY: updatetlevel
  USE array
  USE fields
  USE fourier_grid
  use model, ONLY : tau_e, tau_i, sigma_e, sigma_i, q_e, q_i, lambdaD

  USE prec_const
  IMPLICIT NONE

  INTEGER     :: ikr,ikz, ini,ine, i0j
  REAL(dp)    :: kr, kz
  REAL(dp)    :: k1_, k2_ ! sub kernel factor for recursive build
  REAL(dp)    :: x_e, x_i, alphaD
  COMPLEX(dp) :: gammaD, gammaphi
  COMPLEX(dp) :: sum_kernel2_e,    sum_kernel2_i
  COMPLEX(dp) :: sum_kernel_mom_e, sum_kernel_mom_i

  DO ikr=ikrs,ikre
    DO ikz=ikzs,ikze
      kr = krarray(ikr)
      kz = kzarray(ikz)

      ! Compute electrons sum(Kernel**2) and sum(Kernel * Ne0j)
      x_e    = sqrt(kr**2 + kz**2) * sigma_e * sqrt(2.0*tau_e)/2.0
      sum_kernel2_e    = 0
      sum_kernel_mom_e = 0
      k1_      = moments(1,ikr,ikz,updatetlevel)
      k2_      = 1.0
      if (jmaxe .ge. 0) then
        DO ine=1,jmaxe
          CALL bare(0,ine,i0j)
          k1_ = k1_ * x_e**2/ine
          k2_ = k2_ * x_e**4/ine**2
          sum_kernel_mom_e  = sum_kernel_mom_e  + k1_ * moments(i0j,ikr,ikz,updatetlevel)
          sum_kernel2_e     = sum_kernel2_e     + k2_
        END DO
      endif
      sum_kernel2_e    = sum_kernel2_e    * exp(-x_e**2)
      sum_kernel_mom_e = sum_kernel_mom_e * exp(-x_e**2)**2

      ! Compute ions sum(Kernel**2) and sum(Kernel * Ne0j)
      x_i    = sqrt(kr**2 + kz**2) * sigma_i * sqrt(2.0*tau_i)/2.0
      sum_kernel2_i    = 0
      sum_kernel_mom_i = 0
      k1_      = moments(Nmomi + 1, ikr, ikz, updatetlevel)
      k2_      = 1.0
      if (jmaxi .ge. 0) then
        DO ini=1,jmaxe
          CALL bari(0,ini,i0j)
          k1_ = k1_ * x_i**2/ini
          k2_ = k2_ * x_i**4/ini**2
          sum_kernel_mom_i  = sum_kernel_mom_i  + k1_ * moments(Nmome + i0j,ikr,ikz,updatetlevel)
          sum_kernel2_i     = sum_kernel2_i     + k2_
        END DO
      endif
      sum_kernel2_i    = sum_kernel2_i    * exp(-x_i**2)
      sum_kernel_mom_i = sum_kernel_mom_i * exp(-x_i**2)**2

      ! Assembling the poisson equation
      alphaD   = sqrt(kr**2 + kz**2) * lambdaD**2
      gammaD   = alphaD + q_e**2/tau_e * sum_kernel2_e &
                        + q_i**2/tau_i * sum_kernel2_i

      gammaphi = q_e * sum_kernel_mom_e + q_i * sum_kernel_mom_i
      
      phi(ikr, ikz) =  gammaphi/gammaD
    
    END DO
  END DO

END SUBROUTINE poisson
