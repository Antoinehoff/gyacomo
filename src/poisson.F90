SUBROUTINE poisson
  ! Solve poisson equation to get phi

  USE basic
  USE time_integration, ONLY: updatetlevel
  USE array
  USE fields
  USE grid
  use model, ONLY : tau_e, tau_i, sigma_e, sigma_i, q_e, q_i, lambdaD, DK

  USE prec_const
  IMPLICIT NONE

  INTEGER     :: ini,ine
  REAL(dp)    :: Kne, Kni ! sub kernel factor for recursive build
  REAL(dp)    :: alphaD
  REAL(dp)    :: qe2_taue, qi2_taui ! To avoid redondant computation
  REAL(dp)    :: sum_kernel2_e,    sum_kernel2_i    ! Store sum Kn^2
  COMPLEX(dp) :: sum_kernel_mom_e, sum_kernel_mom_i ! Store sum Kn*Napn
  REAL(dp)    :: gammaD
  COMPLEX(dp) :: gammaD_phi

  ! Execution time start
  CALL cpu_time(t0_poisson)

  !Precompute species dependant factors
  qe2_taue        = (q_e**2)/tau_e ! factor of the gammaD sum
  qi2_taui        = (q_i**2)/tau_i

  DO ikr=ikrs,ikre
    DO ikz=ikzs,ikze

      !!!!!!!!!!!!! Electrons sum(Kernel * Ne0n)  (skm) and sum(Kernel**2) (sk2)
      sum_kernel_mom_e = 0._dp
      sum_kernel2_e    = 0._dp
      ! loop over n only if the max polynomial degree
      DO ine=1,jmaxe+1 ! ine = n+1
        Kne = kernel_e(ine, ikr, ikz)
        sum_kernel_mom_e  = sum_kernel_mom_e  + Kne * moments_e(1, ine, ikr, ikz, updatetlevel)
        sum_kernel2_e     = sum_kernel2_e     + Kne**2 ! ... sum recursively ...
      END DO

      !!!!!!!!!!!!!!!!! Ions sum(Kernel * Ni0n)  (skm) and sum(Kernel**2) (sk2)
        sum_kernel_mom_i = 0
        sum_kernel2_i = 0
        ! loop over n only if the max polynomial degree
        DO ini=1,jmaxi+1
          Kni = kernel_i(ini, ikr, ikz)
          sum_kernel_mom_i  = sum_kernel_mom_i  + Kni * moments_i(1, ini, ikr, ikz, updatetlevel)
          sum_kernel2_i     = sum_kernel2_i     + Kni**2 ! ... sum recursively ...
        END DO

      !!!!!!!!!!!!!!! Assembling the poisson equation !!!!!!!!!!!!!!!!!!!!!!!!!!
      alphaD   = (krarray(ikr)**2 + kzarray(ikz)**2) * lambdaD**2
      gammaD   = alphaD + qe2_taue * (1._dp - sum_kernel2_e) & ! Called Poisson_ in MOLI
                        + qi2_taui * (1._dp - sum_kernel2_i)

      gammaD_phi = q_e * sum_kernel_mom_e + q_i * sum_kernel_mom_i

      phi(ikr, ikz) =  gammaD_phi/gammaD

    END DO
  END DO

! Cancel origin singularity
IF ((ikr_0 .GE. ikrs) .AND. (ikr_0 .LE. ikre)) phi(ikr_0,ikz_0) = 0

! Execution time end
CALL cpu_time(t1_poisson)
tc_poisson = tc_poisson + (t1_poisson - t0_poisson)

END SUBROUTINE poisson
