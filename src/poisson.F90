SUBROUTINE poisson
  ! Solve poisson equation to get phi

  USE time_integration, ONLY: updatetlevel
  USE array
  USE fields
  USE grid
  use model, ONLY : tau_e, tau_i, sigma_e, sigma_i, q_e, q_i, lambdaD

  USE prec_const
  IMPLICIT NONE

  INTEGER     :: ini,ine, i0j
  REAL(dp)    :: ini_dp, ine_dp
  REAL(dp)    :: kr, kz, kperp2
  REAL(dp)    :: Kne, Kni ! sub kernel factor for recursive build
  REAL(dp)    :: b_e2, b_i2, alphaD
  REAL(dp)    :: sigmae2_taue_o2, sigmai2_taui_o2, qe2_taue, qi2_taui ! To avoid redondant computation
  REAL(dp)    :: sum_kernel2_e,    sum_kernel2_i    ! Store sum Kn^2
  COMPLEX(dp) :: sum_kernel_mom_e, sum_kernel_mom_i ! Store sum Kn*Napn
  REAL(dp)    :: gammaD
  COMPLEX(dp) :: gammaD_phi

  !Precompute species dependant factors
  sigmae2_taue_o2 = sigma_e**2 * tau_e/2._dp ! factor of the Kernel argument
  sigmai2_taui_o2 = sigma_i**2 * tau_i/2._dp ! (b_a/2)^2 = (kperp sqrt(2 tau_a) sigma_a/2)^2
                                             !           = kperp^2 tau_a sigma_a^2/2
  qe2_taue        = (q_e**2)/tau_e ! factor of the gammaD sum
  qi2_taui        = (q_i**2)/tau_i

  DO ikr=ikrs,ikre
    DO ikz=ikzs,ikze
      kr     = krarray(ikr)
      kz     = kzarray(ikz)
      kperp2 = kr**2 + kz**2

      !! Electrons sum(Kernel * Ne0n)  (skm) and sum(Kernel**2) (sk2)
      b_e2  =  kperp2 * sigmae2_taue_o2 ! non dim kernel argument (kperp2 sigma_a sqrt(2 tau_a))

      ! Initialization for n = 0 (ine = 1)
      Kne  = exp(-b_e2)
      sum_kernel_mom_e = Kne*moments_e(1, 1, ikr, ikz, updatetlevel)
      sum_kernel2_e    = Kne**2

!write(*,*) 'K0e = ', Kne

      ! loop over n only if the max polynomial degree is 1 or more
      if (jmaxe .GT. 0) then
        DO ine=2,jmaxe+1 ! ine = n+1
          ine_dp = REAL(ine-1,dp)

          Kne  = Kne  * b_e2/ine_dp       ! We update iteratively the kernel functions (to spare factorial computations)

          sum_kernel_mom_e  = sum_kernel_mom_e  + Kne * moments_e(1, ine, ikr, ikz, updatetlevel)
          sum_kernel2_e     = sum_kernel2_e     + Kne**2 ! ... sum recursively ...
!write(*,*) 'K',ine-1,'e = ', Kne
        END DO
      endif

      !! Ions sum(Kernel * Ni0n)  (skm) and sum(Kernel**2) (sk2)
      b_i2  = kperp2 * sigmai2_taui_o2

      ! Initialization for n = 0 (ini = 1)
      Kni  = exp(-b_i2)
      sum_kernel_mom_i = Kni*moments_i(1, 1, ikr, ikz, updatetlevel)
      sum_kernel2_i    = Kni**2

!write(*,*) 'K0i = ', Kni

      ! loop over n only if the max polynomial degree is 1 or more
      if (jmaxi .GT. 0) then
        DO ini=2,jmaxi+1
          ini_dp = REAL(ini-1,dp) ! Real index (0 to jmax)

          Kni  = Kni  * b_i2/ini_dp       ! We update iteratively the kernel functions this time for ions ...

          sum_kernel_mom_i  = sum_kernel_mom_i  + Kni * moments_i(1, ini, ikr, ikz, updatetlevel)
          sum_kernel2_i     = sum_kernel2_i     + Kni**2 ! ... sum recursively ...
!write(*,*) 'K',ini-1,'i = ', Kni
        END DO
      endif
      !!! Assembling the poisson equation
      alphaD   = kperp2 * lambdaD**2
      gammaD   = alphaD + qe2_taue * (1._dp - sum_kernel2_e) & ! Called Poisson_ in MOLI
                        + qi2_taui * (1._dp - sum_kernel2_i)
      gammaD_phi = q_e * sum_kernel_mom_e + q_i * sum_kernel_mom_i

      IF ( (alphaD .EQ. 0._dp) .AND. (gammaD .EQ. 0._dp) ) THEN
        write(*,*) "Warning : 0/0 occuring"
        phi(ikr,ikz) = 0._dp
      ELSE
        phi(ikr, ikz) =  gammaD_phi/gammaD
      ENDIF

    END DO
  END DO

!write(*,*) 'gammaD =',gammaD

END SUBROUTINE poisson
