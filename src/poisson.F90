SUBROUTINE poisson
  ! Solve poisson equation to get phi

  USE time_integration, ONLY: updatetlevel
  USE array
  USE fields
  USE grid
  use model, ONLY : tau_e, tau_i, sigma_e, sigma_i, q_e, q_i, lambdaD, DK

  USE prec_const
  IMPLICIT NONE

  INTEGER     :: ini,ine, i0j
  REAL(dp)    :: ini_dp, ine_dp
  REAL(dp)    :: kr, kz, kperp2
  REAL(dp)    :: Kne, Kni ! sub kernel factor for recursive build
  REAL(dp)    :: be_2, bi_2, alphaD
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
      ! non dim kernel argument (kperp2 sigma_a sqrt(2 tau_a))
      be_2  =  kperp2 * sigmae2_taue_o2
      bi_2  =  kperp2 * sigmai2_taui_o2

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!! Electrons sum(Kernel * Ne0n)  (skm) and sum(Kernel**2) (sk2)
      !! Sum(kernel * Moments_0n)
      ! Initialization for n = 0 (ine = 1)
      Kne  = exp(-be_2)
      sum_kernel_mom_e = Kne*moments_e(1, 1, ikr, ikz, updatetlevel)
      ! loop over n only if the max polynomial degree is 1 or more
      if (jmaxe .GT. 0) then
        DO ine=2,jmaxe+1 ! ine = n+1
          ine_dp = REAL(ine-1,dp)
          ! We update iteratively the kernel functions (to spare factorial computations)
          Kne  = Kne  * be_2/ine_dp
          sum_kernel_mom_e  = sum_kernel_mom_e  + Kne * moments_e(1, ine, ikr, ikz, updatetlevel)
        END DO
      endif
      ! Initialization for n = 0 (ine = 1)
      Kne  = exp(-be_2)
      sum_kernel2_e = Kne**2
      ! loop over n only without caring of jmax since no moment dependency
      ! DO ine=2,10
      if (jmaxe .GT. 0) then
        DO ine=2,jmaxe+1 ! ine = n+1
          ine_dp        = REAL(ine-1,dp)         ! Real index (0 to jmax)
          Kne           = Kne  * be_2/ine_dp     ! update kernel_n
          sum_kernel2_e = sum_kernel2_e + Kne**2 ! ... sum recursively ...
        END DO
      ENDIF

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!! Ions sum(Kernel * Ni0n)  (skm) and sum(Kernel**2) (sk2)
      !! Sum(kernel * Moments_0n)
      ! Initialization for n = 0 (ini = 1)
        Kni  = exp(-bi_2)
        sum_kernel_mom_i = Kni*moments_i(1, 1, ikr, ikz, updatetlevel)
        ! loop over n only if the max polynomial degree is 1 or more
        if (jmaxi .GT. 0) then
          DO ini=2,jmaxi+1
            ini_dp = REAL(ini-1,dp) ! Real index (0 to jmax)
            ! We update iteratively to spare factorial computations
            Kni  = Kni  * bi_2/ini_dp
            sum_kernel_mom_i  = sum_kernel_mom_i  + Kni * moments_i(1, ini, ikr, ikz, updatetlevel)
          END DO
        endif

        ! Initialization for n = 0 (ini = 1)
        Kni  = exp(-bi_2)
        sum_kernel2_i = Kni**2
        ! loop over n only without caring of jmax since no moment dependency
        if (jmaxi .GT. 0) then
          DO ini=2,jmaxi+1
            ini_dp        = REAL(ini-1,dp)         ! Real index (0 to jmax)
            Kni           = Kni  * bi_2/ini_dp     ! update kernel_n
            sum_kernel2_i = sum_kernel2_i + Kni**2 ! ... sum recursively ...
          END DO
        ENDIF

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!! Assembling the poisson equation !!!!!!!!!!!!!!!!!!!!!!!!!!
      alphaD   = kperp2 * lambdaD**2
      gammaD   = alphaD + qe2_taue * (1._dp - sum_kernel2_e) & ! Called Poisson_ in MOLI
                        + qi2_taui * (1._dp - sum_kernel2_i)

      gammaD_phi = q_e * sum_kernel_mom_e + q_i * sum_kernel_mom_i

      IF ( (gammaD .EQ. 0) .AND. (abs(kr)+abs(kz) .NE. 0._dp) ) THEN
        write(*,*) 'Warning gammaD = 0', sum_kernel2_i
      ENDIF

      phi(ikr, ikz) =  gammaD_phi/gammaD

    END DO
  END DO

! Cancel origin singularity
phi(ikr_0,ikz_0) = 0

END SUBROUTINE poisson
