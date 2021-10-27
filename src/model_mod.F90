MODULE model
  ! Module for diagnostic parameters

  USE prec_const
  IMPLICIT NONE
  PRIVATE

  INTEGER,  PUBLIC, PROTECTED ::      CO =  0         ! Collision Operator
  INTEGER,  PUBLIC, PROTECTED ::    CLOS =  0         ! linear truncation method
  INTEGER,  PUBLIC, PROTECTED :: NL_CLOS =  0         ! nonlinear truncation method
  INTEGER,  PUBLIC, PROTECTED ::    KERN =  0         ! Kernel model
  LOGICAL,  PUBLIC, PROTECTED :: NON_LIN =  .true.    ! To turn on non linear bracket term
  REAL(dp), PUBLIC, PROTECTED ::      mu =  0._dp     ! spatial      Hyperdiffusivity coefficient (for num. stability)
  REAL(dp), PUBLIC, PROTECTED ::    mu_p =  0._dp     ! kinetic para hyperdiffusivity coefficient (for num. stability)
  REAL(dp), PUBLIC, PROTECTED ::    mu_j =  0._dp     ! kinetic perp hyperdiffusivity coefficient (for num. stability)
  REAL(dp), PUBLIC, PROTECTED ::      nu =  1._dp     ! Collision frequency
  REAL(dp), PUBLIC, PROTECTED ::   tau_e =  1._dp     ! Temperature
  REAL(dp), PUBLIC, PROTECTED ::   tau_i =  1._dp     !
  REAL(dp), PUBLIC, PROTECTED :: sigma_e =  1._dp     ! Mass
  REAL(dp), PUBLIC, PROTECTED :: sigma_i =  1._dp     !
  REAL(dp), PUBLIC, PROTECTED ::     q_e = -1._dp     ! Charge
  REAL(dp), PUBLIC, PROTECTED ::     q_i =  1._dp     !
  REAL(dp), PUBLIC, PROTECTED ::     K_n =  1._dp     ! Density drive
  REAL(dp), PUBLIC, PROTECTED ::     K_T =  0._dp     ! Temperature drive
  REAL(dp), PUBLIC, PROTECTED ::     K_E =  0._dp     ! Backg. electric field drive
  REAL(dp), PUBLIC, PROTECTED ::   GradB =  1._dp     ! Magnetic gradient
  REAL(dp), PUBLIC, PROTECTED ::   CurvB =  1._dp     ! Magnetic curvature
  REAL(dp), PUBLIC, PROTECTED :: lambdaD =  1._dp     ! Debye length

  REAL(dp), PUBLIC, PROTECTED :: taue_qe         ! factor of the magnetic moment coupling
  REAL(dp), PUBLIC, PROTECTED :: taui_qi         !
  REAL(dp), PUBLIC, PROTECTED :: qi_taui         !
  REAL(dp), PUBLIC, PROTECTED :: qe_taue         !
  REAL(dp), PUBLIC, PROTECTED :: sqrtTaue_qe          ! factor of parallel moment term
  REAL(dp), PUBLIC, PROTECTED :: sqrtTaui_qi          !
  REAL(dp), PUBLIC, PROTECTED :: qe_sigmae_sqrtTaue   ! factor of parallel phi term
  REAL(dp), PUBLIC, PROTECTED :: qi_sigmai_sqrtTaui   !
  REAL(dp), PUBLIC, PROTECTED :: sigmae2_taue_o2      ! factor of the Kernel argument
  REAL(dp), PUBLIC, PROTECTED :: sigmai2_taui_o2      !
  REAL(dp), PUBLIC, PROTECTED :: sqrt_sigmae2_taue_o2      ! factor of the Kernel argument
  REAL(dp), PUBLIC, PROTECTED :: sqrt_sigmai2_taui_o2  
  REAL(dp), PUBLIC, PROTECTED :: nu_e,  nu_i          ! electron-ion, ion-ion collision frequency
  REAL(dp), PUBLIC, PROTECTED :: nu_ee, nu_ie         ! e-e, i-e coll. frequ.
  REAL(dp), PUBLIC, PROTECTED :: qe2_taue, qi2_taui   ! factor of the gammaD sum

  PUBLIC :: model_readinputs, model_outputinputs

CONTAINS

  SUBROUTINE model_readinputs
    !    Read the input parameters
    USE basic, ONLY : lu_in
    USE prec_const
    IMPLICIT NONE

    NAMELIST /MODEL_PAR/ CO, CLOS, NL_CLOS, KERN, NON_LIN, mu, mu_p, mu_j, nu, tau_e, tau_i, sigma_e, sigma_i, &
                         q_e, q_i, K_n, K_T, K_E, GradB, CurvB, lambdaD

    READ(lu_in,model_par)

    taue_qe    = tau_e/q_e ! factor of the magnetic moment coupling
    taui_qi    = tau_i/q_i ! factor of the magnetic moment coupling
    qe_taue    = q_e/tau_e
    qi_taui    = q_i/tau_i
    sqrtTaue_qe     = sqrt(tau_e)/q_e   ! factor of parallel moment term
    sqrtTaui_qi     = sqrt(tau_i)/q_i   ! factor of parallel moment term
    qe_sigmae_sqrtTaue = q_e/sigma_e/SQRT(tau_e) ! factor of parallel phi term
    qi_sigmai_sqrtTaui = q_i/sigma_i/SQRT(tau_i)
    qe2_taue        = (q_e**2)/tau_e ! factor of the gammaD sum
    qi2_taui        = (q_i**2)/tau_i
    sigmae2_taue_o2 = sigma_e**2 * tau_e/2._dp ! factor of the Kernel argument
    sigmai2_taui_o2 = sigma_i**2 * tau_i/2._dp
    sqrt_sigmae2_taue_o2 = SQRT(sigma_e**2 * tau_e/2._dp) ! to avoid multiple SQRT eval
    sqrt_sigmai2_taui_o2 = SQRT(sigma_i**2 * tau_i/2._dp)
    !! We use the ion-ion collision as normalization with definition
    !   nu_ii = 4 sqrt(pi)/3 T_i^(-3/2) m_i^(-1/2) q^4 n_i0 ln(Lambda)
    !
    nu_e            = nu/sigma_e * (tau_e)**(3._dp/2._dp)  ! electron-ion collision frequency (where already multiplied by 0.532)
    nu_i            = nu ! ion-ion collision frequ.
    nu_ee           = nu_e ! e-e coll. frequ.
    nu_ie           = nu_i ! i-e coll. frequ.

    ! Old normalization (MOLI Jorge/Frei)
    ! nu_e            = 0.532_dp*nu ! electron-ion collision frequency (where already multiplied by 0.532)
    ! nu_i            = 0.532_dp*nu*sigma_e*tau_e**(-3._dp/2._dp)/SQRT2 ! ion-ion collision frequ.
    ! nu_ee           = nu_e/SQRT2 ! e-e coll. frequ.
    ! nu_ie           = 0.532_dp*nu*sigma_e**2 ! i-e coll. frequ.

  END SUBROUTINE model_readinputs


  SUBROUTINE model_outputinputs(fidres, str)
    !    Write the input parameters to the results_xx.h5 file

    USE futils, ONLY: attach
    USE prec_const
    IMPLICIT NONE

    INTEGER, INTENT(in) :: fidres
    CHARACTER(len=256), INTENT(in) :: str
    CALL attach(fidres, TRIM(str),        "CO",      CO)
    CALL attach(fidres, TRIM(str),      "CLOS",    CLOS)
    CALL attach(fidres, TRIM(str),      "KERN",    KERN)
    CALL attach(fidres, TRIM(str),   "NON_LIN", NON_LIN)
    CALL attach(fidres, TRIM(str),        "nu",      nu)
    CALL attach(fidres, TRIM(str),        "mu",      mu)
    CALL attach(fidres, TRIM(str),     "tau_e",   tau_e)
    CALL attach(fidres, TRIM(str),     "tau_i",   tau_i)
    CALL attach(fidres, TRIM(str),   "sigma_e", sigma_e)
    CALL attach(fidres, TRIM(str),   "sigma_i", sigma_i)
    CALL attach(fidres, TRIM(str),       "q_e",     q_e)
    CALL attach(fidres, TRIM(str),       "q_i",     q_i)
    CALL attach(fidres, TRIM(str),   "K_n", K_n)
    CALL attach(fidres, TRIM(str),   "K_T", K_T)
    CALL attach(fidres, TRIM(str), "K_E", K_E)
    CALL attach(fidres, TRIM(str),   "lambdaD", lambdaD)
  END SUBROUTINE model_outputinputs

END MODULE model
