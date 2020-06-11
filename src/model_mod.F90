MODULE model
  ! Module for diagnostic parameters

  USE prec_const
  IMPLICIT NONE
  PRIVATE

  REAL(dp), PUBLIC, PROTECTED :: nu = 1._dp     ! Collision frequency

  ! (Numerical) diffusion coefficient
  REAL(dp), PUBLIC, PROTECTED ::  diff_theta = 1._dp
  REAL(dp), PUBLIC, PROTECTED ::  diff_temp = 1._dp
  REAL(dp), PUBLIC, PROTECTED ::  diff_vpar = 1._dp
  REAL(dp), PUBLIC, PROTECTED ::  diff_moments = 1._dp

  ! Fast Fourier Transform to filter out high frequence
  LOGICAL, PUBLIC, PROTECTED :: fft_suppress_high_freq = .false.
  REAL(dp), PUBLIC, PROTECTED ::  fft_sigma = 2._dp


  CHARACTER(len=3), PUBLIC, PROTECTED :: gradient_scheme='fd4'

  LOGICAL, PUBLIC, PROTECTED :: freeze_theta = .false.
  LOGICAL, PUBLIC, PROTECTED :: freeze_temp = .false.
  LOGICAL, PUBLIC, PROTECTED :: freeze_vpar = .false.
  LOGICAL, PUBLIC, PROTECTED :: freeze_moments = .false.
  LOGICAL, PUBLIC, PROTECTED :: freeze_phi = .false.

    PUBLIC :: model_readinputs, model_outputinputs

CONTAINS

  SUBROUTINE model_readinputs
    !    Read the input parameters

    USE basic, ONLY : lu_in
    USE prec_const
    IMPLICIT NONE

    NAMELIST /MODEL_PAR/ nu
    NAMELIST /MODEL_PAR/ diff_theta, diff_temp, diff_vpar, diff_moments
    NAMELIST /MODEL_PAR/ gradient_scheme
    NAMELIST /MODEL_PAR/ freeze_theta, freeze_temp, freeze_vpar, freeze_moments, freeze_phi
    NAMELIST /MODEL_PAR/ fft_suppress_high_freq, fft_sigma

    READ(lu_in,model_par)
    WRITE(*,model_par)

    ! Collision Frequency Normalization ... to match fluid limit
    nu = nu*0.532_dp

  END SUBROUTINE model_readinputs


  SUBROUTINE model_outputinputs(fidres, str)
    !    Write the input parameters to the results_xx.h5 file

    USE futils, ONLY: attach
    USE prec_const
    IMPLICIT NONE
    INTEGER, INTENT(in) :: fidres
    CHARACTER(len=256), INTENT(in) :: str

    CALL attach(fidres, TRIM(str), "nu", nu)
    CALL attach(fidres, TRIM(str), "diff_theta", diff_theta)
    CALL attach(fidres, TRIM(str), "diff_temp", diff_temp)
    CALL attach(fidres, TRIM(str), "diff_vpar", diff_vpar)
    CALL attach(fidres, TRIM(str), "diff_moments", diff_moments)
    CALL attach(fidres, TRIM(str), "gradient_scheme", gradient_scheme)
    CALL attach(fidres, TRIM(str), "freeze_theta", freeze_theta)
    CALL attach(fidres, TRIM(str), "freeze_temp", freeze_temp)
    CALL attach(fidres, TRIM(str), "freeze_vpar", freeze_vpar)
    CALL attach(fidres, TRIM(str), "freeze_moments", freeze_moments)
    CALL attach(fidres, TRIM(str), "freeze_phi", freeze_phi)
    CALL attach(fidres, TRIM(str), "fft_suppress_high_freq", fft_suppress_high_freq)
    CALL attach(fidres, TRIM(str), "fft_sigma", fft_sigma)

  END SUBROUTINE model_outputinputs

END MODULE model
