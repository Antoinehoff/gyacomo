MODULE initial_par
  !   Module for initial parameters

  USE prec_const
  IMPLICIT NONE
  PRIVATE

  ! Initial background level
  REAL(dp), PUBLIC, PROTECTED ::  initback_density=0._dp
  REAL(dp), PUBLIC, PROTECTED ::  initback_temp=0._dp
  REAL(dp), PUBLIC, PROTECTED ::  initback_vpar=0._dp
  REAL(dp), PUBLIC, PROTECTED ::  initback_moments=0._dp

  ! Initial background noise amplitude
  REAL(dp), PUBLIC, PROTECTED ::  initnoise_density=1E-6_dp
  REAL(dp), PUBLIC, PROTECTED ::  initnoise_temp=1E-6_dp
  REAL(dp), PUBLIC, PROTECTED ::  initnoise_vpar=1E-6_dp
  REAL(dp), PUBLIC, PROTECTED ::  initnoise_moments=1E-6_dp

  ! Initialization for random number generator
  INTEGER, PUBLIC, PROTECTED :: iseed=42

  ! Parameters of initial smooth sine profiles
  REAL(dp), PUBLIC, PROTECTED ::  init_nb_oscil_density=2._dp ! Number of oscillations
  REAL(dp), PUBLIC, PROTECTED ::  init_nb_oscil_temp=2._dp 
  REAL(dp), PUBLIC, PROTECTED ::  init_ampli_density=0.1_dp ! Oscillation amplitude
  REAL(dp), PUBLIC, PROTECTED ::  init_ampli_temp=0.1_dp 

  PUBLIC :: initial_outputinputs, initial_readinputs

CONTAINS
  
  
  SUBROUTINE initial_readinputs
    ! Read the input parameters

    USE basic, ONLY : lu_in
    USE prec_const
    IMPLICIT NONE

    NAMELIST /INITIAL_CON/ initback_density, initback_temp, initback_vpar, initback_moments
    NAMELIST /INITIAL_CON/ initnoise_density, initnoise_temp, initnoise_vpar, initnoise_moments
    NAMELIST /INITIAL_CON/ iseed
    NAMELIST /INITIAL_CON/ init_nb_oscil_density, init_nb_oscil_temp, init_ampli_density, init_ampli_temp

    READ(lu_in,initial_con)
    WRITE(*,initial_con)

  END SUBROUTINE initial_readinputs


  SUBROUTINE initial_outputinputs(fidres, str)
    ! Write the input parameters to the results_xx.h5 file

    USE futils, ONLY: attach
    USE prec_const
    IMPLICIT NONE
    INTEGER, INTENT(in) :: fidres
    CHARACTER(len=256), INTENT(in) :: str

    CALL attach(fidres, TRIM(str), "initback_theta", initback_density)
    CALL attach(fidres, TRIM(str), "initback_temp", initback_temp)
    CALL attach(fidres, TRIM(str), "initback_vpar", initback_vpar)
    CALL attach(fidres, TRIM(str), "initback_moments", initback_moments)

    CALL attach(fidres, TRIM(str), "initnoise_theta", initnoise_density)
    CALL attach(fidres, TRIM(str), "initnoise_temp", initnoise_temp)
    CALL attach(fidres, TRIM(str), "initnoise_vpar", initnoise_vpar)
    CALL attach(fidres, TRIM(str), "initnoise_moments", initnoise_moments)

    CALL attach(fidres, TRIM(str), "iseed", iseed)

    CALL attach(fidres, TRIM(str), "init_nb_oscil_theta", init_nb_oscil_density)
    CALL attach(fidres, TRIM(str), "init_nb_oscil_temp", init_nb_oscil_temp)
    CALL attach(fidres, TRIM(str), "init_ampli_theta", init_ampli_density)
    CALL attach(fidres, TRIM(str), "init_ampli_temp", init_ampli_temp)

  END SUBROUTINE initial_outputinputs

END MODULE initial_par
