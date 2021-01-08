MODULE initial_par
  !   Module for initial parameters

  USE prec_const
  IMPLICIT NONE
  PRIVATE

  ! Initial background level
  REAL(dp), PUBLIC, PROTECTED :: initback_moments=0._dp
  ! Initial background level
  LOGICAL,  PUBLIC, PROTECTED :: only_Na00 = .false.
  ! Initial background noise amplitude
  REAL(dp), PUBLIC, PROTECTED :: initnoise_moments=1E-6_dp
  ! Initialization for random number generator
  INTEGER,  PUBLIC, PROTECTED :: iseed=42

  ! Parameters of initial smooth sine profiles
  REAL(dp), PUBLIC, PROTECTED ::  init_nb_oscil_density=2._dp ! Number of oscillations
  REAL(dp), PUBLIC, PROTECTED ::  init_nb_oscil_temp=2._dp
  REAL(dp), PUBLIC, PROTECTED ::  init_ampli_density=0.1_dp ! Oscillation amplitude
  REAL(dp), PUBLIC, PROTECTED ::  init_ampli_temp=0.1_dp

  CHARACTER(len=128), PUBLIC :: selfmat_file  ! COSOlver matrix file names
  CHARACTER(len=128), PUBLIC :: iemat_file  ! COSOlver matrix file names
  CHARACTER(len=128), PUBLIC :: eimat_file  ! COSOlver matrix file names

  PUBLIC :: initial_outputinputs, initial_readinputs

CONTAINS


  SUBROUTINE initial_readinputs
    ! Read the input parameters

    USE basic, ONLY : lu_in, RESTART
    USE prec_const
    IMPLICIT NONE

    NAMELIST /INITIAL_CON/ only_Na00
    NAMELIST /INITIAL_CON/ initback_moments
    NAMELIST /INITIAL_CON/ initnoise_moments
    NAMELIST /INITIAL_CON/ iseed
    NAMELIST /INITIAL_CON/ selfmat_file
    NAMELIST /INITIAL_CON/ iemat_file
    NAMELIST /INITIAL_CON/ eimat_file

    READ(lu_in,initial_con)
    !WRITE(*,initial_con)

  END SUBROUTINE initial_readinputs


  SUBROUTINE initial_outputinputs(fidres, str)
    ! Write the input parameters to the results_xx.h5 file

    USE futils, ONLY: attach
    USE prec_const
    IMPLICIT NONE
    INTEGER, INTENT(in) :: fidres
    CHARACTER(len=256), INTENT(in) :: str

    CALL attach(fidres, TRIM(str), "only_Na00", only_Na00)

    CALL attach(fidres, TRIM(str), "initback_moments", initback_moments)

    CALL attach(fidres, TRIM(str), "initnoise_moments", initnoise_moments)

    CALL attach(fidres, TRIM(str), "iseed", iseed)

  END SUBROUTINE initial_outputinputs

END MODULE initial_par
