MODULE initial_par
  !   Module for initial parameters

  USE prec_const
  IMPLICIT NONE
  PRIVATE

  ! Initialization through a noisy phi
  LOGICAL,  PUBLIC, PROTECTED :: INIT_NOISY_PHI = .false.
  ! Initialization through a zonal flow phi
  INTEGER,  PUBLIC, PROTECTED :: INIT_ZF_PHI    = 0
  ! Initial background level
  REAL(dp), PUBLIC, PROTECTED :: init_background=0._dp
  ! Initial noise amplitude
  REAL(dp), PUBLIC, PROTECTED :: init_noiselvl=1E-6_dp
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

    NAMELIST /INITIAL_CON/ INIT_NOISY_PHI
    NAMELIST /INITIAL_CON/ INIT_ZF_PHI
    NAMELIST /INITIAL_CON/ init_background
    NAMELIST /INITIAL_CON/ init_noiselvl
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

    CALL attach(fidres, TRIM(str), "INIT_NOISY_PHI", INIT_NOISY_PHI)

    CALL attach(fidres, TRIM(str), "INIT_ZF_PHI", INIT_ZF_PHI)

    CALL attach(fidres, TRIM(str), "init_background", init_background)

    CALL attach(fidres, TRIM(str), "init_noiselvl", init_noiselvl)

    CALL attach(fidres, TRIM(str), "iseed", iseed)

  END SUBROUTINE initial_outputinputs

END MODULE initial_par
