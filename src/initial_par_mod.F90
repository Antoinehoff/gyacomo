MODULE initial_par
  !   Module for initial parameters

  USE prec_const
  IMPLICIT NONE
  PRIVATE

  ! Initialization through a noisy phi
  LOGICAL,  PUBLIC, PROTECTED :: INIT_NOISY_PHI = .false.
  ! Initialization through a zonal flow phi
  INTEGER,  PUBLIC, PROTECTED :: INIT_ZF    = 0
  REAL(DP), PUBLIC, PROTECTED :: ZF_AMP     = 1E+3_dp
  ! Act on modes artificially (keep/wipe, zonal, non zonal, entropy mode etc.)
  CHARACTER(len=32),  PUBLIC, PROTECTED :: ACT_ON_MODES = 'nothing'
  ! Wipe turbulence in the restart (=1) or at each step (=2)
  INTEGER,  PUBLIC, PROTECTED :: WIPE_TURB = 0
  ! Init a Gaussian blob density in the middle
  LOGICAL,  PUBLIC, PROTECTED :: INIT_BLOB = .false.
  ! Initial background level
  REAL(dp), PUBLIC, PROTECTED :: init_background=0._dp
  ! Initial noise amplitude
  REAL(dp), PUBLIC, PROTECTED :: init_noiselvl=1E-6_dp
  ! Initialization for random number generator
  INTEGER,  PUBLIC, PROTECTED :: iseed=42

  PUBLIC :: initial_outputinputs, initial_readinputs

CONTAINS


  SUBROUTINE initial_readinputs
    ! Read the input parameters

    USE basic, ONLY : lu_in
    USE prec_const
    IMPLICIT NONE

    NAMELIST /INITIAL_CON/ INIT_NOISY_PHI
    NAMELIST /INITIAL_CON/ INIT_ZF
    NAMELIST /INITIAL_CON/ ACT_ON_MODES
    NAMELIST /INITIAL_CON/ WIPE_TURB
    NAMELIST /INITIAL_CON/ INIT_BLOB
    NAMELIST /INITIAL_CON/ init_background
    NAMELIST /INITIAL_CON/ init_noiselvl
    NAMELIST /INITIAL_CON/ iseed

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

    CALL attach(fidres, TRIM(str), "INIT_ZF", INIT_ZF)

    CALL attach(fidres, TRIM(str), "init_background", init_background)

    CALL attach(fidres, TRIM(str), "init_noiselvl", init_noiselvl)

    CALL attach(fidres, TRIM(str), "iseed", iseed)

  END SUBROUTINE initial_outputinputs

END MODULE initial_par
