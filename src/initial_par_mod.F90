MODULE initial_par
  !   Module for initial parameters

  USE prec_const
  IMPLICIT NONE
  PRIVATE

  ! Initialization option (phi/mom00/allmom/blob)
  CHARACTER(len=32), PUBLIC, PROTECTED :: INIT_OPT = 'phi'
  ! Initialization through a zonal flow phi
  INTEGER,  PUBLIC, PROTECTED :: INIT_ZF    = 0
  REAL(xp), PUBLIC, PROTECTED :: ZF_AMP     = 1E+3_xp
  ! Act on modes artificially (keep/wipe, zonal, non zonal, entropy mode etc.)
  CHARACTER(len=32),  PUBLIC, PROTECTED :: ACT_ON_MODES = 'nothing'
  ! Initial background level
  REAL(xp), PUBLIC, PROTECTED :: init_background=0._xp
  ! Initial noise amplitude
  REAL(xp), PUBLIC, PROTECTED :: init_noiselvl=1E-6_xp
  ! Initialization for random number generator
  INTEGER,  PUBLIC, PROTECTED :: iseed=42

  PUBLIC :: initial_outputinputs, initial_readinputs

CONTAINS


  SUBROUTINE initial_readinputs
    ! Read the input parameters

    USE basic, ONLY : lu_in
    USE prec_const
    IMPLICIT NONE

    NAMELIST /INITIAL_CON/ INIT_OPT
    NAMELIST /INITIAL_CON/ ACT_ON_MODES
    NAMELIST /INITIAL_CON/ init_background
    NAMELIST /INITIAL_CON/ init_noiselvl
    NAMELIST /INITIAL_CON/ iseed

    READ(lu_in,initial_con)
    !WRITE(*,initial_con)

  END SUBROUTINE initial_readinputs


  SUBROUTINE initial_outputinputs(fid)
    ! Write the input parameters to the results_xx.h5 file
    USE futils, ONLY: attach, creatd
    IMPLICIT NONE
    INTEGER, INTENT(in) :: fid
    CHARACTER(len=256)  :: str
    WRITE(str,'(a)') '/data/input/intial'
    CALL creatd(fid, 0,(/0/),TRIM(str),'Initial Input')
    CALL attach(fid, TRIM(str), "INIT_OPT", INIT_OPT)
    CALL attach(fid, TRIM(str), "init_background", init_background)
    CALL attach(fid, TRIM(str), "init_noiselvl", init_noiselvl)
    CALL attach(fid, TRIM(str), "iseed", iseed)
  END SUBROUTINE initial_outputinputs

END MODULE initial_par
