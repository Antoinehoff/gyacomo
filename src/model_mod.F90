MODULE model
  ! Module for diagnostic parameters
  USE prec_const
  IMPLICIT NONE
  PRIVATE
  ! INPUTS
  INTEGER,  PUBLIC, PROTECTED ::    KERN =  0         ! Kernel model
  CHARACTER(len=16), &
            PUBLIC, PROTECTED ::LINEARITY= 'linear'   ! To turn on non linear bracket term
  REAL(xp), PUBLIC, PROTECTED ::    mu_x =  0._xp     ! spatial    x-Hyperdiffusivity coefficient (for num. stability)
  REAL(xp), PUBLIC, PROTECTED ::    mu_y =  0._xp     ! spatial    y-Hyperdiffusivity coefficient (for num. stability)
  INTEGER,  PUBLIC, PROTECTED ::    N_HD =  4         ! order of numerical spatial diffusion
  LOGICAL,  PUBLIC, PROTECTED ::   HDz_h =  .false.    ! to apply z-hyperdiffusion on non adiab part
  REAL(xp), PUBLIC, PROTECTED ::    mu_z =  0._xp     ! spatial    z-Hyperdiffusivity coefficient (for num. stability)
  REAL(xp), PUBLIC, PROTECTED ::    mu_p =  0._xp     ! kinetic para hyperdiffusivity coefficient (for num. stability)
  REAL(xp), PUBLIC, PROTECTED ::    mu_j =  0._xp     ! kinetic perp hyperdiffusivity coefficient (for num. stability)
  CHARACTER(len=16), &
  PUBLIC, PROTECTED ::   HYP_V = 'hypcoll'  ! hyperdiffusion model for velocity space ('none','hypcoll','dvpar4')
  INTEGER,  PUBLIC, PROTECTED ::      Na =  1         ! number of evolved species
  REAL(xp), PUBLIC, PROTECTED ::      nu =  0._xp     ! collision frequency parameter
  REAL(xp), PUBLIC, PROTECTED ::    k_gB =  1._xp     ! Magnetic gradient strength (L_ref/L_gB)
  REAL(xp), PUBLIC, PROTECTED ::    k_cB =  1._xp     ! Magnetic curvature strength (L_ref/L_cB)
  REAL(xp), PUBLIC, PROTECTED :: lambdaD =  0._xp     ! Debye length
  REAL(xp), PUBLIC, PROTECTED ::    beta =  0._xp     ! electron plasma Beta (8piNT_e/B0^2)
  LOGICAL,  PUBLIC            :: ADIAB_E =  .false.   ! adiabatic electron model
  LOGICAL,  PUBLIC            :: ADIAB_I =  .false.   ! adiabatic ion model
  REAL(xp), PUBLIC, PROTECTED ::   tau_i =  1.0       ! electron-ion temperature ratio for ion adiabatic model
  ! Auxiliary variable
  LOGICAL,  PUBLIC, PROTECTED ::      EM =  .false.   ! Electromagnetic effects flag
  LOGICAL,  PUBLIC, PROTECTED ::  MHD_PD =  .false.   ! MHD pressure drift

  ! Removes Landau damping in temperature and higher equation (Ivanov 2022)
  LOGICAL,  PUBLIC, PROTECTED :: RM_LD_T_EQ = .false.
  PUBLIC :: model_readinputs, model_outputinputs

CONTAINS

  SUBROUTINE model_readinputs
    !    Read the input parameters
    USE basic,    ONLY: lu_in, speak
    USE parallel, ONLY: num_procs_p
    USE prec_const
    IMPLICIT NONE

    NAMELIST /MODEL_PAR/ KERN, LINEARITY, RM_LD_T_EQ, &
                         mu_x, mu_y, N_HD, HDz_h, mu_z, mu_p, mu_j, HYP_V, Na,&
                         nu, k_gB, k_cB, lambdaD, MHD_PD, beta, ADIAB_E, ADIAB_I, tau_i

    READ(lu_in,model_par)

    IF (ADIAB_E .AND. ADIAB_I) &
      ERROR STOP '>> ERROR << cannot have both adiab e and adiab i models'

    IF((HYP_V .EQ. 'dvpar4') .AND. (num_procs_p .GT. 1)) &
      ERROR STOP '>> ERROR << dvpar4 velocity dissipation is not compatible with current p parallelization'

    IF(Na .EQ. 1) THEN
      IF((.NOT. ADIAB_E) .AND. (.NOT. ADIAB_I)) ERROR STOP "With one species, ADIAB_E or ADIAB_I must be set to .true. STOP"
      IF(ADIAB_E) THEN
        CALL speak('Adiabatic electron model -> beta = 0')
        beta = 0._xp
      ENDIF
      IF(ADIAB_I) CALL speak('Adiabatic ion model')
    ENDIF

    IF(beta .GT. 0) THEN
      CALL speak('Electromagnetic effects are included')
      EM   = .TRUE.
    ENDIF

  END SUBROUTINE model_readinputs

  SUBROUTINE model_outputinputs(fid)
    ! Write the input parameters to the results_xx.h5 file
    USE futils, ONLY: attach, creatd
    IMPLICIT NONE
    INTEGER, INTENT(in) :: fid
    CHARACTER(len=256)  :: str
    WRITE(str,'(a)') '/data/input/model'
    CALL creatd(fid, 0,(/0/),TRIM(str),'Model Input')
    CALL attach(fid, TRIM(str),      "KERN",    KERN)
    CALL attach(fid, TRIM(str), "LINEARITY", LINEARITY)
    CALL attach(fid, TRIM(str),"RM_LD_T_EQ",RM_LD_T_EQ)
    CALL attach(fid, TRIM(str),      "mu_x",    mu_x)
    CALL attach(fid, TRIM(str),      "mu_y",    mu_y)
    CALL attach(fid, TRIM(str),      "N_HD",    N_HD)
    CALL attach(fid, TRIM(str),      "mu_z",    mu_z)
    CALL attach(fid, TRIM(str),     "HDz_h",   HDz_h)
    CALL attach(fid, TRIM(str),      "mu_p",    mu_p)
    CALL attach(fid, TRIM(str),      "mu_j",    mu_j)
    CALL attach(fid, TRIM(str),     "HYP_V",   HYP_V)
    CALL attach(fid, TRIM(str),        "Na",      Na)
    CALL attach(fid, TRIM(str),        "nu",      nu)
    CALL attach(fid, TRIM(str),      "k_gB",    k_gB)
    CALL attach(fid, TRIM(str),      "k_cB",    k_cB)
    CALL attach(fid, TRIM(str),   "lambdaD", lambdaD)
    CALL attach(fid, TRIM(str),    "MHD_PD",  MHD_PD)
    CALL attach(fid, TRIM(str),      "beta",    beta)
    CALL attach(fid, TRIM(str),   "ADIAB_E", ADIAB_E)
    CALL attach(fid, TRIM(str),   "ADIAB_I", ADIAB_I)
    CALL attach(fid, TRIM(str),     "tau_i",   tau_i)
  END SUBROUTINE model_outputinputs

END MODULE model
