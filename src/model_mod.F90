MODULE model
  ! Module for diagnostic parameters
  USE prec_const
  IMPLICIT NONE
  PRIVATE
  ! INPUTS
  CHARACTER(len=16), &
            PUBLIC, PROTECTED ::LINEARITY= 'linear'   ! To turn on non linear bracket term
  REAL(xp), PUBLIC, PROTECTED ::    mu_x =  0._xp     ! spatial    x-Hyperdiffusivity coefficient (for num. stability)
  REAL(xp), PUBLIC, PROTECTED ::    mu_y =  0._xp     ! spatial    y-Hyperdiffusivity coefficient (for num. stability)
  INTEGER,  PUBLIC, PROTECTED ::    N_HD =  4         ! order of numerical perpendicular spatial diffusion
  INTEGER,  PUBLIC, PROTECTED ::   N_HDz =  4         ! order of numerical parallel spatial diffusion
  LOGICAL,  PUBLIC, PROTECTED ::   HDz_h =  .false.    ! to apply z-hyperdiffusion on non adiab part
  REAL(xp), PUBLIC, PROTECTED ::    mu_z =  0._xp     ! spatial    z-Hyperdiffusivity coefficient (for num. stability)
  REAL(xp), PUBLIC, PROTECTED ::    mu_p =  0._xp     ! kinetic para hyperdiffusivity coefficient (for num. stability)
  REAL(xp), PUBLIC, PROTECTED ::    mu_j =  0._xp     ! kinetic perp hyperdiffusivity coefficient (for num. stability)
  CHARACTER(len=16), &
  PUBLIC, PROTECTED ::   HYP_V = 'hypcoll'  ! hyperdiffusion model for velocity space ('none','hypcoll','dvpar4')
  INTEGER,  PUBLIC, PROTECTED ::      Na =  1         ! number of evolved species
  LOGICAL,  PUBLIC            :: ADIAB_E =  .false.   ! adiabatic electron model
  LOGICAL,  PUBLIC            :: ADIAB_I =  .false.   ! adiabatic ion model
  REAL(xp), PUBLIC, PROTECTED ::   tau_i =  1.0       ! electron-ion temperature ratio for ion adiabatic model
  REAL(xp), PUBLIC, PROTECTED ::     q_i =  1.0       ! ion charge for ion adiabatic model
  REAL(xp), PUBLIC, PROTECTED ::      nu =  0._xp     ! collision frequency parameter
  REAL(xp), PUBLIC, PROTECTED ::    k_gB =  1._xp     ! artificial magnetic gradient tuner  (L_ref/L_gB)
  REAL(xp), PUBLIC, PROTECTED ::    k_cB =  1._xp     ! artificial magnetic curvature tuner (L_ref/L_cB)
  REAL(xp), PUBLIC, PROTECTED ::    k_mB =  1._xp     ! artificial mirror force tuner       (L_ref/L_cB)
  REAL(xp), PUBLIC, PROTECTED ::    k_tB =  1._xp     ! artificial trapping term tuner      (L_ref/L_cB)
  REAL(xp), PUBLIC, PROTECTED ::   k_ldB =  1._xp     ! artificial Landau damping tuner     (L_ref/L_cB)
  REAL(xp), PUBLIC, PROTECTED ::    beta =  0._xp     ! electron plasma Beta (8piNT_e/B0^2)
  REAL(xp), PUBLIC, PROTECTED :: ExBrate =  0._xp     ! ExB background shearing rate (radially constant shear flow)
  INTEGER,  PUBLIC, PROTECTED ::   ikxZF =  0         ! Background zonal mode wavenumber (acts in the nonlinear term)
  REAL(xp), PUBLIC, PROTECTED ::  ZFrate =  0._xp     ! Shearing rate of the background zonal mode
  LOGICAL,  PUBLIC, PROTECTED :: ZF_ONLY =  .false.   ! Cancels the nonlinear term excepts the background ZF contribution
  ! Auxiliary variable
  LOGICAL,  PUBLIC, PROTECTED ::      EM =  .true.    ! Electromagnetic effects flag
  LOGICAL,  PUBLIC, PROTECTED ::  MHD_PD =  .true.    ! MHD pressure drift
  ! Removes Landau damping in temperature and higher equation (Ivanov 2022)
  LOGICAL,  PUBLIC, PROTECTED :: RM_LD_T_EQ = .false.
  ! Flag to force the reality condition symmetry for the kx at ky=0
  LOGICAL,  PUBLIC, PROTECTED :: FORCE_SYMMETRY = .false.
  ! Add or remove the ExB nonlinear correction (Mcmillan 2019)
  LOGICAL,  PUBLIC, PROTECTED :: ExB_NL_CORRECTION = .true.
  CHARACTER(len=16), &
  PUBLIC, PROTECTED           :: KN_MODEL  =  'std'   ! Kernel model
  INTEGER,  PUBLIC, PROTECTED :: ORDER     = 1        ! order for Taylor expansion
  INTEGER,  PUBLIC, PROTECTED :: ORDER_NUM = 2        ! numerator order for Pade approx
  INTEGER,  PUBLIC, PROTECTED :: ORDER_DEN = 4        ! denominator order for Pade approx

  ! Module's routines
  PUBLIC :: model_readinputs, model_outputinputs

CONTAINS

  SUBROUTINE model_readinputs
    !    Read the input parameters
    USE basic,          ONLY: lu_in, speak
    USE parallel,       ONLY: num_procs_p
    USE prec_const,     ONLY: xp
    IMPLICIT NONE

    NAMELIST /MODEL/     LINEARITY, RM_LD_T_EQ, FORCE_SYMMETRY, MHD_PD, &
                         Na, ADIAB_E, ADIAB_I, q_i, tau_i, &
                         mu_x, mu_y, N_HD, N_HDz, HDz_h, mu_z, mu_p, mu_j, HYP_V, &
                         nu, k_gB, k_cB, k_mB, k_tB, k_ldB, &
                         beta, ExBrate, ExB_NL_CORRECTION,&
                         ikxZF, ZFrate, ZF_ONLY, KN_MODEL, ORDER, ORDER_NUM, ORDER_DEN

    READ(lu_in,model)

    IF (ADIAB_E .AND. ADIAB_I) &
      ERROR STOP '>> ERROR << cannot have both adiab e and adiab i models'

    IF((HYP_V .EQ. 'dvpar4') .AND. (num_procs_p .GT. 1)) &
      ERROR STOP '>> ERROR << dvpar4 velocity dissipation is not compatible with current p parallelization'

    IF(Na .EQ. 1) THEN
      IF((.NOT. ADIAB_E) .AND. (.NOT. ADIAB_I)) ERROR STOP "With one species, ADIAB_E or ADIAB_I must be set to .true. STOP"
      IF(ADIAB_E) THEN
        CALL speak('Adiabatic electron model -> beta = 0',1)
        beta = 0._xp
        EM   = .false.
      ENDIF
      IF(ADIAB_I) CALL speak('Adiabatic ion model',1)
    ENDIF

    IF(beta .GT. 0) THEN
      CALL speak('Electromagnetic effects are included',1)
      EM = .TRUE.
    ELSE
      EM = .FALSE.
    ENDIF

    IF(abs(ExBrate) .LT. epsilon(ExBrate))&
      ExB_NL_CORRECTION = .false.

  END SUBROUTINE model_readinputs

  SUBROUTINE model_outputinputs
    ! Write the input parameters to the results_xx.h5 file
    USE h5fortran
    IMPLICIT NONE
    CALL h5write("outputinput.h5", "/model/LINEARITY", LINEARITY)
    CALL h5write("outputinput.h5", "/model/RM_LD_T_EQ", RM_LD_T_EQ)
    CALL h5write("outputinput.h5", "/model/mu_x", mu_x)
    CALL h5write("outputinput.h5", "/model/mu_y", mu_y)
    CALL h5write("outputinput.h5", "/model/N_HD", N_HD)
    CALL h5write("outputinput.h5", "/model/mu_z", mu_z)
    CALL h5write("outputinput.h5", "/model/HDz_h", HDz_h)
    CALL h5write("outputinput.h5", "/model/mu_p", mu_p)
    CALL h5write("outputinput.h5", "/model/mu_j", mu_j)
    CALL h5write("outputinput.h5", "/model/HYP_V", HYP_V)
    CALL h5write("outputinput.h5", "/model/Na", Na)
    CALL h5write("outputinput.h5", "/model/nu", nu)
    CALL h5write("outputinput.h5", "/model/k_gB", k_gB)
    CALL h5write("outputinput.h5", "/model/k_cB", k_cB)
    CALL h5write("outputinput.h5", "/model/MHD_PD", MHD_PD)
    CALL h5write("outputinput.h5", "/model/beta", beta)
    CALL h5write("outputinput.h5", "/model/ExBrate", ExBrate)
    CALL h5write("outputinput.h5", "/model/ikxZF", ikxZF)
    CALL h5write("outputinput.h5", "/model/ZFrate", ZFrate)
    CALL h5write("outputinput.h5", "/model/ADIAB_E", ADIAB_E)
    CALL h5write("outputinput.h5", "/model/ADIAB_I", ADIAB_I)
    CALL h5write("outputinput.h5", "/model/tau_i", tau_i)
    CALL h5write("outputinput.h5", "/model/kern model", KN_MODEL)
    CALL h5write("outputinput.h5", "/model/order den", ORDER_DEN)
    CALL h5write("outputinput.h5", "/model/order num", ORDER_NUM)
    CALL h5write("outputinput.h5", "/model/order", ORDER)
  END SUBROUTINE model_outputinputs

END MODULE model
