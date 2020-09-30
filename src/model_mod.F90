MODULE model
  ! Module for diagnostic parameters

  USE prec_const
  IMPLICIT NONE
  PRIVATE

  INTEGER,  PUBLIC, PROTECTED ::      CO =  0         ! Collision Operator
  LOGICAL,  PUBLIC, PROTECTED :: NON_LIN =  .true.    ! To turn on non linear bracket term
  REAL(dp), PUBLIC, PROTECTED ::      nu =  1._dp     ! Collision frequency
  REAL(dp), PUBLIC, PROTECTED ::   tau_e =  1._dp     ! Temperature
  REAL(dp), PUBLIC, PROTECTED ::   tau_i =  1._dp     !
  REAL(dp), PUBLIC, PROTECTED :: sigma_e =  1._dp     ! Mass
  REAL(dp), PUBLIC, PROTECTED :: sigma_i =  1._dp     !
  REAL(dp), PUBLIC, PROTECTED ::     q_e = -1._dp     ! Charge
  REAL(dp), PUBLIC, PROTECTED ::     q_i =  1._dp     !
  REAL(dp), PUBLIC, PROTECTED ::   eta_n =  1._dp     ! Density gradient
  REAL(dp), PUBLIC, PROTECTED ::   eta_T =  1._dp     ! Temperature gradient
  REAL(dp), PUBLIC, PROTECTED ::   eta_B =  1._dp     ! Magnetic gradient
  REAL(dp), PUBLIC, PROTECTED :: lambdaD =  1._dp     ! Debye length

  PUBLIC :: model_readinputs, model_outputinputs

CONTAINS

  SUBROUTINE model_readinputs
    !    Read the input parameters

    USE basic, ONLY : lu_in
    USE prec_const
    IMPLICIT NONE

    NAMELIST /MODEL_PAR/ CO, NON_LIN, nu, tau_e, tau_i, sigma_e, sigma_i, &
                         q_e, q_i, eta_n, eta_T, eta_B, lambdaD

    READ(lu_in,model_par)
    !WRITE(*,model_par)

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
    CALL attach(fidres, TRIM(str),      "CO",      CO)
    CALL attach(fidres, TRIM(str), "NON_LIN", NON_LIN)
    CALL attach(fidres, TRIM(str),      "nu",      nu)
    CALL attach(fidres, TRIM(str),   "tau_e",   tau_e)
    CALL attach(fidres, TRIM(str),   "tau_i",   tau_i)
    CALL attach(fidres, TRIM(str), "sigma_e", sigma_e)
    CALL attach(fidres, TRIM(str), "sigma_i", sigma_i)
    CALL attach(fidres, TRIM(str),     "q_e",     q_e)
    CALL attach(fidres, TRIM(str),     "q_i",     q_i)
    CALL attach(fidres, TRIM(str),   "eta_n",   eta_n)
    CALL attach(fidres, TRIM(str),   "eta_T",   eta_T)
    CALL attach(fidres, TRIM(str),   "eta_B",   eta_B)
    CALL attach(fidres, TRIM(str), "lambdaD", lambdaD)
  END SUBROUTINE model_outputinputs

END MODULE model
