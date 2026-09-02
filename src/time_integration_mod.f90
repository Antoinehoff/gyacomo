MODULE time_integration

  USE prec_const

  IMPLICIT NONE
  PRIVATE

  INTEGER, PUBLIC, PROTECTED :: ntimelevel=4 ! Total number of time levels required by the numerical scheme
  INTEGER, PUBLIC, PROTECTED :: updatetlevel ! Current time level to be updated

  real(xp),PUBLIC,PROTECTED,DIMENSION(:,:),ALLOCATABLE :: A_E
  real(xp),PUBLIC,PROTECTED,DIMENSION(:,:),ALLOCATABLE :: b_E
  real(xp),PUBLIC,PROTECTED,DIMENSION(:),ALLOCATABLE :: c_E !Coeff for Expl/Implic time integration in case of time dependent RHS (i.e. dy/dt = f(y,t)) see Baptiste Frei CSE Rapport 06/17

  character(len=10),PUBLIC,PROTECTED :: numerical_scheme='RK4'


  PUBLIC :: set_updatetlevel, time_integration_readinputs, time_integration_outputinputs

CONTAINS

  SUBROUTINE time_integration_readinputs
    ! Read the input parameters
    USE prec_const
    USE basic, ONLY : lu_in
    IMPLICIT NONE
    NAMELIST /TIME_INTEGRATION/ numerical_scheme

    READ(lu_in,time_integration)
    CALL set_numerical_scheme
  END SUBROUTINE time_integration_readinputs


  SUBROUTINE time_integration_outputinputs(fid)
    ! Write the input parameters to the results_xx.h5 file
    USE futils, ONLY: attach, creatd
    IMPLICIT NONE
    INTEGER, INTENT(in) :: fid
    CHARACTER(len=256)  :: str
    WRITE(str,'(a)') '/data/input/time_integration'
    CALL creatd(fid, 0,(/0/),TRIM(str),'Time Integration Input')
    CALL attach(fid, TRIM(str), "numerical_scheme", numerical_scheme)
  END SUBROUTINE time_integration_outputinputs

  SUBROUTINE set_numerical_scheme
    ! Initialize Butcher coefficient of set_numerical_scheme
    use parallel, ONLY: my_id
    use basic,    ONLY: speak
    IMPLICIT NONE
    SELECT CASE (numerical_scheme)
    ! Order 1 method
    CASE ('EE')
      CALL EE
    ! Order 2 methods
    CASE ('RK2')
      CALL RK2
    ! Order 3 methods
    CASE ('RK3')
      CALL RK3
    CASE ('SSP_RK3')
      CALL SSP_RK3
    ! Embedded pairs. NOTE: these are integrated at the fixed dt of the BASIC
    ! namelist using the higher order weights b_E(:,1). The step size is not
    ! adapted: the embedded estimate b_E(:,2) is stored but not used.
    CASE ('ODE23')
      CALL ode23
    ! Order 4 methods
    CASE ('RK4')
      CALL RK4
    ! Order 5 methods
    CASE ('DOPRI5')
      CALL DOPRI5
    CASE DEFAULT
       ERROR STOP 'Cannot initialize time integration scheme. Name invalid.'
    END SELECT
    CALL speak("-Time integration with "// numerical_scheme,2)
  END SUBROUTINE set_numerical_scheme

  !!! first order time schemes
  SUBROUTINE EE
    ! Butcher coeff for Euler Explicit scheme
    USE basic
    USE prec_const
    IMPLICIT NONE
    INTEGER,PARAMETER :: nbstep = 1
    CALL allocate_array(c_E,1,nbstep)
    CALL allocate_array(b_E,1,nbstep,1,1)
    CALL allocate_array(A_E,1,nbstep,1,nbstep)
    ntimelevel = 1
    c_E(1)   = 0._xp
    b_E(1,1) = 1._xp
    A_E(1,1) = 0._xp
  END SUBROUTINE EE

  !!! second order time schemes
  SUBROUTINE RK2
    ! Butcher coeff for clasical RK2 (Heun's)
    USE basic
    USE prec_const
    IMPLICIT NONE
    INTEGER,PARAMETER :: nbstep = 2
    CALL allocate_array(c_E,1,nbstep)
    CALL allocate_array(b_E,1,nbstep,1,1)
    CALL allocate_array(A_E,1,nbstep,1,nbstep)
    ntimelevel = 2
    c_E(1)   = 0.0_xp
    c_E(2)   = 1.0_xp
    b_E(1,1) = 1._xp/2._xp
    b_E(2,1) = 1._xp/2._xp
    A_E(2,1) = 1._xp
  END SUBROUTINE RK2

  !!! third order time schemes
  SUBROUTINE RK3
    ! Butcher coeff for classical RK3
    USE basic
    USE prec_const
    IMPLICIT NONE
    INTEGER,PARAMETER :: nbstep = 3
    CALL allocate_array(c_E,1,nbstep)
    CALL allocate_array(b_E,1,nbstep,1,1)
    CALL allocate_array(A_E,1,nbstep,1,nbstep)
    ntimelevel = 3
    c_E(1)   = 0.0_xp
    c_E(2)   = 1.0_xp/2.0_xp
    c_E(3)   = 1.0_xp
    b_E(1,1) = 1._xp/6._xp
    b_E(2,1) = 2._xp/3._xp
    b_E(3,1) = 1._xp/6._xp
    A_E(2,1) = 1.0_xp/2.0_xp
    A_E(3,1) = -1._xp
    A_E(3,2) = 2._xp
  END SUBROUTINE RK3

  SUBROUTINE SSP_RK3
    ! Butcher coeff for strong stability  preserving RK3
    USE basic
    USE prec_const
    IMPLICIT NONE
    INTEGER,PARAMETER :: nbstep = 3
    CALL allocate_array(c_E,1,nbstep)
    CALL allocate_array(b_E,1,nbstep,1,1)
    CALL allocate_array(A_E,1,nbstep,1,nbstep)
    ntimelevel = 3
    c_E(1)   = 0.0_xp
    c_E(2)   = 1.0_xp
    c_E(3)   = 1.0_xp/2.0_xp
    b_E(1,1) = 1._xp/6._xp
    b_E(2,1) = 1._xp/6._xp
    b_E(3,1) = 2._xp/3._xp
    A_E(2,1) = 1._xp
    A_E(3,1) = 1._xp/4._xp
    A_E(3,2) = 1._xp/4._xp
  END SUBROUTINE SSP_RK3

  !!! fourth order time schemes
  SUBROUTINE RK4
    ! Butcher coeff for RK4 (default)
    USE basic
    USE prec_const
    IMPLICIT NONE
    INTEGER,PARAMETER :: nbstep = 4
    CALL allocate_array(c_E,1,nbstep)
    CALL allocate_array(b_E,1,nbstep,1,1)
    CALL allocate_array(A_E,1,nbstep,1,nbstep)
    ntimelevel = 4
    c_E(1)   = 0.0_xp
    c_E(2)   = 1.0_xp/2.0_xp
    c_E(3)   = 1.0_xp/2.0_xp
    c_E(4)   = 1.0_xp
    b_E(1,1) = 1.0_xp/6.0_xp
    b_E(2,1) = 1.0_xp/3.0_xp
    b_E(3,1) = 1.0_xp/3.0_xp
    b_E(4,1) = 1.0_xp/6.0_xp
    A_E(2,1) = 1.0_xp/2.0_xp
    A_E(3,2) = 1.0_xp/2.0_xp
    A_E(4,3) = 1.0_xp
  END SUBROUTINE RK4

  !!! fifth order time schemes
  SUBROUTINE DOPRI5
    ! Butcher coeff for DOPRI5 --> Stiffness detection
    ! DOPRI5 used for stiffness detection.
    ! 5 order method/7 stages
    USE basic
    IMPLICIT NONE
    INTEGER,PARAMETER :: nbstep =7
    CALL allocate_array(c_E,1,nbstep)
    CALL allocate_array(b_E,1,nbstep,1,2)
    CALL allocate_array(A_E,1,nbstep,1,nbstep)

    ntimelevel = 7

    c_E = 0._xp
    c_E(2) = 1.0_xp/5.0_xp
    c_E(3) = 3.0_xp /10.0_xp
    c_E(4) = 4.0_xp/5.0_xp
    c_E(5) = 8.0_xp/9.0_xp
    c_E(6) = 1.0_xp
    c_E(7) = 1.0_xp

    A_E = 0._xp
    A_E(2,1) = 1.0_xp/5.0_xp
    A_E(3,1) = 3.0_xp/40.0_xp
    A_E(3,2) = 9.0_xp/40.0_xp
    A_E(4,1) = 44.0_xp/45.0_xp
    A_E(4,2) = -56.0_xp/15.0_xp
    A_E(4,3) = 32.0_xp/9.0_xp
    A_E(5,1 ) = 19372.0_xp/6561.0_xp
    A_E(5,2) = -25360.0_xp/2187.0_xp
    A_E(5,3) = 64448.0_xp/6561.0_xp
    A_E(5,4) = -212.0_xp/729.0_xp
    A_E(6,1) = 9017.0_xp/3168.0_xp
    A_E(6,2)= -355.0_xp/33.0_xp
    A_E(6,3) = 46732.0_xp/5247.0_xp
    A_E(6,4) = 49.0_xp/176.0_xp
    A_E(6,5) = -5103.0_xp/18656.0_xp
    A_E(7,1) = 35.0_xp/384.0_xp
    A_E(7,3) = 500.0_xp/1113.0_xp
    A_E(7,4) = 125.0_xp/192.0_xp
    A_E(7,5) = -2187.0_xp/6784.0_xp
    A_E(7,6) = 11.0_xp/84.0_xp

    b_E = 0._xp
    b_E(1, 1) = 35._xp/384._xp
    b_E(2, 1) = 0._xp
    b_E(3, 1) = 500._xp/1113._xp
    b_E(4, 1) = 125._xp/192._xp
    b_E(5, 1) = -2187._xp/6784._xp
    b_E(6, 1) = 11._xp/84._xp
    b_E(7, 1) = 0._xp
    b_E(1, 2) = 5179._xp/57600._xp
    b_E(2, 2) = 0._xp
    b_E(3, 2) = 7571._xp/16695._xp
    b_E(4, 2) = 393._xp/640._xp
    b_E(5, 2) = -92097._xp/339200._xp
    b_E(6, 2) = 187._xp/2100._xp
    b_E(7, 2) = 1._xp/40._xp
  END SUBROUTINE DOPRI5

  subroutine set_updatetlevel(new_updatetlevel)
  integer, intent(in) :: new_updatetlevel
  updatetlevel = new_updatetlevel
  end subroutine set_updatetlevel

  !> Sets up Butcher coefficients for the Bogacki-Shampine method (orders 3 and 2)
  subroutine ode23
    use basic
    integer,parameter :: nbstep = 4
    CALL allocate_array(c_E,1,nbstep)
    CALL allocate_array(b_E,1,nbstep, 1, 2)
    CALL allocate_array(A_E,1,nbstep,1,nbstep)

    ntimelevel = 4

    c_E = 0._xp
    c_E(2) = 1.0_xp/2.0_xp
    c_E(3) = 3.0_xp /4.0_xp
    c_E(4) = 1._xp

    A_E = 0._xp
    A_E(2,1) = 1.0_xp/2.0_xp
    A_E(3,2) = 3.0_xp/4.0_xp
    A_E(4,1) = 2.0_xp/9.0_xp
    A_E(4,2) = 1.0_xp/3.0_xp
    A_E(4,3) = 4.0_xp/9.0_xp

    b_E = 0._xp
    b_E(1, 1) = 2._xp/9._xp
    b_E(2, 1) = 1._xp/3._xp
    b_E(3, 1) = 4._xp/9._xp
    b_E(4, 1) = 0._xp
    b_E(1, 2) = 7._xp/24._xp
    b_E(2, 2) = 1._xp/4._xp
    b_E(3, 2) = 1._xp/3._xp
    b_E(4, 2) = 1._xp/8._xp
  end subroutine ode23
  !!-------------------------------------------------------------------------
  !!-------------------------------------------------------------------------

END MODULE time_integration
