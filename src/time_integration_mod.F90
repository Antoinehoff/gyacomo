MODULE time_integration

  USE prec_const
  IMPLICIT NONE
  PRIVATE

  INTEGER, PUBLIC, PROTECTED :: ntimelevel=4 ! Total number of time levels required by the numerical scheme
  INTEGER, PUBLIC, PROTECTED :: updatetlevel ! Current time level to be updated

  real(xp),PUBLIC,PROTECTED,DIMENSION(:,:),ALLOCATABLE :: A_E,A_I
  real(xp),PUBLIC,PROTECTED,DIMENSION(:),ALLOCATABLE :: b_E,b_Es,b_I
  real(xp),PUBLIC,PROTECTED,DIMENSION(:),ALLOCATABLE :: c_E,c_I !Coeff for Expl/Implic time integration in case of time dependent RHS (i.e. dy/dt = f(y,t)) see Baptiste Frei CSE Rapport 06/17

  character(len=10),PUBLIC,PROTECTED :: numerical_scheme='RK4'

  PUBLIC :: set_updatetlevel, time_integration_readinputs, time_integration_outputinputs

CONTAINS

  SUBROUTINE set_updatetlevel(new_updatetlevel)
    INTEGER, INTENT(in) :: new_updatetlevel
    updatetlevel = new_updatetlevel
  END SUBROUTINE set_updatetlevel

  SUBROUTINE time_integration_readinputs
    ! Read the input parameters
    USE prec_const
    USE basic, ONLY : lu_in
    IMPLICIT NONE
    NAMELIST /TIME_INTEGRATION_PAR/ numerical_scheme
    READ(lu_in,time_integration_par)
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
    IMPLICIT NONE
    SELECT CASE (numerical_scheme)
    ! Order 2 methods
    CASE ('RK2')
      CALL RK2
    CASE ('SSPx_RK2')
      CALL SSPx_RK2
    ! Order 3 methods
    CASE ('RK3')
      CALL RK3
    CASE ('SSP_RK3')
      CALL SSP_RK3
    CASE ('SSPx_RK3')
      CALL SSPx_RK3
    CASE ('IMEX_SSP2')
      CALL IMEX_SSP2
    CASE ('ARK2')
      CALL ARK2
    ! Order 4 methods
    CASE ('RK4')
      CALL RK4
    ! Order 5 methods
    CASE ('DOPRI5')
      CALL DOPRI5
    CASE DEFAULT
       IF (my_id .EQ. 0) WRITE(*,*) 'Cannot initialize time integration scheme. Name invalid.'
    END SELECT
    IF (my_id .EQ. 0) WRITE(*,*) " Time integration with ", numerical_scheme
  END SUBROUTINE set_numerical_scheme

  !!! second order time schemes
  SUBROUTINE RK2
    ! Butcher coeff for clasical RK2 (Heun's)
    USE basic
    USE prec_const
    IMPLICIT NONE
    INTEGER,PARAMETER :: nbstep = 2
    CALL allocate_array(c_E,1,nbstep)
    CALL allocate_array(b_E,1,nbstep)
    CALL allocate_array(A_E,1,nbstep,1,nbstep)
    ntimelevel = 2
    c_E(1) = 0.0_xp
    c_E(2) = 1.0_xp
    b_E(1) = 1._xp/2._xp
    b_E(2) = 1._xp/2._xp
    A_E(2,1) = 1._xp
  END SUBROUTINE RK2

  SUBROUTINE SSPx_RK2
    ! DOESNT WORK
    ! Butcher coeff for modified strong stability  preserving RK2
    ! used in GX (Hammett 2022, Mandell et al. 2022)
    USE basic
    USE prec_const
    IMPLICIT NONE
    INTEGER,PARAMETER :: nbstep = 2
    REAL(xp) :: alpha, beta
    alpha = 1._xp/SQRT(2._xp)
    beta  = SQRT(2._xp) - 1._xp
    CALL allocate_array(c_E,1,nbstep)
    CALL allocate_array(b_E,1,nbstep)
    CALL allocate_array(A_E,1,nbstep,1,nbstep)
    ntimelevel = 2
    c_E(1)   = 0.0_xp
    c_E(2)   = 1.0_xp/2.0_xp
    b_E(1)   = alpha*beta/2._xp
    b_E(2)   = alpha/2._xp
    A_E(2,1) = alpha
    ! b_E(1) = 1._xp
    ! b_E(2) = 1._xp/SQRT(2._xp)
    ! A_E(2,1) = 1._xp/SQRT(2._xp)
  END SUBROUTINE SSPx_RK2

  !!! third order time schemes
  SUBROUTINE RK3
    ! Butcher coeff for classical RK3
    USE basic
    USE prec_const
    IMPLICIT NONE
    INTEGER,PARAMETER :: nbstep = 3
    CALL allocate_array(c_E,1,nbstep)
    CALL allocate_array(b_E,1,nbstep)
    CALL allocate_array(A_E,1,nbstep,1,nbstep)
    ntimelevel = 3
    c_E(1)   = 0.0_xp
    c_E(2)   = 1.0_xp/2.0_xp
    c_E(3)   = 1.0_xp
    b_E(1)   = 1._xp/6._xp
    b_E(2)   = 2._xp/3._xp
    b_E(3)   = 1._xp/6._xp
    A_E(2,1) = 1.0_xp/2.0_xp
    A_E(3,1) = -1._xp
    A_E(3,2) = 2._xp
  END SUBROUTINE RK3

  SUBROUTINE SSPx_RK3
    ! DOESNT WORK
    ! Butcher coeff for modified strong stability  preserving RK3
    ! used in GX (Hammett 2022, Mandell et al. 2022)
    USE basic
    USE prec_const
    IMPLICIT NONE
    INTEGER,PARAMETER :: nbstep = 3
    REAL(xp) :: a1, a2, a3, w1, w2, w3
    a1 = (1._xp/6._xp)**(1._xp/3._xp)! (1/6)^(1/3)
    ! a1 = 0.5503212081491044571635029569733887910843_xp ! (1/6)^(1/3)
    a2 = a1
    a3 = a1
    w1 = 0.5_xp*(-1._xp + SQRT( 9._xp - 2._xp * 6._xp**(2._xp/3._xp))) ! (-1 + sqrt(9-2*6^(2/3)))/2
    ! w1 = 0.2739744023885328783052273138309828937054_xp ! (sqrt(9-2*6^(2/3))-1)/2
    w2 = 0.5_xp*(-1._xp + 6._xp**(2._xp/3._xp) - SQRT(9._xp - 2._xp * 6._xp**(2._xp/3._xp))) ! (6^(2/3)-1-sqrt(9-2*6^(2/3)))/2
    ! w2 = 0.3769892220587804931852815570891834795475_xp ! (6^(2/3)-1-sqrt(9-2*6^(2/3)))/2
    w3 = 1._xp/a1 - w2 * (1._xp + w1)
    ! w3 = 1.3368459739528868457369981115334667265415_xp
    CALL allocate_array(c_E,1,nbstep)
    CALL allocate_array(b_E,1,nbstep)
    CALL allocate_array(A_E,1,nbstep,1,nbstep)
    ntimelevel = 3
    c_E(1) = 0.0_xp
    c_E(2) = 1.0_xp/2.0_xp
    c_E(3) = 1.0_xp/2.0_xp
    b_E(1) = a1 * (w1*w2 + w3)
    b_E(2) = a2 * w2
    b_E(3) = a3
    A_E(2,1) = a1
    A_E(3,1) = a1 * w1
    A_E(3,2) = a2
  END SUBROUTINE SSPx_RK3

  SUBROUTINE IMEX_SSP2
    !! Version of Rokhzadi 2017 (An Optimally Stable and Accurate Second-Order
    !   SSP Runge-Kutta IMEX Scheme for Atmospheric Applications)
    USE basic
    USE prec_const
    IMPLICIT NONE
    INTEGER,PARAMETER :: nbstep = 3
    CALL allocate_array(c_E,1,nbstep)
    CALL allocate_array(b_E,1,nbstep)
    CALL allocate_array(A_E,1,nbstep,1,nbstep)
    ntimelevel = 3
    c_E(1)   = 0._xp
    c_E(2)   = 0.711664700366941_xp
    c_E(3)   = 0.994611536833690_xp
    b_E(1)   = 0.398930808264688_xp
    b_E(2)   = 0.345755244189623_xp
    b_E(3)   = 0.255313947545689_xp
    A_E(2,1) = 0.711664700366941_xp
    A_E(3,1) = 0.077338168947683_xp
    A_E(3,2) = 0.917273367886007_xp
  END SUBROUTINE IMEX_SSP2

  SUBROUTINE ARK2
    !! Version of Rokhzadi 2017 (An Optimally Stable and Accurate Second-Order
    !   SSP Runge-Kutta IMEX Scheme for Atmospheric Applications)
    USE basic
    USE prec_const
    IMPLICIT NONE
    INTEGER,PARAMETER :: nbstep = 3
    CALL allocate_array(c_E,1,nbstep)
    CALL allocate_array(b_E,1,nbstep)
    CALL allocate_array(A_E,1,nbstep,1,nbstep)
    ntimelevel = 3
    c_E(1)   = 0._xp
    c_E(2)   = 2._xp*(1._xp - 1._xp/SQRT2)
    c_E(3)   = 1._xp
    b_E(1)   = 1._xp/(2._xp*SQRT2)
    b_E(2)   = 1._xp/(2._xp*SQRT2)
    b_E(3)   = 1._xp - 1._xp/SQRT2
    A_E(2,1) = 2._xp*(1._xp - 1._xp/SQRT2)
    A_E(3,1) = 1._xp - (3._xp + 2._xp*SQRT2)/6._xp
    A_E(3,2) = (3._xp + 2._xp*SQRT2)/6._xp
  END SUBROUTINE ARK2

  SUBROUTINE SSP_RK3
    ! Butcher coeff for strong stability  preserving RK3
    USE basic
    USE prec_const
    IMPLICIT NONE
    INTEGER,PARAMETER :: nbstep = 3
    CALL allocate_array(c_E,1,nbstep)
    CALL allocate_array(b_E,1,nbstep)
    CALL allocate_array(A_E,1,nbstep,1,nbstep)
    ntimelevel = 3
    c_E(1)   = 0.0_xp
    c_E(2)   = 1.0_xp
    c_E(3)   = 1.0_xp/2.0_xp
    b_E(1)   = 1._xp/6._xp
    b_E(2)   = 1._xp/6._xp
    b_E(3)   = 2._xp/3._xp
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
    CALL allocate_array(b_E,1,nbstep)
    CALL allocate_array(A_E,1,nbstep,1,nbstep)
    ntimelevel = 4
    c_E(1)   = 0.0_xp
    c_E(2)   = 1.0_xp/2.0_xp
    c_E(3)   = 1.0_xp/2.0_xp
    c_E(4)   = 1.0_xp
    b_E(1)   = 1.0_xp/6.0_xp
    b_E(2)   = 1.0_xp/3.0_xp
    b_E(3)   = 1.0_xp/3.0_xp
    b_E(4)   = 1.0_xp/6.0_xp
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
    CALL allocate_array(b_E,1,nbstep)
    CALL allocate_array(A_E,1,nbstep,1,nbstep)
    ntimelevel = 7
    c_E(1) = 0._xp
    c_E(2) = 1.0_xp/5.0_xp
    c_E(3) = 3.0_xp /10.0_xp
    c_E(4) = 4.0_xp/5.0_xp
    c_E(5) = 8.0_xp/9.0_xp
    c_E(6) = 1.0_xp
    c_E(7) = 1.0_xp
    A_E(2,1) = 1.0_xp/5.0_xp
    A_E(3,1) = 3.0_xp/40.0_xp
    A_E(3,2) = 9.0_xp/40.0_xp
    A_E(4,1) = 	44.0_xp/45.0_xp
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
    b_E(1) = 35.0_xp/384.0_xp
    b_E(2) = 0._xp
    b_E(3) = 500.0_xp/1113.0_xp
    b_E(4) = 125.0_xp/192.0_xp
    b_E(5) = -2187.0_xp/6784.0_xp
    b_E(6) = 11.0_xp/84.0_xp
    b_E(7) = 0._xp
  END SUBROUTINE DOPRI5

END MODULE time_integration
