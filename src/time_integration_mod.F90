MODULE time_integration

  USE prec_const
  IMPLICIT NONE
  PRIVATE

  INTEGER, PUBLIC, PROTECTED :: ntimelevel=4 ! Total number of time levels required by the numerical scheme
  INTEGER, PUBLIC, PROTECTED :: updatetlevel ! Current time level to be updated

  real(dp),PUBLIC,PROTECTED,DIMENSION(:,:),ALLOCATABLE :: A_E,A_I
  real(dp),PUBLIC,PROTECTED,DIMENSION(:),ALLOCATABLE :: b_E,b_Es,b_I
  real(dp),PUBLIC,PROTECTED,DIMENSION(:),ALLOCATABLE :: c_E,c_I !Coeff for Expl/Implic time integration in case of time dependent RHS (i.e. dy/dt = f(y,t)) see Baptiste Frei CSE Rapport 06/17

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


  SUBROUTINE time_integration_outputinputs(fidres, str)
    ! Write the input parameters to the results_xx.h5 file

    USE prec_const
    USE futils, ONLY: attach
    IMPLICIT NONE
    INTEGER, INTENT(in) :: fidres
    CHARACTER(len=256), INTENT(in) :: str

    CALL attach(fidres, TRIM(str), "numerical_scheme", numerical_scheme)

  END SUBROUTINE time_integration_outputinputs



  SUBROUTINE set_numerical_scheme
    ! Initialize Butcher coefficient of numerical_scheme

    use basic
    IMPLICIT NONE

    SELECT CASE (numerical_scheme)
    CASE ('RK4')
      CALL RK4
    CASE ('DOPRI5')
       CALL DOPRI5
    CASE ('DOPRI5_ADAPT')
      CALL DOPRI5_ADAPT
    CASE DEFAULT
       IF (my_id .EQ. 0) WRITE(*,*) 'Cannot initialize time integration scheme. Name invalid.'
    END SELECT

    IF (my_id .EQ. 0) WRITE(*,*) " Time integration with ", numerical_scheme

  END SUBROUTINE set_numerical_scheme

  SUBROUTINE RK4
    ! Butcher coeff for RK4 (default)

    USE basic
    USE prec_const
    IMPLICIT NONE
    INTEGER,PARAMETER :: nbstep = 4

    CALL allocate_array_dp1(c_E,1,nbstep)
    CALL allocate_array_dp1(b_E,1,nbstep)
    CALL allocate_array_dp2(A_E,1,nbstep,1,nbstep)

    ntimelevel = 4

    c_E(1) = 0.0_dp
    c_E(2) = 1.0_dp/2.0_dp
    c_E(3) = 1.0_dp/2.0_dp
    c_E(4) = 1.0_dp

    b_E(1) = 1.0_dp/6.0_dp
    b_E(2) = 1.0_dp/3.0_dp
    b_E(3) = 1.0_dp/3.0_dp
    b_E(4) = 1.0_dp/6.0_dp

    A_E(2,1) = 1.0_dp/2.0_dp
    A_E(3,2) = 1.0_dp/2.0_dp
    A_E(4,3) = 1.0_dp

  END SUBROUTINE RK4



  SUBROUTINE DOPRI5
    ! Butcher coeff for DOPRI5 --> Stiffness detection
    ! DOPRI5 used for stiffness detection.
    ! 5 order method/7 stages

    USE basic
    IMPLICIT NONE
    INTEGER,PARAMETER :: nbstep =7

    CALL allocate_array_dp1(c_E,1,nbstep)
    CALL allocate_array_dp1(b_E,1,nbstep)
    CALL allocate_array_dp2(A_E,1,nbstep,1,nbstep)

    ntimelevel = 7
    !
    c_E(1) = 0._dp
    c_E(2) = 1.0_dp/5.0_dp
    c_E(3) = 3.0_dp /10.0_dp
    c_E(4) = 4.0_dp/5.0_dp
    c_E(5) = 8.0_dp/9.0_dp
    c_E(6) = 1.0_dp
    c_E(7) = 1.0_dp
    !
    A_E(2,1) = 1.0_dp/5.0_dp
    A_E(3,1) = 3.0_dp/40.0_dp
    A_E(3,2) = 9.0_dp/40.0_dp
    A_E(4,1) = 	44.0_dp/45.0_dp
    A_E(4,2) = -56.0_dp/15.0_dp
    A_E(4,3) = 32.0_dp/9.0_dp
    A_E(5,1 ) = 19372.0_dp/6561.0_dp
    A_E(5,2) = -25360.0_dp/2187.0_dp
    A_E(5,3) = 64448.0_dp/6561.0_dp
    A_E(5,4) = -212.0_dp/729.0_dp
    A_E(6,1) = 9017.0_dp/3168.0_dp
    A_E(6,2)= -355.0_dp/33.0_dp
    A_E(6,3) = 46732.0_dp/5247.0_dp
    A_E(6,4) = 49.0_dp/176.0_dp
    A_E(6,5) = -5103.0_dp/18656.0_dp
    A_E(7,1) = 35.0_dp/384.0_dp
    A_E(7,3) = 500.0_dp/1113.0_dp
    A_E(7,4) = 125.0_dp/192.0_dp
    A_E(7,5) = -2187.0_dp/6784.0_dp
    A_E(7,6) = 11.0_dp/84.0_dp
    !
    b_E(1) = 35.0_dp/384.0_dp
    b_E(2) = 0._dp
    b_E(3) = 500.0_dp/1113.0_dp
    b_E(4) = 125.0_dp/192.0_dp
    b_E(5) = -2187.0_dp/6784.0_dp
    b_E(6) = 11.0_dp/84.0_dp
    b_E(7) = 0._dp
    !
  END SUBROUTINE DOPRI5

  SUBROUTINE DOPRI5_ADAPT
    ! Butcher coeff for DOPRI5 --> Stiffness detection
    ! DOPRI5 used for stiffness detection.
    ! 5 order method/7 stages

    USE basic
    IMPLICIT NONE
    INTEGER,PARAMETER :: nbstep =7

    CALL allocate_array_dp1(c_E,1,nbstep)
    CALL allocate_array_dp1(b_E,1,nbstep)
    CALL allocate_array_dp1(b_Es,1,nbstep)
    CALL allocate_array_dp2(A_E,1,nbstep,1,nbstep)

    ntimelevel = 7
    !
    c_E(1) = 0._dp
    c_E(2) = 1.0_dp/5.0_dp
    c_E(3) = 3.0_dp /10.0_dp
    c_E(4) = 4.0_dp/5.0_dp
    c_E(5) = 8.0_dp/9.0_dp
    c_E(6) = 1.0_dp
    c_E(7) = 1.0_dp
    !
    A_E(2,1) = 1.0_dp/5.0_dp
    A_E(3,1) = 3.0_dp/40.0_dp
    A_E(3,2) = 9.0_dp/40.0_dp
    A_E(4,1) = 	44.0_dp/45.0_dp
    A_E(4,2) = -56.0_dp/15.0_dp
    A_E(4,3) = 32.0_dp/9.0_dp
    A_E(5,1 ) = 19372.0_dp/6561.0_dp
    A_E(5,2) = -25360.0_dp/2187.0_dp
    A_E(5,3) = 64448.0_dp/6561.0_dp
    A_E(5,4) = -212.0_dp/729.0_dp
    A_E(6,1) = 9017.0_dp/3168.0_dp
    A_E(6,2)= -355.0_dp/33.0_dp
    A_E(6,3) = 46732.0_dp/5247.0_dp
    A_E(6,4) = 49.0_dp/176.0_dp
    A_E(6,5) = -5103.0_dp/18656.0_dp
    A_E(7,1) = 35.0_dp/384.0_dp
    A_E(7,3) = 500.0_dp/1113.0_dp
    A_E(7,4) = 125.0_dp/192.0_dp
    A_E(7,5) = -2187.0_dp/6784.0_dp
    A_E(7,6) = 11.0_dp/84.0_dp
    !
    b_E(1) = 35.0_dp/384.0_dp
    b_E(2) = 0._dp
    b_E(3) = 500.0_dp/1113.0_dp
    b_E(4) = 125.0_dp/192.0_dp
    b_E(5) = -2187.0_dp/6784.0_dp
    b_E(6) = 11.0_dp/84.0_dp
    b_E(7) = 0._dp
    !
    b_Es(1) = 5179.0_dp/57600.0_dp
    b_Es(2) = 0._dp
    b_Es(3) = 7571.0_dp/16695.0_dp
    b_Es(4) = 393.0_dp/640.0_dp
    b_Es(5) = -92097.0_dp/339200.0_dp
    b_Es(6) = 187.0_dp/2100.0_dp
    b_Es(7) = 1._dp/40._dp
    !
  END SUBROUTINE DOPRI5_ADAPT

END MODULE time_integration
