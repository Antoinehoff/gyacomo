PROGRAM test_time_integration
  ! Verifies the Butcher tableaux built by time_integration_mod.
  !
  ! A mistyped digit in a tableau does not crash anything: the run simply
  ! integrates at a lower order than advertised, which looks like a physics
  ! result. The tableaux are therefore checked against the algebraic order
  ! conditions of explicit Runge-Kutta methods, so each scheme has to prove the
  ! order its name claims.
  !
  ! The scheme is taken from the command line because allocate_array does not
  ! deallocate first: time_integration_readinputs can only be called once per
  ! process. "make check" runs this binary once per scheme.
  !
  !   usage: test_time_integration [EE|RK2|RK3|SSP_RK3|ODE23|RK4|DOPRI5]
  USE prec_const,       ONLY: xp
  USE basic,            ONLY: lu_in
  USE parallel,         ONLY: ppinit
  USE time_integration, ONLY: time_integration_readinputs, numerical_scheme, &
                              ntimelevel, A_E, b_E, c_E
  USE test_harness
  IMPLICIT NONE

  CHARACTER(len=32) :: scheme
  INTEGER :: nargs, s, i, j
  INTEGER :: want_stages, want_order, want_embedded_order
  LOGICAL :: explicit_ok, consistent_ok

  CALL ppinit
  nargs = COMMAND_ARGUMENT_COUNT()
  IF (nargs .GE. 1) THEN
    CALL GET_COMMAND_ARGUMENT(1, scheme)
  ELSE
    scheme = 'RK4'
  ENDIF

  CALL start_suite('time_integration ['//TRIM(scheme)//']')

  ! What the scheme is expected to deliver. want_embedded_order is 0 when the
  ! scheme has no embedded pair, i.e. b_E has a single column.
  SELECT CASE (TRIM(scheme))
  CASE ('EE');      want_stages = 1; want_order = 1; want_embedded_order = 0
  CASE ('RK2');     want_stages = 2; want_order = 2; want_embedded_order = 0
  CASE ('RK3');     want_stages = 3; want_order = 3; want_embedded_order = 0
  CASE ('SSP_RK3'); want_stages = 3; want_order = 3; want_embedded_order = 0
  CASE ('ODE23');   want_stages = 4; want_order = 3; want_embedded_order = 2
  CASE ('RK4');     want_stages = 4; want_order = 4; want_embedded_order = 0
  CASE ('DOPRI5');  want_stages = 7; want_order = 5; want_embedded_order = 4
  CASE DEFAULT
    WRITE(*,'(A)') 'unknown scheme "'//TRIM(scheme)//'"'
    CALL check('scheme is known to the test', .FALSE.)
    CALL finish_suite
    STOP
  END SELECT

  ! Feed the namelist through a scratch file so nothing is left on disk.
  OPEN(UNIT=lu_in, STATUS='SCRATCH', ACTION='READWRITE', FORM='FORMATTED')
  WRITE(lu_in,'(A)') '&TIME_INTEGRATION'
  WRITE(lu_in,'(A)') "  numerical_scheme = '"//TRIM(scheme)//"'"
  WRITE(lu_in,'(A)') '/'
  REWIND(lu_in)
  CALL time_integration_readinputs
  CLOSE(lu_in)

  CALL section('the requested scheme was selected')
  CALL check('numerical_scheme matches the input', TRIM(numerical_scheme) .EQ. TRIM(scheme))
  CALL check('coefficients allocated', ALLOCATED(A_E) .AND. ALLOCATED(b_E) .AND. ALLOCATED(c_E))

  s = SIZE(c_E)
  CALL section('tableau has the expected shape')
  CALL check_int('number of stages',      s,           want_stages)
  CALL check_int('ntimelevel',            ntimelevel,  want_stages)
  CALL check_int('A_E is square',         SIZE(A_E,1), s)
  CALL check_int('A_E second dimension',  SIZE(A_E,2), s)
  CALL check_int('b_E first dimension',   SIZE(b_E,1), s)
  IF (want_embedded_order .GT. 0) THEN
    CALL check_int('b_E holds an embedded pair', SIZE(b_E,2), 2)
  ELSE
    CALL check_int('b_E holds a single weight set', SIZE(b_E,2), 1)
  ENDIF

  CALL section('the method is explicit')
  ! A strictly lower triangular A is what makes the stages computable in order.
  explicit_ok = .TRUE.
  DO i = 1,s
    DO j = i,s
      IF (A_E(i,j) .NE. 0._xp) explicit_ok = .FALSE.
    ENDDO
  ENDDO
  CALL check('A_E is strictly lower triangular', explicit_ok)

  CALL section('stage times match the stage weights')
  ! c_i = sum_j A_ij, otherwise the stages are evaluated at the wrong times.
  consistent_ok = .TRUE.
  DO i = 1,s
    IF (ABS(c_E(i) - SUM(A_E(i,1:s))) .GT. tol_of(1.0e3_xp)) consistent_ok = .FALSE.
  ENDDO
  CALL check('row sums of A_E equal c_E', consistent_ok)
  CALL check_real('first stage is at the start of the step', c_E(1), 0._xp)

  CALL section('order conditions, main weights')
  CALL check_order(A_E, b_E(:,1), c_E, want_order)

  IF (want_embedded_order .GT. 0) THEN
    CALL section('order conditions, embedded weights')
    CALL check_order(A_E, b_E(:,2), c_E, want_embedded_order)
  ENDIF

  CALL finish_suite

CONTAINS

  !> Check every rooted tree condition of an explicit Runge-Kutta method up to
  !! the requested order. Satisfying all of them is necessary and sufficient for
  !! the method to be of that order.
  SUBROUTINE check_order(a, b, c, ord)
    REAL(xp), INTENT(IN) :: a(:,:), b(:), c(:)
    INTEGER,  INTENT(IN) :: ord
    REAL(xp), DIMENSION(SIZE(c)) :: ac, ac2, ac3, a2c
    REAL(xp) :: tol

    tol = tol_of(1.0e5_xp)
    ac  = MATMUL(a, c)
    ac2 = MATMUL(a, c**2)
    ac3 = MATMUL(a, c**3)
    a2c = MATMUL(a, ac)

    ! order 1
    CALL check_real('sum b = 1', SUM(b), 1._xp, tol)
    IF (ord .LT. 2) RETURN

    ! order 2
    CALL check_real('b.c = 1/2', DOT_PRODUCT(b,c), 1._xp/2._xp, tol)
    IF (ord .LT. 3) RETURN

    ! order 3
    CALL check_real('b.c^2 = 1/3',  DOT_PRODUCT(b,c**2), 1._xp/3._xp, tol)
    CALL check_real('b.(Ac) = 1/6', DOT_PRODUCT(b,ac),   1._xp/6._xp, tol)
    IF (ord .LT. 4) RETURN

    ! order 4
    CALL check_real('b.c^3 = 1/4',      DOT_PRODUCT(b,c**3),  1._xp/4._xp,  tol)
    CALL check_real('b.(c*Ac) = 1/8',   DOT_PRODUCT(b,c*ac),  1._xp/8._xp,  tol)
    CALL check_real('b.(Ac^2) = 1/12',  DOT_PRODUCT(b,ac2),   1._xp/12._xp, tol)
    CALL check_real('b.(A^2c) = 1/24',  DOT_PRODUCT(b,a2c),   1._xp/24._xp, tol)
    IF (ord .LT. 5) RETURN

    ! order 5
    CALL check_real('b.c^4 = 1/5',           DOT_PRODUCT(b,c**4),      1._xp/5._xp,   tol)
    CALL check_real('b.(c^2*Ac) = 1/10',     DOT_PRODUCT(b,c**2*ac),   1._xp/10._xp,  tol)
    CALL check_real('b.(c*Ac^2) = 1/15',     DOT_PRODUCT(b,c*ac2),     1._xp/15._xp,  tol)
    CALL check_real('b.(Ac)^2 = 1/20',       DOT_PRODUCT(b,ac*ac),     1._xp/20._xp,  tol)
    CALL check_real('b.(Ac^3) = 1/20',       DOT_PRODUCT(b,ac3),       1._xp/20._xp,  tol)
    CALL check_real('b.(c*A^2c) = 1/30',     DOT_PRODUCT(b,c*a2c),     1._xp/30._xp,  tol)
    CALL check_real('b.A(c*Ac) = 1/40',      DOT_PRODUCT(b,MATMUL(a,c*ac)),  1._xp/40._xp,  tol)
    CALL check_real('b.A^2c^2 = 1/60',       DOT_PRODUCT(b,MATMUL(a,ac2)),   1._xp/60._xp,  tol)
    CALL check_real('b.A^3c = 1/120',        DOT_PRODUCT(b,MATMUL(a,a2c)),   1._xp/120._xp, tol)
  END SUBROUTINE check_order

END PROGRAM test_time_integration
