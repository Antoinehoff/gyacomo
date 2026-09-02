PROGRAM test_numerics
  ! Exercises the Taylor kernel approximation of numerics_mod.
  !
  ! The gyroaverage kernels are K_n(y) = y^n exp(-y)/n!. When KN_MODEL is set to
  ! 'taylor', evaluate_kernels replaces them by taylor_kernel_n, the truncation
  ! of that expression at a given power of y. Its value is only as good as its
  ! agreement with the closed form, so that is what is measured here, together
  ! with the boundary behaviour the caller relies on (n above the truncation
  ! order is zeroed by evaluate_kernels, n = order must give the leading term).
  USE prec_const, ONLY: xp
  USE numerics,   ONLY: taylor_kernel_n
  USE test_harness
  IMPLICIT NONE

  REAL(xp) :: y, approx, exact, e_lo, e_hi
  INTEGER  :: n, order

  CALL start_suite('numerics')

  !--------------------------------------------------------- degenerate values
  CALL section('taylor_kernel_n at y = 0')
  CALL check_real('K_0(0) = 1', taylor_kernel_n(8, 0, 0._xp), 1._xp)
  DO n = 1,4
    CALL check_real('K_n(0) = 0 for n > 0', taylor_kernel_n(8, n, 0._xp), 0._xp)
  ENDDO

  CALL section('lowest order truncation is the leading term y^n/n!')
  DO n = 0,4
    y = 0.3_xp
    CALL check_real('order = n gives y^n/n!', &
                    taylor_kernel_n(n, n, y), y**n/gamma_int(n))
  ENDDO

  CALL section('a truncation shorter than n is empty')
  ! The DO loop runs from n to order, so order < n yields no terms at all.
  ! evaluate_kernels relies on this never being reached: it zeroes j > ORDER
  ! itself. The value is pinned here so that assumption stays visible.
  CALL check_real('order < n gives 0', taylor_kernel_n(2, 5, 0.3_xp), 0._xp)

  !------------------------------------------------------ agreement with y^n e^-y/n!
  CALL section('converges to the exact kernel for small y')
  DO n = 0,3
    y      = 0.25_xp
    approx = taylor_kernel_n(30, n, y)
    exact  = y**n*EXP(-y)/gamma_int(n)
    CALL check_real('K_n(0.25), order 30', approx, exact, tol_of(1.0e3_xp))
  ENDDO
  DO n = 0,3
    y      = 2.0_xp
    approx = taylor_kernel_n(40, n, y)
    exact  = y**n*EXP(-y)/gamma_int(n)
    CALL check_real('K_n(2.0), order 40', approx, exact, tol_of(1.0e4_xp))
  ENDDO

  CALL section('truncation error shrinks as the order grows')
  ! Not just "it converges": each extra power of y must actually help, which is
  ! what makes ORDER a meaningful input knob.
  y = 0.5_xp
  DO n = 0,3
    exact = y**n*EXP(-y)/gamma_int(n)
    e_lo  = ABS(taylor_kernel_n(n+2, n, y) - exact)
    e_hi  = ABS(taylor_kernel_n(n+6, n, y) - exact)
    CALL check('order n+6 beats order n+2', e_hi .LT. e_lo)
  ENDDO

  CALL section('kernels sum to one, as the exact ones do')
  ! sum_n K_n(y) = 1 exactly. A high order truncation must reproduce it.
  DO order = 20,30,10
    y      = 0.4_xp
    approx = 0._xp
    DO n = 0,order
      approx = approx + taylor_kernel_n(order, n, y)
    ENDDO
    CALL check_real('sum over n of K_n(0.4)', approx, 1._xp, tol_of(1.0e4_xp))
  ENDDO

  CALL finish_suite

CONTAINS

  !> n! as a real, matching the GAMMA(n+1) used inside the module.
  REAL(xp) FUNCTION gamma_int(n)
    INTEGER, INTENT(IN) :: n
    gamma_int = GAMMA(REAL(n,xp) + 1._xp)
  END FUNCTION gamma_int

END PROGRAM test_numerics
