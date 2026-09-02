PROGRAM test_prec_const
  ! Checks the hard-coded literals of prec_const against values recomputed at
  ! run time. Every constant there is typed out digit by digit, so a fat finger
  ! is invisible at compile time and silently wrong everywhere downstream.
  USE prec_const
  USE test_harness
  IMPLICIT NONE
  REAL(xp) :: tol

  CALL start_suite('prec_const')
  tol = tol_of(4._xp) ! these are exact to the last bit up to rounding of the literal

  CALL section('working precision')
#ifdef SINGLE_PRECISION
  CALL check_int('xp is REAL32', xp, sp)
#else
  CALL check_int('xp is REAL64', xp, dp)
#endif
  CALL check('single and double precision differ', sp .NE. dp)

  CALL section('pi and its multiples')
  CALL check_real('PI',     PI,     4._xp*ATAN(1._xp),        tol)
  CALL check_real('PIO2',   PIO2,   2._xp*ATAN(1._xp),        tol)
  CALL check_real('TWOPI',  TWOPI,  8._xp*ATAN(1._xp),        tol)
  CALL check_real('ONEOPI', ONEOPI, 1._xp/(4._xp*ATAN(1._xp)),tol)

  CALL section('square roots')
  CALL check_real('SQRT2',    SQRT2,    SQRT(2._xp),      tol)
  CALL check_real('SQRT3',    SQRT3,    SQRT(3._xp),      tol)
  CALL check_real('SQRT6',    SQRT6,    SQRT(6._xp),      tol)
  CALL check_real('INVSQRT2', INVSQRT2, 1._xp/SQRT(2._xp),tol)
  CALL check_real('SQRT2*SQRT3 = SQRT6', SQRT2*SQRT3, SQRT6, tol)

  CALL section('rational constants')
  CALL check_real('onetwelfth',        onetwelfth,        1._xp/12._xp, tol)
  CALL check_real('onetwentyfourth',   onetwentyfourth,   1._xp/24._xp, tol)
  CALL check_real('onethird',          onethird,          1._xp/3._xp,  tol)
  CALL check_real('twothird',          twothird,          2._xp/3._xp,  tol)
  CALL check_real('onesixth',          onesixth,          1._xp/6._xp,  tol)
  CALL check_real('fivesixths',        fivesixths,        5._xp/6._xp,  tol)
  CALL check_real('sevensixths',       sevensixths,       7._xp/6._xp,  tol)
  CALL check_real('elevensixths',      elevensixths,     11._xp/6._xp,  tol)
  CALL check_real('nineeighths',       nineeighths,       9._xp/8._xp,  tol)
  CALL check_real('onesixteenth',      onesixteenth,      1._xp/16._xp, tol)
  CALL check_real('ninesixteenths',    ninesixteenths,    9._xp/16._xp, tol)
  CALL check_real('thirteentwelfths',  thirteentwelfths, 13._xp/12._xp, tol)

  CALL section('imaginary unit')
  CALL check_cmpx('imagu**2 = -1', imagu*imagu, (-1._xp,0._xp), tol)
  CALL check_real('|imagu| = 1',   ABS(imagu),  1._xp,          tol)

  CALL finish_suite

END PROGRAM test_prec_const
