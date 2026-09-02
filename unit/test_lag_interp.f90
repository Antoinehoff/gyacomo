PROGRAM test_lag_interp
  ! Exercises lagrange_interpolation, the third order interpolator used by the
  ! Miller and circular geometry modules to resample the equilibrium.
  !
  ! A cubic Lagrange rule reproduces any cubic exactly, whatever the node
  ! spacing, and its derivative rule reproduces the cubic's derivative exactly.
  ! That is the property tested here: it pins the stencil selection and the
  ! weights at once, and it holds for the descending grids and the
  ! extrapolation regime that the routine explicitly supports.
  USE prec_const, ONLY: xp
  USE lagrange_interpolation
  USE test_harness
  IMPLICIT NONE

  INTEGER,  PARAMETER :: NIN = 9, NOUT = 17
  ! Coefficients of the reference cubic p(x) = c0 + c1 x + c2 x^2 + c3 x^3
  REAL(xp), PARAMETER :: c0 = 0.7_xp, c1 = -1.3_xp, c2 = 0.45_xp, c3 = 0.2_xp

  REAL(xp) :: xin(NIN), yin(NIN), xout(NOUT), yout(NOUT), yref(NOUT)
  REAL(xp) :: xin_d(NIN), yin_d(NIN), xout_d(NOUT), yout_d(NOUT), yref_d(NOUT)
  COMPLEX(xp) :: zin(NIN), zout(NOUT), zref(NOUT)
  REAL(xp) :: ys, dys, tol
  INTEGER  :: i, j

  ! 2D case
  INTEGER,  PARAMETER :: N1I = 8, N2I = 9, N1O = 5, N2O = 6
  REAL(xp) :: x1i(N1I), x2i(N2I), x1o(N1O), x2o(N2O)
  REAL(xp) :: f2i(N1I,N2I), f2o(N1O,N2O), f2ref(N1O,N2O)

  CALL start_suite('lagrange_interpolation')
  tol = tol_of(1.0e4_xp)

  !------------------------------------------------------- non-uniform ascending grid
  ! Deliberately unevenly spaced: a cubic rule must not need equal spacing.
  DO i = 1,NIN
    xin(i) = -2._xp + 0.5_xp*REAL(i-1,xp) + 0.07_xp*SIN(REAL(i,xp))
    yin(i) = cubic(xin(i))
  ENDDO
  ! Output points spread over the input range and a little beyond it on both
  ! sides, since the routine documents that it extrapolates.
  DO j = 1,NOUT
    xout(j) = xin(1) - 0.4_xp + (xin(NIN)-xin(1)+0.8_xp)*REAL(j-1,xp)/REAL(NOUT-1,xp)
    yref(j) = cubic(xout(j))
  ENDDO

  CALL section('lag3interp reproduces a cubic on an ascending grid')
  CALL lag3interp(yin, xin, NIN, yout, xout, NOUT)
  CALL check_real_arr('interpolated cubic', yout, yref, tol)

  CALL section('lag3interp is exact at the input nodes')
  CALL lag3interp(yin, xin, NIN, yout(1:NIN), xin, NIN)
  CALL check_real_arr('nodes returned unchanged', yout(1:NIN), yin, tol)

  CALL section('lag3deriv reproduces the derivative of a cubic')
  DO j = 1,NOUT
    yref(j) = dcubic(xout(j))
  ENDDO
  CALL lag3deriv(yin, xin, NIN, yout, xout, NOUT)
  CALL check_real_arr('interpolated derivative', yout, yref, tol)

  CALL section('scalar wrappers agree with the array versions')
  CALL lag3interp(yin, xin, NIN, yout, xout, NOUT)
  CALL lag3interp(yin, xin, NIN, ys, xout(5))
  CALL check_real('lag3interp_scalar', ys, yout(5), tol)
  CALL lag3deriv(yin, xin, NIN, yout, xout, NOUT)
  CALL lag3deriv(yin, xin, NIN, dys, xout(5))
  CALL check_real('lag3deriv_scalar', dys, yout(5), tol)

  !------------------------------------------------------- descending grid branch
  ! The routine switches stencil direction on x_in(n) < x_in(1); that branch is
  ! reached by the geometry modules whenever the equilibrium is stored inwards.
  CALL section('descending input grid')
  DO i = 1,NIN
    xin_d(i) = xin(NIN-i+1)
    yin_d(i) = yin(NIN-i+1)
  ENDDO
  DO j = 1,NOUT
    xout_d(j) = xout(NOUT-j+1)
    yref_d(j) = cubic(xout_d(j))
  ENDDO
  CALL lag3interp(yin_d, xin_d, NIN, yout_d, xout_d, NOUT)
  CALL check_real_arr('interpolated cubic, descending', yout_d, yref_d, tol)
  DO j = 1,NOUT
    yref_d(j) = dcubic(xout_d(j))
  ENDDO
  CALL lag3deriv(yin_d, xin_d, NIN, yout_d, xout_d, NOUT)
  CALL check_real_arr('derivative, descending', yout_d, yref_d, tol)

  !------------------------------------------------------- complex variant
  ! lag3interp_complex duplicates the real stencil logic by hand, so it is
  ! checked against the real routine rather than trusted to stay in sync.
  CALL section('lag3interp_complex')
  DO i = 1,NIN
    zin(i) = CMPLX(cubic(xin(i)), dcubic(xin(i)), xp)
  ENDDO
  DO j = 1,NOUT
    zref(j) = CMPLX(cubic(xout(j)), dcubic(xout(j)), xp)
  ENDDO
  CALL lag3interp_complex(zin, xin, NIN, zout, xout, NOUT)
  CALL check_cmpx_arr('complex cubic, real and imaginary parts', zout, zref, tol)

  !------------------------------------------------------- 2D variant
  ! Two successive 1D cubic passes are exact for a tensor product cubic.
  CALL section('lag3interp_2d reproduces a tensor product cubic')
  DO i = 1,N1I
    x1i(i) = -1._xp + 0.4_xp*REAL(i-1,xp)
  ENDDO
  DO j = 1,N2I
    x2i(j) = 0.2_xp + 0.3_xp*REAL(j-1,xp)
  ENDDO
  DO i = 1,N1I
    DO j = 1,N2I
      f2i(i,j) = cubic(x1i(i))*cubic(x2i(j))
    ENDDO
  ENDDO
  DO i = 1,N1O
    x1o(i) = -0.8_xp + 0.5_xp*REAL(i-1,xp)
  ENDDO
  DO j = 1,N2O
    x2o(j) = 0.35_xp + 0.37_xp*REAL(j-1,xp)
  ENDDO
  DO i = 1,N1O
    DO j = 1,N2O
      f2ref(i,j) = cubic(x1o(i))*cubic(x2o(j))
    ENDDO
  ENDDO
  CALL lag3interp_2d(f2i, x1i, N1I, x2i, N2I, f2o, x1o, N1O, x2o, N2O)
  CALL check_real_arr('bicubic interpolation', &
                      RESHAPE(f2o,   (/N1O*N2O/)), &
                      RESHAPE(f2ref, (/N1O*N2O/)), tol)

  CALL finish_suite

CONTAINS

  REAL(xp) FUNCTION cubic(x)
    REAL(xp), INTENT(IN) :: x
    cubic = c0 + c1*x + c2*x*x + c3*x*x*x
  END FUNCTION cubic

  REAL(xp) FUNCTION dcubic(x)
    REAL(xp), INTENT(IN) :: x
    dcubic = c1 + 2._xp*c2*x + 3._xp*c3*x*x
  END FUNCTION dcubic

END PROGRAM test_lag_interp
