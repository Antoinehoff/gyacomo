PROGRAM test_calculus
  ! Exercises the parallel-direction finite difference stencils of calculus_mod.
  !
  ! Each stencil is a polynomial rule of a known degree, so the sharpest test is
  ! exactness: feed it a polynomial of that degree and the answer must be the
  ! analytic derivative to round off. That pins both the weights and the index
  ! offsets, the latter being easy to get wrong given the ghost cell layout
  !
  !   |x|x|_ _ ... _|x|x|      storage index k holds z = z0 + (k-1)*dz
  !    1 2 3       N+2 N+3     interior point iz lives at k = iz + ngz/2
  !
  ! and the half cell shift of the staggered (odd <-> even) variants.
  USE prec_const, ONLY: xp, PI, imagu
  USE calculus,   ONLY: grad_z, grad_z2, grad_z4, grad_z_5D, grad_z4_5D, interp_z
  USE test_harness
  IMPLICIT NONE

  INTEGER,  PARAMETER :: NZ = 12, NGZ = 4
  REAL(xp), PARAMETER :: DZ = 0.5_xp, Z0 = -4._xp
  ! Reference polynomial coefficients, real and imaginary part
  REAL(xp), PARAMETER :: ar(0:5) = (/ 0.30_xp,-0.80_xp, 0.45_xp,-0.25_xp, 0.12_xp,-0.05_xp /)
  REAL(xp), PARAMETER :: ai(0:5) = (/-0.15_xp, 0.60_xp,-0.35_xp, 0.18_xp,-0.09_xp, 0.04_xp /)

  COMPLEX(xp) :: f(NZ+NGZ), g(NZ), ref(NZ)
  COMPLEX(xp) :: f5(1,1,1,1,1,NZ+NGZ), g5(1,1,1,1,1,NZ)
  COMPLEX(xp) :: fs(NZ+2), gs(NZ)
  REAL(xp)    :: invdz, z, err1, err2
  INTEGER     :: iz, deg

  CALL start_suite('calculus')
  invdz = 1._xp/DZ

  !--------------------------------------------------- centered first derivative
  ! dz_usu is the 4th order centered rule: exact up to degree 4.
  CALL section('grad_z, centered: exact on a quartic')
  deg = 4
  CALL fill(f, deg)
  DO iz = 1,NZ
    ref(iz) = dpoly(zc(iz+NGZ/2), deg)
  ENDDO
  CALL grad_z(0, NZ, NGZ, invdz, f, g)
  CALL check_cmpx_arr('d/dz of a quartic', g, ref, tol_of(1.0e4_xp))

  !------------------------------------------- staggered first derivatives (o2e/e2o)
  ! dz_o2e / dz_e2o evaluate the derivative half a cell away from the input
  ! points, on the companion staggered grid: exact up to degree 3.
  CALL section('grad_z, odd to even: exact on a cubic, half cell shifted')
  deg = 3
  CALL fill(f, deg)
  DO iz = 1,NZ
    ref(iz) = dpoly(zc(iz+1) + 0.5_xp*DZ, deg)
  ENDDO
  CALL grad_z(1, NZ, NGZ, invdz, f, g)
  CALL check_cmpx_arr('d/dz of a cubic on the even grid', g, ref, tol_of(1.0e4_xp))

  CALL section('grad_z, even to odd: exact on a cubic, half cell shifted')
  DO iz = 1,NZ
    ref(iz) = dpoly(zc(iz+2) + 0.5_xp*DZ, deg)
  ENDDO
  CALL grad_z(2, NZ, NGZ, invdz, f, g)
  CALL check_cmpx_arr('d/dz of a cubic on the odd grid', g, ref, tol_of(1.0e4_xp))

  !--------------------------------------------------------- second derivative
  ! dz2_usu is the 4th order centered second derivative: exact up to degree 5.
  CALL section('grad_z2: exact on a quintic')
  deg = 5
  CALL fill(f, deg)
  DO iz = 1,NZ
    ref(iz) = d2poly(zc(iz+NGZ/2), deg)
  ENDDO
  CALL grad_z2(NZ, NGZ, invdz, f, g)
  CALL check_cmpx_arr('d2/dz2 of a quintic', g, ref, tol_of(1.0e5_xp))

  !--------------------------------------------------------- fourth derivative
  ! dz4_usu is the 2nd order centered fourth derivative used for the parallel
  ! hyperdiffusion: exact up to degree 5.
  CALL section('grad_z4: exact on a quintic')
  DO iz = 1,NZ
    ref(iz) = d4poly(zc(iz+NGZ/2), deg)
  ENDDO
  CALL grad_z4(NZ, NGZ, invdz, f, g)
  CALL check_cmpx_arr('d4/dz4 of a quintic', g, ref, tol_of(1.0e6_xp))

  !------------------------------------------------------------- interpolation
  CALL section('interp_z, no staggering: strips the ghost cells')
  deg = 5
  CALL fill(f, deg)
  DO iz = 1,NZ
    ref(iz) = f(iz+NGZ/2)
  ENDDO
  CALL interp_z(0, NZ, NGZ, f, g)
  CALL check_cmpx_arr('identity on the interior', g, ref, tol_of(1._xp))

  ! iz_o2e / iz_e2o are the 4th order midpoint rules: exact up to degree 3.
  CALL section('interp_z, odd to even: exact on a cubic')
  deg = 3
  CALL fill(f, deg)
  DO iz = 1,NZ
    ref(iz) = poly(zc(iz+1) + 0.5_xp*DZ, deg)
  ENDDO
  CALL interp_z(1, NZ, NGZ, f, g)
  CALL check_cmpx_arr('cubic onto the even grid', g, ref, tol_of(1.0e4_xp))

  CALL section('interp_z, even to odd: exact on a cubic')
  DO iz = 1,NZ
    ref(iz) = poly(zc(iz+2) + 0.5_xp*DZ, deg)
  ENDDO
  CALL interp_z(2, NZ, NGZ, f, g)
  CALL check_cmpx_arr('cubic onto the odd grid', g, ref, tol_of(1.0e4_xp))

  !------------------------------------------------------- 5D array wrappers
  ! The 5D variants duplicate the 1D stencils for the whole moment hierarchy;
  ! they must agree with the scalar routines point for point.
  CALL section('grad_z_5D and grad_z4_5D match the 1D routines')
  deg = 5
  CALL fill(f, deg)
  f5(1,1,1,1,1,:) = f
  CALL grad_z(0, NZ, NGZ, invdz, f, g)
  CALL grad_z_5D(NZ, NGZ, invdz, f5, g5)
  CALL check_cmpx_arr('grad_z_5D vs grad_z', g5(1,1,1,1,1,:), g, tol_of(1._xp))
  CALL grad_z4(NZ, NGZ, invdz, f, g)
  CALL grad_z4_5D(NZ, NGZ, invdz, f5, g5)
  CALL check_cmpx_arr('grad_z4_5D vs grad_z4', g5(1,1,1,1,1,:), g, tol_of(1._xp))

  !-------------------------------------------------- too few ghost cells guard
  ! With fewer than four ghost cells the stencil does not fit and every routine
  ! must return zero rather than read out of bounds.
  CALL section('fewer than 4 ghost cells returns zero')
  fs = (1._xp, 1._xp)
  gs = (9._xp, 9._xp)
  CALL grad_z(0, NZ, 2, invdz, fs, gs)
  CALL check('grad_z zeroed',    MAXVAL(ABS(gs)) .EQ. 0._xp)
  gs = (9._xp, 9._xp)
  CALL grad_z2(NZ, 2, invdz, fs, gs)
  CALL check('grad_z2 zeroed',   MAXVAL(ABS(gs)) .EQ. 0._xp)
  gs = (9._xp, 9._xp)
  CALL grad_z4(NZ, 2, invdz, fs, gs)
  CALL check('grad_z4 zeroed',   MAXVAL(ABS(gs)) .EQ. 0._xp)

  !------------------------------------------------------------- order of accuracy
  ! On a function the stencil is not exact for, halving dz must divide the error
  ! by 2^4. This is what actually certifies the rule is 4th and not 2nd order.
  CALL section('grad_z is fourth order accurate')
  err1 = fourier_error(16)
  err2 = fourier_error(32)
  CALL check('error decreases when refining', err2 .LT. err1)
  CALL check_real('error ratio is 2^4', err1/err2, 16._xp, 0.15_xp)

  CALL finish_suite

CONTAINS

  !> Position stored at index k of a ghosted z array.
  REAL(xp) FUNCTION zc(k)
    INTEGER, INTENT(IN) :: k
    zc = Z0 + REAL(k-1,xp)*DZ
  END FUNCTION zc

  COMPLEX(xp) FUNCTION poly(x, deg)
    REAL(xp), INTENT(IN) :: x
    INTEGER,  INTENT(IN) :: deg
    INTEGER :: k
    poly = (0._xp,0._xp)
    DO k = 0,deg
      poly = poly + CMPLX(ar(k), ai(k), xp)*x**k
    ENDDO
  END FUNCTION poly

  COMPLEX(xp) FUNCTION dpoly(x, deg)
    REAL(xp), INTENT(IN) :: x
    INTEGER,  INTENT(IN) :: deg
    INTEGER :: k
    dpoly = (0._xp,0._xp)
    DO k = 1,deg
      dpoly = dpoly + REAL(k,xp)*CMPLX(ar(k), ai(k), xp)*x**(k-1)
    ENDDO
  END FUNCTION dpoly

  COMPLEX(xp) FUNCTION d2poly(x, deg)
    REAL(xp), INTENT(IN) :: x
    INTEGER,  INTENT(IN) :: deg
    INTEGER :: k
    d2poly = (0._xp,0._xp)
    DO k = 2,deg
      d2poly = d2poly + REAL(k*(k-1),xp)*CMPLX(ar(k), ai(k), xp)*x**(k-2)
    ENDDO
  END FUNCTION d2poly

  COMPLEX(xp) FUNCTION d4poly(x, deg)
    REAL(xp), INTENT(IN) :: x
    INTEGER,  INTENT(IN) :: deg
    INTEGER :: k
    d4poly = (0._xp,0._xp)
    DO k = 4,deg
      d4poly = d4poly + REAL(k*(k-1)*(k-2)*(k-3),xp)*CMPLX(ar(k), ai(k), xp)*x**(k-4)
    ENDDO
  END FUNCTION d4poly

  !> Sample the reference polynomial over a ghosted array.
  SUBROUTINE fill(arr, deg)
    COMPLEX(xp), INTENT(OUT) :: arr(:)
    INTEGER,     INTENT(IN)  :: deg
    INTEGER :: k
    DO k = 1,SIZE(arr)
      arr(k) = poly(zc(k), deg)
    ENDDO
  END SUBROUTINE fill

  !> Max error of the centered stencil on exp(i z) sampled over one period with
  !! n interior points and periodic ghost cells.
  REAL(xp) FUNCTION fourier_error(n)
    INTEGER, INTENT(IN) :: n
    COMPLEX(xp), ALLOCATABLE :: fp(:), dfp(:)
    REAL(xp) :: h, zz
    INTEGER  :: k
    ALLOCATE(fp(n+NGZ), dfp(n))
    h = 2._xp*PI/REAL(n,xp)
    ! index k holds z = (k-1-ngz/2)*h, so the ghosts wrap around by periodicity
    DO k = 1,n+NGZ
      zz    = REAL(k-1-NGZ/2,xp)*h
      fp(k) = EXP(imagu*zz)
    ENDDO
    CALL grad_z(0, n, NGZ, 1._xp/h, fp, dfp)
    fourier_error = 0._xp
    DO k = 1,n
      zz = REAL(k-1,xp)*h
      fourier_error = MAX(fourier_error, ABS(dfp(k) - imagu*EXP(imagu*zz)))
    ENDDO
    DEALLOCATE(fp, dfp)
  END FUNCTION fourier_error

END PROGRAM test_calculus
