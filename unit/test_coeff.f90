PROGRAM test_coeff
  ! Verifies the Laguerre triple product coefficients of coeff_mod.
  !
  ! ALL2L(n,k,s,0) returns d_nks defined by L_n L_k = sum_s d_nks L_s, that is
  !
  !     d_nks = int_0^inf exp(-x) L_n(x) L_k(x) L_s(x) dx
  !
  ! since the Laguerre polynomials are orthonormal for that weight. These
  ! coefficients are the entire Laguerre content of the nonlinear term (they
  ! populate the dnjs table), and they are computed through a multiple precision
  ! triple sum that no one is going to check by inspection. The reference values
  ! below were obtained by exact rational integration of the same integral.
  USE prec_const, ONLY: xp
  USE coeff,      ONLY: ALL2L
  USE FMZM,       ONLY: TO_DP
  USE test_harness
  IMPLICIT NONE

  INTEGER,  PARAMETER :: NMAX = 3 ! degrees swept by the property checks
  REAL(xp) :: tol, dv, dref
  INTEGER  :: n, k, s
  LOGICAL  :: sym_ok, cutoff_ok

  CALL start_suite('coeff')
  ! FM carries far more digits than xp; the loss is in the final TO_DP.
  tol = tol_of(1.0e4_xp)

  !------------------------------------------------------------- known values
  CALL section('tabulated values of d_nks')
  CALL check_real('d(1,1,0)', dnks(1,1,0),  1._xp, tol)
  CALL check_real('d(1,1,1)', dnks(1,1,1), -2._xp, tol)
  CALL check_real('d(1,1,2)', dnks(1,1,2),  2._xp, tol)
  CALL check_real('d(1,2,1)', dnks(1,2,1),  2._xp, tol)
  CALL check_real('d(1,2,2)', dnks(1,2,2), -4._xp, tol)
  CALL check_real('d(1,2,3)', dnks(1,2,3),  3._xp, tol)
  CALL check_real('d(2,2,1)', dnks(2,2,1), -4._xp, tol)
  CALL check_real('d(2,2,2)', dnks(2,2,2), 10._xp, tol)
  CALL check_real('d(2,2,3)', dnks(2,2,3),-12._xp, tol)
  CALL check_real('d(2,2,4)', dnks(2,2,4),  6._xp, tol)
  CALL check_real('d(2,3,4)', dnks(2,3,4),-24._xp, tol)
  CALL check_real('d(2,3,5)', dnks(2,3,5), 10._xp, tol)

  !------------------------------------------------------------- L_0 is the unit
  CALL section('multiplying by L_0 is the identity')
  ! L_0 = 1, so d_0ks must be the Kronecker delta on (k,s).
  DO k = 0,NMAX
    DO s = 0,NMAX
      dref = 0._xp
      IF (k .EQ. s) dref = 1._xp
      CALL check_real('d(0,k,s) = delta_ks', dnks(0,k,s), dref, tol)
    ENDDO
  ENDDO

  !------------------------------------------------------------- orthonormality
  CALL section('projection on L_0 recovers orthonormality')
  ! d_nk0 = int exp(-x) L_n L_k = delta_nk.
  DO n = 0,NMAX
    DO k = 0,NMAX
      dref = 0._xp
      IF (n .EQ. k) dref = 1._xp
      CALL check_real('d(n,k,0) = delta_nk', dnks(n,k,0), dref, tol)
    ENDDO
  ENDDO

  !------------------------------------------------------------- degree cutoff
  CALL section('coefficients vanish above degree n+k')
  ! L_n L_k has degree n+k, so it is orthogonal to every higher L_s.
  cutoff_ok = .TRUE.
  DO n = 0,2
    DO k = 0,2
      DO s = n+k+1, n+k+2
        IF (ABS(dnks(n,k,s)) .GT. tol) cutoff_ok = .FALSE.
      ENDDO
    ENDDO
  ENDDO
  CALL check('d(n,k,s) = 0 for s > n+k', cutoff_ok)

  !------------------------------------------------------------- symmetry
  CALL section('d_nks is symmetric in its three indices')
  ! build_dnjs_table fills five of the six permutations by assignment rather
  ! than by evaluation, so the symmetry it assumes has to actually hold.
  sym_ok = .TRUE.
  DO n = 0,NMAX
    DO k = n,NMAX
      DO s = k,NMAX
        dv = dnks(n,k,s)
        IF (ABS(dnks(n,s,k) - dv) .GT. tol) sym_ok = .FALSE.
        IF (ABS(dnks(k,n,s) - dv) .GT. tol) sym_ok = .FALSE.
        IF (ABS(dnks(k,s,n) - dv) .GT. tol) sym_ok = .FALSE.
        IF (ABS(dnks(s,n,k) - dv) .GT. tol) sym_ok = .FALSE.
        IF (ABS(dnks(s,k,n) - dv) .GT. tol) sym_ok = .FALSE.
      ENDDO
    ENDDO
  ENDDO
  CALL check('all six permutations agree', sym_ok)

  CALL finish_suite

CONTAINS

  !> d_nks as a working precision real, the way numerics_mod consumes it.
  REAL(xp) FUNCTION dnks(n, k, s)
    INTEGER, INTENT(IN) :: n, k, s
    dnks = REAL(TO_DP(ALL2L(n, k, s, 0)), xp)
  END FUNCTION dnks

END PROGRAM test_coeff
