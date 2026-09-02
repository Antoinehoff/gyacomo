PROGRAM test_utility
  ! Exercises the helpers of utility_mod.
  !
  ! decomp1D is the single point where the code decides which slice of p, ky or
  ! z each rank owns. If it ever drops or duplicates an index the run does not
  ! crash, it silently integrates the wrong domain, so the invariants below are
  ! checked exhaustively over a range of (n, numprocs) rather than spot-checked.
  USE prec_const, ONLY: xp
  USE utility,    ONLY: decomp1D, is_nan, is_inf
  USE test_harness
  USE, INTRINSIC :: ieee_arithmetic
  IMPLICIT NONE

  INTEGER, PARAMETER :: NMAX = 64, PMAX = 9
  INTEGER :: n, np, id, is, ie
  INTEGER :: owner(NMAX), nloc(0:PMAX-1)
  LOGICAL :: covered_ok, balanced_ok, contiguous_ok, ordered_ok
  INTEGER :: prev_ie, sz, szmin, szmax, i
  REAL    :: nan_val, inf_val

  CALL start_suite('utility')

  !-------------------------------------------------------------------- decomp1D
  ! Every index of 1..n must be owned by exactly one rank.
  CALL section('decomp1D: partition covers 1..n exactly once')
  covered_ok = .TRUE.
  DO n = 1,NMAX
    DO np = 1,PMAX
      owner(1:n) = 0
      DO id = 0,np-1
        CALL decomp1D(n, np, id, is, ie)
        DO i = is,ie
          IF (i .LT. 1 .OR. i .GT. n) THEN
            covered_ok = .FALSE.
          ELSE
            owner(i) = owner(i) + 1
          ENDIF
        ENDDO
      ENDDO
      IF (ANY(owner(1:n) .NE. 1)) covered_ok = .FALSE.
    ENDDO
  ENDDO
  CALL check('no index lost or owned twice, n<=64 and np<=9', covered_ok)

  ! Slices must follow each other without gap and in rank order.
  CALL section('decomp1D: slices are contiguous and ordered by rank')
  contiguous_ok = .TRUE.
  ordered_ok    = .TRUE.
  DO n = 1,NMAX
    DO np = 1,PMAX
      prev_ie = 0
      DO id = 0,np-1
        CALL decomp1D(n, np, id, is, ie)
        IF (ie .GE. is) THEN ! non-empty slice
          IF (is .NE. prev_ie+1) contiguous_ok = .FALSE.
          IF (is .LE. prev_ie)   ordered_ok    = .FALSE.
          prev_ie = ie
        ENDIF
      ENDDO
      IF (prev_ie .NE. n) contiguous_ok = .FALSE.
    ENDDO
  ENDDO
  CALL check('slices tile [1,n] back to back', contiguous_ok)
  CALL check('slice starts increase with rank', ordered_ok)

  ! Load balance: no rank may carry more than one extra point.
  CALL section('decomp1D: load is balanced to within one point')
  balanced_ok = .TRUE.
  DO n = 1,NMAX
    DO np = 1,PMAX
      DO id = 0,np-1
        CALL decomp1D(n, np, id, is, ie)
        nloc(id) = MAX(0, ie-is+1)
      ENDDO
      szmin = MINVAL(nloc(0:np-1))
      szmax = MAXVAL(nloc(0:np-1))
      IF (szmax-szmin .GT. 1) balanced_ok = .FALSE.
      IF (SUM(nloc(0:np-1)) .NE. n) balanced_ok = .FALSE.
    ENDDO
  ENDDO
  CALL check('max slice - min slice <= 1', balanced_ok)

  CALL section('decomp1D: known cases')
  CALL decomp1D(10, 1, 0, is, ie)
  CALL check_int('n=10 np=1 id=0 start', is, 1)
  CALL check_int('n=10 np=1 id=0 end',   ie, 10)
  ! 10 over 4 -> 3,3,2,2
  CALL decomp1D(10, 4, 0, is, ie); sz = ie-is+1
  CALL check_int('n=10 np=4 id=0 size', sz, 3)
  CALL decomp1D(10, 4, 1, is, ie); sz = ie-is+1
  CALL check_int('n=10 np=4 id=1 size', sz, 3)
  CALL decomp1D(10, 4, 2, is, ie); sz = ie-is+1
  CALL check_int('n=10 np=4 id=2 size', sz, 2)
  CALL decomp1D(10, 4, 3, is, ie)
  CALL check_int('n=10 np=4 id=3 size', ie-is+1, 2)
  CALL check_int('last rank ends on n',  ie, 10)
  ! More ranks than points: the surplus ranks get an empty slice (ie < is).
  CALL decomp1D(3, 5, 4, is, ie)
  CALL check('n=3 np=5 id=4 slice is empty', ie .LT. is)
  CALL decomp1D(3, 5, 0, is, ie)
  CALL check_int('n=3 np=5 id=0 owns one point', ie-is+1, 1)

  !------------------------------------------------------------ is_nan / is_inf
  ! Note: these take a default REAL, not REAL(xp).
  CALL section('is_nan')
  nan_val = ieee_value(1.0, ieee_quiet_nan)
  CALL check('NaN is detected',            is_nan(nan_val, ''))
  CALL check('a finite value is not NaN', .NOT. is_nan(1.0,  ''))
  CALL check('zero is not NaN',           .NOT. is_nan(0.0,  ''))

  CALL section('is_inf')
  inf_val = ieee_value(1.0, ieee_positive_inf)
  CALL check('+Inf is detected',           is_inf( inf_val, ''))
  CALL check('-Inf is detected',           is_inf(-inf_val, ''))
  CALL check('a finite value is not Inf', .NOT. is_inf(1.0,  ''))
  CALL check('zero is not Inf',           .NOT. is_inf(0.0,  ''))

  CALL finish_suite

END PROGRAM test_utility
