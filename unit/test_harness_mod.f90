MODULE test_harness
  ! Minimal assertion and reporting helpers shared by the unit test programs.
  !
  ! A test program calls start_suite once, groups its assertions with section,
  ! and ends with finish_suite. finish_suite finalizes MPI and stops with a
  ! non-zero code as soon as one assertion failed, so that "make check" can
  ! simply look at the exit status of each binary.
  !
  ! Tolerances are expressed as multiples of EPSILON(1._xp) so that the same
  ! test source is meaningful in both the double and the single precision build.
  USE prec_const, ONLY: xp
  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER :: MAX_REPORTED = 12 ! stop echoing individual failures past this

  INTEGER,           SAVE :: n_check    = 0
  INTEGER,           SAVE :: n_fail     = 0
  INTEGER,           SAVE :: n_shown    = 0
  INTEGER,           SAVE :: sect_check = 0
  INTEGER,           SAVE :: sect_fail  = 0
  CHARACTER(len=64), SAVE :: suite_name = 'unnamed'
  CHARACTER(len=64), SAVE :: sect_name  = ''
  LOGICAL,           SAVE :: in_section = .FALSE.

  PUBLIC :: start_suite, finish_suite, section
  PUBLIC :: check, check_int, check_real, check_real_arr, check_cmpx, check_cmpx_arr
  PUBLIC :: tol_of

CONTAINS

  !> Tolerance expressed as a multiple of the working precision epsilon.
  REAL(xp) FUNCTION tol_of(factor)
    REAL(xp), INTENT(IN) :: factor
    tol_of = factor*EPSILON(1._xp)
  END FUNCTION tol_of

  !> Open a test suite. MPI is initialized here unless the caller already did it
  !! (a test that needs the Gyacomo communicators calls ppinit before this).
  SUBROUTINE start_suite(name)
    USE mpi
    CHARACTER(len=*), INTENT(IN) :: name
    LOGICAL :: already
    INTEGER :: ierr
    CALL MPI_INITIALIZED(already, ierr)
    IF (.NOT. already) CALL MPI_INIT(ierr)
    suite_name = name
    n_check = 0; n_fail = 0; n_shown = 0
    WRITE(*,'(A)') '=== '//TRIM(suite_name)//' ==='
  END SUBROUTINE start_suite

  !> Close the current section (if any) and report its tally.
  SUBROUTINE close_section
    IF (.NOT. in_section) RETURN
    IF (sect_fail .EQ. 0) THEN
      WRITE(*,'(A,A,A,I0,A)') '  [ ok ] ', TRIM(sect_name), '  (', sect_check, ' checks)'
    ELSE
      WRITE(*,'(A,A,A,I0,A,I0,A)') '  [FAIL] ', TRIM(sect_name), '  (', &
                                   sect_fail, '/', sect_check, ' failed)'
    ENDIF
    in_section = .FALSE.
  END SUBROUTINE close_section

  !> Start a named group of assertions.
  SUBROUTINE section(name)
    CHARACTER(len=*), INTENT(IN) :: name
    CALL close_section
    sect_name  = name
    sect_check = 0
    sect_fail  = 0
    in_section = .TRUE.
  END SUBROUTINE section

  !> Report the suite tally, finalize MPI and stop 1 if anything failed.
  SUBROUTINE finish_suite
    USE mpi
    INTEGER :: ierr
    LOGICAL :: already
    CALL close_section
    IF (n_fail .EQ. 0) THEN
      WRITE(*,'(A,I0,A)') 'PASS  '//TRIM(suite_name)//'  (', n_check, ' checks)'
    ELSE
      WRITE(*,'(A,I0,A,I0,A)') 'FAIL  '//TRIM(suite_name)//'  (', n_fail, '/', n_check, ' failed)'
    ENDIF
    WRITE(*,'(A)') ''
    FLUSH(6)
    CALL MPI_INITIALIZED(already, ierr)
    IF (already) CALL MPI_FINALIZE(ierr)
    IF (n_fail .NE. 0) STOP 1
  END SUBROUTINE finish_suite

  !> Book-keeping common to every assertion.
  SUBROUTINE record(label, ok, detail)
    CHARACTER(len=*), INTENT(IN) :: label
    LOGICAL,          INTENT(IN) :: ok
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: detail
    n_check    = n_check    + 1
    sect_check = sect_check + 1
    IF (ok) RETURN
    n_fail    = n_fail    + 1
    sect_fail = sect_fail + 1
    n_shown   = n_shown   + 1
    IF (n_shown .LE. MAX_REPORTED) THEN
      IF (PRESENT(detail)) THEN
        WRITE(*,'(A)') '    x '//TRIM(label)//': '//TRIM(detail)
      ELSE
        WRITE(*,'(A)') '    x '//TRIM(label)
      ENDIF
    ELSEIF (n_shown .EQ. MAX_REPORTED+1) THEN
      WRITE(*,'(A)') '    x ... (further failures not listed)'
    ENDIF
  END SUBROUTINE record

  !> Assert that a condition holds.
  SUBROUTINE check(label, cond)
    CHARACTER(len=*), INTENT(IN) :: label
    LOGICAL,          INTENT(IN) :: cond
    CALL record(label, cond)
  END SUBROUTINE check

  !> Assert equality of two integers.
  SUBROUTINE check_int(label, got, want)
    CHARACTER(len=*), INTENT(IN) :: label
    INTEGER,          INTENT(IN) :: got, want
    CHARACTER(len=128) :: msg
    WRITE(msg,'(A,I0,A,I0)') 'got ', got, ', expected ', want
    CALL record(label, got .EQ. want, msg)
  END SUBROUTINE check_int

  !> Assert that two reals agree within tol, measured relative to max(1,|want|).
  SUBROUTINE check_real(label, got, want, tol)
    CHARACTER(len=*), INTENT(IN) :: label
    REAL(xp),         INTENT(IN) :: got, want
    REAL(xp),         INTENT(IN), OPTIONAL :: tol
    REAL(xp) :: t, err
    CHARACTER(len=160) :: msg
    t = tol_of(1.0e3_xp); IF (PRESENT(tol)) t = tol
    err = ABS(got-want)/MAX(1._xp, ABS(want))
    WRITE(msg,'(A,ES14.7,A,ES14.7,A,ES9.2,A,ES9.2)') &
      'got ', got, ', expected ', want, ', rel.err ', err, ' > tol ', t
    CALL record(label, err .LE. t, msg)
  END SUBROUTINE check_real

  !> Assert that two real arrays agree; the error is normalized by the array scale.
  SUBROUTINE check_real_arr(label, got, want, tol)
    CHARACTER(len=*), INTENT(IN) :: label
    REAL(xp), DIMENSION(:), INTENT(IN) :: got, want
    REAL(xp),               INTENT(IN), OPTIONAL :: tol
    REAL(xp) :: t, err, scal
    CHARACTER(len=160) :: msg
    t = tol_of(1.0e3_xp); IF (PRESENT(tol)) t = tol
    IF (SIZE(got) .NE. SIZE(want)) THEN
      CALL record(label, .FALSE., 'size mismatch')
      RETURN
    ENDIF
    scal = MAX(1._xp, MAXVAL(ABS(want)))
    err  = MAXVAL(ABS(got-want))/scal
    WRITE(msg,'(A,ES9.2,A,ES9.2)') 'max rel.err ', err, ' > tol ', t
    CALL record(label, err .LE. t, msg)
  END SUBROUTINE check_real_arr

  !> Assert that two complex numbers agree within tol.
  SUBROUTINE check_cmpx(label, got, want, tol)
    CHARACTER(len=*), INTENT(IN) :: label
    COMPLEX(xp),      INTENT(IN) :: got, want
    REAL(xp),         INTENT(IN), OPTIONAL :: tol
    REAL(xp) :: t, err
    CHARACTER(len=160) :: msg
    t = tol_of(1.0e3_xp); IF (PRESENT(tol)) t = tol
    err = ABS(got-want)/MAX(1._xp, ABS(want))
    WRITE(msg,'(A,2ES14.7,A,2ES14.7,A,ES9.2,A,ES9.2)') &
      'got ', got, ', expected ', want, ', rel.err ', err, ' > tol ', t
    CALL record(label, err .LE. t, msg)
  END SUBROUTINE check_cmpx

  !> Assert that two complex arrays agree; the error is normalized by the array scale.
  SUBROUTINE check_cmpx_arr(label, got, want, tol)
    CHARACTER(len=*), INTENT(IN) :: label
    COMPLEX(xp), DIMENSION(:), INTENT(IN) :: got, want
    REAL(xp),                  INTENT(IN), OPTIONAL :: tol
    REAL(xp) :: t, err, scal
    CHARACTER(len=160) :: msg
    t = tol_of(1.0e3_xp); IF (PRESENT(tol)) t = tol
    IF (SIZE(got) .NE. SIZE(want)) THEN
      CALL record(label, .FALSE., 'size mismatch')
      RETURN
    ENDIF
    scal = MAX(1._xp, MAXVAL(ABS(want)))
    err  = MAXVAL(ABS(got-want))/scal
    WRITE(msg,'(A,ES9.2,A,ES9.2)') 'max rel.err ', err, ' > tol ', t
    CALL record(label, err .LE. t, msg)
  END SUBROUTINE check_cmpx_arr

END MODULE test_harness
