MODULE utility
  IMPLICIT NONE
  PUBLIC :: is_nan, is_inf!. checkfield, checkelem
CONTAINS

  FUNCTION is_nan(x,str) RESULT(isn)
    USE basic,            ONLY: cstep
    USE time_integration, ONLY: updatetlevel
    USE prec_const,       ONLY: xp, stdout
    IMPLICIT NONE

    real, INTENT(IN) :: x
    CHARACTER(LEN=*), INTENT(IN) :: str
    LOGICAL :: isn

    isn=.FALSE.
    IF(x .NE. x) THEN
       isn=.TRUE.
    END IF
    IF((isn).AND.(str.NE.'')) THEN
       WRITE(*,'(a20,a25,i6.6,a20,i1)') str,' = NaN  at timestep',cstep, ' and substep',updatetlevel
       CALL FLUSH(stdout)
    END IF
  END FUNCTION is_nan


  FUNCTION is_inf(x,str) RESULT(isi)
    USE prec_const,       ONLY: xp, stdout
    IMPLICIT NONE

    real, INTENT(IN) :: x
    CHARACTER(LEN=*), INTENT(IN) :: str
    LOGICAL :: isi

    isi=.FALSE.
    IF(x+1.0== x) THEN
       isi=.TRUE.
    END IF
    IF((isi).AND.(str.NE.'')) THEN
       !WRITE(*,'(a20,a25,i6.6,a20,i1)') str,' = Inf at timestep',cstep, ' and substep',updatetlevel
       CALL FLUSH(stdout)
    END IF

  END FUNCTION is_inf

  ! FUNCTION checkfield(n1,n2,n3,field,str) RESULT(mlend)
  !   use prec_const, ONLY: xp
  !   IMPLICIT NONE
  !   !! BUG found (or feature?)
  !   ! if one put the commented first line (commented) instead of the second one,
  !   ! no error will be risen by the compiler even if the rank of the array is not matching (should be 3D!)
  !   ! COMPLEX(xp), DIMENSION(ikys:ikye,ikxs:ikxe), INTENT(IN) :: field
  !   INTEGER, INTENT(in) :: n1,n2,n3
  !   COMPLEX(xp), DIMENSION(n1,n2,n3), INTENT(IN) :: field
  !   CHARACTER(LEN=*), INTENT(IN) :: str
  !   LOGICAL :: mlend
  !   COMPLEX(xp) :: sumfield

  !   sumfield=SUM(field)

  !   mlend= is_nan( REAL(sumfield),str).OR.is_inf( REAL(sumfield),str) &
  !     .OR. is_nan(AIMAG(sumfield),str).OR.is_inf(AIMAG(sumfield),str)
  ! END FUNCTION checkfield

  ! FUNCTION checkelem(elem,str) RESULT(mlend)
  !   use prec_const, ONLY: xp
  !   IMPLICIT NONE
  !   COMPLEX(xp), INTENT(IN) :: elem
  !   CHARACTER(LEN=*), INTENT(IN) :: str
  !   LOGICAL :: mlend

  !   mlend= is_nan( REAL(elem),str).OR.is_inf( REAL(elem),str) &
  !     .OR. is_nan(AIMAG(elem),str).OR.is_inf(AIMAG(elem),str)
  ! END FUNCTION checkelem
END MODULE utility
