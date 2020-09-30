MODULE utility

  USE basic

  use prec_const
  IMPLICIT NONE

CONTAINS

  FUNCTION is_nan(x,str) RESULT(isn)

    USE time_integration

    use prec_const
    IMPLICIT NONE

    real(dp), INTENT(IN) :: x
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

    USE time_integration

    use prec_const
    IMPLICIT NONE

    real(dp), INTENT(IN) :: x
    CHARACTER(LEN=*), INTENT(IN) :: str
    LOGICAL :: isi

    isi=.FALSE.
    IF(x+1.0== x) THEN
       isi=.TRUE.
    END IF
    IF((isi).AND.(str.NE.'')) THEN
       WRITE(*,'(a20,a25,i6.6,a20,i1)') str,' = Inf at timestep',cstep, ' and substep',updatetlevel
       CALL FLUSH(stdout)
    END IF

  END FUNCTION is_inf

  FUNCTION checkfield(field,str) RESULT(mlend)

    USE grid

    use prec_const
    IMPLICIT NONE

    COMPLEX(dp), DIMENSION(ikrs:ikre,ikzs:ikze), INTENT(IN) :: field
    CHARACTER(LEN=*), INTENT(IN) :: str
    LOGICAL :: mlend
    COMPLEX(dp) :: sumfield

    sumfield=SUM(field)

    mlend= is_nan( REAL(sumfield),str).OR.is_inf( REAL(sumfield),str) &
      .OR. is_nan(AIMAG(sumfield),str).OR.is_inf(AIMAG(sumfield),str)
  END FUNCTION checkfield

END MODULE utility


