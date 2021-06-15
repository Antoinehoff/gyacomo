MODULE utility

  USE basic

  use prec_const
  IMPLICIT NONE
  PUBLIC :: manual_2D_bcast

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
       !WRITE(*,'(a20,a25,i6.6,a20,i1)') str,' = Inf at timestep',cstep, ' and substep',updatetlevel
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


    !!!!! This is a manual way to do MPI_BCAST !!!!!!!!!!!
  SUBROUTINE manual_2D_bcast(field_)
    USE grid
    IMPLICIT NONE
    COMPLEX(dp), INTENT(INOUT) :: field_(ikrs:ikre,ikzs:ikze)
    COMPLEX(dp) :: buffer(ikrs:ikre,ikzs:ikze)
    INTEGER     :: i_, root, world_rank, world_size
    root = 0;
    CALL MPI_COMM_RANK(comm_p,world_rank,ierr)
    CALL MPI_COMM_SIZE(comm_p,world_size,ierr)
    IF (world_size .GT. 1) THEN
      !! Broadcast phi to the other processes on the same k range (communicator along p)
      IF (world_rank .EQ. root) THEN
        ! Fill the buffer
        DO ikr = ikrs,ikre
          DO ikz = ikzs,ikze
            buffer(ikr,ikz) = field_(ikr,ikz)
          ENDDO
        ENDDO
        ! Send it to all the other processes
        DO i_ = 0,num_procs_p-1
          IF (i_ .NE. world_rank) &
          CALL MPI_SEND(buffer, local_nkr * nkz , MPI_DOUBLE_COMPLEX, i_, 0, comm_p, ierr)
        ENDDO
      ELSE
        ! Recieve buffer from root
        CALL MPI_RECV(buffer, local_nkr * nkz , MPI_DOUBLE_COMPLEX, root, 0, comm_p, MPI_STATUS_IGNORE, ierr)
        ! Write it in phi
        DO ikr = ikrs,ikre
          DO ikz = ikzs,ikze
            field_(ikr,ikz) = buffer(ikr,ikz)
          ENDDO
        ENDDO
      ENDIF
    ENDIF
  END SUBROUTINE manual_2D_bcast

END MODULE utility
