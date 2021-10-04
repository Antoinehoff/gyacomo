SUBROUTINE tesend
  !   Test for run completion

  USE basic

  use prec_const
  IMPLICIT NONE
  LOGICAL :: mlend, mlexist
  REAL    :: tnow
  INTEGER :: ncheck_stop = 100
  CHARACTER(len=*), PARAMETER :: stop_file = 'mystop'

  !________________________________________________________________________________
  !                   1.  Some processors had set nlend
  CALL mpi_allreduce(nlend, mlend, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
  IF( mlend ) THEN
    nlend   = .TRUE.
    crashed = .TRUE.
    IF (my_id .EQ. 0) WRITE(*,'(/a)') 'rhs are NaN/Inf'
    IF (my_id .EQ. 0) WRITE(*,*) 'Run terminated at cstep=',cstep
    RETURN
  END IF

  !________________________________________________________________________________
  !                   2.  Test on NRUN
  nlend = step .GT. nrun
  IF ( nlend ) THEN
     WRITE(*,'(/a)') 'NRUN steps done'
     RETURN
  END IF


  !________________________________________________________________________________
  !                   3.  Test on TMAX
  nlend = time .GT. tmax
  IF ( nlend ) THEN
     IF (my_id .EQ. 0) WRITE(*,'(/a)') 'TMAX reached'
     RETURN
  END IF
  !
  !

  !________________________________________________________________________________
  !                   4.  Test on run time
  CALL cpu_time(tnow)
  mlend = (1.2*(tnow-start)) .GT. maxruntime


  CALL mpi_allreduce(mlend, nlend, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
  IF ( nlend ) THEN
     IF(my_id.EQ.0) WRITE(*,'(/a)') 'Max run time reached'
     RETURN
  END IF
  !________________________________________________________________________________
  !                   5.  NRUN modified through "stop file"
  !
  IF( (my_id .EQ. 0) .AND. (MOD(cstep, ncheck_stop) == 0) ) THEN
     INQUIRE(file=stop_file, exist=mlexist)
     IF( mlexist ) THEN
        OPEN(lu_stop, file=stop_file)
        mlend = mlexist ! Send stop status asa the file exists
        WRITE(*,'(/a,i4,a)') 'Stop file found -> finishing..'
        CLOSE(lu_stop, status='delete')
     END IF
  END IF
  CALL mpi_allreduce(mlend, nlend, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
  !
  RETURN
  !
END SUBROUTINE tesend
