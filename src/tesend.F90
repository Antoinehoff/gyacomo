SUBROUTINE tesend
  !   Test for run completion

  USE basic
  USE mpi
  USE parallel, ONLY: my_id
  IMPLICIT NONE
  LOGICAL :: mlend, mlexist
  REAL    :: tnow
  INTEGER :: ncheck_stop = 100, ierr
  CHARACTER(len=*), PARAMETER :: stop_file = 'mystop'

  !________________________________________________________________________________
  !                   1.  Some processors had set nlend
  CALL mpi_allreduce(nlend, mlend, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
  IF( mlend ) THEN
    nlend   = .TRUE.
    crashed = .TRUE.
    CALL speak('rhs are NaN/Inf')
    CALL speak('Run terminated at cstep='//str(cstep))
    RETURN
  END IF

  !________________________________________________________________________________
  !                   2.  Test on NRUN
  nlend = step .GE. nrun
  IF ( nlend ) THEN
     CALL speak('NRUN steps done')
     RETURN
  END IF


  !________________________________________________________________________________
  !                   3.  Test on TMAX
  nlend = time .GE. tmax
  IF ( nlend ) THEN
     CALL speak('TMAX reached')
     RETURN
  END IF
  !
  !

  !________________________________________________________________________________
  !                   4.  Test on run time
  CALL cpu_time(tnow)
  mlend = (1.1*(tnow-chrono_runt%tstart)) .GT. maxruntime


  CALL mpi_allreduce(mlend, nlend, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
  IF ( nlend ) THEN
     CALL speak('Max run time reached')
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
        CALL speak('Stop file found -> finishing..')
        CLOSE(lu_stop, status='delete')
     END IF
  END IF
  CALL mpi_allreduce(mlend, nlend, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
  !
  RETURN
  !
END SUBROUTINE tesend
