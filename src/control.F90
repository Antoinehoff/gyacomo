SUBROUTINE control
  !   Control the run

  use basic
  use prec_const
  IMPLICIT NONE
  REAL(dp) :: t_init_diag_0, t_init_diag_1

  CALL cpu_time(start)
  !________________________________________________________________________________
  !              1.   Prologue
  !                   1.1     Initialize the parallel environment
  CALL ppinit
  IF (my_id .EQ. 0) WRITE(*,'(a/)') 'MPI initialized'


  CALL daytim('Start at ')

  !                   1.2     Define data specific to run
  IF (my_id .EQ. 0) WRITE(*,*) 'Load basic data...'
  CALL basic_data
  ! CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  IF (my_id .EQ. 0) WRITE(*,'(a/)') '...basic data loaded.'

  !                   1.3   Read input parameters from input file
  IF (my_id .EQ. 0) WRITE(*,*) 'Read input parameters...'
  CALL readinputs
  ! CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  IF (my_id .EQ. 0) WRITE(*,'(a/)') '...input parameters read'

  !                   1.4     Set auxiliary values (allocate arrays, set grid, ...)
  IF (my_id .EQ. 0) WRITE(*,*) 'Calculate auxval...'
  CALL auxval
  ! CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  IF (my_id .EQ. 0) WRITE(*,'(a/)') '...auxval calculated'

  !                   1.5     Initial conditions
  IF (my_id .EQ. 0) WRITE(*,*) 'Create initial state...'
  CALL inital
  ! CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  IF (my_id .EQ. 0) WRITE(*,'(a/)') '...initial state created'

  !                   1.6     Initial diagnostics
  IF (my_id .EQ. 0) WRITE(*,*) 'Initial diagnostics...'
  CALL cpu_time(t_init_diag_0) ! Measure the time of the init diag
  CALL diagnose(0)
  CALL cpu_time(t_init_diag_1)
  ! CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  IF (my_id .EQ. 0) THEN
    WRITE(*,'(a)') '...initial diagnostics done'
    WRITE(*,'(a,F6.3,a/)') '(',t_init_diag_1-t_init_diag_0,'[s])'
  ENDIF

  CALL FLUSH(stdout)
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  !________________________________________________________________________________

  IF (my_id .EQ. 0) WRITE(*,*) 'Time integration loop..'
  !________________________________________________________________________________
  !              2.   Main loop
  DO
     CALL cpu_time(t0_step) ! Measuring time

     step  = step  + 1
     cstep = cstep + 1
     CALL stepon

     time  = time  + dt

     CALL tesend
     IF( nlend ) EXIT ! exit do loop

     CALL diagnose(step)

    CALL cpu_time(t1_step);
    tc_step = tc_step + (t1_step - t0_step)

  END DO

  IF (my_id .EQ. 0) WRITE(*,'(a/)') '...time integration done'
  !________________________________________________________________________________
  !              9.   Epilogue

  CALL diagnose(-1)
  CALL endrun

  IF (my_id .EQ. 0) CALL daytim('Done at ')

  CALL ppexit

END SUBROUTINE control
