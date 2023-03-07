SUBROUTINE control
  !   Control the run

  use basic,      ONLY: str,daytim,speak,basic_data,start,t0_step,t1_step,tc_step,&
                        nlend,step,increase_step,increase_time,increase_cstep
  use prec_const, ONLY: dp, stdout
  USE parallel,   ONLY: ppinit
  USE mpi
  IMPLICIT NONE
  REAL(dp) :: t_init_diag_0, t_init_diag_1
  INTEGER  :: ierr
  CALL cpu_time(start)
  !________________________________________________________________________________
  !              1.   Prologue
  !                   1.1     Initialize the parallel environment
  CALL ppinit
  CALL speak('MPI initialized')


  CALL daytim('Start at ')

  !                   1.2     Define data specific to run
  CALL speak( 'Load basic data...')
  CALL basic_data
  ! CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  CALL speak('...basic data loaded.')

  !                   1.3   Read input parameters from input file
  CALL speak('Read input parameters...')
  CALL readinputs
  ! CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  CALL speak('...input parameters read')

  !                   1.4     Set auxiliary values (allocate arrays, set grid, ...)
  CALL speak('Calculate auxval...')
  CALL auxval
  ! CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  CALL speak('...auxval calculated')

  !                   1.5     Initial conditions
  CALL speak( 'Create initial state...')
  CALL inital
  ! CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  CALL speak('...initial state created')

  !                   1.6     Initial diagnostics
  CALL speak( 'Initial diagnostics...')
  CALL cpu_time(t_init_diag_0) ! Measure the time of the init diag
  CALL diagnose(0)
  CALL cpu_time(t_init_diag_1)
  ! CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  CALL speak('...initial diagnostics done')
  CALL speak('('//str(t_init_diag_1-t_init_diag_0)//'[s])')

  CALL FLUSH(stdout)
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  !________________________________________________________________________________

  CALL speak( 'Time integration loop..')
  !________________________________________________________________________________
  !              2.   Main loop
  DO
     CALL cpu_time(t0_step) ! Measuring time

     CALL increase_step
     CALL increase_cstep
     CALL stepon

     CALL increase_time

     CALL tesend
     IF( nlend ) EXIT ! exit do loop

     CALL diagnose(step)

    CALL cpu_time(t1_step);
    tc_step = tc_step + (t1_step - t0_step)

  END DO

  CALL speak('...time integration done')
  !________________________________________________________________________________
  !              9.   Epilogue

  CALL diagnose(-1)
  CALL endrun

  CALL daytim('Done at ')

  CALL ppexit

END SUBROUTINE control
