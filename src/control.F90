SUBROUTINE control
  !   Control the run

  use basic,      ONLY: str,daytim,speak,basic_data,&
                        nlend,step,increase_step,increase_time,increase_cstep,&
                        chrono_runt,chrono_step, chrono_diag, start_chrono, stop_chrono
  use prec_const,  ONLY: xp, stdout
  USE parallel,    ONLY: ppinit
  USE initial,     ONLY: initialize
  USE mpi,         ONLY: MPI_COMM_WORLD
  USE diagnostics, ONLY: diagnose
  IMPLICIT NONE
  REAL(xp) :: t_init_diag_0, t_init_diag_1
  INTEGER  :: ierr
  ! start the chronometer for the total runtime
  CALL start_chrono(chrono_runt)
  !________________________________________________________________________________
  !              1.   Prologue
  !                   1.1     Initialize the parallel environment
  CALL ppinit
  CALL speak('MPI initialized')
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)

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
  CALL initialize
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
    CALL start_chrono(chrono_step) ! Measuring time per step
    
    ! Test if the stopping requirements are met (update nlend)
    CALL tesend
    IF( nlend ) EXIT ! exit do loop

    ! Increment steps and csteps (private in basic module)
    CALL increase_step
    CALL increase_cstep

    ! Do a full RK step (e.g. 4 substeps for ERK4)
    CALL stepon

    ! Increment time (private in basic module)
    CALL increase_time

    ! Periodic diagnostics
    CALL start_chrono(chrono_diag)
      CALL diagnose(step)
    CALL stop_chrono(chrono_diag)

    CALL stop_chrono(chrono_step)

  END DO

  CALL speak('...time integration done')
  !________________________________________________________________________________
  !              9.   Epilogue
  ! Stop total run chronometer (must be done before the last diagnostic)
  CALL stop_chrono(chrono_runt)
  ! last diagnostic
  CALL diagnose(-1)
  ! end the run
  CALL endrun
  ! display final time
  CALL daytim('Done at ')
  ! close mpi environement
  CALL ppexit

END SUBROUTINE control
