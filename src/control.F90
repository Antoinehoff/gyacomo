SUBROUTINE control
  !   Control the run

  use basic,          ONLY: str,daytim,speak,basic_data,&
                        nlend,step,increase_step,increase_time,increase_cstep,&
                        chrono_runt,chrono_step, chrono_diag, chrono_ExBs,&
                        start_chrono, stop_chrono, day_and_time_str, show_title
  use prec_const,     ONLY: xp, stdout
  USE parallel,       ONLY: ppinit
  USE initial,        ONLY: initialize
  USE mpi,            ONLY: MPI_COMM_WORLD
  USE diagnostics,    ONLY: diagnose
  USE ExB_shear_flow, ONLY: Update_ExB_shear_flow
  IMPLICIT NONE
  REAL(xp) :: t_init_diag_0, t_init_diag_1
  INTEGER  :: ierr
  ! start the chronometer for the total runtime
  CALL start_chrono(chrono_runt)
  !________________________________________________________________________________
  !           ]   1.   Prologue

  !                   1.0     Initialize the parallel environment
  CALL ppinit
  CALL speak('] MPI initialized',2)
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)

  !                   1.1     Title and start time
  CALL show_title
  CALL speak('Start time: '//day_and_time_str(),0)
  
  !                   1.2     Define data specific to run
  CALL speak('Load basic data [',2)
  CALL basic_data
  CALL speak('] basic data loaded.',2)
  
  !                   1.3   Read input parameters from input file
  CALL speak('Read input parameters [',2)
  CALL readinputs
  CALL speak('] input parameters read',2)

  !                   1.4     Set auxiliary values (allocate arrays, set grid, ...)
  CALL speak('Setting auxiliary values [',2)
  CALL auxval
  CALL speak('] auxval set ',2)
  
  !                   1.5     Initial conditions
  CALL speak( 'Create initial state [',2)
  CALL initialize
  ! CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  CALL speak('] initial state created',2)
  
  !                   1.6     Initial diagnostics
  CALL speak( 'Initial diagnostics [',2)
  CALL cpu_time(t_init_diag_0) ! Measure the time of the init diag
  CALL diagnose(0)
  CALL cpu_time(t_init_diag_1)
  ! CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  CALL speak('] initial diagnostics done',2)
  CALL speak('('//str(t_init_diag_1-t_init_diag_0)//'[s])',2)

  CALL FLUSH(stdout)
  CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  !________________________________________________________________________________

  CALL speak( 'Time integration loop [',2)
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

    ! Update the ExB shear flow for the next step
    ! This call includes :
    !  - the ExB shear value (s(ky)) update for the next time step
    !  - the kx grid update
    !  - the ExB NL correction factor update (exp(+/- ixkySdts))
    !  - (optional) the kernel, poisson op. and ampere op update
    CALL start_chrono(chrono_ExBs)
      CALL Update_ExB_shear_flow(-1)
    CALL stop_chrono(chrono_ExBs)

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

  CALL speak('] time integration done',2)
  !________________________________________________________________________________
  !              9.   Epilogue
  ! Stop total run chronometer (must be done before the last diagnostic)
  CALL stop_chrono(chrono_runt)
  ! last diagnostic
  CALL diagnose(-1)
  ! end the run
  CALL endrun
  ! display final time
  ! CALL daytim('Run terminated at ')
  CALL speak('Finish time: '//day_and_time_str(),0)
  ! close mpi environement
  CALL ppexit

END SUBROUTINE control
