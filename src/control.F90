SUBROUTINE control
  !   Control the run

  use basic

  use prec_const
  IMPLICIT NONE

  !________________________________________________________________________________
  !              1.   Prologue
  !                   1.1     Initialize the parallel environment
  WRITE(*,*) 'Initialize MPI...'
  CALL ppinit
  WRITE(*,*) '...MPI initialized'


  call daytim('Start at ')
  
  !                   1.2     Define data specific to run
  WRITE(*,*) 'Load basic data...'
  CALL basic_data
  WRITE(*,*) '...basic data loaded.'

  
  !                   1.3   Read input parameters from input file
  WRITE(*,*) 'Read input parameters...'
  CALL readinputs
  WRITE(*,*) '...input parameters read'

  !                   1.4     Set auxiliary values (allocate arrays, set laplace operator, ...)
  WRITE(*,*) 'Calculate auxval...'
  CALL auxval
  WRITE(*,*) '...auxval calculated'
  !
  !                   1.5     Initial conditions
  WRITE(*,*) 'Create initial state...'
     CALL inital
  WRITE(*,*) '...initial state created'

  !                   1.6     Initial diagnostics
  WRITE(*,*) 'Initial diagnostics...'
  CALL diagnose(0)
  WRITE(*,*) '...initial diagnostics done'
  !
  CALL FLUSH(stdout)
  
  !________________________________________________________________________________
  !              2.   Main loop
  DO
     step = step+1
     cstep = cstep+1
     CALL stepon
     time = time + dt

     CALL tesend
     CALL diagnose(step)
     IF( nlend ) EXIT ! exit do loop

     ! CALL write_restart ! if want to write a restart file every so often (in case of crash)
  END DO
  !________________________________________________________________________________
  !              9.   Epilogue

  CALL diagnose(-1)
  CALL endrun

  CALL daytim('Done at ')
  CALL ppexit

END SUBROUTINE control
