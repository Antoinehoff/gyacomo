SUBROUTINE control
  !   Control the run

  use basic
  use prec_const
  IMPLICIT NONE

  call cpu_time(start)
  !________________________________________________________________________________
  !              1.   Prologue
  !                   1.1     Initialize the parallel environment
  WRITE(*,*) 'Initialize MPI...'
  CALL ppinit
  WRITE(*,'(a/)') '...MPI initialized'


  call daytim('Start at ')

  !                   1.2     Define data specific to run
  WRITE(*,*) 'Load basic data...'
  CALL basic_data
  WRITE(*,'(a/)') '...basic data loaded.'


  !                   1.3   Read input parameters from input file
  WRITE(*,*) 'Read input parameters...'
  CALL readinputs
  WRITE(*,'(a/)') '...input parameters read'

  !                   1.4     Set auxiliary values (allocate arrays, set laplace operator, ...)
  WRITE(*,*) 'Calculate auxval...'
  CALL auxval
  WRITE(*,'(a/)') '...auxval calculated'
  !
  !                   1.5     Initial conditions
  WRITE(*,*) 'Create initial state...'
  CALL inital
  WRITE(*,'(a/)') '...initial state created'

  !                   1.6     Initial diagnostics
  WRITE(*,*) 'Initial diagnostics...'
  CALL diagnose(0)
  WRITE(*,'(a/)') '...initial diagnostics done'
  !
  CALL FLUSH(stdout)
  !________________________________________________________________________________

  WRITE(*,*) 'Time integration loop..'
  !________________________________________________________________________________
  !              2.   Main loop
  DO
     step  = step  + 1
     cstep = cstep + 1

     CALL stepon
     time  = time  + dt

     CALL tesend
     IF( nlend ) EXIT ! exit do loop

     CALL diagnose(step)

  END DO
  WRITE(*,'(a/)') '...time integration done'
  !________________________________________________________________________________
  !              9.   Epilogue

  CALL diagnose(-1)
  CALL endrun

  CALL daytim('Done at ')
  CALL ppexit

END SUBROUTINE control
