SUBROUTINE ppinit
  !   Parallel environment

  USE basic
  use prec_const
  IMPLICIT NONE

  INTEGER :: version_prov=-1
  ! Variables for cartesian domain decomposition
  INTEGER, PARAMETER :: ndims=2 ! p and kr
  INTEGER, DIMENSION(ndims) :: dims=0, coords=0, coords_L=0, coords_R=0
  LOGICAL :: periods(ndims) = .FALSE., reorder=.FALSE.
  CHARACTER(len=32) :: str
  INTEGER :: nargs, i, l
  INTEGER :: nghb_L, nghb_R ! Left and right neighbors along p
  INTEGER :: source

  CALL MPI_INIT(ierr)
  ! CALL MPI_INIT_THREAD(MPI_THREAD_SINGLE,version_prov,ierr)
  ! CALL MPI_INIT_THREAD(MPI_THREAD_FUNNELED,version_prov,ierr)

  CALL MPI_COMM_RANK (MPI_COMM_WORLD,     my_id, ierr)
  CALL MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)
  
  ! nargs = COMMAND_ARGUMENT_COUNT()
  ! IF( nargs .NE. 0 .AND. nargs .NE. ndims ) THEN
  !    IF(my_id .EQ. 0) WRITE(*, '(a,i4,a)') 'Number of arguments not equal to NDIMS =', ndims, '!'
  !    CALL MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
  ! END IF
  ! !
  ! IF( nargs .NE. 0 ) THEN
  !    DO i=1,nargs
  !       CALL GET_COMMAND_ARGUMENT(i, str, l, ierr)
  !       READ(str(1:l),'(i3)')  dims(i)
  !       ncp = dims(1) ! Number of processes along p
  !       ncr = dims(2)  ! Number of processes along kr
  !    END DO
  !    IF( PRODUCT(dims) .NE. num_procs ) THEN
  !     IF(my_id .EQ. 0) WRITE(*, '(a,i4,a,i4)') 'Product of dims: ', PRODUCT(dims), " is not consistent WITH NPROCS=",num_procs
  !       CALL MPI_ABORT(MPI_COMM_WORLD, -2, ierr)
  !    END IF
  ! ELSE
  !    CALL MPI_DIMS_CREATE(num_procs, ndims, dims, ierr)
  ! END IF
  ! !
  ! !periodicity in p
  ! periods(1)=.FALSE.
  ! !periodicity in kr
  ! periods(2)=.FALSE.

  ! CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periods, reorder, comm0, ierr)

  ! CALL MPI_COMM_RANK(comm0, me_0,  ierr)

  ! CALL MPI_CART_COORDS(comm0,me_0,ndims,coords,ierr)

  ! DO i=0,num_procs-1
  !   CALL mpi_barrier(MPI_COMM_WORLD, ierr)
  !   IF (my_id .EQ. i) THEN
  !     WRITE(*,*) 'The coords of process ',me_0,' are: ', coords

  !     IF( coords(1) .GT. 0 ) THEN
  !       CALL MPI_CART_SHIFT(comm0, 0, -1, source , nghb_L , ierr)
  !       CALL MPI_CART_COORDS(comm0,nghb_L,ndims,coords_L,ierr)
  !       WRITE(*,*) coords_L,' are L from ', coords
  !     ELSE
  !       nghb_L = MPI_PROC_NULL
  !     ENDIF


  !     IF( coords(1) .LT. dims(1)-1 ) THEN
  !       CALL MPI_CART_SHIFT(comm0, 0, +1, source , nghb_R , ierr)
  !       CALL MPI_CART_COORDS(comm0,nghb_R,ndims,coords_R,ierr)
  !       WRITE(*,*) coords_R,' are R from ', coords
  !     ELSE
  !       nghb_R = MPI_PROC_NULL
  !     ENDIF
  !   ENDIF
  ! ENDDO

  ! !
  ! !  Partitions 2-dim cartesian topology of comm0 into 1-dim cartesian subgrids
  ! !
  ! CALL MPI_CART_SUB (comm0, (/.TRUE.,.FALSE./), commp, ierr)
  ! CALL MPI_CART_SUB (comm0, (/.FALSE.,.TRUE./), commr, ierr)

  ! CALL MPI_COMM_RANK(commp, me_p, ierr)
  ! CALL MPI_COMM_RANK(commr, me_r, ierr)

END SUBROUTINE ppinit
