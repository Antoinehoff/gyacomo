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

  CALL MPI_INIT(ierr)

  CALL MPI_COMM_RANK (MPI_COMM_WORLD,     my_id, ierr)
  CALL MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)

  nargs = COMMAND_ARGUMENT_COUNT()
  IF( nargs .NE. 0 .AND. nargs .NE. ndims ) THEN
     IF(my_id .EQ. 0) WRITE(*, '(a,i4,a)') 'Number of arguments not equal to NDIMS =', ndims, '!'
     CALL MPI_ABORT(MPI_COMM_WORLD, -1, ierr)
  END IF
  !
  IF( nargs .NE. 0 ) THEN
     DO i=1,nargs
        CALL GET_COMMAND_ARGUMENT(i, str, l, ierr)
        READ(str(1:l),'(i3)')  dims(i)
     END DO
     IF( PRODUCT(dims) .NE. num_procs ) THEN
      IF(my_id .EQ. 0) WRITE(*, '(a,i4,a,i4)') 'Product of dims: ', PRODUCT(dims), " is not consistent WITH NPROCS=",num_procs
        CALL MPI_ABORT(MPI_COMM_WORLD, -2, ierr)
     END IF
  ELSE
    !  CALL MPI_DIMS_CREATE(num_procs, ndims, dims, ierr)
    dims(1) = 1
    dims(2) = num_procs
  END IF

  num_procs_p = dims(1) ! Number of processes along p
  num_procs_kr = dims(2) ! Number of processes along kr

  !
  !periodicity in p
  periods(1)=.FALSE.
  !periodicity in kr
  periods(2)=.FALSE.

  CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periods, reorder, comm0, ierr)

  CALL MPI_COMM_RANK(comm0, rank_0,  ierr)

  CALL MPI_CART_COORDS(comm0,rank_0,ndims,coords,ierr)

  !
  !  Partitions 2-dim cartesian topology of comm0 into 1-dim cartesian subgrids
  !
  CALL MPI_CART_SUB (comm0, (/.TRUE.,.FALSE./), commp, ierr)
  CALL MPI_CART_SUB (comm0, (/.FALSE.,.TRUE./), commr, ierr)
  ! Find id inside the sub communicators
  CALL MPI_COMM_RANK(commp, rank_p, ierr)
  CALL MPI_COMM_RANK(commr, rank_r, ierr)
  ! Find neighbours
  CALL MPI_CART_SHIFT(comm0, 0, 1, nbr_L, nbr_R, ierr)
  CALL MPI_CART_SHIFT(comm0, 1, 1, nbr_B, nbr_T, ierr)

END SUBROUTINE ppinit
