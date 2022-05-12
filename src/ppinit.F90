SUBROUTINE ppinit
  !   Parallel environment

  USE basic
  use prec_const
  IMPLICIT NONE

  INTEGER :: version_prov=-1
  ! Variables for cartesian domain decomposition
  INTEGER, PARAMETER :: ndims=3 ! p, kx and z
  INTEGER, DIMENSION(ndims) :: dims=0, coords=0, coords_L=0, coords_R=0
  LOGICAL :: periods(ndims) = .FALSE., reorder=.FALSE.
  CHARACTER(len=32) :: str
  INTEGER :: nargs, i, l

  CALL MPI_INIT(ierr)

  CALL MPI_COMM_RANK (MPI_COMM_WORLD,     my_id, ierr)
  CALL MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)

  nargs = COMMAND_ARGUMENT_COUNT()
  !
  IF( nargs .GT. 1 ) THEN
     DO i=1,ndims
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
    dims(3) = 1
  END IF

  num_procs_p  = dims(1) ! Number of processes along p
  num_procs_ky = dims(2) ! Number of processes along kx
  num_procs_z  = dims(3) ! Number of processes along z

  !
  !periodicity in p
  periods(1)=.FALSE.
  !periodicity in ky
  periods(2)=.FALSE.
  !periodicity in z
  periods(3)=.TRUE.

  CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periods, reorder, comm0, ierr)

  CALL MPI_COMM_RANK(comm0, rank_0,  ierr)

  CALL MPI_CART_COORDS(comm0,rank_0,ndims,coords,ierr)

  !
  !  Partitions 3-dim cartesian topology of comm0 into 1-dim cartesian subgrids
  !
  CALL MPI_CART_SUB (comm0, (/.TRUE.,.FALSE.,.FALSE./),  comm_p, ierr)
  CALL MPI_CART_SUB (comm0, (/.FALSE.,.TRUE.,.FALSE./), comm_ky, ierr)
  CALL MPI_CART_SUB (comm0, (/.FALSE.,.FALSE.,.TRUE./),  comm_z, ierr)
  ! Find id inside the 1d-sub communicators
  CALL MPI_COMM_RANK(comm_p,  rank_p,  ierr)
  CALL MPI_COMM_RANK(comm_ky, rank_ky, ierr)
  CALL MPI_COMM_RANK(comm_z,  rank_z,  ierr)

  ! 2D communicator
  CALL MPI_CART_SUB (comm0, (/.TRUE.,.FALSE.,.TRUE./),  comm_pz,  ierr)
  CALL MPI_CART_SUB (comm0, (/.FALSE.,.TRUE.,.TRUE./),  comm_kyz, ierr)

  ! Find neighbours
  CALL MPI_CART_SHIFT(comm0, 0, 1, nbr_L, nbr_R, ierr) !left   right neighbours
  CALL MPI_CART_SHIFT(comm0, 1, 1, nbr_B, nbr_T, ierr) !bottom top   neighbours
  CALL MPI_CART_SHIFT(comm0, 2, 1, nbr_D, nbr_U, ierr) !down   up    neighbours

END SUBROUTINE ppinit
