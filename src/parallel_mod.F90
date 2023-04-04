MODULE parallel
  use prec_const, ONLY : xp, mpi_xp_c
  USE mpi
  IMPLICIT NONE
  ! Auxiliary variables
  INTEGER, PUBLIC, PROTECTED  :: comm0                 ! Default communicator with a topology
  INTEGER, PUBLIC, PROTECTED  :: group0                ! Default group with a topology
  INTEGER, PUBLIC, PROTECTED  :: rank_0                ! Ranks in comm0
  ! Communicators for 1-dim cartesian subgrids of comm0
  INTEGER, PUBLIC, PROTECTED  :: comm_p, comm_ky, comm_z
  INTEGER, PUBLIC, PROTECTED  :: rank_p, rank_ky, rank_z! Ranks
  INTEGER, PUBLIC, PROTECTED  :: comm_pz,  rank_pz      ! 2D comm for N_a(p,j,z) ouptut (mspfile)
  INTEGER, PUBLIC, PROTECTED  :: comm_kyz, rank_kyz     ! 2D comm for N_a(p,j,z) ouptut (mspfile)
  INTEGER, PUBLIC, PROTECTED  :: comm_ky0, rank_ky0     ! comm along ky with p=0
  INTEGER, PUBLIC, PROTECTED  :: comm_z0,  rank_z0      ! comm along z  with p=0
  INTEGER, PUBLIC, PROTECTED  :: group_ky0, group_z0
  INTEGER, PUBLIC, PROTECTED  :: ierr                  ! flag for MPI error
  INTEGER, PUBLIC, PROTECTED  :: my_id                 ! Rank in COMM_WORLD
  INTEGER, PUBLIC, PROTECTED  :: num_procs             ! number of MPI processes
  INTEGER, PUBLIC, PROTECTED  :: num_procs_p           ! Number of processes in p
  INTEGER, PUBLIC, PROTECTED  :: num_procs_ky          ! Number of processes in r
  INTEGER, PUBLIC, PROTECTED  :: num_procs_z           ! Number of processes in z
  INTEGER, PUBLIC, PROTECTED  :: num_procs_pz          ! Number of processes in pz comm
  INTEGER, PUBLIC, PROTECTED  :: num_procs_kyz         ! Number of processes in kyz comm
  INTEGER, PUBLIC, PROTECTED  :: nbr_L, nbr_R          ! Left and right neighbours (along p)
  INTEGER, PUBLIC, PROTECTED  :: nbr_T, nbr_B          ! Top and bottom neighbours (along kx)
  INTEGER, PUBLIC, PROTECTED  :: nbr_U, nbr_D          ! Upstream and downstream neighbours (along z)


  ! recieve and displacement counts for gatherv
  INTEGER, DIMENSION(:), ALLOCATABLE :: rcv_p, dsp_p
  INTEGER, DIMENSION(:), ALLOCATABLE :: rcv_y, dsp_y
  INTEGER, DIMENSION(:), ALLOCATABLE :: rcv_zy, dsp_zy
  INTEGER, DIMENSION(:), ALLOCATABLE :: rcv_zp,  dsp_zp
  INTEGER, DIMENSION(:), ALLOCATABLE :: rcv_yp,  dsp_yp
  INTEGER, DIMENSION(:), ALLOCATABLE :: rcv_zyp, dsp_zyp

  PUBLIC :: ppinit, manual_0D_bcast, manual_3D_bcast, init_parallel_var, &
            gather_xyz, gather_pjz, gather_pjxyz, exchange_ghosts_1D

CONTAINS

  SUBROUTINE ppinit
    ! Init the parallel environment
    IMPLICIT NONE
    ! Variables for cartesian domain decomposition
    INTEGER, PARAMETER :: ndims=3 ! p, kx and z
    INTEGER, DIMENSION(ndims) :: dims=0, coords=0
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
      dims(1) = 1
      dims(2) = num_procs
      dims(3) = 1
    END IF
    num_procs_p  = dims(1) ! Number of processes along p
    num_procs_ky = dims(2) ! Number of processes along kx
    num_procs_z  = dims(3) ! Number of processes along z
    !
    !periodiciyt in p
    periods(1)=.FALSE.
    !periodiciyt in ky
    periods(2)=.FALSE.
    !periodiciyt in z
    periods(3)=.TRUE.
    CALL MPI_CART_CREATE(MPI_COMM_WORLD, ndims, dims, periods, reorder, comm0, ierr)
    CALL MPI_COMM_GROUP(comm0,group0, ierr)
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
    ! Count the number of processes in 2D comms
    CALL MPI_COMM_SIZE(comm_pz, num_procs_pz, ierr)
    CALL MPI_COMM_SIZE(comm_kyz,num_procs_kyz,ierr)
    ! Find id inside the 1d-sub communicators
    CALL MPI_COMM_RANK(comm_pz,  rank_pz,   ierr)
    CALL MPI_COMM_RANK(comm_kyz, rank_kyz,  ierr)
    ! Find neighbours
    CALL MPI_CART_SHIFT(comm0, 0, 1, nbr_L, nbr_R, ierr) !left   right neighbours
    CALL MPI_CART_SHIFT(comm0, 1, 1, nbr_B, nbr_T, ierr) !bottom top   neighbours
    CALL MPI_CART_SHIFT(comm0, 2, 1, nbr_D, nbr_U, ierr) !down   up    neighbours
  END SUBROUTINE ppinit

  SUBROUTINE init_parallel_var(np_loc,np_tot,nky_loc,nky_tot,nz_loc)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: np_loc,np_tot,nky_loc,nky_tot,nz_loc
    INTEGER :: i_
    !     !! P reduction at constant x,y,z,j
    ! ALLOCATE(rcv_p(0:num_procs_p-1),dsp_p(0:num_procs_p-1)) !Displacement sizes for balance diagnostic
    ! ! all processes share their local number of points
    ! CALL MPI_ALLGATHER(np_loc,1,MPI_INTEGER,rcv_p,1,MPI_INTEGER,comm_p,ierr)
    ! ! the displacement array can be build from each np_loc as
    ! dsp_p(0)=0
    ! DO i_=1,num_procs_p-1
    !    dsp_p(i_) =dsp_p(i_-1) + rcv_p(i_-1)
    ! END DO
    ! !!!!!! XYZ gather variables
    ! !! Y reduction at constant x and z
    ! ! number of points recieved and displacement for the y reduction
    ! ALLOCATE(rcv_y(0:num_procs_ky-1),dsp_y(0:num_procs_ky-1)) !Displacement sizes for balance diagnostic
    ! ! all processes share their local number of points
    ! CALL MPI_ALLGATHER(nky_loc,1,MPI_INTEGER,rcv_y,1,MPI_INTEGER,comm_ky,ierr)
    ! ! the displacement array can be build from each np_loc as
    ! dsp_y(0)=0
    ! DO i_=1,num_procs_ky-1
    !    dsp_y(i_) =dsp_y(i_-1) + rcv_y(i_-1)
    ! END DO
    ! !! Z reduction for full slices of y data but constant x
    ! ! number of points recieved and displacement for the z reduction
    ! ALLOCATE(rcv_zy(0:num_procs_z-1),dsp_zy(0:num_procs_z-1)) !Displacement sizes for balance diagnostic
    ! ! all processes share their local number of points
    ! CALL MPI_ALLGATHER(nz_loc*nky_tot,1,MPI_INTEGER,rcv_zy,1,MPI_INTEGER,comm_z,ierr)
    ! ! the displacement array can be build from each np_loc as
    ! dsp_zy(0)=0
    ! DO i_=1,num_procs_z-1
    !    dsp_zy(i_) =dsp_zy(i_-1) + rcv_zy(i_-1)
    ! END DO
    ! !!!!! PJZ gather variables
    ! !! P reduction at constant j and z is already done in module GRID
    ! !! Z reduction for full slices of p data but constant j
    ! ! number of points recieved and displacement for the z reduction
    ! ALLOCATE(rcv_zp(0:num_procs_z-1),dsp_zp(0:num_procs_z-1)) !Displacement sizes for balance diagnostic
    ! ! all processes share their local number of points
    ! CALL MPI_ALLGATHER(nz_loc*np_tot,1,MPI_INTEGER,rcv_zp,1,MPI_INTEGER,comm_z,ierr)
    ! ! the displacement array can be build from each np_loc as
    ! dsp_zp(0)=0
    ! DO i_=1,num_procs_z-1
    !    dsp_zp(i_) =dsp_zp(i_-1) + rcv_zp(i_-1)
    ! END DO
    ! !!!!! PJXYZ gather variables
    ! !! Y reduction for full slices of p data but constant j
    ! ! number of points recieved and displacement for the y reduction
    ! ALLOCATE(rcv_yp(0:num_procs_ky-1),dsp_yp(0:num_procs_ky-1)) !Displacement sizes for balance diagnostic
    ! ! all processes share their local number of points
    ! CALL MPI_ALLGATHER(nky_loc*np_tot,1,MPI_INTEGER,rcv_yp,1,MPI_INTEGER,comm_ky,ierr)
    ! ! the displacement array can be build from each np_loc as
    ! dsp_yp(0)=0
    ! DO i_=1,num_procs_ky-1
    !    dsp_yp(i_) =dsp_yp(i_-1) + rcv_yp(i_-1)
    ! END DO
    ! !! Z reduction for full slices of py data but constant j
    ! ! number of points recieved and displacement for the z reduction
    ! ALLOCATE(rcv_zyp(0:num_procs_z-1),dsp_zyp(0:num_procs_z-1)) !Displacement sizes for balance diagnostic
    ! ! all processes share their local number of points
    ! CALL MPI_ALLGATHER(nz_loc*np_tot*nky_tot,1,MPI_INTEGER,rcv_zyp,1,MPI_INTEGER,comm_z,ierr)
    ! ! the displacement array can be build from each np_loc as
    ! dsp_zyp(0)=0
    ! DO i_=1,num_procs_z-1
    !    dsp_zyp(i_) =dsp_zyp(i_-1) + rcv_zyp(i_-1)
    ! END DO
    ! P reduction at constant x,y,z,j
    ALLOCATE(rcv_p(num_procs_p),dsp_p(num_procs_p)) !Displacement sizes for balance diagnostic
    CALL MPI_ALLGATHER(np_loc,1,MPI_INTEGER,rcv_p,1,MPI_INTEGER,comm_p,ierr)
    dsp_p(1)=0
    DO i_=2,num_procs_p
       dsp_p(i_) =dsp_p(i_-1) + rcv_p(i_-1)
    END DO
    !!!!!! XYZ gather variables
    ALLOCATE(rcv_y(num_procs_ky),dsp_y(num_procs_ky)) !Displacement sizes for balance diagnostic
    CALL MPI_ALLGATHER(nky_loc,1,MPI_INTEGER,rcv_y,1,MPI_INTEGER,comm_ky,ierr)
    dsp_y(1)=0
    DO i_=2,num_procs_ky
       dsp_y(i_) =dsp_y(i_-1) + rcv_y(i_-1)
    END DO
    !! Z reduction for full slices of y data but constant x
    ALLOCATE(rcv_zy(num_procs_z),dsp_zy(num_procs_z)) !Displacement sizes for balance diagnostic
    CALL MPI_ALLGATHER(nz_loc*nky_tot,1,MPI_INTEGER,rcv_zy,1,MPI_INTEGER,comm_z,ierr)
    dsp_zy(1)=0
    DO i_=2,num_procs_z
       dsp_zy(i_) =dsp_zy(i_-1) + rcv_zy(i_-1)
    END DO
    !!!!! PJZ gather variables
    ALLOCATE(rcv_zp(num_procs_z),dsp_zp(num_procs_z)) !Displacement sizes for balance diagnostic
    CALL MPI_ALLGATHER(nz_loc*np_tot,1,MPI_INTEGER,rcv_zp,1,MPI_INTEGER,comm_z,ierr)
    dsp_zp(1)=0
    DO i_=2,num_procs_z
       dsp_zp(i_) =dsp_zp(i_-1) + rcv_zp(i_-1)
    END DO
    !!!!! PJXYZ gather variables
    !! Y reduction for full slices of p data but constant j
    ALLOCATE(rcv_yp(num_procs_ky),dsp_yp(num_procs_ky)) !Displacement sizes for balance diagnostic
    CALL MPI_ALLGATHER(nky_loc*np_tot,1,MPI_INTEGER,rcv_yp,1,MPI_INTEGER,comm_ky,ierr)
    dsp_yp(1)=0
    DO i_=2,num_procs_ky
       dsp_yp(i_) =dsp_yp(i_-1) + rcv_yp(i_-1)
    END DO
    !! Z reduction for full slices of py data but constant j
    ALLOCATE(rcv_zyp(num_procs_z),dsp_zyp(num_procs_z)) !Displacement sizes for balance diagnostic
    CALL MPI_ALLGATHER(nz_loc*np_tot*nky_tot,1,MPI_INTEGER,rcv_zyp,1,MPI_INTEGER,comm_z,ierr)
    dsp_zyp(1)=0
    DO i_=2,num_procs_z
       dsp_zyp(i_) =dsp_zyp(i_-1) + rcv_zyp(i_-1)
    END DO
  END SUBROUTINE init_parallel_var

  SUBROUTINE parallel_ouptutinputs(fid)
    !
    !    Write the input parameters to the results_xx.h5 file
    !
    USE futils, ONLY: attach, creatd
    IMPLICIT NONE
    INTEGER, INTENT(in) :: fid
    CHARACTER(len=256)  :: str
    WRITE(str,'(a)') '/data/input/parallel'
    CALL creatd(fid, 0,(/0/),TRIM(str),'Parallel Input')
    CALL attach(fid, TRIM(str),       "Nproc",   num_procs)
    CALL attach(fid, TRIM(str),       "Np_p" , num_procs_p)
    CALL attach(fid, TRIM(str),       "Np_kx",num_procs_ky)
    CALL attach(fid, TRIM(str),        "Np_z", num_procs_z)
  END SUBROUTINE parallel_ouptutinputs

  !!!! Gather a field in spatial coordinates on rank 0 !!!!!
  SUBROUTINE gather_xyz(field_loc,field_tot,nky_loc,nky_tot,nkx_tot,nz_loc,nz_tot)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nky_loc,nky_tot,nkx_tot,nz_loc,nz_tot
    COMPLEX(xp), DIMENSION(:,:,:), INTENT(IN)  :: field_loc
    COMPLEX(xp), DIMENSION(:,:,:), INTENT(OUT) :: field_tot
    COMPLEX(xp), DIMENSION(nky_tot,nz_loc) :: buffer_yt_zl !full  y, local z
    COMPLEX(xp), DIMENSION(nky_tot,nz_tot) :: buffer_yt_zt !full  y, full  z
    COMPLEX(xp), DIMENSION(nky_loc):: buffer_yl_zc !local y, constant z
    COMPLEX(xp), DIMENSION(nky_tot):: buffer_yt_zc !full  y, constant z
    INTEGER :: snd_y, snd_z, root_p, root_z, root_ky, ix, iz

    snd_y  = nky_loc    ! Number of points to send along y (per z)
    snd_z  = nky_tot*nz_loc ! Number of points to send along z (full y)
    root_p = 0; root_z = 0; root_ky = 0
    IF(rank_p .EQ. root_p) THEN
      DO ix = 1,nkx_tot
        DO iz = 1,nz_loc
          ! fill a buffer to contain a slice of data at constant kx and z
          buffer_yl_zc(1:nky_loc) = field_loc(1:nky_loc,ix,iz)
          CALL MPI_GATHERV(buffer_yl_zc, snd_y,        mpi_xp_c, &
                           buffer_yt_zc, rcv_y, dsp_y, mpi_xp_c, &
                           root_ky, comm_ky, ierr)
          buffer_yt_zl(1:nky_tot,iz) = buffer_yt_zc(1:nky_tot)
        ENDDO
        ! send the full line on y contained by root_ky
        IF(rank_ky .EQ. root_ky) THEN
          CALL MPI_GATHERV(buffer_yt_zl, snd_z,          mpi_xp_c, &
                           buffer_yt_zt, rcv_zy, dsp_zy, mpi_xp_c, &
                           root_z, comm_z, ierr)
        ENDIF
        ! ID 0 (the one who output) rebuild the whole array
        IF(my_id .EQ. 0) &
          field_tot(1:nky_tot,ix,1:nz_tot) = buffer_yt_zt(1:nky_tot,1:nz_tot)
      ENDDO
    ENDIF
  END SUBROUTINE gather_xyz
  

  !!!!! Gather a field in kinetic + z coordinates on rank 0 !!!!!
  SUBROUTINE gather_pjz(field_loc,field_tot,na_tot,np_loc,np_tot,nj_tot,nz_loc,nz_tot)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: na_tot,np_loc,np_tot, nj_tot, nz_loc,nz_tot
    COMPLEX(xp), DIMENSION(:,:,:,:), INTENT(IN)  :: field_loc
    COMPLEX(xp), DIMENSION(:,:,:,:), INTENT(OUT) :: field_tot
    COMPLEX(xp), DIMENSION(np_tot,nz_loc) :: buffer_pt_zl !full  p, local z
    COMPLEX(xp), DIMENSION(np_tot,nz_tot) :: buffer_pt_zt !full  p, full  z
    COMPLEX(xp), DIMENSION(np_loc) :: buffer_pl_zc !local p, constant z
    COMPLEX(xp), DIMENSION(np_tot) :: buffer_pt_zc !full  p, constant z
    INTEGER :: snd_p, snd_z, root_p, root_z, root_ky, ij, iz, ia

    snd_p  = np_loc    ! Number of points to send along p (per z)
    snd_z  = np_tot*nz_loc ! Number of points to send along z (full p)

    root_p = 0; root_z = 0; root_ky = 0
    IF(rank_ky .EQ. root_ky) THEN
      DO ia= 1,na_tot
        DO ij = 1,nj_tot
          DO iz = 1,nz_loc
            ! fill a buffer to contain a slice of data at constant j and z
            buffer_pl_zc(1:np_loc) = field_loc(ia,1:np_loc,ij,iz)
            CALL MPI_GATHERV(buffer_pl_zc, snd_p,        mpi_xp_c, &
                            buffer_pt_zc, rcv_p, dsp_p, mpi_xp_c, &
                            root_p, comm_p, ierr)
            buffer_pt_zl(1:np_tot,iz) = buffer_pt_zc(1:np_tot)
          ENDDO
          ! send the full line on y contained by root_p
          IF(rank_p .EQ. root_p) THEN
            CALL MPI_GATHERV(buffer_pt_zl, snd_z,         mpi_xp_c, &
                            buffer_pt_zt, rcv_zp, dsp_zp, mpi_xp_c, &
                            root_z, comm_z, ierr)
          ENDIF
          ! ID 0 (the one who output) rebuild the whole array
          IF(my_id .EQ. 0) &
            field_tot(ia,1:np_tot,ij,1:nz_tot) = buffer_pt_zt(1:np_tot,1:nz_tot)
        ENDDO
      ENDDO
    ENDIF
  END SUBROUTINE gather_pjz

  !!!!! Gather a field in kinetic + spatial coordinates on rank 0 !!!!!
  !!!!! Gather a field in spatial coordinates on rank 0 !!!!!
  SUBROUTINE gather_pjxyz(field_loc,field_tot,na_tot,np_loc,np_tot,nj_tot,nky_loc,nky_tot,nkx_tot,nz_loc,nz_tot)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: na_tot,np_loc,np_tot,nj_tot,nky_loc,nky_tot,nkx_tot,nz_loc,nz_tot
    COMPLEX(xp), DIMENSION(:,:,:,:,:,:), INTENT(IN)  :: field_loc
    COMPLEX(xp), DIMENSION(:,:,:,:,:,:), INTENT(OUT) :: field_tot
    COMPLEX(xp), DIMENSION(np_tot,nky_tot,nz_loc) :: buffer_pt_yt_zl ! full p,     full y, local    z
    COMPLEX(xp), DIMENSION(np_tot,nky_tot,nz_tot) :: buffer_pt_yt_zt ! full p,     full y, full     z
    COMPLEX(xp), DIMENSION(np_tot,nky_loc) :: buffer_pt_yl_zc     ! full p,    local y, constant z
    COMPLEX(xp), DIMENSION(np_tot,nky_tot) :: buffer_pt_yt_zc     ! full p,     full y, constant z
    COMPLEX(xp), DIMENSION(np_loc)       :: buffer_pl_cy_zc     !local p, constant y, constant z
    COMPLEX(xp), DIMENSION(np_tot)       :: buffer_pt_cy_zc     ! full p, constant y, constant z
    INTEGER :: snd_p, snd_y, snd_z, root_p, root_z, root_ky, iy, ix, iz, ij, ia
    snd_p  = np_loc                ! Number of points to send along p (per z,y)
    snd_y  = np_tot*nky_loc        ! Number of points to send along y (per z, full p)
    snd_z  = np_tot*nky_tot*nz_loc ! Number of points to send along z (full y, full p)
    root_p = 0; root_z = 0; root_ky = 0
    a: DO ia= 1,na_tot
      j: DO ij = 1,nj_tot
        x: DO ix = 1,nkx_tot
          z: DO iz = 1,nz_loc
            y: DO iy = 1,nky_loc
              ! fill a buffer to contain a slice of p data at constant j, ky, kx and z
              buffer_pl_cy_zc(1:np_loc) = field_loc(ia,1:np_loc,ij,iy,ix,iz)
              CALL MPI_GATHERV(buffer_pl_cy_zc, snd_p,       mpi_xp_c, &
                              buffer_pt_cy_zc, rcv_p, dsp_p, mpi_xp_c, &
                              root_p, comm_p, ierr)
              buffer_pt_yl_zc(1:np_tot,iy) = buffer_pt_cy_zc(1:np_tot)
            ENDDO y
            ! send the full line on p contained by root_p
            IF(rank_p .EQ. 0) THEN
              CALL MPI_GATHERV(buffer_pt_yl_zc, snd_y,         mpi_xp_c, &
                              buffer_pt_yt_zc, rcv_yp, dsp_yp, mpi_xp_c, &
                              root_ky, comm_ky, ierr)
              buffer_pt_yt_zl(1:np_tot,1:nky_tot,iz) = buffer_pt_yt_zc(1:np_tot,1:nky_tot)
            ENDIF
          ENDDO z
          ! send the full line on y contained by root_kyas
          IF(rank_ky .EQ. 0) THEN
            CALL MPI_GATHERV(buffer_pt_yt_zl, snd_z,           mpi_xp_c, &
                            buffer_pt_yt_zt, rcv_zyp, dsp_zyp, mpi_xp_c, &
                            root_z, comm_z, ierr)
          ENDIF
          ! ID 0 (the one who ouptut) rebuild the whole array
          IF(my_id .EQ. 0) &
            field_tot(ia,1:np_tot,ij,1:nky_tot,ix,1:nz_tot) = buffer_pt_yt_zt(1:np_tot,1:nky_tot,1:nz_tot)
        ENDDO x
      ENDDO j
    ENDDO a
  END SUBROUTINE gather_pjxyz

  !!!!! This is a manual way to do MPI_BCAST !!!!!!!!!!!
  SUBROUTINE manual_3D_bcast(field_,n1,n2,n3)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n1,n2,n3
    COMPLEX(xp), DIMENSION(:,:,:), INTENT(INOUT) :: field_
    COMPLEX(xp) :: buffer(n1,n2,n3)
    INTEGER     :: i_, root, world_rank, world_size, count, i1,i2,i3
    root = 0;
    count = n1*n2*n3;

    CALL MPI_COMM_RANK(comm_p,world_rank,ierr)
    CALL MPI_COMM_SIZE(comm_p,world_size,ierr)
    IF (world_size .GT. 1) THEN
      !! Broadcast phi to the other processes on the same k range (communicator along p)
      IF (world_rank .EQ. root) THEN
        ! Fill the buffer
        DO i3 = 1,n3
          DO i2 = 1,n2
            DO i1 = 1,n1
                buffer(i1,i2,i3) = field_(i1,i2,i3)
              ENDDO
          ENDDO
        ENDDO
        ! Send it to all the other processes
        DO i_ = 0,num_procs_p-1
          IF (i_ .NE. world_rank) &
          CALL MPI_SEND(buffer, count, mpi_xp_c, i_, 0, comm_p, ierr)
        ENDDO
      ELSE
        ! Recieve buffer from root
        CALL MPI_RECV(buffer, count, mpi_xp_c, root, 0, comm_p, MPI_STATUS_IGNORE, ierr)
        ! Write it in phi
        DO i3 = 1,n3
          DO i2 = 1,n2
            DO i1 = 1,n1
              field_(i1,i2,i3) = buffer(i1,i2,i3)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ENDIF
  END SUBROUTINE manual_3D_bcast

  !!!!! This is a manual way to do MPI_BCAST !!!!!!!!!!!
  SUBROUTINE manual_0D_bcast(v)
    IMPLICIT NONE
    COMPLEX(xp), INTENT(INOUT) :: v
    COMPLEX(xp) :: buffer
    INTEGER     :: i_, root, world_rank, world_size, count
    root = 0;
    count = 1;

    CALL MPI_COMM_RANK(comm_z,world_rank,ierr)
    CALL MPI_COMM_SIZE(comm_z,world_size,ierr)
    IF (world_size .GT. 1) THEN
      !! Broadcast phi to the other processes on the same k range (communicator along p)
      IF (world_rank .EQ. root) THEN
        ! Fill the buffer
        buffer = v
        ! Send it to all the other processes
        DO i_ = 0,num_procs_z-1
          IF (i_ .NE. world_rank) &
          CALL MPI_SEND(buffer, count, mpi_xp_c, i_, 0, comm_z, ierr)
        ENDDO
      ELSE
        ! Recieve buffer from root
        CALL MPI_RECV(buffer, count, mpi_xp_c, root, 0, comm_z, MPI_STATUS_IGNORE, ierr)
        ! Write it in phi
        v = buffer
      ENDIF
    ENDIF
  END SUBROUTINE manual_0D_bcast

  ! Routine that exchange ghosts on one dimension
  SUBROUTINE exchange_ghosts_1D(f,nbr_L,nbr_R,np,ng)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: np,ng, nbr_L, nbr_R
    COMPLEX(xp), DIMENSION(:), INTENT(INOUT) :: f
    INTEGER :: ierr, first, last, ig
    first = 1 + ng/2
    last  = np + ng/2
    !!!!!!!!!!! Send ghost to right neighbour !!!!!!!!!!!!!!!!!!!!!!
    DO ig = 1,ng/2
      CALL mpi_sendrecv(f(last-(ig-1)), 1, mpi_xp_c, nbr_R, 14+ig, &
                           f(first-ig), 1, mpi_xp_c, nbr_L, 14+ig, &
                        comm0, MPI_STATUS_IGNORE, ierr)
    ENDDO
    !!!!!!!!!!! Send ghost to left neighbour !!!!!!!!!!!!!!!!!!!!!!
    DO ig = 1,ng/2
    CALL mpi_sendrecv(f(first+(ig-1)), 1, mpi_xp_c, nbr_L, 16+ig, &
                           f(last+ig), 1, mpi_xp_c, nbr_R, 16+ig, &
                      comm0, MPI_STATUS_IGNORE, ierr)
    ENDDO
  END SUBROUTINE exchange_ghosts_1D

END MODULE parallel
