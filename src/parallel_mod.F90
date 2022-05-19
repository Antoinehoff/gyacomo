MODULE parallel
  USE basic
  USE grid, ONLY: Nkx,Nky,Nz, ikys,ikye, izs,ize, local_nky, local_nz
  use prec_const
  IMPLICIT NONE
  ! recieve and displacement counts for gatherv
  INTEGER, DIMENSION(:), ALLOCATABLE :: rcv_z, rcv_y, dsp_z, dsp_y

  PUBLIC :: manual_0D_bcast, manual_3D_bcast, init_parallel_var, &
            gather_xyz!, gather_pjz, gather_pjxyz

CONTAINS

  SUBROUTINE init_parallel_var
    INTEGER :: i_
    !! Y reduction at constant x and z
    ! number of points recieved and displacement for the y reduction
    ALLOCATE(rcv_y(0:num_procs_ky-1),dsp_y(0:num_procs_ky-1)) !Displacement sizes for balance diagnostic
    ! all processes share their local number of points
    CALL MPI_ALLGATHER(local_nky,1,MPI_INTEGER,rcv_y,1,MPI_INTEGER,comm_ky,ierr)
    ! the displacement array can be build from each local_np as
    dsp_y(0)=0
    DO i_=1,num_procs_ky-1
       dsp_y(i_) =dsp_y(i_-1) + rcv_y(i_-1)
    END DO
    print*, rcv_y, dsp_y
    !! Z reduction for full slices of y data but constant x
    ! number of points recieved and displacement for the z reduction
    ALLOCATE(rcv_z(0:num_procs_z-1),dsp_z(0:num_procs_z-1)) !Displacement sizes for balance diagnostic
    ! all processes share their local number of points
    CALL MPI_ALLGATHER(local_nz*Nky,1,MPI_INTEGER,rcv_z,1,MPI_INTEGER,comm_z,ierr)
    ! the displacement array can be build from each local_np as
    dsp_z(0)=0
    DO i_=1,num_procs_z-1
       dsp_z(i_) =dsp_z(i_-1) + rcv_z(i_-1)
    END DO
    print*, rcv_z, dsp_z
  END SUBROUTINE init_parallel_var

  !!!!! Gather a field in spatial coordinates on rank 0 !!!!!
  SUBROUTINE gather_xyz(field_sub,field_full)
    COMPLEX(dp), DIMENSION(ikys:ikye, 1:Nkx, izs:ize), INTENT(IN)    :: field_sub
    COMPLEX(dp), DIMENSION(   1:Nky,  1:Nkx,   1:Nz),  INTENT(INOUT) :: field_full
    COMPLEX(dp), DIMENSION(ikys:ikye)           :: buffer_ly_cz !local y, constant z
    COMPLEX(dp), DIMENSION(   1:Nky )           :: buffer_fy_cz !full  y, constant z
    COMPLEX(dp), DIMENSION(   1:Nky,  izs:ize ) :: buffer_fy_lz !full  y, local z
    COMPLEX(dp), DIMENSION(   1:Nky,    1:Nz  ) :: buffer_fy_fz !full  y, full  z
    INTEGER :: send_count, root_p, root_z, root_ky, ix, iz

    root_p = 0; root_z = 0; root_ky = 0
    IF(rank_p .EQ. root_p) THEN
      DO ix = 1,Nkx
        DO iz = izs,ize
          ! fill a buffer to contain a slice of data at constant kx and z
          print*,iz
          buffer_ly_cz(ikys:ikye) = field_sub(ikys:ikye,ix,iz)
          CALL MPI_GATHERV(buffer_ly_cz,local_nky,MPI_DOUBLE_COMPLEX, &
                           buffer_fy_cz, rcv_y, dsp_y, MPI_DOUBLE_COMPLEX, &
                           root_ky, comm_ky)
          buffer_fy_lz(1:Nky,iz) = buffer_fy_cz(1:Nky)
        ENDDO

        ! send the full line on y contained by root_ky
        IF(rank_ky .EQ. 0) THEN
          CALL MPI_GATHERV(buffer_fy_lz,Nky*local_nz,MPI_DOUBLE_COMPLEX, &
                           buffer_fy_fz, rcv_z, dsp_z, MPI_DOUBLE_COMPLEX, &
                           root_z, comm_z)
        ENDIF
        ! ID 0 (the one who output) rebuild the whole array
        IF(my_id .EQ. 0) &
          field_full(1:Nky,ix,1:Nz) = buffer_fy_fz(1:Nky,1:Nz)
      ENDDO
    ENDIF
  END SUBROUTINE gather_xyz

  !!!!! Gather a field in kinetic + z coordinates on rank 0 !!!!!

  !!!!! Gather a field in kinetic + spatial coordinates on rank 0 !!!!!

  !!!!! This is a manual way to do MPI_BCAST !!!!!!!!!!!
  SUBROUTINE manual_3D_bcast(field_)
    USE grid
    IMPLICIT NONE
    COMPLEX(dp), INTENT(INOUT) :: field_(ikys:ikye,ikxs:ikxe,izs:ize)
    COMPLEX(dp) :: buffer(ikys:ikye,ikxs:ikxe,izs:ize)
    INTEGER     :: i_, root, world_rank, world_size, count
    root = 0;
    count = (ikye-ikys+1) * (ikxe-ikxs+1) * (ize-izs+1);

    CALL MPI_COMM_RANK(comm_p,world_rank,ierr)
    CALL MPI_COMM_SIZE(comm_p,world_size,ierr)
    IF (world_size .GT. 1) THEN
      !! Broadcast phi to the other processes on the same k range (communicator along p)
      IF (world_rank .EQ. root) THEN
        ! Fill the buffer
        DO iz = izs,ize
          DO ikx = ikxs,ikxe
            DO iky = ikys,ikye
                buffer(iky,ikx,iz) = field_(iky,ikx,iz)
              ENDDO
          ENDDO
        ENDDO
        ! Send it to all the other processes
        DO i_ = 0,num_procs_p-1
          IF (i_ .NE. world_rank) &
          CALL MPI_SEND(buffer, count, MPI_DOUBLE_COMPLEX, i_, 0, comm_p, ierr)
        ENDDO
      ELSE
        ! Recieve buffer from root
        CALL MPI_RECV(buffer, count, MPI_DOUBLE_COMPLEX, root, 0, comm_p, MPI_STATUS_IGNORE, ierr)
        ! Write it in phi
        DO iz = izs,ize
          DO ikx = ikxs,ikxe
            DO iky = ikys,ikye
              field_(iky,ikx,iz) = buffer(iky,ikx,iz)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ENDIF
  END SUBROUTINE manual_3D_bcast

  !!!!! This is a manual way to do MPI_BCAST !!!!!!!!!!!
  SUBROUTINE manual_0D_bcast(v)
    USE grid
    IMPLICIT NONE
    COMPLEX(dp), INTENT(INOUT) :: v
    COMPLEX(dp) :: buffer
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
          CALL MPI_SEND(buffer, count, MPI_DOUBLE_COMPLEX, i_, 0, comm_z, ierr)
        ENDDO
      ELSE
        ! Recieve buffer from root
        CALL MPI_RECV(buffer, count, MPI_DOUBLE_COMPLEX, root, 0, comm_z, MPI_STATUS_IGNORE, ierr)
        ! Write it in phi
        v = buffer
      ENDIF
    ENDIF
  END SUBROUTINE manual_0D_bcast


END MODULE parallel
