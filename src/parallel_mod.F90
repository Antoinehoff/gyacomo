MODULE parallel
  USE basic
  USE grid, ONLY: Nkx,Nky,Nz, ikys,ikye, izs,ize, local_nky, local_nz, &
                  local_np_e, local_np_i, total_np_e, total_np_i, &
                  pmaxi, pmaxe, ips_e, ipe_e, ips_i, ipe_i, &
                  jmaxi, jmaxe, ijs_e, ije_e, ijs_i, ije_i, &
                  rcv_p_e, rcv_p_i, dsp_p_e, dsp_p_i
  use prec_const
  IMPLICIT NONE
  ! recieve and displacement counts for gatherv
  INTEGER, DIMENSION(:), ALLOCATABLE :: rcv_y, dsp_y, rcv_zy, dsp_zy
  INTEGER, DIMENSION(:), ALLOCATABLE :: rcv_zp_e, dsp_zp_e
  INTEGER, DIMENSION(:), ALLOCATABLE :: rcv_zp_i, dsp_zp_i

  PUBLIC :: manual_0D_bcast, manual_3D_bcast, init_parallel_var, &
            gather_xyz, gather_pjz_e, gather_pjz_i!, gather_pjxyz

CONTAINS

  SUBROUTINE init_parallel_var
    INTEGER :: i_
    !!!!!! XYZ gather variables
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
    !! Z reduction for full slices of y data but constant x
    ! number of points recieved and displacement for the z reduction
    ALLOCATE(rcv_zy(0:num_procs_z-1),dsp_zy(0:num_procs_z-1)) !Displacement sizes for balance diagnostic
    ! all processes share their local number of points
    CALL MPI_ALLGATHER(local_nz*Nky,1,MPI_INTEGER,rcv_zy,1,MPI_INTEGER,comm_z,ierr)
    ! the displacement array can be build from each local_np as
    dsp_zy(0)=0
    DO i_=1,num_procs_z-1
       dsp_zy(i_) =dsp_zy(i_-1) + rcv_zy(i_-1)
    END DO
    print*, rcv_zy, dsp_zy

    !!!!! PJZ gather variables
    ! ELECTRONS
    ! P variables are already defined in grid module (see rcv_p_a, dsp_p_a)
    !! Z reduction for full slices of p data but constant j
    ! number of points recieved and displacement for the z reduction
    ALLOCATE(rcv_zp_e(0:num_procs_z-1),dsp_zp_e(0:num_procs_z-1)) !Displacement sizes for balance diagnostic
    ! all processes share their local number of points
    CALL MPI_ALLGATHER(local_nz*total_np_e,1,MPI_INTEGER,rcv_zp_e,1,MPI_INTEGER,comm_z,ierr)
    ! the displacement array can be build from each local_np as
    dsp_zp_e(0)=0
    DO i_=1,num_procs_z-1
       dsp_zp_e(i_) =dsp_zp_e(i_-1) + dsp_zp_e(i_-1)
    END DO
    ! IONS
    ! P variables are already defined in grid module (see rcv_p_a, dsp_p_a)
    !! Z reduction for full slices of p data but constant j
    ! number of points recieved and displacement for the z reduction
    ALLOCATE(rcv_zp_i(0:num_procs_z-1),dsp_zp_i(0:num_procs_z-1)) !Displacement sizes for balance diagnostic
    ! all processes share their local number of points
    CALL MPI_ALLGATHER(local_nz*total_np_i,1,MPI_INTEGER,rcv_zp_i,1,MPI_INTEGER,comm_z,ierr)
    ! the displacement array can be build from each local_np as
    dsp_zp_i(0)=0
    DO i_=1,num_procs_z-1
       dsp_zp_i(i_) =dsp_zp_i(i_-1) + dsp_zp_i(i_-1)
    END DO

  END SUBROUTINE init_parallel_var

  !!!!! Gather a field in spatial coordinates on rank 0 !!!!!
  SUBROUTINE gather_xyz(field_sub,field_full)
    COMPLEX(dp), DIMENSION(ikys:ikye, 1:Nkx, izs:ize), INTENT(IN)    :: field_sub
    COMPLEX(dp), DIMENSION(   1:Nky,  1:Nkx,   1:Nz),  INTENT(INOUT) :: field_full
    COMPLEX(dp), DIMENSION(ikys:ikye)           :: buffer_ly_cz !local y, constant z
    COMPLEX(dp), DIMENSION(   1:Nky )           :: buffer_fy_cz !full  y, constant z
    COMPLEX(dp), DIMENSION(   1:Nky,  izs:ize ) :: buffer_fy_lz !full  y, local z
    COMPLEX(dp), DIMENSION(   1:Nky,    1:Nz  ) :: buffer_fy_fz !full  y, full  z
    INTEGER :: snd_y, snd_z, root_p, root_z, root_ky, ix, iz

    snd_y  = local_nky    ! Number of points to send along y (per z)
    snd_z  = Nky*local_nz ! Number of points to send along z (full y)

    root_p = 0; root_z = 0; root_ky = 0
    IF(rank_p .EQ. root_p) THEN
      DO ix = 1,Nkx
        DO iz = izs,ize
          ! fill a buffer to contain a slice of data at constant kx and z
          buffer_ly_cz(ikys:ikye) = field_sub(ikys:ikye,ix,iz)
          CALL MPI_GATHERV(buffer_ly_cz, snd_y,        MPI_DOUBLE_COMPLEX, &
                           buffer_fy_cz, rcv_y, dsp_y, MPI_DOUBLE_COMPLEX, &
                           root_ky, comm_ky, ierr)
          buffer_fy_lz(1:Nky,iz) = buffer_fy_cz(1:Nky)
        ENDDO

        ! send the full line on y contained by root_kyas
        IF(rank_ky .EQ. 0) THEN
          CALL MPI_GATHERV(buffer_fy_lz, snd_z,        MPI_DOUBLE_COMPLEX, &
                           buffer_fy_fz, rcv_zy, dsp_zy, MPI_DOUBLE_COMPLEX, &
                           root_z, comm_z, ierr)
        ENDIF
        ! ID 0 (the one who output) rebuild the whole array
        IF(my_id .EQ. 0) &
          field_full(1:Nky,ix,1:Nz) = buffer_fy_fz(1:Nky,1:Nz)
      ENDDO
    ENDIF
  END SUBROUTINE gather_xyz

  !!!!! Gather a field in kinetic + z coordinates on rank 0 !!!!!
  SUBROUTINE gather_pjz_e(field_sub,field_full)
    COMPLEX(dp), DIMENSION(ips_e:ipe_e, ijs_e:ije_e, izs:ize), INTENT(IN)    :: field_sub
    COMPLEX(dp), DIMENSION(   1:pmaxe+1,  1:jmaxe+1,   1:Nz),  INTENT(INOUT) :: field_full
    COMPLEX(dp), DIMENSION(ips_e:ipe_e)         :: buffer_lp_cz !local p, constant z
    COMPLEX(dp), DIMENSION(   1:pmaxe+1 )       :: buffer_fp_cz !full  p, constant z
    COMPLEX(dp), DIMENSION(   1:pmaxe+1,  izs:ize ) :: buffer_fp_lz !full  p, local z
    COMPLEX(dp), DIMENSION(   1:pmaxe+1,    1:Nz  ) :: buffer_fp_fz !full  p, full  z
    INTEGER :: snd_p, snd_z, root_p, root_z, root_ky, ij, iz, Npe

    Npe    = pmaxe+1      ! total number of hermite moments
    snd_p  = local_np_e   ! Number of points to send along y (per z)
    snd_z  = Npe*local_nz ! Number of points to send along z (full y)

    root_p = 0; root_z = 0; root_ky = 0
    IF(rank_ky .EQ. root_ky) THEN
      DO ij = 1,jmaxe+1
        DO iz = izs,ize
          ! fill a buffer to contain a slice of data at constant kx and z
          buffer_lp_cz(ips_e:ipe_e) = field_sub(ips_e:ipe_e,ij,iz)
          CALL MPI_GATHERV(buffer_lp_cz, snd_p,           MPI_DOUBLE_COMPLEX, &
                           buffer_fp_cz, rcv_p_e, dsp_p_e, MPI_DOUBLE_COMPLEX, &
                           root_p, comm_p, ierr)
          buffer_fp_lz(1:Npe,iz) = buffer_fp_cz(1:Npe)
        ENDDO

        ! send the full line on y contained by root_kyas
        IF(rank_p .EQ. 0) THEN
          CALL MPI_GATHERV(buffer_fp_lz, snd_z,              MPI_DOUBLE_COMPLEX, &
                           buffer_fp_fz, rcv_zp_e, dsp_zp_e, MPI_DOUBLE_COMPLEX, &
                           root_z, comm_z, ierr)
        ENDIF
        ! ID 0 (the one who output) rebuild the whole array
        IF(my_id .EQ. 0) &
          field_full(1:Npe,ij,1:Nz) = buffer_fp_fz(1:Npe,1:Nz)
      ENDDO
    ENDIF
  END SUBROUTINE gather_pjz_e

  SUBROUTINE gather_pjz_i(field_sub,field_full)
    COMPLEX(dp), DIMENSION(ips_i:ipe_i, ijs_i:ije_i, izs:ize), INTENT(IN)    :: field_sub
    COMPLEX(dp), DIMENSION(   1:pmaxi+1,  1:jmaxi+1,   1:Nz),  INTENT(INOUT) :: field_full
    COMPLEX(dp), DIMENSION(ips_i:ipe_i)         :: buffer_lp_cz !local p, constant z
    COMPLEX(dp), DIMENSION(   1:pmaxi+1 )       :: buffer_fp_cz !full  p, constant z
    COMPLEX(dp), DIMENSION(   1:pmaxi+1,  izs:ize ) :: buffer_fp_lz !full  p, local z
    COMPLEX(dp), DIMENSION(   1:pmaxi+1,    1:Nz  ) :: buffer_fp_fz !full  p, full  z
    INTEGER :: snd_p, snd_z, root_p, root_z, root_ky, ij, iz, Npi

    Npi    = pmaxi+1      ! total number of hermite moments
    snd_p  = local_np_i   ! Number of points to send along y (per z)
    snd_z  = Npi*local_nz ! Number of points to send along z (full y)

    root_p = 0; root_z = 0; root_ky = 0
    IF(rank_ky .EQ. root_ky) THEN
      DO ij = 1,jmaxi+1
        DO iz = izs,ize
          ! fill a buffer to contain a slice of data at constant kx and z
          buffer_lp_cz(ips_i:ipe_i) = field_sub(ips_i:ipe_i,ij,iz)
          CALL MPI_GATHERV(buffer_lp_cz, snd_p,            MPI_DOUBLE_COMPLEX, &
                           buffer_fp_cz, rcv_p_i, dsp_p_i, MPI_DOUBLE_COMPLEX, &
                           root_p, comm_p, ierr)
          buffer_fp_lz(1:Npi,iz) = buffer_fp_cz(1:Npi)
        ENDDO

        ! send the full line on y contained by root_kyas
        IF(rank_p .EQ. 0) THEN
          CALL MPI_GATHERV(buffer_fp_lz, snd_z,              MPI_DOUBLE_COMPLEX, &
                           buffer_fp_fz, rcv_zp_i, dsp_zp_i, MPI_DOUBLE_COMPLEX, &
                           root_z, comm_z, ierr)
        ENDIF
        ! ID 0 (the one who output) rebuild the whole array
        IF(my_id .EQ. 0) &
          field_full(1:Npi,ij,1:Nz) = buffer_fp_fz(1:Npi,1:Nz)
      ENDDO
    ENDIF
  END SUBROUTINE gather_pjz_i

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
