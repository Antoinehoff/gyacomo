module ghosts
USE mpi
USE prec_const, ONLY: xp, mpi_xp_c
IMPLICIT NONE

INTEGER :: status(MPI_STATUS_SIZE), source, dest, count, ipg

PUBLIC :: update_ghosts_moments, update_ghosts_EM

CONTAINS

SUBROUTINE update_ghosts_moments
  USE grid,     ONLY: total_nz
  USE parallel, ONLY: num_procs_p
  IMPLICIT NONE
  IF (num_procs_p .GT. 1) THEN ! Do it only if we share the p
    CALL update_ghosts_p_mom
  ENDIF
  IF(total_nz .GT. 1) THEN
    CALL update_ghosts_z_mom
  ENDIF
END SUBROUTINE update_ghosts_moments

SUBROUTINE update_ghosts_EM
  USE model,  ONLY :  beta
  USE grid,   ONLY: total_nz
  USE fields, ONLY: phi, psi
  IMPLICIT NONE
  IF(total_nz .GT. 1) THEN
    CALL update_ghosts_z_3D(phi)
    IF(beta .GT. 0._xp) &
      CALL update_ghosts_z_3D(psi)
  ENDIF
END SUBROUTINE update_ghosts_EM


!Communicate p+1, p+2 moments to left neighboor and p-1, p-2 moments to right one
! [a b|C D|e f] : proc n has moments a to f where a,b,e,f are ghosts
!
!proc 0: [0 1 2 3 4|5 6]
!               V V ^ ^
!proc 1:       [3 4|5 6 7 8|9 10]
!                       V V ^  ^
!proc 2:               [7 8|9 10 11 12|13 14]
!                                 V  V  ^  ^
!proc 3:                        [11 12|13 14 15 16|17 18]
!                                                   ^  ^
!Closure by zero truncation :                       0  0

!Communicate p+1, p+2 moments to left neighboor and p-1, p-2 moments to right one
SUBROUTINE update_ghosts_p_mom
  USE time_integration, ONLY: updatetlevel
  USE fields,   ONLY: moments
  USE parallel, ONLY: nbr_R,nbr_L,comm0, exchange_ghosts_1D
  USE grid,     ONLY: local_na,local_np,local_nj,local_nky,local_nkx,local_nz,&
                              ngp,ngj,ngz
  IMPLICIT NONE
  INTEGER :: ierr, first, last, count
  first = 1 + ngp/2
  last  = local_np + ngp/2

  count = (ngp/2)*local_na*(local_nj+ngj)*local_nky*local_nkx*(local_nz+ngz) ! Number of elements to send
  !!!!!! Send to the right, receive from the left
  CALL mpi_sendrecv(moments(:,(last-(ngp/2-1)):(last),:,:,:,:,updatetlevel), count, mpi_xp_c, nbr_R, 14, &
                    moments(:,(first-ngp/2):(first-1),:,:,:,:,updatetlevel), count, mpi_xp_c, nbr_L, 14, &
                    comm0, status, ierr)
  !!!!!!! Send to the left, receive from the right
  CALL mpi_sendrecv(moments(:,(first):(first+(ngp/2-1)),:,:,:,:,updatetlevel), count, mpi_xp_c, nbr_L, 16, &
                    moments(:,(last+1):(last+ngp/2)    ,:,:,:,:,updatetlevel), count, mpi_xp_c, nbr_R, 16, &
                    comm0, status, ierr)
END SUBROUTINE update_ghosts_p_mom

!Communicate z+1, z+2 moments to left neighboor and z-1, z-2 moments to right one
! [a b|C D|e f] : proc n has moments a to f where a,b,e,f are ghosts
!
!proc 0: [0 1 2 3 4|5 6]
!               V V ^ ^
!proc 1:       [3 4|5 6 7 8|9 10]
!                       V V ^  ^
!proc 2:               [7 8|9 10 11 12|13 14]
!                                 V  V  ^  ^
!proc 3:                        [11 12|13 14 15 16|17 18]
!                                                   ^  ^
!Periodic:                                          0  1

SUBROUTINE update_ghosts_z_mom
  USE geometry,         ONLY: ikx_zBC_L, ikx_zBC_R, pb_phase_L, pb_phase_R
  USE time_integration, ONLY: updatetlevel
  USE parallel, ONLY: comm0,nbr_U,nbr_D,num_procs_z
  USE fields,   ONLY: moments
  USE grid,     ONLY: local_na,local_np,local_nj,local_nky,local_nkx,local_nz,&
                      ngp,ngj,ngz
  IMPLICIT NONE
  !! buffer for data exchanges
  COMPLEX(xp),DIMENSION(local_na,local_np+ngp,local_nj+ngj,local_nky,local_nkx,-Ngz/2:Ngz/2) :: buff_pjxy_zBC
  INTEGER :: ikxBC_L, ikxBC_R, ikx, iky, first, last, ig, ierr
  first = 1 + ngz/2
  last  = local_nz + ngz/2
  count = local_na*(local_np+ngp)*(local_nj+ngj)*local_nky*local_nkx
  IF (num_procs_z .GT. 1) THEN
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    !!!!!!!!!!! Send ghost to up neighbour !!!!!!!!!!!!!!!!!!!!!!
    ! Send the last local moment to fill the -1 neighbour ghost
    DO ig=1,ngz/2
      CALL mpi_sendrecv(moments(:,:,:,:,:,last-(ig-1),updatetlevel),count,mpi_xp_c,nbr_U,24+ig, & ! Send to Up the last
                                       buff_pjxy_zBC(:,:,:,:,:,-ig),count,mpi_xp_c,nbr_D,24+ig, & ! Recieve from Down the first-1
                        comm0, status, ierr)
    ENDDO
    !!!!!!!!!!! Send ghost to down neighbour !!!!!!!!!!!!!!!!!!!!!!
    DO ig=1,ngz/2
      CALL mpi_sendrecv(moments(:,:,:,:,:,first+(ig-1),updatetlevel),count,mpi_xp_c,nbr_D,26+ig, & ! Send to Up the last
                                         buff_pjxy_zBC(:,:,:,:,:,ig),count,mpi_xp_c,nbr_U,26+ig, & ! Recieve from Down the first-1
                        comm0, status, ierr)
    ENDDO
  ELSE !No parallel (just copy)
    DO ig=1,ngz/2
      buff_pjxy_zBC(:,:,:,:,:,-ig) = moments(:,:,:,:,:, last-(ig-1),updatetlevel)
      buff_pjxy_zBC(:,:,:,:,:, ig) = moments(:,:,:,:,:,first+(ig-1),updatetlevel)
    ENDDO
  ENDIF
  DO iky = 1,local_nky
    DO ikx = 1,local_nkx
      ikxBC_L = ikx_zBC_L(iky,ikx);
      ! Exchanging the modes that have a periodic pair (from sheared BC)
      IF (ikxBC_L .NE. -99) THEN
        ! Fill the lower ghosts cells with the buffer value (upper cells from LEFT process)
        DO ig=1,ngz/2
          moments(:,:,:,iky,ikx,first-ig,updatetlevel) = pb_phase_L(iky)*buff_pjxy_zBC(:,:,:,iky,ikxBC_L,-ig)
        ENDDO
      ELSE
        DO ig=1,ngz/2
          moments(:,:,:,iky,ikx,first-ig,updatetlevel) = 0._xp
        ENDDO
      ENDIF
      ikxBC_R = ikx_zBC_R(iky,ikx);
      ! Exchanging the modes that have a periodic pair (from sheared BC)
      IF (ikxBC_R .NE. -99) THEN
        ! Fill the upper ghosts cells with the buffer value (lower cells from RIGHT process)
        DO ig=1,ngz/2
          moments(:,:,:,iky,ikx,last+ig,updatetlevel) = pb_phase_R(iky)*buff_pjxy_zBC(:,:,:,iky,ikxBC_R,ig)
        ENDDO
      ELSE
        DO ig=1,ngz/2
          moments(:,:,:,iky,ikx,last+ig,updatetlevel) = 0._xp
        ENDDO
      ENDIF
    ENDDO
  ENDDO
END SUBROUTINE update_ghosts_z_mom

SUBROUTINE update_ghosts_z_3D(field)
  USE geometry, ONLY: ikx_zBC_L, ikx_zBC_R, pb_phase_L, pb_phase_R
  USE parallel, ONLY: nbr_U,nbr_D,comm0,num_procs_z
  USE grid,     ONLY: local_nky,local_nkx,local_nz,ngz
  IMPLICIT NONE
  !! buffer for data exchanges
  COMPLEX(xp),DIMENSION(local_nky,local_nkx,-ngz/2:ngz/2) :: buff_xy_zBC
  COMPLEX(xp), INTENT(INOUT) :: field(local_nky,local_nkx,local_nz+ngz)
  INTEGER :: ikxBC_L, ikxBC_R, ikx, iky, first, last, ig, ierr
  first = 1 + ngz/2
  last  = local_nz + ngz/2
  count = local_nky * local_nkx
  IF (num_procs_z .GT. 1) THEN
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      !!!!!!!!!!! Send ghost to up neighbour !!!!!!!!!!!!!!!!!!!!!!
      DO ig = 1,ngz/2
      CALL mpi_sendrecv(     field(:,:,last-(ig-1)), count, mpi_xp_c, nbr_U, 30+ig, & ! Send to Up the last
                       buff_xy_zBC(:,:,-ig),         count, mpi_xp_c, nbr_D, 30+ig, & ! Receive from Down the first-1
                        comm0, status, ierr)
      ENDDO
      !!!!!!!!!!! Send ghost to down neighbour !!!!!!!!!!!!!!!!!!!!!!
      DO ig = 1,ngz/2
      CALL mpi_sendrecv(     field(:,:,first+(ig-1)), count, mpi_xp_c, nbr_D, 32+ig, & ! Send to Down the first
                       buff_xy_zBC(:,:,ig),           count, mpi_xp_c, nbr_U, 32+ig, & ! Recieve from Up the last+1
                        comm0, status, ierr)
      ENDDO
   ELSE
     ! no parallelization so just copy last cell into first ghosts and vice versa
     DO ig = 1,ngz/2
       buff_xy_zBC(:,:,-ig) = field(:,:,last-(ig-1))
       buff_xy_zBC(:,:, ig) = field(:,:,first+(ig-1))
     ENDDO
   ENDIF
  DO iky = 1,local_nky
    DO ikx = 1,local_nkx
      ikxBC_L = ikx_zBC_L(iky,ikx);
      IF (ikxBC_L .NE. -99) THEN ! Exchanging the modes that have a periodic pair (a)
        DO ig = 1,ngz/2
          field(iky,ikx,first-ig) = pb_phase_L(iky)*buff_xy_zBC(iky,ikxBC_L,-ig)
        ENDDO
      ELSE
        DO ig = 1,ngz/2
          field(iky,ikx,first-ig) = 0._xp
        ENDDO
      ENDIF
      ikxBC_R = ikx_zBC_R(iky,ikx);
      IF (ikxBC_R .NE. -99) THEN ! Exchanging the modes that have a periodic pair (a)
        ! last+1 gets first
        DO ig = 1,ngz/2
          field(iky,ikx,last+ig) = pb_phase_R(iky)*buff_xy_zBC(iky,ikxBC_R,ig)
        ENDDO
      ELSE
        DO ig = 1,ngz/2
          field(iky,ikx,last+ig) = 0._xp
        ENDDO
      ENDIF
    ENDDO
  ENDDO
END SUBROUTINE update_ghosts_z_3D

END MODULE ghosts
