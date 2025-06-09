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
                              ngp_o2,ngj,ngz
  IMPLICIT NONE
  !! buffer for data exchanges
  ! COMPLEX(xp),DIMENSION(local_na,-ngp_o2:ngp_o2,local_nj+ngj,local_nky,local_nkx,local_nz+ngz) :: buff_jxyz
  COMPLEX(xp),DIMENSION(local_na,local_nj+ngj,local_nky,local_nkx,local_nz+ngz) :: buff_jxyz_in, buff_jxyz_out
  INTEGER :: ierr, first, last, count, ig
  first = 1 + ngp_o2
  last  = local_np + ngp_o2
#if 0
  count = (ngp_o2)*local_na*(local_nj+ngj)*local_nky*local_nkx*(local_nz+ngz) ! Number of elements to send
  !!!!!! Send to the right, receive from the left
  CALL mpi_sendrecv(moments(:,(last-(ngp_o2-1)):(last),:,:,:,:,updatetlevel), count, mpi_xp_c, nbr_R, 14, &
                    moments(:,(first-ngp_o2):(first-1),:,:,:,:,updatetlevel), count, mpi_xp_c, nbr_L, 14, &
                    comm0, status, ierr)
  !!!!!!! Send to the left, receive from the right
  CALL mpi_sendrecv(moments(:,(first):(first+(ngp_o2-1)),:,:,:,:,updatetlevel), count, mpi_xp_c, nbr_L, 16, &
                    moments(:,(last+1):(last+ngp_o2)    ,:,:,:,:,updatetlevel), count, mpi_xp_c, nbr_R, 16, &
                    comm0, status, ierr)
#else
  !! Alternative version (explicit copy here, EXPERIMENTAL)
  !! Send the last local moment to fill the -1 neighbour ghost
  count = local_na*(local_nj+ngj)*local_nky*local_nkx*(local_nz+ngz) ! Number of elements to send
  DO ig=1,ngp_o2
    buff_jxyz_out = moments(:,last-(ig-1),:,:,:,:,updatetlevel)
    CALL mpi_sendrecv(buff_jxyz_out,count,mpi_xp_c,nbr_R,14+ig, & ! Send to Right the last
                      buff_jxyz_in ,count,mpi_xp_c,nbr_L,14+ig, & ! Recieve from Left the first-1
                      comm0, status, ierr)
    moments(:,first-ig,:,:,:,:,updatetlevel) = buff_jxyz_in
  ENDDO
  ! !!!!!!! Send to the left, receive from the right
  DO ig=1,ngp_o2
    buff_jxyz_out = moments(:,first+(ig-1),:,:,:,:,updatetlevel)
    CALL mpi_sendrecv(buff_jxyz_out,count,mpi_xp_c,nbr_L,16+ig, & ! Send to Left the first
                      buff_jxyz_in ,count,mpi_xp_c,nbr_R,16+ig, & ! Recieve from Right the last+1
                      comm0, status, ierr)
    moments(:,last+ig,:,:,:,:,updatetlevel) = buff_jxyz_in
  ENDDO
#endif
END SUBROUTINE update_ghosts_p_mom

! |t =     0.01| Px  = -0.608E-19| Qx  =  0.152E-01|
! |t =     0.02| Px  = -0.336E-19| Qx  =  0.303E-01|
! |t =     0.03| Px  =  0.683E-19| Qx  =  0.454E-01|
! |t =     0.04| Px  = -0.122E-18| Qx  =  0.605E-01|
! |t =     0.05| Px  =  0.322E-19| Qx  =  0.756E-01|
! |t =     0.06| Px  = -0.905E-19| Qx  =  0.906E-01|
! |t =     0.07| Px  =  0.833E-19| Qx  =  0.106    |
! |t =     0.08| Px  =  0.276E-18| Qx  =  0.120    |
! |t =     0.09| Px  = -0.889E-19| Qx  =  0.135    |
! |t =     0.10| Px  =  0.284E-18| Qx  =  0.150    |


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
                      ngp, ngj, ngz_o2
  IMPLICIT NONE
  !! buffer for data exchanges
  COMPLEX(xp),DIMENSION(local_na,local_np+ngp,local_nj+ngj,local_nky,local_nkx,-ngz_o2:ngz_o2) :: buff_pjxy_zBC
  INTEGER :: ikxBC_L, ikxBC_R, ikx, iky, first, last, ig, ierr
  first = 1 + ngz_o2
  last  = local_nz + ngz_o2
  count = local_na*(local_np+ngp)*(local_nj+ngj)*local_nky*local_nkx
  IF (num_procs_z .GT. 1) THEN
    !!!!!!!!!!! Send ghost to up neighbour !!!!!!!!!!!!!!!!!!!!!!
    ! Send the last local moment to fill the -1 neighbour ghost
    DO ig=1,ngz_o2
      CALL mpi_sendrecv(moments(:,:,:,:,:,last-(ig-1),updatetlevel),count,mpi_xp_c,nbr_U,24+ig, & ! Send to Up the last
                                       buff_pjxy_zBC(:,:,:,:,:,-ig),count,mpi_xp_c,nbr_D,24+ig, & ! Recieve from Down the first-1
                        comm0, status, ierr)
    ENDDO
    !!!!!!!!!!! Send ghost to down neighbour !!!!!!!!!!!!!!!!!!!!!!
    DO ig=1,ngz_o2
      CALL mpi_sendrecv(moments(:,:,:,:,:,first+(ig-1),updatetlevel),count,mpi_xp_c,nbr_D,26+ig, & ! Send to Down the first
                                         buff_pjxy_zBC(:,:,:,:,:,ig),count,mpi_xp_c,nbr_U,26+ig, & ! Recieve from Up the last+1
                        comm0, status, ierr)
    ENDDO
  ELSE !No parallel (just copy)
    DO ig=1,ngz_o2
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
        DO ig=1,ngz_o2
          moments(:,:,:,iky,ikx,first-ig,updatetlevel) = pb_phase_L(iky)*buff_pjxy_zBC(:,:,:,iky,ikxBC_L,-ig)
        ENDDO
      ELSE
        DO ig=1,ngz_o2
          moments(:,:,:,iky,ikx,first-ig,updatetlevel) = 0._xp
        ENDDO
      ENDIF
      ikxBC_R = ikx_zBC_R(iky,ikx);
      ! Exchanging the modes that have a periodic pair (from sheared BC)
      IF (ikxBC_R .NE. -99) THEN
        ! Fill the upper ghosts cells with the buffer value (lower cells from RIGHT process)
        DO ig=1,ngz_o2
          moments(:,:,:,iky,ikx,last+ig,updatetlevel) = pb_phase_R(iky)*buff_pjxy_zBC(:,:,:,iky,ikxBC_R,ig)
        ENDDO
      ELSE
        DO ig=1,ngz_o2
          moments(:,:,:,iky,ikx,last+ig,updatetlevel) = 0._xp
        ENDDO
      ENDIF
    ENDDO
  ENDDO
END SUBROUTINE update_ghosts_z_mom

SUBROUTINE update_ghosts_z_3D(field)
  USE geometry, ONLY: ikx_zBC_L, ikx_zBC_R, pb_phase_L, pb_phase_R
  USE parallel, ONLY: nbr_U,nbr_D,comm0,num_procs_z
  USE grid,     ONLY: local_nky,local_nkx,local_nz,ngz,ngz_o2
  IMPLICIT NONE
  !! buffer for data exchanges, the last dimension is indexing the z ghost cells
  ! Example in the full periodic case
  ! (down) |x|x|a|b|...|c|d|x|x| (UP) array along z with old ghost cells
  !                  V
  !              |a|b|c|d|            buffer
  !               1 2 3 4             buffer indices
  !                  V
  !        |c|d|a|b|...|c|d|a|b|      array along z with update ghost cells
  !         3 4             1 2       buffer indices
  COMPLEX(xp),DIMENSION(local_nky,local_nkx,ngz) :: buff_xy_zBC
  COMPLEX(xp), INTENT(INOUT) :: field(local_nky,local_nkx,local_nz+ngz)
  INTEGER :: ikxBC_L, ikxBC_R, ikx, iky, first, last, ig, ierr
  first = 1 + ngz_o2
  last  = local_nz + ngz_o2
  count = local_nky * local_nkx
  IF (num_procs_z .GT. 1) THEN
      !!!!!!!!!!! Send ghost to up neighbour !!!!!!!!!!!!!!!!!!!!!!
      DO ig = 1,ngz_o2
      CALL mpi_sendrecv(     field(:,:,last-(ig-1)), count, mpi_xp_c, nbr_U, 30+ig, & ! Send to Up the last
                       buff_xy_zBC(:,:, ngz-(ig-1)), count, mpi_xp_c, nbr_D, 30+ig, & ! Receive from Down the first-1 (idx 3,4)
                        comm0, status, ierr)
      ENDDO
      !!!!!!!!!!! Send ghost to down neighbour !!!!!!!!!!!!!!!!!!!!!!
      DO ig = 1,ngz_o2
      CALL mpi_sendrecv(     field(:,:,first+(ig-1)), count, mpi_xp_c, nbr_D, 32+ig, & ! Send to Down the first
                       buff_xy_zBC(:,:,ig),           count, mpi_xp_c, nbr_U, 32+ig, & ! Recieve from Up the last+1 (idx 1 2)
                        comm0, status, ierr)
      ENDDO
   ELSE
     ! no parallelization so just copy last cell into first ghosts and vice versa
     DO ig = 1,ngz_o2
       buff_xy_zBC(:,:, ngz-(ig-1)) = field(:,:,last -(ig-1))
       buff_xy_zBC(:,:, ig)         = field(:,:,first+(ig-1))
     ENDDO
   ENDIF
  DO iky = 1,local_nky
    DO ikx = 1,local_nkx
      ikxBC_L = ikx_zBC_L(iky,ikx)
      IF (ikxBC_L .NE. -99) THEN ! Exchanging the modes that have a periodic pair (a)
        DO ig = 1,ngz_o2
          field(iky,ikx,first-ig) = pb_phase_L(iky)*buff_xy_zBC(iky,ikxBC_L,ngz-(ig-1))
        ENDDO
      ELSE
        DO ig = 1,ngz_o2
          field(iky,ikx,first-ig) = 0._xp
        ENDDO
      ENDIF
      ikxBC_R = ikx_zBC_R(iky,ikx)
      IF (ikxBC_R .NE. -99) THEN ! Exchanging the modes that have a periodic pair (a)
        ! last+1 gets first
        DO ig = 1,ngz_o2
          field(iky,ikx,last+ig) = pb_phase_R(iky)*buff_xy_zBC(iky,ikxBC_R,ig)
        ENDDO
      ELSE
        DO ig = 1,ngz_o2
          field(iky,ikx,last+ig) = 0._xp
        ENDDO
      ENDIF
    ENDDO
  ENDDO
END SUBROUTINE update_ghosts_z_3D

END MODULE ghosts
