module ghosts
USE basic
USE grid
USE time_integration
USE model, ONLY : KIN_E, beta
USE geometry, ONLY : SHEARED, ikx_zBC_L, ikx_zBC_R
IMPLICIT NONE

INTEGER :: status(MPI_STATUS_SIZE), source, dest, count, ipg

PUBLIC :: update_ghosts_moments, update_ghosts_EM

CONTAINS

SUBROUTINE update_ghosts_moments
  CALL cpu_time(t0_ghost)

  IF (num_procs_p .GT. 1) THEN ! Do it only if we share the p
      IF(KIN_E)&
      CALL update_ghosts_p_e
      CALL update_ghosts_p_i
  ENDIF

  IF(Nz .GT. 1) THEN
    IF(KIN_E) &
    CALL update_ghosts_z_e
    CALL update_ghosts_z_i
  ENDIF

  tc_ghost = tc_ghost + (t1_ghost - t0_ghost)
END SUBROUTINE update_ghosts_moments

SUBROUTINE update_ghosts_EM
  CALL cpu_time(t0_ghost)

  IF(Nz .GT. 1) THEN
    CALL update_ghosts_z_phi

    IF(beta .GT. 0._dp) &
      CALL update_ghosts_z_psi
  ENDIF

  tc_ghost = tc_ghost + (t1_ghost - t0_ghost)
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
SUBROUTINE update_ghosts_p_e
    USE fields, ONLY : moments_e
    IMPLICIT NONE

    count = (ijge_e-ijgs_e+1)*(ikye-ikys+1)*(ikxe-ikxs+1)*(izge-izgs+1)

    !!!!!!!!!!! Send ghost to right neighbour !!!!!!!!!!!!!!!!!!!!!!
    ! Send the last local moment to fill the -1 neighbour ghost
    CALL mpi_sendrecv(moments_e(ipe_e  ,:,:,:,:,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_R, 10, & ! Send to right
                      moments_e(ips_e-1,:,:,:,:,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_L, 10, & ! Recieve from left
                      comm0, status, ierr)
    IF (deltape .EQ. 1) & ! If we have odd Hermite degrees we need a 2nd order stencil
    CALL mpi_sendrecv(moments_e(ipe_e-1,:,:,:,:,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_R, 11, & ! Send to right
                      moments_e(ips_e-2,:,:,:,:,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_L, 11, & ! Recieve from left
                      comm0, status, ierr)

    !!!!!!!!!!! Send ghost to left neighbour !!!!!!!!!!!!!!!!!!!!!!
    CALL mpi_sendrecv(moments_e(ips_e  ,:,:,:,:,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_L, 12, & ! Send to left
                      moments_e(ipe_e+1,:,:,:,:,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_R, 12, & ! Recieve from right
                      comm0, status, ierr)
    IF (deltape .EQ. 1) & ! If we have odd Hermite degrees we need a 2nd order stencil
    CALL mpi_sendrecv(moments_e(ips_e+1,:,:,:,:,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_L, 13, & ! Send to left
                      moments_e(ipe_e+2,:,:,:,:,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_R, 13, & ! Recieve from right
                      comm0, status, ierr)

END SUBROUTINE update_ghosts_p_e

!Communicate p+1, p+2 moments to left neighboor and p-1, p-2 moments to right one
SUBROUTINE update_ghosts_p_i
  USE fields, ONLY : moments_i
    IMPLICIT NONE

    count = (ijge_i-ijgs_i+1)*(ikye-ikys+1)*(ikxe-ikxs+1)*(izge-izgs+1) ! Number of elements sent

    !!!!!!!!!!! Send ghost to right neighbour !!!!!!!!!!!!!!!!!!!!!!
    CALL mpi_sendrecv(moments_i(ipe_i  ,:,:,:,:,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_R, 14, &
                      moments_i(ips_i-1,:,:,:,:,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_L, 14, &
                      comm0, status, ierr)
    IF (deltapi .EQ. 1) & ! If we have odd Hermite degrees we need a 2nd order stencil
    CALL mpi_sendrecv(moments_i(ipe_i-1,:,:,:,:,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_R, 15, &
                      moments_i(ips_i-2,:,:,:,:,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_L, 15, &
                      comm0, status, ierr)

    !!!!!!!!!!! Send ghost to left neighbour !!!!!!!!!!!!!!!!!!!!!!
    CALL mpi_sendrecv(moments_i(ips_i  ,:,:,:,:,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_L, 16, &
                      moments_i(ipe_i+1,:,:,:,:,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_R, 16, &
                      comm0, status, ierr)
    IF (deltapi .EQ. 1) & ! If we have odd Hermite degrees we need a 2nd order stencil
    CALL mpi_sendrecv(moments_i(ips_i+1,:,:,:,:,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_L, 17, &
                      moments_i(ipe_i+2,:,:,:,:,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_R, 17, &
                      comm0, status, ierr)

END SUBROUTINE update_ghosts_p_i

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
SUBROUTINE update_ghosts_z_e
  USE parallel, ONLY : buff_pjxy_zBC_e
  USE fields, ONLY : moments_e
  IMPLICIT NONE
  INTEGER :: ikxBC_L, ikxBC_R
  IF(Nz .GT. 1) THEN
    IF (num_procs_z .GT. 1) THEN
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      count = (ipge_e-ipgs_e+1)*(ijge_e-ijgs_e+1)*(ikye-ikys+1)*(ikxe-ikxs+1)

      !!!!!!!!!!! Send ghost to up neighbour !!!!!!!!!!!!!!!!!!!!!!
      ! Send the last local moment to fill the -1 neighbour ghost
      CALL mpi_sendrecv(moments_e(:,:,:,:,ize  ,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_U, 20, & ! Send to Up the last
                                  buff_pjxy_zBC_e(:,:,:,:,-1), count, MPI_DOUBLE_COMPLEX, nbr_D, 20, & ! Recieve from Down the first-1
                        comm0, status, ierr)
      CALL mpi_sendrecv(moments_e(:,:,:,:,ize-1,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_U, 21, & ! Send to Up the last
                                  buff_pjxy_zBC_e(:,:,:,:,-2), count, MPI_DOUBLE_COMPLEX, nbr_D, 21, & ! Recieve from Down the first-1
                        comm0, status, ierr)
      !!!!!!!!!!! Send ghost to down neighbour !!!!!!!!!!!!!!!!!!!!!!
      CALL mpi_sendrecv(moments_e(:,:,:,:,izs  ,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_D, 22, & ! Send to Up the last
                                  buff_pjxy_zBC_e(:,:,:,:,+1), count, MPI_DOUBLE_COMPLEX, nbr_U, 22, & ! Recieve from Down the first-1
                        comm0, status, ierr)
      CALL mpi_sendrecv(moments_e(:,:,:,:,izs+1,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_D, 23, & ! Send to Up the last
                                  buff_pjxy_zBC_e(:,:,:,:,+2), count, MPI_DOUBLE_COMPLEX, nbr_U, 23, & ! Recieve from Down the first-1
                        comm0, status, ierr)
    ELSE !No parallel (copy)
      buff_pjxy_zBC_e(:,:,:,:,-1) = moments_e(:,:,:,:,ize  ,updatetlevel)
      buff_pjxy_zBC_e(:,:,:,:,-2) = moments_e(:,:,:,:,ize-1,updatetlevel)
      buff_pjxy_zBC_e(:,:,:,:,+1) = moments_e(:,:,:,:,izs  ,updatetlevel)
      buff_pjxy_zBC_e(:,:,:,:,+2) = moments_e(:,:,:,:,izs+1,updatetlevel)
    ENDIF
    DO iky = ikys,ikye
      DO ikx = ikxs,ikxe
        ikxBC_L = ikx_zBC_L(iky,ikx);
        IF (ikxBC_L .NE. -99) THEN ! Exchanging the modes that have a periodic pair (a)
          ! first-1 gets last
          moments_e(:,:,iky,ikx,izs-1,updatetlevel) = buff_pjxy_zBC_e(:,:,iky,ikxBC_L,-1)
          ! first-2 gets last-1
          moments_e(:,:,iky,ikx,izs-2,updatetlevel) = buff_pjxy_zBC_e(:,:,iky,ikxBC_L,-2)
        ELSE
          moments_e(:,:,iky,ikx,izs-1,updatetlevel) = 0._dp
          moments_e(:,:,iky,ikx,izs-2,updatetlevel) = 0._dp
        ENDIF
        ikxBC_R = ikx_zBC_R(iky,ikx);
        IF (ikxBC_R .NE. -99) THEN ! Exchanging the modes that have a periodic pair (a)
          ! last+1 gets first
          moments_e(:,:,iky,ikx,ize+1,updatetlevel) = buff_pjxy_zBC_e(:,:,iky,ikxBC_R,+1)
          ! last+2 gets first+1
          moments_e(:,:,iky,ikx,ize+2,updatetlevel) = buff_pjxy_zBC_e(:,:,iky,ikxBC_R,+2)
        ELSE
          moments_e(:,:,iky,ikx,ize+1,updatetlevel) = 0._dp
          moments_e(:,:,iky,ikx,ize+2,updatetlevel) = 0._dp
        ENDIF
      ENDDO
    ENDDO
  ENDIF
END SUBROUTINE update_ghosts_z_e

SUBROUTINE update_ghosts_z_i
  USE parallel, ONLY : buff_pjxy_zBC_i
  USE fields, ONLY : moments_i
  IMPLICIT NONE
  INTEGER :: ikxBC_L, ikxBC_R
  IF(Nz .GT. 1) THEN
    IF (num_procs_z .GT. 1) THEN
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
      count = (ipge_i-ipgs_i+1)*(ijge_i-ijgs_i+1)*(ikye-ikys+1)*(ikxe-ikxs+1)

      !!!!!!!!!!! Send ghost to up neighbour !!!!!!!!!!!!!!!!!!!!!!
      ! Send the last local moment to fill the -1 neighbour ghost
      CALL mpi_sendrecv(moments_i(:,:,:,:,ize  ,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_U, 24, & ! Send to Up the last
                                  buff_pjxy_zBC_i(:,:,:,:,-1), count, MPI_DOUBLE_COMPLEX, nbr_D, 24, & ! Recieve from Down the first-1
                        comm0, status, ierr)
      CALL mpi_sendrecv(moments_i(:,:,:,:,ize-1,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_U, 25, & ! Send to Up the last
                                  buff_pjxy_zBC_i(:,:,:,:,-2), count, MPI_DOUBLE_COMPLEX, nbr_D, 25, & ! Recieve from Down the first-1
                        comm0, status, ierr)
      !!!!!!!!!!! Send ghost to down neighbour !!!!!!!!!!!!!!!!!!!!!!
      CALL mpi_sendrecv(moments_i(:,:,:,:,izs  ,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_D, 26, & ! Send to Up the last
                                  buff_pjxy_zBC_i(:,:,:,:,+1), count, MPI_DOUBLE_COMPLEX, nbr_U, 26, & ! Recieve from Down the first-1
                        comm0, status, ierr)
      CALL mpi_sendrecv(moments_i(:,:,:,:,izs+1,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_D, 27, & ! Send to Up the last
                                  buff_pjxy_zBC_i(:,:,:,:,+2), count, MPI_DOUBLE_COMPLEX, nbr_U, 27, & ! Recieve from Down the first-1
                        comm0, status, ierr)
    ELSE !No parallel (copy)
      buff_pjxy_zBC_i(:,:,:,:,-1) = moments_i(:,:,:,:,ize  ,updatetlevel)
      buff_pjxy_zBC_i(:,:,:,:,-2) = moments_i(:,:,:,:,ize-1,updatetlevel)
      buff_pjxy_zBC_i(:,:,:,:,+1) = moments_i(:,:,:,:,izs  ,updatetlevel)
      buff_pjxy_zBC_i(:,:,:,:,+2) = moments_i(:,:,:,:,izs+1,updatetlevel)
    ENDIF
    DO iky = ikys,ikye
      DO ikx = ikxs,ikxe
        ikxBC_L = ikx_zBC_L(iky,ikx);
        IF (ikxBC_L .NE. -99) THEN ! Exchanging the modes that have a periodic pair (a)
          ! first-1 gets last
          moments_i(:,:,iky,ikx,izs-1,updatetlevel) = buff_pjxy_zBC_i(:,:,iky,ikxBC_L,-1)
          ! first-2 gets last-1
          moments_i(:,:,iky,ikx,izs-2,updatetlevel) = buff_pjxy_zBC_i(:,:,iky,ikxBC_L,-2)
        ELSE
          moments_i(:,:,iky,ikx,izs-1,updatetlevel) = 0._dp
          moments_i(:,:,iky,ikx,izs-2,updatetlevel) = 0._dp
        ENDIF
        ikxBC_R = ikx_zBC_R(iky,ikx);
        IF (ikxBC_R .NE. -99) THEN ! Exchanging the modes that have a periodic pair (a)
          ! last+1 gets first
          moments_i(:,:,iky,ikx,ize+1,updatetlevel) = buff_pjxy_zBC_i(:,:,iky,ikxBC_R,+1)
          ! last+2 gets first+1
          moments_i(:,:,iky,ikx,ize+2,updatetlevel) = buff_pjxy_zBC_i(:,:,iky,ikxBC_R,+2)
        ELSE
          moments_i(:,:,iky,ikx,ize+1,updatetlevel) = 0._dp
          moments_i(:,:,iky,ikx,ize+2,updatetlevel) = 0._dp
        ENDIF
      ENDDO
    ENDDO
  ENDIF
END SUBROUTINE update_ghosts_z_i

SUBROUTINE update_ghosts_z_phi
  USE parallel, ONLY : buff_xy_zBC
  USE fields,   ONLY : phi
  IMPLICIT NONE
  INTEGER :: ikxBC_L, ikxBC_R
  CALL cpu_time(t1_ghost)
  IF(Nz .GT. 1) THEN
    IF (num_procs_z .GT. 1) THEN
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        count = (ikye-ikys+1) * (ikxe-ikxs+1)
        !!!!!!!!!!! Send ghost to up neighbour !!!!!!!!!!!!!!!!!!!!!!
        CALL mpi_sendrecv(     phi(:,:,ize  ), count, MPI_DOUBLE_COMPLEX, nbr_U, 30, & ! Send to Up the last
                          buff_xy_zBC(:,:,-1), count, MPI_DOUBLE_COMPLEX, nbr_D, 30, & ! Receive from Down the first-1
                          comm0, status, ierr)

        CALL mpi_sendrecv(     phi(:,:,ize-1), count, MPI_DOUBLE_COMPLEX, nbr_U, 31, & ! Send to Up the last
                          buff_xy_zBC(:,:,-2), count, MPI_DOUBLE_COMPLEX, nbr_D, 31, & ! Receive from Down the first-2
                          comm0, status, ierr)

        !!!!!!!!!!! Send ghost to down neighbour !!!!!!!!!!!!!!!!!!!!!!
        CALL mpi_sendrecv(     phi(:,:,izs  ), count, MPI_DOUBLE_COMPLEX, nbr_D, 32, & ! Send to Down the first
                          buff_xy_zBC(:,:,+1), count, MPI_DOUBLE_COMPLEX, nbr_U, 32, & ! Recieve from Up the last+1
                          comm0, status, ierr)

        CALL mpi_sendrecv(     phi(:,:,izs+1), count, MPI_DOUBLE_COMPLEX, nbr_D, 33, & ! Send to Down the first
                          buff_xy_zBC(:,:,+2), count, MPI_DOUBLE_COMPLEX, nbr_U, 33, & ! Recieve from Up the last+2
                          comm0, status, ierr)
     ELSE
       buff_xy_zBC(:,:,-1) = phi(:,:,ize  )
       buff_xy_zBC(:,:,-2) = phi(:,:,ize-1)
       buff_xy_zBC(:,:,+1) = phi(:,:,izs  )
       buff_xy_zBC(:,:,+2) = phi(:,:,izs+1)
     ENDIF
    DO iky = ikys,ikye
      DO ikx = ikxs,ikxe
        ikxBC_L = ikx_zBC_L(iky,ikx);
        IF (ikxBC_L .NE. -99) THEN ! Exchanging the modes that have a periodic pair (a)
          ! first-1 gets last
          phi(iky,ikx,izs-1) = buff_xy_zBC(iky,ikxBC_L,-1)
          ! first-2 gets last-1
          phi(iky,ikx,izs-2) = buff_xy_zBC(iky,ikxBC_L,-2)
        ELSE
          phi(iky,ikx,izs-1) = 0._dp
          phi(iky,ikx,izs-2) = 0._dp
        ENDIF
        ikxBC_R = ikx_zBC_R(iky,ikx);
        IF (ikxBC_R .NE. -99) THEN ! Exchanging the modes that have a periodic pair (a)
          ! last+1 gets first
          phi(iky,ikx,ize+1) = buff_xy_zBC(iky,ikxBC_R,+1)
          ! last+2 gets first+1
          phi(iky,ikx,ize+2) = buff_xy_zBC(iky,ikxBC_R,+2)
        ELSE
          phi(iky,ikx,ize+1) = 0._dp
          phi(iky,ikx,ize+2) = 0._dp
        ENDIF
      ENDDO
    ENDDO
  ENDIF
  CALL cpu_time(t1_ghost)
  tc_ghost = tc_ghost + (t1_ghost - t0_ghost)
END SUBROUTINE update_ghosts_z_phi

SUBROUTINE update_ghosts_z_psi
  USE parallel, ONLY : buff_xy_zBC
  USE fields, ONLY : psi
  IMPLICIT NONE
  INTEGER :: ikxBC_L, ikxBC_R
  CALL cpu_time(t1_ghost)
  IF(Nz .GT. 1) THEN
    IF (num_procs_z .GT. 1) THEN
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        count = (ikye-ikys+1) * (ikxe-ikxs+1)
        !!!!!!!!!!! Send ghost to up neighbour !!!!!!!!!!!!!!!!!!!!!!
        CALL mpi_sendrecv(     psi(:,:,ize  ), count, MPI_DOUBLE_COMPLEX, nbr_U, 30, & ! Send to Up the last
                          buff_xy_zBC(:,:,-1), count, MPI_DOUBLE_COMPLEX, nbr_D, 30, & ! Receive from Down the first-1
                          comm0, status, ierr)

        CALL mpi_sendrecv(     psi(:,:,ize-1), count, MPI_DOUBLE_COMPLEX, nbr_U, 31, & ! Send to Up the last
                          buff_xy_zBC(:,:,-2), count, MPI_DOUBLE_COMPLEX, nbr_D, 31, & ! Receive from Down the first-2
                          comm0, status, ierr)

        !!!!!!!!!!! Send ghost to down neighbour !!!!!!!!!!!!!!!!!!!!!!
        CALL mpi_sendrecv(     psi(:,:,izs  ), count, MPI_DOUBLE_COMPLEX, nbr_D, 32, & ! Send to Down the first
                          buff_xy_zBC(:,:,+1), count, MPI_DOUBLE_COMPLEX, nbr_U, 32, & ! Recieve from Up the last+1
                          comm0, status, ierr)

        CALL mpi_sendrecv(     psi(:,:,izs+1), count, MPI_DOUBLE_COMPLEX, nbr_D, 33, & ! Send to Down the first
                          buff_xy_zBC(:,:,+2), count, MPI_DOUBLE_COMPLEX, nbr_U, 33, & ! Recieve from Up the last+2
                          comm0, status, ierr)
     ELSE
       buff_xy_zBC(:,:,-1) = psi(:,:,ize  )
       buff_xy_zBC(:,:,-2) = psi(:,:,ize-1)
       buff_xy_zBC(:,:,+1) = psi(:,:,izs  )
       buff_xy_zBC(:,:,+2) = psi(:,:,izs+1)
     ENDIF
    DO iky = ikys,ikye
      DO ikx = ikxs,ikxe
        ikxBC_L = ikx_zBC_L(iky,ikx);
        IF (ikxBC_L .NE. -99) THEN ! Exchanging the modes that have a periodic pair (a)
          ! first-1 gets last
          psi(iky,ikx,izs-1) = buff_xy_zBC(iky,ikxBC_L,-1)
          ! first-2 gets last-1
          psi(iky,ikx,izs-2) = buff_xy_zBC(iky,ikxBC_L,-2)
        ELSE
          psi(iky,ikx,izs-1) = 0._dp
          psi(iky,ikx,izs-2) = 0._dp
        ENDIF
        ikxBC_R = ikx_zBC_R(iky,ikx);
        IF (ikxBC_R .NE. -99) THEN ! Exchanging the modes that have a periodic pair (a)
          ! last+1 gets first
          psi(iky,ikx,ize+1) = buff_xy_zBC(iky,ikxBC_R,+1)
          ! last+2 gets first+1
          psi(iky,ikx,ize+2) = buff_xy_zBC(iky,ikxBC_R,+2)
        ELSE
          psi(iky,ikx,ize+1) = 0._dp
          psi(iky,ikx,ize+2) = 0._dp
        ENDIF
      ENDDO
    ENDDO
  ENDIF
  CALL cpu_time(t1_ghost)
  tc_ghost = tc_ghost + (t1_ghost - t0_ghost)
END SUBROUTINE update_ghosts_z_psi


END MODULE ghosts
