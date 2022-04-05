module ghosts
USE basic
USE fields, ONLY : moments_e, moments_i, phi
USE grid
USE time_integration
USE model, ONLY : KIN_E

IMPLICIT NONE

INTEGER :: status(MPI_STATUS_SIZE), source, dest, count, ipg

PUBLIC :: update_ghosts_p_moments, update_ghosts_z_phi, update_ghosts_z_moments

CONTAINS

SUBROUTINE update_ghosts_p_moments
    CALL cpu_time(t0_ghost)

    IF (num_procs_p .GT. 1) THEN ! Do it only if we share the p
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        IF(KIN_E) CALL update_ghosts_p_e

        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        CALL update_ghosts_p_i
    ENDIF

    CALL update_ghosts_z_moments
    CALL cpu_time(t1_ghost)
    tc_ghost = tc_ghost + (t1_ghost - t0_ghost)
END SUBROUTINE update_ghosts_p_moments


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

    IMPLICIT NONE
    count = (ijge_e-ijgs_e+1)*(ikxe-ikxs+1)*(ikye-ikys+1)*(izge-izgs+1)

    !!!!!!!!!!! Send ghost to right neighbour !!!!!!!!!!!!!!!!!!!!!!
    ! Send the last local moment to fill the -1 neighbour ghost
    CALL mpi_sendrecv(moments_e(ipe_e  ,ijgs_e:ijge_e,ikxs:ikxe,ikys:ikye,izgs:izge,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_R, 10, & ! Send to right
                      moments_e(ips_e-1,ijgs_e:ijge_e,ikxs:ikxe,ikys:ikye,izgs:izge,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_L, 10, & ! Recieve from left
                      comm0, status, ierr)
    IF (deltape .EQ. 1) & ! If we have odd Hermite degrees we need a 2nd order stencil
    CALL mpi_sendrecv(moments_e(ipe_e-1,ijgs_e:ijge_e,ikxs:ikxe,ikys:ikye,izgs:izge,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_R, 11, & ! Send to right
                      moments_e(ips_e-2,ijgs_e:ijge_e,ikxs:ikxe,ikys:ikye,izgs:izge,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_L, 11, & ! Recieve from left
                      comm0, status, ierr)

    !!!!!!!!!!! Send ghost to left neighbour !!!!!!!!!!!!!!!!!!!!!!
    CALL mpi_sendrecv(moments_e(ips_e  ,ijgs_e:ijge_e,ikxs:ikxe,ikys:ikye,izgs:izge,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_L, 12, & ! Send to left
                      moments_e(ipe_e+1,ijgs_e:ijge_e,ikxs:ikxe,ikys:ikye,izgs:izge,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_R, 12, & ! Recieve from right
                      comm0, status, ierr)
    IF (deltape .EQ. 1) & ! If we have odd Hermite degrees we need a 2nd order stencil
    CALL mpi_sendrecv(moments_e(ips_e+1,ijgs_e:ijge_e,ikxs:ikxe,ikys:ikye,izgs:izge,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_L, 13, & ! Send to left
                      moments_e(ipe_e+2,ijgs_e:ijge_e,ikxs:ikxe,ikys:ikye,izgs:izge,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_R, 13, & ! Recieve from right
                      comm0, status, ierr)

END SUBROUTINE update_ghosts_p_e

!Communicate p+1, p+2 moments to left neighboor and p-1, p-2 moments to right one
SUBROUTINE update_ghosts_p_i

    IMPLICIT NONE

    count = (ijge_i-ijgs_i+1)*(ikxe-ikxs+1)*(ikye-ikys+1)*(izge-izgs+1) ! Number of elements sent

    !!!!!!!!!!! Send ghost to right neighbour !!!!!!!!!!!!!!!!!!!!!!
    CALL mpi_sendrecv(moments_i(ipe_i  ,ijgs_i:ijge_i,ikxs:ikxe,ikys:ikye,izgs:izge,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_R, 14, &
                      moments_i(ips_i-1,ijgs_i:ijge_i,ikxs:ikxe,ikys:ikye,izgs:izge,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_L, 14, &
                      comm0, status, ierr)
    IF (deltapi .EQ. 1) & ! If we have odd Hermite degrees we need a 2nd order stencil
    CALL mpi_sendrecv(moments_i(ipe_i-1,ijgs_i:ijge_i,ikxs:ikxe,ikys:ikye,izgs:izge,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_R, 15, &
                      moments_i(ips_i-2,ijgs_i:ijge_i,ikxs:ikxe,ikys:ikye,izgs:izge,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_L, 15, &
                      comm0, status, ierr)

    !!!!!!!!!!! Send ghost to left neighbour !!!!!!!!!!!!!!!!!!!!!!
    CALL mpi_cart_shift(comm0, 0, -1, source , dest , ierr)
    CALL mpi_sendrecv(moments_i(ips_i  ,ijgs_i:ijge_i,ikxs:ikxe,ikys:ikye,izgs:izge,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_L, 16, &
                      moments_i(ipe_i+1,ijgs_i:ijge_i,ikxs:ikxe,ikys:ikye,izgs:izge,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_R, 16, &
                      comm0, status, ierr)
    IF (deltapi .EQ. 1) & ! If we have odd Hermite degrees we need a 2nd order stencil
    CALL mpi_sendrecv(moments_i(ips_i+1,ijgs_i:ijge_i,ikxs:ikxe,ikys:ikye,izgs:izge,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_L, 17, &
                      moments_i(ipe_i+2,ijgs_i:ijge_i,ikxs:ikxe,ikys:ikye,izgs:izge,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_R, 17, &
                      comm0, status, ierr)

END SUBROUTINE update_ghosts_p_i

SUBROUTINE update_ghosts_z_moments
  IMPLICIT NONE

    IF(KIN_E) &
    CALL update_ghosts_z_e
    CALL update_ghosts_z_i

END SUBROUTINE update_ghosts_z_moments

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

  IMPLICIT NONE

  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  IF (num_procs_z .GT. 1) THEN
    count = (ipge_e-ipgs_e+1)*(ijge_e-ijgs_e+1)*(ikxe-ikxs+1)*(ikye-ikys+1)

    !!!!!!!!!!! Send ghost to right neighbour !!!!!!!!!!!!!!!!!!!!!!
    ! Send the last local moment to fill the -1 neighbour ghost
    CALL mpi_sendrecv(moments_e(ipgs_e:ipge_e,ijgs_e:ijge_e,ikxs:ikxe,ikys:ikye,ize  ,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_U, 20, & ! Send to right
                      moments_e(ipgs_e:ipge_e,ijgs_e:ijge_e,ikxs:ikxe,ikys:ikye,ize-1,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_D, 20, & ! Recieve from left
                      comm0, status, ierr)
    CALL mpi_sendrecv(moments_e(ipgs_e:ipge_e,ijgs_e:ijge_e,ikxs:ikxe,ikys:ikye,ize-1,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_U, 21, & ! Send to right
                      moments_e(ipgs_e:ipge_e,ijgs_e:ijge_e,ikxs:ikxe,ikys:ikye,ize-2,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_D, 21, & ! Recieve from left
                      comm0, status, ierr)

    !!!!!!!!!!! Send ghost to left neighbour !!!!!!!!!!!!!!!!!!!!!!
    CALL mpi_sendrecv(moments_e(ipgs_e:ipge_e,ijgs_e:ijge_e,ikxs:ikxe,ikys:ikye,izs  ,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_D, 22, & ! Send to left
                      moments_e(ipgs_e:ipge_e,ijgs_e:ijge_e,ikxs:ikxe,ikys:ikye,izs+1,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_U, 22, & ! Recieve from right
                      comm0, status, ierr)
    CALL mpi_sendrecv(moments_e(ipgs_e:ipge_e,ijgs_e:ijge_e,ikxs:ikxe,ikys:ikye,izs+1,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_D, 23, & ! Send to left
                      moments_e(ipgs_e:ipge_e,ijgs_e:ijge_e,ikxs:ikxe,ikys:ikye,izs+2,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_U, 23, & ! Recieve from right
                      comm0, status, ierr)
  ELSE ! still need to perform periodic boundary conditions
      moments_e(:,:,:,:,   0,updatetlevel) = moments_e(:,:,:,:,  Nz,updatetlevel)
      moments_e(:,:,:,:,  -1,updatetlevel) = moments_e(:,:,:,:,Nz-1,updatetlevel)
      moments_e(:,:,:,:,Nz+1,updatetlevel) = moments_e(:,:,:,:,   1,updatetlevel)
      moments_e(:,:,:,:,Nz+2,updatetlevel) = moments_e(:,:,:,:,   2,updatetlevel)
  ENDIF

END SUBROUTINE update_ghosts_z_e

SUBROUTINE update_ghosts_z_i

  IMPLICIT NONE
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  IF (num_procs_z .GT. 1) THEN
    count = (ipge_i-ipgs_i+1)*(ijge_i-ijgs_i+1)*(ikxe-ikxs+1)*(ikye-ikys+1)

    !!!!!!!!!!! Send ghost to right neighbour !!!!!!!!!!!!!!!!!!!!!!
    ! Send the last local moment to fill the -1 neighbour ghost
    CALL mpi_sendrecv(moments_i(ipgs_i:ipge_i,ijgs_i:ijge_i,ikxs:ikxe,ikys:ikye,ize  ,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_U, 30, & ! Send to right
                      moments_i(ipgs_i:ipge_i,ijgs_i:ijge_i,ikxs:ikxe,ikys:ikye,ize-1,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_D, 30, & ! Recieve from left
                      comm0, status, ierr)
    CALL mpi_sendrecv(moments_i(ipgs_i:ipge_i,ijgs_i:ijge_i,ikxs:ikxe,ikys:ikye,ize-1,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_U, 31, & ! Send to right
                      moments_i(ipgs_i:ipge_i,ijgs_i:ijge_i,ikxs:ikxe,ikys:ikye,ize-2,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_D, 31, & ! Recieve from left
                      comm0, status, ierr)

    !!!!!!!!!!! Send ghost to left neighbour !!!!!!!!!!!!!!!!!!!!!!
    CALL mpi_sendrecv(moments_i(ipgs_i:ipge_i,ijgs_i:ijge_i,ikxs:ikxe,ikys:ikye,izs  ,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_D, 32, & ! Send to left
                      moments_i(ipgs_i:ipge_i,ijgs_i:ijge_i,ikxs:ikxe,ikys:ikye,izs+1,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_U, 32, & ! Recieve from right
                      comm0, status, ierr)
    CALL mpi_sendrecv(moments_i(ipgs_i:ipge_i,ijgs_i:ijge_i,ikxs:ikxe,ikys:ikye,izs+1,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_D, 33, & ! Send to left
                      moments_i(ipgs_i:ipge_i,ijgs_i:ijge_i,ikxs:ikxe,ikys:ikye,izs+2,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_U, 33, & ! Recieve from right
                      comm0, status, ierr)
  ELSE ! still need to perform periodic boundary conditions
    moments_i(:,:,:,:,   0,updatetlevel) = moments_i(:,:,:,:,  Nz,updatetlevel)
    moments_i(:,:,:,:,  -1,updatetlevel) = moments_i(:,:,:,:,Nz-1,updatetlevel)
    moments_i(:,:,:,:,Nz+1,updatetlevel) = moments_i(:,:,:,:,   1,updatetlevel)
    moments_i(:,:,:,:,Nz+2,updatetlevel) = moments_i(:,:,:,:,   2,updatetlevel)
  ENDIF
END SUBROUTINE update_ghosts_z_i

SUBROUTINE update_ghosts_z_phi

  IMPLICIT NONE
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  IF (num_procs_z .GT. 1) THEN
    count = (ikxe-ikxs+1)*(ikye-ikys+1)

    !!!!!!!!!!! Send ghost to right neighbour !!!!!!!!!!!!!!!!!!!!!!
    ! Send the last local moment to fill the -1 neighbour ghost
    CALL mpi_sendrecv(phi(ikxs:ikxe,ikys:ikye,  ize), count, MPI_DOUBLE_COMPLEX, nbr_U, 40, & ! Send to right
                      phi(ikxs:ikxe,ikys:ikye,ize-1), count, MPI_DOUBLE_COMPLEX, nbr_D, 40, & ! Recieve from left
                      comm0, status, ierr)
    CALL mpi_sendrecv(phi(ikxs:ikxe,ikys:ikye,ize-1), count, MPI_DOUBLE_COMPLEX, nbr_U, 41, & ! Send to right
                      phi(ikxs:ikxe,ikys:ikye,ize-2), count, MPI_DOUBLE_COMPLEX, nbr_D, 41, & ! Recieve from left
                      comm0, status, ierr)
                      !!!!!!!!!!! Send ghost to left neighbour !!!!!!!!!!!!!!!!!!!!!!

    CALL mpi_sendrecv(phi(ikxs:ikxe,ikys:ikye,  izs), count, MPI_DOUBLE_COMPLEX, nbr_U, 42, & ! Send to right
                      phi(ikxs:ikxe,ikys:ikye,izs+1), count, MPI_DOUBLE_COMPLEX, nbr_D, 42, & ! Recieve from left
                      comm0, status, ierr)
    CALL mpi_sendrecv(phi(ikxs:ikxe,ikys:ikye,izs+1), count, MPI_DOUBLE_COMPLEX, nbr_U, 43, & ! Send to right
                      phi(ikxs:ikxe,ikys:ikye,izs+2), count, MPI_DOUBLE_COMPLEX, nbr_D, 43, & ! Recieve from left
                      comm0, status, ierr)

  ELSE ! still need to perform periodic boundary conditions
    phi(:,:,   0) = phi(:,:,   Nz)
    phi(:,:,  -1) = phi(:,:, Nz-1)
    phi(:,:,Nz+1) = phi(:,:,    1)
    phi(:,:,Nz+2) = phi(:,:,    2)
  ENDIF
END SUBROUTINE update_ghosts_z_phi

END MODULE ghosts
