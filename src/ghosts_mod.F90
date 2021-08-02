module ghosts
USE basic
USE fields, ONLY : moments_e, moments_i
USE grid
USE time_integration

IMPLICIT NONE

INTEGER :: status(MPI_STATUS_SIZE), source, dest, count, ipg

PUBLIC :: update_ghosts

CONTAINS

SUBROUTINE update_ghosts
    CALL cpu_time(t0_ghost)

    IF (num_procs_p .GT. 1) THEN ! Do it only if we share the p
        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        CALL update_ghosts_p_e

        CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
        CALL update_ghosts_p_i
    ENDIF

    CALL cpu_time(t1_ghost)
    tc_ghost = tc_ghost + (t1_ghost - t0_ghost)
END SUBROUTINE update_ghosts


!Communicate p+1, p+2 moments to left neighboor and p-1, p-2 moments to right one
SUBROUTINE update_ghosts_p_e

    IMPLICIT NONE
    count = (ijeg_e-ijsg_e+1)*(ikxe-ikxs+1)*(ikye-ikys+1)*(ize-izs+1)

    !!!!!!!!!!! Send ghost to right neighbour !!!!!!!!!!!!!!!!!!!!!!
    CALL mpi_sendrecv(moments_e(ipe_e  ,ijsg_e:ijeg_e,ikxs:ikxe,ikys:ikye,izs:ize,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_R, 10, &
                      moments_e(ips_e-1,ijsg_e:ijeg_e,ikxs:ikxe,ikys:ikye,izs:ize,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_L, 10, &
                      comm0, status, ierr)
    CALL mpi_sendrecv(moments_e(ipe_e-1,ijsg_e:ijeg_e,ikxs:ikxe,ikys:ikye,izs:ize,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_R, 11, &
                      moments_e(ips_e-2,ijsg_e:ijeg_e,ikxs:ikxe,ikys:ikye,izs:ize,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_L, 11, &
                      comm0, status, ierr)

    !!!!!!!!!!! Send ghost to left neighbour !!!!!!!!!!!!!!!!!!!!!!
    CALL mpi_sendrecv(moments_e(ips_e  ,ijsg_e:ijeg_e,ikxs:ikxe,ikys:ikye,izs:ize,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_L, 12, &
                      moments_e(ipe_e+1,ijsg_e:ijeg_e,ikxs:ikxe,ikys:ikye,izs:ize,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_R, 12, &
                      comm0, status, ierr)
    CALL mpi_sendrecv(moments_e(ips_e+1,ijsg_e:ijeg_e,ikxs:ikxe,ikys:ikye,izs:ize,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_L, 13, &
                      moments_e(ipe_e+2,ijsg_e:ijeg_e,ikxs:ikxe,ikys:ikye,izs:ize,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_R, 13, &
                      comm0, status, ierr)

END SUBROUTINE update_ghosts_p_e

!Communicate p+1, p+2 moments to left neighboor and p-1, p-2 moments to right one
SUBROUTINE update_ghosts_p_i

    IMPLICIT NONE

    count = (ijeg_i-ijsg_i+1)*(ikxe-ikxs+1)*(ikye-ikys+1)*(ize-izs+1) ! Number of elements sent

    !!!!!!!!!!! Send ghost to right neighbour !!!!!!!!!!!!!!!!!!!!!!
    CALL mpi_sendrecv(moments_i(ipe_i  ,ijsg_i:ijeg_i,ikxs:ikxe,ikys:ikye,izs:ize,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_R, 14, &
                      moments_i(ips_i-1,ijsg_i:ijeg_i,ikxs:ikxe,ikys:ikye,izs:ize,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_L, 14, &
                      comm0, status, ierr)
    CALL mpi_sendrecv(moments_i(ipe_i-1,ijsg_i:ijeg_i,ikxs:ikxe,ikys:ikye,izs:ize,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_R, 15, &
                      moments_i(ips_i-2,ijsg_i:ijeg_i,ikxs:ikxe,ikys:ikye,izs:ize,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_L, 15, &
                      comm0, status, ierr)

    !!!!!!!!!!! Send ghost to left neighbour !!!!!!!!!!!!!!!!!!!!!!
    CALL mpi_cart_shift(comm0, 0, -1, source , dest , ierr)
    CALL mpi_sendrecv(moments_i(ips_i  ,ijsg_i:ijeg_i,ikxs:ikxe,ikys:ikye,izs:ize,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_L, 16, &
                      moments_i(ipe_i+1,ijsg_i:ijeg_i,ikxs:ikxe,ikys:ikye,izs:ize,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_R, 16, &
                      comm0, status, ierr)
    CALL mpi_sendrecv(moments_i(ips_i+1,ijsg_i:ijeg_i,ikxs:ikxe,ikys:ikye,izs:ize,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_L, 17, &
                      moments_i(ipe_i+2,ijsg_i:ijeg_i,ikxs:ikxe,ikys:ikye,izs:ize,updatetlevel), count, MPI_DOUBLE_COMPLEX, nbr_R, 17, &
                      comm0, status, ierr)

END SUBROUTINE update_ghosts_p_i

END MODULE ghosts
