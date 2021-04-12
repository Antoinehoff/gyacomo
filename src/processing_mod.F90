MODULE processing
    ! contains the Hermite-Laguerre collision operators. Solved using COSOlver.
    USE basic
    USE prec_const
    USE grid
    implicit none

    REAL(dp), PUBLIC, PROTECTED :: pflux_ri, gflux_ri
    
    PUBLIC :: compute_radial_ion_transport
    
CONTAINS

SUBROUTINE compute_radial_ion_transport
    USE fields,           ONLY : moments_i, phi
    USE array,            ONLY : kernel_i
    USE time_integration, ONLY : updatetlevel
    IMPLICIT NONE
    COMPLEX(dp) :: pflux_local, gflux_local
    REAL(dp)    :: kz_, buffer(1:2)
    INTEGER     :: i_, world_rank, world_size, root

    pflux_local = 0._dp
    gflux_local = 0._dp
    IF(ips_i .EQ. 1) THEN
        ! Loop to compute gamma_kr = sum_kz sum_j -i k_z Kernel_j Ni00 * phi
        DO ikz = ikzs,ikze
            kz_ = kzarray(ikz)
            DO ikr = ikrs,ikre
                gflux_local = gflux_local - &
                    imagu * kz_ * moments_i(1,1,ikr,ikz,updatetlevel) * CONJG(phi(ikr,ikz))
                DO ij = ijs_i, ije_i
                    pflux_local = pflux_local - &
                        imagu * kz_ * kernel_i(ij,ikr,ikz) * moments_i(1,ij,ikr,ikz,updatetlevel) * CONJG(phi(ikr,ikz))
                ENDDO
            ENDDO
        ENDDO

        buffer(1) = REAL(gflux_local)
        buffer(2) = REAL(pflux_local)

        root = 0
        !Gather manually among the rank_p=0 processes and perform the sum
        gflux_ri = 0
        pflux_ri = 0
        IF (num_procs_kr .GT. 1) THEN
            !! Everyone sends its local_sum to root = 0
            IF (rank_kr .NE. root) THEN
                CALL MPI_SEND(buffer, 2 , MPI_DOUBLE_PRECISION, root, 1234, comm_kr, ierr)
            ELSE
                ! Recieve from all the other processes
                DO i_ = 0,num_procs_kr-1
                    IF (i_ .NE. rank_kr) &
                        CALL MPI_RECV(buffer, 2 , MPI_DOUBLE_PRECISION, i_, 1234, comm_kr, MPI_STATUS_IGNORE, ierr)
                        gflux_ri = gflux_ri + buffer(1)
                        pflux_ri = pflux_ri + buffer(2)
                ENDDO
            ENDIF
        ENDIF
    ENDIF


END SUBROUTINE compute_radial_ion_transport

END MODULE processing