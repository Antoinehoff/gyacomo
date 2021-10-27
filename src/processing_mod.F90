MODULE processing
    ! contains the Hermite-Laguerre collision operators. Solved using COSOlver.
    USE basic
    USE prec_const
    USE grid
    USE utility
    implicit none

    REAL(dp), PUBLIC, PROTECTED :: pflux_ri, gflux_ri
    REAL(dp), PUBLIC, PROTECTED :: hflux_x

    PUBLIC :: compute_density, compute_temperature, compute_nadiab_moments
    PUBLIC :: compute_radial_ion_transport, compute_radial_heatflux
CONTAINS

! 1D diagnostic to compute the average radial particle transport <n_i v_ExB>_r
SUBROUTINE compute_radial_ion_transport
    USE fields,           ONLY : moments_i, phi
    USE array,            ONLY : kernel_i
    USE time_integration, ONLY : updatetlevel
    IMPLICIT NONE
    COMPLEX(dp) :: pflux_local, gflux_local
    REAL(dp)    :: ky_, buffer(1:2)
    INTEGER     :: i_, world_rank, world_size, root

    pflux_local = 0._dp ! particle flux
    gflux_local = 0._dp ! gyrocenter flux
    IF(ips_i .EQ. 1) THEN
        ! Loop to compute gamma_kx = sum_ky sum_j -i k_z Kernel_j Ni00 * phi
        DO iz = izs,ize
          DO iky = ikys,ikye
              ky_ = kyarray(iky)
              DO ikx = ikxs,ikxe
                  gflux_local = gflux_local + &
                      imagu * ky_ * moments_i(ip0_i,1,ikx,iky,iz,updatetlevel) * CONJG(phi(ikx,iky,iz))
                  DO ij = ijs_i, ije_i
                      pflux_local = pflux_local + &
                          imagu * ky_ * kernel_i(ij,ikx,iky,iz) * moments_i(ip0_i,ij,ikx,iky,iz,updatetlevel) * CONJG(phi(ikx,iky,iz))
                  ENDDO
              ENDDO
          ENDDO
        ENDDO
        gflux_local = gflux_local/Nz ! Average over parallel planes
        pflux_local = pflux_local/Nz

        buffer(1) = REAL(gflux_local)
        buffer(2) = REAL(pflux_local)
        root = 0
        !Gather manually among the rank_p=0 processes and perform the sum
        gflux_ri = 0
        pflux_ri = 0
        IF (num_procs_kx .GT. 1) THEN
            !! Everyone sends its local_sum to root = 0
            IF (rank_kx .NE. root) THEN
                CALL MPI_SEND(buffer, 2 , MPI_DOUBLE_PRECISION, root, 1234, comm_kx, ierr)
            ELSE
                ! Recieve from all the other processes
                DO i_ = 0,num_procs_kx-1
                    IF (i_ .NE. rank_kx) &
                        CALL MPI_RECV(buffer, 2 , MPI_DOUBLE_PRECISION, i_, 1234, comm_kx, MPI_STATUS_IGNORE, ierr)
                        gflux_ri = gflux_ri + buffer(1)
                        pflux_ri = pflux_ri + buffer(2)
                ENDDO
            ENDIF
        ENDIF
    ENDIF
END SUBROUTINE compute_radial_ion_transport

! 1D diagnostic to compute the average radial particle transport <n_i v_ExB>_r
SUBROUTINE compute_radial_heatflux
    USE fields,           ONLY : moments_i, moments_e, phi
    USE array,            ONLY : dens_e, dens_i, kernel_e, kernel_i
    USE time_integration, ONLY : updatetlevel
    USE model, ONLY : q_e, q_i, tau_e, tau_i
    IMPLICIT NONE
    COMPLEX(dp) :: hflux_local
    REAL(dp)    :: ky_, buffer(1:2), j_dp
    INTEGER     :: i_, world_rank, world_size, root

    hflux_local = 0._dp ! particle flux
    IF(ips_i .EQ. 1 .AND. ips_e .EQ. 1) THEN
        ! Loop to compute gamma_kx = sum_ky sum_j -i k_z Kernel_j Ni00 * phi
        DO iz = izs,ize
          DO iky = ikys,ikye
            ky_ = kyarray(iky)
            DO ikx = ikxs,ikxe
              DO ij = ijs_i, ije_i
                j_dp = REAL(ij-1,dp)
                hflux_local = hflux_local - imagu*ky_*CONJG(phi(ikx,iky,iz))&
                 *(2._dp/3._dp * (2._dp*j_dp*kernel_i(ij,ikx,iky,iz) - (j_dp+1)*kernel_i(ij+1,ikx,iky,iz) - j_dp*kernel_i(ij-1,ikx,iky,iz))&
                 * (moments_i(ip0_i,ij,ikx,iky,iz,updatetlevel)+q_i/tau_i*kernel_i(ij,ikx,iky,iz)*phi(ikx,iky,iz)) &
                + SQRT2/3._dp * kernel_i(ij,ikx,iky,iz) * moments_i(ip2_i,ij,ikx,iky,iz,updatetlevel))
              ENDDO
              DO ij = ijs_e, ije_e
                j_dp = REAL(ij-1,dp)
                hflux_local = hflux_local - imagu*ky_*CONJG(phi(ikx,iky,iz))&
                 *(2._dp/3._dp *(2._dp*j_dp * kernel_e(  ij,ikx,iky,iz) &
                                  -(j_dp+1)  * kernel_e(ij+1,ikx,iky,iz) &
                                  -j_dp      * kernel_e(ij-1,ikx,iky,iz))&
                               *(moments_e(ip0_e,ij,ikx,iky,iz,updatetlevel)&
                                  +q_e/tau_e * kernel_e(  ij,ikx,iky,iz) * phi(ikx,iky,iz)) &
                  +SQRT2/3._dp * kernel_e(ij,ikx,iky,iz) *                                                                                                     moments_e(ip2_e,ij,ikx,iky,iz,updatetlevel))
              ENDDO
            ENDDO
          ENDDO
        ENDDO
        hflux_local = hflux_local/Nz

        buffer(2) = REAL(hflux_local)
        root = 0
        !Gather manually among the rank_p=0 processes and perform the sum
        hflux_x = 0
        IF (num_procs_kx .GT. 1) THEN
            !! Everyone sends its local_sum to root = 0
            IF (rank_kx .NE. root) THEN
                CALL MPI_SEND(buffer, 2 , MPI_DOUBLE_PRECISION, root, 1234, comm_kx, ierr)
            ELSE
                ! Recieve from all the other processes
                DO i_ = 0,num_procs_kx-1
                    IF (i_ .NE. rank_kx) &
                        CALL MPI_RECV(buffer, 2 , MPI_DOUBLE_PRECISION, i_, 1234, comm_kx, MPI_STATUS_IGNORE, ierr)
                        hflux_x = hflux_x + buffer(2)
                ENDDO
            ENDIF
        ENDIF
    ENDIF
END SUBROUTINE compute_radial_heatflux

SUBROUTINE compute_nadiab_moments
  ! evaluate the non-adiabatique ion moments
  !
  ! n_{pi} = N^{pj} + kernel(j) /tau_i phi delta_p0
  !
  USE fields,           ONLY : moments_i, moments_e, phi
  USE array,            ONLY : kernel_e, kernel_i, nadiab_moments_e, nadiab_moments_i
  USE time_integration, ONLY : updatetlevel
  USE model,            ONLY : qe_taue, qi_taui
  implicit none

  ! Add non-adiabatique term
  DO ip=ipsg_e,ipeg_e
    IF(parray_e(ip) .EQ. 0) THEN
      DO ij=ijsg_e,ijeg_e
        nadiab_moments_e(ip,ij,ikxs:ikxe,ikys:ikye,izs:ize)&
         = moments_e(ip,ij,ikxs:ikxe,ikys:ikye,izs:ize,updatetlevel) &
           + qe_taue*kernel_e(ij,ikxs:ikxe,ikys:ikye,izs:ize)*phi(ikx,iky,iz)
      ENDDO
    ELSE
      nadiab_moments_e(ip,ijsg_e:ijeg_e,ikxs:ikxe,ikys:ikye,izs:ize) &
      = moments_e(ip,ijsg_e:ijeg_e,ikxs:ikxe,ikys:ikye,izs:ize,updatetlevel)
    ENDIF
  ENDDO
  ! Add non-adiabatique term
  DO ip=ipsg_i,ipeg_i
    IF(parray_i(ip) .EQ. 0) THEN
      DO ij=ijsg_i,ijeg_i
        nadiab_moments_i(ip,ij,ikxs:ikxe,ikys:ikye,izs:ize)&
         = moments_i(ip,ij,ikxs:ikxe,ikys:ikye,izs:ize,updatetlevel) &
           + qi_taui*kernel_i(ij,ikxs:ikxe,ikys:ikye,izs:ize)*phi(ikx,iky,iz)
      ENDDO
    ELSE
      nadiab_moments_i(ip,ijsg_i:ijeg_i,ikxs:ikxe,ikys:ikye,izs:ize) &
      = moments_i(ip,ijsg_i:ijeg_i,ikxs:ikxe,ikys:ikye,izs:ize,updatetlevel)
    ENDIF
  ENDDO
  !
END SUBROUTINE compute_nadiab_moments

! Compute the 2D particle density for electron and ions (sum over Laguerre)
SUBROUTINE compute_density
  USE fields,           ONLY : moments_i, moments_e, phi
  USE array,            ONLY : dens_e, dens_i, kernel_e, kernel_i
  USE time_integration, ONLY : updatetlevel
  USE model, ONLY : q_e, q_i, tau_e, tau_i
  IMPLICIT NONE

  IF( (ips_i .EQ. 1) .AND. (ips_e .EQ. 1) ) THEN
      ! Loop to compute dens_i = sum_j kernel_j Ni0j at each k
      DO iky = ikys,ikye
        DO ikx = ikxs,ikxe
          DO iz = izs,ize
            ! electron density
            dens_e(ikx,iky,iz) = 0._dp
            DO ij = ijs_e, ije_e
                dens_e(ikx,iky,iz) = dens_e(ikx,iky,iz) + kernel_e(ij,ikx,iky,iz) * &
                 (moments_e(ip0_e,ij,ikx,iky,iz,updatetlevel)+q_e/tau_e*kernel_e(ij,ikx,iky,iz)*phi(ikx,iky,iz))
            ENDDO
            ! ion density
            dens_i(ikx,iky,iz) = 0._dp
            DO ij = ijs_i, ije_i
                dens_i(ikx,iky,iz) = dens_i(ikx,iky,iz) + kernel_i(ij,ikx,iky,iz) * &
                (moments_i(ip0_i,ij,ikx,iky,iz,updatetlevel)+q_i/tau_i*kernel_i(ij,ikx,iky,iz)*phi(ikx,iky,iz))
            ENDDO
          ENDDO
        ENDDO
      ENDDO
  ENDIF
  CALL manual_3D_bcast(dens_e(ikxs:ikxe,ikys:ikye,izs:ize))
  CALL manual_3D_bcast(dens_i(ikxs:ikxe,ikys:ikye,izs:ize))
END SUBROUTINE compute_density

! Compute the 2D particle temperature for electron and ions (sum over Laguerre)
SUBROUTINE compute_temperature
  USE fields,           ONLY : moments_i, moments_e, phi
  USE array,            ONLY : temp_e, temp_i, kernel_e, kernel_i
  USE time_integration, ONLY : updatetlevel
  USE model, ONLY : q_e, q_i, tau_e, tau_i
  IMPLICIT NONE
  REAL(dp)    :: j_dp
  COMPLEX(dp) :: Tperp, Tpar

  IF( ((ips_i .EQ. 1) .AND. (ips_e .EQ. 1)) ) THEN
      ! Loop to compute T = 1/3*(Tpar + 2Tperp)
      DO iky = ikys,ikye
        DO ikx = ikxs,ikxe
          DO iz = izs,ize
            ! electron temperature
            temp_e(ikx,iky,iz) = 0._dp
            DO ij = ijs_e, ije_e
              j_dp = REAL(ij-1,dp)
              temp_e(ikx,iky,iz) = temp_e(ikx,iky,iz) + &
                2._dp/3._dp * (2._dp*j_dp*kernel_e(ij,ikx,iky,iz) - (j_dp+1)*kernel_e(ij+1,ikx,iky,iz) - j_dp*kernel_e(ij-1,ikx,iky,iz))&
                 * (moments_e(ip0_e,ij,ikx,iky,iz,updatetlevel)+q_e/tau_e*kernel_e(ij,ikx,iky,iz)*phi(ikx,iky,iz)) &
                + SQRT2/3._dp * kernel_e(ij,ikx,iky,iz) * moments_e(ip2_e,ij,ikx,iky,iz,updatetlevel)
            ENDDO

            ! ion temperature
            temp_i(ikx,iky,iz) = 0._dp
            DO ij = ijs_i, ije_i
              j_dp = REAL(ij-1,dp)
              temp_i(ikx,iky,iz) = temp_i(ikx,iky,iz) + &
                2._dp/3._dp * (2._dp*j_dp*kernel_i(ij,ikx,iky,iz) - (j_dp+1)*kernel_i(ij+1,ikx,iky,iz) - j_dp*kernel_i(ij-1,ikx,iky,iz))&
                 * (moments_i(ip0_i,ij,ikx,iky,iz,updatetlevel)+q_i/tau_i*kernel_i(ij,ikx,iky,iz)*phi(ikx,iky,iz)) &
                + SQRT2/3._dp * kernel_i(ij,ikx,iky,iz) * moments_i(ip2_i,ij,ikx,iky,iz,updatetlevel)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
  ENDIF
  CALL manual_3D_bcast(temp_e(ikxs:ikxe,ikys:ikye,izs:ize))
  CALL manual_3D_bcast(temp_i(ikxs:ikxe,ikys:ikye,izs:ize))
END SUBROUTINE compute_temperature

END MODULE processing
