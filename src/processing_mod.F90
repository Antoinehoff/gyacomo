MODULE processing
    ! contains the Hermite-Laguerre collision operators. Solved using COSOlver.
    USE basic
    USE prec_const
    USE grid
    USE utility
    implicit none

    REAL(dp), PUBLIC, PROTECTED :: pflux_ri, gflux_ri

    PUBLIC :: compute_radial_ion_transport, compute_density, compute_temperature

CONTAINS

! 1D diagnostic to compute the average radial particle transport <n_i v_ExB>_r
SUBROUTINE compute_radial_ion_transport
    USE fields,           ONLY : moments_i, phi
    USE array,            ONLY : kernel_i
    USE time_integration, ONLY : updatetlevel
    IMPLICIT NONE
    COMPLEX(dp) :: pflux_local, gflux_local
    REAL(dp)    :: kz_, buffer(1:2)
    INTEGER     :: i_, world_rank, world_size, root

    pflux_local = 0._dp ! particle flux
    gflux_local = 0._dp ! gyrocenter flux
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

! Compute the 2D particle density for electron and ions (sum over Laguerre)
SUBROUTINE compute_density
  USE fields,           ONLY : moments_i, moments_e
  USE array,            ONLY : dens_e, dens_i, kernel_e, kernel_i
  USE time_integration, ONLY : updatetlevel
  IMPLICIT NONE

  IF( (ips_i .EQ. 1) .AND. (ips_e .EQ. 1) ) THEN
      ! Loop to compute dens_i = sum_j kernel_j Ni0j at each k
      DO ikz = ikzs,ikze
        DO ikr = ikrs,ikre
          ! electron density
          dens_e(ikr,ikz) = 0._dp
          DO ij = ijs_e, ije_e
              dens_e(ikr,ikz) = dens_e(ikr,ikz) + &
                    kernel_e(ij,ikr,ikz) * moments_e(1,ij,ikr,ikz,updatetlevel)
          ENDDO
          ! ion density
          dens_i(ikr,ikz) = 0._dp
          DO ij = ijs_i, ije_i
              dens_i(ikr,ikz) = dens_i(ikr,ikz) + &
                    kernel_i(ij,ikr,ikz) * moments_i(1,ij,ikr,ikz,updatetlevel)
          ENDDO
        ENDDO
      ENDDO
  ENDIF
  CALL manual_2D_bcast(dens_e(ikrs:ikre,ikzs:ikze))
  CALL manual_2D_bcast(dens_i(ikrs:ikre,ikzs:ikze))
END SUBROUTINE compute_density

! Compute the 2D particle temperature for electron and ions (sum over Laguerre)
SUBROUTINE compute_temperature
  USE fields,           ONLY : moments_i, moments_e
  USE array,            ONLY : temp_e, temp_i, kernel_e, kernel_i
  USE time_integration, ONLY : updatetlevel
  IMPLICIT NONE
  REAL(dp)    :: j_dp
  COMPLEX(dp) :: Tperp, Tpar

  IF( ((ips_i .EQ. 1) .AND. (ips_e .EQ. 1)) ) THEN
      ! Loop to compute T = 1/3*(Tpar + 2Tperp)
      DO ikz = ikzs,ikze
        DO ikr = ikrs,ikre
          ! electron temperature
          Tpar = 0._dp; Tperp = 0._dp
          DO ij = ijs_e, ije_e
            j_dp = REAL(ij-1,dp)
            Tpar = Tpar + kernel_e(ij,ikr,ikz)* &
              (SQRT2*moments_e(3,ij,ikr,ikz,updatetlevel) + moments_e(1,ij,ikr,ikz,updatetlevel))
            Tperp = Tperp + moments_e(1,ij,ikr,ikz,updatetlevel)*&
              ((2._dp*j_dp+1)*kernel_e(ij,ikr,ikz) - (j_dp+1)*kernel_e(ij+1,ikr,ikz) - j_dp*kernel_e(ij-1,ikr,ikz))
          ENDDO
          temp_e(ikr,ikz) = (Tpar + 2._dp*Tperp)/3._dp

          ! ion temperature
          Tpar = 0._dp; Tperp = 0._dp
          DO ij = ijs_i, ije_i
            j_dp = REAL(ij-1,dp)
            Tpar = Tpar + kernel_i(ij,ikr,ikz)* &
              (SQRT2*moments_i(3,ij,ikr,ikz,updatetlevel) + moments_i(1,ij,ikr,ikz,updatetlevel))
            Tperp = Tperp + moments_i(1,ij,ikr,ikz,updatetlevel)*&
              ((2._dp*j_dp+1)*kernel_i(ij,ikr,ikz) - (j_dp+1)*kernel_i(ij+1,ikr,ikz) - j_dp*kernel_i(ij-1,ikr,ikz))
          ENDDO
          temp_i(ikr,ikz) = (Tpar + 2._dp*Tperp)/3._dp
        ENDDO
      ENDDO
  ENDIF
  CALL manual_2D_bcast(temp_e(ikrs:ikre,ikzs:ikze))
  CALL manual_2D_bcast(temp_i(ikrs:ikre,ikzs:ikze))
END SUBROUTINE compute_temperature

END MODULE processing
