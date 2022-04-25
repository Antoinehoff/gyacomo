MODULE processing
    ! contains the Hermite-Laguerre collision operators. Solved using COSOlver.
    USE basic
    USE prec_const
    USE grid
    USE geometry
    USE utility
    USE calculus
    implicit none

    REAL(dp), PUBLIC, PROTECTED :: pflux_ri, gflux_ri
    REAL(dp), PUBLIC, PROTECTED :: hflux_x

    PUBLIC :: compute_nadiab_moments
    PUBLIC :: compute_density, compute_upar, compute_uperp
    PUBLIC :: compute_Tpar, compute_Tperp, compute_fluid_moments
    PUBLIC :: compute_radial_ion_transport, compute_radial_heatflux
CONTAINS

! 1D diagnostic to compute the average radial particle transport <n_i v_ExB>_r
SUBROUTINE compute_radial_ion_transport
    USE fields,           ONLY : moments_i, phi
    USE array,            ONLY : kernel_i
    USE geometry,         ONLY : Jacobian
    USE time_integration, ONLY : updatetlevel
    IMPLICIT NONE
    COMPLEX(dp) :: pflux_local, gflux_local, integral
    REAL(dp)    :: ky_, buffer(1:2)
    INTEGER     :: i_, world_rank, world_size, root
    COMPLEX(dp), DIMENSION(izgs:izge) :: integrant

    pflux_local = 0._dp ! particle flux
    gflux_local = 0._dp ! gyrocenter flux
    integrant   = 0._dp ! auxiliary variable for z integration
    IF(ips_i .EQ. 1) THEN
      !! Gyro center flux (drift kinetic)
      DO iky = ikys,ikye
          ky_ = kyarray(iky)
          DO ikx = ikxs,ikxe
            DO iz = izgs,izge
              integrant(iz) = Jacobian(iz,0)*imagu * ky_ * moments_i(ip0_i,1,ikx,iky,iz,updatetlevel) * CONJG(phi(ikx,iky,iz))
            ENDDO
            ! Integrate over z
            call simpson_rule_z(integrant,integral)
            ! add contribution
            gflux_local = gflux_local + integral*iInt_Jacobian
        ENDDO
      ENDDO
      !! Particle flux
      DO iky = ikys,ikye
          ky_ = kyarray(iky)
          DO ikx = ikxs,ikxe
            integrant   = 0._dp ! auxiliary variable for z integration
              DO ij = ijs_i, ije_i
                DO iz = izgs,izge
                  integrant(iz) = integrant(iz) + &
                      Jacobian(iz,0)*imagu * ky_ * kernel_i(ij,ikx,iky,iz,0) &
                      *moments_i(ip0_i,ij,ikx,iky,iz,updatetlevel) * CONJG(phi(ikx,iky,iz))
                ENDDO
              ENDDO
              ! Integrate over z
              call simpson_rule_z(integrant,integral)
              ! add contribution
              pflux_local = pflux_local + integral*iInt_Jacobian
        ENDDO
      ENDDO
    ENDIF
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
    ELSE
      gflux_ri = gflux_local
      pflux_ri = pflux_local
    ENDIF
    ! if(my_id .eq. 0) write(*,*) 'pflux_ri = ',pflux_ri
END SUBROUTINE compute_radial_ion_transport

! 1D diagnostic to compute the average radial particle transport <n_i v_ExB>_r
SUBROUTINE compute_radial_heatflux
    USE fields,           ONLY : moments_i, moments_e, phi
    USE array,            ONLY : kernel_e, kernel_i
    USE geometry,         ONLY : Jacobian
    USE time_integration, ONLY : updatetlevel
    USE model, ONLY : qe_taue, qi_taui, KIN_E
    IMPLICIT NONE
    COMPLEX(dp) :: hflux_local, integral
    REAL(dp)    :: ky_, buffer(1:2), j_dp
    INTEGER     :: i_, world_rank, world_size, root
    COMPLEX(dp), DIMENSION(izgs:izge) :: integrant        ! charge density q_a n_a

    hflux_local = 0._dp ! particle flux
    IF(ips_i .EQ. 1 .AND. ips_e .EQ. 1) THEN
        ! Loop to compute gamma_kx = sum_ky sum_j -i k_z Kernel_j Ni00 * phi
        DO iky = ikys,ikye
          ky_ = kyarray(iky)
          DO ikx = ikxs,ikxe
            integrant = 0._dp
            DO ij = ijs_i, ije_i
              DO iz = izgs,izge
              integrant(iz) = integrant(iz) + Jacobian(iz,0)*imagu*ky_*CONJG(phi(ikx,iky,iz))&
               *(twothird * (   2._dp*j_dp  * kernel_i(ij  ,ikx,iky,iz,0) &
                                - (j_dp+1)  * kernel_i(ij+1,ikx,iky,iz,0) &
                                -    j_dp   * kernel_i(ij-1,ikx,iky,iz,0))&
               * (moments_i(ip0_i,ij,ikx,iky,iz,updatetlevel)+qi_taui*kernel_i(ij,ikx,iky,iz,0)*phi(ikx,iky,iz)) &
              + SQRT2*onethird * kernel_i(ij,ikx,iky,iz,0) * moments_i(ip2_i,ij,ikx,iky,iz,updatetlevel))
              ENDDO
            ENDDO
            IF(KIN_E) THEN
            DO ij = ijs_e, ije_e
              DO iz = izgs,izge
              integrant(iz) = integrant(iz) + Jacobian(iz,0)*imagu*ky_*CONJG(phi(ikx,iky,iz))&
               *(twothird * (   2._dp*j_dp  * kernel_e(ij  ,ikx,iky,iz,0) &
                                - (j_dp+1)  * kernel_e(ij+1,ikx,iky,iz,0) &
                                -    j_dp   * kernel_e(ij-1,ikx,iky,iz,0))&
               * (moments_e(ip0_e,ij,ikx,iky,iz,updatetlevel)+qe_taue*kernel_e(ij,ikx,iky,iz,0)*phi(ikx,iky,iz)) &
              + SQRT2*onethird * kernel_e(ij,ikx,iky,iz,0) * moments_e(ip2_e,ij,ikx,iky,iz,updatetlevel))
              ENDDO
            ENDDO
            ENDIF
            ! Integrate over z
            call simpson_rule_z(integrant,integral)
            hflux_local = hflux_local + integral*iInt_Jacobian
          ENDDO
        ENDDO
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
        ELSE
          hflux_x = hflux_local
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
  USE model,            ONLY : qe_taue, qi_taui, KIN_E
  implicit none

  ! Electron non adiab moments
  xloop: DO ikx = ikxs,ikxe
  yloop: DO iky = ikys,ikye
  zloop: DO iz  = izgs,izge
    IF(KIN_E) THEN
      DO ip=ipgs_e,ipge_e
        IF(parray_e(ip) .EQ. 0) THEN
          DO ij=ijgs_e,ijge_e
            nadiab_moments_e(ip,ij,ikx,iky,iz) = moments_e(ip,ij,ikx,iky,iz,updatetlevel) &
                                      + qe_taue*kernel_e(ij,ikx,iky,iz,0)*phi(ikx,iky,iz)
          ENDDO
        ELSE
          DO ij=ijgs_e,ijge_e
            nadiab_moments_e(ip,ij,ikx,iky,iz) = moments_e(ip,ij,ikx,iky,iz,updatetlevel)
          ENDDO
        ENDIF
      ENDDO
    ENDIF
    ! Ions non adiab moments
    DO ip=ipgs_i,ipge_i
      IF(parray_i(ip) .EQ. 0) THEN
        DO ij=ijgs_i,ijge_i
          nadiab_moments_i(ip,ij,ikx,iky,iz) = moments_i(ip,ij,ikx,iky,iz,updatetlevel) &
                                    + qi_taui*kernel_i(ij,ikx,iky,iz,0)*phi(ikx,iky,iz)
        ENDDO
      ELSE
        DO ij=ijgs_i,ijge_i
          nadiab_moments_i(ip,ij,ikx,iky,iz) = moments_i(ip,ij,ikx,iky,iz,updatetlevel)
        ENDDO
      ENDIF
    ENDDO
  ENDDO zloop
  ENDDO yloop
  ENDDO xloop
  !
END SUBROUTINE compute_nadiab_moments

! Compute the 2D particle density for electron and ions (sum over Laguerre)
SUBROUTINE compute_density
  USE array, ONLY : dens_e, dens_i, kernel_e, kernel_i, nadiab_moments_e, nadiab_moments_i
  USE model, ONLY : KIN_E
  IMPLICIT NONE
  COMPLEX(dp) :: dens

  IF ( CONTAINS_ip0_e .AND. CONTAINS_ip0_i ) THEN
      ! Loop to compute dens_i = sum_j kernel_j Ni0j at each k
      DO iky = ikys,ikye
        DO ikx = ikxs,ikxe
          DO iz = izs,ize
            IF(KIN_E) THEN
            ! electron density
            dens = 0._dp
            DO ij = ijs_e, ije_e
                dens = dens + kernel_e(ij,ikx,iky,iz,0) * nadiab_moments_e(ip0_e,ij,ikx,iky,iz)
            ENDDO
            dens_e(ikx,iky,iz) = dens
            ENDIF
            ! ion density
            dens = 0._dp
            DO ij = ijs_i, ije_i
                dens = dens + kernel_i(ij,ikx,iky,iz,0) * nadiab_moments_i(ip0_e,ij,ikx,iky,iz)
            ENDDO
            dens_i(ikx,iky,iz) = dens
          ENDDO
        ENDDO
      ENDDO
  ENDIF
  ! IF(KIN_E)&
  ! CALL manual_3D_bcast(dens_e(ikxs:ikxe,ikys:ikye,izs:ize))
  ! CALL manual_3D_bcast(dens_i(ikxs:ikxe,ikys:ikye,izs:ize))
END SUBROUTINE compute_density

! Compute the 2D particle fluid perp velocity for electron and ions (sum over Laguerre)
SUBROUTINE compute_uperp
  USE array, ONLY : uper_e, uper_i, kernel_e, kernel_i, nadiab_moments_e, nadiab_moments_i
  USE model, ONLY : KIN_E
  IMPLICIT NONE
  COMPLEX(dp) :: uperp

  IF ( CONTAINS_ip0_e .AND. CONTAINS_ip0_i ) THEN
      DO iky = ikys,ikye
        DO ikx = ikxs,ikxe
          DO iz = izs,ize
            IF(KIN_E) THEN
            ! electron
            uperp = 0._dp
            DO ij = ijs_e, ije_e
                uperp = uperp + kernel_e(ij,ikx,iky,iz,0) *&
                 0.5_dp*(nadiab_moments_e(ip0_e,ij,ikx,iky,iz) - nadiab_moments_e(ip0_e,ij-1,ikx,iky,iz))
            ENDDO
            uper_e(ikx,iky,iz) = uperp
            ENDIF
            ! ion
            uperp = 0._dp
            DO ij = ijs_i, ije_i
              uperp = uperp + kernel_i(ij,ikx,iky,iz,0) *&
               0.5_dp*(nadiab_moments_i(ip0_i,ij,ikx,iky,iz) - nadiab_moments_i(ip0_i,ij-1,ikx,iky,iz))
             ENDDO
            uper_i(ikx,iky,iz) = uperp
          ENDDO
        ENDDO
      ENDDO
  ENDIF
  ! IF(KIN_E)&
  ! CALL manual_3D_bcast(uper_e(ikxs:ikxe,ikys:ikye,izs:ize))
  ! CALL manual_3D_bcast(uper_i(ikxs:ikxe,ikys:ikye,izs:ize))
END SUBROUTINE compute_uperp

! Compute the 2D particle fluid par velocity for electron and ions (sum over Laguerre)
SUBROUTINE compute_upar
  USE array, ONLY : upar_e, upar_i, kernel_e, kernel_i, nadiab_moments_e, nadiab_moments_i
  USE model, ONLY : KIN_E
  IMPLICIT NONE
  COMPLEX(dp) :: upar

  IF ( CONTAINS_ip1_e .AND. CONTAINS_ip1_i ) THEN
      DO iky = ikys,ikye
        DO ikx = ikxs,ikxe
          DO iz = izs,ize
            IF(KIN_E) THEN
            ! electron
            upar = 0._dp
            DO ij = ijs_e, ije_e
              upar = upar + kernel_e(ij,ikx,iky,iz,1)*nadiab_moments_e(ip1_e,ij,ikx,iky,iz)
            ENDDO
            upar_e(ikx,iky,iz) = upar
            ENDIF
            ! ion
            upar = 0._dp
            DO ij = ijs_i, ije_i
              upar = upar + kernel_i(ij,ikx,iky,iz,1)*nadiab_moments_i(ip1_i,ij,ikx,iky,iz)
             ENDDO
            upar_i(ikx,iky,iz) = upar
          ENDDO
        ENDDO
      ENDDO
  ELSE
    IF(KIN_E)&
    upar_e = 0
    upar_i = 0
  ENDIF
  ! IF(KIN_E)&
  ! CALL manual_3D_bcast(upar_e(ikxs:ikxe,ikys:ikye,izs:ize))
  ! CALL manual_3D_bcast(upar_i(ikxs:ikxe,ikys:ikye,izs:ize))
END SUBROUTINE compute_upar

! Compute the 2D particle temperature for electron and ions (sum over Laguerre)
SUBROUTINE compute_tperp
  USE array, ONLY : Tper_e, Tper_i, kernel_e, kernel_i, nadiab_moments_e, nadiab_moments_i
  USE model, ONLY : KIN_E
  IMPLICIT NONE
  REAL(dp)    :: j_dp
  COMPLEX(dp) :: Tperp

  IF ( CONTAINS_ip0_e .AND. CONTAINS_ip0_i .AND. &
       CONTAINS_ip2_e .AND. CONTAINS_ip2_i ) THEN
      ! Loop to compute T = 1/3*(Tpar + 2Tperp)
      DO iky = ikys,ikye
        DO ikx = ikxs,ikxe
          DO iz = izs,ize
            ! electron temperature
            IF(KIN_E) THEN
            Tperp  = 0._dp
            DO ij = ijs_e, ije_e
              j_dp = REAL(ij-1,dp)
              Tperp = Tperp + kernel_e(ij,ikx,iky,iz,0)*&
                  ((2_dp*j_dp+1)*nadiab_moments_e(ip0_e,ij  ,ikx,iky,iz)&
                  -j_dp         *nadiab_moments_e(ip0_e,ij-1,ikx,iky,iz)&
                  -j_dp+1       *nadiab_moments_e(ip0_e,ij+1,ikx,iky,iz))
            ENDDO
            Tper_e(ikx,iky,iz) = Tperp
            ENDIF
            ! ion temperature
            Tperp = 0._dp
            DO ij = ijs_i, ije_i
              j_dp = REAL(ij-1,dp)
              Tperp = Tperp + kernel_i(ij,ikx,iky,iz,0)*&
                  ((2_dp*j_dp+1)*nadiab_moments_i(ip0_i,ij  ,ikx,iky,iz)&
                  -j_dp         *nadiab_moments_i(ip0_i,ij-1,ikx,iky,iz)&
                  -j_dp+1       *nadiab_moments_i(ip0_i,ij+1,ikx,iky,iz))
            ENDDO
            Tper_i(ikx,iky,iz) = Tperp
          ENDDO
        ENDDO
      ENDDO
  ENDIF
  ! IF(KIN_E)&
  ! CALL manual_3D_bcast(Tper_e(ikxs:ikxe,ikys:ikye,izs:ize))
  ! CALL manual_3D_bcast(Tper_i(ikxs:ikxe,ikys:ikye,izs:ize))
END SUBROUTINE compute_Tperp

! Compute the 2D particle temperature for electron and ions (sum over Laguerre)
SUBROUTINE compute_Tpar
  USE array, ONLY : Tpar_e, Tpar_i, kernel_e, kernel_i, nadiab_moments_e, nadiab_moments_i
  USE model, ONLY : KIN_E
  IMPLICIT NONE
  REAL(dp)    :: j_dp
  COMPLEX(dp) :: tpar

  IF ( CONTAINS_ip0_e .AND. CONTAINS_ip0_i .AND. &
       CONTAINS_ip2_e .AND. CONTAINS_ip2_i ) THEN
      ! Loop to compute T = 1/3*(Tpar + 2Tperp)
      DO iky = ikys,ikye
        DO ikx = ikxs,ikxe
          DO iz = izs,ize
            ! electron temperature
            IF(KIN_E) THEN
            Tpar  = 0._dp
            DO ij = ijs_e, ije_e
              j_dp = REAL(ij-1,dp)
              Tpar  = Tpar + kernel_e(ij,ikx,iky,iz,0)*&
               (SQRT2 * nadiab_moments_e(ip2_e,ij,ikx,iky,iz) &
                      + nadiab_moments_e(ip0_e,ij,ikx,iky,iz))
            ENDDO
            Tpar_e(ikx,iky,iz) = Tpar
            ENDIF
            ! ion temperature
            Tpar  = 0._dp
            DO ij = ijs_i, ije_i
              j_dp = REAL(ij-1,dp)
              Tpar  = Tpar + kernel_i(ij,ikx,iky,iz,0)*&
               (SQRT2 * nadiab_moments_i(ip2_i,ij,ikx,iky,iz) &
                      + nadiab_moments_i(ip0_i,ij,ikx,iky,iz))
            ENDDO
            Tpar_i(ikx,iky,iz) = Tpar
          ENDDO
        ENDDO
      ENDDO
  ENDIF
  ! IF(KIN_E)&
  ! CALL manual_3D_bcast(Tpar_e(ikxs:ikxe,ikys:ikye,izs:ize))
  ! CALL manual_3D_bcast(Tpar_i(ikxs:ikxe,ikys:ikye,izs:ize))
END SUBROUTINE compute_Tpar

! Compute the 2D particle fluid moments for electron and ions (sum over Laguerre)
SUBROUTINE compute_fluid_moments
  USE array, ONLY : dens_e, Tpar_e, Tper_e, dens_i, Tpar_i, Tper_i, temp_e, temp_i
  USE model, ONLY : KIN_E
  IMPLICIT NONE
  CALL compute_density
  CALL compute_upar
  CALL compute_uperp
  CALL compute_Tpar
  CALL compute_Tperp
  ! Temperature
  temp_e = (Tpar_e - 2._dp * Tper_e)/3._dp - dens_e
  temp_i = (Tpar_i - 2._dp * Tper_i)/3._dp - dens_i
END SUBROUTINE compute_fluid_moments

! Compute the 2D particle temperature for electron and ions (sum over Laguerre)
! SUBROUTINE compute_temperature
!   USE array, ONLY : temp_e, temp_i, kernel_e, kernel_i, nadiab_moments_e, nadiab_moments_i
!   USE model, ONLY : KIN_E
!   IMPLICIT NONE
!   REAL(dp)    :: j_dp
!   COMPLEX(dp) :: Tperp, Tpar, dens
!
!   IF ( CONTAINS_ip0_e .AND. CONTAINS_ip0_i .AND. &
!        CONTAINS_ip2_e .AND. CONTAINS_ip2_i ) THEN
!       ! Loop to compute T = 1/3*(Tpar + 2Tperp)
!       DO iky = ikys,ikye
!         DO ikx = ikxs,ikxe
!           DO iz = izs,ize
!             ! electron temperature
!             IF(KIN_E) THEN
!             dens  = 0._dp
!             Tpar  = 0._dp
!             Tperp = 0._dp
!             DO ij = ijs_e, ije_e
!               j_dp = REAL(ij-1,dp)
!               temp_e(ikx,iky,iz) = temp_e(ikx,iky,iz) + &
!                 2._dp/3._dp * (2._dp*j_dp*kernel_e(ij,ikx,iky,iz,0) - (j_dp+1)*kernel_e(ij+1,ikx,iky,iz,0) - j_dp*kernel_e(ij-1,ikx,iky,iz,0))&
!                  * (moments_e(ip0_e,ij,ikx,iky,iz,updatetlevel)+q_e/tau_e*kernel_e(ij,ikx,iky,iz,0)*phi(ikx,iky,iz)) &
!                 + SQRT2/3._dp * kernel_e(ij,ikx,iky,iz,0) * moments_e(ip2_e,ij,ikx,iky,iz,updatetlevel)
!             ENDDO
!             temp_e(ikx,iky,iz) = (Tpar - 2._dp * Tperp)/3._dp - dens
!             ENDIF
!             ! ion temperature
!             dens  = 0._dp
!             Tpar  = 0._dp
!             Tperp = 0._dp
!             DO ij = ijs_i, ije_i
!               j_dp = REAL(ij-1,dp)
!               temp_i(ikx,iky,iz) = temp_i(ikx,iky,iz) + &
!                 2._dp/3._dp * (2._dp*j_dp*kernel_i(ij,ikx,iky,iz,0) - (j_dp+1)*kernel_i(ij+1,ikx,iky,iz,0) - j_dp*kernel_i(ij-1,ikx,iky,iz,0))&
!                  * (moments_i(ip0_i,ij,ikx,iky,iz,updatetlevel)+q_i/tau_i*kernel_i(ij,ikx,iky,iz,0)*phi(ikx,iky,iz)) &
!                 + SQRT2/3._dp * kernel_i(ij,ikx,iky,iz,0) * moments_i(ip2_i,ij,ikx,iky,iz,updatetlevel)
!             ENDDO
!             temp_i(ikx,iky,iz) = (Tpar - 2._dp * Tperp)/3._dp - dens
!           ENDDO
!         ENDDO
!       ENDDO
!   ENDIF
!   IF(KIN_E) CALL manual_3D_bcast(temp_e(ikxs:ikxe,ikys:ikye,izs:ize))
!   CALL manual_3D_bcast(temp_i(ikxs:ikxe,ikys:ikye,izs:ize))
! END SUBROUTINE compute_temperature


END MODULE processing
