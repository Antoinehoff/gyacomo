MODULE processing
    USE basic
    USE prec_const
    USE grid
    implicit none

    REAL(dp), PUBLIC, PROTECTED :: pflux_ri, gflux_ri, pflux_re, gflux_re
    REAL(dp), PUBLIC, PROTECTED :: hflux_xi, hflux_xe

    PUBLIC :: compute_nadiab_moments_z_gradients_and_interp
    PUBLIC :: compute_density, compute_upar, compute_uperp
    PUBLIC :: compute_Tpar, compute_Tperp, compute_fluid_moments
    PUBLIC :: compute_radial_ion_transport, compute_radial_electron_transport
    PUBLIC :: compute_radial_ion_heatflux, compute_radial_electron_heatflux
    PUBLIC :: compute_Napjz_spectrum
CONTAINS

! 1D diagnostic to compute the average radial particle transport <n_i v_ExB_x>_xyz
SUBROUTINE compute_radial_ion_transport
  USE fields,           ONLY : moments_i, phi, psi
  USE array,            ONLY : kernel_i
  USE geometry,         ONLY : Jacobian, iInt_Jacobian
  USE time_integration, ONLY : updatetlevel
  USE calculus,         ONLY : simpson_rule_z
  USE model,            ONLY : sqrt_tau_o_sigma_i, EM
  IMPLICIT NONE
  COMPLEX(dp) :: pflux_local, gflux_local, integral
  REAL(dp)    :: ky_, buffer(1:2)
  INTEGER     :: i_, world_rank, world_size, root
  COMPLEX(dp), DIMENSION(izgs:izge) :: integrant

  pflux_local = 0._dp ! particle flux
  gflux_local = 0._dp ! gyrocenter flux
  integrant   = 0._dp ! auxiliary variable for z integration
  !!---------- Gyro center flux (drift kinetic) ------------
  ! Electrostatic part
  IF(CONTAINS_ip0_i) THEN
    DO iky = ikys,ikye
        ky_ = kyarray(iky)
        DO ikx = ikxs,ikxe
          integrant(izgs:izge) = integrant(izgs:izge) &
           +moments_i(ip0_i,ij0_i,iky,ikx,izgs:izge,updatetlevel)&
            *imagu*ky_*CONJG(phi(iky,ikx,izgs:izge))
      ENDDO
    ENDDO
  ENDIF
  ! Electromagnetic part
  IF( EM .AND. CONTAINS_ip1_i ) THEN
    DO iky = ikys,ikye
        ky_ = kyarray(iky)
        DO ikx = ikxs,ikxe
          integrant(izgs:izge) = integrant(izgs:izge)&
            -sqrt_tau_o_sigma_i*moments_i(ip1_i,ij0_i,iky,ikx,izgs:izge,updatetlevel)&
             *imagu*ky_*CONJG(psi(iky,ikx,izgs:izge))
      ENDDO
    ENDDO
  ENDIF
  ! Integrate over z
  integrant(izgs:izge) = Jacobian(izgs:izge,0)*integrant(izgs:izge)
  call simpson_rule_z(integrant,integral)
  ! Get process local gyrocenter flux
  gflux_local = integral*iInt_Jacobian

  !
  integrant   = 0._dp ! reset auxiliary variable
  !!---------- Particle flux (gyrokinetic) ------------
  ! Electrostatic part
  IF(CONTAINS_ip0_i) THEN
    DO iky = ikys,ikye
        ky_ = kyarray(iky)
        DO ikx = ikxs,ikxe
          DO ij = ijs_i, ije_i
            integrant(izgs:izge) = integrant(izgs:izge)&
              +moments_i(ip0_i,ij,iky,ikx,izgs:izge,updatetlevel)&
              *imagu*ky_*kernel_i(ij,iky,ikx,izgs:izge,0)*CONJG(phi(iky,ikx,izgs:izge))
          ENDDO
      ENDDO
    ENDDO
  ENDIF
  ! Electromagnetic part
  IF( EM .AND. CONTAINS_ip1_i ) THEN
    DO iky = ikys,ikye
      ky_ = kyarray(iky)
      DO ikx = ikxs,ikxe
        integrant   = 0._dp ! auxiliary variable for z integration
        DO ij = ijs_i, ije_i
          integrant(izgs:izge) = integrant(izgs:izge)&
            -sqrt_tau_o_sigma_i*moments_i(ip1_i,ij,iky,ikx,izgs:izge,updatetlevel)&
            *imagu*ky_*kernel_i(ij,iky,ikx,izgs:izge,0)*CONJG(psi(iky,ikx,izgs:izge))
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  ! Integrate over z
  integrant(izgs:izge) = Jacobian(izgs:izge,0)*integrant(izgs:izge)
  call simpson_rule_z(integrant,integral)
  ! Get process local particle flux
  pflux_local = integral*iInt_Jacobian

  !!!!---------- Sum over all processes ----------
  buffer(1) = 2._dp*REAL(gflux_local)
  buffer(2) = 2._dp*REAL(pflux_local)
  root = 0
  !Gather manually among the rank_p=0 processes and perform the sum
  gflux_ri = 0
  pflux_ri = 0
  IF (num_procs_ky .GT. 1) THEN
      !! Everyone sends its local_sum to root = 0
      IF (rank_ky .NE. root) THEN
          CALL MPI_SEND(buffer, 2 , MPI_DOUBLE_PRECISION, root, 1234, comm_ky, ierr)
      ELSE
          ! Recieve from all the other processes
          DO i_ = 0,num_procs_ky-1
              IF (i_ .NE. rank_ky) &
                  CALL MPI_RECV(buffer, 2 , MPI_DOUBLE_PRECISION, i_, 1234, comm_ky, MPI_STATUS_IGNORE, ierr)
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

! 1D diagnostic to compute the average radial particle transport <n_e v_ExB_x>_xyz
SUBROUTINE compute_radial_electron_transport
  USE fields,           ONLY : moments_e, phi, psi
  USE array,            ONLY : kernel_e
  USE geometry,         ONLY : Jacobian, iInt_Jacobian
  USE time_integration, ONLY : updatetlevel
  USE calculus,         ONLY : simpson_rule_z
  USE model,            ONLY : sqrt_tau_o_sigma_e, EM
  IMPLICIT NONE
  COMPLEX(dp) :: pflux_local, gflux_local, integral
  REAL(dp)    :: ky_, buffer(1:2)
  INTEGER     :: i_, world_rank, world_size, root
  COMPLEX(dp), DIMENSION(izgs:izge) :: integrant

  pflux_local = 0._dp ! particle flux
  gflux_local = 0._dp ! gyrocenter flux
  integrant   = 0._dp ! auxiliary variable for z integration
  !!---------- Gyro center flux (drift kinetic) ------------
  ! Electrostatic part
  IF(CONTAINS_ip0_e) THEN
    DO iky = ikys,ikye
        ky_ = kyarray(iky)
        DO ikx = ikxs,ikxe
          integrant(izgs:izge) = integrant(izgs:izge) &
           +moments_e(ip0_e,ij0_e,iky,ikx,izgs:izge,updatetlevel)&
            *imagu*ky_*CONJG(phi(iky,ikx,izgs:izge))
      ENDDO
    ENDDO
  ENDIF
  ! Electromagnetic part
  IF( EM .AND. CONTAINS_ip1_e ) THEN
    DO iky = ikys,ikye
        ky_ = kyarray(iky)
        DO ikx = ikxs,ikxe
          integrant(izgs:izge) = integrant(izgs:izge)&
            -sqrt_tau_o_sigma_e*moments_e(ip1_e,ij0_e,iky,ikx,izgs:izge,updatetlevel)&
             *imagu*ky_*CONJG(psi(iky,ikx,izgs:izge))
      ENDDO
    ENDDO
  ENDIF
  ! Integrate over z
  integrant(izgs:izge) = Jacobian(izgs:izge,0)*integrant(izgs:izge)
  call simpson_rule_z(integrant,integral)
  ! Get process local gyrocenter flux
  gflux_local = integral*iInt_Jacobian
  !
  integrant   = 0._dp ! reset auxiliary variable
  !!---------- Particle flux (gyrokinetic) ------------
  ! Electrostatic part
  IF(CONTAINS_ip0_e) THEN
    DO iky = ikys,ikye
        ky_ = kyarray(iky)
        DO ikx = ikxs,ikxe
          DO ij = ijs_e, ije_e
            integrant(izgs:izge) = integrant(izgs:izge)&
              +moments_e(ip0_e,ij,iky,ikx,izgs:izge,updatetlevel)&
              *imagu*ky_*kernel_e(ij,iky,ikx,izgs:izge,0)*CONJG(phi(iky,ikx,izgs:izge))
          ENDDO
      ENDDO
    ENDDO
  ENDIF
  ! Electromagnetic part
  IF( EM .AND. CONTAINS_ip1_e ) THEN
    DO iky = ikys,ikye
      ky_ = kyarray(iky)
      DO ikx = ikxs,ikxe
        integrant   = 0._dp ! auxiliary variable for z integration
        DO ij = ijs_e, ije_e
          integrant(izgs:izge) = integrant(izgs:izge)&
            -sqrt_tau_o_sigma_e*moments_e(ip1_e,ij,iky,ikx,izgs:izge,updatetlevel)&
            *imagu*ky_*kernel_e(ij,iky,ikx,izgs:izge,0)*CONJG(psi(iky,ikx,izgs:izge))
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  ! Integrate over z
  integrant(izgs:izge) = Jacobian(izgs:izge,0)*integrant(izgs:izge)
  call simpson_rule_z(integrant,integral)
  ! Get process local particle flux
  pflux_local = integral*iInt_Jacobian

  !!!!---------- Sum over all processes ----------
  buffer(1) = 2._dp*REAL(gflux_local)
  buffer(2) = 2._dp*REAL(pflux_local)
  root = 0
  !Gather manually among the rank_p=0 processes and perform the sum
  gflux_re = 0
  pflux_re = 0
  IF (num_procs_ky .GT. 1) THEN
      !! Everyone sends its local_sum to root = 0
      IF (rank_ky .NE. root) THEN
          CALL MPI_SEND(buffer, 2 , MPI_DOUBLE_PRECISION, root, 1234, comm_ky, ierr)
      ELSE
          ! Recieve from all the other processes
          DO i_ = 0,num_procs_ky-1
              IF (i_ .NE. rank_ky) &
                  CALL MPI_RECV(buffer, 2 , MPI_DOUBLE_PRECISION, i_, 1234, comm_ky, MPI_STATUS_IGNORE, ierr)
                  gflux_re = gflux_re + buffer(1)
                  pflux_re = pflux_re + buffer(2)
          ENDDO
      ENDIF
  ELSE
    gflux_re = gflux_local
    pflux_re = pflux_local
  ENDIF
END SUBROUTINE compute_radial_electron_transport

! 1D diagnostic to compute the average radial particle transport <T_i v_ExB_x>_xyz
SUBROUTINE compute_radial_ion_heatflux
  USE fields,           ONLY : moments_i, phi, psi
  USE array,            ONLY : kernel_i, HF_phi_correction_operator
  USE geometry,         ONLY : Jacobian, iInt_Jacobian
  USE time_integration, ONLY : updatetlevel
  USE model,            ONLY : qi_taui, KIN_E, tau_i, EM
  USE calculus,         ONLY : simpson_rule_z
  USE model,            ONLY : sqrt_tau_o_sigma_i
  IMPLICIT NONE
  COMPLEX(dp) :: hflux_local, integral
  REAL(dp)    :: ky_, buffer(1:2), n_dp
  INTEGER     :: i_, world_rank, world_size, root, in
  COMPLEX(dp), DIMENSION(izgs:izge) :: integrant        ! charge density q_a n_a

  hflux_local = 0._dp ! heat flux
  integrant   = 0._dp ! z integration auxiliary variable
  !!----------------ELECTROSTATIC CONTRIBUTION---------------------------
  IF(CONTAINS_ip0_i .AND. CONTAINS_ip2_i) THEN
    ! Loop to compute gamma_kx = sum_ky sum_j -i k_z Kernel_j Ni00 * phi
    DO iky = ikys,ikye
    ky_ = kyarray(iky)
    DO ikx = ikxs,ikxe
      DO in = ijs_i, ije_i
        n_dp = jarray_i(in)
        integrant(izgs:izge) = integrant(izgs:izge) + tau_i*imagu*ky_*CONJG(phi(iky,ikx,izgs:izge))&
         *kernel_i(in,iky,ikx,izgs:izge,0)*(&
                       0.5_dp*SQRT2*moments_i(ip2_i,in  ,iky,ikx,izgs:izge,updatetlevel)&
             +(2._dp*n_dp + 1.5_dp)*moments_i(ip0_i,in  ,iky,ikx,izgs:izge,updatetlevel)&
                      -(n_dp+1._dp)*moments_i(ip0_i,in+1,iky,ikx,izgs:izge,updatetlevel)&
                              -n_dp*moments_i(ip0_i,in-1,iky,ikx,izgs:izge,updatetlevel))
      ENDDO
    ENDDO
    ENDDO
  ENDIF
  IF(EM .AND. CONTAINS_ip1_i .AND. CONTAINS_ip3_i) THEN
    !!----------------ELECTROMAGNETIC CONTRIBUTION--------------------
    DO iky = ikys,ikye
    ky_ = kyarray(iky)
    DO ikx = ikxs,ikxe
      DO in = ijs_i, ije_i
        n_dp = jarray_i(in)
        integrant(izgs:izge) = integrant(izgs:izge) &
         +tau_i*sqrt_tau_o_sigma_i*imagu*ky_*CONJG(psi(iky,ikx,izgs:izge))&
           *kernel_i(in,iky,ikx,izgs:izge,0)*(&
                   0.5_dp*SQRT2*SQRT3*moments_i(ip3_i,in  ,iky,ikx,izgs:izge,updatetlevel)&
                        +1.5_dp*CONJG(moments_i(ip1_i,in  ,iky,ikx,izgs:izge,updatetlevel))& !?????
                  +(2._dp*n_dp+1._dp)*moments_i(ip1_i,in  ,iky,ikx,izgs:izge,updatetlevel)&
                        -(n_dp+1._dp)*moments_i(ip1_i,in+1,iky,ikx,izgs:izge,updatetlevel)&
                                -n_dp*moments_i(ip1_i,in-1,iky,ikx,izgs:izge,updatetlevel))
      ENDDO
    ENDDO
    ENDDO
  ENDIF
  ! Add polarisation contribution
  integrant(izs:ize) = integrant(izs:ize) + tau_i*imagu*ky_&
  *CONJG(phi(iky,ikx,izs:ize))*phi(iky,ikx,izs:ize) * HF_phi_correction_operator(iky,ikx,izs:ize)
  ! Integrate over z
  integrant(izgs:izge) = Jacobian(izgs:izge,0)*integrant(izgs:izge)
  call simpson_rule_z(integrant,integral)
  hflux_local = hflux_local + integral*iInt_Jacobian
  ! Double it for taking into account the other half plane
  buffer(2) = 2._dp*REAL(hflux_local)
  root = 0
  !Gather manually among the rank_p=0 processes and perform the sum
  hflux_xi = 0
  IF (num_procs_ky .GT. 1) THEN
      !! Everyone sends its local_sum to root = 0
      IF (rank_ky .NE. root) THEN
          CALL MPI_SEND(buffer, 2 , MPI_DOUBLE_PRECISION, root, 1234, comm_ky, ierr)
      ELSE
          ! Recieve from all the other processes
          DO i_ = 0,num_procs_ky-1
              IF (i_ .NE. rank_ky) &
                  CALL MPI_RECV(buffer, 2 , MPI_DOUBLE_PRECISION, i_, 1234, comm_ky, MPI_STATUS_IGNORE, ierr)
                  hflux_xi = hflux_xi + buffer(2)
          ENDDO
      ENDIF
  ELSE
    hflux_xi = hflux_local
  ENDIF
END SUBROUTINE compute_radial_ion_heatflux


! 1D diagnostic to compute the average radial particle transport <T_e v_ExB_x>_xyz
SUBROUTINE compute_radial_electron_heatflux
  USE fields,           ONLY : moments_e, phi, psi
  USE array,            ONLY : kernel_e, HF_phi_correction_operator
  USE geometry,         ONLY : Jacobian, iInt_Jacobian
  USE time_integration, ONLY : updatetlevel
  USE model,            ONLY : qi_taui, KIN_E, tau_e, EM
  USE calculus,         ONLY : simpson_rule_z
  USE model,            ONLY : sqrt_tau_o_sigma_e
  IMPLICIT NONE
  COMPLEX(dp) :: hflux_local, integral
  REAL(dp)    :: ky_, buffer(1:2), n_dp
  INTEGER     :: i_, world_rank, world_size, root, in
  COMPLEX(dp), DIMENSION(izgs:izge) :: integrant        ! charge density q_a n_a

  hflux_local = 0._dp ! heat flux
  integrant   = 0._dp ! z integration auxiliary variable
  !!----------------ELECTROSTATIC CONTRIBUTION---------------------------
  IF(CONTAINS_ip0_e .AND. CONTAINS_ip2_e) THEN
    ! Loop to compute gamma_kx = sum_ky sum_j -i k_z Kernel_j Ni00 * phi
    DO iky = ikys,ikye
    ky_ = kyarray(iky)
    DO ikx = ikxs,ikxe
      DO in = ijs_e, ije_e
        n_dp = jarray_e(in)
        integrant(izgs:izge) = integrant(izgs:izge) + tau_e*imagu*ky_*CONJG(phi(iky,ikx,izgs:izge))&
         *kernel_e(in,iky,ikx,izgs:izge,0)*(&
                       0.5_dp*SQRT2*moments_e(ip2_e,in  ,iky,ikx,izgs:izge,updatetlevel)&
             +(2._dp*n_dp + 1.5_dp)*moments_e(ip0_e,in  ,iky,ikx,izgs:izge,updatetlevel)&
                      -(n_dp+1._dp)*moments_e(ip0_e,in+1,iky,ikx,izgs:izge,updatetlevel)&
                              -n_dp*moments_e(ip0_e,in-1,iky,ikx,izgs:izge,updatetlevel))
      ENDDO
    ENDDO
    ENDDO
  ENDIF
  IF(EM .AND. CONTAINS_ip1_e .AND. CONTAINS_ip3_e) THEN
    !!----------------ELECTROMAGNETIC CONTRIBUTION--------------------
    DO iky = ikys,ikye
    ky_ = kyarray(iky)
    DO ikx = ikxs,ikxe
      DO in = ijs_e, ije_e
        n_dp = jarray_e(in)
        integrant(izgs:izge) = integrant(izgs:izge) &
         +tau_e*sqrt_tau_o_sigma_e*imagu*ky_*CONJG(psi(iky,ikx,izgs:izge))&
           *kernel_e(in,iky,ikx,izgs:izge,0)*(&
                   0.5_dp*SQRT2*SQRT3*moments_e(ip3_e,in  ,iky,ikx,izgs:izge,updatetlevel)&
                        +1.5_dp*CONJG(moments_e(ip1_e,in  ,iky,ikx,izgs:izge,updatetlevel))& !?????
                  +(2._dp*n_dp+1._dp)*moments_e(ip1_e,in  ,iky,ikx,izgs:izge,updatetlevel)&
                        -(n_dp+1._dp)*moments_e(ip1_e,in+1,iky,ikx,izgs:izge,updatetlevel)&
                                -n_dp*moments_e(ip1_e,in-1,iky,ikx,izgs:izge,updatetlevel))
      ENDDO
    ENDDO
    ENDDO
  ENDIF
  ! Add polarisation contribution
  integrant(izs:ize) = integrant(izs:ize) + tau_e*imagu*ky_&
  *CONJG(phi(iky,ikx,izs:ize))*phi(iky,ikx,izs:ize) * HF_phi_correction_operator(iky,ikx,izs:ize)
  ! Integrate over z
  integrant(izgs:izge) = Jacobian(izgs:izge,0)*integrant(izgs:izge)
  call simpson_rule_z(integrant,integral)
  hflux_local = hflux_local + integral*iInt_Jacobian
  ! Double it for taking into account the other half plane
  buffer(2) = 2._dp*REAL(hflux_local)
  root = 0
  !Gather manually among the rank_p=0 processes and perform the sum
  hflux_xi = 0
  IF (num_procs_ky .GT. 1) THEN
      !! Everyone sends its local_sum to root = 0
      IF (rank_ky .NE. root) THEN
          CALL MPI_SEND(buffer, 2 , MPI_DOUBLE_PRECISION, root, 1234, comm_ky, ierr)
      ELSE
          ! Recieve from all the other processes
          DO i_ = 0,num_procs_ky-1
              IF (i_ .NE. rank_ky) &
                  CALL MPI_RECV(buffer, 2 , MPI_DOUBLE_PRECISION, i_, 1234, comm_ky, MPI_STATUS_IGNORE, ierr)
                  hflux_xi = hflux_xi + buffer(2)
          ENDDO
      ENDIF
  ELSE
    hflux_xi = hflux_local
  ENDIF
END SUBROUTINE compute_radial_electron_heatflux



SUBROUTINE compute_nadiab_moments_z_gradients_and_interp
  ! evaluate the non-adiabatique ion moments
  !
  ! n_{pi} = N^{pj} + kernel(j) /tau_i phi delta_p0
  !
  USE fields,           ONLY : moments_i, moments_e, phi, psi
  USE array,            ONLY : kernel_e, kernel_i, nadiab_moments_e, nadiab_moments_i, &
                               ddz_nepj, ddzND_nepj, interp_nepj,&
                               ddz_nipj, ddzND_nipj, interp_nipj
  USE time_integration, ONLY : updatetlevel
  USE model,            ONLY : qe_taue, qi_taui,q_o_sqrt_tau_sigma_e, q_o_sqrt_tau_sigma_i, &
                               KIN_E, CLOS, beta
  USE calculus,         ONLY : grad_z, grad_z2, grad_z4, interp_z
  IMPLICIT NONE
  INTEGER :: eo, p_int, j_int
  CALL cpu_time(t0_process)

  ! Electron non adiab moments

    IF(KIN_E) THEN
      DO ip=ipgs_e,ipge_e
        IF(parray_e(ip) .EQ. 0) THEN
          DO ij=ijgs_e,ijge_e
            nadiab_moments_e(ip,ij,ikys:ikye,ikxs:ikxe,izgs:izge) = moments_e(ip,ij,ikys:ikye,ikxs:ikxe,izgs:izge,updatetlevel) &
                                      + qe_taue*kernel_e(ij,ikys:ikye,ikxs:ikxe,izgs:izge,0)*phi(ikys:ikye,ikxs:ikxe,izgs:izge)
          ENDDO
        ELSEIF( (parray_e(ip) .EQ. 1) .AND. (beta .GT. 0) ) THEN
          DO ij=ijgs_e,ijge_e
            nadiab_moments_e(ip,ij,ikys:ikye,ikxs:ikxe,izgs:izge) = moments_e(ip,ij,ikys:ikye,ikxs:ikxe,izgs:izge,updatetlevel) &
                                      - q_o_sqrt_tau_sigma_e*kernel_e(ij,ikys:ikye,ikxs:ikxe,izgs:izge,0)*psi(ikys:ikye,ikxs:ikxe,izgs:izge)
          ENDDO
        ELSE
          DO ij=ijgs_e,ijge_e
            nadiab_moments_e(ip,ij,ikys:ikye,ikxs:ikxe,izgs:izge) = moments_e(ip,ij,ikys:ikye,ikxs:ikxe,izgs:izge,updatetlevel)
          ENDDO
        ENDIF
      ENDDO
    ENDIF
    ! Ions non adiab moments
    DO ip=ipgs_i,ipge_i
      IF(parray_i(ip) .EQ. 0) THEN
        DO ij=ijgs_i,ijge_i
          nadiab_moments_i(ip,ij,ikys:ikye,ikxs:ikxe,izgs:izge) = moments_i(ip,ij,ikys:ikye,ikxs:ikxe,izgs:izge,updatetlevel) &
                                    + qi_taui*kernel_i(ij,ikys:ikye,ikxs:ikxe,izgs:izge,0)*phi(ikys:ikye,ikxs:ikxe,izgs:izge)
        ENDDO
      ELSEIF( (parray_i(ip) .EQ. 1) .AND. (beta .GT. 0) ) THEN
        DO ij=ijgs_i,ijge_i
          nadiab_moments_i(ip,ij,ikys:ikye,ikxs:ikxe,izgs:izge) = moments_i(ip,ij,ikys:ikye,ikxs:ikxe,izgs:izge,updatetlevel) &
                                    - q_o_sqrt_tau_sigma_i*kernel_i(ij,ikys:ikye,ikxs:ikxe,izgs:izge,0)*psi(ikys:ikye,ikxs:ikxe,izgs:izge)
        ENDDO
      ELSE
        DO ij=ijgs_i,ijge_i
          nadiab_moments_i(ip,ij,ikys:ikye,ikxs:ikxe,izgs:izge) = moments_i(ip,ij,ikys:ikye,ikxs:ikxe,izgs:izge,updatetlevel)
        ENDDO
      ENDIF
    ENDDO

    !! Ensure to kill the moments too high if closue option is set to 1
    IF(CLOS .EQ. 1) THEN
      IF(KIN_E) THEN
        DO ip=ipgs_e,ipge_e
          p_int = parray_e(ip)
            DO ij=ijgs_e,ijge_e
              j_int = jarray_e(ij)
              IF(p_int+2*j_int .GT. dmaxe) &
              nadiab_moments_e(ip,ij,:,:,:) = 0._dp
            ENDDO
        ENDDO
      ENDIF
      DO ip=ipgs_i,ipge_i
        p_int = parray_i(ip)
          DO ij=ijgs_i,ijge_i
            j_int = jarray_i(ij)
            IF(p_int+2*j_int .GT. dmaxi) &
            nadiab_moments_i(ip,ij,:,:,:) = 0._dp
          ENDDO
      ENDDO
    ENDIF

 !------------- INTERP AND GRADIENTS ALONG Z ----------------------------------

  IF (KIN_E) THEN
  DO ikx = ikxs,ikxe
    DO iky = ikys,ikye
      DO ij = ijgs_e,ijge_e
        DO ip = ipgs_e,ipge_e
          p_int = parray_e(ip)
          eo    = MODULO(p_int,2) ! Indicates if we are on even or odd z grid
          ! Compute z derivatives
          CALL   grad_z(eo,nadiab_moments_e(ip,ij,iky,ikx,izgs:izge),    ddz_nepj(ip,ij,iky,ikx,izs:ize))
          ! Parallel hyperdiffusion
          CALL  grad_z4(moments_e(ip,ij,iky,ikx,izgs:izge,updatetlevel),ddzND_nepj(ip,ij,iky,ikx,izs:ize))
          ! CALL  grad_z2(nadiab_moments_e(ip,ij,iky,ikx,izgs:izge),ddzND_nepj(ip,ij,iky,ikx,izs:ize))
          ! Compute even odd grids interpolation
          CALL interp_z(eo,nadiab_moments_e(ip,ij,iky,ikx,izgs:izge), interp_nepj(ip,ij,iky,ikx,izs:ize))
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  ENDIF

  DO ikx = ikxs,ikxe
    DO iky = ikys,ikye
      DO ij = ijgs_i,ijge_i
        DO ip = ipgs_i,ipge_i
          p_int = parray_i(ip)
          eo    = MODULO(p_int,2) ! Indicates if we are on even or odd z grid
          ! Compute z first derivative
          CALL   grad_z(eo,nadiab_moments_i(ip,ij,iky,ikx,izgs:izge),    ddz_nipj(ip,ij,iky,ikx,izs:ize))
          ! Parallel numerical diffusion
          CALL  grad_z4(moments_i(ip,ij,iky,ikx,izgs:izge,updatetlevel),ddzND_nipj(ip,ij,iky,ikx,izs:ize))
          ! CALL  grad_z2(nadiab_moments_i(ip,ij,iky,ikx,izgs:izge),ddzND_nipj(ip,ij,iky,ikx,izs:ize))
          ! Compute even odd grids interpolation
          CALL interp_z(eo,nadiab_moments_i(ip,ij,iky,ikx,izgs:izge), interp_nipj(ip,ij,iky,ikx,izs:ize))
        ENDDO
      ENDDO
    ENDDO
  ENDDO

  ! Execution time end
  CALL cpu_time(t1_process)
  tc_process = tc_process + (t1_process - t0_process)
END SUBROUTINE compute_nadiab_moments_z_gradients_and_interp

SUBROUTINE compute_Napjz_spectrum
  USE fields, ONLY : moments_e, moments_i
  USE model,  ONLY : KIN_E
  USE array,  ONLY : Nipjz, Nepjz
  USE time_integration, ONLY : updatetlevel
  IMPLICIT NONE
  REAL(dp), DIMENSION(ips_e:ipe_e,ijs_e:ije_e,izs:ize) :: local_sum_e,global_sum_e, buffer_e
  REAL(dp), DIMENSION(ips_i:ipe_i,ijs_i:ije_i,izs:ize) :: local_sum_i,global_sum_i, buffer_i
  INTEGER  :: i_, world_rank, world_size, root, count
  root = 0
  ! Electron moments spectrum
  IF (KIN_E) THEN
    ! build local sum
    local_sum_e = 0._dp
    DO ikx = ikxs,ikxe
      DO iky = ikys,ikye
        local_sum_e(:,:,:)  = local_sum_e(:,:,:)  + &
        moments_e(:,:,iky,ikx,:,updatetlevel) * CONJG(moments_e(:,:,iky,ikx,:,updatetlevel))
      ENDDO
    ENDDO
    ! sum reduction
    buffer_e     = local_sum_e
    global_sum_e = 0._dp
    count = (ipe_e-ips_e+1)*(ije_e-ijs_e+1)*(ize-izs+1)
    IF (num_procs_ky .GT. 1) THEN
        !! Everyone sends its local_sum to root = 0
        IF (rank_ky .NE. root) THEN
            CALL MPI_SEND(buffer_e, count , MPI_DOUBLE_PRECISION, root, 1234, comm_ky, ierr)
        ELSE
            ! Recieve from all the other processes
            DO i_ = 0,num_procs_ky-1
                IF (i_ .NE. rank_ky) &
                    CALL MPI_RECV(buffer_e, count , MPI_DOUBLE_PRECISION, i_, 1234, comm_ky, MPI_STATUS_IGNORE, ierr)
                    global_sum_e = global_sum_e + buffer_e
            ENDDO
        ENDIF
    ELSE
      global_sum_e = local_sum_e
    ENDIF
    Nepjz = global_sum_e
  ENDIF
  ! Ion moment spectrum
  ! build local sum
  local_sum_i = 0._dp
  DO ikx = ikxs,ikxe
    DO iky = ikys,ikye
      local_sum_i(:,:,:)  = local_sum_i(:,:,:)  + &
      moments_i(:,:,iky,ikx,:,updatetlevel) * CONJG(moments_i(:,:,iky,ikx,:,updatetlevel))
    ENDDO
  ENDDO
  ! sum reduction
  buffer_i     = local_sum_i
  global_sum_i = 0._dp
  count = (ipe_i-ips_i+1)*(ije_i-ijs_i+1)*(ize-izs+1)
  IF (num_procs_ky .GT. 1) THEN
      !! Everyone sends its local_sum to root = 0
      IF (rank_ky .NE. root) THEN
          CALL MPI_SEND(buffer_i, count , MPI_DOUBLE_PRECISION, root, 1234, comm_ky, ierr)
      ELSE
          ! Recieve from all the other processes
          DO i_ = 0,num_procs_ky-1
              IF (i_ .NE. rank_ky) &
                  CALL MPI_RECV(buffer_i, count , MPI_DOUBLE_PRECISION, i_, 1234, comm_ky, MPI_STATUS_IGNORE, ierr)
                  global_sum_i = global_sum_i + buffer_i
          ENDDO
      ENDIF
  ELSE
    global_sum_i = local_sum_i
  ENDIF
  Nipjz = global_sum_i
END SUBROUTINE compute_Napjz_spectrum

!_____________________________________________________________________________!
!!!!! FLUID MOMENTS COMPUTATIONS !!!!!
! Compute the 2D particle density for electron and ions (sum over Laguerre)
SUBROUTINE compute_density
  USE array, ONLY : dens_e, dens_i, kernel_e, kernel_i, nadiab_moments_e, nadiab_moments_i
  USE model, ONLY : KIN_E
  IMPLICIT NONE
  COMPLEX(dp) :: dens

  IF ( CONTAINS_ip0_e .AND. CONTAINS_ip0_i ) THEN
      ! Loop to compute dens_i = sum_j kernel_j Ni0j at each k
      DO iz = izs,ize
        DO iky = ikys,ikye
          DO ikx = ikxs,ikxe

            IF(KIN_E) THEN
            ! electron density
            dens = 0._dp
            DO ij = ijs_e, ije_e
                dens = dens + kernel_e(ij,iky,ikx,iz,0) * nadiab_moments_e(ip0_e,ij,iky,ikx,iz)
            ENDDO
            dens_e(iky,ikx,iz) = dens
            ENDIF
            ! ion density
            dens = 0._dp
            DO ij = ijs_i, ije_i
                dens = dens + kernel_i(ij,iky,ikx,iz,0) * nadiab_moments_i(ip0_e,ij,iky,ikx,iz)
            ENDDO
            dens_i(iky,ikx,iz) = dens
          ENDDO
        ENDDO
      ENDDO
  ENDIF
END SUBROUTINE compute_density

! Compute the 2D particle fluid perp velocity for electron and ions (sum over Laguerre)
SUBROUTINE compute_uperp
  USE array, ONLY : uper_e, uper_i, kernel_e, kernel_i, nadiab_moments_e, nadiab_moments_i
  USE model, ONLY : KIN_E
  IMPLICIT NONE
  COMPLEX(dp) :: uperp

  IF ( CONTAINS_ip0_e .AND. CONTAINS_ip0_i ) THEN
      DO iz = izs,ize
        DO iky = ikys,ikye
          DO ikx = ikxs,ikxe

            IF(KIN_E) THEN
            ! electron
            uperp = 0._dp
            DO ij = ijs_e, ije_e
                uperp = uperp + kernel_e(ij,iky,ikx,iz,0) *&
                 0.5_dp*(nadiab_moments_e(ip0_e,ij,iky,ikx,iz) - nadiab_moments_e(ip0_e,ij-1,iky,ikx,iz))
            ENDDO
            uper_e(iky,ikx,iz) = uperp
            ENDIF
            ! ion
            uperp = 0._dp
            DO ij = ijs_i, ije_i
              uperp = uperp + kernel_i(ij,iky,ikx,iz,0) *&
               0.5_dp*(nadiab_moments_i(ip0_i,ij,iky,ikx,iz) - nadiab_moments_i(ip0_i,ij-1,iky,ikx,iz))
             ENDDO
            uper_i(iky,ikx,iz) = uperp
          ENDDO
        ENDDO
      ENDDO
  ENDIF
END SUBROUTINE compute_uperp

! Compute the 2D particle fluid par velocity for electron and ions (sum over Laguerre)
SUBROUTINE compute_upar
  USE array, ONLY : upar_e, upar_i, kernel_e, kernel_i, nadiab_moments_e, nadiab_moments_i
  USE model, ONLY : KIN_E
  IMPLICIT NONE
  COMPLEX(dp) :: upar

  IF ( CONTAINS_ip1_e .AND. CONTAINS_ip1_i ) THEN
    DO iz = izs,ize
      DO iky = ikys,ikye
        DO ikx = ikxs,ikxe
          IF(KIN_E) THEN
          ! electron
          upar = 0._dp
          DO ij = ijs_e, ije_e
            upar = upar + kernel_e(ij,iky,ikx,iz,1)*nadiab_moments_e(ip1_e,ij,iky,ikx,iz)
          ENDDO
          upar_e(iky,ikx,iz) = upar
          ENDIF
          ! ion
          upar = 0._dp
          DO ij = ijs_i, ije_i
            upar = upar + kernel_i(ij,iky,ikx,iz,1)*nadiab_moments_i(ip1_i,ij,iky,ikx,iz)
           ENDDO
          upar_i(iky,ikx,iz) = upar
        ENDDO
      ENDDO
    ENDDO
  ELSE
    IF(KIN_E)&
    upar_e = 0
    upar_i = 0
  ENDIF
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
      DO iz = izs,ize
        DO iky = ikys,ikye
          DO ikx = ikxs,ikxe
            ! electron temperature
            IF(KIN_E) THEN
            Tperp  = 0._dp
            DO ij = ijs_e, ije_e
              j_dp = REAL(ij-1,dp)
              Tperp = Tperp + kernel_e(ij,iky,ikx,iz,0)*&
                  ((2_dp*j_dp+1)*nadiab_moments_e(ip0_e,ij  ,iky,ikx,iz)&
                  -j_dp         *nadiab_moments_e(ip0_e,ij-1,iky,ikx,iz)&
                  -j_dp+1       *nadiab_moments_e(ip0_e,ij+1,iky,ikx,iz))
            ENDDO
            Tper_e(iky,ikx,iz) = Tperp
            ENDIF
            ! ion temperature
            Tperp = 0._dp
            DO ij = ijs_i, ije_i
              j_dp = REAL(ij-1,dp)
              Tperp = Tperp + kernel_i(ij,iky,ikx,iz,0)*&
                  ((2_dp*j_dp+1)*nadiab_moments_i(ip0_i,ij  ,iky,ikx,iz)&
                  -j_dp         *nadiab_moments_i(ip0_i,ij-1,iky,ikx,iz)&
                  -j_dp+1       *nadiab_moments_i(ip0_i,ij+1,iky,ikx,iz))
            ENDDO
            Tper_i(iky,ikx,iz) = Tperp
          ENDDO
        ENDDO
      ENDDO
  ENDIF
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
      DO iz = izs,ize
        DO iky = ikys,ikye
          DO ikx = ikxs,ikxe
            ! electron temperature
            IF(KIN_E) THEN
            Tpar  = 0._dp
            DO ij = ijs_e, ije_e
              j_dp = REAL(ij-1,dp)
              Tpar  = Tpar + kernel_e(ij,iky,ikx,iz,0)*&
               (SQRT2 * nadiab_moments_e(ip2_e,ij,iky,ikx,iz) &
                      + nadiab_moments_e(ip0_e,ij,iky,ikx,iz))
            ENDDO
            Tpar_e(iky,ikx,iz) = Tpar
            ENDIF
            ! ion temperature
            Tpar  = 0._dp
            DO ij = ijs_i, ije_i
              j_dp = REAL(ij-1,dp)
              Tpar  = Tpar + kernel_i(ij,iky,ikx,iz,0)*&
               (SQRT2 * nadiab_moments_i(ip2_i,ij,iky,ikx,iz) &
                      + nadiab_moments_i(ip0_i,ij,iky,ikx,iz))
            ENDDO
            Tpar_i(iky,ikx,iz) = Tpar
          ENDDO
        ENDDO
      ENDDO
  ENDIF
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

END MODULE processing
