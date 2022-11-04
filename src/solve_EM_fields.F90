SUBROUTINE solve_EM_fields
  ! Solve Poisson and Ampere equations
  USE model, ONLY : beta
  USE basic
  USE prec_const
  IMPLICIT NONE

  CALL poisson

  IF(beta .GT. 0._dp) &
  CALL ampere

CONTAINS

  SUBROUTINE poisson
    ! Solve poisson equation to get phi
    USE time_integration, ONLY: updatetlevel
    USE array,            ONLY: kernel_e, kernel_i, inv_poisson_op, inv_pol_ion
    USE fields,           ONLY: phi, moments_e, moments_i
    USE grid
    USE calculus,         ONLY : simpson_rule_z
    USE parallel,         ONLY : manual_3D_bcast
    use model,            ONLY : q_e, q_i, lambdaD, KIN_E, sigma_e, sigma_i
    USE processing,       ONLY : compute_density
    USE geometry,         ONLY : iInt_Jacobian, Jacobian
    IMPLICIT NONE

    INTEGER     :: ini,ine
    COMPLEX(dp) :: fsa_phi, intf_   ! current flux averaged phi
    COMPLEX(dp), DIMENSION(izs:ize) :: rho_i, rho_e, integrant  ! charge density q_a n_a and aux var

    ! Execution time start
    CALL cpu_time(t0_poisson)
    !! Poisson can be solved only for process containing p=0
    IF ( SOLVE_POISSON ) THEN
      kxloop: DO ikx = ikxs,ikxe
        kyloop: DO iky = ikys,ikye
          phi(iky,ikx,izs:ize)  = 0._dp
          !!!! Compute ion particle charge density q_i n_i
          rho_i = 0._dp
          DO ini=1,jmaxi+1
            rho_i(izs:ize) = rho_i(izs:ize) &
             +q_i*kernel_i(ini,iky,ikx,izs:ize,0)*moments_i(ip0_i,ini,iky,ikx,izs:ize,updatetlevel)
          END DO
          !!!! Compute electron particle charge density q_e n_e
          rho_e = 0._dp
          IF (KIN_E) THEN ! Kinetic electrons
          DO ine=1,jmaxe+1
            rho_e(izs:ize) = rho_e(izs:ize) &
             +q_e*kernel_e(ine,iky,ikx,izs:ize,0)*moments_e(ip0_e,ine,iky,ikx,izs:ize,updatetlevel)
          END DO
          ELSE  ! Adiabatic electrons
            ! Adiabatic charge density (linked to flux surface averaged phi)
            ! We compute the flux surface average solving a flux surface averaged
            ! Poisson equation, i.e.
            ! [qi^2(1-sum_j K_j^2)/tau_i] <phi>_psi = <q_i n_i >_psi
            !       inv_pol_ion^-1         fsa_phi  = simpson(Jacobian rho_i ) * iInt_Jacobian
            fsa_phi = 0._dp
            IF(kyarray(iky).EQ.0._dp) THEN ! take ky=0 mode (y-average)
              ! Prepare integrant for z-average
              integrant(izs:ize) = Jacobian(izs:ize,0)*rho_i(izs:ize)*inv_pol_ion(iky,ikx,izs:ize)
              call simpson_rule_z(integrant(izs:ize),intf_) ! get the flux averaged phi
              fsa_phi = intf_ * iInt_Jacobian !Normalize by 1/int(Jxyz)dz
            ENDIF
            rho_e(izs:ize) = fsa_phi
          ENDIF
          !!!!!!!!!!!!!!! Inverting the poisson equation !!!!!!!!!!!!!!!!!!!!!!!!!!
          phi(iky,ikx,izs:ize) = (rho_e(izs:ize) + rho_i(izs:ize))*inv_poisson_op(iky,ikx,izs:ize)
        END DO kyloop
      END DO kxloop
      ! Cancel origin singularity
      IF (contains_kx0 .AND. contains_ky0) phi(iky_0,ikx_0,:) = 0._dp
    ENDIF

    ! Transfer phi to all the others process along p
    CALL manual_3D_bcast(phi(ikys:ikye,ikxs:ikxe,izs:ize))

    ! Execution time end
    CALL cpu_time(t1_poisson)
    tc_poisson = tc_poisson + (t1_poisson - t0_poisson)
  END SUBROUTINE poisson

  SUBROUTINE ampere
    ! Solve ampere equation to get psi
    USE time_integration, ONLY: updatetlevel
    USE array,            ONLY: kernel_e, kernel_i, inv_ampere_op
    USE fields,           ONLY: moments_i, moments_e, psi
    USE grid
    USE parallel,         ONLY : manual_3D_bcast
    use model,            ONLY : sqrt_tau_o_sigma_e, sqrt_tau_o_sigma_i, q_e, q_i, beta, KIN_E
    IMPLICIT NONE

    INTEGER     :: ini,ine

    ! Execution time start
    CALL cpu_time(t0_poisson)
    !! Ampere can be solved only with beta > 0 and for process containing p=1 moments
    IF ( SOLVE_AMPERE ) THEN
      psi(ikys:ikye,ikxs:ikxe,izs:ize)  = 0._dp
      !!!! ion particle current density contribution "q_i u_i"
      DO ini=1,jmaxi+1
        psi(ikys:ikye,ikxs:ikxe,izs:ize) = psi(ikys:ikye,ikxs:ikxe,izs:ize) &
         +q_i*sqrt_tau_o_sigma_i*kernel_i(ini,ikys:ikye,ikxs:ikxe,izs:ize,0)*moments_i(ip1_i,ini,ikys:ikye,ikxs:ikxe,izs:ize,updatetlevel)
      END DO
      !!!! electron particle current density contribution "q_e u_e"
      DO ine=1,jmaxe+1
        psi(ikys:ikye,ikxs:ikxe,izs:ize) = psi(ikys:ikye,ikxs:ikxe,izs:ize) &
          +q_e*sqrt_tau_o_sigma_e*kernel_e(ine,ikys:ikye,ikxs:ikxe,izs:ize,0)*moments_e(ip1_e,ine,ikys:ikye,ikxs:ikxe,izs:ize,updatetlevel)
      END DO
      !!!!!!!!!!!!!!! Inverting the poisson equation !!!!!!!!!!!!!!!!!!!!!!!!!!
      psi(ikys:ikye,ikxs:ikxe,izs:ize) = beta*psi(ikys:ikye,ikxs:ikxe,izs:ize)*inv_ampere_op(ikys:ikye,ikxs:ikxe,izs:ize)

      ! Cancel origin singularity
      IF (contains_kx0 .AND. contains_ky0) psi(iky_0,ikx_0,:) = 0._dp
    ENDIF

    ! Transfer phi to all the others process along p
    CALL manual_3D_bcast(psi(ikys:ikye,ikxs:ikxe,izs:ize))

    ! Execution time end
    CALL cpu_time(t1_poisson)
    tc_poisson = tc_poisson + (t1_poisson - t0_poisson)
  END SUBROUTINE ampere

END SUBROUTINE solve_EM_fields
