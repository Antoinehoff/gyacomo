SUBROUTINE solve_EM_fields
  ! Solve Poisson and Ampere equations
  USE model, ONLY : beta
  USE basic
  USE prec_const
  IMPLICIT NONE

  CALL poisson

  IF(beta .GT. 0._xp) &
  CALL ampere

CONTAINS

  SUBROUTINE poisson
    ! Solve poisson equation to get phi
    USE time_integration, ONLY: updatetlevel
    USE array,            ONLY: kernel, inv_poisson_op, inv_pol_ion
    USE fields,           ONLY: phi, moments
    USE grid,             ONLY: local_na, local_nky, local_nkx, local_nz,ngz, SOLVE_POISSON,&
                                kyarray, contains_kx0, contains_ky0,ikx0,iky0, zweights_SR, ieven,&
                                ip0, total_nj, ngj
    USE calculus,         ONLY: simpson_rule_z
    USE parallel,         ONLY: manual_3D_bcast
    USE model,            ONLY: lambdaD, ADIAB_E
    use species,          ONLY: q
    USE processing,       ONLY: compute_density
    USE geometry,         ONLY: iInt_Jacobian, Jacobian
    IMPLICIT NONE

    INTEGER     :: in, ia, ikx, iky, iz, izi, ini
    COMPLEX(xp) :: fsa_phi, intf_, rhtot   ! current flux averaged phi
    COMPLEX(xp), DIMENSION(local_nz) :: rho, integrant  ! charge density q_a n_a and aux var
    rhtot = 0
    !! Poisson can be solved only for process containng p=0
    IF ( SOLVE_POISSON ) THEN
        x:DO ikx = 1,local_nkx
          y:DO iky = 1,local_nky
            !!!!!!!!!!!!!!! Compute particle charge density q_a n_a for each evolved species
            DO iz = 1,local_nz
              izi = iz+ngz/2
              rho(iz) = 0._xp
              DO in = 1,total_nj
                ini = in+ngj/2
                DO ia = 1,local_na
                  rho(iz) = rho(iz) + q(ia)*kernel(ia,ini,iky,ikx,izi,ieven)&
                              *moments(ia,ip0,ini,iky,ikx,izi,updatetlevel)
                END DO
              END DO
            END DO
            !!!!!!!!!!!!!!! adiabatic electron contribution if asked
            ! Adiabatic charge density (linked to flux surface averaged phi)
            ! We compute the flux surface average solving a flux surface averaged
            ! Poisson equation, i.e.
            ! [qi^2(1-sum_j K_j^2)/tau_i] <phi>_psi = <q_i n_i >_psi
            !       inv_pol_ion^-1         fsa_phi  = simpson(Jacobian rho_i ) * iInt_Jacobian
            IF (ADIAB_E) THEN
              fsa_phi = 0._xp
              IF(kyarray(iky).EQ.0._xp) THEN ! take ky=0 mode (y-average)
                ! Prepare integrant for z-average
                DO iz = 1,local_nz
                  izi = iz+ngz/2
                  integrant(iz) = Jacobian(izi,ieven)*rho(iz)*inv_pol_ion(iky,ikx,iz)
                ENDDO
                call simpson_rule_z(local_nz,zweights_SR,integrant,intf_) ! get the flux averaged phi
                fsa_phi = intf_ * iInt_Jacobian !Normalize by 1/int(Jxyz)dz
              ENDIF
              rho = rho + fsa_phi
            ENDIF
            !!!!!!!!!!!!!!! adiabatic ions ?
            ! IF (ADIAB_I) THEN
            ! ENDIF
            !!!!!!!!!!!!!!! Inverting the poisson equation
            DO iz = 1,local_nz
              phi(iky,ikx,iz+ngz/2) = inv_poisson_op(iky,ikx,iz)*rho(iz)
            ENDDO
            rhtot = rhtot + sum(real(rho))
          ENDDO y
        ENDDO x
      ! Cancel origin singularity
      IF (contains_kx0 .AND. contains_ky0) phi(iky0,ikx0,:) = 0._xp
    ENDIF
    ! Transfer phi to all the others process along p
    CALL manual_3D_bcast(phi,local_nky,local_nkx,local_nz+ngz)
  END SUBROUTINE poisson

  SUBROUTINE ampere
    ! Solve ampere equation to get psi
    USE time_integration, ONLY: updatetlevel
    USE array,            ONLY: kernel, inv_ampere_op
    USE fields,           ONLY: moments, psi
    USE grid,             ONLY: local_na, local_nky, local_nkx, local_nz,ngz, SOLVE_AMPERE,&
                                contains_kx0, contains_ky0,ikx0,iky0, iodd,&
                                ip1, total_nj, ngj
    USE parallel,         ONLY: manual_3D_bcast
    use species,          ONLY: sqrt_tau_o_sigma, q
    use model,            ONLY: beta
    IMPLICIT NONE
    COMPLEX(xp) :: j_a ! current density
    INTEGER     :: in, ia, ikx, iky, iz, ini, izi
    !! Ampere can be solved only with beta > 0 and for process containng p=1 moments
    IF ( SOLVE_AMPERE ) THEN
      z:DO iz = 1,local_nz
      izi = iz+ngz/2
        x:DO ikx = 1,local_nkx
          y:DO iky = 1,local_nky
          !!!!!!!!!!!!!!! compute current density contribution "j_a = q_a u_a" for each species
          j_a = 0._xp
          n:DO in=1,total_nj
          ini = in+ngj/2
            a:DO ia = 1,local_na
            j_a = j_a &
              +q(ia)*sqrt_tau_o_sigma(ia)*kernel(ia,ini,iky,ikx,izi,iodd)*moments(ia,ip1,ini,iky,ikx,izi,updatetlevel)
            ENDDO a
          ENDDO n
          !!!!!!!!!!!!!!! Inverting the Ampere equation
          psi(iky,ikx,iz+ngz/2) = beta*inv_ampere_op(iky,ikx,iz)*j_a
          ENDDO y
        ENDDO x
      ENDDO z
    ENDIF
    ! Cancel origin singularity
    IF (contains_kx0 .AND. contains_ky0) psi(iky0,ikx0,:) = 0._xp
    ! Transfer phi to all the others process along p
    CALL manual_3D_bcast(psi,local_nky,local_nkx,local_nz+ngz)
  END SUBROUTINE ampere

END SUBROUTINE solve_EM_fields
