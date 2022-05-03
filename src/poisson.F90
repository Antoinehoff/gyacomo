SUBROUTINE poisson
  ! Solve poisson equation to get phi

  USE basic
  USE time_integration, ONLY: updatetlevel
  USE array
  USE fields
  USE grid
  USE calculus, ONLY : simpson_rule_z
  USE utility,  ONLY : manual_3D_bcast
  use model, ONLY : qe2_taue, qi2_taui, q_e, q_i, lambdaD, KIN_E
  USE processing, ONLY : compute_density
  USE prec_const
  USE geometry, ONLY : iInt_Jacobian, Jacobian
  IMPLICIT NONE

  INTEGER     :: ini,ine, i_, root_bcast
  COMPLEX(dp) :: fa_phi, intf_   ! current flux averaged phi
  INTEGER     :: count !! mpi integer to broadcast the electric potential at the end
  COMPLEX(dp), DIMENSION(izs:ize) :: rho_i, rho_e     ! charge density q_a n_a

  ! Execution time start
  CALL cpu_time(t0_poisson)
  !! Poisson can be solved only for process containing p=0
  IF ( (ips_e .EQ. ip0_e) .AND. (ips_i .EQ. ip0_i) ) THEN


    kxloop: DO ikx = ikxs,ikxe
      kyloop: DO iky = ikys,ikye
        !!!! Compute ion gyro charge density
        rho_i = 0._dp
        DO ini=1,jmaxi+1
          DO iz = izs,ize
          rho_i(iz) = rho_i(iz) &
           +q_i*kernel_i(ini,iky,ikx,iz,0)*moments_i(ip0_i,ini,iky,ikx,iz,updatetlevel)
          ENDDO
        END DO

        !!!! Compute electron gyro charge density
        rho_e = 0._dp
        IF (KIN_E) THEN ! Kinetic electrons
          DO ine=1,jmaxe+1
            DO iz = izs,ize
            rho_e(iz) = rho_e(iz) &
             +q_e*kernel_e(ine,iky,ikx,iz,0)*moments_e(ip0_e,ine,iky,ikx,iz,updatetlevel)
            ENDDO
          END DO
        ELSE  ! Adiabatic electrons
          ! Adiabatic charge density (linked to flux averaged phi)
          fa_phi = 0._dp
          IF(kyarray(iky).EQ.0._dp) THEN ! take ky=0 mode (average)
            DO ini=1,jmaxi+1 ! some over FLR contributions
                rho_e(izs:ize) = Jacobian(izs:ize,0)*moments_i(ip0_i,ini,iky,ikx,izs:ize,updatetlevel)&
                           *kernel_i(ini,iky,ikx,izs:ize,0)*(inv_poisson_op(iky,ikx,izs:ize)-1._dp)
                call simpson_rule_z(rho_e(izs:ize),intf_) ! integrate over z
                fa_phi = fa_phi + intf_
            ENDDO
            rho_e(izs:ize) = fa_phi*iInt_Jacobian !Normalize by 1/int(Jxyz)dz
          ENDIF
        ENDIF

        !!!!!!!!!!!!!!! Inverting the poisson equation !!!!!!!!!!!!!!!!!!!!!!!!!!
        DO iz = izs,ize
        phi(iky, ikx, iz) =  (rho_e(iz) + rho_i(iz))*inv_poisson_op(iky,ikx,iz)
        ENDDO
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
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE poisson
