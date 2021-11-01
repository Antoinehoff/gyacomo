SUBROUTINE poisson
  ! Solve poisson equation to get phi

  USE basic
  USE time_integration, ONLY: updatetlevel
  USE array
  USE fields
  USE grid
  USE utility
  use model, ONLY : qe2_taue, qi2_taui, q_e, q_i, lambdaD, KIN_E
  USE processing, ONLY : compute_density
  USE prec_const
  USE geometry, ONLY : iInt_Jacobian, Jacobian
  IMPLICIT NONE

  INTEGER     :: ini,ine, i_, root_bcast
  COMPLEX(dp) :: fa_phi, intf_   ! current flux averaged phi
  INTEGER     :: count !! mpi integer to broadcast the electric potential at the end
  COMPLEX(dp) :: buffer(ikxs:ikxe,ikys:ikye)
  COMPLEX(dp), DIMENSION(izs:ize) :: rho_i, rho_e     ! charge density q_a n_a
  COMPLEX(dp), DIMENSION(izs:ize) :: integrant        ! charge density q_a n_a

  !! Poisson can be solved only for process containing ip=1
  IF ( (ips_e .EQ. ip0_e) .AND. (ips_i .EQ. ip0_i) ) THEN
    ! Execution time start
    CALL cpu_time(t0_poisson)

    kxloop: DO ikx = ikxs,ikxe
      kyloop: DO iky = ikys,ikye
        !!!! Compute ion gyro charge density
        rho_i = 0._dp
        DO ini=1,jmaxi+1
          rho_i = rho_i &
           +q_i*kernel_i(ini,ikx,iky,:)*moments_i(ip0_i,ini,ikx,iky,:,updatetlevel)
        END DO

        !!!! Compute electron gyro charge density
        rho_e = 0._dp
        IF (KIN_E) THEN ! Kinetic electrons
          DO ine=1,jmaxe+1
            rho_e = rho_e &
             +q_e*kernel_e(ine,ikx,iky,:)*moments_e(ip0_e,ine, ikx,iky,:,updatetlevel)
          END DO
        ELSE  ! Adiabatic electrons
          ! Adiabatic charge density (linked to flux averaged phi)
          fa_phi = 0._dp
          IF(kyarray(iky).EQ.0._dp) THEN
            DO ini=1,jmaxi+1
                rho_e(:) = Jacobian(:)*moments_i(ip0_i,ini,ikx,iky,:,updatetlevel)&
                           *kernel_i(ini,ikx,iky,:)*(inv_poisson_op(ikx,iky,:)-1._dp)
                call simpson_rule_z(rho_e,intf_)
                fa_phi = fa_phi + intf_
            ENDDO
            rho_e = fa_phi*iInt_Jacobian !Normalize by 1/int(Jxyz)dz
          ENDIF
        ENDIF

        !!!!!!!!!!!!!!! Inverting the poisson equation !!!!!!!!!!!!!!!!!!!!!!!!!!
        phi(ikx, iky, :) =  (rho_e + rho_i)*inv_poisson_op(ikx,iky,:)
      END DO kyloop
    END DO kxloop

    ! Cancel origin singularity
    IF (contains_kx0 .AND. contains_ky0) phi(ikx_0,iky_0,:) = 0._dp

  ENDIF

  ! Transfer phi to all the others process along p
  CALL manual_3D_bcast(phi(ikxs:ikxe,ikys:ikye,izs:ize))

  ! Execution time end
  CALL cpu_time(t1_poisson)
  tc_poisson = tc_poisson + (t1_poisson - t0_poisson)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE poisson
