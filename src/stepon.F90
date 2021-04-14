SUBROUTINE stepon
  !   Advance one time step, (num_step=4 for Runge Kutta 4 scheme)
  USE advance_field_routine, ONLY: advance_time_level, advance_field
  USE array , ONLY: moments_rhs_e, moments_rhs_i, Sepj, Sipj
  USE basic
  USE closure
  USE collision, ONLY : compute_TColl
  USE fields, ONLY: moments_e, moments_i, phi
  USE ghosts
  USE grid
  USE model
  use prec_const
  USE time_integration
  USE utility, ONLY: checkfield

  IMPLICIT NONE

  INTEGER :: num_step
  LOGICAL :: mlend

   DO num_step=1,ntimelevel ! eg RK4 compute successively k1, k2, k3, k4
   !----- BEFORE: All fields are updated for step = n

      ! Compute right hand side from current fields
      ! N_rhs(N_n,phi_n, S_n, Tcoll_n)
      CALL moments_eq_rhs_e
      CALL moments_eq_rhs_i

      ! ---- step n -> n+1 transition
      ! Advance from updatetlevel to updatetlevel+1 (according to num. scheme)
      CALL advance_time_level
      ! Update moments with the hierarchy RHS (step by step)
      ! N_n+1 = N_n + N_rhs(n)
      CALL advance_moments
      ! Closure enforcement of N_n+1
      CALL apply_closure_model
      ! Exchanges the ghosts values of N_n+1
      CALL update_ghosts

      ! Update collision C_n+1 = C(N_n+1)
      CALL compute_TColl

      ! Update electrostatic potential phi_n = phi(N_n+1)
      CALL poisson

      ! Update nonlinear term S_n -> S_n+1(phi_n+1,N_n+1)
      IF ( NON_LIN ) THEN
        CALL compute_Sapj
      ENDIF

      !-  Check before next step
      CALL checkfield_all()
      IF( nlend ) EXIT ! exit do loop

      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

   !----- AFTER: All fields are updated for step = n+1
   END DO

   CONTAINS

      SUBROUTINE advance_moments
        ! Execution time start
        CALL cpu_time(t0_adv_field)

        DO ip=ips_e,ipe_e
          DO ij=ijs_e,ije_e
            CALL advance_field(moments_e(ip,ij,:,:,:),moments_rhs_e(ip,ij,:,:,:))
          ENDDO
        ENDDO
        DO ip=ips_i,ipe_i
          DO ij=ijs_i,ije_i
            CALL advance_field(moments_i(ip,ij,:,:,:),moments_rhs_i(ip,ij,:,:,:))
          ENDDO
        ENDDO

        ! Execution time end
        CALL cpu_time(t1_adv_field)
        tc_adv_field = tc_adv_field + (t1_adv_field - t0_adv_field)
      END SUBROUTINE advance_moments


      SUBROUTINE checkfield_all ! Check all the fields for inf or nan
        ! Execution time start
        CALL cpu_time(t0_checkfield)

        IF ( (ikrs .EQ. 1) .AND. (NON_LIN) ) CALL enforce_symetry() ! Enforcing symmetry on kr = 0

        mlend=.FALSE.
        IF(.NOT.nlend) THEN
           mlend=mlend .or. checkfield(phi,' phi')
           DO ip=ips_e,ipe_e
             DO ij=ijs_e,ije_e
              mlend=mlend .or. checkfield(moments_e(ip,ij,:,:,updatetlevel),' moments_e')
             ENDDO
           ENDDO
           DO ip=ips_i,ipe_i
             DO ij=ijs_i,ije_i
              mlend=mlend .or. checkfield(moments_i(ip,ij,:,:,updatetlevel),' moments_i')
             ENDDO
           ENDDO
           CALL MPI_ALLREDUCE(mlend, nlend, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
        ENDIF

        ! Execution time end
        CALL cpu_time(t1_checkfield)
        tc_checkfield = tc_checkfield + (t1_checkfield - t0_checkfield)
      END SUBROUTINE checkfield_all

      SUBROUTINE anti_aliasing
        DO ip=ips_e,ipe_e
          DO ij=ijs_e,ije_e
            DO ikr=ikrs,ikre
              DO ikz=ikzs,ikze
                moments_e( ip,ij,ikr,ikz,:) = AA_r(ikr)* AA_z(ikz) * moments_e( ip,ij,ikr,ikz,:)
              END DO
            END DO
          END DO
        END DO
        DO ip=ips_i,ipe_i
          DO ij=ijs_i,ije_i
            DO ikr=ikrs,ikre
              DO ikz=ikzs,ikze
                moments_i( ip,ij,ikr,ikz,:) = AA_r(ikr)* AA_z(ikz) * moments_i( ip,ij,ikr,ikz,:)
              END DO
            END DO
          END DO
        END DO
      END SUBROUTINE anti_aliasing

      SUBROUTINE enforce_symetry ! Force X(k) = X(N-k)* complex conjugate symmetry
        IF ( contains_kr0 ) THEN
          ! Electron moments
          DO ip=ips_e,ipe_e
            DO ij=ijs_e,ije_e
              DO ikz=2,Nkz/2 !symmetry at kr = 0
                moments_e( ip,ij,ikr_0,ikz, :) = CONJG(moments_e( ip,ij,ikr_0,Nkz+2-ikz, :))
              END DO
              ! must be real at origin
              moments_e(ip,ij, ikr_0,ikz_0, :) = REAL(moments_e(ip,ij, ikr_0,ikz_0, :))
            END DO
          END DO
          ! Ion moments
          DO ip=ips_i,ipe_i
            DO ij=ijs_i,ije_i
              DO ikz=2,Nkz/2 !symmetry at kr = 0
                moments_i( ip,ij,ikr_0,ikz, :) = CONJG(moments_i( ip,ij,ikr_0,Nkz+2-ikz, :))
              END DO
              ! must be real at origin and top right
              moments_i(ip,ij, ikr_0,ikz_0, :) = REAL(moments_i(ip,ij, ikr_0,ikz_0, :))
            END DO
          END DO
          ! Phi
          DO ikz=2,Nkz/2 !symmetry at kr = 0
            phi(ikr_0,ikz) = phi(ikr_0,Nkz+2-ikz)
          END DO
          ! must be real at origin
          phi(ikr_0,ikz_0) = REAL(phi(ikr_0,ikz_0))
        ENDIF
      END SUBROUTINE enforce_symetry

END SUBROUTINE stepon
