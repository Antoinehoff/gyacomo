SUBROUTINE stepon
  !   Advance one time step, (num_step=4 for Runge Kutta 4 scheme)
  USE advance_field_routine, ONLY: advance_time_level, advance_field, advance_moments
  USE basic
  USE closure
  USE ghosts, ONLY: update_ghosts
  USE grid
  USE model, ONLY : LINEARITY, KIN_E
  use prec_const
  USE time_integration
  USE numerics, ONLY: play_with_modes
  USE utility, ONLY: checkfield
  IMPLICIT NONE

  INTEGER :: num_step
  LOGICAL :: mlend

   DO num_step=1,ntimelevel ! eg RK4 compute successively k1, k2, k3, k4
   !----- BEFORE: All fields+ghosts are updated for step = n
      ! Compute right hand side from current fields
      ! N_rhs(N_n, nadia_n, phi_n, S_n, Tcoll_n)
      CALL assemble_RHS

      ! ---- step n -> n+1 transition
      ! Advance from updatetlevel to updatetlevel+1 (according to num. scheme)
      CALL advance_time_level
      ! ----

      ! Update moments with the hierarchy RHS (step by step)
      ! N_n+1 = N_n + N_rhs(n)
      CALL advance_moments
      ! Closure enforcement of moments
      CALL apply_closure_model
      ! Exchanges the ghosts values of N_n+1
      CALL update_ghosts

      ! Update electrostatic potential phi_n = phi(N_n+1)
      CALL poisson
      CALL update_ghosts_z_phi

      ! Numerical experiments
      ! Store or cancel/maintain zonal modes artificially
      ! CALL play_with_modes

      !-  Check before next step
      CALL checkfield_all()
      IF( nlend ) EXIT ! exit do loop

      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
   !----- AFTER: All fields are updated for step = n+1
   END DO

   CONTAINS
     !!!! Basic structure to simplify stepon
     SUBROUTINE assemble_RHS
       USE moments_eq_rhs, ONLY: compute_moments_eq_rhs
       USE collision,      ONLY: compute_TColl
       USE nonlinear,      ONLY: compute_Sapj
       USE processing,     ONLY: compute_nadiab_moments_z_gradients_and_interp
       IMPLICIT NONE
         ! compute auxiliary non adiabatic moments array, gradients and interp
         CALL compute_nadiab_moments_z_gradients_and_interp
         ! compute nonlinear term ("if linear" is included inside)
         CALL compute_Sapj
         ! compute collision term ("if coll, if nu >0" is included inside)
         CALL compute_TColl
         ! compute the moments equation rhs
         CALL compute_moments_eq_rhs
     END SUBROUTINE assemble_RHS

      SUBROUTINE checkfield_all ! Check all the fields for inf or nan
        ! Execution time start
        CALL cpu_time(t0_checkfield)

        IF(LINEARITY .NE. 'linear') CALL anti_aliasing   ! ensure 0 mode for 2/3 rule
        IF(LINEARITY .NE. 'linear') CALL enforce_symmetry ! Enforcing symmetry on kx = 0

        mlend=.FALSE.
        IF(.NOT.nlend) THEN
           mlend=mlend .or. checkfield(phi,' phi')
           IF(KIN_E) THEN
           DO ij=ijgs_e,ijge_e
             DO ip=ipgs_e,ipge_e
              mlend=mlend .or. checkfield(moments_e(ip,ij,:,:,:,updatetlevel),' moments_e')
             ENDDO
           ENDDO
           ENDIF
           DO ij=ijgs_i,ijge_i
             DO ip=ipgs_i,ipge_i
              mlend=mlend .or. checkfield(moments_i(ip,ij,:,:,:,updatetlevel),' moments_i')
             ENDDO
           ENDDO
           CALL MPI_ALLREDUCE(mlend, nlend, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
        ENDIF

        ! Execution time end
        CALL cpu_time(t1_checkfield)
        tc_checkfield = tc_checkfield + (t1_checkfield - t0_checkfield)
      END SUBROUTINE checkfield_all

      SUBROUTINE anti_aliasing
        IF(KIN_E)THEN
        DO iz=izgs,izge
          DO ikx=ikxs,ikxe
            DO iky=ikys,ikye
              DO ij=ijgs_e,ijge_e
                DO ip=ipgs_e,ipge_e
                  moments_e( ip,ij,iky,ikx,iz,:) = AA_x(ikx)* AA_y(iky) * moments_e( ip,ij,iky,ikx,iz,:)
                END DO
              END DO
            END DO
          END DO
        END DO
        ENDIF
        DO iz=izgs,izge
          DO ikx=ikxs,ikxe
            DO iky=ikys,ikye
              DO ij=ijs_i,ije_i
                DO ip=ipgs_i,ipge_i
                  moments_i( ip,ij,iky,ikx,iz,:) = AA_x(ikx)* AA_y(iky) * moments_i( ip,ij,iky,ikx,iz,:)
                END DO
              END DO
            END DO
          END DO
        END DO
      END SUBROUTINE anti_aliasing

      SUBROUTINE enforce_symmetry ! Force X(k) = X(N-k)* complex conjugate symmetry
        IF ( contains_ky0 ) THEN
          ! Electron moments
          IF(KIN_E) THEN
            DO iz=izs,ize
              DO ij=ijgs_e,ijge_e
                DO ip=ipgs_e,ipge_e
                  DO ikx=2,Nkx/2 !symmetry at ky = 0
                    moments_e( ip,ij,iky_0,ikx,iz, :) = CONJG(moments_e( ip,ij,iky_0,Nkx+2-iky,iz, :))
                  END DO
                ! must be real at origin
                moments_e(ip,ij, iky_0,ikx_0,iz, :) = REAL(moments_e(ip,ij, iky_0,ikx_0,iz, :))
              END DO
            END DO
          END DO
          ENDIF
          ! Ion moments
          DO iz=izs,ize
            DO ij=ijgs_i,ijge_i
              DO ip=ipgs_i,ipge_i
                DO ikx=2,Nkx/2 !symmetry at ky = 0
                  moments_i( ip,ij,iky_0,ikx,iz, :) = CONJG(moments_i( ip,ij,iky_0,Nkx+2-ikx,iz, :))
                END DO
                ! must be real at origin and top right
                moments_i(ip,ij, iky_0,ikx_0,iz, :) = REAL(moments_i(ip,ij, iky_0,ikx_0,iz, :))
              END DO
            END DO
          END DO
          ! Phi
          DO ikx=2,Nkx/2 !symmetry at ky = 0
            phi(iky_0,ikx,izs:ize) = phi(iky_0,Nkx+2-ikx,izs:ize)
          END DO
          ! must be real at origin
          phi(iky_0,ikx_0,izs:ize) = REAL(phi(iky_0,ikx_0,izs:ize))
        ENDIF
      END SUBROUTINE enforce_symmetry

END SUBROUTINE stepon
