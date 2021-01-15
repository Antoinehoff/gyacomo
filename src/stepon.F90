SUBROUTINE stepon
  !   Advance one time step, (num_step=4 for Runge Kutta 4 scheme)

  USE basic
  USE time_integration
  USE fields, ONLY: moments_e, moments_i, phi
  USE array , ONLY: moments_rhs_e, moments_rhs_i, Sepj, Sipj
  USE grid
  USE advance_field_routine, ONLY: advance_time_level, advance_field
  USE model
  USE closure
  USE utility, ONLY: checkfield
  use prec_const
  IMPLICIT NONE

  INTEGER :: num_step
  LOGICAL :: mlend

   DO num_step=1,ntimelevel ! eg RK4 compute successively k1, k2, k3, k4
      ! Compute right hand side of moments hierarchy equation
      CALL moments_eq_rhs_e
      CALL moments_eq_rhs_i

      ! Advance from updatetlevel to updatetlevel+1 (according to num. scheme)
      CALL advance_time_level

      ! Update the moments with the hierarchy RHS (step by step)

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

      ! Closure enforcement
      CALL apply_closure_model

      ! Execution time end
      CALL cpu_time(t1_adv_field)
      tc_adv_field = tc_adv_field + (t1_adv_field - t0_adv_field)

      ! Update electrostatic potential
      CALL poisson

      ! Update nonlinear term
      IF ( NON_LIN ) THEN
        CALL compute_Sapj
      ENDIF

      !(The two routines above are called in inital for t=0)
      CALL checkfield_all()
      IF( nlend ) EXIT ! exit do loop

   END DO

   CONTAINS

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

      SUBROUTINE enforce_symetry
        DO ip=ips_e,ipe_e
          DO ij=ijs_e,ije_e
            DO ikz=2,Nkz/2 !symmetry at kr = 0
              moments_e( ip,ij,1,ikz, :) = CONJG(moments_e( ip,ij,1,Nkz+2-ikz, :))
            END DO
          END DO
        END DO
        DO ip=ips_i,ipe_i
          DO ij=ijs_i,ije_i
            DO ikz=2,Nkz/2 !symmetry at kr = 0
              moments_i( ip,ij,1,ikz, :) = CONJG(moments_i( ip,ij,1,Nkz+2-ikz, :))
            END DO
          END DO
        END DO
        DO ikz=2,Nkz/2 !symmetry at kr = 0
          phi(1,ikz) = phi(1,Nkz+2-ikz)
        END DO
      END SUBROUTINE enforce_symetry

END SUBROUTINE stepon
