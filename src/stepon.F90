SUBROUTINE stepon
  !   Advance one time step, (num_step=4 for Runge Kutta 4 scheme)

  USE basic 
  USE time_integration
  USE fields, ONLY: moments_e, moments_i, phi
  USE array , ONLY: moments_rhs_e, moments_rhs_i
  USE fourier_grid
  USE advance_field_routine, ONLY: advance_time_level, advance_field
  USE model
  USE utility, ONLY: checkfield

  use prec_const
  IMPLICIT NONE

  INTEGER :: num_step, ip,ij


   DO num_step=1,ntimelevel ! eg RK4 compute successively k1, k2, k3, k4

      ! Compute right hand side
      !WRITE (*,*) 'Compute right hand side .. nstep = ', num_step
      CALL moments_eq_rhs

      !WRITE (*,*) 'Advance time level .. nstep = ', num_step
      CALL advance_time_level ! Advance from updatetlevel to updatetlevel+1

      !WRITE (*,*) 'Advance field .. nstep = ', num_step

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

      ! Solving Poisson equation
      !WRITE (*,*) 'Solving Poisson equation .. nstep = ', num_step
      CALL poisson
      CALL checkfield_all()

   END DO

   CONTAINS

      SUBROUTINE checkfield_all ! Check all the fields for inf or nan
        IF(.NOT.nlend) THEN
           nlend=nlend .or. checkfield(phi,' phi')
           DO ip=ips_e,ipe_e
             DO ij=ijs_e,ije_e
              nlend=nlend .or. checkfield(moments_e(ip,ij,:,:,updatetlevel),' moments_e')
             ENDDO
           ENDDO

           DO ip=ips_i,ipe_i
             DO ij=ijs_i,ije_i
              nlend=nlend .or. checkfield(moments_i(ip,ij,:,:,updatetlevel),' moments_i')
             ENDDO
           ENDDO
        ENDIF
      END SUBROUTINE checkfield_all

END SUBROUTINE stepon
