SUBROUTINE stepon
  !   Advance one time step, (num_step=4 for Runge Kutta 4 scheme)

  USE basic 
  USE time_integration
  USE fields, ONLY: moments, phi
  USE array , ONLY: moments_rhs
  USE fourier_grid
  USE advance_field_routine, ONLY: advance_time_level, advance_field
  USE model
  USE utility, ONLY: checkfield

  use prec_const
  IMPLICIT NONE

  INTEGER :: num_step, ipj


   DO num_step=1,ntimelevel ! eg RK4 compute successively k1, k2, k3, k4

      ! Compute right hand side
      CALL moments_eq_rhs

      CALL advance_time_level ! Advance from updatetlevel to updatetlevel+1

      do ipj=ipjs,ipje
            CALL advance_field(moments(ipj,:,:,:),moments_rhs(ipj,:,:,:))
      enddo

      ! Solving Poisson equation
      CALL poisson
      CALL checkfield_all()

   END DO

   CONTAINS

      SUBROUTINE checkfield_all ! Check all the fields for inf or nan
        IF(.NOT.nlend) THEN
           nlend=nlend .or. checkfield(phi,' phi')
           do ipj=ipjs,ipje
             nlend=nlend .or. checkfield(moments(ipj,:,:,updatetlevel),' moments')
           enddo
        ENDIF
      END SUBROUTINE checkfield_all

END SUBROUTINE stepon
