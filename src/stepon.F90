SUBROUTINE stepon
  !   Advance one time step, (num_step=4 for Runge Kutta 4 scheme)

  USE basic 
  USE time_integration
  USE fields, ONLY: theta, temp, vpar, moments, phi
  USE array , ONLY: theta_rhs, temp_rhs, vpar_rhs, moments_rhs
  USE space_grid
  USE advance_field_routine, ONLY: advance_time_level, advance_field
  USE fft
  USE model
  USE utility, ONLY: checkfield

  use prec_const
  IMPLICIT NONE

  INTEGER :: num_step, ip


   DO num_step=1,ntimelevel ! eg RK4 compute successively k1, k2, k3, k4

      CALL evaluation_auxfield_total ! Compute gradients, etc.

      ! Compute right hand side of derivative the fields
      if (.not. freeze_theta)   CALL theta_eq_rhs
      if (.not. freeze_temp)    CALL temp_eq_rhs
      if (.not. freeze_vpar)    CALL vpar_eq_rhs
      if (.not. freeze_moments) CALL moments_eq_rhs

      CALL advance_time_level ! Advance from updatetlevel to updatetlevel+1

      if (.not. freeze_theta)     CALL advance_field(theta,theta_rhs)
      if (.not. freeze_temp)      CALL advance_field(temp,temp_rhs)
      if (.not. freeze_vpar)      CALL advance_field(vpar,vpar_rhs)
      do ip=ips,ipe
         do ij=ijs,ije
            if (.not. freeze_moments) CALL advance_field(moments(ip,ij,:,:,:),moments_rhs(ip,ij,:,:,:))
         enddo
      enddo

      if (.not. freeze_phi) CALL poisson
      CALL checkfield_all()
   END DO


   if ( fft_suppress_high_freq ) then ! should have updatetlevel = 1, only suppress high frequences after the runge kutta step
      call suppress_high_freq(theta(:,updatetlevel))
      call suppress_high_freq(temp(:,updatetlevel))
      call suppress_high_freq(vpar(:,updatetlevel))
      do ip=ips,ipe
         call suppress_high_freq(moments(ip,:,updatetlevel))
      enddo
   endif

   CONTAINS

      SUBROUTINE checkfield_all ! Check all the fields for inf or nan
        IF(.NOT.nlend) THEN
           nlend=nlend .or.checkfield(theta(:,updatetlevel),' theta')
           nlend=nlend .or. checkfield(temp(:,updatetlevel),' temp')
           nlend=nlend .or. checkfield(vpar(:,updatetlevel),' vpar')
           nlend=nlend .or. checkfield(phi,' phi')
           do ip=ips,ipe
             nlend=nlend .or. checkfield(moments(ip,:,updatetlevel),' moments')
           enddo
        ENDIF
      END SUBROUTINE checkfield_all

END SUBROUTINE stepon
