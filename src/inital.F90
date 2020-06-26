SUBROUTINE inital
  
  USE basic
  USE model, ONLY : CO
  USE prec_const
  implicit none

  WRITE(*,'(a/)') '=== Set initial conditions ==='

  CALL init_profiles

  IF (CO .EQ. -1) THEN
    WRITE(*,'(a/)') '=== Load Full Coulomb matrix ==='

    CALL load_FC_mat
  ENDIF
  !

END SUBROUTINE inital


SUBROUTINE init_profiles
  !   Set initial conditions

  USE basic
  USE fourier_grid
  USE fields
  USE initial_par
  USE time_integration

  USE prec_const
  IMPLICIT NONE

  INTEGER :: ip,ij, ikr,ikz
  INTEGER, DIMENSION(12) :: iseedarr
  REAL(dp) :: noise

  ! Seed random number generator
  iseedarr(:)=iseed
  CALL RANDOM_SEED(PUT=iseedarr)

  CALL set_updatetlevel(1)
  
  DO ikr=ikrs,ikre  
    DO ikz=ikzs,ikze

      DO ip=ips_e,ipe_e
        DO ij=ijs_e,ije_e
          CALL RANDOM_NUMBER(noise)
          !moments_e( ip,ij, ikr,ikz, :) = initback_moments + initnoise_moments*(noise-0.5_dp)
        END DO
      END DO

      DO ip=ips_i,ipe_i
        DO ij=ijs_i,ije_i
          CALL RANDOM_NUMBER(noise)
          !moments_i( ip,ij, ikr,ikz, :) = initback_moments + initnoise_moments*(noise-0.5_dp)
        END DO
      END DO
      
      ! Poke initialization on only Ne00 and Ni00
      moments_e( 1,1, ikr,ikz, :) = initback_moments + initnoise_moments*(noise-0.5_dp)
      moments_i( 1,1, ikr,ikz, :) = initback_moments + initnoise_moments*(noise-0.5_dp)
      
    END DO
  END DO
  
  CALL poisson ! To set phi

END SUBROUTINE init_profiles

SUBROUTINE load_FC_mat
  USE basic
  USE array
  IMPLICIT NONE

END SUBROUTINE load_FC_mat