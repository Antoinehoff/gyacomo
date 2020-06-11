SUBROUTINE inital
  
  use basic

  use prec_const
  implicit none

  WRITE(*,'(a/)') '=== Set initial conditions ==='

  CALL init_profiles

END SUBROUTINE inital


SUBROUTINE init_profiles
  !   Set initial conditions

  USE basic
  USE space_grid
  USE fields
  use initial_par
  USE time_integration
  USE model, ONLY: gradient_scheme

  use prec_const
  IMPLICIT NONE

  INTEGER :: ip, ij, ikr, ikz
  integer, dimension(12) :: iseedarr
  real(dp) :: noiser, noisez

  ! Seed random number generator
  iseedarr(:)=iseed
  CALL RANDOM_SEED(PUT=iseedarr)


  CALL set_updatetlevel(1)

  
  DO ikr=ikrs,ikre  
    DO ikz=ikzs,ikze
      
      CALL RANDOM_NUMBER(noiser)
      CALL RANDOM_NUMBER(noisez)
      theta(ikr,ikz,:) = log( initback_density &
          + init_ampli_density * sin(init_nb_oscil_density*TWOPI/(krmax-krmin)*rarray(ikr)) + initnoise_density*(noiser-0.5_dp) & ! sine profile with noise
          + init_ampli_density * sin(init_nb_oscil_density*TWOPI/(kzmax-kzmin)*zarray(ikz)) + initnoise_density*(noisez-0.5_dp) ) ! sine profile with noise
      CALL RANDOM_NUMBER(noiser)
      CALL RANDOM_NUMBER(noisez)
      temp(ikr,ikz,:) = log( initback_temp + init_ampli_temp * sin(init_nb_oscil_temp*TWOPI/(zmax-zmin)*zarray(iz)) + initnoise_temp*(noise-0.5_dp) ) ! sine profile with noise

      CALL RANDOM_NUMBER(noiser)
      CALL RANDOM_NUMBER(noisez)
      vpar(ikr,ikz,:) = initback_vpar + initnoise_vpar*(noisez-0.5_dp)

      DO ip=ips,ipe
        DO ij=ijs,ije
          CALL RANDOM_NUMBER(noiser)
          CALL RANDOM_NUMBER(noisez)
          moments(ip,ij,ikr,ikz,:) = initback_moments + initnoise_moments*(noisez-0.5_dp)
        END DO
      END DO

    END DO
  END DO
  
  call poisson ! To set phi

END SUBROUTINE init_profiles

