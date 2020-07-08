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

SUBROUTINE load_FC_mat ! Works only for a partiular file for now with P,J <= 15,6
  USE basic
  USE fourier_grid
  USE initial_par
  USE array
  use model
  USE futils, ONLY : openf, getarr, closef
  IMPLICIT NONE

  INTEGER :: ip,ij, ip2,ij2
  INTEGER :: fid1, fid2, fid3, fid4

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!! Electron matrices !!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  WRITE(*,*) 'Load self electron collision matrix..'
  ! get the self electron colision matrix
  CALL openf(selfmat_file,fid1, 'r', 'D');
  CALL getarr(fid1, '/Caapj/Ceepj', Ceepj) ! get array (moli format)
  CALL closef(fid1)

  ! get the Test and Back field electron ion collision matrix
  WRITE(*,*) 'Load test + field electron collision matrix..'
  CALL openf(eimat_file,fid2, 'r', 'D');
  CALL getarr(fid2, '/Ceipj/CeipjT', CeipjT)
  CALL getarr(fid2, '/Ceipj/CeipjF', CeipjF)
  CALL closef(fid2)

  write(25,*) eimat_file
  write(25,*) 'Ceepj(3,5)'
  write(25,*) Ceepj(3,5)
  write(25,*) 'CeipjT(12,12)'
  write(25,*) CeipjT(12,12)
  write(25,*) 'CeipjF(12,2)'
  write(25,*) CeipjF(12,2)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!! Ion matrices !!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  WRITE(*,*) 'Load self ion collision matrix..'
  ! get the self electron colision matrix
  CALL openf(selfmat_file, fid3, 'r', 'D');
  IF ( (pmaxe .EQ. pmaxi) .AND. (jmaxe .EQ. jmaxi) ) THEN ! if same degrees ion and electron matrices are alike so load Ceepj
    CALL getarr(fid3, '/Caapj/Ceepj', Ciipj) ! get array (moli format)
  ELSE
    CALL getarr(fid3, '/Caapj/Ciipj', Ciipj) ! get array (moli format)
  ENDIF
  CALL closef(fid3)
  ! get the Test and Back field ion electron collision matrix
  WRITE(*,*) 'Load test + field ion collision matrix..'
  CALL openf(iemat_file, fid4, 'r', 'D');
  CALL getarr(fid4, '/Ciepj/CiepjT', CiepjT)
  CALL getarr(fid4, '/Ciepj/CiepjF', CiepjF)
  CALL closef(fid4)

  write(25,*) iemat_file
  write(25,*) 'Ciipj(3,5)'
  write(25,*) Ciipj(3,5)
  write(25,*) 'CiepjT(10,10)'
  write(25,*) CiepjT(10,10)
  write(25,*) 'CiepjF(10,2)'
  write(25,*) CiepjF(10,2)

END SUBROUTINE load_FC_mat
