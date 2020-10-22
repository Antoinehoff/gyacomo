!******************************************************************************!
!!!!!! initialize the moments and load/build coeff tables
!******************************************************************************!
SUBROUTINE inital

  USE basic
  USE model, ONLY : CO, NON_LIN
  USE prec_const
  USE coeff
  USE time_integration
  USE array, ONLY : Sepj,Sipj

  implicit none

  CALL set_updatetlevel(1)
  !!!!!! Set the moments arrays Nepj, Nipj !!!!!!
  write(*,*) 'Init moments'
  IF ( RESTART ) THEN
    CALL load_cp
  ELSE
    CALL init_moments
  ENDIF
  !!!!!! Set phi !!!!!!
  write(*,*) 'Init phi'
  CALL poisson

  !!!!!! Set Sepj, Sipj and dnjs coeff table !!!!!!
  IF ( NON_LIN .OR. (A0KH .NE. 0)) THEN;
    write(*,*) 'Init Sapj'
    CALL compute_Sapj
    WRITE(*,*) 'Building Dnjs table'
    CALL build_dnjs_table
  ENDIF
  !!!!!! Load the full coulomb collision operator coefficients !!!!!!
  IF (CO .EQ. -1) THEN
    WRITE(*,*) '=== Load Full Coulomb matrix ==='
    CALL load_FC_mat
  ENDIF
END SUBROUTINE inital
!******************************************************************************!

!******************************************************************************!
!!!!!!! Initialize the moments randomly
!******************************************************************************!
SUBROUTINE init_moments
  USE basic
  USE grid
  USE fields
  USE prec_const
  USE utility, ONLY: checkfield
  USE initial_par
  IMPLICIT NONE

  REAL(dp) :: noise
  REAL(dp) :: kr, kz, sigma, gain, kz_shift
  INTEGER, DIMENSION(12) :: iseedarr

  ! Seed random number generator
  iseedarr(:)=iseed
  CALL RANDOM_SEED(PUT=iseedarr)

  IF ( only_Na00 ) THEN ! Spike initialization on density only

    DO ikr=ikrs,ikre
      DO ikz=ikzs,ikze
        moments_e( 1,1, ikr,ikz, :) = initback_moments + initnoise_moments*(noise-0.5_dp)
        moments_i( 1,1, ikr,ikz, :) = initback_moments + initnoise_moments*(noise-0.5_dp)
      END DO
    END DO

  ELSE
    sigma    = 5._dp     ! Gaussian sigma
    gain     = 0.5_dp    ! Gaussian mean
    kz_shift = 0.5_dp    ! Gaussian z shift
    !**** Gaussian initialization (Hakim 2017) *********************************
    ! moments_i = 0; moments_e = 0;
    !   DO ikr=ikrs,ikre
    !     kr = krarray(ikr)
    !     DO ikz=ikzs,ikze
    !       kz = kzarray(ikz)
    !       moments_i( 1,1, ikr,ikz, :) = gain*sigma/SQRT2 * EXP(-(kr**2+(kz-kz_shift)**2)*sigma**2/4._dp)
    !       moments_e( 1,1, ikr,ikz, :) = gain*sigma/SQRT2 * EXP(-(kr**2+(kz-kz_shift)**2)*sigma**2/4._dp)
    !     END DO
    !   END DO

    !**** Broad noise initialization *******************************************
    DO ip=ips_e,ipe_e
      DO ij=ijs_e,ije_e
        DO ikr=ikrs,ikre
          DO ikz=ikzs,ikze
            CALL RANDOM_NUMBER(noise)
            moments_e( ip,ij,     ikr,    ikz, :) = initback_moments + initnoise_moments*(noise-0.5_dp)
            ! moments_e( ip,ij, Nkz-ikz,    ikr, :) = moments_e( ip,ij,     ikr,    ikz, :) ! Symmetry
            ! moments_e( ip,ij,     ikz,Nkr-ikr, :) = moments_e( ip,ij,     ikr,    ikz, :) ! Symmetry
            ! moments_e( ip,ij, Nkz-ikz,Nkr-ikr, :) = moments_e( ip,ij,     ikr,    ikz, :) ! Symmetry
          END DO
        END DO
      END DO
    END DO
    DO ip=ips_i,ipe_i
      DO ij=ijs_i,ije_i
        DO ikr=ikrs,ikre
          DO ikz=ikzs,ikze
            CALL RANDOM_NUMBER(noise)
            moments_i( ip,ij,     ikr,    ikz, :) = initback_moments + initnoise_moments*(noise-0.5_dp)
            ! moments_i( ip,ij, Nkz-ikz,    ikr, :) = moments_i( ip,ij,     ikr,    ikz, :) ! Symmetry
            ! moments_i( ip,ij,     ikz,Nkr-ikr, :) = moments_i( ip,ij,     ikr,    ikz, :) ! Symmetry
            ! moments_i( ip,ij, Nkz-ikz,Nkr-ikr, :) = moments_i( ip,ij,     ikr,    ikz, :) ! Symmetry
          END DO
        END DO
      END DO
    END DO

    !**** Sinusoidal phi initialization for Kelvin-Helmholtz *******************
  !   phi(1,FLOOR(0.5_dp/deltakr)) = 1._dp ! Trigger only mode kr = 1
  !   ! moments_e( :,:, 1,FLOOR(0.5_dp/deltakr), :) = 1._dp
  !   ! moments_i( :,:, 1,FLOOR(0.5_dp/deltakr), :) = 1._dp
  !
  !   DO ikr=ikrs,ikre
  !     DO ikz=ikzs,ikze
  !       CALL RANDOM_NUMBER(noise)
  !       moments_e( :,:, ikr,ikz, :) = moments_e( :,:, ikr,ikz, :) &
  !       + initnoise_moments*(noise-0.5_dp) ! adding noise
  !
  !       CALL RANDOM_NUMBER(noise)
  !       moments_i( :,:, ikr,ikz, :) = moments_i( :,:, ikr,ikz, :) &
  !       + initnoise_moments*(noise-0.5_dp) ! adding noise
  !
  !       CALL RANDOM_NUMBER(noise)
  !       phi(ikr,ikz) = phi(ikr,ikz) + initnoise_moments*(noise-0.5_dp) ! adding noise
  !     ENDDO
  !   ENDDO
  !   CALL moments_eq_rhs
  ENDIF
END SUBROUTINE init_moments
!******************************************************************************!

!******************************************************************************!
!!!!!!! Load moments from a previous save
!******************************************************************************!
SUBROUTINE load_cp
  USE basic
  USE futils,          ONLY: openf, closef, getarr, getatt
  USE grid
  USE fields
  USE diagnostics_par
  USE time_integration
  IMPLICIT NONE

  WRITE(rstfile,'(a,a1,i2.2,a3)') TRIM(rstfile0),'_',job2load,'.h5'

  WRITE(*,'(3x,a)') "Resume from previous run"

  CALL openf(rstfile, fidrst)
  CALL getatt(fidrst, '/Basic', 'cstep', cstep)
  CALL getatt(fidrst, '/Basic', 'time', time)
  CALL getatt(fidrst, '/Basic', 'jobnum', jobnum)
  jobnum = jobnum+1
  CALL getatt(fidrst, '/Basic', 'iframe2d',iframe2d)
  CALL getatt(fidrst, '/Basic', 'iframe5d',iframe5d)
  iframe2d = iframe2d-1; iframe5d = iframe5d-1

  ! Read state of system from restart file
  CALL getarr(fidrst, '/Basic/moments_e', moments_e(ips_e:ipe_e,ijs_e:ije_e,ikrs:ikre,ikzs:ikze,1),ionode=0)
  CALL getarr(fidrst, '/Basic/moments_i', moments_i(ips_i:ipe_i,ijs_i:ije_i,ikrs:ikre,ikzs:ikze,1),ionode=0)
  CALL closef(fidrst)

  WRITE(*,'(3x,a)') "Reading from restart file "//TRIM(rstfile)//" completed!"

END SUBROUTINE load_cp
!******************************************************************************!

!******************************************************************************!
!!!!!!! Build the Laguerre-Laguerre coupling coefficient table
!******************************************************************************!
SUBROUTINE build_dnjs_table
  USE basic
  USE array, Only : dnjs
  USE grid, Only : jmaxe, jmaxi
  USE coeff
  IMPLICIT NONE

  INTEGER :: in, ij, is, J
  INTEGER :: n_, j_, s_

  J = max(jmaxe,jmaxi)

  DO in = 1,J+1 ! Nested dependent loops to make benefit from dnjs symmetry
    n_ = in - 1
    DO ij = in,J+1
      j_ = ij - 1
      DO is = ij,J+1
        s_ = is - 1

        dnjs(in,ij,is) = TO_DP(ALL2L(n_,j_,s_,0))
        ! By symmetry
        dnjs(in,is,ij) = dnjs(in,ij,is)
        dnjs(ij,in,is) = dnjs(in,ij,is)
        dnjs(ij,is,in) = dnjs(in,ij,is)
        dnjs(is,ij,in) = dnjs(in,ij,is)
        dnjs(is,in,ij) = dnjs(in,ij,is)
      ENDDO
    ENDDO
  ENDDO
END SUBROUTINE
!******************************************************************************!

!******************************************************************************!
!!!!!!! Load the Full coulomb coefficient table from COSOlver results
!******************************************************************************!
SUBROUTINE load_FC_mat ! Works only for a partiular file for now with P,J <= 15,6
  USE basic
  USE grid
  USE initial_par
  USE array
  use model
  USE futils, ONLY : openf, getarr, closef
  IMPLICIT NONE

  INTEGER :: ip2,ij2
  INTEGER :: fid1, fid2, fid3, fid4

  !!!!!!!!!!!! Electron matrices !!!!!!!!!!!!
  ! get the self electron colision matrix
  CALL openf(selfmat_file,fid1, 'r', 'D');
  CALL getarr(fid1, '/Caapj/Ceepj', Ceepj) ! get array (moli format)
  CALL closef(fid1)
  ! get the Test and Back field electron ion collision matrix
  CALL openf(eimat_file,fid2, 'r', 'D');
  CALL getarr(fid2, '/Ceipj/CeipjT', CeipjT)
  CALL getarr(fid2, '/Ceipj/CeipjF', CeipjF)
  CALL closef(fid2)

  !!!!!!!!!!!!!!! Ion matrices !!!!!!!!!!!!!!
  ! get the self electron colision matrix
  CALL openf(selfmat_file, fid3, 'r', 'D');
  IF ( (pmaxe .EQ. pmaxi) .AND. (jmaxe .EQ. jmaxi) ) THEN ! if same degrees ion and electron matrices are alike so load Ceepj
    CALL getarr(fid3, '/Caapj/Ceepj', Ciipj) ! get array (moli format)
  ELSE
    CALL getarr(fid3, '/Caapj/Ciipj', Ciipj) ! get array (moli format)
  ENDIF
  CALL closef(fid3)
  ! get the Test and Back field ion electron collision matrix
  CALL openf(iemat_file, fid4, 'r', 'D');
  CALL getarr(fid4, '/Ciepj/CiepjT', CiepjT)
  CALL getarr(fid4, '/Ciepj/CiepjF', CiepjF)
  CALL closef(fid4)
END SUBROUTINE load_FC_mat
!******************************************************************************!
