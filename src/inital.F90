!******************************************************************************!
!!!!!! initialize the moments and load/build coeff tables
!******************************************************************************!
SUBROUTINE inital

  USE basic
  USE model, ONLY : CO, NON_LIN
  USE prec_const
  USE coeff
  USE array, ONLY : Sepj,Sipj

  implicit none
  !!!!!! Set the moments arrays Nepj, Nipj !!!!!!
  IF ( RESTART ) THEN
    CALL init_moments
  ELSE
    CALL load_moments
  ENDIF
  !!!!!! Set phi !!!!!!
  CALL poisson
  !!!!!! Set Sepj, Sipj and dnjs coeff table !!!!!!
  IF ( NON_LIN ) THEN;
    CALL compute_Sapj(Sepj,Sipj)
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

  IF ( only_Na00 ) THEN ! Spike initialization on density only

    DO ikr=ikrs,ikre
      DO ikz=ikzs,ikze
        moments_e( 1,1, ikr,ikz, :) = initback_moments + initnoise_moments*(noise-0.5_dp)
        moments_i( 1,1, ikr,ikz, :) = initback_moments + initnoise_moments*(noise-0.5_dp)
      END DO
    END DO

  ELSE ! Broad noise initialization

    DO ikr=ikrs,ikre
      DO ikz=ikzs,ikze
        DO ip=ips_e,ipe_e
          DO ij=ijs_e,ije_e
            CALL RANDOM_NUMBER(noise)
            moments_e( ip,ij,     ikr,    ikz, :) = initback_moments + initnoise_moments*(noise-0.5_dp)
            moments_e( ip,ij, Nkz-ikz,    ikr, :) = moments_e( ip,ij,     ikr,    ikz, :) ! Symmetry
            moments_e( ip,ij,     ikz,Nkr-ikr, :) = moments_e( ip,ij,     ikr,    ikz, :) ! Symmetry
            moments_e( ip,ij, Nkz-ikz,Nkr-ikr, :) = moments_e( ip,ij,     ikr,    ikz, :) ! Symmetry
          END DO
        END DO
        DO ip=ips_i,ipe_i
          DO ij=ijs_i,ije_i
            CALL RANDOM_NUMBER(noise)
            moments_i( ip,ij,     ikr,    ikz, :) = initback_moments + initnoise_moments*(noise-0.5_dp)
            moments_i( ip,ij, Nkz-ikz,    ikr, :) = moments_i( ip,ij,     ikr,    ikz, :) ! Symmetry
            moments_i( ip,ij,     ikz,Nkr-ikr, :) = moments_i( ip,ij,     ikr,    ikz, :) ! Symmetry
            moments_i( ip,ij, Nkz-ikz,Nkr-ikr, :) = moments_i( ip,ij,     ikr,    ikz, :) ! Symmetry
          END DO
        END DO
      END DO
    END DO
  ENDIF
END SUBROUTINE init_moments
!******************************************************************************!

!******************************************************************************!
!!!!!!! Load moments from a previous save
!******************************************************************************!
SUBROUTINE load_moments
!  USE basic
!  USE initial_par
!  USE fields, ONLY : moments_e, moments_i
!  USE futils, ONLY : openf, getarr, closef
!  IMPLICIT NONE
!
!  INTEGER :: fid
!
!  CALL openf(backup_file,fid, 'r', 'D');
!  CALL getarr(fid, '/moments_e', moments_e) ! Nepj
!  CALL getarr(fid, '/moments_i', moments_i) ! Nipj
!  CALL getarr(fid, '/time', time) ! time
!  CALL closef(fid)
END SUBROUTINE
!******************************************************************************!

!******************************************************************************!
!!!!!!! Build the Laguerre-Laguerre coupling coefficient table
!******************************************************************************!
SUBROUTINE build_dnjs_table
  USE basic
  USE array, Only : dnjs
  USE fourier_grid, Only : jmaxe, jmaxi
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
  USE fourier_grid
  USE initial_par
  USE array
  use model
  USE futils, ONLY : openf, getarr, closef
  IMPLICIT NONE

  INTEGER :: ip,ij, ip2,ij2
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
