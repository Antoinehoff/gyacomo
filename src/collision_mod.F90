module collision
! contains the Hermite-Laguerre collision operators. Solved using COSOlver.

use prec_const
use fields
use array
use grid
use basic
use futils
use initial_par
use model

implicit none

PUBLIC :: load_FC_mat

CONTAINS

!******************************************************************************!
SUBROUTINE LenardBernsteinDK

END SUBROUTINE LenardBernsteinDK

!******************************************************************************!
SUBROUTINE DoughertyDK

END SUBROUTINE DoughertyDK

!******************************************************************************!
SUBROUTINE FullCoulombDK

END SUBROUTINE FullCoulombDK


!******************************************************************************!
!!!!!!! Load the Full coulomb coefficient table from COSOlver results
!******************************************************************************!
SUBROUTINE load_FC_mat ! Load a sub matrix from iCa files (works for pmaxa,jmaxa<=P_full,J_full)
  IMPLICIT NONE

  INTEGER :: irow_sub, irow_full, icol_sub, icol_full
  INTEGER :: fid1, fid2, fid3, fid4

  INTEGER :: ip_e, ij_e, il_e, ik_e
  INTEGER :: pdime, jdime
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: Ceepj_full, CeipjT_full
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: CeipjF_full
  INTEGER :: ip_i, ij_i, il_i, ik_i
  INTEGER :: pdimi, jdimi
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: Ciipj_full, CiepjT_full
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: CiepjF_full

  !!!!!!!!!!!! Electron matrices !!!!!!!!!!!!
  IF (my_id .EQ. 0) WRITE(*,*) 'Load ee FC mat...'
  ! get the self electron colision matrix
  CALL openf(selfmat_file,fid1, 'r', 'D',mpicomm=MPI_COMM_WORLD)

  CALL getatt(fid1,'/Caapj/Ceepj/','Pmaxe',pdime)
  CALL getatt(fid1,'/Caapj/Ceepj/','Jmaxe',jdime)

  IF ( ((pdime .LT. pmaxe) .OR. (jdime .LT. jmaxe)) .AND. (my_id .EQ. 0)) WRITE(*,*) '!! FC Matrix too small !!'

  CALL allocate_array(  Ceepj_full, 1,(pdime+1)*(jdime+1), 1,(pdime+1)*(jdime+1))
  CALL getarr(fid1, '/Caapj/Ceepj', Ceepj_full) ! get array (moli format)

  ! Fill sub array with only usefull polynmials degree
  DO ip_e = 0,pmaxe ! Loop over rows
  DO ij_e = 0,jmaxe
        irow_sub  = (jmaxe +1)*ip_e + ij_e +1
        irow_full = (jdime +1)*ip_e + ij_e +1
        DO il_e = 0,pmaxe ! Loop over columns
        DO ik_e = 0,jmaxe
              icol_sub  = (jmaxe +1)*il_e + ik_e +1
              icol_full = (jdime +1)*il_e + ik_e +1
              Ceepj (irow_sub,icol_sub) = Ceepj_full (irow_full,icol_full)
        ENDDO
        ENDDO
  ENDDO
  ENDDO

  CALL closef(fid1)
  DEALLOCATE(Ceepj_full)

  IF (my_id .EQ. 0) WRITE(*,*) 'Load ei FC mat...'
  ! get the Test and Back field electron ion collision matrix
  CALL openf(eimat_file,fid2, 'r', 'D');

  CALL getatt(fid2,'/Ceipj/CeipjT/','Pmaxi',pdimi)
  CALL getatt(fid2,'/Ceipj/CeipjT/','Jmaxi',jdimi)
  IF ( (pdimi .LT. pmaxi) .OR. (jdimi .LT. jmaxi) ) WRITE(*,*) '!! FC Matrix too small !!'

  CALL allocate_array( CeipjT_full, 1,(pdime+1)*(jdime+1), 1,(pdime+1)*(jdime+1))
  CALL allocate_array( CeipjF_full, 1,(pdime+1)*(jdime+1), 1,(pdimi+1)*(jdimi+1))

  CALL getarr(fid2, '/Ceipj/CeipjT', CeipjT_full)
  CALL getarr(fid2, '/Ceipj/CeipjF', CeipjF_full)

  ! Fill sub array with only usefull polynmials degree
  DO ip_e = 0,pmaxe ! Loop over rows
  DO ij_e = 0,jmaxe
        irow_sub  = (jmaxe +1)*ip_e + ij_e +1
        irow_full = (jdime +1)*ip_e + ij_e +1
        DO il_e = 0,pmaxe ! Loop over columns
        DO ik_e = 0,jmaxe
              icol_sub  = (jmaxe +1)*il_e + ik_e +1
              icol_full = (jdime +1)*il_e + ik_e +1
              CeipjT(irow_sub,icol_sub) = CeipjT_full(irow_full,icol_full)
        ENDDO
        ENDDO
        DO il_i = 0,pmaxi ! Loop over columns
        DO ik_i = 0,jmaxi
              icol_sub  = (Jmaxi +1)*il_i + ik_i +1
              icol_full = (jdimi +1)*il_i + ik_i +1
              CeipjF(irow_sub,icol_sub) = CeipjF_full(irow_full,icol_full)
        ENDDO
        ENDDO
  ENDDO
  ENDDO

  CALL closef(fid2)
  DEALLOCATE(CeipjF_full)
  DEALLOCATE(CeipjT_full)

  !!!!!!!!!!!!!!! Ion matrices !!!!!!!!!!!!!!
  IF (my_id .EQ. 0) WRITE(*,*) 'Load ii FC mat...'
  ! get the self electron colision matrix
  CALL openf(selfmat_file, fid3, 'r', 'D',mpicomm=MPI_COMM_WORLD);

  CALL allocate_array(  Ciipj_full, 1,(pdimi+1)*(jdimi+1), 1,(pdimi+1)*(jdimi+1))

  IF ( (pmaxe .EQ. pmaxi) .AND. (jmaxe .EQ. jmaxi) ) THEN ! if same degrees ion and electron matrices are alike so load Ceepj
    CALL getarr(fid3, '/Caapj/Ceepj', Ciipj_full) ! get array (moli format)
  ELSE
    CALL getarr(fid3, '/Caapj/Ciipj', Ciipj_full) ! get array (moli format)
  ENDIF

  ! Fill sub array with only usefull polynmials degree
  DO ip_i = 0,Pmaxi ! Loop over rows
  DO ij_i = 0,Jmaxi
        irow_sub  = (Jmaxi +1)*ip_i + ij_i +1
        irow_full = (jdimi +1)*ip_i + ij_i +1
        DO il_i = 0,Pmaxi ! Loop over columns
        DO ik_i = 0,Jmaxi
              icol_sub  = (Jmaxi +1)*il_i + ik_i +1
              icol_full = (jdimi +1)*il_i + ik_i +1
              Ciipj (irow_sub,icol_sub) = Ciipj_full (irow_full,icol_full)
        ENDDO
        ENDDO
  ENDDO
  ENDDO

  CALL closef(fid3)
  DEALLOCATE(Ciipj_full)

  ! get the Test and Back field electron ion collision matrix
  CALL openf(iemat_file,fid4, 'r', 'D');

  CALL allocate_array( CiepjT_full, 1,(pdimi+1)*(jdimi+1), 1,(pdimi+1)*(jdimi+1))
  CALL allocate_array( CiepjF_full, 1,(pdimi+1)*(jdimi+1), 1,(pdime+1)*(jdime+1))

  IF (my_id .EQ. 0) WRITE(*,*) 'Load ie FC mat...'
  CALL getarr(fid4, '/Ciepj/CiepjT', CiepjT_full)
  CALL getarr(fid4, '/Ciepj/CiepjF', CiepjF_full)

  ! Fill sub array with only usefull polynmials degree
  DO ip_i = 0,Pmaxi ! Loop over rows
  DO ij_i = 0,Jmaxi
        irow_sub  = (Jmaxi +1)*ip_i + ij_i +1
        irow_full = (jdimi +1)*ip_i + ij_i +1
        DO il_i = 0,Pmaxi ! Loop over columns
        DO ik_i = 0,Jmaxi
              icol_sub  = (Jmaxi +1)*il_i + ik_i +1
              icol_full = (jdimi +1)*il_i + ik_i +1
              CiepjT(irow_sub,icol_sub) = CiepjT_full(irow_full,icol_full)
        ENDDO
        ENDDO
        DO il_e = 0,pmaxe ! Loop over columns
        DO ik_e = 0,jmaxe
              icol_sub  = (jmaxe +1)*il_e + ik_e +1
              icol_full = (jdime +1)*il_e + ik_e +1
              CiepjF(irow_sub,icol_sub) = CiepjF_full(irow_full,icol_full)
        ENDDO
        ENDDO
  ENDDO
  ENDDO

  CALL closef(fid4)
  DEALLOCATE(CiepjF_full)
  DEALLOCATE(CiepjT_full)

END SUBROUTINE load_FC_mat
!******************************************************************************!

end module collision
