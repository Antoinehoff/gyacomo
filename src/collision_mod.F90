module collision
! contains the Hermite-Laguerre collision operators. Solved using COSOlver.

implicit none

PUBLIC :: load_FC_mat
PUBLIC :: DoughertyGK_e, DoughertyGK_i

CONTAINS

  !******************************************************************************!
  SUBROUTINE LenardBernsteinDK

  END SUBROUTINE LenardBernsteinDK

  !******************************************************************************!
  !! Doughtery gyrokinetic collision operator for electrons
  SUBROUTINE DoughertyGK_e(ip,ij,ikr,ikz,TColl)
    USE fields,           ONLY: moments_e, phi
    USE array,            ONLY: Kernel_e
    USE basic
    USE grid,             ONLY: jmaxe, parray_e, jarray_e, krarray, kzarray
    USE prec_const
    USE time_integration, ONLY: updatetlevel
    USE model,            ONLY: sigmae2_taue_o2, tau_e, nu_e
    IMPLICIT NONE
    INTEGER, INTENT(IN)      :: ip,ij,ikr,ikz
    COMPLEX(dp), INTENT(INOUT) :: TColl

    COMPLEX(dp) :: n_,upar_,uperp_,Tpar_, Tperp_
    COMPLEX(dp) :: Dpj, Ppj, T_, ibe_
    COMPLEX(dp) :: adiab_moment_0j, adiab_moment_0jp1, adiab_moment_0jm1
    COMPLEX(dp) :: adiab_moment_pj
    INTEGER     :: in_
    REAL        :: n_dp, j_dp, p_dp, be_2

    !** Auxiliary variables **
    p_dp  = REAL(parray_e(ij),dp)
    j_dp  = REAL(jarray_e(ij),dp)
    be_2  =  (krarray(ikr)**2 + kzarray(ikz)**2) * sigmae2_taue_o2
    ibe_  = imagu*2._dp*SQRT(be_2)

    !** build fluid moments used for collison operator **
    n_     = 0._dp
    upar_  = 0._dp; uperp_ = 0._dp
    Tpar_  = 0._dp; Tperp_ = 0._dp
    DO in_ = 1,jmaxe+1
      n_dp = REAL(in_-1,dp)

      ! Adiabatic moments (only different from moments when p=0)
      adiab_moment_0j = moments_e(1,in_,ikr,ikz,updatetlevel) + kernel_e(in_,ikr,ikz)*phi(ikr,ikz)/tau_e
      adiab_moment_0jp1 = 0._dp; adiab_moment_0jm1 = 0._dp
      IF (n_dp+1 .LE. jmaxe)& ! Truncation
        adiab_moment_0jp1 = moments_e(1,in_+1,ikr,ikz,updatetlevel) + kernel_e(in_+1,ikr,ikz)*phi(ikr,ikz)/tau_e
      IF (n_dp-1 .GE. 0)& ! Truncation
        adiab_moment_0jm1 = moments_e(1,in_-1,ikr,ikz,updatetlevel) + kernel_e(in_-1,ikr,ikz)*phi(ikr,ikz)/tau_e

      ! Density
      n_     = n_     + Kernel_e(in_,ikr,ikz) * adiab_moment_0j
      ! Parallel velocity
      upar_  = upar_  + Kernel_e(in_,ikr,ikz) * moments_e(2,in_,ikr,ikz,updatetlevel)
      ! Perpendicular velocity
      uperp_ = uperp_ + ibe_*0.5_dp*Kernel_e(in_,ikr,ikz) * (adiab_moment_0j -adiab_moment_0jp1)
      ! Parallel temperature
      Tpar_  = Tpar_  + Kernel_e(in_,ikr,ikz) * (SQRT2*moments_e(3,in_,ikr,ikz,updatetlevel) + adiab_moment_0j)
      ! Perpendicular temperature
      Tperp_ = Tperp_ + Kernel_e(in_,ikr,ikz) * ((2._dp*n_dp+1._dp)*adiab_moment_0j - n_dp*adiab_moment_0jm1 - (n_dp+1)*adiab_moment_0jp1)

    ENDDO
    T_    = (Tpar_ + 2._dp*Tperp_)/3._dp - n_

    !** Assembling collison operator **
    ! Velocity-space diffusion
    TColl = -(2._dp*p_dp + j_dp + be_2)*moments_e(ip,ij,ikr,ikz,updatetlevel)
    IF( p_dp .EQ. 0 ) & ! Get adiabatic moment
      TColl = TColl - (2._dp*p_dp + j_dp + be_2) * Kernel_e(ij,ikr,ikz)*phi(ikr,ikz)/tau_e

    ! Add energy restoring term
    IF( p_dp .eq. 0 ) THEN ! kronecker p0
      TColl = TColl + T_* 4._dp*j_dp * Kernel_e(ij,ikr,ikz)
      IF( j_dp+1 .LE. jmaxe ) & ! Truncation
        TColl = TColl - T_*2._dp*(j_dp + 1._dp)*Kernel_e(ij+1,ikr,ikz)
      IF( j_dp-1 .GE. 0 )  & ! Truncation
        TColl = TColl - T_*2._dp*j_dp*Kernel_e(ij-1,ikr,ikz)
    ENDIF
    IF( p_dp .eq. 2 ) & ! kronecker p2
      TColl = TColl + T_*SQRT2*Kernel_e(ij,ikr,ikz)

    !Add momentum restoring term
    IF( p_dp .eq. 1 ) & ! kronecker p1
      TColl = TColl + upar_*Kernel_e(ij,ikr,ikz)
    IF( p_dp .eq. 0 ) & ! kronecker p0
      TColl = TColl + uperp_*ibe_*( (j_dp + 1._dp)*Kernel_e(ij,ikr,ikz) - j_dp*Kernel_e(ij-1,ikr,ikz))

    ! Multiply by electron-ion collision coefficient
    TColl = nu_e * TColl

  END SUBROUTINE DoughertyGK_e

  !******************************************************************************!
  !! Doughtery gyrokinetic collision operator for ions
  SUBROUTINE DoughertyGK_i(ip,ij,ikr,ikz,TColl)
    USE fields,           ONLY: moments_i, phi
    USE array,            ONLY: Kernel_i
    USE basic
    USE grid,             ONLY: jmaxi, parray_i, jarray_i, krarray, kzarray
    USE prec_const
    USE time_integration, ONLY: updatetlevel
    USE model,            ONLY: sigmai2_taui_o2, tau_i, nu_i
    IMPLICIT NONE
    INTEGER, INTENT(IN)      :: ip,ij,ikr,ikz
    COMPLEX(dp), INTENT(INOUT) :: TColl

    COMPLEX(dp) :: n_,upar_,uperp_,Tpar_, Tperp_
    COMPLEX(dp) :: Dpj, Ppj, T_, ibi_
    COMPLEX(dp) :: adiab_moment_0j, adiab_moment_0jp1, adiab_moment_0jm1
    COMPLEX(dp) :: adiab_moment_pj
    INTEGER     :: in_
    REAL        :: n_dp, j_dp, p_dp, bi_2

    !** Auxiliary variables **
    p_dp  = REAL(parray_i(ij),dp)
    j_dp  = REAL(jarray_i(ij),dp)
    bi_2  =  (krarray(ikr)**2 + kzarray(ikz)**2) * sigmai2_taui_o2
    ibi_  = imagu*2._dp*SQRT(bi_2)

    !** build fluid moments used for collison operator **
    n_     = 0._dp
    upar_  = 0._dp; uperp_ = 0._dp
    Tpar_  = 0._dp; Tperp_ = 0._dp
    DO in_ = 1,jmaxi+1
      n_dp = REAL(in_-1,dp)

      ! Adiabatic moments (only different from moments when j=0)
      adiab_moment_0j = moments_i(1,in_,ikr,ikz,updatetlevel) + Kernel_i(in_  ,ikr,ikz)*phi(ikr,ikz)/tau_i
      adiab_moment_0jp1 = 0._dp; adiab_moment_0jm1 = 0._dp
      IF (n_dp+1 .LE. jmaxi)&
       adiab_moment_0jp1 = moments_i(1,in_+1,ikr,ikz,updatetlevel) + Kernel_i(in_+1,ikr,ikz)*phi(ikr,ikz)/tau_i
      IF (n_dp-1 .GE. 0)&
       adiab_moment_0jm1 = moments_i(1,in_-1,ikr,ikz,updatetlevel) + Kernel_i(in_-1,ikr,ikz)*phi(ikr,ikz)/tau_i

      ! Density
      n_     = n_     + Kernel_i(in_,ikr,ikz) * adiab_moment_0j
      ! Parallel velocity
      upar_  = upar_  + Kernel_i(in_,ikr,ikz) * moments_i(2,in_,ikr,ikz,updatetlevel)
      ! Perpendicular velocity
      uperp_ = uperp_ + ibi_*0.5_dp*Kernel_i(in_,ikr,ikz) * (adiab_moment_0j -adiab_moment_0jp1)
      ! Parallel temperature
      Tpar_  = Tpar_  + Kernel_i(in_,ikr,ikz) * (SQRT2*moments_i(3,in_,ikr,ikz,updatetlevel) + adiab_moment_0j)
      ! Perpendicular temperature
      Tperp_ = Tperp_ + Kernel_i(in_,ikr,ikz) * ((2._dp*n_dp+1._dp)*adiab_moment_0j - n_dp*adiab_moment_0jm1 - (n_dp+1)*adiab_moment_0jp1)

    ENDDO
    T_    = (Tpar_ + 2._dp*Tperp_)/3._dp - n_

    !** Assembling collison operator **
    ! Velocity-space diffusion
    TColl = -(2._dp*p_dp + j_dp + bi_2)*moments_i(ip,ij,ikr,ikz,updatetlevel)
    IF( p_dp .EQ. 0 ) & ! Get adiabatic moment
      TColl = TColl - (2._dp*p_dp + j_dp + bi_2) * Kernel_i(ij,ikr,ikz)*phi(ikr,ikz)/tau_i

    ! Add energy restoring term
    IF( p_dp .EQ. 0 ) THEN
      TColl = TColl + T_* 4._dp*j_dp * Kernel_i(ij,ikr,ikz)
      IF( j_dp+1 .LE. jmaxi ) &
        TColl = TColl - T_*2._dp*(j_dp + 1._dp)*Kernel_i(ij+1,ikr,ikz)
      IF( j_dp-1 .GE. 0 )     &
        TColl = TColl - T_*2._dp*j_dp*Kernel_i(ij-1,ikr,ikz)
    ENDIF
    IF( p_dp .eq. 2 ) &
      TColl = TColl + T_*SQRT2*Kernel_i(ij,ikr,ikz)

    ! Add momentum restoring term
    IF( p_dp .eq. 1 ) &
      TColl = TColl + upar_*Kernel_i(ij,ikr,ikz)
    IF( p_dp .eq. 0 ) &
      TColl = TColl + uperp_*ibi_*( (j_dp + 1._dp)*Kernel_i(ij,ikr,ikz) - j_dp*Kernel_i(ij-1,ikr,ikz))
    ! Multiply by ion-ion collision coefficient
    TColl = nu_i * TColl

  END SUBROUTINE DoughertyGK_i

  !******************************************************************************!
  !!!!!!! Compute ion Full Coulomb collision operator
  !******************************************************************************!
  SUBROUTINE FullCoulombDK_e(p_int,j_int,ikr,ikz,TColl)
    USE basic
    USE fields,           ONLY: moments_i, moments_e
    USE array,            ONLY: Ceepj, CeipjT, CeipjF
    USE grid,             ONLY: pmaxe,jmaxe, pmaxi,jmaxi, parray_e, parray_i, &
                                jarray_e, jarray_i, bari, bare
    USE prec_const
    USE time_integration, ONLY: updatetlevel
    USE model,            ONLY: nu_e, nu_ee
    IMPLICIT NONE

    INTEGER,     INTENT(IN)    :: p_int,j_int ,ikr,ikz
    COMPLEX(dp), INTENT(INOUT) :: TColl

    INTEGER     :: ip2,ij2, p2_int,j2_int

    TColl = 0._dp ! Initialization

    ploopee: DO ip2 = 1,pmaxe+1 ! sum the electron-self and electron-ion test terms
      p2_int = parray_e(ip2)
      jloopee: DO ij2 = 1,jmaxe+1
        j2_int = jarray_e(ij2)
        TColl = TColl + moments_e(ip2,ij2,ikr,ikz,updatetlevel) &
           *( nu_e  * CeipjT(bare(p_int,j_int), bare(p2_int,j2_int)) &
             +nu_ee * Ceepj (bare(p_int,j_int), bare(p2_int,j2_int)))
      ENDDO jloopee
    ENDDO ploopee
    ploopei: DO ip2 = 1,pmaxi+1 ! sum the electron-ion field terms
      p2_int = parray_i(ip2)
      jloopei: DO ij2 = 1,jmaxi+1
        j2_int = jarray_i(ij2)
        TColl = TColl + moments_i(ip2,ij2,ikr,ikz,updatetlevel) &
          *(nu_e * CeipjF(bare(p_int,j_int), bari(p2_int,j2_int)))
      END DO jloopei
    ENDDO ploopei

  END SUBROUTINE FullCoulombDK_e

  !******************************************************************************!
  !!!!!!! Compute ion Full Coulomb collision operator
  !******************************************************************************!
  SUBROUTINE FullCoulombDK_i(p_int,j_int,ikr,ikz,TColl)
    USE basic
    USE fields,           ONLY: moments_i, moments_e
    USE array,            ONLY: Ciipj, CiepjT, CiepjF
    USE grid,             ONLY: pmaxe,jmaxe, pmaxi,jmaxi, parray_e, parray_i, &
                                jarray_e, jarray_i, bari, bare
    USE prec_const
    USE time_integration, ONLY: updatetlevel
    USE model,            ONLY: nu_i, nu_ie
    IMPLICIT NONE

    INTEGER,     INTENT(IN)    :: p_int,j_int,ikr,ikz
    COMPLEX(dp), INTENT(INOUT) :: TColl

    INTEGER     :: ip2,ij2, p2_int,j2_int

    TColl = 0._dp ! Initialization

    ploopii: DO ip2 = 1,pmaxi+1 ! sum the ion-self and ion-electron test terms
      p2_int = parray_i(ip2)
      jloopii: DO ij2 = 1,jmaxi+1
        j2_int = jarray_i(ij2)
        TColl = TColl + moments_i(ip2,ij2,ikr,ikz,updatetlevel) &
            *( nu_ie * CiepjT(bari(p_int,j_int), bari(p2_int,j2_int)) &
              +nu_i  * Ciipj (bari(p_int,j_int), bari(p2_int,j2_int)))
      ENDDO jloopii
    ENDDO ploopii

    ploopie: DO ip2 = 1,pmaxe+1 ! sum the ion-electron field terms
      p2_int = parray_e(ip2)
      jloopie: DO ij2 = 1,jmaxe+1
        j2_int = jarray_e(ij2)
        TColl = TColl + moments_e(ip2,ij2,ikr,ikz,updatetlevel) &
          *(nu_ie * CiepjF(bari(p_int,j_int), bare(p2_int,j2_int)))
      ENDDO jloopie
    ENDDO ploopie

  END SUBROUTINE FullCoulombDK_i


  !******************************************************************************!
  !!!!!!! Load the Full coulomb coefficient table from COSOlver results
  !******************************************************************************!
  SUBROUTINE load_FC_mat ! Load a sub matrix from iCa files (works for pmaxa,jmaxa<=P_full,J_full)
    use prec_const
    use fields
    use array
    use grid
    use basic
    use futils
    use initial_par
    use model
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
