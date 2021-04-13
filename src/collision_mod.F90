module collision
! contains the Hermite-Laguerre collision operators. Solved using COSOlver.
USE fields
USE array
USE basic
USE grid
USE prec_const
USE time_integration
USE model
USE utility
IMPLICIT NONE

PUBLIC :: compute_TColl
PUBLIC :: DoughertyGK_e, DoughertyGK_i
PUBLIC :: load_COSOlver_mat
PUBLIC :: apply_COSOlver_mat_e, apply_COSOlver_mat_i

CONTAINS


  !******************************************************************************!
  !! Doughtery gyrokinetic collision operator for electrons
  !******************************************************************************!
  SUBROUTINE DoughertyGK_e(ip_,ij_,ikr_,ikz_,TColl_)
    IMPLICIT NONE
    INTEGER,     INTENT(IN)    :: ip_,ij_,ikr_,ikz_
    COMPLEX(dp), INTENT(INOUT) :: TColl_

    COMPLEX(dp) :: n_,upar_,uperp_,Tpar_, Tperp_
    COMPLEX(dp) :: Dpj, Ppj, T_, ibe_
    COMPLEX(dp) :: nadiab_moment_0j, nadiab_moment_0jp1, nadiab_moment_0jm1
    INTEGER     :: in_
    REAL(dp)    :: n_dp, j_dp, p_dp, be_2, q_e_tau_e

    !** Auxiliary variables **
    p_dp      = REAL(parray_e(ip_),dp)
    j_dp      = REAL(jarray_e(ij_),dp)
    be_2      = (krarray(ikr_)**2 + kzarray(ikz_)**2) * sigmae2_taue_o2
    ibe_      = imagu*2._dp*SQRT(be_2)
    q_e_tau_e = q_e/tau_e

    !** Assembling collison operator **
    ! Velocity-space diffusion (similar to Lenhard Bernstein)
    TColl_ = -(p_dp + 2._dp*j_dp + be_2)*moments_e(ip_,ij_,ikr_,ikz_,updatetlevel)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Non zero term for p = 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF( p_dp .EQ. 0 ) THEN ! Kronecker p0
      ! Get adiabatic moment
      TColl_ = TColl_ - (p_dp + 2._dp*j_dp + be_2) * q_e_tau_e * Kernel_e(ij_,ikr_,ikz_)*phi(ikr_,ikz_)
        !** build required fluid moments **
        n_     = 0._dp
        upar_  = 0._dp; uperp_ = 0._dp
        Tpar_  = 0._dp; Tperp_ = 0._dp
        DO in_ = 1,jmaxe+1
          n_dp = REAL(in_-1,dp)
          ! Nonadiabatic moments (only different from moments when p=0)
          nadiab_moment_0j   = moments_e(1,in_  ,ikr_,ikz_,updatetlevel) + q_e_tau_e * kernel_e(in_  ,ikr_,ikz_)*phi(ikr_,ikz_)
          nadiab_moment_0jp1 = moments_e(1,in_+1,ikr_,ikz_,updatetlevel) + q_e_tau_e * kernel_e(in_+1,ikr_,ikz_)*phi(ikr_,ikz_)
          nadiab_moment_0jm1 = moments_e(1,in_-1,ikr_,ikz_,updatetlevel) + q_e_tau_e * kernel_e(in_-1,ikr_,ikz_)*phi(ikr_,ikz_)
          ! Density
          n_     = n_     + Kernel_e(in_,ikr_,ikz_) * nadiab_moment_0j
          ! Perpendicular velocity
          uperp_ = uperp_ + ibe_*0.5_dp*Kernel_e(in_,ikr_,ikz_) * (nadiab_moment_0j - nadiab_moment_0jp1)
          ! Parallel temperature
          Tpar_  = Tpar_  + Kernel_e(in_,ikr_,ikz_) * (SQRT2*moments_e(3,in_,ikr_,ikz_,updatetlevel) + nadiab_moment_0j)
          ! Perpendicular temperature
          Tperp_ = Tperp_ + Kernel_e(in_,ikr_,ikz_) * ((2._dp*n_dp+1._dp)* nadiab_moment_0j   &
                                                     -            n_dp * nadiab_moment_0jm1 &
                                                     -         (n_dp+1)* nadiab_moment_0jp1)
        ENDDO
        T_  = (Tpar_ + 2._dp*Tperp_)/3._dp - n_
      ! Add energy restoring term
      TColl_ = TColl_ + T_* 4._dp *  j_dp          * Kernel_e(ij_  ,ikr_,ikz_)
      TColl_ = TColl_ - T_* 2._dp * (j_dp + 1._dp) * Kernel_e(ij_+1,ikr_,ikz_)
      TColl_ = TColl_ - T_* 2._dp *  j_dp          * Kernel_e(ij_-1,ikr_,ikz_)
      TColl_ = TColl_ + uperp_*ibe_*((j_dp + 1._dp)* Kernel_e(ij_  ,ikr_,ikz_) &
                                   -j_dp         * Kernel_e(ij_-1,ikr_,ikz_))
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Non zero term for p = 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSEIF( p_dp .eq. 1 ) THEN ! kronecker p1
        !** build required fluid moments **
        upar_  = 0._dp
        DO in_ = 1,jmaxe+1
          ! Parallel velocity
           upar_  = upar_  + Kernel_e(in_,ikr_,ikz_) * moments_e(2,in_,ikr_,ikz_,updatetlevel)
        ENDDO
      TColl_ = TColl_ + upar_*Kernel_e(ij_,ikr_,ikz_)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Non zero term for p = 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSEIF( p_dp .eq. 2 ) THEN ! kronecker p2
        !** build required fluid moments **
        n_     = 0._dp
        Tpar_  = 0._dp; Tperp_ = 0._dp
        DO in_ = 1,jmaxe+1
          n_dp = REAL(in_-1,dp)
          ! Nonadiabatic moments (only different from moments when p=0)
          nadiab_moment_0j   = moments_e(1,in_  ,ikr_,ikz_,updatetlevel) + q_e_tau_e * kernel_e(in_  ,ikr_,ikz_)*phi(ikr_,ikz_)/tau_e
          nadiab_moment_0jp1 = moments_e(1,in_+1,ikr_,ikz_,updatetlevel) + q_e_tau_e * kernel_e(in_+1,ikr_,ikz_)*phi(ikr_,ikz_)/tau_e
          nadiab_moment_0jm1 = moments_e(1,in_-1,ikr_,ikz_,updatetlevel) + q_e_tau_e * kernel_e(in_-1,ikr_,ikz_)*phi(ikr_,ikz_)/tau_e
          ! Density
          n_     = n_     + Kernel_e(in_,ikr_,ikz_) * nadiab_moment_0j
          ! Parallel temperature
          Tpar_  = Tpar_  + Kernel_e(in_,ikr_,ikz_) * (SQRT2*moments_e(3,in_,ikr_,ikz_,updatetlevel) + nadiab_moment_0j)
          ! Perpendicular temperature
          Tperp_ = Tperp_ + Kernel_e(in_,ikr_,ikz_) * ((2._dp*n_dp+1._dp)*nadiab_moment_0j - n_dp*nadiab_moment_0jm1 - (n_dp+1)*nadiab_moment_0jp1)
        ENDDO
        T_    = (Tpar_ + 2._dp*Tperp_)/3._dp - n_
      TColl_ = TColl_ + T_*SQRT2*Kernel_e(ij_,ikr_,ikz_)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ENDIF
    ! Multiply by electron-ion collision coefficient
    TColl_ = nu_e * TColl_

  END SUBROUTINE DoughertyGK_e

  !******************************************************************************!
  !! Doughtery gyrokinetic collision operator for electrons
  SUBROUTINE DoughertyGK_i(ip_,ij_,ikr_,ikz_,TColl_)
    IMPLICIT NONE
    INTEGER,     INTENT(IN)    :: ip_,ij_,ikr_,ikz_
    COMPLEX(dp), INTENT(INOUT) :: TColl_

    COMPLEX(dp) :: n_,upar_,uperp_,Tpar_, Tperp_
    COMPLEX(dp) :: Dpj, Ppj, T_, ibi_
    COMPLEX(dp) :: nadiab_moment_0j, nadiab_moment_0jp1, nadiab_moment_0jm1
    INTEGER     :: in_
    REAL(dp)    :: n_dp, j_dp, p_dp, bi_2, q_i_tau_i

    !** Auxiliary variables **
    p_dp      = REAL(parray_i(ip_),dp)
    j_dp      = REAL(jarray_i(ij_),dp)
    bi_2      =  (krarray(ikr_)**2 + kzarray(ikz_)**2) * sigmai2_taui_o2
    ibi_      = imagu*2._dp*SQRT(bi_2)
    q_i_tau_i = q_i/tau_i

    !** Assembling collison operator **
    ! Velocity-space diffusion (similar to Lenhard Bernstein)
    TColl_ = -(p_dp + 2._dp*j_dp + bi_2)*moments_i(ip_,ij_,ikr_,ikz_,updatetlevel)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Non zero term for p = 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF( p_dp .EQ. 0 ) THEN ! Kronecker p0
      ! Get adiabatic moment
      TColl_ = TColl_ - (p_dp + 2._dp*j_dp + bi_2) * q_i_tau_i * Kernel_i(ij_,ikr_,ikz_)*phi(ikr_,ikz_)
        !** build required fluid moments **
        n_     = 0._dp
        upar_  = 0._dp; uperp_ = 0._dp
        Tpar_  = 0._dp; Tperp_ = 0._dp
        DO in_ = 1,jmaxi+1
          n_dp = REAL(in_-1,dp)
          ! Nonadiabatic moments (only different from moments when p=0)
          nadiab_moment_0j   = moments_i(1,in_  ,ikr_,ikz_,updatetlevel) + q_i_tau_i * kernel_i(in_  ,ikr_,ikz_)*phi(ikr_,ikz_)
          nadiab_moment_0jp1 = moments_i(1,in_+1,ikr_,ikz_,updatetlevel) + q_i_tau_i * kernel_i(in_+1,ikr_,ikz_)*phi(ikr_,ikz_)
          nadiab_moment_0jm1 = moments_i(1,in_-1,ikr_,ikz_,updatetlevel) + q_i_tau_i * kernel_i(in_-1,ikr_,ikz_)*phi(ikr_,ikz_)
          ! Density
          n_     = n_     + Kernel_i(in_,ikr_,ikz_) * nadiab_moment_0j
          ! Perpendicular velocity
          uperp_ = uperp_ + ibi_*0.5_dp*Kernel_i(in_,ikr_,ikz_) * (nadiab_moment_0j -nadiab_moment_0jp1)
          ! Parallel temperature
          Tpar_  = Tpar_  + Kernel_i(in_,ikr_,ikz_) * (SQRT2*moments_i(3,in_,ikr_,ikz_,updatetlevel) + nadiab_moment_0j)
          ! Perpendicular temperature
          Tperp_ = Tperp_ + Kernel_i(in_,ikr_,ikz_) * ((2._dp*n_dp+1._dp)* nadiab_moment_0j   &
                                                     -            n_dp * nadiab_moment_0jm1 &
                                                     -         (n_dp+1)* nadiab_moment_0jp1)
        ENDDO
        T_  = (Tpar_ + 2._dp*Tperp_)/3._dp - n_
      ! Add energy restoring term
      TColl_ = TColl_ + T_* 4._dp *  j_dp          * Kernel_i(ij_  ,ikr_,ikz_)
      TColl_ = TColl_ - T_* 2._dp * (j_dp + 1._dp) * Kernel_i(ij_+1,ikr_,ikz_)
      TColl_ = TColl_ - T_* 2._dp *  j_dp          * Kernel_i(ij_-1,ikr_,ikz_)
      TColl_ = TColl_ + uperp_*ibi_*((j_dp + 1._dp)* Kernel_i(ij_  ,ikr_,ikz_) &
                                   -j_dp         * Kernel_i(ij_-1,ikr_,ikz_))
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Non zero term for p = 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSEIF( p_dp .eq. 1 ) THEN ! kronecker p1
        !** build required fluid moments **
        upar_  = 0._dp
        DO in_ = 1,jmaxi+1
          ! Parallel velocity
           upar_  = upar_  + Kernel_i(in_,ikr_,ikz_) * moments_i(2,in_,ikr_,ikz_,updatetlevel)
        ENDDO
      TColl_ = TColl_ + upar_*Kernel_i(ij_,ikr_,ikz_)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Non zero term for p = 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSEIF( p_dp .eq. 2 ) THEN ! kronecker p2
        !** build required fluid moments **
        n_     = 0._dp
        Tpar_  = 0._dp; Tperp_ = 0._dp
        DO in_ = 1,jmaxi+1
          n_dp = REAL(in_-1,dp)
          ! Nonadiabatic moments (only different from moments when p=0)
          nadiab_moment_0j   = moments_i(1,in_  ,ikr_,ikz_,updatetlevel) + q_i_tau_i * kernel_i(in_  ,ikr_,ikz_)*phi(ikr_,ikz_)
          nadiab_moment_0jp1 = moments_i(1,in_+1,ikr_,ikz_,updatetlevel) + q_i_tau_i * kernel_i(in_+1,ikr_,ikz_)*phi(ikr_,ikz_)
          nadiab_moment_0jm1 = moments_i(1,in_-1,ikr_,ikz_,updatetlevel) + q_i_tau_i * kernel_i(in_-1,ikr_,ikz_)*phi(ikr_,ikz_)
          ! Density
          n_     = n_     + Kernel_i(in_,ikr_,ikz_) * nadiab_moment_0j
          ! Parallel temperature
          Tpar_  = Tpar_  + Kernel_i(in_,ikr_,ikz_) * (SQRT2*moments_i(3,in_,ikr_,ikz_,updatetlevel) + nadiab_moment_0j)
          ! Perpendicular temperature
          Tperp_ = Tperp_ + Kernel_i(in_,ikr_,ikz_) * ((2._dp*n_dp+1._dp)*nadiab_moment_0j - n_dp*nadiab_moment_0jm1 - (n_dp+1)*nadiab_moment_0jp1)
        ENDDO
        T_    = (Tpar_ + 2._dp*Tperp_)/3._dp - n_
      TColl_ = TColl_ + T_*SQRT2*Kernel_i(ij_,ikr_,ikz_)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ENDIF
    ! Multiply by ion-ion collision coefficient
    TColl_ = nu_i * TColl_

  END SUBROUTINE DoughertyGK_i

  !******************************************************************************!
  !! compute the collision terms in a (Np x Nj x Nkr x Nkz) matrix all at once
  !******************************************************************************!
  SUBROUTINE compute_TColl
    IMPLICIT NONE
    COMPLEX(dp), DIMENSION(1:pmaxe+1)   :: local_sum_e, buffer_e, total_sum_e
    COMPLEX(dp), DIMENSION(ips_e:ipe_e) :: TColl_distr_e
    COMPLEX(dp), DIMENSION(1:pmaxi+1)   :: local_sum_i, buffer_i, total_sum_i
    COMPLEX(dp), DIMENSION(ips_i:ipe_i) :: TColl_distr_i
    COMPLEX(dp) :: TColl
    INTEGER :: ikrs_C, ikre_C, ikzs_C, ikze_C
    IF (ABS(CO) .GE. 2) THEN !compute only if COSOlver matrices are used

      DO ikr = ikrs,ikre
        DO ikz = ikzs,ikze
          ! Electrons
          DO ij = 1,Jmaxe+1
            ! Loop over all p to compute sub collision term
            DO ip = 1,Pmaxe+1
              CALL apply_COSOlver_mat_e(ip,ij,ikr,ikz,TColl)
              local_sum_e(ip) = TColl
            ENDDO
            ! Sum up all the sub collision terms on root 0
            CALL MPI_REDUCE(local_sum_e, buffer_e, pmaxe+1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, comm_p, ierr)
            ! IF(rank_p .EQ. 0) THEN
            !   DO ip = 1,Pmaxe+1
            !     total_sum_e(ip) = buffer_e(ip)
            !   ENDDO
            ! ENDIF
            ! distribute the sum over the process among p
            CALL MPI_SCATTERV(buffer_e, counts_np_e, displs_np_e, MPI_DOUBLE_COMPLEX,&
                              TColl_distr_e, local_np_e, MPI_DOUBLE_COMPLEX,&
                              0, comm_p, ierr)
            DO ip = ips_e,ipe_e
              TColl_e(ip,ij,ikr,ikz) = TColl_distr_e(ip)
            ENDDO
          ENDDO
          ! Ions
          DO ij = 1,Jmaxi+1
            DO ip = 1,Pmaxi+1
              CALL apply_COSOlver_mat_i(ip,ij,ikr,ikz,TColl)
              local_sum_i(ip) = TColl
            ENDDO
            ! Reduce the local_sums to root = 0
            CALL MPI_REDUCE(local_sum_i, buffer_i, pmaxi+1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, comm_p, ierr)
            ! IF(rank_p .EQ. 0) THEN
            !   DO ip = 1,Pmaxi+1
            !     total_sum_i(ip) = buffer_i(ip)
            !   ENDDO
            ! ENDIF
            ! buffer_e contains the entire collision term along p, scatter it among
            ! the other processes (use of scatterv since Pmax/Np is not an integer)
            CALL MPI_SCATTERV(buffer_i, counts_np_i, displs_np_i, MPI_DOUBLE_COMPLEX,&
                              TColl_distr_i, local_np_i, MPI_DOUBLE_COMPLEX, &
                              0, comm_p, ierr)
            DO ip = ips_i,ipe_i
              TColl_i(ip,ij,ikr,ikz) = TColl_distr_i(ip)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    ! IF(cstep .EQ. 1) THEN
    !   write(*,*) rank_p, ': local_sum = ', local_sum_e
    !   write(*,*) rank_p, ': buffer = ', buffer_e
    !   write(*,*) rank_p, ': dist_sum = ', TColl_distr_e
    ! ENDIF

  END SUBROUTINE compute_TColl

  !******************************************************************************!
  !!!!!!! Compute ion collision term
  !******************************************************************************!
  SUBROUTINE apply_COSOlver_mat_e(ip_,ij_,ikr_,ikz_,TColl_)
    IMPLICIT NONE

    INTEGER,     INTENT(IN)  :: ip_, ij_ ,ikr_, ikz_
    COMPLEX(dp), INTENT(OUT) :: TColl_

    INTEGER     :: ip2,ij2, p_int,j_int, p2_int,j2_int, ikr_C, ikz_C
    p_int = ip_-1; j_int = ij_-1

    IF (CO .GT. 0) THEN ! GK operator (k-dependant)
      ikr_C = ikr_; ikz_C = ikz_
    ELSEIF (CO .LT. 0) THEN ! DK operator (only one mat for every k)
      ikr_C = 1;   ikz_C = 1
    ENDIF

    TColl_ = 0._dp ! Initialization of the local sum

    ! sum the electron-self and electron-ion test terms
    ploopee: DO ip2 = ips_e,ipe_e
      p2_int = parray_e(ip2)
      jloopee: DO ij2 = ijs_e,ije_e
        j2_int = jarray_e(ij2)
        TColl_ = TColl_ + moments_e(ip2,ij2,ikr_,ikz_,updatetlevel) &
           *( nu_e  * CeipjT(bare(p_int,j_int), bare(p2_int,j2_int),ikr_C, ikz_C) &
             +nu_ee * Ceepj (bare(p_int,j_int), bare(p2_int,j2_int),ikr_C, ikz_C))
      ENDDO jloopee
    ENDDO ploopee

    ! sum the electron-ion field terms
    ploopei: DO ip2 = ips_i,ipe_i
      p2_int = parray_i(ip2)
      jloopei: DO ij2 = ijs_i,ije_i
        j2_int = jarray_i(ij2)
        TColl_ = TColl_ + moments_i(ip2,ij2,ikr_,ikz_,updatetlevel) &
          *(nu_e * CeipjF(bare(p_int,j_int), bari(p2_int,j2_int),ikr_C, ikz_C))
      END DO jloopei
    ENDDO ploopei

  END SUBROUTINE apply_COSOlver_mat_e

  !******************************************************************************!
  !!!!!!! Compute ion collision term
  !******************************************************************************!
  SUBROUTINE apply_COSOlver_mat_i(ip_,ij_,ikr_,ikz_,TColl_)
    IMPLICIT NONE

    INTEGER,     INTENT(IN)    :: ip_, ij_ ,ikr_, ikz_
    COMPLEX(dp), INTENT(OUT)   :: TColl_

    INTEGER     :: ip2,ij2, p_int,j_int, p2_int,j2_int, ikr_C, ikz_C
    p_int = ip_-1; j_int = ij_-1

    IF (CO .GT. 0) THEN ! GK operator (k-dependant)
      ikr_C = ikr_; ikz_C = ikz_
    ELSEIF (CO .LT. 0) THEN ! DK operator (only one mat for every k)
      ikr_C = 1;   ikz_C = 1
    ENDIF

    TColl_ = 0._dp ! Initialization
    ! sum the ion-self and ion-electron test terms
    ploopii: DO ip2 = ips_i,ipe_i
      p2_int = parray_i(ip2)
      jloopii: DO ij2 = ijs_i,ije_i
        j2_int = jarray_i(ij2)
        TColl_ = TColl_ + moments_i(ip2,ij2,ikr_,ikz_,updatetlevel) &
            *( nu_ie * CiepjT(bari(p_int,j_int), bari(p2_int,j2_int), ikr_C, ikz_C) &
              +nu_i  * Ciipj (bari(p_int,j_int), bari(p2_int,j2_int), ikr_C, ikz_C))
      ENDDO jloopii
    ENDDO ploopii

    ploopie: DO ip2 = ips_e,ipe_e ! sum the ion-electron field terms
      p2_int = parray_e(ip2)
      jloopie: DO ij2 = ijs_e,ije_e
        j2_int = jarray_e(ij2)
        TColl_ = TColl_ + moments_e(ip2,ij2,ikr_,ikz_,updatetlevel) &
          *(nu_ie * CiepjF(bari(p_int,j_int), bare(p2_int,j2_int), ikr_C, ikz_C))
      ENDDO jloopie
    ENDDO ploopie

  END SUBROUTINE apply_COSOlver_mat_i


  !******************************************************************************!
  !!!!!!! Load the collision matrix coefficient table from COSOlver results
  !******************************************************************************!
  SUBROUTINE load_COSOlver_mat ! Load a sub matrix from iCa files (works for pmaxa,jmaxa<=P_full,J_full)
    use futils
    use initial_par
    IMPLICIT NONE
    ! Indices for row and columns of the COSOlver matrix (4D compressed 2D matrices)
    INTEGER :: irow_sub, irow_full, icol_sub, icol_full
    INTEGER :: fid1, fid2, fid3, fid4 ! file indexation

    INTEGER :: ip_e, ij_e, il_e, ik_e, ikps_C, ikpe_C                  ! indices for electrons loops
    INTEGER :: pdime, jdime                                            ! dimensions of the COSOlver matrices
    REAL(dp), DIMENSION(:,:),   ALLOCATABLE :: Ceepj_full, CeipjT_full ! To load the entire matrix
    REAL(dp), DIMENSION(:,:),   ALLOCATABLE :: CeipjF_full             ! ''
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: Ceepj__kp, CeipjT_kp    ! To store the coeff that will be used along kperp
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: CeipjF_kp               ! ''
    INTEGER :: ip_i, ij_i, il_i, ik_i                                  ! same for ions
    INTEGER :: pdimi, jdimi                                            ! .
    REAL(dp), DIMENSION(:,:),   ALLOCATABLE :: Ciipj_full, CiepjT_full ! .
    REAL(dp), DIMENSION(:,:),   ALLOCATABLE :: CiepjF_full             ! .
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: Ciipj__kp, CiepjT_kp    ! .
    REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: CiepjF_kp               ! .
    INTEGER  :: NFLR

    INTEGER  :: ikp_next, ikp_prev
    REAL(dp) ::  kp_next,  kp_prev, kperp, zerotoone

    CHARACTER(len=256) :: mat_filename, kperp_string, NFLR_string

    !! Some terminal info
    IF (CO .EQ. 2) THEN
      IF (my_id .EQ. 0) WRITE(*,*) '=== Load GK Sugama matrix ==='
    ELSEIF(CO .EQ. 3) THEN
      IF (my_id .EQ. 0) WRITE(*,*) '=== Load GK Full Coulomb matrix ==='
    ELSEIF(CO .EQ. -2) THEN
      IF (my_id .EQ. 0) WRITE(*,*) '=== Load DK Sugama matrix ==='
    ELSEIF(CO .EQ. -3) THEN
      IF (my_id .EQ. 0) WRITE(*,*) '=== Load DK Full Coulomb matrix ==='
    ENDIF

    IF (CO .GT. 0) THEN ! GK operator (k-dependant)
      ikps_C = ikps; ikpe_C = ikpe
    ELSEIF (CO .LT. 0) THEN ! DK operator (only one mat for every k)
      ikps_C = 1; ikpe_C = 1
    ENDIF

    CALL allocate_array(  Ceepj__kp, 1,(pmaxe+1)*(jmaxe+1), 1,(pmaxe+1)*(jmaxe+1), ikps_C,ikpe_C)
    CALL allocate_array(  CeipjT_kp, 1,(pmaxe+1)*(jmaxe+1), 1,(pmaxe+1)*(jmaxe+1), ikps_C,ikpe_C)
    CALL allocate_array(  CeipjF_kp, 1,(pmaxe+1)*(jmaxe+1), 1,(pmaxe+1)*(jmaxe+1), ikps_C,ikpe_C)
    CALL allocate_array(  Ciipj__kp, 1,(pmaxe+1)*(jmaxe+1), 1,(pmaxe+1)*(jmaxe+1), ikps_C,ikpe_C)
    CALL allocate_array(  CiepjT_kp, 1,(pmaxe+1)*(jmaxe+1), 1,(pmaxe+1)*(jmaxe+1), ikps_C,ikpe_C)
    CALL allocate_array(  CiepjF_kp, 1,(pmaxe+1)*(jmaxe+1), 1,(pmaxe+1)*(jmaxe+1), ikps_C,ikpe_C)

    DO ikp = ikps_C,ikpe_C ! Loop over everz kperp values
      write(kperp_string,'(f6.4)') kparray(ikp)
      NFLR = MIN(25,MAX(5,CEILING(kparray(ikp)**2)))
      write(NFLR_string,'(i2.1)') NFLR
      !!!!!!!!!!!! Electron matrices !!!!!!!!!!!!
      ! get the self electron colision matrix
      IF (CO .GT. 0) THEN
        WRITE(mat_filename,'(a,a,a,a,a,a)') TRIM(selfmat_file),&
            'NFLR_',TRIM(ADJUSTL(NFLR_string)),'_kperp_',TRIM(ADJUSTL(kperp_string)),'.h5'
      ELSE
        mat_filename = selfmat_file
      ENDIF

      WRITE(*,*) 'loading : ', mat_filename

      CALL openf(mat_filename,fid1, 'r', 'D', mpicomm=comm_p);
      CALL getatt(fid1,'/Caapj/Ceepj/','Pmaxe',pdime)
      CALL getatt(fid1,'/Caapj/Ceepj/','Jmaxe',jdime)
      IF ( ((pdime .LT. pmaxe) .OR. (jdime .LT. jmaxe)) .AND. (my_id .EQ. 0)) WRITE(*,*) '!! P,J Matrix too small !!'

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
                  Ceepj__kp (irow_sub,icol_sub,ikp) = Ceepj_full (irow_full,icol_full)
            ENDDO
            ENDDO
      ENDDO
      ENDDO

      CALL closef(fid1)
      DEALLOCATE(Ceepj_full)

      ! get the Test and Back field electron ion collision matrix
      IF (CO .GT. 0) THEN
        WRITE(mat_filename,'(a,a,a,a,a,a,a,a)') TRIM(eimat_file),&
            'NFLRe_',TRIM(ADJUSTL(NFLR_string)),'_NFLRi_',TRIM(ADJUSTL(NFLR_string)),&
            '_kperp_',TRIM(ADJUSTL(kperp_string)),'.h5'
          ELSE
        mat_filename = eimat_file
      ENDIF

      ! WRITE(*,*) 'loading : ', mat_filename

      CALL openf(mat_filename,fid2, 'r', 'D', mpicomm=comm_p);

      CALL getatt(fid2,'/Ceipj/CeipjT/','Pmaxi',pdimi)
      CALL getatt(fid2,'/Ceipj/CeipjT/','Jmaxi',jdimi)
      IF ( (pdimi .LT. pmaxi) .OR. (jdimi .LT. jmaxi) ) WRITE(*,*) '!! Sugama Matrix too small !!'

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
                  CeipjT_kp(irow_sub,icol_sub,ikp) = CeipjT_full(irow_full,icol_full)
            ENDDO
            ENDDO
            DO il_i = 0,pmaxi ! Loop over columns
            DO ik_i = 0,jmaxi
                  icol_sub  = (Jmaxi +1)*il_i + ik_i +1
                  icol_full = (jdimi +1)*il_i + ik_i +1
                  CeipjF_kp(irow_sub,icol_sub,ikp) = CeipjF_full(irow_full,icol_full)
            ENDDO
            ENDDO
      ENDDO
      ENDDO

      CALL closef(fid2)
      DEALLOCATE(CeipjF_full)
      DEALLOCATE(CeipjT_full)

      !!!!!!!!!!!!!!! Ion matrices !!!!!!!!!!!!!!
      ! get the self electron colision matrix
      IF (CO .GT. 0) THEN
        WRITE(mat_filename,'(a,a,a,a,a,a)') TRIM(selfmat_file),&
            'NFLR_',TRIM(ADJUSTL(NFLR_string)),'_kperp_',TRIM(ADJUSTL(kperp_string)),'.h5'
      ELSE
        mat_filename = selfmat_file
      ENDIF

      ! WRITE(*,*) 'loading : ', mat_filename

      CALL openf(mat_filename, fid3, 'r', 'D', mpicomm=comm_p);

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
                  Ciipj__kp (irow_sub,icol_sub,ikp) = Ciipj_full (irow_full,icol_full)
            ENDDO
            ENDDO
      ENDDO
      ENDDO

      CALL closef(fid3)
      DEALLOCATE(Ciipj_full)

      ! get the Test and Back field electron ion collision matrix
      IF (CO .GT. 0) THEN
        WRITE(mat_filename,'(a,a,a,a,a,a,a,a)') TRIM(iemat_file),&
            'NFLRe_',TRIM(ADJUSTL(NFLR_string)),'_NFLRi_',TRIM(ADJUSTL(NFLR_string)),&
            '_kperp_',TRIM(ADJUSTL(kperp_string)),'.h5'
      ELSE
        mat_filename = iemat_file
      ENDIF

      ! write(*,*) 'loading : ', mat_filename

      CALL openf(mat_filename,fid4, 'r', 'D', mpicomm=comm_p);

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
                  CiepjT_kp(irow_sub,icol_sub,ikp) = CiepjT_full(irow_full,icol_full)
            ENDDO
            ENDDO
            DO il_e = 0,pmaxe ! Loop over columns
            DO ik_e = 0,jmaxe
                  icol_sub  = (jmaxe +1)*il_e + ik_e +1
                  icol_full = (jdime +1)*il_e + ik_e +1
                  CiepjF_kp(irow_sub,icol_sub,ikp) = CiepjF_full(irow_full,icol_full)
            ENDDO
            ENDDO
      ENDDO
      ENDDO

      CALL closef(fid4)
      DEALLOCATE(CiepjF_full)
      DEALLOCATE(CiepjT_full)
    ENDDO

    IF (CO .GT. 0) THEN ! Interpolation of the kperp matrix values on kr kz grid
      IF (my_id .EQ. 0 ) WRITE(*,*) '...Interpolation from kperp to kr,kz...'
      DO ikr = ikrs,ikre
        DO ikz = ikzs,ikze
          kperp = SQRT(krarray(ikr)**2+kzarray(ikz)**2)
          ! Find correspondings previous and next kperp values of matrices grid
          ikp_prev  = INT(  FLOOR(kperp/deltakr))+1
          ikp_next  = INT(  FLOOR(kperp/deltakr))+2
          ! 0->1 variable for linear interp
          zerotoone = (kperp - kparray(ikp_prev))/(kparray(ikp_next) - kparray(ikp_prev))
          ! Linear interpolation between previous and next kperp matrix values
          Ceepj (:,:,ikr,ikz) = (Ceepj__kp(:,:,ikp_next) - Ceepj__kp(:,:,ikp_prev))*zerotoone + Ceepj__kp(:,:,ikp_prev)
          CeipjT(:,:,ikr,ikz) = (CeipjT_kp(:,:,ikp_next) - CeipjT_kp(:,:,ikp_prev))*zerotoone + CeipjT_kp(:,:,ikp_prev)
          CeipjF(:,:,ikr,ikz) = (CeipjF_kp(:,:,ikp_next) - CeipjF_kp(:,:,ikp_prev))*zerotoone + CeipjF_kp(:,:,ikp_prev)
          Ciipj (:,:,ikr,ikz) = (Ciipj__kp(:,:,ikp_next) - Ciipj__kp(:,:,ikp_prev))*zerotoone + Ciipj__kp(:,:,ikp_prev)
          CiepjT(:,:,ikr,ikz) = (CiepjT_kp(:,:,ikp_next) - CiepjT_kp(:,:,ikp_prev))*zerotoone + CiepjT_kp(:,:,ikp_prev)
          CiepjF(:,:,ikr,ikz) = (CiepjF_kp(:,:,ikp_next) - CiepjF_kp(:,:,ikp_prev))*zerotoone + CiepjF_kp(:,:,ikp_prev)
        ENDDO
      ENDDO
    ELSE ! DK -> No kperp dep, copy simply to final collision matrices
      Ceepj (:,:,1,1) = Ceepj__kp (:,:,1)
      CeipjT(:,:,1,1) = CeipjT_kp(:,:,1)
      CeipjF(:,:,1,1) = CeipjF_kp(:,:,1)
      Ciipj (:,:,1,1) = Ciipj__kp (:,:,1)
      CiepjT(:,:,1,1) = CiepjT_kp(:,:,1)
      CiepjF(:,:,1,1) = CiepjF_kp(:,:,1)
    ENDIF
    ! Deallocate auxiliary variables
    DEALLOCATE (Ceepj__kp ); DEALLOCATE (CeipjT_kp); DEALLOCATE (CeipjF_kp)
    DEALLOCATE (Ciipj__kp ); DEALLOCATE (CiepjT_kp); DEALLOCATE (CiepjF_kp)

    IF (my_id .EQ. 0) WRITE(*,*) '============DONE==========='

  END SUBROUTINE load_COSOlver_mat
  !******************************************************************************!

end module collision
