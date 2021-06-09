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
PUBLIC :: DoughertyGK_e, DoughertyGK_i!, DoughertyGK
PUBLIC :: load_COSOlver_mat
PUBLIC :: apply_COSOlver_mat_e, apply_COSOlver_mat_i

CONTAINS

  ! !******************************************************************************!
  ! !! Doughtery gyrokinetic collision operator for electrons
  ! !******************************************************************************!
  ! SUBROUTINE DoughertyGK(ip_,ij_,ikr_,ikz_,TColl_,specie_)
  !   IMPLICIT NONE
  !   INTEGER,     INTENT(IN)        :: ip_,ij_,ikr_,ikz_
  !   CHARACTER(len = 1), INTENT(IN) :: specie_
  !   COMPLEX(dp), INTENT(INOUT)     :: TColl_
  !
  !   COMPLEX(dp) :: n_,upar_,uperp_,Tpar_, Tperp_
  !   COMPLEX(dp) :: Dpj, Ppj, T_
  !   COMPLEX(dp) :: nadiab_moment_0j
  !   COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE :: moments_
  !   REAL(dp),    DIMENSION(:),   ALLOCATABLE :: kernel_
  !   REAL(dp)    :: Knp0, Knp1, Knm1
  !   INTEGER     :: in_, jmax_
  !   REAL(dp)    :: n_dp, j_dp, p_dp, b_, bo2_2_, q_tau_, nu_
  !
  !   !** Auxiliary variables **
  !   !! If electrons !!
  !   IF ( specie_ .EQ. 'e' ) THEN
  !     p_dp      = REAL(parray_e(ip_),dp)
  !     j_dp      = REAL(jarray_e(ij_),dp)
  !     jmax_     = jmaxe
  !     bo2_2_    = (krarray(ikr_)**2 + kzarray(ikz_)**2) * sigmae2_taue_o2 ! this is (be/2)^2
  !     b_        = 2_dp*SQRT(bo2_2_) ! this is be
  !     q_tau_   = q_e/tau_e
  !     nu_       = nu_ee
  !     CALL allocate_array(moments_, ips_e,ipe_e, ijs_e,ije_e)
  !     moments_(ips_e:ipe_e,ijs_e:ije_e) = moments_e(ips_e:ipe_e,ijs_e:ije_e,ikr_,ikz_,updatetlevel)
  !     CALL allocate_array(kernel_, ijs_e,ije_e)
  !     kernel_(ijsg_e:ijeg_e)    = kernel_e(ijsg_e:ijeg_e,ikr_,ikz_)
  !   !! If ions !!
  !   ELSEIF ( specie_ .EQ. 'i') THEN
  !     p_dp      = REAL(parray_i(ip_),dp)
  !     j_dp      = REAL(jarray_i(ij_),dp)
  !     jmax_     = jmaxi
  !     bo2_2_    = (krarray(ikr_)**2 + kzarray(ikz_)**2) * sigmai2_taui_o2 ! this is (bi/2)^2
  !     b_        = 2_dp*SQRT(bo2_2_) ! this is be
  !     q_tau_   = q_i/tau_i
  !     nu_       = nu_i
  !     CALL allocate_array(moments_, ips_i,ipe_i, ijs_i,ije_i)
  !     moments_(ips_i:ipe_i,ijs_i:ije_i) = moments_i(ips_i:ipe_i,ijs_i:ije_i,ikr_,ikz_,updatetlevel)
  !     CALL allocate_array(kernel_, ijsg_i,ijeg_i)
  !     kernel_(ijsg_i:ijeg_i)    = kernel_i(ijsg_i:ijeg_i,ikr_,ikz_)
  !   ENDIF
  !
  !   !** Assembling collison operator **
  !   ! Velocity-space diffusion (similar to Lenhard Bernstein)
  !   ! -nuee (p + 2j + b^2/2) Nepj
  !   TColl_ = -(p_dp + 2._dp*j_dp + 2._dp*bo2_2_)*moments_(ip_,ij_)
  !
  !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Non zero term for p = 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   IF( p_dp .EQ. 0 ) THEN ! Kronecker p0
  !     ! Get adiabatic moment
  !     TColl_ = TColl_ - (p_dp + 2._dp*j_dp + 2._dp*bo2_2_) * q_tau_ * kernel_(ij_)*phi(ikr_,ikz_)
  !       !** build required fluid moments **
  !       n_     = 0._dp
  !       upar_  = 0._dp; uperp_ = 0._dp
  !       Tpar_  = 0._dp; Tperp_ = 0._dp
  !       DO in_ = 1,jmaxe+1
  !         n_dp = REAL(in_-1,dp)
  !         ! Store the kernels for sparing readings
  !         Knp0 =  kernel_(in_)
  !         Knp1 =  kernel_(in_+1)
  !         Knm1 =  kernel_(in_-1)
  !         ! Nonadiabatic moments (only different from moments when p=0)
  !         nadiab_moment_0j   = moments_(1,in_) + q_tau_ * Knp0 *phi(ikr_,ikz_)
  !         ! Density
  !         n_     = n_     + Knp0 * nadiab_moment_0j
  !         ! Perpendicular velocity
  !         uperp_ = uperp_ + b_*0.5_dp*(Knp0 - Knm1) * nadiab_moment_0j
  !         ! Parallel temperature
  !         Tpar_  = Tpar_  + Knp0 * (SQRT2*moments_(3,in_) + nadiab_moment_0j)
  !         ! Perpendicular temperature
  !         Tperp_ = Tperp_ + ((2._dp*n_dp+1._dp)*Knp0 - (n_dp+1._dp)*Knp1 - n_dp*Knm1)*nadiab_moment_0j
  !       ENDDO
  !     T_  = (Tpar_ + 2._dp*Tperp_)/3._dp - n_
  !     ! Add energy restoring term
  !     TColl_ = TColl_ + T_* 4._dp *  j_dp          * kernel_(ij_)
  !     TColl_ = TColl_ - T_* 2._dp * (j_dp + 1._dp) * kernel_(ij_+1)
  !     TColl_ = TColl_ - T_* 2._dp *  j_dp          * kernel_(ij_-1)
  !     TColl_ = TColl_ + uperp_*b_*  (kernel_(ij_) - kernel_(ij_-1))
  !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Non zero term for p = 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   ELSEIF( p_dp .eq. 1 ) THEN ! kronecker p1
  !     !** build required fluid moments **
  !     upar_  = 0._dp
  !     DO in_ = 1,jmax_+1
  !       ! Parallel velocity
  !        upar_  = upar_  + Kernel_(in_) * moments_(2,in_)
  !     ENDDO
  !     TColl_ = TColl_ + upar_*Kernel_(ij_)
  !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Non zero term for p = 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !   ELSEIF( p_dp .eq. 2 ) THEN ! kronecker p2
  !     !** build required fluid moments **
  !     n_     = 0._dp
  !     upar_  = 0._dp; uperp_ = 0._dp
  !     Tpar_  = 0._dp; Tperp_ = 0._dp
  !     DO in_ = 1,jmaxe+1
  !       n_dp = REAL(in_-1,dp)
  !       ! Store the kernels for sparing readings
  !       Knp0 =  kernel_(in_)
  !       Knp1 =  kernel_(in_+1)
  !       Knm1 =  kernel_(in_-1)
  !       ! Nonadiabatic moments (only different from moments when p=0)
  !       nadiab_moment_0j   = moments_(1,in_) + q_tau_*Knp0*phi(ikr_,ikz_)
  !       ! Density
  !       n_     = n_     + Knp0 * nadiab_moment_0j
  !       ! Parallel temperature
  !       Tpar_  = Tpar_  + Knp0 * (SQRT2*moments_(3,in_) + nadiab_moment_0j)
  !       ! Perpendicular temperature
  !       Tperp_ = Tperp_ + ((2._dp*n_dp+1._dp)*Knp0 - (n_dp+1) * Knp1 - n_dp * Knm1)*nadiab_moment_0j
  !     ENDDO
  !     T_  = (Tpar_ + 2._dp*Tperp_)/3._dp - n_
  !     TColl_ = TColl_ + T_*SQRT2*kernel_(ij_)
  !   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !   ENDIF
  !   ! Multiply by specieslike collision coefficient
  !   TColl_ = nu_ * TColl_
  !
  ! END SUBROUTINE DoughertyGK

  !******************************************************************************!
  !! Doughtery gyrokinetic collision operator for electrons
  !******************************************************************************!
  SUBROUTINE DoughertyGK_e(ip_,ij_,ikr_,ikz_,TColl_)
    IMPLICIT NONE
    INTEGER,     INTENT(IN)    :: ip_,ij_,ikr_,ikz_
    COMPLEX(dp), INTENT(INOUT) :: TColl_

    COMPLEX(dp) :: n_,upar_,uperp_,Tpar_, Tperp_
    COMPLEX(dp) :: Dpj, Ppj, T_, be_
    COMPLEX(dp) :: nadiab_moment_0j
    REAL(dp)    :: Knp0, Knp1, Knm1
    INTEGER     :: in_
    REAL(dp)    :: n_dp, j_dp, p_dp, be_2, q_e_tau_e

    !** Auxiliary variables **
    p_dp      = REAL(parray_e(ip_),dp)
    j_dp      = REAL(jarray_e(ij_),dp)
    be_2      = (krarray(ikr_)**2 + kzarray(ikz_)**2) * sigmae2_taue_o2 ! this is (be/2)^2
    ! ibe_       = imagu*2._dp*SQRT(be_2)
    be_       = 2_dp*SQRT(be_2) ! this is be
    q_e_tau_e = q_e/tau_e

    !** Assembling collison operator **
    ! Velocity-space diffusion (similar to Lenhard Bernstein)
    ! -nuee (p + 2j + b^2/2) Nepj
    TColl_ = -(p_dp + 2._dp*j_dp + 2._dp*be_2)*moments_e(ip_,ij_,ikr_,ikz_,updatetlevel)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Non zero term for p = 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF( p_dp .EQ. 0 ) THEN ! Kronecker p0
      ! Get adiabatic moment
      TColl_ = TColl_ - (p_dp + 2._dp*j_dp + 2._dp*be_2) * q_e_tau_e * Kernel_e(ij_,ikr_,ikz_)*phi(ikr_,ikz_)
        !** build required fluid moments **
        n_     = 0._dp
        upar_  = 0._dp; uperp_ = 0._dp
        Tpar_  = 0._dp; Tperp_ = 0._dp
        DO in_ = 1,jmaxe+1
          n_dp = REAL(in_-1,dp)
          ! Store the kernels for sparing readings
          Knp0 =  Kernel_e(in_,ikr_,ikz_)
          Knp1 =  Kernel_e(in_+1,ikr_,ikz_)
          Knm1 =  Kernel_e(in_-1,ikr_,ikz_)
          ! Nonadiabatic moments (only different from moments when p=0)
          nadiab_moment_0j   = moments_e(1,in_  ,ikr_,ikz_,updatetlevel) + q_e_tau_e * Knp0 *phi(ikr_,ikz_)
          ! Density
          n_     = n_     + Knp0 * nadiab_moment_0j
          ! Perpendicular velocity
          uperp_ = uperp_ + be_*0.5_dp*(Knp0 - Knm1) * nadiab_moment_0j
          ! Parallel temperature
          Tpar_  = Tpar_  + Knp0 * (SQRT2*moments_e(3,in_,ikr_,ikz_,updatetlevel) + nadiab_moment_0j)
          ! Perpendicular temperature
          Tperp_ = Tperp_ + ((2._dp*n_dp+1._dp)*Knp0 - (n_dp+1._dp)*Knp1 - n_dp*Knm1)*nadiab_moment_0j
        ENDDO
      T_  = (Tpar_ + 2._dp*Tperp_)/3._dp - n_
      ! Add energy restoring term
      TColl_ = TColl_ + T_* 4._dp *  j_dp          * Kernel_e(ij_,ikr_,ikz_)
      TColl_ = TColl_ - T_* 2._dp * (j_dp + 1._dp) * Kernel_e(ij_+1,ikr_,ikz_)
      TColl_ = TColl_ - T_* 2._dp *  j_dp          * Kernel_e(ij_-1,ikr_,ikz_)
      TColl_ = TColl_ + uperp_*be_* (Kernel_e(ij_,ikr_,ikz_) - Kernel_e(ij_-1,ikr_,ikz_))
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
      upar_  = 0._dp; uperp_ = 0._dp
      Tpar_  = 0._dp; Tperp_ = 0._dp
      DO in_ = 1,jmaxe+1
        n_dp = REAL(in_-1,dp)
        ! Store the kernels for sparing readings
        Knp0 =  Kernel_e(in_,ikr_,ikz_)
        Knp1 =  Kernel_e(in_+1,ikr_,ikz_)
        Knm1 =  Kernel_e(in_-1,ikr_,ikz_)
        ! Nonadiabatic moments (only different from moments when p=0)
        nadiab_moment_0j   = moments_e(1,in_  ,ikr_,ikz_,updatetlevel) + q_e_tau_e*Knp0*phi(ikr_,ikz_)
        ! Density
        n_     = n_     + Knp0 * nadiab_moment_0j
        ! Parallel temperature
        Tpar_  = Tpar_  + Knp0 * (SQRT2*moments_e(3,in_,ikr_,ikz_,updatetlevel) + nadiab_moment_0j)
        ! Perpendicular temperature
        Tperp_ = Tperp_ + ((2._dp*n_dp+1._dp)*Knp0 - (n_dp+1) * Knp1 - n_dp * Knm1)*nadiab_moment_0j
      ENDDO
      T_  = (Tpar_ + 2._dp*Tperp_)/3._dp - n_
      TColl_ = TColl_ + T_*SQRT2*Kernel_e(ij_,ikr_,ikz_)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ENDIF
    ! Multiply by electron-electron collision coefficient
    TColl_ = nu_ee * TColl_

  END SUBROUTINE DoughertyGK_e
  !******************************************************************************!
  !! Doughtery gyrokinetic collision operator for ions
  !******************************************************************************!
  SUBROUTINE DoughertyGK_i(ip_,ij_,ikr_,ikz_,TColl_)
    IMPLICIT NONE
    INTEGER,     INTENT(IN)    :: ip_,ij_,ikr_,ikz_
    COMPLEX(dp), INTENT(INOUT) :: TColl_

    COMPLEX(dp) :: n_,upar_,uperp_,Tpar_, Tperp_
    COMPLEX(dp) :: Dpj, Ppj, T_, bi_
    COMPLEX(dp) :: nadiab_moment_0j
    REAL(dp)    :: Knp0, Knp1, Knm1
    INTEGER     :: in_
    REAL(dp)    :: n_dp, j_dp, p_dp, bi_2, q_i_tau_i

    !** Auxiliary variables **
    p_dp      = REAL(parray_i(ip_),dp)
    j_dp      = REAL(jarray_i(ij_),dp)
    bi_2      = (krarray(ikr_)**2 + kzarray(ikz_)**2) * sigmai2_taui_o2 ! this is (bi/2)^2
    bi_       = 2_dp*SQRT(bi_2) ! this is be
    q_i_tau_i = q_i/tau_i

    !** Assembling collison operator **
    ! Velocity-space diffusion (similar to Lenhard Bernstein)
    ! -nui (p + 2j + b^2/2) Nipj
    TColl_ = -(p_dp + 2._dp*j_dp + 2._dp*bi_2)*moments_i(ip_,ij_,ikr_,ikz_,updatetlevel)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Non zero term for p = 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF( p_dp .EQ. 0 ) THEN ! Kronecker p0
      ! Get adiabatic moment
      TColl_ = TColl_ - (p_dp + 2._dp*j_dp + 2._dp*bi_2) * q_i_tau_i * Kernel_i(ij_,ikr_,ikz_)*phi(ikr_,ikz_)
        !** build required fluid moments **
        n_     = 0._dp
        upar_  = 0._dp; uperp_ = 0._dp
        Tpar_  = 0._dp; Tperp_ = 0._dp
        DO in_ = 1,jmaxi+1
          n_dp = REAL(in_-1,dp)
          ! Store the kernels for sparing readings
          Knp0 =  Kernel_i(in_,ikr_,ikz_)
          Knp1 =  Kernel_i(in_+1,ikr_,ikz_)
          Knm1 =  Kernel_i(in_-1,ikr_,ikz_)
          ! Nonadiabatic moments (only different from moments when p=0)
          nadiab_moment_0j   = moments_i(1,in_  ,ikr_,ikz_,updatetlevel) + q_i_tau_i * Kernel_i(in_  ,ikr_,ikz_)*phi(ikr_,ikz_)
          ! Density
          n_     = n_     + Knp0 * nadiab_moment_0j
          ! Perpendicular velocity
          ! uperp_ = uperp_ + ibi_*0.5_dp*Kernel_i(in_,ikr_,ikz_) * (nadiab_moment_0j - nadiab_moment_0jp1)
          ! uperp_ = uperp_ + b_i*0.5_dp*Kernel_i(in_,ikr_,ikz_) * (nadiab_moment_0j - nadiab_moment_0jp1)
          uperp_ = uperp_ + bi_*0.5_dp*(Knp0 - Knm1) * nadiab_moment_0j
          ! Parallel temperature
          Tpar_  = Tpar_  + Knp0 * (SQRT2*moments_i(3,in_,ikr_,ikz_,updatetlevel) + nadiab_moment_0j)
          ! Perpendicular temperature
          ! Tperp_ = Tperp_ + Kernel_i(in_,ikr_,ikz_) * ((2._dp*n_dp+1._dp)* nadiab_moment_0j   &
          !                                            -            n_dp * nadiab_moment_0jm1 &
          !                                            -         (n_dp+1)* nadiab_moment_0jp1)
          Tperp_ = Tperp_ + ((2._dp*n_dp+1._dp)*Knp0 - (n_dp+1) * Knp1 - n_dp * Knm1)*nadiab_moment_0j
        ENDDO
        T_  = (Tpar_ + 2._dp*Tperp_)/3._dp - n_
      ! Add energy restoring term
      TColl_ = TColl_ + T_* 4._dp *  j_dp          * Kernel_i(ij_,ikr_,ikz_)
      TColl_ = TColl_ - T_* 2._dp * (j_dp + 1._dp) * Kernel_i(ij_+1,ikr_,ikz_)
      TColl_ = TColl_ - T_* 2._dp *  j_dp          * Kernel_i(ij_-1,ikr_,ikz_)
      TColl_ = TColl_ + uperp_*bi_* (Kernel_i(ij_,ikr_,ikz_) - Kernel_i(ij_-1,ikr_,ikz_))
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
      upar_  = 0._dp; uperp_ = 0._dp
      Tpar_  = 0._dp; Tperp_ = 0._dp
      DO in_ = 1,jmaxi+1
        n_dp = REAL(in_-1,dp)
        ! Store the kernels for sparing readings
        Knp0 =  Kernel_i(in_,ikr_,ikz_)
        Knp1 =  Kernel_i(in_+1,ikr_,ikz_)
        Knm1 =  Kernel_i(in_-1,ikr_,ikz_)
        ! Nonadiabatic moments (only different from moments when p=0)
        nadiab_moment_0j   = moments_i(1,in_  ,ikr_,ikz_,updatetlevel) + q_i_tau_i * Kernel_i(in_  ,ikr_,ikz_)*phi(ikr_,ikz_)
        ! Density
        n_     = n_     + Knp0 * nadiab_moment_0j
        ! Parallel temperature
        Tpar_  = Tpar_  + Knp0 * (SQRT2*moments_i(3,in_,ikr_,ikz_,updatetlevel) + nadiab_moment_0j)
        ! Perpendicular temperature
        Tperp_ = Tperp_ + ((2._dp*n_dp+1._dp)*Knp0 - (n_dp+1) * Knp1 - n_dp * Knm1)*nadiab_moment_0j
      ENDDO
      T_  = (Tpar_ + 2._dp*Tperp_)/3._dp - n_
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

    ! Execution time start
    CALL cpu_time(t0_coll)

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
            IF (num_procs_p .GT. 1) THEN
              ! Sum up all the sub collision terms on root 0
              CALL MPI_REDUCE(local_sum_e, buffer_e, pmaxe+1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, comm_p, ierr)
              ! distribute the sum over the process among p
              CALL MPI_SCATTERV(buffer_e, counts_np_e, displs_np_e, MPI_DOUBLE_COMPLEX,&
                                TColl_distr_e, local_np_e, MPI_DOUBLE_COMPLEX,&
                                0, comm_p, ierr)
            ELSE
              TColl_distr_e = local_sum_e
            ENDIF
            ! Write in output variable
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
            IF (num_procs_p .GT. 1) THEN
              ! Reduce the local_sums to root = 0
              CALL MPI_REDUCE(local_sum_i, buffer_i, pmaxi+1, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, comm_p, ierr)
              ! buffer contains the entire collision term along p, we scatter it between
              ! the other processes (use of scatterv since Pmax/Np is not an integer)
              CALL MPI_SCATTERV(buffer_i, counts_np_i, displs_np_i, MPI_DOUBLE_COMPLEX,&
                                TColl_distr_i, local_np_i, MPI_DOUBLE_COMPLEX, &
                                0, comm_p, ierr)
            ELSE
              TColl_distr_i = local_sum_i
            ENDIF
            ! Write in output variable
            DO ip = ips_i,ipe_i
              TColl_i(ip,ij,ikr,ikz) = TColl_distr_i(ip)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDIF

    ! Execution time end
    CALL cpu_time(t1_coll)
    tc_coll = tc_coll + (t1_coll - t0_coll)
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
      INTEGER :: fid ! file indexation

      INTEGER :: ip_e, ij_e, il_e, ik_e, ikps_C, ikpe_C                  ! indices for electrons loops
      REAL(dp), DIMENSION(2) :: dims_e
      INTEGER :: pdime, jdime                                            ! dimensions of the COSOlver matrices
      REAL(dp), DIMENSION(:,:),   ALLOCATABLE :: Ceepj_full, CeipjT_full ! To load the entire matrix
      REAL(dp), DIMENSION(:,:),   ALLOCATABLE :: CeipjF_full             ! ''
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: Ceepj__kp, CeipjT_kp    ! To store the coeff that will be used along kperp
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: CeipjF_kp               ! ''
      INTEGER :: ip_i, ij_i, il_i, ik_i                                  ! same for ions
      INTEGER,  DIMENSION(2) :: dims_i
      INTEGER :: pdimi, jdimi                                            ! dimensions of the COSOlver matrices
      REAL(dp), DIMENSION(:,:),   ALLOCATABLE :: Ciipj_full, CiepjT_full ! .
      REAL(dp), DIMENSION(:,:),   ALLOCATABLE :: CiepjF_full             ! .
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: Ciipj__kp, CiepjT_kp    ! .
      REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: CiepjF_kp               ! .
      INTEGER  :: NFLR

      REAL(dp), DIMENSION(:),     ALLOCATABLE :: kp_grid_mat             ! kperp grid of the matrices
      INTEGER  :: ikp_next, ikp_prev, nkp_mat, ikp_mat
      REAL(dp) ::  kp_next,  kp_prev, kperp_sim, kperp_mat, zerotoone

      CHARACTER(len=128) :: var_name, kperp_string, ikp_string

      LOGICAL     :: CO_AA_ONLY = .false. ! Flag to remove ei ie collision

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

      ! Opening the compiled cosolver matrices results
      if(my_id.EQ.0)write(*,*) mat_file
      CALL openf(mat_file,fid, 'r', 'D', mpicomm=comm_p);

      ! Get matrices dimensions (polynomials degrees and kperp grid)
      CALL getarr(fid, '/dims_e', dims_e) ! Get the electron polynomial degrees
      pdime = dims_e(1); jdime = dims_e(2);
      CALL getarr(fid, '/dims_i', dims_i) ! Get the ion      polynomial degrees
      pdimi = dims_i(1); jdimi = dims_i(2);
      IF ( ((pdime .LT. pmaxe) .OR. (jdime .LT. jmaxe)) .AND. (my_id .EQ. 0)) WRITE(*,*) '!! Pe,Je Matrix too small !!'
      IF ( ((pdimi .LT. pmaxi) .OR. (jdimi .LT. jmaxi)) .AND. (my_id .EQ. 0)) WRITE(*,*) '!! Pi,Ji Matrix too small !!'

      CALL getsize(fid, '/coordkperp', nkp_mat) ! Get the dimension kperp grid of the matrices
      CALL allocate_array(kp_grid_mat, 1,nkp_mat)
      CALL getarr(fid, '/coordkperp', kp_grid_mat)
      IF (NON_LIN) THEN ! check that we have enough kperps mat
        IF ( (kp_grid_mat(nkp_mat) .LT. two_third_kpmax) .AND. (my_id .EQ. 0)) WRITE(*,*) '!! Matrix kperp grid too small !!'
      ELSE
        IF ( (kp_grid_mat(nkp_mat) .LT. kp_max) .AND. (my_id .EQ. 0)) WRITE(*,*) '!! Matrix kperp grid too small !!'
      ENDIF

      IF (CO .GT. 0) THEN ! GK operator (k-dependant)
        ikps_C = 1; ikpe_C = nkp_mat
      ELSEIF (CO .LT. 0) THEN ! DK operator (only one mat for all k)
        ikps_C = 1; ikpe_C = 1
      ENDIF

      CALL allocate_array(  Ceepj__kp, 1,(pmaxe+1)*(jmaxe+1), 1,(pmaxe+1)*(jmaxe+1), ikps_C,ikpe_C)
      CALL allocate_array(  CeipjT_kp, 1,(pmaxe+1)*(jmaxe+1), 1,(pmaxe+1)*(jmaxe+1), ikps_C,ikpe_C)
      CALL allocate_array(  CeipjF_kp, 1,(pmaxe+1)*(jmaxe+1), 1,(pmaxe+1)*(jmaxe+1), ikps_C,ikpe_C)
      CALL allocate_array(  Ciipj__kp, 1,(pmaxe+1)*(jmaxe+1), 1,(pmaxe+1)*(jmaxe+1), ikps_C,ikpe_C)
      CALL allocate_array(  CiepjT_kp, 1,(pmaxe+1)*(jmaxe+1), 1,(pmaxe+1)*(jmaxe+1), ikps_C,ikpe_C)
      CALL allocate_array(  CiepjF_kp, 1,(pmaxe+1)*(jmaxe+1), 1,(pmaxe+1)*(jmaxe+1), ikps_C,ikpe_C)

      DO ikp = ikps_C,ikpe_C ! Loop over everz kperp values
        ! we put zeros if kp>2/3kpmax because thoses frequenvies are filtered through AA
        IF(kp_grid_mat(ikp) .GT. two_third_kpmax .AND. NON_LIN) THEN
          CiepjT_kp(:,:,ikp) = 0._dp
          CiepjF_kp(:,:,ikp) = 0._dp
          CeipjT_kp(:,:,ikp) = 0._dp
          CeipjF_kp(:,:,ikp) = 0._dp
          Ceepj__kp(:,:,ikp) = 0._dp
          Ciipj__kp(:,:,ikp) = 0._dp
        ELSE
          ! Kperp value in string format
          IF (CO .GT. 0) THEN
            write(ikp_string,'(i5.5)') ikp-1
          ELSE
            write(ikp_string,'(i5.5)') 0
          ENDIF
          !!!!!!!!!!!! Electron matrices !!!!!!!!!!!!
          ! get the self electron colision matrix
          ! Allocate space for storing full collision matrix
          CALL allocate_array(  Ceepj_full, 1,(pdime+1)*(jdime+1), 1,(pdime+1)*(jdime+1))
          ! Naming of the array to load (kperp dependant)
          WRITE(var_name,'(a,a)') TRIM(ADJUSTL(ikp_string)),'/Caapj/Ceepj'
          CALL getarr(fid, var_name, Ceepj_full) ! get array (moli format)
          ! Fill sub array with the usefull polynmial degrees only
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
          DEALLOCATE(Ceepj_full)

          ! Get test and field e-i collision matrices
          CALL allocate_array( CeipjT_full, 1,(pdime+1)*(jdime+1), 1,(pdime+1)*(jdime+1))
          CALL allocate_array( CeipjF_full, 1,(pdime+1)*(jdime+1), 1,(pdimi+1)*(jdimi+1))
          WRITE(var_name,'(a,a)') TRIM(ADJUSTL(ikp_string)),'/Ceipj/CeipjT'
          CALL getarr(fid, var_name, CeipjT_full)
          WRITE(var_name,'(a,a)') TRIM(ADJUSTL(ikp_string)),'/Ceipj/CeipjF'
          CALL getarr(fid, var_name, CeipjF_full)
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
          DEALLOCATE(CeipjF_full)
          DEALLOCATE(CeipjT_full)

          !!!!!!!!!!!!!!! Ion matrices !!!!!!!!!!!!!!
          ! get the self electron colision matrix
          CALL allocate_array(  Ciipj_full, 1,(pdimi+1)*(jdimi+1), 1,(pdimi+1)*(jdimi+1))
          WRITE(var_name,'(a,a,a)') TRIM(ADJUSTL(ikp_string)),'/Caapj/Ciipj'
          CALL getarr(fid, var_name, Ciipj_full) ! get array (moli format)
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
          DEALLOCATE(Ciipj_full)

          ! get the Test and Back field electron ion collision matrix
          CALL allocate_array( CiepjT_full, 1,(pdimi+1)*(jdimi+1), 1,(pdimi+1)*(jdimi+1))
          CALL allocate_array( CiepjF_full, 1,(pdimi+1)*(jdimi+1), 1,(pdime+1)*(jdime+1))
          WRITE(var_name,'(a,a,a)') TRIM(ADJUSTL(ikp_string)),'/Ciepj/CiepjT'
          CALL getarr(fid, var_name, CiepjT_full)
          WRITE(var_name,'(a,a,a)') TRIM(ADJUSTL(ikp_string)),'/Ciepj/CiepjF'
          CALL getarr(fid, var_name, CiepjF_full)
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
          DEALLOCATE(CiepjF_full)
          DEALLOCATE(CiepjT_full)
        ENDIF
      ENDDO
      CALL closef(fid)

      IF (CO .GT. 0) THEN ! Interpolation of the kperp matrix values on kr kz grid
        IF (my_id .EQ. 0 ) WRITE(*,*) '...Interpolation from matrices kperp to simulation kr,kz...'
        DO ikr = ikrs,ikre
          DO ikz = ikzs,ikze
            kperp_sim = SQRT(krarray(ikr)**2+kzarray(ikz)**2) ! current simulation kperp

            ! Find the interval in kp grid mat where kperp_sim is contained
            ! Loop over the whole kp mat grid to find the smallest kperp that is
            ! larger than the current kperp_sim (brute force...)
            DO ikp=1,nkp_mat
              ikp_mat   = ikp ! the first indice of the interval (k0)
              kperp_mat = kp_grid_mat(ikp)
              IF(kperp_mat .GT. kperp_sim) EXIT
            ENDDO
            ! interval boundaries
            ikp_prev  = ikp_mat - 1 !index of k0
            ikp_next  = ikp_mat     !index of k1
            ! write(*,*) kp_grid_mat(ikp_prev), '<', kperp_sim, '<', kp_grid_mat(ikp_next)

            ! 0->1 variable for linear interp, i.e. zero2one = (k-k0)/(k1-k0)
            zerotoone = (kperp_sim - kp_grid_mat(ikp_prev))/(kp_grid_mat(ikp_next) - kp_grid_mat(ikp_prev))

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

      IF( CO_AA_ONLY ) THEN
        CeipjF = 0._dp;
        CeipjT = 0._dp;
        CiepjF = 0._dp;
        CiepjT = 0._dp;
      ENDIF

      IF (my_id .EQ. 0) WRITE(*,*) '============DONE==========='

    END SUBROUTINE load_COSOlver_mat
    !******************************************************************************!

end module collision
