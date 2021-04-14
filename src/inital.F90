!******************************************************************************!
!!!!!! initialize the moments and load/build coeff tables
!******************************************************************************!
SUBROUTINE inital

  USE basic
  USE model, ONLY : CO, NON_LIN
  USE initial_par
  USE prec_const
  USE coeff
  USE time_integration
  USE array, ONLY : Sepj,Sipj
  USE collision
  USE closure
  USE ghosts
  USE restarts

  implicit none

  CALL set_updatetlevel(1)

  IF (my_id .EQ. 0) WRITE(*,*) 'Evaluate kernels'
  CALL evaluate_kernels

  !!!!!! Set the moments arrays Nepj, Nipj and phi!!!!!!
  IF ( RESTART ) THEN
    IF (my_id .EQ. 0) WRITE(*,*) 'Load moments'
    CALL load_moments ! get N_0

    IF (my_id .EQ. 0) WRITE(*,*) 'Init phi with Poisson'
    CALL poisson ! compute phi_0=phi(N_0)

  ELSE
    IF (INIT_NOISY_PHI) THEN
      IF (my_id .EQ. 0) WRITE(*,*) 'Init noisy phi'
      CALL init_phi ! init noisy phi_0, N_0 = 0
    ELSE
      IF (my_id .EQ. 0) WRITE(*,*) 'Init noisy moments and ghosts'
      CALL init_moments ! init noisy N_0
      IF (my_id .EQ. 0) WRITE(*,*) 'Init phi with Poisson'
      CALL poisson ! get phi_0 = phi(N_0)
    ENDIF

  ENDIF

  IF (my_id .EQ. 0) WRITE(*,*) 'Apply closure'
  CALL apply_closure_model

  IF (my_id .EQ. 0) WRITE(*,*) 'Ghosts communication'
  CALL update_ghosts

  !!!!!! Set Sepj, Sipj and dnjs coeff table !!!!!!
  IF ( NON_LIN ) THEN;
    IF (my_id .EQ. 0) WRITE(*,*) 'Building Dnjs table'
    CALL build_dnjs_table

    IF (my_id .EQ. 0) WRITE(*,*) 'Init Sapj'
    CALL compute_Sapj ! compute S_0 = S(phi_0,N_0)
  ENDIF

  !!!!!! Load the COSOlver collision operator coefficients !!!!!!
  IF (ABS(CO) .GT. 1) THEN
    CALL load_COSOlver_mat
    ! Compute collision
    CALL compute_TColl ! compute C_0 = C(N_0)
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
  USE model, ONLY : NON_LIN
  IMPLICIT NONE

  REAL(dp) :: noise
  REAL(dp) :: kr, kz, sigma, gain, kz_shift
  INTEGER, DIMENSION(12) :: iseedarr

  ! Seed random number generator
  iseedarr(:)=iseed
  CALL RANDOM_SEED(PUT=iseedarr+my_id)

    !**** Broad noise initialization *******************************************
    DO ip=ips_e,ipe_e
      DO ij=ijs_e,ije_e

        DO ikr=ikrs,ikre
          DO ikz=ikzs,ikze
            CALL RANDOM_NUMBER(noise)
            moments_e( ip,ij, ikr,ikz, :) = (initback_moments + initnoise_moments*(noise-0.5_dp))
          END DO
        END DO

        IF ( contains_kr0 ) THEN
          DO ikz=2,Nkz/2 !symmetry at kr = 0
            moments_e( ip,ij,ikr_0,ikz, :) = moments_e( ip,ij,ikr_0,Nkz+2-ikz, :)
          END DO
        ENDIF

      END DO
    END DO

    DO ip=ips_i,ipe_i
      DO ij=ijs_i,ije_i

        DO ikr=ikrs,ikre
          DO ikz=ikzs,ikze
            CALL RANDOM_NUMBER(noise)
            moments_i( ip,ij, ikr,ikz, :) = (initback_moments + initnoise_moments*(noise-0.5_dp))
          END DO
        END DO

        IF ( contains_kr0 ) THEN
          DO ikz=2,Nkz/2 !symmetry at kr = 0
            moments_i( ip,ij,ikr_0,ikz, :) = moments_i( ip,ij,ikr_0,Nkz+2-ikz, :)
          END DO
        ENDIF

      END DO
    END DO

    ! Putting to zero modes that are not in the 2/3 Orszag rule
    IF (NON_LIN) THEN
      DO ip=ips_e,ipe_e
      DO ij=ijs_e,ije_e
      DO ikr=ikrs,ikre
      DO ikz=ikzs,ikze
        moments_e( ip,ij,ikr,ikz, :) = moments_e( ip,ij,ikr,ikz, :)*AA_r(ikr)*AA_z(ikz)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      DO ip=ips_i,ipe_i
      DO ij=ijs_i,ije_i
      DO ikr=ikrs,ikre
      DO ikz=ikzs,ikze
        moments_i( ip,ij,ikr,ikz, :) = moments_i( ip,ij,ikr,ikz, :)*AA_r(ikr)*AA_z(ikz)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDIF
END SUBROUTINE init_moments
!******************************************************************************!


!******************************************************************************!
!!!!!!! Initialize a noisy ES potential and cancel the moments
!******************************************************************************!
SUBROUTINE init_phi
  USE basic
  USE grid
  USE fields
  USE prec_const
  USE initial_par
  IMPLICIT NONE

  REAL(dp) :: noise
  REAL(dp) :: kr, kz, sigma, gain, kz_shift
  INTEGER, DIMENSION(12) :: iseedarr

  IF (INIT_NOISY_PHI) THEN

    ! Seed random number generator
    iseedarr(:)=iseed
    CALL RANDOM_SEED(PUT=iseedarr+my_id)

      !**** noise initialization *******************************************

      DO ikr=ikrs,ikre
        DO ikz=ikzs,ikze
          CALL RANDOM_NUMBER(noise)
          phi(ikr,ikz) = (initback_moments + initnoise_moments*(noise-0.5_dp))*AA_r(ikr)*AA_z(ikz)
        END DO
      END DO

      !symmetry at kr = 0 to keep real inverse transform
      IF ( contains_kr0 ) THEN
        DO ikz=2,Nkz/2
          phi(ikr_0,ikz) = phi(ikr_0,Nkz+2-ikz)
        END DO
        phi(ikr_0,Nz/2) = REAL(phi(ikr_0,Nz/2)) !origin must be real
      ENDIF

      !**** Cancel previous moments initialization
      DO ip=ips_e,ipe_e
        DO ij=ijs_e,ije_e
          DO ikr=ikrs,ikre
            DO ikz=ikzs,ikze
              moments_e( ip,ij,ikr,ikz, :) = 0._dp
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      DO ip=ips_i,ipe_i
        DO ij=ijs_i,ije_i
          DO ikr=ikrs,ikre
            DO ikz=ikzs,ikze
              moments_i( ip,ij,ikr,ikz, :) = 0._dp
            ENDDO
          ENDDO
        ENDDO
      ENDDO

    ELSE ! we compute phi from noisy moments and poisson

      CALL poisson
    ENDIF

END SUBROUTINE init_phi
!******************************************************************************!

!******************************************************************************!
!!!!!!! Build the Laguerre-Laguerre coupling coefficient table for nonlin
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
END SUBROUTINE build_dnjs_table
!******************************************************************************!

!******************************************************************************!
!!!!!!! Evaluate the kernels once for all
!******************************************************************************!
SUBROUTINE evaluate_kernels
  USE basic
  USE array, Only : kernel_e, kernel_i
  USE grid
  use model, ONLY : tau_e, tau_i, sigma_e, sigma_i, q_e, q_i, lambdaD, CLOS
  IMPLICIT NONE

  REAL(dp)    :: factj, j_dp, j_int
  REAL(dp)    :: sigmae2_taue_o2, sigmai2_taui_o2
  REAL(dp)    :: be_2, bi_2, alphaD
  REAL(dp)    :: kr, kz, kperp2

  !!!!! Electron kernels !!!!!
  !Precompute species dependant factors
  sigmae2_taue_o2 = sigma_e**2 * tau_e/2._dp ! factor of the Kernel argument

  factj = 1.0 ! Start of the recursive factorial
  DO ij = 1, jmaxe+1
    j_int = jarray_e(ij)
    j_dp = REAL(j_int,dp) ! REAL of degree

    ! Recursive factorial
    IF (j_dp .GT. 0) THEN
      factj = factj * j_dp
    ELSE
      factj = 1._dp
    ENDIF

    DO ikr = ikrs,ikre
      kr     = krarray(ikr)
      DO ikz = ikzs,ikze
        kz    = kzarray(ikz)

        be_2  =  (kr**2 + kz**2) * sigmae2_taue_o2

        kernel_e(ij, ikr, ikz) = be_2**j_int * exp(-be_2)/factj

      ENDDO
    ENDDO
  ENDDO
  ! Kernels closure
  DO ikr = ikrs,ikre
    kr     = krarray(ikr)
    DO ikz = ikzs,ikze
      kz    = kzarray(ikz)
      be_2  =  (kr**2 + kz**2) * sigmae2_taue_o2
      kernel_e(ijeg_e,ikr,ikz) = be_2/(real(ijeg_e,dp))*kernel_e(ije_e,ikr,ikz)
    ENDDO
  ENDDO

  !!!!! Ion kernels !!!!!
  sigmai2_taui_o2 = sigma_i**2 * tau_i/2._dp ! (b_a/2)^2 = (kperp sqrt(2 tau_a) sigma_a/2)^2

  factj = 1.0 ! Start of the recursive factorial
  DO ij = 1, jmaxi+1
    j_int = jarray_e(ij)
    j_dp = REAL(j_int,dp) ! REAL of degree

    ! Recursive factorial
    IF (j_dp .GT. 0) THEN
      factj = factj * j_dp
    ELSE
      factj = 1._dp
    ENDIF

    DO ikr = ikrs,ikre
      kr     = krarray(ikr)
      DO ikz = ikzs,ikze
        kz    = kzarray(ikz)

        bi_2  =  (kr**2 + kz**2) * sigmai2_taui_o2

        kernel_i(ij, ikr, ikz) = bi_2**j_int * exp(-bi_2)/factj

      ENDDO
    ENDDO
  ENDDO
  ! Kernels closure
  DO ikr = ikrs,ikre
    kr     = krarray(ikr)
    DO ikz = ikzs,ikze
      kz    = kzarray(ikz)
      bi_2  =  (kr**2 + kz**2) * sigmai2_taui_o2
      kernel_i(ijeg_i,ikr,ikz) = bi_2/(real(ijeg_i,dp))*kernel_e(ije_i,ikr,ikz)
    ENDDO
  ENDDO
END SUBROUTINE evaluate_kernels
!******************************************************************************!
