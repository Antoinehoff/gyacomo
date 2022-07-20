!******************************************************************************!
!!!!!! initialize the moments and load/build coeff tables
!******************************************************************************!
SUBROUTINE inital

  USE basic,            ONLY: my_id
  USE initial_par,      ONLY: INIT_OPT
  USE time_integration, ONLY: set_updatetlevel
  USE collision,        ONLY: load_COSOlver_mat, cosolver_coll
  USE closure,          ONLY: apply_closure_model
  USE ghosts,           ONLY: update_ghosts_moments, update_ghosts_EM
  USE restarts,         ONLY: load_moments, job2load
  USE numerics,         ONLY: play_with_modes, save_EM_ZF_modes
  USE processing,       ONLY: compute_fluid_moments
  USE model,            ONLY: KIN_E, LINEARITY
  USE nonlinear,        ONLY: compute_Sapj, nonlinear_init
  IMPLICIT NONE

  CALL set_updatetlevel(1)

  !!!!!! Set the moments arrays Nepj, Nipj and phi!!!!!!
  ! through loading a previous state
  IF ( job2load .GE. 0 ) THEN
    IF (my_id .EQ. 0) WRITE(*,*) 'Load moments'
    CALL load_moments ! get N_0
    CALL update_ghosts_moments
    CALL solve_EM_fields ! compute phi_0=phi(N_0)
    CALL update_ghosts_EM
  ! through initialization
  ELSE
    SELECT CASE (INIT_OPT)
    ! set phi with noise and set moments to 0
    CASE ('phi')
      IF (my_id .EQ. 0) WRITE(*,*) 'Init noisy phi'
      CALL init_phi
      CALL update_ghosts_EM
    ! set moments_00 (GC density) with noise and compute phi afterwards
    CASE('mom00')
      IF (my_id .EQ. 0) WRITE(*,*) 'Init noisy gyrocenter density'
      CALL init_gyrodens ! init only gyrocenter density
      CALL update_ghosts_moments
      CALL solve_EM_fields
      CALL update_ghosts_EM
    ! init all moments randomly (unadvised)
    CASE('allmom')
      IF (my_id .EQ. 0) WRITE(*,*) 'Init noisy moments'
      CALL init_moments ! init all moments
      CALL update_ghosts_moments
      CALL solve_EM_fields
      CALL update_ghosts_EM
    ! init a gaussian blob in gyrodens
    CASE('blob')
      IF (my_id .EQ. 0) WRITE(*,*) '--init a blob'
      CALL initialize_blob
      CALL update_ghosts_moments
      CALL solve_EM_fields
      CALL update_ghosts_EM
    ! init moments 00 with a power law similarly to GENE
    CASE('ppj')
      IF (my_id .EQ. 0) WRITE(*,*) 'ppj init ~ GENE'
      call init_ppj
      CALL update_ghosts_moments
      CALL solve_EM_fields
      CALL update_ghosts_EM
    END SELECT
  ENDIF
  ! closure of j>J, p>P and j<0, p<0 moments
  IF (my_id .EQ. 0) WRITE(*,*) 'Apply closure'
  CALL apply_closure_model
  ! ghosts for p parallelization
  IF (my_id .EQ. 0) WRITE(*,*) 'Ghosts communication'
  CALL update_ghosts_moments
  CALL update_ghosts_EM
  !! End of phi and moments initialization

  ! Save (kx,0) and (0,ky) modes for num exp
  CALL save_EM_ZF_modes
  ! Freeze/Wipe some selected modes (entropy,zonal,turbulent)
  CALL play_with_modes

  ! Load the COSOlver collision operator coefficients
  IF(cosolver_coll) &
  CALL load_COSOlver_mat

  !! Preparing auxiliary arrays at initial state
  ! particle density, fluid velocity and temperature (used in diagnose)
  IF (my_id .EQ. 0) WRITE(*,*) 'Computing fluid moments'
  CALL compute_fluid_moments

  ! init auxval for nonlinear
  CALL nonlinear_init
  ! compute nonlinear for t=0 diagnostic
  CALL compute_Sapj ! compute S_0 = S(phi_0,N_0)

END SUBROUTINE inital
!******************************************************************************!

!******************************************************************************!
!!!!!!! Initialize all the moments randomly
!******************************************************************************!
SUBROUTINE init_moments
  USE basic
  USE grid
  USE initial_par,ONLY: iseed, init_noiselvl, init_background
  USE fields,     ONLY: moments_e, moments_i
  USE prec_const, ONLY: dp
  USE utility,    ONLY: checkfield
  USE model,      ONLY : LINEARITY, KIN_E
  IMPLICIT NONE

  REAL(dp) :: noise
  INTEGER, DIMENSION(12) :: iseedarr

  ! Seed random number generator
  iseedarr(:)=iseed
  CALL RANDOM_SEED(PUT=iseedarr+my_id)

    !**** Broad noise initialization *******************************************
    ! Electron init
    IF(KIN_E) THEN
    DO ip=ips_e,ipe_e
      DO ij=ijs_e,ije_e
        DO ikx=ikxs,ikxe
          DO iky=ikys,ikye
            DO iz=izs,ize
              CALL RANDOM_NUMBER(noise)
              moments_e(ip,ij,iky,ikx,iz,:) = (init_background + init_noiselvl*(noise-0.5_dp))
            END DO
          END DO
        END DO
        IF ( contains_ky0 ) THEN
          DO iky=2,Nky/2 !symmetry at ky = 0 for all z
            moments_e(ip,ij,iky_0,ikx,:,:) = moments_e( ip,ij,iky_0,Nkx+2-ikx,:, :)
          END DO
        ENDIF
      END DO
    END DO
    ENDIF
    ! Ion init
    DO ip=ips_i,ipe_i
      DO ij=ijs_i,ije_i

        DO ikx=ikxs,ikxe
          DO iky=ikys,ikye
            DO iz=izs,ize
              CALL RANDOM_NUMBER(noise)
              moments_i(ip,ij,iky,ikx,iz,:) = (init_background + init_noiselvl*(noise-0.5_dp))
            END DO
          END DO
        END DO

        IF ( contains_ky0 ) THEN
          DO ikx=2,Nkx/2 !symmetry at ky = 0 for all z
            moments_i( ip,ij,iky_0,ikx,:,:) = moments_i( ip,ij,iky_0,Nkx+2-ikx,:,:)
          END DO
        ENDIF

      END DO
    END DO

    ! Putting to zero modes that are not in the 2/3 Orszag rule
    IF (LINEARITY .NE. 'linear') THEN
      DO ikx=ikxs,ikxe
      DO iky=ikys,ikye
      DO iz=izs,ize
        IF(KIN_E) THEN
        DO ip=ips_e,ipe_e
        DO ij=ijs_e,ije_e
          moments_e( ip,ij,iky,ikx,iz, :) = moments_e( ip,ij,iky,ikx,iz, :)*AA_x(ikx)*AA_y(iky)
        ENDDO
        ENDDO
        ENDIF
        DO ip=ips_i,ipe_i
        DO ij=ijs_i,ije_i
          moments_i( ip,ij,iky,ikx,iz, :) = moments_i( ip,ij,iky,ikx,iz, :)*AA_x(ikx)*AA_y(iky)
        ENDDO
        ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDIF
END SUBROUTINE init_moments
!******************************************************************************!

!******************************************************************************!
!!!!!!! Initialize the gyrocenter density randomly
!******************************************************************************!
SUBROUTINE init_gyrodens
  USE basic
  USE grid
  USE fields
  USE prec_const
  USE utility, ONLY: checkfield
  USE initial_par
  USE model, ONLY: KIN_E, LINEARITY
  IMPLICIT NONE

  REAL(dp) :: noise
  REAL(dp) :: kx, ky, sigma, gain, ky_shift
  INTEGER, DIMENSION(12) :: iseedarr

  ! Seed random number generator
  iseedarr(:)=iseed
  CALL RANDOM_SEED(PUT=iseedarr+my_id)

    !**** Broad noise initialization *******************************************
    IF(KIN_E) THEN
    DO ip=ips_e,ipe_e
      DO ij=ijs_e,ije_e
        DO ikx=ikxs,ikxe
          DO iky=ikys,ikye
            DO iz=izs,ize
              CALL RANDOM_NUMBER(noise)
              IF ( (ip .EQ. 1) .AND. (ij .EQ. 1) ) THEN
                moments_e(ip,ij,iky,ikx,iz,:) = (init_background + init_noiselvl*(noise-0.5_dp))
              ELSE
                moments_e(ip,ij,iky,ikx,iz,:) = 0._dp
              ENDIF
            END DO
          END DO
        END DO
        IF ( contains_ky0 ) THEN
          DO ikx=2,Nkx/2 !symmetry at ky = 0 for all z
            moments_e(ip,ij,iky_0,ikx,:,:) = moments_e( ip,ij,iky_0,Nkx+2-ikx,:, :)
          END DO
        ENDIF
      END DO
    END DO
    ENDIF

    DO ip=ips_i,ipe_i
      DO ij=ijs_i,ije_i
        DO ikx=ikxs,ikxe
          DO iky=ikys,ikye
            DO iz=izs,ize
              CALL RANDOM_NUMBER(noise)
              IF ( (ip .EQ. 1) .AND. (ij .EQ. 1) ) THEN
                moments_i(ip,ij,iky,ikx,iz,:) = (init_background + init_noiselvl*(noise-0.5_dp))
              ELSE
                moments_i(ip,ij,iky,ikx,iz,:) = 0._dp
              ENDIF
            END DO
          END DO
        END DO
        IF ( contains_ky0 ) THEN
          DO iky=2,Nky/2 !symmetry at ky = 0 for all z
            moments_i( ip,ij,iky_0,ikx,:,:) = moments_i( ip,ij,iky_0,Nkx+2-ikx,:,:)
          END DO
        ENDIF
      END DO
    END DO

    ! Putting to zero modes that are not in the 2/3 Orszag rule
    IF (LINEARITY .NE. 'linear') THEN
      DO ikx=ikxs,ikxe
      DO iky=ikys,ikye
      DO iz=izs,ize
        IF(KIN_E) THEN
        DO ip=ips_e,ipe_e
        DO ij=ijs_e,ije_e
          moments_e( ip,ij,iky,ikx,iz, :) = moments_e( ip,ij,iky,ikx,iz, :)*AA_x(ikx)*AA_y(iky)
        ENDDO
        ENDDO
        ENDIF
        DO ip=ips_i,ipe_i
        DO ij=ijs_i,ije_i
          moments_i( ip,ij,iky,ikx,iz, :) = moments_i( ip,ij,iky,ikx,iz, :)*AA_x(ikx)*AA_y(iky)
        ENDDO
        ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDIF
END SUBROUTINE init_gyrodens
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
  USE model, ONLY: KIN_E, LINEARITY
  IMPLICIT NONE

  REAL(dp) :: noise
  REAL(dp) :: kx, ky, kp, sigma, gain, ky_shift
  INTEGER, DIMENSION(12) :: iseedarr

  ! Seed random number generator
  iseedarr(:)=iseed
  CALL RANDOM_SEED(PUT=iseedarr+my_id)

    !**** noise initialization *******************************************

    DO ikx=ikxs,ikxe
      DO iky=ikys,ikye
        DO iz=izs,ize
          CALL RANDOM_NUMBER(noise)
          phi(iky,ikx,iz) = (init_background + init_noiselvl*(noise-0.5_dp))!*AA_x(ikx)*AA_y(iky)
        ENDDO
      END DO
    END DO

    !symmetry at ky = 0 to keep real inverse transform
    IF ( contains_ky0 ) THEN
      DO ikx=2,Nkx/2
        phi(iky_0,ikx,izs:ize) = phi(iky_0,Nkx+2-ikx,izs:ize)
      END DO
      phi(iky_0,ikx_0,izs:ize) = REAL(phi(iky_0,ikx_0,izs:ize)) !origin must be real
    ENDIF

    !**** ensure no previous moments initialization
    IF(KIN_E) moments_e = 0._dp
    moments_i = 0._dp

    !**** Zonal Flow initialization *******************************************
    ! put a mode at ikx = mode number + 1, symmetry is already included since kx>=0
    IF(INIT_ZF .GT. 0) THEN
      IF (my_id .EQ. 0) WRITE(*,*) 'Init ZF phi'
      IF( (INIT_ZF+1 .GT. ikxs) .AND. (INIT_ZF+1 .LT. ikxe) ) THEN
        DO iz = izs,ize
          phi(iky_0,INIT_ZF+1,iz) = ZF_AMP*(2._dp*PI)**2/deltakx/deltaky/2._dp * COS((iz-1)/Nz*2._dp*PI)
          moments_i(1,1,iky_0,INIT_ZF+1,iz,:) = kxarray(INIT_ZF+1)**2*phi(iky_0,INIT_ZF+1,iz)* COS((iz-1)/Nz*2._dp*PI)
          IF(KIN_E) moments_e(1,1,iky_0,INIT_ZF+1,iz,:) = 0._dp
        ENDDO
      ENDIF
    ENDIF

END SUBROUTINE init_phi
!******************************************************************************!

!******************************************************************************!
!******************************************************************************!
!!!!!!! Initialize an ionic Gaussian blob on top of the preexisting modes
!******************************************************************************!
SUBROUTINE initialize_blob
  USE fields
  USE grid
  USE model, ONLY: sigmai2_taui_o2, KIN_E, LINEARITY
  IMPLICIT NONE
  REAL(dp) ::kx, ky, sigma, gain
  sigma = 0.5_dp
  gain  = 5e2_dp

  DO ikx=ikxs,ikxe
    kx = kxarray(ikx)
  DO iky=ikys,ikye
    ky = kyarray(iky)
  DO iz=izs,ize
    DO ip=ips_i,ipe_i
    DO ij=ijs_i,ije_i
      IF( (iky .NE. iky_0) .AND. (ip .EQ. 1) .AND. (ij .EQ. 1)) THEN
        moments_i( ip,ij,iky,ikx,iz, :) = moments_i( ip,ij,iky,ikx,iz, :) &
        + gain*sigma/SQRT2 * exp(-(kx**2+ky**2)*sigma**2/4._dp) &
          * AA_x(ikx)*AA_y(iky)!&
          ! * exp(sigmai2_taui_o2*(kx**2+ky**2))
      ENDIF
    ENDDO
    ENDDO
    IF(KIN_E) THEN
    DO ip=ips_e,ipe_e
    DO ij=ijs_e,ije_e
      IF( (iky .NE. iky_0) .AND. (ip .EQ. 1) .AND. (ij .EQ. 1)) THEN
        moments_e( ip,ij,iky,ikx,iz, :) = moments_e( ip,ij,iky,ikx,iz, :) &
        + gain*sigma/SQRT2 * exp(-(kx**2+ky**2)*sigma**2/4._dp) &
          * AA_x(ikx)*AA_y(iky)!&
          ! * exp(sigmai2_taui_o2*(kx**2+ky**2))
      ENDIF
    ENDDO
    ENDDO
    ENDIF
  ENDDO
  ENDDO
  ENDDO
END SUBROUTINE initialize_blob
!******************************************************************************!


!******************************************************************************!
!!!!!!! Initialize the gyrocenter in a ppj gene manner (power law)
!******************************************************************************!
SUBROUTINE init_ppj
  USE basic
  USE grid
  USE fields
  USE prec_const
  USE utility, ONLY: checkfield
  USE initial_par
  USE model, ONLY: KIN_E, LINEARITY
  USE geometry, ONLY: Jacobian, iInt_Jacobian

  IMPLICIT NONE

  REAL(dp) :: noise
  REAL(dp) :: kx, ky, sigma_z, amp, ky_shift, z
  INTEGER, DIMENSION(12) :: iseedarr

  sigma_z = pi/4._dp
  amp = 1.0_dp

    !**** Broad noise initialization *******************************************
    ! Electrons
    IF (KIN_E) THEN
    DO ip=ips_e,ipe_e
      DO ij=ijs_e,ije_e
        IF ( (ip .EQ. 1) .AND. (ij .EQ. 1) ) THEN
          DO ikx=ikxs,ikxe
            kx = kxarray(ikx)
            DO iky=ikys,ikye
              ky = kyarray(iky)
              DO iz=izs,ize
                z = zarray(iz,0)
                IF (kx .EQ. 0) THEN
                  IF(ky .EQ. 0) THEN
                    moments_e(ip,ij,iky,ikx,iz,:) = 0._dp
                  ELSE
                    moments_e(ip,ij,iky,ikx,iz,:) = 0.5_dp * ky_min/(ABS(ky)+ky_min)
                  ENDIF
                ELSE
                  IF(ky .GT. 0) THEN
                    moments_e(ip,ij,iky,ikx,iz,:) = (kx_min/(ABS(kx)+kx_min))*(ky_min/(ABS(ky)+ky_min))
                  ELSE
                    moments_e(ip,ij,iky,ikx,iz,:) = 0.5_dp*amp*(kx_min/(ABS(kx)+kx_min))
                  ENDIF
                ENDIF
                CALL RANDOM_NUMBER(noise)
                ! z-dep and noise
                moments_e(ip,ij,iky,ikx,iz,:) = moments_e(ip,ij,iky,ikx,iz,:) * &
                (Jacobian(iz,0)*iInt_Jacobian)**2 !+ (init_background + init_noiselvl*(noise-0.5_dp))
              END DO
            END DO
          END DO

          IF ( contains_ky0 ) THEN
            DO ikx=2,Nkx/2 !symmetry at kx = 0 for all z
              moments_e(ip,ij,iky_0,ikx,:,:) = moments_e( ip,ij,iky_0,Nkx+2-ikx,:, :)
            END DO
          ENDIF
        ELSE
          moments_e(ip,ij,:,:,:,:) = 0._dp
        ENDIF
      END DO
    END DO
    ENDIF

    ! Ions
    DO ip=ips_i,ipe_i
      DO ij=ijs_i,ije_i
        IF ( (ip .EQ. 1) .AND. (ij .EQ. 1) ) THEN
          DO ikx=ikxs,ikxe
            kx = kxarray(ikx)
            DO iky=ikys,ikye
              ky = kyarray(iky)
              DO iz=izs,ize
                z = zarray(iz,0)
                IF (kx .EQ. 0) THEN
                  IF(ky .EQ. 0) THEN
                    moments_i(ip,ij,iky,ikx,iz,:) = 0._dp
                  ELSE
                    moments_i(ip,ij,iky,ikx,iz,:) = 0.5_dp * ky_min/(ABS(ky)+ky_min)
                  ENDIF
                ELSE
                  IF(ky .GT. 0) THEN
                    moments_i(ip,ij,iky,ikx,iz,:) = (kx_min/(ABS(kx)+kx_min))*(ky_min/(ABS(ky)+ky_min))
                  ELSE
                    moments_i(ip,ij,iky,ikx,iz,:) = 0.5_dp*amp*(kx_min/(ABS(kx)+kx_min))
                  ENDIF
                ENDIF
                CALL RANDOM_NUMBER(noise)
                ! z-dep and noise
                moments_i(ip,ij,iky,ikx,iz,:) = moments_i(ip,ij,iky,ikx,iz,:) * &
                (Jacobian(iz,0)*iInt_Jacobian)**2 !+ (init_background + init_noiselvl*(noise-0.5_dp))
              END DO
            END DO
          END DO

          IF ( contains_ky0 ) THEN
            DO ikx=2,Nkx/2 !symmetry at kx = 0 for all z
              moments_i(ip,ij,iky_0,ikx,:,:) = moments_i( ip,ij,iky_0,Nkx+2-ikx,:, :)
            END DO
          ENDIF
        ELSE
          moments_i(ip,ij,:,:,:,:) = 0._dp
        ENDIF
      END DO
    END DO

    ! Putting to zero modes that are not in the 2/3 Orszag rule
    IF (LINEARITY .NE. 'linear') THEN
      DO ikx=ikxs,ikxe
      DO iky=ikys,ikye
      DO iz=izs,ize
        DO ip=ips_e,ipe_e
        DO ij=ijs_e,ije_e
          moments_e( ip,ij,iky,ikx,iz, :) = moments_e( ip,ij,iky,ikx,iz, :)*AA_x(ikx)*AA_y(iky)
        ENDDO
        ENDDO
        DO ip=ips_i,ipe_i
        DO ij=ijs_i,ije_i
          moments_i( ip,ij,iky,ikx,iz, :) = moments_i( ip,ij,iky,ikx,iz, :)*AA_x(ikx)*AA_y(iky)
        ENDDO
        ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDIF

END SUBROUTINE init_ppj
!******************************************************************************!
