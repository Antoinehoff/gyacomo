MODULE initial
  !   Module for initial parameters

  USE prec_const
  IMPLICIT NONE
  PRIVATE

  ! Initialization option (phi/mom00/allmom/blob)
  CHARACTER(len=32), PUBLIC, PROTECTED :: INIT_OPT = 'phi'
  ! Act on modes artificially (keep/wipe, zonal, non zonal, entropy mode etc.)
  ! REMOVED FEATURE (only here for retrocompatibility)
  CHARACTER(len=32),  PUBLIC, PROTECTED :: ACT_ON_MODES = 'nothing'
  ! Initial amplitude (for specific modes or image init.)
  REAL(xp), PUBLIC, PROTECTED :: init_amp        = 1._xp
  ! Initial background level (in addition to noise init)
  REAL(xp), PUBLIC, PROTECTED :: init_background = 0._xp
  ! Initial noise amplitude (for noise init)
  REAL(xp), PUBLIC, PROTECTED :: init_noiselvl   = 1E-6_xp
  ! Initialization for random number generator
  INTEGER,  PUBLIC, PROTECTED :: iseed=42

  PUBLIC :: initial_outputinputs, initial_readinputs, initialize

CONTAINS


  SUBROUTINE initial_readinputs
    ! Read the input parameters
    USE basic, ONLY : lu_in
    USE prec_const
    IMPLICIT NONE

    NAMELIST /INITIAL/ INIT_OPT,ACT_ON_MODES,&
                       init_amp,init_background,init_noiselvl,iseed
    READ(lu_in,initial)
  END SUBROUTINE initial_readinputs

  !******************************************************************************!
  !!!!!! initialize the moments and load/build coeff tables
  !******************************************************************************!
  SUBROUTINE initialize
    USE basic,            ONLY: speak
    USE time_integration, ONLY: set_updatetlevel
    USE collision,        ONLY: init_collision
    USE closure,          ONLY: apply_closure_model
    USE ghosts,           ONLY: update_ghosts_moments, update_ghosts_EM
    USE restarts,         ONLY: load_moments, job2load
    USE processing,       ONLY: compute_fluid_moments, compute_radial_heatflux, compute_radial_transport
    USE nonlinear,        ONLY: compute_Sapj, nonlinear_init
    IMPLICIT NONE
    CALL set_updatetlevel(1)
    !!!!!! Set the moments arrays Nepj, Nipj and phi!!!!!!
    ! through loading a previous state
    IF ( job2load .GE. 0 ) THEN
      CALL speak('Load moments')
      CALL load_moments ! get N_0
      CALL apply_closure_model
      CALL update_ghosts_moments
      CALL solve_EM_fields ! compute phi_0=phi(N_0)
      CALL update_ghosts_EM
    ! through initialization
    ELSE
      SELECT CASE (INIT_OPT)
      ! set phi with noise and set moments to 0
      CASE ('phi')
        CALL speak('Init noisy phi')
        CALL init_phi
        CALL update_ghosts_EM
      CASE ('phi_ppj')
        CALL speak('Init noisy phi')
        CALL init_phi_ppj
        CALL update_ghosts_EM
      ! set moments_00 (GC density) with noise and compute phi afterwards
      CASE('mom00')
        CALL speak('Init noisy gyrocenter density')
        CALL init_gyrodens ! init only gyrocenter density
        CALL update_ghosts_moments
        CALL solve_EM_fields
        CALL update_ghosts_EM
      ! set moments_00 (GC density) with a single kx0,ky0 mode
      ! kx0 is setup but by the noise value and ky0 by the background value
      CASE('mom00_single_mode')
        CALL speak('Init noisy gyrocenter density')
        CALL init_single_mode ! init single mode in gyrocenter density
        CALL update_ghosts_moments
        CALL solve_EM_fields
        CALL update_ghosts_EM
      ! init all moments randomly (unadvised)
      CASE('allmom')
        CALL speak('Init noisy moments')
        CALL init_moments ! init all moments
        CALL update_ghosts_moments
        CALL solve_EM_fields
        CALL update_ghosts_EM
      ! init a gaussian blob in gyrodens
      CASE('blob')
        CALL speak('--init a blob')
        CALL initialize_blob
        CALL update_ghosts_moments
        CALL solve_EM_fields
        CALL update_ghosts_EM
      ! init moments 00 with a power law similarly to GENE
      CASE('ppj')
        CALL speak('ppj init ~ GENE')
        call init_ppj
        CALL update_ghosts_moments
        CALL solve_EM_fields
        CALL update_ghosts_EM
      CASE('ricci')
        CALL speak('Init Ricci')
        CALL init_ricci ! init only gyrocenter density
        CALL update_ghosts_moments
        CALL solve_EM_fields
        CALL update_ghosts_EM
    END SELECT
    ENDIF
    ! closure of j>J, p>P and j<0, p<0 moments
    CALL speak('Apply closure')
    CALL apply_closure_model
    ! ghosts for p parallelization
    CALL speak('Ghosts communication')
    CALL update_ghosts_moments
    CALL update_ghosts_EM
    !! End of phi and moments initialization
    ! Init collision operator
    CALL init_collision
    !! Preparing auxiliary arrays at initial state
    ! particle density, fluid velocity and temperature (used in diagnose)
    CALL speak('Computing fluid moments and transport')
    CALL compute_fluid_moments
    CALL compute_radial_transport
    CALL compute_radial_heatflux
    ! init auxval for nonlinear
    CALL nonlinear_init
    ! compute nonlinear for t=0 diagnostic
    CALL compute_Sapj ! compute S_0 = S(phi_0,N_0)
  END SUBROUTINE initialize
  !******************************************************************************!

  !******************************************************************************!
  !!!!!!! Initialize all the moments randomly
  !******************************************************************************!
  SUBROUTINE init_moments
    USE grid,       ONLY: local_na, local_np, local_nj, total_nkx, local_nky, local_nz,&
                          ngp, ngj, ngz, iky0, contains_ky0, AA_x, AA_y
    USE fields,     ONLY: moments
    USE prec_const, ONLY: xp
    USE model,      ONLY: LINEARITY
    USE parallel,   ONLY: my_id
    IMPLICIT NONE

    REAL(xp) :: noise
    INTEGER, DIMENSION(12) :: iseedarr
    INTEGER  :: ia,ip,ij,ikx,iky,iz, ipi,iji,izi

    ! Seed random number generator
    iseedarr(:)=iseed
    CALL RANDOM_SEED(PUT=iseedarr+my_id)

      !**** Broad noise initialization *******************************************
    DO ia=1,local_na
      DO ip=1,local_np
        ipi = ip+ngp/2
        DO ij=1,local_nj
          iji = ij+ngj/2
          DO ikx=1,total_nkx
            DO iky=1,local_nky
              DO iz=1,local_nz
                izi = iz+ngz/2
                CALL RANDOM_NUMBER(noise)
                moments(ia,ipi,iji,iky,ikx,izi,:) = (init_background + init_noiselvl*(noise-0.5_xp))
              END DO
            END DO
          END DO
          IF ( contains_ky0 ) THEN
            DO ikx=2,total_nkx/2 !symmetry at ky = 0 for all z
              moments(ia,ipi,iji,iky0,ikx,:,:) = moments(ia,ipi,iji,iky0,total_nkx+2-ikx,:,:)
            END DO
          ENDIF
        END DO
      END DO
      ! Putting to zero modes that are not in the 2/3 Orszag rule
      IF (LINEARITY .NE. 'linear') THEN
        DO ikx=1,total_nkx
          DO iky=1,local_nky
            DO iz=1,local_nz
              izi = iz+ngz/2
              DO ip=1,local_np
                ipi = ip+ngp/2
                DO ij=1,local_nj
                  iji = ij+ngj/2
                  moments(ia,ipi,iji,iky,ikx,izi, :) = moments(ia, ipi,iji,iky,ikx,izi, :)*AA_x(ikx)*AA_y(iky)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  END SUBROUTINE init_moments
  !******************************************************************************!

  !******************************************************************************!
  !!!!!!! Initialize the gyrocenter density randomly
  !******************************************************************************!
  SUBROUTINE init_gyrodens
    USE grid,       ONLY: local_na, local_np, local_nj, total_nkx, local_nky, local_nz,&
                          ngp, ngj, ngz, iky0, parray, jarray, contains_ky0, AA_x, AA_y
    USE fields,     ONLY: moments
    USE prec_const, ONLY: xp
    USE model,      ONLY: LINEARITY
    USE parallel,   ONLY: my_id
    IMPLICIT NONE

    REAL(xp) :: noise
    INTEGER  :: ia,ip,ij,ikx,iky,iz
    INTEGER, DIMENSION(12) :: iseedarr

    ! Seed random number generator
    iseedarr(:)=iseed
    CALL RANDOM_SEED(PUT=iseedarr+my_id)
    moments = 0._xp
      !**** Broad noise initialization *******************************************
    DO ia=1,local_na
      DO ip=1+ngp/2,local_np+ngp/2
        DO ij=1+ngj/2,local_nj+ngj/2
          DO ikx=1,total_nkx
            DO iky=1,local_nky
              DO iz=1+ngz/2,local_nz+ngz/2
                CALL RANDOM_NUMBER(noise)
                IF ( (parray(ip) .EQ. 0) .AND. (jarray(ij) .EQ. 0) ) THEN
                  moments(ia,ip,ij,iky,ikx,iz,:) = (init_background + init_noiselvl*(noise-0.5_xp))
                ELSE
                  moments(ia,ip,ij,iky,ikx,iz,:) = 0._xp
                ENDIF
              END DO
            END DO
          END DO
          IF ( contains_ky0 ) THEN
            DO ikx=2,total_nkx/2 !symmetry at ky = 0 for all z
              moments(ia, ip,ij,iky0,ikx,:,:) = moments(ia, ip,ij,iky0,total_nkx+2-ikx,:,:)
            END DO
          ENDIF
        END DO
      END DO
      ! Putting to zero modes that are not in the 2/3 Orszag rule
      IF (LINEARITY .NE. 'linear') THEN
        DO ikx=1,total_nkx
          DO iky=1,local_nky
            DO iz=1,local_nz+ngz
              DO ip=1,local_np+ngp
                DO ij=1,local_nj+ngj
                  moments(ia, ip,ij,iky,ikx,iz, :) = moments(ia, ip,ij,iky,ikx,iz, :)*AA_x(ikx)*AA_y(iky)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  END SUBROUTINE init_gyrodens
  !******************************************************************************!

  !******************************************************************************!
  !!!!!!! Initialize a single mode in the gyrocenter density
  !******************************************************************************!
  SUBROUTINE init_single_mode
    USE fields,     ONLY: moments
    USE prec_const, ONLY: xp
    USE parallel,   ONLY: my_id
    IMPLICIT NONE
    moments   = 0._xp
    IF (my_id .EQ. 0) THEN
      moments(:,:,:,2,2,:,:) = init_amp
      moments(:,:,:,3,3,:,:) = init_amp
      moments(:,:,:,4,4,:,:) = init_amp
      moments(:,:,:,5,5,:,:) = init_amp
    ENDIF
  END SUBROUTINE init_single_mode
  !******************************************************************************!

  !******************************************************************************!
  !!!!!!! Initialize a noisy ES potential and cancel the moments
  !******************************************************************************!
  SUBROUTINE init_phi
    USE grid,       ONLY: total_nkx, local_nky, local_nz,&
                          ngz, iky0, ikx0, contains_ky0, AA_x, AA_y
    USE fields,     ONLY: phi, moments
    USE prec_const, ONLY: xp
    USE model,      ONLY: LINEARITY
    USE parallel,   ONLY: my_id
    IMPLICIT NONE
    REAL(xp) :: noise
    INTEGER, DIMENSION(12) :: iseedarr
    INTEGER :: iky,ikx,iz
    ! Seed random number generator
    iseedarr(:)=iseed
    CALL RANDOM_SEED(PUT=iseedarr+my_id)
    !**** noise initialization *******************************************
    DO ikx=1,total_nkx
      DO iky=1,local_nky
        DO iz=1,local_nz+ngz
          CALL RANDOM_NUMBER(noise)
          phi(iky,ikx,iz) = (init_background + init_noiselvl*(noise-0.5_xp))!*AA_x(ikx)*AA_y(iky)
        ENDDO
      END DO
    END DO
    !symmetry at ky = 0 to keep real inverse transform
    IF ( contains_ky0 ) THEN
      DO iz=1,local_nz+ngz
        DO ikx=2,total_nkx/2
          phi(iky0,ikx,iz) = phi(iky0,total_nkx+2-ikx,iz)
        ENDDO
      phi(iky0,ikx0,iz) = REAL(phi(iky0,ikx0,iz),xp) !origin must be real
    END DO
    ENDIF
    ! Putting to zero modes that are not in the 2/3 Orszag rule
    IF (LINEARITY .NE. 'linear') THEN
      DO ikx=1,total_nkx
        DO iky=1,local_nky
          DO iz=1,local_nz+ngz
            phi(iky,ikx,iz) = phi(iky,ikx,iz)*AA_x(ikx)*AA_y(iky)
          ENDDO
        ENDDO
      ENDDO
    ENDIF
    !**** ensure no previous moments initialization
    moments = 0._xp
  END SUBROUTINE init_phi
  !******************************************************************************!

  !******************************************************************************!
  !!!!!!! Initialize a ppj ES potential and cancel the moments
  !******************************************************************************!
  SUBROUTINE init_phi_ppj
    USE grid,       ONLY: total_nkx, local_nky, local_nz, AA_x, AA_y,&
                          ngz, iky0, ikx0, contains_ky0, ieven, kxarray, kyarray, zarray, deltakx
    USE fields,     ONLY: phi, moments
    USE prec_const, ONLY: xp
    USE model,      ONLY: LINEARITY
    USE geometry,   ONLY: Jacobian, iInt_Jacobian
    IMPLICIT NONE
    REAL(xp) :: kx, ky, z, amp
    INTEGER  :: ikx, iky, iz
    amp = 1.0_xp
      !**** ppj initialization *******************************************
        DO iky=1,local_nky
          ky = kyarray(iky)
          DO ikx=1,total_nkx
            kx = kxarray(iky,ikx)
            DO iz=1,local_nz+ngz
              z = zarray(iz,ieven)
              IF (ky .NE. 0) THEN
                phi(iky,ikx,iz) = 0._xp
              ELSE
                phi(iky,ikx,iz) = 0.5_xp*amp*(deltakx/(ABS(kx)+deltakx))
              ENDIF
              ! z-dep and noise
              phi(iky,ikx,iz) = phi(iky,ikx,iz) * &
              (Jacobian(iz,ieven)*iInt_Jacobian)**2
            END DO
          END DO
        END DO
      !symmetry at ky = 0 to keep real inverse transform
      IF ( contains_ky0 ) THEN
        DO iz=1,local_nz+ngz
          DO ikx=2,total_nkx/2
            phi(iky0,ikx,iz) = phi(iky0,total_nkx+2-ikx,iz)
          ENDDO
          phi(iky0,ikx0,iz) = REAL(phi(iky0,ikx0,iz),xp) !origin must be real
        END DO
      ENDIF
      ! Putting to zero modes that are not in the 2/3 Orszag rule
      IF (LINEARITY .NE. 'linear') THEN
        DO ikx=1,total_nkx
          DO iky=1,local_nky
            DO iz=1,local_nz+ngz
              phi(iky,ikx,iz) = phi(iky,ikx,iz)*AA_x(ikx)*AA_y(iky)
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      !**** ensure no previous moments initialization
      moments = 0._xp
  END SUBROUTINE init_phi_ppj
  !******************************************************************************!

  !******************************************************************************!
  !******************************************************************************!
  !!!!!!! Initialize an ionic Gaussian blob on top of the preexisting modes
  !******************************************************************************!
  SUBROUTINE initialize_blob
    USE grid,       ONLY: local_na, local_np, local_nj, total_nkx, local_nky, local_nz, total_nz,&
                          AA_x, AA_y, parray, jarray,&
                          ngp,ngj,ngz, iky0, ieven, kxarray, kyarray, zarray
    USE fields,     ONLY: moments
    USE prec_const, ONLY: xp
    USE model,      ONLY: LINEARITY
    USE geometry,   ONLY: Jacobian, iInt_Jacobian, shear
    IMPLICIT NONE
    REAL(xp) ::kx, ky, z, sigma_x, sigma_y, gain
    INTEGER :: ia,iky,ikx,iz,ip,ij, p, j
    sigma_y = 0.5
    sigma_x = sigma_y
    gain  = 1.0
    ! One can increase the gain if we run 3D sim
    IF(total_nz .GT. 1) THEN
      sigma_y = 1.0
      sigma_x = sigma_y
      gain = 10.0
    ENDIF

    DO ia=1,local_na
      DO iky=1,local_nky
        ky = kyarray(iky)
        DO iz=1+ngz/2,local_nz+ngz/2
          z  = zarray(iz,ieven)
          DO ikx=1,total_nkx
            kx = kxarray(iky,ikx) + z*shear*ky
            DO ip=1+ngp/2,local_np+ngp/2
              p = parray(ip)
              DO ij=1+ngj/2,local_nj+ngj/2
                j = jarray(ij)
                IF( (iky .NE. iky0) .AND. (p .EQ. 0) .AND. (j .EQ. 0)) THEN
                  moments(ia,ip,ij,iky,ikx,iz, :) = moments(ia,ip,ij,iky,ikx,iz, :) &
                  + gain * exp(-((kx/sigma_x)**2+(ky/sigma_y)**2)) &
                    * AA_x(ikx)*AA_y(iky)* &
                    (Jacobian(iz,ieven)*iInt_Jacobian)**2!&
                    ! * exp(sigmai2_taui_o2*(kx**2+ky**2))
                ENDIF
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      ! Putting to zero modes that are not in the 2/3 Orszag rule
      IF (LINEARITY .NE. 'linear') THEN
        DO ikx=1,total_nkx
          DO iky=1,local_nky
            DO iz=1,local_nz+ngz
              DO ip=1,local_np+ngp
                DO ij=1,local_nj+ngj
                  moments(ia, ip,ij,iky,ikx,iz, :) = moments(ia, ip,ij,iky,ikx,iz, :)*AA_x(ikx)*AA_y(iky)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  END SUBROUTINE initialize_blob
  !******************************************************************************!
  !******************************************************************************!
  !!!!!!! Initialize the gyrocenter in a ppj gene manner (power law)
  !******************************************************************************!
  SUBROUTINE init_ppj
    USE grid,       ONLY: local_na, local_np, local_nj, total_nkx, local_nky, local_nz,&
                          AA_x, AA_y, deltakx, deltaky,contains_ky0,&
                          ngp,ngj,ngz, iky0, ieven, kxarray, kyarray, zarray
    USE fields,     ONLY: moments
    USE prec_const, ONLY: xp, pi
    USE model,      ONLY: LINEARITY
    USE geometry,   ONLY: Jacobian, iInt_Jacobian
    IMPLICIT NONE
    REAL(xp) :: kx, ky, sigma_z, z
    INTEGER :: ia,iky,ikx,iz,ip,ij
    sigma_z = pi/4._xp
      !**** Broad noise initialization *******************************************
      DO ia=1,local_na
        DO ip=1,local_np+ngp
          DO ij=1,local_nj+ngj
            IF ( (ip .EQ. 1) .AND. (ij .EQ. 1) ) THEN
              DO iky=1,local_nky
                ky = kyarray(iky)
                DO ikx=1,total_nkx
                  kx = kxarray(iky,ikx)
                  DO iz=1,local_nz+ngz
                    z = zarray(iz,ieven)
                    IF (kx .EQ. 0) THEN
                      IF(ky .EQ. 0) THEN
                        moments(ia,ip,ij,iky,ikx,iz,:) = 0._xp
                      ELSE
                        moments(ia,ip,ij,iky,ikx,iz,:) = 0.5_xp * deltaky/(ABS(ky)+deltaky)
                      ENDIF
                    ELSE
                      IF(ky .GT. 0) THEN
                        moments(ia,ip,ij,iky,ikx,iz,:) = (deltakx/(ABS(kx)+deltakx))*(deltaky/(ABS(ky)+deltaky))
                      ELSE
                        moments(ia,ip,ij,iky,ikx,iz,:) = 0.5_xp*(deltakx/(ABS(kx)+deltakx))
                      ENDIF
                    ENDIF
                    ! z-dep and noise
                    moments(ia,ip,ij,iky,ikx,iz,:) = moments(ia,ip,ij,iky,ikx,iz,:) * &
                    (Jacobian(iz,ieven)*iInt_Jacobian)**2
                    ! divide by kernel_0 to adjust to particle density (n = kernel_0 N00)
                    ! moments(ia,ip,ij,iky,ikx,iz,:) = moments(ia,ip,ij,iky,ikx,iz,:)/kernel(ia,ij,iky,ikx,iz,0)
                  END DO
                END DO
              END DO
              IF ( contains_ky0 ) THEN
                DO ikx=2,total_nkx/2 !symmetry at kx = 0 for all z
                  moments(ia,ip,ij,iky0,ikx,:,:) = moments(ia, ip,ij,iky0,total_nkx+2-ikx,:, :)
                END DO
              ENDIF
          ELSE
            moments(ia,ip,ij,:,:,:,:) = 0._xp
          ENDIF
        END DO
      END DO
      ! Putting to zero modes that are not in the 2/3 Orszag rule
      IF (LINEARITY .NE. 'linear') THEN
        DO ikx=1,total_nkx
        DO iky=1,local_nky
        DO iz=1,local_nz+ngz
          DO ip=1,local_np+ngp
          DO ij=1,local_nj+ngj
            moments(ia, ip,ij,iky,ikx,iz, :) = moments(ia, ip,ij,iky,ikx,iz, :)*AA_x(ikx)*AA_y(iky)
          ENDDO
          ENDDO
        ENDDO
        ENDDO
        ENDDO
      ENDIF
    ENDDO
  END SUBROUTINE init_ppj
  !******************************************************************************!

  !******************************************************************************!
  !!!!!!! Initialize Ricci density
  ! We take the picture of Paolo Ricci and use it as input for the initialization
  ! of the gyrodensity. The picture is expected to be 186x96 fourier modes stored
  ! in the gyacomo/Gallery directory, in two 2D matrices (real.txt and imag.txt)  
  ! with only spaces as separators. 
  ! One can scale the amplitude of the modes with the init_background parameter 
  ! (1e-5 approx is needed to calm down the system if it is in nonlinear evolution)
  !******************************************************************************!
  SUBROUTINE init_ricci
    USE basic,      ONLY: maindir
    USE grid,       ONLY: local_na, local_np, local_nj, total_nkx, local_nky, local_nz,&
                          local_nkx_offset, local_nky_offset, kxarray, kyarray, &
                          ngp, ngj, ngz, parray, jarray,&
                          deltakx, deltaky,&
                          AA_x, AA_y
    USE fields,     ONLY: moments
    USE prec_const, ONLY: xp, imagu
    USE model,      ONLY: LINEARITY
    IMPLICIT NONE
    REAL(xp), DIMENSION(186,94) :: ricci_mat_real, ricci_mat_imag
    REAL(xp) :: scaling
    INTEGER  :: ia,ip,ij,ikx,iky,iz, LPFx, LPFy
    CHARACTER(256) ::  filename

    ! open data file
    ricci_mat_real = 0; ricci_mat_imag = 0
    filename = TRIM(maindir) // '/Gallery/fourier_ricci_real.txt'
    OPEN(unit = 1 , file = filename)
    READ(1,*) ricci_mat_real
    CLOSE(1)
    filename = TRIM(maindir) // '/Gallery/fourier_ricci_imag.txt'
    OPEN(unit = 2 , file = filename)
    READ(2,*) ricci_mat_imag
    CLOSE(2)
    scaling = init_amp*1e-5
    moments = 0._xp
      !**** initialization *******************************************
    DO ia=1,local_na
      DO ikx=1,total_nkx
        DO iky=1,local_nky
          DO ip=1+ngp/2,local_np+ngp/2
            DO ij=1+ngj/2,local_nj+ngj/2
              DO iz=1+ngz/2,local_nz+ngz/2
                IF((ikx+local_nkx_offset .LE. 186) .AND. (iky+local_nky_offset .LE. 94)) THEN
                  !IF (.TRUE.) THEN
                  IF ( (parray(ip) .EQ. 0) .AND. (jarray(ij) .EQ. 0) ) THEN
                    moments(ia,ip,ij,iky,ikx,iz,:) = scaling*(ricci_mat_real(ikx+local_nkx_offset,iky+local_nky_offset)&
                    - imagu*ricci_mat_imag(ikx+local_nkx_offset,iky+local_nky_offset))
                  ELSE
                    moments(ia,ip,ij,iky,ikx,iz,:) = 0._xp
                  ENDIF
                ELSE
                  moments(ia,ip,ij,iky,ikx,iz,:) = 0._xp
                ENDIF
              END DO
            END DO
          END DO
        END DO
      END DO
      ! Putting to zero modes that are not in the 2/3 Orszag rule
      IF (LINEARITY .NE. 'linear') THEN
        DO ikx=1,total_nkx
          DO iky=1,local_nky
            DO iz=1,local_nz+ngz
              DO ip=1,local_np+ngp
                DO ij=1,local_nj+ngj
                  moments(ia, ip,ij,iky,ikx,iz, :) = moments(ia, ip,ij,iky,ikx,iz, :)*AA_x(ikx)*AA_y(iky)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDIF
    ENDDO
    !! Play with some filtering
    LPFx = 0
    LPFy = 0
      DO ikx=1,total_nkx
        DO iky=1,local_nky
          IF ( (abs(kxarray(ikx,iky)) .LT. LPFx*deltakx ) .AND. (abs(kyarray(iky)) .LT. LPFy*deltaky ) )&
            moments(:,:,:,iky,ikx,:,:) = 0._xp
        ENDDO
      ENDDO
  END SUBROUTINE init_ricci
  !******************************************************************************!


  SUBROUTINE initial_outputinputs(fid)
    ! Write the input parameters to the results_xx.h5 file
    USE futils, ONLY: attach, creatd
    IMPLICIT NONE
    INTEGER, INTENT(in) :: fid
    CHARACTER(len=256)  :: str
    WRITE(str,'(a)') '/data/input/intial'
    CALL creatd(fid, 0,(/0/),TRIM(str),'Initial Input')
    CALL attach(fid, TRIM(str), "INIT_OPT", INIT_OPT)
    CALL attach(fid, TRIM(str), "init_background", init_background)
    CALL attach(fid, TRIM(str), "init_noiselvl", init_noiselvl)
    CALL attach(fid, TRIM(str), "iseed", iseed)
  END SUBROUTINE initial_outputinputs

END MODULE initial
