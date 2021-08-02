!******************************************************************************!
!!!!!! initialize the moments and load/build coeff tables
!******************************************************************************!
SUBROUTINE inital

  USE basic
  USE model, ONLY : CO, NON_LIN
  USE initial_par
  USE prec_const
  USE time_integration
  USE array, ONLY : Sepj,Sipj
  USE collision
  USE closure
  USE ghosts
  USE restarts
  IMPLICIT NONE

  CALL set_updatetlevel(1)

  !!!!!! Set the moments arrays Nepj, Nipj and phi!!!!!!
  ! through loading a previou state
  IF ( RESTART ) THEN
    IF (my_id .EQ. 0) WRITE(*,*) 'Load moments'
    CALL load_moments ! get N_0
    CALL poisson ! compute phi_0=phi(N_0)
  ! through initialization
  ELSE
    ! set phi with noise and set moments to 0
    IF (INIT_NOISY_PHI) THEN
      CALL init_phi

    ! set moments_00 (GC density) with noise and compute phi afterwards
    ELSE
      IF (my_id .EQ. 0) WRITE(*,*) 'Init noisy moments'
      CALL init_moments ! init noisy N_0
      CALL poisson ! get phi_0 = phi(N_0)
    ENDIF
  ENDIF

  ! Option for wiping the turbulence and check growth of secondary inst.
  IF ( WIPE_TURB ) THEN
    IF (my_id .EQ. 0) WRITE(*,*) '-Wiping turbulence'
    CALL wipe_turbulence
  ENDIF
  ! Option for initializing a gaussian blob on the zonal profile
  IF ( INIT_BLOB ) THEN
    IF (my_id .EQ. 0) WRITE(*,*) '--init a blob'
    CALL put_blob
  ENDIF

  IF (my_id .EQ. 0) WRITE(*,*) 'Apply closure'
  CALL apply_closure_model

  IF (my_id .EQ. 0) WRITE(*,*) 'Ghosts communication'
  CALL update_ghosts

  !!!!!! Set Sepj, Sipj and dnjs coeff table !!!!!!
  IF ( NON_LIN ) THEN;
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
  REAL(dp) :: kx, ky, sigma, gain, ky_shift
  INTEGER, DIMENSION(12) :: iseedarr

  ! Seed random number generator
  iseedarr(:)=iseed
  CALL RANDOM_SEED(PUT=iseedarr+my_id)

    !**** Broad noise initialization *******************************************
    DO ip=ips_e,ipe_e
      DO ij=ijs_e,ije_e

        DO ikx=ikxs,ikxe
          DO iky=ikys,ikye
            DO iz=izs,ize
              CALL RANDOM_NUMBER(noise)
              moments_e(ip,ij,ikx,iky,iz,:) = (init_background + init_noiselvl*(noise-0.5_dp))
            END DO
          END DO
        END DO

        IF ( contains_kx0 ) THEN
          DO iky=2,Nky/2 !symmetry at kx = 0 for all z
            moments_e(ip,ij,ikx_0,iky,:,:) = moments_e( ip,ij,ikx_0,Nky+2-iky,:, :)
          END DO
        ENDIF

      END DO
    END DO

    DO ip=ips_i,ipe_i
      DO ij=ijs_i,ije_i

        DO ikx=ikxs,ikxe
          DO iky=ikys,ikye
            DO iz=izs,ize
              CALL RANDOM_NUMBER(noise)
              moments_i(ip,ij,ikx,iky,iz,:) = (init_background + init_noiselvl*(noise-0.5_dp))
            END DO
          END DO
        END DO

        IF ( contains_kx0 ) THEN
          DO iky=2,Nky/2 !symmetry at kx = 0 for all z
            moments_i( ip,ij,ikx_0,iky,:,:) = moments_i( ip,ij,ikx_0,Nky+2-iky,:,:)
          END DO
        ENDIF

      END DO
    END DO

    ! Putting to zero modes that are not in the 2/3 Orszag rule
    IF (NON_LIN) THEN
      DO ikx=ikxs,ikxe
      DO iky=ikys,ikye
      DO iz=izs,ize
        DO ip=ips_e,ipe_e
        DO ij=ijs_e,ije_e
          moments_e( ip,ij,ikx,iky,iz, :) = moments_e( ip,ij,ikx,iky,iz, :)*AA_x(ikx)*AA_y(iky)
        ENDDO
        ENDDO
        DO ip=ips_i,ipe_i
        DO ij=ijs_i,ije_i
          moments_i( ip,ij,ikx,iky,iz, :) = moments_i( ip,ij,ikx,iky,iz, :)*AA_x(ikx)*AA_y(iky)
        ENDDO
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
  REAL(dp) :: kx, ky, sigma, gain, ky_shift
  INTEGER, DIMENSION(12) :: iseedarr

  IF (INIT_NOISY_PHI) THEN
    IF (my_id .EQ. 0) WRITE(*,*) 'Init noisy phi'
    ! Seed random number generator
    iseedarr(:)=iseed
    CALL RANDOM_SEED(PUT=iseedarr+my_id)

      !**** noise initialization *******************************************

      DO ikx=ikxs,ikxe
        DO iky=ikys,ikye
          DO iz=izs,ize
            CALL RANDOM_NUMBER(noise)
            phi(ikx,iky,iz) = (init_background + init_noiselvl*(noise-0.5_dp))*AA_x(ikx)*AA_y(iky)
          ENDDO
        END DO
      END DO

      !symmetry at kx = 0 to keep real inverse transform
      IF ( contains_kx0 ) THEN
        DO iky=2,Nky/2
          phi(ikx_0,iky,:) = phi(ikx_0,Nky+2-iky,:)
        END DO
        phi(ikx_0,Ny/2,:) = REAL(phi(ikx_0,Ny/2,:)) !origin must be real
      ENDIF

      !**** ensure no previous moments initialization
      moments_e = 0._dp; moments_i = 0._dp

      !**** Zonal Flow initialization *******************************************
      ! put a mode at ikx = mode number + 1, symmetry is already included since kx>=0
      IF(INIT_ZF .GT. 0) THEN
      IF (my_id .EQ. 0) WRITE(*,*) 'Init ZF phi'
        IF( (INIT_ZF+1 .GT. ikxs) .AND. (INIT_ZF+1 .LT. ikxe) ) THEN
          DO iz = izs,ize
            phi(INIT_ZF+1,iky_0,iz) = ZF_AMP*(2._dp*PI)**2/deltakx/deltaky/2._dp * COS((iz-1)/Nz*2._dp*PI)
            moments_i(1,1,INIT_ZF+1,iky_0,iz,:) = kxarray(INIT_ZF+1)**2*phi(INIT_ZF+1,iky_0,iz)* COS((iz-1)/Nz*2._dp*PI)
            moments_e(1,1,INIT_ZF+1,iky_0,iz,:) = 0._dp
          ENDDO
        ENDIF
      ENDIF
    ELSE ! we compute phi from noisy moments and poisson
      CALL poisson
    ENDIF

END SUBROUTINE init_phi
!******************************************************************************!

!******************************************************************************!
!!!!!!! Remove all ky!=0 modes to conserve only zonal modes in a restart
!******************************************************************************!
SUBROUTINE wipe_turbulence
  USE fields
  USE grid
  IMPLICIT NONE
  DO ikx=ikxs,ikxe
  DO iky=ikys,ikye
  DO iz=izs,ize
    DO ip=ips_e,ipe_e
    DO ij=ijs_e,ije_e
      IF( iky .NE. iky_0) THEN
        moments_e( ip,ij,ikx,iky,iz, :) = 0e-3_dp*moments_e( ip,ij,ikx,iky,iz, :)
      ELSE
        moments_e( ip,ij,ikx,iky,iz, :) = 1e+0_dp*moments_e( ip,ij,ikx,iky,iz, :)
      ENDIF
    ENDDO
    ENDDO
    DO ip=ips_i,ipe_i
    DO ij=ijs_i,ije_i
      IF( iky .NE. iky_0) THEN
        moments_i( ip,ij,ikx,iky,iz, :) = 0e-3_dp*moments_i( ip,ij,ikx,iky,iz, :)
      ELSE
        moments_i( ip,ij,ikx,iky,iz, :) = 1e+0_dp*moments_i( ip,ij,ikx,iky,iz, :)
      ENDIF
    ENDDO
    ENDDO
  ENDDO
  ENDDO
  ENDDO
END SUBROUTINE
!******************************************************************************!
!******************************************************************************!
!!!!!!! Initialize an ionic Gaussian blob on top of the preexisting modes
!******************************************************************************!
SUBROUTINE put_blob
  USE fields
  USE grid
  USE model, ONLY: sigmai2_taui_o2
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
        moments_i( ip,ij,ikx,iky,iz, :) = moments_i( ip,ij,ikx,iky,iz, :) &
        + gain*sigma/SQRT2 * exp(-(kx**2+ky**2)*sigma**2/4._dp) &
          * AA_x(ikx)*AA_y(iky)!&
          ! * exp(sigmai2_taui_o2*(kx**2+ky**2))
      ENDIF
    ENDDO
    ENDDO
  ENDDO
  ENDDO
  ENDDO
END SUBROUTINE put_blob
!******************************************************************************!
