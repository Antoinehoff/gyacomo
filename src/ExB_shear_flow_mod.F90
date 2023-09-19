MODULE ExB_shear_flow
    ! This module contains the necessary tools to implement ExB shearing flow effects.
    ! The algorithm is taken from the presentation of Hammett et al. 2006 (APS) and
    ! it the one used in GS2.
    USE prec_const, ONLY: xp, imagu, pi

    IMPLICIT NONE
    ! Variables
    REAL(xp),   PUBLIC, PROTECTED :: gamma_E = 0._xp     ! ExB background shearing rate \gamma_E
    REAL(xp),   PUBLIC, PROTECTED :: t0, inv_t0 = 0._xp  ! charact. shear time
    REAL(xp),   DIMENSION(:),   ALLOCATABLE, PUBLIC, PROTECTED :: sky_ExB      ! shift of the kx modes, kx* = kx + s(ky)
    REAL(xp),   DIMENSION(:),   ALLOCATABLE, PUBLIC, PROTECTED :: sky_ExB_full ! full ky version
    REAL(xp),   DIMENSION(:),   ALLOCATABLE, PUBLIC, PROTECTED :: dkx_ExB      ! correction to obtain the exact kx mode
    INTEGER,    DIMENSION(:),   ALLOCATABLE, PUBLIC, PROTECTED :: jump_ExB     ! jump to do to shift the kx grids
    LOGICAL,    DIMENSION(:),   ALLOCATABLE, PUBLIC, PROTECTED :: shiftnow_ExB ! Indicates if there is a line to shift
    COMPLEX(xp),DIMENSION(:,:), ALLOCATABLE, PUBLIC, PROTECTED :: ExB_NL_factor! factor for nonlinear term
    COMPLEX(xp),DIMENSION(:,:), ALLOCATABLE, PUBLIC, PROTECTED :: inv_ExB_NL_factor
    LOGICAL,  PUBLIC, PROTECTED ::     ExB =  .false.   ! presence of ExB background shearing rate
    ! Routines
    PUBLIC :: Setup_ExB_shear_flow, Array_shift_ExB_shear_flow, Update_ExB_shear_flow

CONTAINS

    ! Setup the variables for the ExB shear
    SUBROUTINE Setup_ExB_shear_flow(ExBrate)
        USE grid,     ONLY: Nx, local_nky, total_nky, local_nx, Ny, deltakx, deltaky,&
                            kyarray, kyarray_full
        USE geometry, ONLY: Cyq0_x0
        USE basic,    ONLY: dt
        IMPLICIT NONE
        INTEGER :: iky
        REAL(xp), INTENT(IN) :: ExBrate

        ! Setup the ExB shearing rate and aux var
        gamma_E = -ExBrate*Cyq0_x0
        IF(abs(gamma_E) .GT. EPSILON(gamma_E)) THEN
            ExB    = .TRUE.
            t0     = deltakx/deltaky/gamma_E
            inv_t0 = 1._xp/t0
        ELSE ! avoid 1/0 division (t0 is killed anyway in this case)
            ExB    = .FALSE.
            t0     = 0._xp
            inv_t0 = 0._xp
        ENDIF

        ! Setup the ExB shift array
        ALLOCATE(sky_ExB(local_nky))
        ! consider no initial shift (maybe changed if restart)
        sky_ExB = 0._xp
        ! Midpoint init
        DO iky = 1,local_nky
            sky_ExB(iky) = sky_ExB(iky) !+ 0.5_xp*kyarray(iky)*gamma_E*dt
        ENDDO

        ALLOCATE(sky_ExB_full(total_nky+1))
        ! consider no initial shift (maybe changed if restart)
        sky_ExB_full = 0._xp
        ! Midpoint init
        DO iky = 1,total_nky+1
            sky_ExB_full(iky) = sky_ExB_full(iky) !+ 0.5_xp*REAL(iky-1,xp)*deltaky*gamma_E*dt
        ENDDO

        ! Setup the kx correction array
        ALLOCATE(dkx_ExB(local_nky))
        dkx_ExB = 0._xp

        ! Setup the jump array
        ALLOCATE(jump_ExB(local_nky))
        jump_ExB = 0

        ! Setup the shifting flag array
        ALLOCATE(shiftnow_ExB(local_nky))
        shiftnow_ExB = .FALSE.

        ! Setup nonlinear factor
        ALLOCATE(    ExB_NL_factor(Nx,local_nky))
        ALLOCATE(inv_ExB_NL_factor(Ny/2+2,local_nx))
            ExB_NL_factor = 1._xp
        inv_ExB_NL_factor = 1._xp

    END SUBROUTINE Setup_ExB_shear_flow

    ! update the ExB shear value for the next time step
    SUBROUTINE Update_ExB_shear_flow(step_number)
        USE basic,      ONLY: dt,time
        USE grid,       ONLY: local_nky, total_nky, kyarray, inv_dkx, kyarray_full, update_grids, deltaky
        USE geometry,   ONLY: gxx,gxy,gyy,inv_hatB2, evaluate_magn_curv
        USE numerics,   ONLY: evaluate_EM_op, evaluate_kernels
        USE model,      ONLY: LINEARITY
        USE time_integration, ONLY: c_E
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: step_number
        ! local var
        INTEGER :: iky
        REAL(xp):: tnext

        ! tnext = time + c_E(step_number)*dt
        tnext = time + 0._xp*dt
        ! do nothing if no ExB
        IF(ExB) THEN
            ! Update new shear value
            DO iky = 1,local_nky
                !! This must be done incrementely to be able to pull it back
                !  when a grid shift occurs                
                sky_ExB(iky)      = sky_ExB(iky) - kyarray(iky)*gamma_E*dt
                jump_ExB(iky)     = NINT(sky_ExB(iky)*inv_dkx)
                ! If the jump is 1 or more for a given ky, we flag the index
                ! in shiftnow_ExB and will use it in Shift_fields to avoid
                ! zero-shiftings that may be majoritary.
                shiftnow_ExB(iky) = (abs(jump_ExB(iky)) .GT. 0)
            ENDDO
            ! Update the full skyExB array too
            DO iky = 1,total_nky+1
                sky_ExB_full(iky) = sky_ExB_full(iky) - REAL(iky-1,xp)*deltaky*gamma_E*dt
            ENDDO
            ! Shift the arrays if the shear value sky is too high
            CALL Array_shift_ExB_shear_flow

            ! We update the operators and grids
            !   update the grids  
            CALL update_grids(sky_ExB,gxx,gxy,gyy,inv_hatB2)
            !   update the EM op., the kernels and the curvature op.
            CALL evaluate_kernels
            CALL evaluate_EM_op
            CALL evaluate_magn_curv
            !   update the ExB nonlinear factor...
            IF(LINEARITY .EQ. 'nonlinear') &
             CALL update_nonlinear_ExB_factors
        ENDIF
    END SUBROUTINE Update_ExB_shear_flow

    ! According to the current ExB shear value we update
    ! the fields by imposing a shift on kx
    SUBROUTINE Array_shift_ExB_shear_flow
        USE grid,       ONLY: local_nky, total_nky, update_grids, &
            total_nkx, deltakx, kx_min, kx_max, kxarray0, inv_dkx
        USE prec_const, ONLY: PI
        USE fields,     ONLY: moments, phi, psi
        USE numerics,   ONLY: evaluate_EM_op, evaluate_kernels
        IMPLICIT NONE
        ! local var
        INTEGER :: iky, ikx, ikx_s, i_, loopstart, loopend, increment, jump_
        IF(ExB) THEN
            ! shift all local fields and correct the local shift value
            DO iky = 1,local_Nky
                IF(shiftnow_ExB(iky)) THEN
                    ! We shift the array from left to right or right to left according to the jump
                    ! This avoids to make copy
                    IF(jump_ExB(iky) .GT. 0) THEN
                        loopstart = 1
                        loopend   = total_nkx
                        increment = 1
                    ELSE
                        loopstart = total_nkx
                        loopend   = 1
                        increment = -1
                    ENDIF
                    !loop to go through the array in a monotonic kx order
                    ! Recall: the kx array is organized as
                    !   6   7   8   1   2   3   4   5 (indices ikx)
                    !  -3  -2  -1   0   1   2   3   4 (values in dkx)
                    ! so to go along the array in a monotonic way one must travel as
                    ! 67812345 or 54321678
                    DO i_ = loopstart, loopend, increment
                        IF (i_ .LT. total_nkx/2) THEN ! go to the negative kx region
                            ikx = i_ + total_nkx/2 + 1
                        ELSE ! positive
                            ikx = i_ - total_nkx/2 + 1
                        ENDIF
                        ikx_s = ikx + jump_ExB(iky)
                        ! adjust the shift accordingly
                        IF (ikx_s .LE. 0) &
                            ikx_s = ikx_s + total_nkx
                        IF (ikx_s .GT. total_nkx) &
                            ikx_s = ikx_s - total_nkx
                        ! Then we test if the shifted modes are still in contained in our resolution
                        IF ( (kxarray0(ikx)+jump_ExB(iky)*deltakx .LE. kx_max) .AND. &
                             (kxarray0(ikx)+jump_ExB(iky)*deltakx .GE. kx_min)) THEN
                                moments(:,:,:,iky,ikx,:,:) = moments(:,:,:,iky,ikx_s,:,:)
                                phi(iky,ikx,:)             = phi(iky,ikx_s,:)
                                psi(iky,ikx,:)             = psi(iky,ikx_s,:)
                        ELSE ! if it is not, it is lost (~dissipation for high modes)
                            moments(:,:,:,iky,ikx,:,:) = 0._xp
                            phi(iky,ikx,:)             = 0._xp
                            psi(iky,ikx,:)             = 0._xp
                        ENDIF
                    ENDDO
                    ! correct the shift value s(ky) for this row
                    sky_ExB(iky) = sky_ExB(iky) - jump_ExB(iky)*deltakx
                    ! reset the flag
                    shiftnow_ExB(iky) = .FALSE.
                ENDIF
            ENDDO
        ENDIF
        ! Check the global shift values
        DO iky = 1,total_nky+1
            jump_ = NINT(sky_ExB_full(iky)*inv_dkx)
            IF (ABS(jump_) .GT. 0) &
                sky_ExB_full(iky) = sky_ExB_full(iky) - jump_*deltakx
        ENDDO

    END SUBROUTINE Array_shift_ExB_shear_flow

    SUBROUTINE Update_nonlinear_ExB_factors
        USE grid,  ONLY: local_nky, local_nky_offset, xarray, Nx, Ny, local_nx, deltakx,&
                         local_nx_offset, ikyarray, inv_ikyarray, deltaky, update_grids
        USE basic, ONLY: time, dt
        IMPLICIT NONE
        INTEGER :: iky, ix
        REAL(xp):: dt_ExB, J_xp, inv_J, xval, tnext
        tnext = time + dt
        DO iky = 1,local_nky ! WARNING: Local indices ky loop
            ! for readability
            ! J_xp  = ikyarray(iky+local_nky_offset)
            ! inv_J = inv_ikyarray(iky+local_nky_offset)
            J_xp   = REAL(iky-1,xp)
            IF(J_xp .GT. 0._xp) THEN
                inv_J  = 1._xp/J_xp
            ELSE
                inv_J  = 0._xp
            ENDIF
            ! compute dt factor
            dt_ExB = (tnext - t0*inv_J*ANINT(J_xp*tnext*inv_t0,xp))
            DO ix = 1,Nx
                xval = 2._xp*pi/deltakx*REAL(ix-1,xp)/REAL(Nx,xp)!xarray(ix)
                ! assemble the ExB nonlin factor
                ! ExB_NL_factor(ix,iky) = EXP(-imagu*x*gamma_E*J_xp*deltaky*dt_ExB)
                ExB_NL_factor(ix,iky) = EXP(-imagu*sky_ExB(iky)*xval)               
                ! ExB_NL_factor(ix,iky) = EXP(-imagu*sky_ExB_full(iky+local_nky_offset)*xval)               
            ENDDO
        ENDDO
        ! ... and the inverse
        DO iky = 1,Ny/2+2 ! WARNING: Global indices ky loop
            ! for readability
            J_xp   = REAL(iky-1,xp)
            IF(J_xp .GT. 0._xp) THEN
                inv_J  = 1._xp/J_xp
            ELSE
                inv_J  = 0._xp
            ENDIF
            ! compute dt factor
            dt_ExB = (tnext - t0*inv_J*ANINT(J_xp*tnext*inv_t0,xp))
            DO ix = 1,local_nx
                xval = 2._xp*pi/deltakx*REAL(ix-1,xp)/REAL(Nx,xp)!xarray(ix+local_nx_offset)
                ! assemble the inverse ExB nonlin factor
                ! inv_ExB_NL_factor(iky,ix) = EXP(imagu*x*gamma_E*J_xp*deltaky*dt_ExB)
                inv_ExB_NL_factor(iky,ix) = EXP(imagu*sky_ExB_full(iky)*xval)
            ENDDO
        ENDDO
    END SUBROUTINE Update_nonlinear_ExB_factors

END MODULE ExB_shear_flow