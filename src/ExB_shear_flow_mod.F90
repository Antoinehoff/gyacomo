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
    INTEGER,    DIMENSION(:),   ALLOCATABLE, PUBLIC, PROTECTED :: jump_ExB     ! jump to do to shift the kx grids
    LOGICAL,    DIMENSION(:),   ALLOCATABLE, PUBLIC, PROTECTED :: shiftnow_ExB ! Indicates if there is a line to shift
    COMPLEX(xp),DIMENSION(:,:), ALLOCATABLE, PUBLIC, PROTECTED :: ExB_NL_factor! factor for nonlinear term
    COMPLEX(xp),DIMENSION(:,:), ALLOCATABLE, PUBLIC, PROTECTED :: inv_ExB_NL_factor
    LOGICAL,  PUBLIC, PROTECTED :: ExB =  .false.         ! presence of ExB background shearing rate
    REAL(xp), PUBLIC, PROTECTED :: maximal_kx, minimal_kx ! max and min kx evolved for array shifting
    ! Routines
    PUBLIC :: Setup_ExB_shear_flow, Array_shift_ExB_shear_flow, Update_ExB_shear_flow

CONTAINS

    ! Setup the variables for the ExB shear
    SUBROUTINE Setup_ExB_shear_flow(ExBrate)
        USE grid,     ONLY: Nx, local_nky, total_nky, local_nx, Ny, deltakx, deltaky,&
                            kx_max, kx_min !kyarray, kyarray_full
        USE geometry, ONLY: Cyq0_x0, C_y
        USE basic,    ONLY: speak
        USE model,    ONLY: LINEARITY
        IMPLICIT NONE
        INTEGER :: iky
        REAL(xp), INTENT(IN) :: ExBrate

        ! Setup the ExB shearing rate and aux var
        ! In GENE, there is a minus sign here...
        gamma_E = ExBrate*C_y*abs(Cyq0_x0/C_y)
        IF(abs(gamma_E) .GT. EPSILON(gamma_E)) THEN
            CALL speak('-ExB background flow detected-',2)
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

        ! Setup the jump array
        ALLOCATE(jump_ExB(local_nky))
        jump_ExB = 0

        ! Setup the shifting flag array
        ALLOCATE(shiftnow_ExB(local_nky))
        shiftnow_ExB = .FALSE.

        ! Setup nonlinear factor (McMillan 2019)
        ALLOCATE(    ExB_NL_factor(Nx,local_nky))
        ALLOCATE(inv_ExB_NL_factor(Ny/2+1,local_nx))
            ExB_NL_factor = 1._xp
        inv_ExB_NL_factor = 1._xp

        ! Setup maximal evolved modes for the shift conditions
        SELECT CASE (LINEARITY)
        CASE('linear') ! If linear we just use the max kx
            maximal_kx = kx_max
            minimal_kx = kx_min
        CASE('nonlinear') ! If NL, 2/3 rule
            maximal_kx = 2._xp/3._xp*kx_max-deltakx
            minimal_kx = 2._xp/3._xp*kx_min+deltakx
        END SELECT           

    END SUBROUTINE Setup_ExB_shear_flow

    ! update the ExB shear value for the next time step
    SUBROUTINE Update_ExB_shear_flow(step_number)
        USE basic,      ONLY: dt!,time
        USE grid,       ONLY: local_nky, total_nky, kyarray, inv_dkx, update_grids, deltaky!,kyarray_full
        USE geometry,   ONLY: gxx,gxy,gyy,inv_hatB2, evaluate_magn_curv
        USE numerics,   ONLY: evaluate_EM_op, evaluate_kernels
        ! USE model,      ONLY: LINEARITY
        USE time_integration, ONLY: c_E!, ntimelevel
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: step_number
        ! local var
        REAL(xp):: dt_sub
        INTEGER :: iky
        ! do nothing if no ExB
        IF(ExB) THEN
            IF (step_number .LT. 0) THEN ! not updated each substeps (call from control)
                dt_sub = dt
            ELSEIF(step_number .GT. 1) THEN ! updated each substeps
                dt_sub = (c_E(step_number)-c_E(step_number-1))*dt
            ELSE ! first step is at t=0
                dt_sub = 0._xp
            ENDIF
            ! Do nothing if dt is 0
            IF(dt_sub .GT. epsilon(dt_sub)) THEN 
                ! Update new shear value
                DO iky = 1,local_nky
                    !! This must be done incrementely to be able to pull it back
                    !  when a grid shift occurs                
                    sky_ExB(iky)      = sky_ExB(iky) - kyarray(iky)*gamma_E*dt_sub
                    jump_ExB(iky)     = NINT(sky_ExB(iky)*inv_dkx)
                    ! If the jump is 1 or more for a given ky, we flag the index
                    ! in shiftnow_ExB and will use it in Shift_fields to avoid
                    ! zero-shiftings that may be majoritary.
                    shiftnow_ExB(iky) = (abs(jump_ExB(iky)) .GT. 0)
                ENDDO
                ! Update the full skyExB array too
                DO iky = 1,total_nky+1
                    sky_ExB_full(iky) = sky_ExB_full(iky) - REAL(iky-1,xp)*deltaky*gamma_E*dt_sub
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
                ! IF(LINEARITY .EQ. 'nonlinear') &
                ! CALL Update_nonlinear_ExB_factors(dt_sub)
            ENDIF
        ENDIF
    END SUBROUTINE Update_ExB_shear_flow

    ! According to the current ExB shear value we update
    ! the fields by imposing a shift on kx
    SUBROUTINE Array_shift_ExB_shear_flow
        USE grid,       ONLY: local_nky, total_nky, update_grids, &
            total_nkx, deltakx, kxarray0, inv_dkx!,kx_min, kx_max
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
                    IF(jump_ExB(iky) .LT. 0) THEN
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
                    ! to monotonically travel across the kx array the indices write
                    ! 67812345 for positive shift or 54321678 for negative shift
                    DO i_ = loopstart, loopend, increment
                        ! We shift our index since kx is stored in a [0 kmax]U[kmin -dk] fashion
                        IF (i_ .LT. total_nkx/2) THEN ! go to the negative kx region
                            ikx = i_ + total_nkx/2 + 1
                        ELSE ! positive
                            ikx = i_ - total_nkx/2 + 1
                        ENDIF
                        ikx_s = ikx - jump_ExB(iky)
                        ! adjust the shift accordingly
                        IF (ikx_s .LE. 0) &
                            ikx_s = ikx_s + total_nkx
                        IF (ikx_s .GT. total_nkx) &
                            ikx_s = ikx_s - total_nkx
                        ! Then we test if the shifted modes are still contained in our resolution
                        IF ( (kxarray0(ikx)+jump_ExB(iky)*deltakx .LE. maximal_kx) .AND. &
                             (kxarray0(ikx)+jump_ExB(iky)*deltakx .GE. minimal_kx)) THEN
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

    SUBROUTINE Update_nonlinear_ExB_factors(dt_sub)
        USE grid,  ONLY: local_nky, local_nky_offset, Nx, Ny, local_nx, deltakx,&
                         local_nx_offset, deltaky, update_grids!,xarray, ikyarray, inv_ikyarray
        USE basic, ONLY: time!, dt
        ! USE time_integration, ONLY: c_E
        IMPLICIT NONE
        REAL(xp), INTENT(IN) :: dt_sub        
        INTEGER :: iky, ix
        REAL(xp):: dt_ExB, J_xp, inv_J, &
            I_xp, deltax, Lx, xval, tnow, kxExB, dkx_ExB, v_ExB
        ! aux val
        tnow   = time + dt_sub
        Lx     = 2._xp*pi/deltakx
        deltax = Lx/REAL(Nx,xp)

        DO iky = 1,local_nky ! WARNING: Local indices ky loop
            ! for readability
            ! J_xp  = ikyarray(iky+local_nky_offset)
            ! inv_J = inv_ikyarray(iky+local_nky_offset)
            J_xp   = REAL(iky+local_nky_offset-1,xp)
            IF(J_xp .GT. 0._xp) THEN
                inv_J  = 1._xp/J_xp
            ELSE
                inv_J  = 0._xp
            ENDIF
            ! compute dt factor
            dt_ExB = (tnow - t0*inv_J*ANINT(J_xp*tnow*inv_t0,xp))
            v_ExB  = -gamma_E*J_xp*deltaky
            dkx_ExB= v_ExB*tnow - ANINT(v_ExB*tnow/deltakx,xp)*deltakx
            DO ix = 1,Nx
                I_xp  = REAL(ix-1,xp)
                xval  = I_xp*deltax
                kxExB = gamma_E*J_xp*deltaky*dt_ExB
                ! kxExB = sky_ExB(iky) ! GENE
                kxExB = dkx_ExB ! home made
                ExB_NL_factor(ix,iky) = EXP(-imagu*kxExB*xval)
                ! ExB_NL_factor(ix,iky) = 1._xp + (imagu*kxExB*xval)      
            ENDDO
        ENDDO
        ! ... and the inverse
        DO iky = 1,Ny/2+1 ! WARNING: Global indices ky loop
            ! for readability
            J_xp   = REAL(iky-1,xp)
            IF(J_xp .GT. 0._xp) THEN
                inv_J  = 1._xp/J_xp
            ELSE
                inv_J  = 0._xp
            ENDIF
            ! compute dt factor
            dt_ExB = tnow - t0*inv_J*ANINT(J_xp*tnow*inv_t0,xp)
            v_ExB  = -gamma_E*J_xp*deltaky
            dkx_ExB= v_ExB*tnow - ANINT(v_ExB*tnow/deltakx,xp)*deltakx
            DO ix = 1,local_nx
                I_xp  = REAL(ix+local_nx_offset-1,xp)
                xval = I_xp*deltax
                ! kxExB = gamma_E*J_xp*deltaky*dt_ExB
                ! kxExB = sky_ExB_full(iky) ! GENE
                kxExB = dkx_ExB ! home made
                ! assemble the inverse ExB nonlin factor
                inv_ExB_NL_factor(iky,ix) = EXP(imagu*kxExB*xval)
                ! inv_ExB_NL_factor(iky,ix) = 1._xp + (imagu*kxExB*xval)
            ENDDO
        ENDDO
 END SUBROUTINE Update_nonlinear_ExB_factors
END MODULE ExB_shear_flow