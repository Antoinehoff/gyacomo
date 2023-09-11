MODULE ExB_shear_flow
    ! This module contains the necessary tools to implement ExB shearing flow effects.
    ! The algorithm is taken from the presentation of Hammett et al. 2006 (APS) and
    ! it the one used in GS2.
    USE prec_const, ONLY: xp, imagu

    IMPLICIT NONE
    ! Variables
    REAL(xp),   PUBLIC, PROTECTED :: gamma_E = 0._xp     ! ExB background shearing rate \gamma_E
    REAL(xp),   PUBLIC, PROTECTED :: t0, inv_t0 = 0._xp  ! charact. shear time
    REAL(xp),   DIMENSION(:),   ALLOCATABLE, PUBLIC, PROTECTED :: sky_ExB      ! shift of the kx modes, kx* = kx + s(ky)
    INTEGER,    DIMENSION(:),   ALLOCATABLE, PUBLIC, PROTECTED :: jump_ExB     ! jump to do to shift the kx grids
    LOGICAL,    DIMENSION(:),   ALLOCATABLE, PUBLIC, PROTECTED :: shiftnow_ExB ! Indicates if there is a line to shift
    COMPLEX(xp),DIMENSION(:,:), ALLOCATABLE, PUBLIC, PROTECTED :: ExB_NL_factor! factor for nonlinear term
    COMPLEX(xp),DIMENSION(:,:), ALLOCATABLE, PUBLIC, PROTECTED :: inv_ExB_NL_factor
    ! Routines
    PUBLIC :: Setup_ExB_shear_flow, Apply_ExB_shear_flow, Update_ExB_shear_flow

CONTAINS

    ! Setup the variables for the ExB shear
    SUBROUTINE Setup_ExB_shear_flow
        USE grid,  ONLY : total_nkx, local_nky, deltakx, deltaky
        USE model, ONLY : ExBrate  
        IMPLICIT NONE

        ! Setup the ExB shift
        ALLOCATE(sky_ExB(local_nky))
        sky_ExB = 0._xp

        ! Setup the jump array and shifting flag
        ALLOCATE(jump_ExB(local_nky))
        jump_ExB     = 0
        ALLOCATE(shiftnow_ExB(local_nky))
        shiftnow_ExB = .FALSE.

        ! Setup nonlinear factor
        ALLOCATE(    ExB_NL_factor(total_nkx,local_nky))
        ALLOCATE(inv_ExB_NL_factor(total_nkx,local_nky))
            ExB_NL_factor = 1._xp
        inv_ExB_NL_factor = 1._xp
        IF(ExBrate .NE. 0) THEN
            t0     = deltakx/deltaky/ExBrate
            inv_t0 = 1._xp/t0
        ELSE ! avoid 1/0 division (t0 is killed anyway in this case)
            t0     = 0._xp
            inv_t0 = 0._xp
        ENDIF

    END SUBROUTINE Setup_ExB_shear_flow

    ! Update according to the current ExB shear value
    ! -the grids
    ! -the spatial operators 
    ! -the fields by imposing a shift on kx
    SUBROUTINE Apply_ExB_shear_flow
        USE basic,      ONLY: chrono_ExBs, start_chrono, stop_chrono
        USE grid,       ONLY: local_nky, kxarray, update_grids, &
            total_nkx, deltakx, kx_min, kx_max
        USE prec_const, ONLY: PI
        USE geometry,   ONLY: gxx,gxy,inv_hatB2, evaluate_magn_curv
        USE numerics,   ONLY: evaluate_EM_op, evaluate_kernels
        USE fields,     ONLY: moments, phi, psi
        IMPLICIT NONE
        ! local var
        INTEGER :: iky, ikx, ikx_s, i_, loopstart, loopend, increment

        CALL start_chrono(chrono_ExBs)

        ! shift all fields and correct the shift value
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
                    ! We test if the shifted modes are still in contained in our resolution
                    ! IF( (kxarray(iky,ikx)-sky_ExB(iky) .GE. kx_min) .AND. (kxarray(iky,ikx)-sky_ExB(iky) .LE. kx_max) ) THEN
                    IF ( ((ikx_s .GT. 0 ) .AND. (ikx_s .LE. total_nkx )) .AND. &
                         (((ikx   .LE. (total_nkx/2+1)) .AND. (ikx_s .LE. (total_nkx/2+1))) .OR. &
                          ((ikx   .GT. (total_nkx/2+1)) .AND. (ikx_s .GT. (total_nkx/2+1)))) ) THEN
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
            ENDIF
        ENDDO
        ! After shifting and correction, we update the operators and grids
        !   update the grids  
        CALL update_grids(sky_ExB,gxx,gxy,inv_hatB2)

        !   update the EM op., the kernels and the curvature op.
        CALL evaluate_kernels
        CALL evaluate_EM_op
        CALL evaluate_magn_curv

        CALL stop_chrono(chrono_ExBs)
    END SUBROUTINE Apply_ExB_shear_flow

    ! update the ExB shear value for the next time step
    SUBROUTINE Update_ExB_shear_flow
        USE basic,      ONLY: dt, time, chrono_ExBs, start_chrono, stop_chrono
        USE grid,       ONLY: local_nky, kyarray, inv_dkx, xarray,&
                              local_nkx, ikyarray, inv_ikyarray, deltakx, deltaky, deltax
        USE model,      ONLY: ExBrate
        IMPLICIT NONE
        ! local var
        INTEGER :: iky, ix
        REAL(xp):: dtExBshear
        CALL start_chrono(chrono_ExBs)
        ! update the ExB shift, jumps and flags
        shiftnow_ExB = .FALSE.
        DO iky = 1,local_Nky
            sky_ExB(iky)      = sky_ExB(iky) - kyarray(iky)*ExBrate*dt
            jump_ExB(iky)     = NINT(sky_ExB(iky)*inv_dkx)
            ! If the jump is 1 or more for a given ky, we flag the index
            ! in shiftnow_ExB and will use it in Shift_fields to avoid
            ! zero-shiftings that may be majoritary.
            shiftnow_ExB(iky) = (abs(jump_ExB(iky)) .GT. 0)
            ! Update the ExB nonlinear factor
            dtExBshear = time - t0*inv_ikyarray(iky)*ANINT(ikyarray(iky)*time*inv_t0,xp)
            DO ix = 1,local_nkx
                ExB_NL_factor(ix,iky) = EXP(-imagu*xarray(ix)*ExBrate*ikyarray(iky)*dtExBshear)
            inv_ExB_NL_factor(ix,iky) = 1._xp/ExB_NL_factor(ix,iky)
            ENDDO
        ENDDO
        CALL stop_chrono(chrono_ExBs)
    END SUBROUTINE Update_ExB_shear_flow
END MODULE ExB_shear_flow