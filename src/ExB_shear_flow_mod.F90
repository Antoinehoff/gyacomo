MODULE ExB_shear_flow
    ! This module contains the necessary tools to implement ExB shearing flow effects.
    ! The algorithm is taken from the presentation of Hammett et al. 2006 (APS) and
    ! it the one used in GS2.
    USE prec_const, ONLY: xp

    IMPLICIT NONE
    ! Variables
    REAL(xp), PUBLIC, PROTECTED :: gamma_E = 0._xp ! ExB background shearing rate \gamma_E
    REAL(xp), DIMENSION(:), ALLOCATABLE, PUBLIC, PROTECTED :: sky_ExB      ! shift of the kx modes, kx* = kx + s(ky)
    INTEGER,  DIMENSION(:), ALLOCATABLE, PUBLIC, PROTECTED :: jump_ExB     ! jump to do to shift the kx grids
    LOGICAL,  DIMENSION(:), ALLOCATABLE, PUBLIC, PROTECTED :: shiftnow_ExB ! Indicates if there is a line to shift

    ! Routines
    PUBLIC :: Setup_ExB_shear_flow, Apply_ExB_shear_flow, Update_ExB_shear_flow

CONTAINS

    ! Setup the variables for the ExB shear
    SUBROUTINE Setup_ExB_shear_flow
        USE grid,       ONLY : local_nky
        IMPLICIT NONE

        ! Setup the ExB shift
        ALLOCATE(sky_ExB(local_nky))
        sky_ExB = 0._xp

        ! Setup the jump array and shifting flag
        ALLOCATE(jump_ExB(local_nky))
        jump_ExB     = 0
        ALLOCATE(shiftnow_ExB(local_nky))
        shiftnow_ExB = .FALSE.
    END SUBROUTINE Setup_ExB_shear_flow

    ! Update according to the current ExB shear value
    ! -the grids
    ! -the spatial operators 
    ! -the fields by imposing a shift on kx
    SUBROUTINE Apply_ExB_shear_flow
        USE basic,      ONLY: chrono_ExBs, start_chrono, stop_chrono
        USE grid,       ONLY: local_nky, kxarray, update_grids, &
            local_nkx, deltakx, kx_min, kx_max
        USE prec_const, ONLY: PI
        USE geometry,   ONLY: gxx,gxy,inv_hatB2, evaluate_magn_curv
        USE numerics,   ONLY: evaluate_EM_op, evaluate_kernels
        USE fields,     ONLY: moments, phi, psi
        IMPLICIT NONE
        ! local var
        INTEGER :: iky, ikx, ikx_s

        CALL start_chrono(chrono_ExBs)

        ! shift all fields and correct the shift value
        DO iky = 1,local_Nky
            IF(shiftnow_ExB(iky)) THEN
                print*, "SHIFT ARRAYS"
                ! shift all fields
                DO ikx = 1,local_nkx
                    ikx_s = ikx + jump_ExB(iky)
                    IF( (kxarray(iky,ikx) .GE. kx_min) .AND. (kxarray(iky,ikx) .LE. kx_max) ) THEN
                        moments(:,:,:,iky,ikx,:,:) = moments(:,:,:,iky,ikx_s,:,:)
                        phi(iky,ikx,:)             = phi(iky,ikx_s,:)
                        psi(iky,ikx,:)             = psi(iky,ikx_s,:)
                    ELSE
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
        USE basic,      ONLY: dt, chrono_ExBs, start_chrono, stop_chrono
        USE grid,       ONLY: local_nky, kyarray, inv_dkx
        USE model,      ONLY: ExBrate
        IMPLICIT NONE
        ! local var
        INTEGER :: iky
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
        ENDDO
        CALL stop_chrono(chrono_ExBs)
    END SUBROUTINE Update_ExB_shear_flow
END MODULE ExB_shear_flow