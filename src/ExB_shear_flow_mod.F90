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
    PUBLIC :: Setup_ExB_shear_flow, Update_ExB_shear_flow

CONTAINS

    ! Setup the variables for the ExB shear
    SUBROUTINE Setup_ExB_shear_flow(g_E)
        USE basic,      ONLY : time
        USE grid,       ONLY : local_nky, kyarray, Lx
        USE prec_const, ONLY : PI
        IMPLICIT NONE
        REAL(xp), INTENT(IN) :: g_E ! Input shearing rate (comes from model module)
        ! local var
        INTEGER :: iky
        REAL(xp):: inv_dkx

        ! Store the shearing rate in this module
        gamma_E = g_E

        ! Setup the ExB shift
        ALLOCATE(sky_ExB(local_nky))
        sky_ExB = 0._xp

        ! Setup the jump array and shifting flag
        ALLOCATE(shiftnow_ExB(local_nky))
        shiftnow_ExB = .FALSE.
        ALLOCATE(jump_ExB(local_nky))
        jump_ExB     = 0
    END SUBROUTINE Setup_ExB_shear_flow

    ! Update according to the current ExB shear value
    ! -the grids
    ! -the spatial operators 
    ! -the fields by imposing a shift on kx
    ! Then update the ExB shear value for the next time step
    SUBROUTINE Update_ExB_shear_flow
        USE basic,      ONLY : dt
        USE grid,       ONLY : local_nky, kyarray, update_grids
        USE prec_const, ONLY : PI
        USE geometry,   ONLY : gxx,gxy,gyy,inv_hatB2, evaluate_magn_curv
        USE numerics,   ONLY : evaluate_EM_op, evaluate_kernels
        IMPLICIT NONE
        ! local var
        INTEGER :: iky
        REAL(xp):: inv_dkx

        ! update the grids
        CALL update_grids(sky_ExB(iky),gxx,gxy,gyy,inv_hatB2)

        ! update the EM operators and the kernels
        CALL evaluate_kernels
        CALL evaluate_EM_op

        ! shift all fields and correct the shift value
        CALL Shift_fields

        ! update the ExB shift
        shiftnow_ExB = .FALSE.
        DO iky = 1,local_Nky
            sky_ExB(iky)      = sky_ExB(iky) - kyarray(iky)*gamma_E*dt
            jump_ExB(iky)     = NINT(sky_ExB(iky)*inv_dkx)
            shiftnow_ExB(iky) = (abs(jump_ExB(iky)) .GT. 0)
        ENDDO
        CONTAINS

        SUBROUTINE Shift_fields
            USE grid,  ONLY: local_nky, kxarray, kyarray, local_nkx, deltakx, &
                             kx_min, kx_max
            USE fields,ONLY: moments, phi, psi
            IMPLICIT NONE
            ! local var
            INTEGER :: iky, ikx, ikx_s
            REAL(xp):: inv_dkx
            DO iky = 1,local_Nky
                IF(shiftnow_ExB(iky)) THEN
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
        END SUBROUTINE Shift_fields
    END SUBROUTINE Update_ExB_shear_flow


END MODULE ExB_shear_flow