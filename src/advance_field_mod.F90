MODULE advance_field_routine

USE prec_const
implicit none

CONTAINS

  SUBROUTINE advance_time_level

    USE basic
    USE time_integration
    use prec_const
    IMPLICIT NONE
    CALL set_updatetlevel(mod(updatetlevel,ntimelevel)+1)
  END SUBROUTINE advance_time_level

  SUBROUTINE advance_moments

    USE basic,  ONLY: t0_adv_field,t1_adv_field,tc_adv_field, dt
    USE grid,   ONLY:local_na,local_np,local_nj,local_nky,local_nkx,local_nz,&
                     ngp, ngj, ngz, dmax, parray, jarray, dmax
    USE model,  ONLY: CLOS
    use fields, ONLY: moments
    use array,  ONLY: moments_rhs
    USE time_integration, ONLY: updatetlevel, A_E, b_E, ntimelevel
    IMPLICIT NONE
    INTEGER :: ia, ip, ij, ikx,iky,iz, istage, ipi, iji, izi
    CALL cpu_time(t0_adv_field)
    SELECT CASE (updatetlevel)
    CASE(1)
      DO istage=1,ntimelevel
      DO iz    =1,local_nz
        izi = iz+ngz/2
      DO ikx   =1,local_nkx
      DO iky   =1,local_nky
      DO ij    =1,local_nj
        iji = ij+ngj/2
      DO ip    =1,local_np
        ipi = ip+ngp/2
      DO ia    =1,local_na
        IF((CLOS .NE. 1) .OR. (parray(ipi)+2*jarray(iji) .LE. dmax))&
        moments(ia,ipi,iji,iky,ikx,izi,1) = moments(ia,ipi,iji,iky,ikx,izi,1) &
                + dt*b_E(istage)*moments_rhs(ia,ip,ij,iky,ikx,iz,istage)
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
    CASE DEFAULT
      moments(:,:,:,:,:,:,updatetlevel) = moments(:,:,:,:,:,:,1);
      DO istage=1,ntimelevel-1
      DO iz    =1,local_nz
        izi = iz+ngz/2
      DO ikx   =1,local_nkx
      DO iky   =1,local_nky
      DO ij    =1,local_nj
        iji = ij+ngj/2
      DO ip    =1,local_np
        ipi = ip+ngp/2
      DO ia    =1,local_na
        IF((CLOS .NE. 1) .OR. (parray(ipi)+2*jarray(iji) .LE. dmax))&
        moments(ia,ipi,iji,iky,ikx,izi,updatetlevel) = moments(ia,ipi,iji,iky,ikx,izi,updatetlevel) + &
                          dt*A_E(updatetlevel,istage)*moments_rhs(ia,ip,ij,iky,ikx,iz,istage)
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
    END SELECT
    ! Execution time end
    CALL cpu_time(t1_adv_field)
    tc_adv_field = tc_adv_field + (t1_adv_field - t0_adv_field)
  END SUBROUTINE advance_moments
END MODULE advance_field_routine
