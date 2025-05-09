MODULE advance_field_routine

USE prec_const
implicit none

CONTAINS

  SUBROUTINE advance_time_level
    USE time_integration, ONLY :set_updatetlevel, updatetlevel, ntimelevel
    IMPLICIT NONE
    CALL set_updatetlevel(mod(updatetlevel,ntimelevel)+1)
  END SUBROUTINE advance_time_level

  SUBROUTINE advance_moments

    USE basic,  ONLY: dt
    USE grid,   ONLY:local_na,local_np,local_nj,local_nky,local_nkx,local_nz,&
                     ngp_o2, ngj_o2, ngz_o2
    USE closure,ONLY: evolve_mom
    use fields, ONLY: moments
    use array,  ONLY: moments_rhs
    USE time_integration, ONLY: updatetlevel, A_E, b_E, ntimelevel
    IMPLICIT NONE
    INTEGER :: ia, ip, ij, ikx,iky,iz, istage, ipi, iji, izi
    SELECT CASE (updatetlevel)
    CASE(1)
      DO istage=1,ntimelevel
      DO iz    =1,local_nz
        izi = iz+ngz_o2
      DO ikx   =1,local_nkx
      DO iky   =1,local_nky
      DO ij    =1,local_nj
        iji = ij+ngj_o2
      DO ip    =1,local_np
        ipi = ip+ngp_o2
      DO ia    =1,local_na
        IF( evolve_mom(ipi,iji) )&
        moments(ia,ipi,iji,iky,ikx,izi,1) = moments(ia,ipi,iji,iky,ikx,izi,1) &
               + dt*b_E(istage,1)*moments_rhs(ia,ip,ij,iky,ikx,iz,istage)
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
        izi = iz+ngz_o2
      DO ikx   =1,local_nkx
      DO iky   =1,local_nky
      DO ij    =1,local_nj
        iji = ij+ngj_o2
      DO ip    =1,local_np
        ipi = ip+ngp_o2
      DO ia    =1,local_na
        IF( evolve_mom(ipi,iji) )&
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
  END SUBROUTINE advance_moments
END MODULE advance_field_routine
