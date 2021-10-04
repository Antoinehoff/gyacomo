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

    USE basic
    USE time_integration
    USE grid
    use prec_const
    use fields, ONLY: moments_e,     moments_i
    use array,  ONLY: moments_rhs_e, moments_rhs_i
    IMPLICIT NONE

    DO ip=ips_i,ipe_i
      DO ij=ijs_i,ije_i
        CALL advance_field(moments_i(ip,ij,:,:,:,:), moments_rhs_i(ip,ij,:,:,:,:))
      ENDDO
    ENDDO
    DO ip=ips_e,ipe_e
      DO ij=ijs_e,ije_e
        CALL advance_field(moments_e(ip,ij,:,:,:,:), moments_rhs_e(ip,ij,:,:,:,:))
      ENDDO
    ENDDO
    ! Execution time end
    CALL cpu_time(t1_adv_field)
    tc_adv_field = tc_adv_field + (t1_adv_field - t0_adv_field)
  END SUBROUTINE advance_moments


  SUBROUTINE advance_field( f, f_rhs )

    USE basic
    USE time_integration
    USE array
    USE grid
    use prec_const
    use initial_par, ONLY: WIPE_ZF
    IMPLICIT NONE

    COMPLEX(dp), DIMENSION ( ikxs:ikxe, ikys:ikye, izs:ize, ntimelevel ) :: f
    COMPLEX(dp), DIMENSION ( ikxs:ikxe, ikys:ikye, izs:ize, ntimelevel ) :: f_rhs
    INTEGER :: istage

    SELECT CASE (updatetlevel)
    CASE(1)
      DO iky=ikys,ikye
          DO ikx=ikxs,ikxe
            DO iz=izs,ize
              DO istage=1,ntimelevel
                f(ikx,iky,iz,1) = f(ikx,iky,iz,1) + dt*b_E(istage)*f_rhs(ikx,iky,iz,istage)
              END DO
            END DO
          END DO
      END DO
    CASE DEFAULT
      DO iky=ikys,ikye
          DO ikx=ikxs,ikxe
            DO iz=izs,ize
              f(ikx,iky,iz,updatetlevel) = f(ikx,iky,iz,1);
              DO istage=1,updatetlevel-1
                f(ikx,iky,iz,updatetlevel) = f(ikx,iky,iz,updatetlevel) + &
                                      dt*A_E(updatetlevel,istage)*f_rhs(ikx,iky,iz,istage)
              END DO
            END DO
          END DO
      END DO
    END SELECT
  END SUBROUTINE advance_field

END MODULE advance_field_routine
