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
    USE model,  ONLY: CLOS, KIN_E
    use fields, ONLY: moments_e,     moments_i
    use array,  ONLY: moments_rhs_e, moments_rhs_i
    IMPLICIT NONE
    INTEGER :: p_int, j_int

    CALL cpu_time(t0_adv_field)
    IF(KIN_E) THEN
      DO ip=ips_e,ipe_e
        p_int = parray_e(ip)
        DO ij=ijs_e,ije_e
          j_int = jarray_e(ij)
          IF((CLOS .NE. 1) .OR. (p_int+2*j_int .LE. dmaxe))&
          CALL advance_field(moments_e(ip,ij,ikys:ikye,ikxs:ikxe,izs:ize,:), moments_rhs_e(ip,ij,ikys:ikye,ikxs:ikxe,izs:ize,:))
        ENDDO
      ENDDO
    ENDIF
    DO ip=ips_i,ipe_i
      p_int = parray_i(ip)
      DO ij=ijs_i,ije_i
        j_int = jarray_i(ij)
        IF((CLOS .NE. 1) .OR. (p_int+2*j_int .LE. dmaxi))&
        CALL advance_field(moments_i(ip,ij,ikys:ikye,ikxs:ikxe,izs:ize,:), moments_rhs_i(ip,ij,ikys:ikye,ikxs:ikxe,izs:ize,:))
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
    IMPLICIT NONE

    COMPLEX(dp), DIMENSION ( ikys:ikye, ikxs:ikxe, izs:ize, ntimelevel ) :: f
    COMPLEX(dp), DIMENSION ( ikys:ikye, ikxs:ikxe, izs:ize, ntimelevel ) :: f_rhs
    INTEGER :: istage

    SELECT CASE (updatetlevel)
    CASE(1)
        DO istage=1,ntimelevel
          f(ikys:ikye,ikxs:ikxe,izs:ize,1) = f(ikys:ikye,ikxs:ikxe,izs:ize,1) &
                   + dt*b_E(istage)*f_rhs(ikys:ikye,ikxs:ikxe,izs:ize,istage)
        END DO
    CASE DEFAULT
        f(ikys:ikye,ikxs:ikxe,izs:ize,updatetlevel) = f(ikys:ikye,ikxs:ikxe,izs:ize,1);
        DO istage=1,updatetlevel-1
          f(ikys:ikye,ikxs:ikxe,izs:ize,updatetlevel) = f(ikys:ikye,ikxs:ikxe,izs:ize,updatetlevel) + &
                            dt*A_E(updatetlevel,istage)*f_rhs(ikys:ikye,ikxs:ikxe,izs:ize,istage)
        END DO
    END SELECT
  END SUBROUTINE advance_field

END MODULE advance_field_routine
