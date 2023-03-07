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
    USE model,  ONLY: CLOS
    use fields, ONLY: moments
    use array,  ONLY: moments_rhs
    IMPLICIT NONE
    INTEGER :: p_int, j_int, ia, ip, ij

    CALL cpu_time(t0_adv_field)
    DO ia=ias,iae
      DO ip=ips,ipe
        p_int = parray(ip)
        DO ij=ijs,ije
          j_int = jarray(ij)
          IF((CLOS .NE. 1) .OR. (p_int+2*j_int .LE. dmax))&
          CALL advance_field(moments(ia,ip,ij,ikys:ikye,ikxs:ikxe,izs:ize,:), moments_rhs(ia,ip,ij,ikys:ikye,ikxs:ikxe,izs:ize,:))
        ENDDO
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
