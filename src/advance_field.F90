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

    INTEGER :: istage
    ! Execution time start
    CALL cpu_time(t0_adv_field)
    SELECT CASE (updatetlevel)
    CASE(1)
      DO istage=1,ntimelevel
        moments_e(ips_e:ipe_e,ijs_e:ije_e,:,:,1) = moments_e(ips_e:ipe_e,ijs_e:ije_e,:,:,1) &
                    + dt*b_E(istage)*moments_rhs_e(ips_e:ipe_e,ijs_e:ije_e,:,:,istage)
      END DO
      ! Advance ions
      DO istage=1,ntimelevel
        moments_i(ips_i:ipe_i,ijs_i:ije_i,:,:,1) = moments_i(ips_i:ipe_i,ijs_i:ije_i,:,:,1) &
                    + dt*b_E(istage)*moments_rhs_i(ips_i:ipe_i,ijs_i:ije_i,:,:,istage)
      END DO

    CASE DEFAULT
      ! Advance electrons
      moments_e(ips_e:ipe_e,ijs_e:ije_e,:,:,updatetlevel) = moments_e(ips_e:ipe_e,ijs_e:ije_e,:,:,1);
      DO istage=1,updatetlevel-1
        moments_e(ips_e:ipe_e,ijs_e:ije_e,:,:,updatetlevel) = moments_e(ips_e:ipe_e,ijs_e:ije_e,:,:,updatetlevel) &
              + dt*A_E(updatetlevel,istage)*moments_rhs_e(ips_e:ipe_e,ijs_e:ije_e,:,:,istage)
      END DO
      ! Advance ions
      moments_i(ips_i:ipe_i,ijs_i:ije_i,:,:,updatetlevel) = moments_i(ips_i:ipe_i,ijs_i:ije_i,:,:,1);
      DO istage=1,updatetlevel-1
        moments_i(ips_i:ipe_i,ijs_i:ije_i,:,:,updatetlevel) = moments_i(ips_i:ipe_i,ijs_i:ije_i,:,:,updatetlevel) &
              + dt*A_E(updatetlevel,istage)*moments_rhs_i(ips_i:ipe_i,ijs_i:ije_i,:,:,istage)
      END DO
    END SELECT
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

    COMPLEX(dp), DIMENSION ( ikrs:ikre, ikzs:ikze, ntimelevel ) :: f
    COMPLEX(dp), DIMENSION ( ikrs:ikre, ikzs:ikze, ntimelevel ) :: f_rhs
    INTEGER :: istage

    SELECT CASE (updatetlevel)
    CASE(1)
      DO ikz=ikzs,ikze
        DO ikr=ikrs,ikre
          DO istage=1,ntimelevel
            f(ikr,ikz,1) = f(ikr,ikz,1) + dt*b_E(istage)*f_rhs(ikr,ikz,istage)
          END DO
        END DO
      END DO
    CASE DEFAULT
      DO ikz=ikzs,ikze
        DO ikr=ikrs,ikre
          f(ikr,ikz,updatetlevel) = f(ikr,ikz,1);
          DO istage=1,updatetlevel-1
            f(ikr,ikz,updatetlevel) = f(ikr,ikz,updatetlevel) + &
                                  dt*A_E(updatetlevel,istage)*f_rhs(ikr,ikz,istage)
          END DO
        END DO
      END DO
    END SELECT
  END SUBROUTINE advance_field

END MODULE advance_field_routine
