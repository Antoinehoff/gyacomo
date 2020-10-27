MODULE advance_field_routine

USE prec_const
implicit none

CONTAINS

  SUBROUTINE advance_time_level

    USE basic
    USE time_integration
    use prec_const
    IMPLICIT NONE

    call set_updatetlevel(mod(updatetlevel,ntimelevel)+1)

  END SUBROUTINE advance_time_level



  SUBROUTINE advance_field( f, f_rhs )

    USE basic
    USE time_integration
    USE array
    USE grid
    use prec_const
    IMPLICIT NONE

    COMPLEX(dp), DIMENSION ( ikrs:ikre, ikzs:ikze, ntimelevel ) :: f
    COMPLEX(dp), DIMENSION ( ikrs:ikre, ikzs:ikze, ntimelevel ) :: f_rhs
    REAL(dp)    :: error
    INTEGER     :: istage

    SELECT CASE (updatetlevel)

    CASE(1)
      SELECT CASE (numerical_scheme)

      CASE ('DOPRI5_ADAPT')
        error = 0._dp
        DO ikz=ikzs,ikze
          DO ikr=ikrs,ikre
            fs = f(ikr,ikz,1)
            DO istage=1,ntimelevel
              f(ikr,ikz,1) = f(ikr,ikz,1) + dt*b_E(istage)*f_rhs(ikr,ikz,istage)
              fs = fs + dt*b_Es(istage)*f_rhs(ikr,ikz,istage)
            END DO
            IF ( ABS(f(ikr,ikz,1) - fs) .GT. error ) THEN
              error = ABS(f(ikr,ikz,1) - fs)
            ENDIF
          END DO
        END DO
        IF (error > TOL)
      CASE DEFAULT
        DO ikz=ikzs,ikze
          DO ikr=ikrs,ikre
            fs = f(ikr,ikz,1)
            DO istage=1,ntimelevel
            f(ikr,ikz,1) = f(ikr,ikz,1) + dt*b_E(istage)*f_rhs(ikr,ikz,istage)
            END DO
          END DO
        END DO
      END SELECT

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
