SUBROUTINE poisson
  ! Solve poisson equation to get phi

  USE basic
  USE time_integration, ONLY: updatetlevel
  USE array
  USE fields
  USE grid
  USE utility
  use model, ONLY : qe2_taue, qi2_taui, q_e, q_i, lambdaD

  USE prec_const
  IMPLICIT NONE

  INTEGER     :: ini,ine, i_, root_bcast
  REAL(dp)    :: Kne, Kni ! sub kernel factor for recursive build
  REAL(dp)    :: polarisation    ! sum_a(Z_a^2/tau_a (1-sum_n kernel_na^2))
  COMPLEX(dp) :: q_density       ! charge density sum_a q_a n_a
  REAL(dp)    :: gammaD
  COMPLEX(dp) :: gammaD_phi
  INTEGER     :: count !! mpi integer to broadcast the electric potential at the end
  COMPLEX(dp) :: buffer(ikxs:ikxe,ikys:ikye)

  !! Poisson can be solved only for process containing ip=1
  IF ( (ips_e .EQ. ip0_e) .AND. (ips_i .EQ. ip0_i) ) THEN

    ! Execution time start
    CALL cpu_time(t0_poisson)

    kxloop: DO ikx = ikxs,ikxe
      kyloop: DO iky = ikys,ikye
        zloop: DO iz = izs,ize

          q_density      = 0._dp
          polarisation   = 0._dp
          !!!!!!!!!!!!! Electron contribution
          ! loop over n only if the max polynomial degree
          DO ine=1,jmaxe+1 ! ine = n+1
            Kne           = kernel_e(ine,ikx,iky,iz)
            q_density     = q_density     + q_e*Kne*moments_e(ip0_e, ine, ikx, iky, iz, updatetlevel)
            polarisation  = polarisation  + qe2_taue*Kne**2 ! ... sum recursively ...
          END DO

          !!!!!!!!!!!!!!!!! Ions contribution
          ! loop over n only if the max polynomial degree
          DO ini=1,jmaxi+1
            Kni           = kernel_i(ini,ikx,iky,iz)
            q_density     = q_density     + q_i*Kni*moments_i(ip0_i, ini, ikx, iky, iz, updatetlevel)
            polarisation  = polarisation  + qi2_taui*Kni**2 ! ... sum recursively ...
          END DO

          !!!!!!!!!!!!!!! Assembling the poisson equation !!!!!!!!!!!!!!!!!!!!!!!!!!
          phi(ikx, iky, iz) =  q_density/(qe2_taue + qi2_taui - polarisation)

        END DO zloop
      END DO kyloop
    END DO kxloop

    ! Cancel origin singularity
    IF (contains_kx0 .AND. contains_ky0) phi(ikx_0,iky_0,izs:ize) = 0._dp

  ENDIF

  ! Transfer phi to all the others process along p
  CALL manual_3D_bcast(phi(ikxs:ikxe,ikys:ikye,izs:ize))

  ! Execution time end
  CALL cpu_time(t1_poisson)
  tc_poisson = tc_poisson + (t1_poisson - t0_poisson)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE poisson
