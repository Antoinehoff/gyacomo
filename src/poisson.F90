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
  REAL(dp)    :: alphaD
  REAL(dp)    :: sum_kernel2_e,    sum_kernel2_i    ! Store sum Kn^2
  COMPLEX(dp) :: sum_kernel_mom_e, sum_kernel_mom_i ! Store sum Kn*Napn
  REAL(dp)    :: gammaD
  COMPLEX(dp) :: gammaD_phi
  INTEGER     :: count !! mpi integer to broadcast the electric potential at the end
  COMPLEX(dp) :: buffer(ikxs:ikxe,ikys:ikye)

  !! Poisson can be solved only for process containing ip=1
  IF ( (ips_e .EQ. 1) .AND. (ips_i .EQ. 1) ) THEN

    ! Execution time start
    CALL cpu_time(t0_poisson)

    kxloop: DO ikx = ikxs,ikxe
      kyloop: DO iky = ikys,ikye
        zloop: DO iz = izs,ize

          !!!!!!!!!!!!! Electrons sum(Kernel * Ne0n)  (skm) and sum(Kernel**2) (sk2)
          sum_kernel_mom_e = 0._dp
          sum_kernel2_e    = 0._dp
          ! loop over n only if the max polynomial degree
          DO ine=1,jmaxe+1 ! ine = n+1
            Kne = kernel_e(ine, ikx, iky)
            sum_kernel_mom_e  = sum_kernel_mom_e  + Kne * moments_e(1, ine, ikx, iky, iz, updatetlevel)
            sum_kernel2_e     = sum_kernel2_e     + Kne**2 ! ... sum recursively ...
          END DO

          !!!!!!!!!!!!!!!!! Ions sum(Kernel * Ni0n)  (skm) and sum(Kernel**2) (sk2)
            sum_kernel_mom_i = 0._dp
            sum_kernel2_i    = 0._dp
            ! loop over n only if the max polynomial degree
            DO ini=1,jmaxi+1
              Kni = kernel_i(ini, ikx, iky)
              sum_kernel_mom_i  = sum_kernel_mom_i  + Kni * moments_i(1, ini, ikx, iky, iz, updatetlevel)
              sum_kernel2_i     = sum_kernel2_i     + Kni**2 ! ... sum recursively ...
            END DO

          !!!!!!!!!!!!!!! Assembling the poisson equation !!!!!!!!!!!!!!!!!!!!!!!!!!
          alphaD   = (kxarray(ikx)**2 + kyarray(iky)**2) * lambdaD**2
          gammaD   = alphaD + qe2_taue * (1._dp - sum_kernel2_e) & ! Called Poisson_ in MOLI
                            + qi2_taui * (1._dp - sum_kernel2_i)

          gammaD_phi = q_e * sum_kernel_mom_e + q_i * sum_kernel_mom_i

          phi(ikx, iky, iz) =  gammaD_phi/gammaD

        END DO zloop
      END DO kyloop
    END DO kxloop

    ! Cancel origin singularity
    IF ((ikx_0 .GE. ikxs) .AND. (ikx_0 .LE. ikxe)) phi(ikx_0,iky_0,izs:ize) = 0._dp

  ENDIF

  ! Transfer phi to all the others process along p
  CALL manual_3D_bcast(phi(ikxs:ikxe,ikys:ikye,izs:ize))

  ! Execution time end
  CALL cpu_time(t1_poisson)
  tc_poisson = tc_poisson + (t1_poisson - t0_poisson)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE poisson
