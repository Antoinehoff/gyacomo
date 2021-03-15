SUBROUTINE poisson
  ! Solve poisson equation to get phi

  USE basic
  USE time_integration, ONLY: updatetlevel
  USE array
  USE fields
  USE grid
  use model, ONLY : qe2_taue, qi2_taui, q_e, q_i, lambdaD

  USE prec_const
  IMPLICIT NONE

  INTEGER     :: ini,ine, i_, world_rank, world_size, root_bcast
  REAL(dp)    :: Kne, Kni ! sub kernel factor for recursive build
  REAL(dp)    :: alphaD
  REAL(dp)    :: sum_kernel2_e,    sum_kernel2_i    ! Store sum Kn^2
  COMPLEX(dp) :: sum_kernel_mom_e, sum_kernel_mom_i ! Store sum Kn*Napn
  REAL(dp)    :: gammaD
  COMPLEX(dp) :: gammaD_phi
  INTEGER     :: count !! mpi integer to broadcast the electric potential at the end
  COMPLEX(dp) :: buffer(ikrs:ikre,ikzs:ikze)

  !! Poisson can be solved only for process containing ip=1
  IF ( (ips_e .EQ. 1) .AND. (ips_i .EQ. 1) ) THEN

    ! Execution time start
    CALL cpu_time(t0_poisson)

    DO ikr=ikrs,ikre
      DO ikz=ikzs,ikze

        !!!!!!!!!!!!! Electrons sum(Kernel * Ne0n)  (skm) and sum(Kernel**2) (sk2)
        sum_kernel_mom_e = 0._dp
        sum_kernel2_e    = 0._dp
        ! loop over n only if the max polynomial degree
        DO ine=1,jmaxe+1 ! ine = n+1
          Kne = kernel_e(ine, ikr, ikz)
          sum_kernel_mom_e  = sum_kernel_mom_e  + Kne * moments_e(1, ine, ikr, ikz, updatetlevel)
          sum_kernel2_e     = sum_kernel2_e     + Kne**2 ! ... sum recursively ...
        END DO

        !!!!!!!!!!!!!!!!! Ions sum(Kernel * Ni0n)  (skm) and sum(Kernel**2) (sk2)
          sum_kernel_mom_i = 0
          sum_kernel2_i = 0
          ! loop over n only if the max polynomial degree
          DO ini=1,jmaxi+1
            Kni = kernel_i(ini, ikr, ikz)
            sum_kernel_mom_i  = sum_kernel_mom_i  + Kni * moments_i(1, ini, ikr, ikz, updatetlevel)
            sum_kernel2_i     = sum_kernel2_i     + Kni**2 ! ... sum recursively ...
          END DO

        !!!!!!!!!!!!!!! Assembling the poisson equation !!!!!!!!!!!!!!!!!!!!!!!!!!
        alphaD   = (krarray(ikr)**2 + kzarray(ikz)**2) * lambdaD**2
        gammaD   = alphaD + qe2_taue * (1._dp - sum_kernel2_e) & ! Called Poisson_ in MOLI
                          + qi2_taui * (1._dp - sum_kernel2_i)

        gammaD_phi = q_e * sum_kernel_mom_e + q_i * sum_kernel_mom_i

        phi(ikr, ikz) =  gammaD_phi/gammaD

      END DO
    END DO

    ! Cancel origin singularity
    IF ((ikr_0 .GE. ikrs) .AND. (ikr_0 .LE. ikre)) phi(ikr_0,ikz_0) = 0

    ! Execution time end
    CALL cpu_time(t1_poisson)
    tc_poisson = tc_poisson + (t1_poisson - t0_poisson)

  ENDIF

  CALL cpu_time(t0_comm)

  root_bcast = 0 ! Proc zero computes phi for every p

  !!!!! This is a manual way to do MPI_BCAST !!!!!!!!!!!
  CALL MPI_COMM_RANK(commp,world_rank,ierr)
  CALL MPI_COMM_SIZE(commp,world_size,ierr)

  IF (world_size .GT. 1) THEN
    !! Broadcast phi to the other processes on the same k range (communicator along p)
    IF (world_rank .EQ. root_bcast) THEN
      ! Fill the buffer
      DO ikr = ikrs,ikre
        DO ikz = ikzs,ikze
          buffer(ikr,ikz) = phi(ikr,ikz)
        ENDDO
      ENDDO
      ! Send it to all the other processes
      DO i_ = 0,num_procs_p-1
        IF (i_ .NE. world_rank) &
        CALL MPI_SEND(buffer, local_nkr * nkz , MPI_DOUBLE_COMPLEX, i_, 0, commp, ierr)
      ENDDO
    ELSE
      ! Recieve buffer from root
      CALL MPI_RECV(buffer, local_nkr * nkz , MPI_DOUBLE_COMPLEX, root_bcast, 0, commp, MPI_STATUS_IGNORE, ierr)
      ! Write it in phi
      DO ikr = ikrs,ikre
        DO ikz = ikzs,ikze
          phi(ikr,ikz) = buffer(ikr,ikz)
        ENDDO
      ENDDO
    ENDIF
  ENDIF
  
  CALL cpu_time(t1_comm)
  tc_comm = tc_comm + (t1_comm - t0_comm)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE poisson
