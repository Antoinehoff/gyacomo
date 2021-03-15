SUBROUTINE compute_Sapj
  ! This routine is meant to compute the non linear term for each specie and degree
  !! In real space Sapj ~ b*(grad(phi) x grad(g)) which in moments in fourier becomes
  !! Sapj = Sum_n (ikr Kn phi)#(ikz Sum_s d_njs Naps) - (ikz Kn phi)#(ikr Sum_s d_njs Naps)
  !! where # denotes the convolution.
  USE array, ONLY : dnjs, Sepj, Sipj, kernel_i, kernel_e
  USE basic
  USE fourier
  USE fields!, ONLY : phi, moments_e, moments_i
  USE grid
  USE model
  USE prec_const
  USE time_integration!, ONLY : updatetlevel
  IMPLICIT NONE
  INCLUDE 'fftw3-mpi.f03'

  COMPLEX(dp), DIMENSION(ikrs:ikre,ikzs:ikze) :: Fr_cmpx, Gz_cmpx
  COMPLEX(dp), DIMENSION(ikrs:ikre,ikzs:ikze) :: Fz_cmpx, Gr_cmpx, F_conv_G
  REAL(dp),    DIMENSION(irs:ire,izs:ize)     :: fr_real, gz_real
  REAL(dp),    DIMENSION(irs:ire,izs:ize)     :: fz_real, gr_real, f_times_g

  INTEGER :: in, is
  REAL(dp):: kr, kz, kerneln
  LOGICAL :: COMPUTE_ONLY_ODD_P = .true.
  ! Execution time start
  CALL cpu_time(t0_Sapj)

  !!!!!!!!!!!!!!!!!!!! ELECTRON non linear term computation (Sepj)!!!!!!!!!!
  ploope: DO ip = ips_e,ipe_e ! Loop over Hermite moments

    ! we check if poly degree is even (eq to index is odd) to spare computation
    IF (MODULO(ip,2) .EQ. 1 .OR. (.NOT. COMPUTE_ONLY_ODD_P)) THEN
    jloope: DO ij = ijs_e, ije_e ! Loop over Laguerre moments

      real_data_c = 0._dp ! initialize sum over real nonlinear term

      nloope: DO in = 1,jmaxe+1 ! Loop over laguerre for the sum

        krloope: DO ikr = ikrs,ikre ! Loop over kr
          kzloope: DO ikz = ikzs,ikze ! Loop over kz
            kr     = krarray(ikr)
            kz     = kzarray(ikz)
            kerneln = kernel_e(in, ikr, ikz)

            ! First convolution terms
            Fr_cmpx(ikr,ikz) = imagu*kr* phi(ikr,ikz) * kerneln
            Fz_cmpx(ikr,ikz) = imagu*kz* phi(ikr,ikz) * kerneln
            ! Second convolution terms
            Gz_cmpx(ikr,ikz) = 0._dp ! initialization of the sum
            Gr_cmpx(ikr,ikz) = 0._dp ! initialization of the sum
            DO is = 1, MIN( in+ij-1, jmaxe+1 ) ! sum truncation on number of moments
              Gz_cmpx(ikr,ikz) = Gz_cmpx(ikr,ikz) + &
                dnjs(in,ij,is) * moments_e(ip,is,ikr,ikz,updatetlevel)
              Gr_cmpx(ikr,ikz) = Gr_cmpx(ikr,ikz) + &
                dnjs(in,ij,is) * moments_e(ip,is,ikr,ikz,updatetlevel)
            ENDDO
            Gz_cmpx(ikr,ikz) = imagu*kz*Gz_cmpx(ikr,ikz)
            Gr_cmpx(ikr,ikz) = imagu*kr*Gr_cmpx(ikr,ikz)
          ENDDO kzloope
        ENDDO krloope

        ! First term drphi x dzf
        DO ikr = ikrs, ikre
          DO ikz = ikzs, ikze
            cmpx_data_f(ikz,ikr-local_nkr_offset) = Fr_cmpx(ikr,ikz)*AA_r(ikr)*AA_z(ikz) !Anti aliasing filter
            cmpx_data_g(ikz,ikr-local_nkr_offset) = Gz_cmpx(ikr,ikz)*AA_r(ikr)*AA_z(ikz) !Anti aliasing filter
          ENDDO
        ENDDO

        call fftw_mpi_execute_dft_c2r(planb, cmpx_data_f, real_data_f)
        call fftw_mpi_execute_dft_c2r(planb, cmpx_data_g, real_data_g)

        real_data_c = real_data_c + real_data_f/Nz/Nr  * real_data_g/Nz/Nr

        ! Second term -dzphi x drf
        DO ikr = ikrs, ikre
          DO ikz = ikzs, ikze
            cmpx_data_f(ikz,ikr-local_nkr_offset) = Fz_cmpx(ikr,ikz)*AA_r(ikr)*AA_z(ikz) !Anti aliasing filter
            cmpx_data_g(ikz,ikr-local_nkr_offset) = Gr_cmpx(ikr,ikz)*AA_r(ikr)*AA_z(ikz) !Anti aliasing filter
          ENDDO
        ENDDO

        call fftw_mpi_execute_dft_c2r(planb, cmpx_data_f, real_data_f)
        call fftw_mpi_execute_dft_c2r(planb, cmpx_data_g, real_data_g)

        real_data_c = real_data_c - real_data_f/Nz/Nr  * real_data_g/Nz/Nr

      ENDDO nloope

      ! Put the real nonlinear product into k-space
      call fftw_mpi_execute_dft_r2c(planf, real_data_c, cmpx_data_c)

      ! Retrieve convolution in input format
      DO ikr = ikrs, ikre
        DO ikz = ikzs, ikze
          Sepj(ip,ij,ikr,ikz) = cmpx_data_c(ikz,ikr-local_nkr_offset)*AA_r(ikr)*AA_z(ikz) !Anti aliasing filter
        ENDDO
      ENDDO
    ENDDO jloope

    ELSE
      ! Cancel the non lin term if we are dealing with odd Hermite degree
      Sepj(ip,:,:,:) = 0._dp
    ENDIF
  ENDDO ploope
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!! ION non linear term computation (Sipj)!!!!!!!!!!
  ploopi: DO ip = ips_i,ipe_i ! Loop over Hermite moments

    ! we check if poly degree is even (eq to index is odd) to spare computation
    IF (MODULO(ip,2) .EQ. 1 .OR. (.NOT. COMPUTE_ONLY_ODD_P)) THEN

    jloopi: DO ij = ijs_i, ije_i ! Loop over Laguerre moments
      real_data_c = 0._dp ! initialize sum over real nonlinear term

      nloopi: DO in = 1,jmaxi+1 ! Loop over laguerre for the sum

        krloopi: DO ikr = ikrs,ikre ! Loop over kr
          kzloopi: DO ikz = ikzs,ikze ! Loop over kz
            kr      = krarray(ikr)
            kz      = kzarray(ikz)
            kerneln = kernel_i(in, ikr, ikz)

            ! First convolution terms
            Fr_cmpx(ikr,ikz) = imagu*kr* phi(ikr,ikz) * kerneln
            Fz_cmpx(ikr,ikz) = imagu*kz* phi(ikr,ikz) * kerneln
            ! Second convolution terms
            Gz_cmpx(ikr,ikz) = 0._dp ! initialization of the sum
            Gr_cmpx(ikr,ikz) = 0._dp ! initialization of the sum
            DO is = 1, MIN( in+ij-1, jmaxi+1 ) ! sum truncation on number of moments
              Gz_cmpx(ikr,ikz) = Gz_cmpx(ikr,ikz) + &
                dnjs(in,ij,is) * moments_i(ip,is,ikr,ikz,updatetlevel)
              Gr_cmpx(ikr,ikz) = Gr_cmpx(ikr,ikz) + &
                dnjs(in,ij,is) * moments_i(ip,is,ikr,ikz,updatetlevel)
            ENDDO
            Gz_cmpx(ikr,ikz) = imagu*kz*Gz_cmpx(ikr,ikz)
            Gr_cmpx(ikr,ikz) = imagu*kr*Gr_cmpx(ikr,ikz)
          ENDDO kzloopi
        ENDDO krloopi

        ! First term drphi x dzf
        DO ikr = ikrs, ikre
          DO ikz = ikzs, ikze
            cmpx_data_f(ikz,ikr-local_nkr_offset) = Fr_cmpx(ikr,ikz)*AA_r(ikr)*AA_z(ikz)
            cmpx_data_g(ikz,ikr-local_nkr_offset) = Gz_cmpx(ikr,ikz)*AA_r(ikr)*AA_z(ikz)
          ENDDO
        ENDDO

        call fftw_mpi_execute_dft_c2r(planb, cmpx_data_f, real_data_f)
        call fftw_mpi_execute_dft_c2r(planb, cmpx_data_g, real_data_g)

        real_data_c = real_data_c + real_data_f/Nz/Nr  * real_data_g/Nz/Nr

        ! Second term -dzphi x drf
        DO ikr = ikrs, ikre
          DO ikz = ikzs, ikze
            cmpx_data_f(ikz,ikr-local_nkr_offset) = Fz_cmpx(ikr,ikz)*AA_r(ikr)*AA_z(ikz)
            cmpx_data_g(ikz,ikr-local_nkr_offset) = Gr_cmpx(ikr,ikz)*AA_r(ikr)*AA_z(ikz)
          ENDDO
        ENDDO

        call fftw_mpi_execute_dft_c2r(planb, cmpx_data_f, real_data_f)
        call fftw_mpi_execute_dft_c2r(planb, cmpx_data_g, real_data_g)

        real_data_c = real_data_c - real_data_f/Nz/Nr  * real_data_g/Nz/Nr

      ENDDO nloopi

      ! Put the real nonlinear product into k-space
      call fftw_mpi_execute_dft_r2c(planf, real_data_c, cmpx_data_c)

      ! Retrieve convolution in input format
      DO ikr = ikrs, ikre
        DO ikz = ikzs, ikze
          Sipj(ip,ij,ikr,ikz) = cmpx_data_c(ikz,ikr-local_nkr_offset)*AA_r(ikr)*AA_z(ikz)
        ENDDO
      ENDDO

    ENDDO jloopi
    ELSE
      ! Cancel the non lin term if we are dealing with odd Hermite degree
      Sipj(ip,:,:,:) = 0._dp
    ENDIF
  ENDDO ploopi
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Execution time END
  CALL cpu_time(t1_Sapj)
  tc_Sapj = tc_Sapj + (t1_Sapj - t0_Sapj)

END SUBROUTINE compute_Sapj
