SUBROUTINE compute_Sapj
  ! This routine is meant to compute the non linear term for each specie and degree
  !! In real space Sapj ~ b*(grad(phi) x grad(g)) which in moments in fourier becomes
  !! Sapj = Sum_n (ikx Kn phi)#(iky Sum_s d_njs Naps) - (iky Kn phi)#(ikx Sum_s d_njs Naps)
  !! where # denotes the convolution.
  USE array, ONLY : dnjs, Sepj, Sipj, kernel_i, kernel_e,&
                    moments_e_ZF, moments_i_ZF, phi_ZF
  USE initial_par, ONLY : WIPE_ZF
  USE basic
  USE fourier
  USE fields, ONLY : phi, moments_e, moments_i
  USE grid
  USE model
  USE prec_const
  USE time_integration, ONLY : updatetlevel
  IMPLICIT NONE
  INCLUDE 'fftw3-mpi.f03'

  COMPLEX(dp), DIMENSION(ikxs:ikxe,ikys:ikye) :: Fx_cmpx, Gy_cmpx
  COMPLEX(dp), DIMENSION(ikxs:ikxe,ikys:ikye) :: Fy_cmpx, Gx_cmpx, F_conv_G
  REAL(dp),    DIMENSION(ixs:ixe,iys:iye)     :: fr_real, gz_real
  REAL(dp),    DIMENSION(ixs:ixe,iys:iye)     :: fz_real, gr_real, f_times_g

  INTEGER :: in, is, p_int, j_int
  INTEGER :: nmax, smax ! Upper bound of the sums
  REAL(dp):: kx, ky, kerneln

  ! Execution time start
  CALL cpu_time(t0_Sapj)

  IF(NON_LIN .EQ. 1) THEN
    CALL compute_non_linear
  ELSEIF(NON_LIN .EQ. 0.5) THEN
    CALL compute_semi_linear
  ELSE
    Sepj = 0._dp; Sipj = 0._dp
  ENDIF

  ! Execution time END
  CALL cpu_time(t1_Sapj)
  tc_Sapj = tc_Sapj + (t1_Sapj - t0_Sapj)

CONTAINS

SUBROUTINE compute_non_linear
  IMPLICIT NONE
  !!!!!!!!!!!!!!!!!!!! ELECTRON non linear term computation (Sepj)!!!!!!!!!!
  IF(KIN_E) THEN
  zloope: DO iz = izs,ize
  ploope: DO ip = ips_e,ipe_e ! Loop over Hermite moments
    eo = MODULO(parray_e(ip),2)
    jloope: DO ij = ijs_e, ije_e ! Loop over Laguerre moments
      j_int=jarray_e(ij)
      real_data_c = 0._dp ! initialize sum over real nonlinear term
      ! Set non linear sum truncation
      IF (NL_CLOS .EQ. -2) THEN
        nmax = Jmaxe
      ELSEIF (NL_CLOS .EQ. -1) THEN
        nmax = Jmaxe-(ij-1)
      ELSE
        nmax = NL_CLOS
      ENDIF

      nloope: DO in = 1,nmax+1 ! Loop over laguerre for the sum

        kxloope: DO ikx = ikxs,ikxe ! Loop over kx
          kyloope: DO iky = ikys,ikye ! Loop over ky
            kx     = kxarray(ikx)
            ky     = kyarray(iky)
            kerneln = kernel_e(in, ikx, iky, iz, eo)

            ! First convolution terms
            Fx_cmpx(ikx,iky) = imagu*kx* phi(ikx,iky,iz) * kerneln
            Fy_cmpx(ikx,iky) = imagu*ky* phi(ikx,iky,iz) * kerneln
            ! Second convolution terms
            Gy_cmpx(ikx,iky) = 0._dp ! initialization of the sum
            Gx_cmpx(ikx,iky) = 0._dp ! initialization of the sum

            smax = MIN( (in-1)+(ij-1), jmaxe );
            DO is = 1, smax+1 ! sum truncation on number of moments
              Gy_cmpx(ikx,iky) = Gy_cmpx(ikx,iky) + &
                dnjs(in,ij,is) * moments_e(ip,is,ikx,iky,iz,updatetlevel)
              Gx_cmpx(ikx,iky) = Gx_cmpx(ikx,iky) + &
                dnjs(in,ij,is) * moments_e(ip,is,ikx,iky,iz,updatetlevel)
            ENDDO
            Gy_cmpx(ikx,iky) = imagu*ky*Gy_cmpx(ikx,iky)
            Gx_cmpx(ikx,iky) = imagu*kx*Gx_cmpx(ikx,iky)
          ENDDO kyloope
        ENDDO kxloope

        ! First term drphi x dzf
        DO ikx = ikxs, ikxe
          DO iky = ikys, ikye
            cmpx_data_f(iky,ikx-local_nkx_offset) = Fx_cmpx(ikx,iky)*AA_x(ikx)*AA_y(iky) !Anti aliasing filter
            cmpx_data_g(iky,ikx-local_nkx_offset) = Gy_cmpx(ikx,iky)*AA_x(ikx)*AA_y(iky) !Anti aliasing filter
          ENDDO
        ENDDO

        call fftw_mpi_execute_dft_c2r(planb, cmpx_data_f, real_data_f)
        call fftw_mpi_execute_dft_c2r(planb, cmpx_data_g, real_data_g)

        real_data_c = real_data_c + real_data_f/Ny/Nx  * real_data_g/Ny/Nx

        ! Second term -dzphi x drf
        DO ikx = ikxs, ikxe
          DO iky = ikys, ikye
            cmpx_data_f(iky,ikx-local_nkx_offset) = Fy_cmpx(ikx,iky)*AA_x(ikx)*AA_y(iky) !Anti aliasing filter
            cmpx_data_g(iky,ikx-local_nkx_offset) = Gx_cmpx(ikx,iky)*AA_x(ikx)*AA_y(iky) !Anti aliasing filter
          ENDDO
        ENDDO

        call fftw_mpi_execute_dft_c2r(planb, cmpx_data_f, real_data_f)
        call fftw_mpi_execute_dft_c2r(planb, cmpx_data_g, real_data_g)

        real_data_c = real_data_c - real_data_f/Ny/Nx  * real_data_g/Ny/Nx

      ENDDO nloope

      ! Put the real nonlinear product into k-space
      call fftw_mpi_execute_dft_r2c(planf, real_data_c, cmpx_data_c)

      ! Retrieve convolution in input format
      DO ikx = ikxs, ikxe
        DO iky = ikys, ikye
          Sepj(ip,ij,ikx,iky,iz) = cmpx_data_c(iky,ikx-local_nkx_offset)*AA_x(ikx)*AA_y(iky) !Anti aliasing filter
        ENDDO
      ENDDO
    ENDDO jloope
  ENDDO ploope
ENDDO zloope
ENDIF
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!! ION non linear term computation (Sipj)!!!!!!!!!!
zloopi: DO iz = izs,ize
  ploopi: DO ip = ips_i,ipe_i ! Loop over Hermite moments
    eo = MODULO(parray_i(ip),2)
    ! we check if poly degree is even (eq to index is odd) to spare computation
    !EVEN_P_i: IF (.TRUE. .OR. (MODULO(ip,2) .EQ. 1) .OR. (.NOT. COMPUTE_ONLY_EVEN_P)) THEN
    jloopi: DO ij = ijs_i, ije_i ! Loop over Laguerre moments
      j_int=jarray_i(ij)
      real_data_c = 0._dp ! initialize sum over real nonlinear term
      ! Set non linear sum truncation
      IF (NL_CLOS .EQ. -2) THEN
        nmax = Jmaxi
      ELSEIF (NL_CLOS .EQ. -1) THEN
        nmax = Jmaxi-(ij-1)
      ELSE
        nmax = NL_CLOS
      ENDIF

      nloopi: DO in = 1,nmax+1 ! Loop over laguerre for the sum

        kxloopi: DO ikx = ikxs,ikxe ! Loop over kx
          kyloopi: DO iky = ikys,ikye ! Loop over ky
            kx      = kxarray(ikx)
            ky      = kyarray(iky)
            kerneln = kernel_i(in, ikx, iky, iz, eo)

            ! First convolution terms
            Fx_cmpx(ikx,iky) = imagu*kx* phi(ikx,iky,iz) * kerneln
            Fy_cmpx(ikx,iky) = imagu*ky* phi(ikx,iky,iz) * kerneln
            ! Second convolution terms
            Gy_cmpx(ikx,iky) = 0._dp ! initialization of the sum
            Gx_cmpx(ikx,iky) = 0._dp ! initialization of the sum

            smax = MIN( (in-1)+(ij-1), jmaxi );
            DO is = 1, smax+1 ! sum truncation on number of moments
              Gy_cmpx(ikx,iky) = Gy_cmpx(ikx,iky) + &
                dnjs(in,ij,is) * moments_i(ip,is,ikx,iky,iz,updatetlevel)
              Gx_cmpx(ikx,iky) = Gx_cmpx(ikx,iky) + &
                dnjs(in,ij,is) * moments_i(ip,is,ikx,iky,iz,updatetlevel)
            ENDDO
            Gy_cmpx(ikx,iky) = imagu*ky*Gy_cmpx(ikx,iky)
            Gx_cmpx(ikx,iky) = imagu*kx*Gx_cmpx(ikx,iky)
          ENDDO kyloopi
        ENDDO kxloopi

        ! First term drphi x dzf
        DO ikx = ikxs, ikxe
          DO iky = ikys, ikye
            cmpx_data_f(iky,ikx-local_nkx_offset) = Fx_cmpx(ikx,iky)*AA_x(ikx)*AA_y(iky)
            cmpx_data_g(iky,ikx-local_nkx_offset) = Gy_cmpx(ikx,iky)*AA_x(ikx)*AA_y(iky)
          ENDDO
        ENDDO

        call fftw_mpi_execute_dft_c2r(planb, cmpx_data_f, real_data_f)
        call fftw_mpi_execute_dft_c2r(planb, cmpx_data_g, real_data_g)

        real_data_c = real_data_c + real_data_f/Ny/Nx  * real_data_g/Ny/Nx

        ! Second term -dzphi x drf
        DO ikx = ikxs, ikxe
          DO iky = ikys, ikye
            cmpx_data_f(iky,ikx-local_nkx_offset) = Fy_cmpx(ikx,iky)*AA_x(ikx)*AA_y(iky)
            cmpx_data_g(iky,ikx-local_nkx_offset) = Gx_cmpx(ikx,iky)*AA_x(ikx)*AA_y(iky)
          ENDDO
        ENDDO

        call fftw_mpi_execute_dft_c2r(planb, cmpx_data_f, real_data_f)
        call fftw_mpi_execute_dft_c2r(planb, cmpx_data_g, real_data_g)

        real_data_c = real_data_c - real_data_f/Ny/Nx  * real_data_g/Ny/Nx

      ENDDO nloopi

      ! Put the real nonlinear product into k-space
      call fftw_mpi_execute_dft_r2c(planf, real_data_c, cmpx_data_c)

      ! Retrieve convolution in input format
      DO ikx = ikxs, ikxe
        DO iky = ikys, ikye
          Sipj(ip,ij,ikx,iky,iz) = cmpx_data_c(iky,ikx-local_nkx_offset)*AA_x(ikx)*AA_y(iky)
        ENDDO
      ENDDO
    ENDDO jloopi
  ENDDO ploopi
ENDDO zloopi
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE compute_non_linear
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Semi linear computation : Only NZ-ZF convolutions are kept
SUBROUTINE compute_semi_linear
  IMPLICIT NONE
  ! Update the ZF modes if not frozen
  IF (WIPE_ZF .NE. -2) THEN
  moments_e_ZF(ips_e:ipe_e,ijs_e:ije_e,ikxs:ikxe,izs:ize) = moments_e(ips_e:ipe_e,ijs_e:ije_e,ikxs:ikxe,iky_0,izs:ize,updatetlevel)
  moments_i_ZF(ips_i:ipe_i,ijs_i:ije_i,ikxs:ikxe,izs:ize) = moments_i(ips_i:ipe_i,ijs_i:ije_i,ikxs:ikxe,iky_0,izs:ize,updatetlevel)
  phi_ZF(ikxs:ikxe,izs:ize) = phi(ikxs:ikxe,iky_0,izs:ize)
  ENDIF

  !!!!!!!!!!!!!!!!!!!! ELECTRON semi linear term computation (Sepj)!!!!!!!!!!
  IF(KIN_E) THEN
  zloope: DO iz = izs,ize
  ploope: DO ip = ips_e,ipe_e ! Loop over Hermite moments
    eo = MODULO(parray_e(ip),2)
    jloope: DO ij = ijs_e, ije_e ! Loop over Laguerre moments
      j_int=jarray_e(ij)
      real_data_c = 0._dp ! initialize sum over real nonlinear term
      ! Set non linear sum truncation
      IF (NL_CLOS .EQ. -2) THEN
        nmax = Jmaxe
      ELSEIF (NL_CLOS .EQ. -1) THEN
        nmax = Jmaxe-(ij-1)
      ELSE
        nmax = NL_CLOS
      ENDIF

      nloope: DO in = 1,nmax+1 ! Loop over laguerre for the sum
        ! Build the terms to convolve
        kxloope: DO ikx = ikxs,ikxe ! Loop over kx
          kyloope: DO iky = ikys,ikye ! Loop over ky
            kx      = kxarray(ikx)
            ky      = kyarray(iky)
            kerneln = kernel_e(in, ikx, iky, iz, eo)

            ! Zonal terms (=0 for all ky not 0)
            Fx_cmpx(ikx,iky) = 0._dp
            Gx_cmpx(ikx,iky) = 0._dp
            IF(iky .EQ. iky_0) THEN
              Fx_cmpx(ikx,iky) = imagu*kx* phi_ZF(ikx,iz) * kerneln
              smax = MIN( (in-1)+(ij-1), jmaxe );
              DO is = 1, smax+1 ! sum truncation on number of moments
                Gx_cmpx(ikx,iky) = Gx_cmpx(ikx,iky) + &
                  dnjs(in,ij,is) * moments_e_ZF(ip,is,ikx,iz)
              ENDDO
              Gx_cmpx(ikx,iky) = imagu*kx*Gx_cmpx(ikx,iky)
            ENDIF

            ! NZ terms
            Fy_cmpx(ikx,iky) = imagu*ky* phi(ikx,iky,iz) * kerneln
            Gy_cmpx(ikx,iky) = 0._dp ! initialization of the sum
            smax = MIN( (in-1)+(ij-1), jmaxe );
            DO is = 1, smax+1 ! sum truncation on number of moments
              Gy_cmpx(ikx,iky) = Gy_cmpx(ikx,iky) + &
                dnjs(in,ij,is) * moments_e(ip,is,ikx,iky,iz,updatetlevel)
            ENDDO
            Gy_cmpx(ikx,iky) = imagu*ky*Gy_cmpx(ikx,iky)

          ENDDO kyloope
        ENDDO kxloope

        ! First term ddx(phi_ZF) x ddy(f)
        DO ikx = ikxs, ikxe
          DO iky = ikys, ikye
            cmpx_data_f(iky,ikx-local_nkx_offset) = Fx_cmpx(ikx,iky)*AA_x(ikx)*AA_y(iky) !Anti aliasing filter
            cmpx_data_g(iky,ikx-local_nkx_offset) = Gy_cmpx(ikx,iky)*AA_x(ikx)*AA_y(iky) !Anti aliasing filter
          ENDDO
        ENDDO

        call fftw_mpi_execute_dft_c2r(planb, cmpx_data_f, real_data_f)
        call fftw_mpi_execute_dft_c2r(planb, cmpx_data_g, real_data_g)

        real_data_c = real_data_c + real_data_f/Ny/Nx  * real_data_g/Ny/Nx

        ! Second term -dyphi x dxf_ZF
        DO ikx = ikxs, ikxe
          DO iky = ikys, ikye
            cmpx_data_f(iky,ikx-local_nkx_offset) = Fy_cmpx(ikx,iky)*AA_x(ikx)*AA_y(iky) !Anti aliasing filter
            cmpx_data_g(iky,ikx-local_nkx_offset) = Gx_cmpx(ikx,iky)*AA_x(ikx)*AA_y(iky) !Anti aliasing filter
          ENDDO
        ENDDO

        call fftw_mpi_execute_dft_c2r(planb, cmpx_data_f, real_data_f)
        call fftw_mpi_execute_dft_c2r(planb, cmpx_data_g, real_data_g)

        real_data_c = real_data_c - real_data_f/Ny/Nx  * real_data_g/Ny/Nx

      ENDDO nloope

      ! Put the real nonlinear product into k-space
      call fftw_mpi_execute_dft_r2c(planf, real_data_c, cmpx_data_c)

      ! Retrieve convolution in input format
      DO ikx = ikxs, ikxe
        DO iky = ikys, ikye
          Sepj(ip,ij,ikx,iky,iz) = cmpx_data_c(iky,ikx-local_nkx_offset)*AA_x(ikx)*AA_y(iky) !Anti aliasing filter
        ENDDO
      ENDDO
    ENDDO jloope
  ENDDO ploope
ENDDO zloope
ENDIF
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!! ION non linear term computation (Sipj)!!!!!!!!!!
zloopi: DO iz = izs,ize
  ploopi: DO ip = ips_i,ipe_i ! Loop over Hermite moments
    eo = MODULO(parray_i(ip),2)
    ! we check if poly degree is even (eq to index is odd) to spare computation
    !EVEN_P_i: IF (.TRUE. .OR. (MODULO(ip,2) .EQ. 1) .OR. (.NOT. COMPUTE_ONLY_EVEN_P)) THEN
    jloopi: DO ij = ijs_i, ije_i ! Loop over Laguerre moments
      j_int=jarray_i(ij)
      real_data_c = 0._dp ! initialize sum over real nonlinear term
      ! Set non linear sum truncation
      IF (NL_CLOS .EQ. -2) THEN
        nmax = Jmaxi
      ELSEIF (NL_CLOS .EQ. -1) THEN
        nmax = Jmaxi-(ij-1)
      ELSE
        nmax = NL_CLOS
      ENDIF

      nloopi: DO in = 1,nmax+1 ! Loop over laguerre for the sum

        kxloopi: DO ikx = ikxs,ikxe ! Loop over kx
          kyloopi: DO iky = ikys,ikye ! Loop over ky
            ! Zonal terms (=0 for all ky not 0)
            Fx_cmpx(ikx,iky) = 0._dp
            Gx_cmpx(ikx,iky) = 0._dp
            IF(iky .EQ. iky_0) THEN
              Fx_cmpx(ikx,iky) = imagu*kx* phi_ZF(ikx,iz) * kerneln
              smax = MIN( (in-1)+(ij-1), jmaxe );
              DO is = 1, smax+1 ! sum truncation on number of moments
                Gx_cmpx(ikx,iky) = Gx_cmpx(ikx,iky) + &
                  dnjs(in,ij,is) * moments_i_ZF(ip,is,ikx,iz)
              ENDDO
              Gx_cmpx(ikx,iky) = imagu*kx*Gx_cmpx(ikx,iky)
            ENDIF

            ! NZ terms
            Fy_cmpx(ikx,iky) = imagu*ky* phi(ikx,iky,iz) * kerneln
            Gy_cmpx(ikx,iky) = 0._dp ! initialization of the sum
            smax = MIN( (in-1)+(ij-1), jmaxi );
            DO is = 1, smax+1 ! sum truncation on number of moments
              Gy_cmpx(ikx,iky) = Gy_cmpx(ikx,iky) + &
                dnjs(in,ij,is) * moments_i(ip,is,ikx,iky,iz,updatetlevel)
            ENDDO
            Gy_cmpx(ikx,iky) = imagu*ky*Gy_cmpx(ikx,iky)
          ENDDO kyloopi
        ENDDO kxloopi

        ! First term drphi x dzf
        DO ikx = ikxs, ikxe
          DO iky = ikys, ikye
            cmpx_data_f(iky,ikx-local_nkx_offset) = Fx_cmpx(ikx,iky)*AA_x(ikx)*AA_y(iky)
            cmpx_data_g(iky,ikx-local_nkx_offset) = Gy_cmpx(ikx,iky)*AA_x(ikx)*AA_y(iky)
          ENDDO
        ENDDO

        call fftw_mpi_execute_dft_c2r(planb, cmpx_data_f, real_data_f)
        call fftw_mpi_execute_dft_c2r(planb, cmpx_data_g, real_data_g)

        real_data_c = real_data_c + real_data_f/Ny/Nx  * real_data_g/Ny/Nx

        ! Second term -dzphi x drf
        DO ikx = ikxs, ikxe
          DO iky = ikys, ikye
            cmpx_data_f(iky,ikx-local_nkx_offset) = Fy_cmpx(ikx,iky)*AA_x(ikx)*AA_y(iky)
            cmpx_data_g(iky,ikx-local_nkx_offset) = Gx_cmpx(ikx,iky)*AA_x(ikx)*AA_y(iky)
          ENDDO
        ENDDO

        call fftw_mpi_execute_dft_c2r(planb, cmpx_data_f, real_data_f)
        call fftw_mpi_execute_dft_c2r(planb, cmpx_data_g, real_data_g)

        real_data_c = real_data_c - real_data_f/Ny/Nx  * real_data_g/Ny/Nx

      ENDDO nloopi

      ! Put the real nonlinear product into k-space
      call fftw_mpi_execute_dft_r2c(planf, real_data_c, cmpx_data_c)

      ! Retrieve convolution in input format
      DO ikx = ikxs, ikxe
        DO iky = ikys, ikye
          Sipj(ip,ij,ikx,iky,iz) = cmpx_data_c(iky,ikx-local_nkx_offset)*AA_x(ikx)*AA_y(iky)
        ENDDO
      ENDDO
    ENDDO jloopi
  ENDDO ploopi
ENDDO zloopi
END SUBROUTINE compute_semi_linear

END SUBROUTINE compute_Sapj
