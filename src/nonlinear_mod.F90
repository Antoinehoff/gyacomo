MODULE nonlinear
  USE array, ONLY : dnjs, Sepj, Sipj, kernel_i, kernel_e,&
                    moments_e_ZF, moments_i_ZF, phi_ZF
  USE initial_par, ONLY : ACT_ON_MODES
  USE basic
  USE fourier
  USE fields, ONLY : phi, psi, moments_e, moments_i
  USE grid
  USE model
  USE prec_const
  USE time_integration, ONLY : updatetlevel
  IMPLICIT NONE
  INCLUDE 'fftw3-mpi.f03'

  COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE :: F_cmpx, G_cmpx
  COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE :: Fx_cmpx, Gy_cmpx
  COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE :: Fy_cmpx, Gx_cmpx, F_conv_G

  INTEGER :: in, is, p_int, j_int
  INTEGER :: nmax, smax ! Upper bound of the sums
  REAL(dp):: kx, ky, kerneln, sqrt_p, sqrt_pp1
  PUBLIC :: compute_Sapj, nonlinear_init

CONTAINS

SUBROUTINE nonlinear_init
  ALLOCATE( F_cmpx(ikys:ikye,ikxs:ikxe))
  ALLOCATE( G_cmpx(ikys:ikye,ikxs:ikxe))

  ALLOCATE(Fx_cmpx(ikys:ikye,ikxs:ikxe))
  ALLOCATE(Gy_cmpx(ikys:ikye,ikxs:ikxe))
  ALLOCATE(Fy_cmpx(ikys:ikye,ikxs:ikxe))
  ALLOCATE(Gx_cmpx(ikys:ikye,ikxs:ikxe))

  ALLOCATE(F_conv_G(ikys:ikye,ikxs:ikxe))
END SUBROUTINE nonlinear_init

SUBROUTINE compute_Sapj
  ! This routine is meant to compute the non linear term for each specie and degree
  !! In real space Sapj ~ b*(grad(phi) x grad(g)) which in moments in fourier becomes
  !! Sapj = Sum_n (ikx Kn phi)#(iky Sum_s d_njs Naps) - (iky Kn phi)#(ikx Sum_s d_njs Naps)
  !! where # denotes the convolution.

  ! Execution time start
  CALL cpu_time(t0_Sapj)

  SELECT CASE(LINEARITY)
    CASE ('nonlinear')
      CALL compute_nonlinear
    CASE ('ZF_semilin')
      CALL compute_semi_linear_ZF
    CASE ('NZ_semilin')
      CALL compute_semi_linear_NZ
    CASE ('linear')
      Sepj = 0._dp; Sipj = 0._dp
    CASE DEFAULT
      IF(my_id.EQ.0) write(*,*) '/!\ Linearity not recognized /!\'
      stop
  END SELECT

  ! Execution time END
  CALL cpu_time(t1_Sapj)
  tc_Sapj = tc_Sapj + (t1_Sapj - t0_Sapj)

END SUBROUTINE compute_Sapj

SUBROUTINE compute_nonlinear
  IMPLICIT NONE
  !!!!!!!!!!!!!!!!!!!! ELECTRON non linear term computation (Sepj)!!!!!!!!!!
  IF(KIN_E) THEN
  zloope: DO iz = izs,ize

    ploope: DO ip = ips_e,ipe_e ! Loop over Hermite moments
      eo     = MODULO(parray_e(ip),2)
      p_int  = parray_e(ip)
      sqrt_p = SQRT(REAL(p_int,dp)); sqrt_pp1 = SQRT(REAL(p_int,dp)+1._dp);

      jloope: DO ij = ijs_e, ije_e ! Loop over Laguerre moments
        j_int=jarray_e(ij)
        IF((CLOS .NE. 1) .OR. (p_int+2*j_int .LE. dmaxe)) THEN !compute
          ! Set non linear sum truncation
          IF (NL_CLOS .EQ. -2) THEN
            nmax = Jmaxe
          ELSEIF (NL_CLOS .EQ. -1) THEN
            nmax = Jmaxe-j_int
          ELSE
            nmax = min(NL_CLOS,Jmaxe-j_int)
          ENDIF
          bracket_sum_r = 0._dp ! initialize sum over real nonlinear term

          nloope: DO in = 1,nmax+1 ! Loop over laguerre for the sum
!-----------!! ELECTROSTATIC CONTRIBUTION {Sum_s dnjs Naps, Kernel phi}
            ! First convolution terms
            F_cmpx(ikys:ikye,ikxs:ikxe) = phi(ikys:ikye,ikxs:ikxe,iz) * kernel_e(in, ikys:ikye,ikxs:ikxe, iz, eo)
            ! Second convolution terms
            G_cmpx(ikys:ikye,ikxs:ikxe) = 0._dp ! initialization of the sum
            smax = MIN( (in-1)+(ij-1), Jmaxe );
            DO is = 1, smax+1 ! sum truncation on number of moments
              G_cmpx(ikys:ikye,ikxs:ikxe)  = G_cmpx(ikys:ikye,ikxs:ikxe) + &
                dnjs(in,ij,is) * moments_e(ip,is,ikys:ikye,ikxs:ikxe,iz,updatetlevel)
            ENDDO
            !/!\ this function add its result to bracket_sum_r (hard to read sorry) /!\
            CALL poisson_bracket_and_sum(F_cmpx,G_cmpx)

!-----------!! ELECTROMAGNETIC CONTRIBUTION -sqrt(tau)/sigma*{Sum_s dnjs [sqrt(p+1)Nap+1s + sqrt(p)Nap-1s], Kernel psi}
            IF(EM) THEN
            ! First convolution terms
            F_cmpx(ikys:ikye,ikxs:ikxe) = -sqrt_tau_o_sigma_e * psi(ikys:ikye,ikxs:ikxe,iz) * kernel_e(in, ikys:ikye,ikxs:ikxe, iz, eo)
            ! Second convolution terms
            G_cmpx(ikys:ikye,ikxs:ikxe) = 0._dp ! initialization of the sum
            smax = MIN( (in-1)+(ij-1), Jmaxe );
            DO is = 1, smax+1 ! sum truncation on number of moments
              G_cmpx(ikys:ikye,ikxs:ikxe)  = G_cmpx(ikys:ikye,ikxs:ikxe) + &
                dnjs(in,ij,is) * (sqrt_pp1*moments_e(ip+1,is,ikys:ikye,ikxs:ikxe,iz,updatetlevel)&
                                 +sqrt_p  *moments_e(ip-1,is,ikys:ikye,ikxs:ikxe,iz,updatetlevel))
            ENDDO
            !/!\ this function add its result to bracket_sum_r (hard to read sorry) /!\
            CALL poisson_bracket_and_sum(F_cmpx,G_cmpx)
            ENDIF
          ENDDO nloope

!---------! Put back the real nonlinear product into k-space
          call fftw_mpi_execute_dft_r2c(planf, bracket_sum_r, bracket_sum_c)
          ! Retrieve convolution in input format
          DO iky = ikys, ikye
            Sepj(ip,ij,iky,ikxs:ikxe,iz) = bracket_sum_c(ikxs:ikxe,iky-local_nky_offset)*AA_x(ikxs:ikxe)*AA_y(iky) !Anti aliasing filter
          ENDDO
        ELSE
          Sepj(ip,ij,:,:,iz) = 0._dp
        ENDIF
      ENDDO jloope
    ENDDO ploope
  ENDDO zloope
ENDIF
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!! ION non linear term computation (Sipj)!!!!!!!!!!
  zloopi: DO iz = izs,ize
    ploopi: DO ip = ips_i,ipe_i ! Loop over Hermite moments
      p_int = parray_i(ip)
      eo = MODULO(parray_i(ip),2)
      jloopi: DO ij = ijs_i, ije_i ! Loop over Laguerre moments
        j_int=jarray_i(ij)
        IF((CLOS .NE. 1) .OR. (p_int+2*j_int .LE. dmaxi)) THEN !compute for every moments except for closure 1
          ! Set non linear sum truncation
          IF (NL_CLOS .EQ. -2) THEN
            nmax = Jmaxi
          ELSEIF (NL_CLOS .EQ. -1) THEN
            nmax = Jmaxi-j_int
          ELSE
            nmax = min(NL_CLOS,Jmaxi-j_int)
          ENDIF
          bracket_sum_r = 0._dp ! initialize sum over real nonlinear term
          nloopi: DO in = 1,nmax+1 ! Loop over laguerre for the sum
!-----------!! ELECTROSTATIC CONTRIBUTION
            ! First convolution terms
            F_cmpx(ikys:ikye,ikxs:ikxe) = phi(ikys:ikye,ikxs:ikxe,iz) * kernel_i(in, ikys:ikye,ikxs:ikxe, iz, eo)
            ! Second convolution terms
            G_cmpx(ikys:ikye,ikxs:ikxe) = 0._dp ! initialization of the sum
            smax = MIN( (in-1)+(ij-1), jmaxi );
            DO is = 1, smax+1 ! sum truncation on number of moments
              G_cmpx(ikys:ikye,ikxs:ikxe) = G_cmpx(ikys:ikye,ikxs:ikxe) + &
                dnjs(in,ij,is) * moments_i(ip,is,ikys:ikye,ikxs:ikxe,iz,updatetlevel)
            ENDDO
            !/!\ this function add its result to bracket_sum_r (hard to read sorry) /!\
            CALL poisson_bracket_and_sum(F_cmpx,G_cmpx)
!-----------!! ELECTROMAGNETIC CONTRIBUTION -sqrt(tau)/sigma*{Sum_s dnjs [sqrt(p+1)Nap+1s + sqrt(p)Nap-1s], Kernel psi}
            IF(EM) THEN
            ! First convolution terms
            F_cmpx(ikys:ikye,ikxs:ikxe) = -sqrt_tau_o_sigma_i * psi(ikys:ikye,ikxs:ikxe,iz) * kernel_i(in, ikys:ikye,ikxs:ikxe, iz, eo)
            ! Second convolution terms
            G_cmpx(ikys:ikye,ikxs:ikxe) = 0._dp ! initialization of the sum
            smax = MIN( (in-1)+(ij-1), Jmaxi );
            DO is = 1, smax+1 ! sum truncation on number of moments
              G_cmpx(ikys:ikye,ikxs:ikxe)  = G_cmpx(ikys:ikye,ikxs:ikxe) + &
                dnjs(in,ij,is) * (sqrt_pp1*moments_i(ip+1,is,ikys:ikye,ikxs:ikxe,iz,updatetlevel)&
                                 +sqrt_p  *moments_i(ip-1,is,ikys:ikye,ikxs:ikxe,iz,updatetlevel))
            ENDDO
            !/!\ this function add its result to bracket_sum_r (hard to read sorry) /!\
            CALL poisson_bracket_and_sum(F_cmpx,G_cmpx)
            ENDIF
          ENDDO nloopi
          ! Put the real nonlinear product into k-space
          call fftw_mpi_execute_dft_r2c(planf, bracket_sum_r, bracket_sum_c)
          ! Retrieve convolution in input format
          DO iky = ikys, ikye
            Sipj(ip,ij,iky,ikxs:ikxe,iz) = bracket_sum_c(ikxs:ikxe,iky-local_nky_offset)*AA_x(ikxs:ikxe)*AA_y(iky)
          ENDDO
        ELSE
          Sipj(ip,ij,:,:,iz) = 0._dp
        ENDIF
      ENDDO jloopi
    ENDDO ploopi
  ENDDO zloopi
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE compute_nonlinear
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Semi linear computation : Only NZ-ZF convolutions are kept
SUBROUTINE compute_semi_linear_ZF
  IMPLICIT NONE
  !!!!!!!!!!!!!!!!!!!! ELECTRON semi linear term computation (Sepj)!!!!!!!!!!
  IF(KIN_E) THEN
  zloope: DO iz = izs,ize
  ploope: DO ip = ips_e,ipe_e ! Loop over Hermite moments
    eo = MODULO(parray_e(ip),2)
    jloope: DO ij = ijs_e, ije_e ! Loop over Laguerre moments
      j_int=jarray_e(ij)
      ! Set non linear sum truncation
      IF (NL_CLOS .EQ. -2) THEN
        nmax = Jmaxe
      ELSEIF (NL_CLOS .EQ. -1) THEN
        nmax = Jmaxe-(ij-1)
      ELSE
        nmax = NL_CLOS
      ENDIF
      bracket_sum_r = 0._dp ! initialize sum over real nonlinear term
      nloope: DO in = 1,nmax+1 ! Loop over laguerre for the sum
        ! Build the terms to convolve
        kxloope: DO ikx = ikxs,ikxe ! Loop over kx
          kyloope: DO iky = ikys,ikye ! Loop over ky
            kx      = kxarray(ikx)
            ky      = kyarray(iky)
            kerneln = kernel_e(in, ikx, iky, iz, eo)
            ! Zonal terms (=0 for all ky not 0)
            Fx_cmpx(iky,ikx) = 0._dp
            Gx_cmpx(iky,ikx) = 0._dp
            IF(iky .EQ. iky_0) THEN
              Fx_cmpx(iky,ikx) = imagu*kx* phi(iky,ikx,iz) * kerneln
              smax = MIN( (in-1)+(ij-1), jmaxe );
              DO is = 1, smax+1 ! sum truncation on number of moments
                Gx_cmpx(iky,ikx) = Gx_cmpx(iky,ikx) + &
                  dnjs(in,ij,is) * moments_e(ip,is,iky,ikx,iz,updatetlevel)
              ENDDO
              Gx_cmpx(iky,ikx) = imagu*kx*Gx_cmpx(iky,ikx)
            ENDIF
            ! NZ terms
            Fy_cmpx(iky,ikx) = imagu*ky* phi(iky,ikx,iz) * kerneln
            Gy_cmpx(iky,ikx) = 0._dp ! initialization of the sum
            smax = MIN( (in-1)+(ij-1), jmaxe );
            DO is = 1, smax+1 ! sum truncation on number of moments
              Gy_cmpx(iky,ikx) = Gy_cmpx(iky,ikx) + &
                dnjs(in,ij,is) * moments_e(ip,is,iky,ikx,iz,updatetlevel)
            ENDDO
            Gy_cmpx(iky,ikx) = imagu*ky*Gy_cmpx(iky,ikx)
          ENDDO kyloope
        ENDDO kxloope
        ! First term df/dx x dg/dy
        CALL convolve_and_add(Fx_cmpx,Gy_cmpx)
        ! Second term -df/dy x dg/dx
        CALL convolve_and_add(-Fy_cmpx,Gx_cmpx)
      ENDDO nloope
      ! Put the real nonlinear product into k-space
      call fftw_mpi_execute_dft_r2c(planf, bracket_sum_r, bracket_sum_c)
      ! Retrieve convolution in input format
      DO ikx = ikxs, ikxe
        DO iky = ikys, ikye
          Sepj(ip,ij,iky,ikx,iz) = bracket_sum_c(ikx,iky-local_nky_offset)*AA_x(ikx)*AA_y(iky) !Anti aliasing filter
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
    jloopi: DO ij = ijs_i, ije_i ! Loop over Laguerre moments
      j_int=jarray_i(ij)
      ! Set non linear sum truncation
      IF (NL_CLOS .EQ. -2) THEN
        nmax = Jmaxi
      ELSEIF (NL_CLOS .EQ. -1) THEN
        nmax = Jmaxi-(ij-1)
      ELSE
        nmax = NL_CLOS
      ENDIF
      bracket_sum_r = 0._dp ! initialize sum over real nonlinear term
      nloopi: DO in = 1,nmax+1 ! Loop over laguerre for the sum
        kxloopi: DO ikx = ikxs,ikxe ! Loop over kx
          kyloopi: DO iky = ikys,ikye ! Loop over ky
            ! Zonal terms (=0 for all ky not 0)
            Fx_cmpx(iky,ikx) = 0._dp
            Gx_cmpx(iky,ikx) = 0._dp
            IF(iky .EQ. iky_0) THEN
              Fx_cmpx(iky,ikx) = imagu*kx* phi(iky,ikx,iz) * kerneln
              smax = MIN( (in-1)+(ij-1), jmaxi );
              DO is = 1, smax+1 ! sum truncation on number of moments
                Gx_cmpx(iky,ikx) = Gx_cmpx(iky,ikx) + &
                  dnjs(in,ij,is) * moments_i(ip,is,iky,ikx,iz,updatetlevel)
              ENDDO
              Gx_cmpx(iky,ikx) = imagu*kx*Gx_cmpx(iky,ikx)
            ENDIF

            ! NZ terms
            Fy_cmpx(iky,ikx) = imagu*ky* phi(iky,ikx,iz) * kerneln
            Gy_cmpx(iky,ikx) = 0._dp ! initialization of the sum
            smax = MIN( (in-1)+(ij-1), jmaxi );
            DO is = 1, smax+1 ! sum truncation on number of moments
              Gy_cmpx(iky,ikx) = Gy_cmpx(iky,ikx) + &
                dnjs(in,ij,is) * moments_i(ip,is,iky,ikx,iz,updatetlevel)
            ENDDO
            Gy_cmpx(iky,ikx) = imagu*ky*Gy_cmpx(iky,ikx)
          ENDDO kyloopi
        ENDDO kxloopi
        ! First term drphi x dzf
        CALL convolve_and_add(Fy_cmpx,Gx_cmpx)
        ! Second term -dzphi x drf
        CALL convolve_and_add(Fy_cmpx,Gx_cmpx)
      ENDDO nloopi
      ! Put the real nonlinear product into k-space
      call fftw_mpi_execute_dft_r2c(planf, bracket_sum_r, bracket_sum_c)
      ! Retrieve convolution in input format
      DO ikx = ikxs, ikxe
        DO iky = ikys, ikye
          Sipj(ip,ij,iky,ikx,iz) = bracket_sum_c(ikx,iky-local_nky_offset)*AA_x(ikx)*AA_y(iky)
        ENDDO
      ENDDO
    ENDDO jloopi
  ENDDO ploopi
ENDDO zloopi
END SUBROUTINE compute_semi_linear_ZF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Semi linear computation : Only kx=0*all convolutions are kept
SUBROUTINE compute_semi_linear_NZ
  IMPLICIT NONE
  !!!!!!!!!!!!!!!!!!!! ELECTRON semi linear term computation (Sepj)!!!!!!!!!!
  IF(KIN_E) THEN
  zloope: DO iz = izs,ize
  ploope: DO ip = ips_e,ipe_e ! Loop over Hermite moments
    eo = MODULO(parray_e(ip),2)
    jloope: DO ij = ijs_e, ije_e ! Loop over Laguerre moments
      j_int=jarray_e(ij)
      ! Set non linear sum truncation
      IF (NL_CLOS .EQ. -2) THEN
        nmax = Jmaxe
      ELSEIF (NL_CLOS .EQ. -1) THEN
        nmax = Jmaxe-(ij-1)
      ELSE
        nmax = NL_CLOS
      ENDIF
      bracket_sum_r = 0._dp ! initialize sum over real nonlinear term
      nloope: DO in = 1,nmax+1 ! Loop over laguerre for the sum
        ! Build the terms to convolve
        kxloope: DO ikx = ikxs,ikxe ! Loop over kx
          kyloope: DO iky = ikys,ikye ! Loop over ky
            kx      = kxarray(ikx)
            ky      = kyarray(iky)
            kerneln = kernel_e(in, ikx, iky, iz, eo)
            ! All terms
            Fx_cmpx(iky,ikx) = imagu*kx* phi(iky,ikx,iz) * kerneln
            smax = MIN( (in-1)+(ij-1), jmaxe );
            DO is = 1, smax+1 ! sum truncation on number of moments
              Gx_cmpx(iky,ikx) = Gx_cmpx(iky,ikx) + &
                dnjs(in,ij,is) * moments_e(ip,is,iky,ikx,iz,updatetlevel)
            ENDDO
            Gx_cmpx(iky,ikx) = imagu*kx*Gx_cmpx(iky,ikx)
            ! Kx = 0 terms
            Fy_cmpx(iky,ikx) = 0._dp
            Gy_cmpx(iky,ikx) = 0._dp
            IF (ikx .EQ. ikx_0) THEN
              Fy_cmpx(iky,ikx) = imagu*ky* phi(iky,ikx,iz) * kerneln
              Gy_cmpx(iky,ikx) = 0._dp ! initialization of the sum
              smax = MIN( (in-1)+(ij-1), jmaxe );
              DO is = 1, smax+1 ! sum truncation on number of moments
                Gy_cmpx(iky,ikx) = Gy_cmpx(iky,ikx) + &
                  dnjs(in,ij,is) * moments_e(ip,is,iky,ikx,iz,updatetlevel)
              ENDDO
              Gy_cmpx(iky,ikx) = imagu*ky*Gy_cmpx(iky,ikx)
            ENDIF
          ENDDO kyloope
        ENDDO kxloope
        ! First term df/dx x dg/dy
        CALL convolve_and_add(Fx_cmpx,Gy_cmpx)
        ! Second term -df/dy x dg/dx
        CALL convolve_and_add(-Fy_cmpx,Gx_cmpx)
      ENDDO nloope
      ! Put the real nonlinear product into k-space
      call fftw_mpi_execute_dft_r2c(planf, bracket_sum_r, bracket_sum_c)
      ! Retrieve convolution in input format
      DO ikx = ikxs, ikxe
        DO iky = ikys, ikye
          Sepj(ip,ij,iky,ikx,iz) = bracket_sum_c(ikx,iky-local_nky_offset)*AA_x(ikx)*AA_y(iky) !Anti aliasing filter
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
    jloopi: DO ij = ijs_i, ije_i ! Loop over Laguerre moments
      j_int=jarray_i(ij)
      ! Set non linear sum truncation
      IF (NL_CLOS .EQ. -2) THEN
        nmax = Jmaxi
      ELSEIF (NL_CLOS .EQ. -1) THEN
        nmax = Jmaxi-(ij-1)
      ELSE
        nmax = NL_CLOS
      ENDIF
      bracket_sum_r = 0._dp ! initialize sum over real nonlinear term
      nloopi: DO in = 1,nmax+1 ! Loop over laguerre for the sum
        kxloopi: DO ikx = ikxs,ikxe ! Loop over kx
          kyloopi: DO iky = ikys,ikye ! Loop over ky
            ! Zonal terms (=0 for all ky not 0)
            Fx_cmpx(iky,ikx) = imagu*kx* phi(iky,ikx,iz) * kerneln
            smax = MIN( (in-1)+(ij-1), jmaxi );
            DO is = 1, smax+1 ! sum truncation on number of moments
              Gx_cmpx(iky,ikx) = Gx_cmpx(iky,ikx) + &
                dnjs(in,ij,is) * moments_i(ip,is,iky,ikx,iz,updatetlevel)
            ENDDO
            Gx_cmpx(iky,ikx) = imagu*kx*Gx_cmpx(iky,ikx)

            ! Kx = 0 terms
            Fy_cmpx(iky,ikx) = 0._dp
            Gy_cmpx(iky,ikx) = 0._dp
            IF (ikx .EQ. ikx_0) THEN
              Fy_cmpx(iky,ikx) = imagu*ky* phi(iky,ikx,iz) * kerneln
              Gy_cmpx(iky,ikx) = 0._dp ! initialization of the sum
              smax = MIN( (in-1)+(ij-1), jmaxi );
              DO is = 1, smax+1 ! sum truncation on number of moments
                Gy_cmpx(iky,ikx) = Gy_cmpx(iky,ikx) + &
                  dnjs(in,ij,is) * moments_i(ip,is,iky,ikx,iz,updatetlevel)
              ENDDO
              Gy_cmpx(iky,ikx) = imagu*ky*Gy_cmpx(iky,ikx)
            ENDIF
          ENDDO kyloopi
        ENDDO kxloopi
        ! First term drphi x dzf
        CALL convolve_and_add(Fy_cmpx,Gx_cmpx)
        ! Second term -dzphi x drf
        CALL convolve_and_add(Fy_cmpx,Gx_cmpx)
      ENDDO nloopi
      ! Put the real nonlinear product into k-space
      call fftw_mpi_execute_dft_r2c(planf, bracket_sum_r, bracket_sum_c)
      ! Retrieve convolution in input format
      DO ikx = ikxs, ikxe
        DO iky = ikys, ikye
          Sipj(ip,ij,iky,ikx,iz) = bracket_sum_c(ikx,iky-local_nky_offset)*AA_x(ikx)*AA_y(iky)
        ENDDO
      ENDDO
    ENDDO jloopi
  ENDDO ploopi
ENDDO zloopi
END SUBROUTINE compute_semi_linear_NZ

END MODULE nonlinear
