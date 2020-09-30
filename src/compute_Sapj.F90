SUBROUTINE compute_Sapj(Sepj, Sipj)

  USE array, ONLY : dnjs
  USE basic
  USE fourier
  USE fields!, ONLY : phi, moments_e, moments_i
  USE grid
  USE model
  USE prec_const
  USE time_integration!, ONLY : updatetlevel
  IMPLICIT NONE

  COMPLEX(dp), DIMENSION(Nkr,Nkz) :: F_, G_, CONV
  COMPLEX(dp), DIMENSION(ips_e:ipe_e, ijs_e:ije_e, Nkr, Nkz) :: Sepj
  COMPLEX(dp), DIMENSION(ips_i:ipe_i, ijs_i:ije_i, Nkr, Nkz) :: Sipj
  INTEGER :: in, is
  REAL(dp):: kr, kz, kernel, be_2, bi_2, factn
  REAL(dp):: sigmae2_taue_o2, sigmai2_taui_o2

  !!!!!!!!!!!!!!!!!!!! ELECTRON non linear term computation (Sepj)!!!!!!!!!!
  sigmae2_taue_o2 = sigma_e**2 * tau_e/2._dp ! factor of the Kernel argument
  ploope: DO ip = ips_e,ipe_e ! Loop over Hermite moments
    jloope: DO ij = ijs_e, ije_e ! Loop over Laguerre moments
      Sepj(ip,ij,:,:)  = 0._dp
      factn = 1

      nloope: DO in = 1,jmaxe+1 ! Loop over laguerre for the sum

        krloope: DO ikr = 1,Nkr ! Loop over kr
          kzloope: DO ikz = 1,Nkz ! Loop over kz
            kr     = krarray(ikr)
            kz     = kzarray(ikz)
            be_2   = sigmae2_taue_o2 * (kr**2 + kz**2)
            kernel = be_2**(in-1)/factn * EXP(-be_2)

            ! First convolution term
            F_(ikr,ikz) = (kz - kr)  * phi(ikr,ikz) * kernel
            ! Second convolution term
            G_(ikr,ikz) = 0._dp ! initialization of the sum
            DO is = 1, MIN( in+ij-1, jmaxe+1 ) ! sum truncation on number of moments
              G_(ikr,ikz) = G_(ikr,ikz) + &
               dnjs(in,ij,is) * moments_e(ip,is,ikr,ikz,updatetlevel)
            ENDDO
            G_(ikr,ikz) = (kz - kr) * G_(ikr,ikz)
          ENDDO kzloope
        ENDDO krloope

        CALL convolve_2D_F2F( F_, G_, CONV ) ! Convolve and go back to Fourier space
        !CALL convolve_2D_F2R( F_, G_, CONV ) ! .. or convolve and keep the results in real space**

        Sepj(ip,ij,:,:) = Sepj(ip,ij,:,:) + CONV ! Add it to Sepj (Real space here)

        IF ( in + 1 .LE. jmaxe+1 ) THEN
          factn = real(in,dp) * factn ! factorial(n+1)
        ENDIF
      ENDDO nloope

      !CALL forward_FFT(Sepj(ip,ij,:,:)) !**then put the sum back to fourier space (spares n FFT)

    ENDDO jloope
  ENDDO ploope
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !!!!!!!!!!!!!!!!!!!! ION non linear term computation (Sipj)!!!!!!!!!!
  sigmai2_taui_o2 = sigma_i**2 * tau_i/2._dp

  ploopi: DO ip = ips_i,ipe_i ! Loop over Hermite moments
    jloopi: DO ij = ijs_i, ije_i ! Loop over Laguerre moments
      Sipj(ip,ij,:,:)  = 0._dp
      factn = 1

      nloopi: DO in = 1,jmaxi+1 ! Loop over laguerre for the sum

        krloopi: DO ikr = 1,Nkr ! Loop over kr
          kzloopi: DO ikz = 1,Nkz ! Loop over kz
            kr     = krarray(ikr)
            kz     = kzarray(ikz)
            bi_2   = sigmai2_taui_o2 * (kr**2 + kz**2)
            kernel = bi_2**(in-1)/factn * EXP(-bi_2)

            F_(ikr,ikz) = (kz - kr)  * phi(ikr,ikz) * kernel

            G_(ikr,ikz) = 0._dp ! initialization of the sum
            DO is = 1, MIN( in+ij-1, jmaxi+1 )
              G_(ikr,ikz) = G_(ikr,ikz) + &
               dnjs(in,ij,is) * moments_i(ip,is,ikr,ikz,updatetlevel)
            ENDDO
            G_(ikr,ikz) = (kz - kr) * G_(ikr,ikz)
          ENDDO kzloopi
        ENDDO krloopi

        CALL convolve_2D_F2F( F_, G_, CONV ) ! Convolve and back to Fourier
        !CALL convolve_2D_F2R( F_, G_, CONV ) ! or convolve and keep the results in real space**

        Sipj(ip,ij,:,:) = Sipj(ip,ij,:,:) + CONV ! Add it to Sipj (Real space here)

        IF ( in + 1 .LE. jmaxi+1 ) THEN
          factn = real(in,dp) * factn ! factorial(n+1)
        ENDIF

      ENDDO nloopi

      !CALL forward_FFT(Sipj(ip,ij,:,:)) ! **then put the sum back to fourier space (spares n FFT)

    ENDDO jloopi
  ENDDO ploopi
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE compute_Sapj
