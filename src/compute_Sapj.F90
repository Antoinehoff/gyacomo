SUBROUTINE compute_Sapj

  USE array, ONLY : dnjs, Sepj, Sipj
  USE basic
  USE fourier
  USE fields!, ONLY : phi, moments_e, moments_i
  USE grid
  USE model
  USE prec_const
  USE time_integration!, ONLY : updatetlevel
  IMPLICIT NONE

  COMPLEX(dp), DIMENSION(Nkr,Nkz) :: F_, G_, CONV
  ! COMPLEX(dp), DIMENSION(ips_e:ipe_e, ijs_e:ije_e, Nkr, Nkz) :: Sepj
  ! COMPLEX(dp), DIMENSION(ips_i:ipe_i, ijs_i:ije_i, Nkr, Nkz) :: Sipj
  INTEGER :: in, is
  REAL(dp):: kr, kz, kerneln, be_2, bi_2, factn
  REAL(dp):: sigmae2_taue_o2, sigmai2_taui_o2

  !!!!!!!!!!!!!!!!!!!! ELECTRON non linear term computation (Sepj)!!!!!!!!!!
  sigmae2_taue_o2 = sigma_e**2 * tau_e/2._dp ! factor of the kerneln argument
  ploope: DO ip = ips_e,ipe_e ! Loop over Hermite moments
    jloope: DO ij = ijs_e, ije_e ! Loop over Laguerre moments
      Sepj(ip,ij,:,:)  = 0._dp
      factn = 1

      nloope: DO in = 1,jmaxe+1 ! Loop over laguerre for the sum

        krloope1: DO ikr = 1,Nkr ! Loop over kr
          kzloope1: DO ikz = 1,Nkz ! Loop over kz
            kr     = krarray(ikr)
            kz     = kzarray(ikz)
            be_2   = sigmae2_taue_o2 * (kr**2 + kz**2)
            kerneln = be_2**(in-1)/factn * EXP(-be_2)

            ! First convolution term
            F_(ikr,ikz) = imagu*kz* phi(ikr,ikz) * kerneln
            ! Second convolution term
            G_(ikr,ikz) = 0._dp ! initialization of the sum
            DO is = 1, MIN( in+ij-1, jmaxe+1 ) ! sum truncation on number of moments
              G_(ikr,ikz) = G_(ikr,ikz) + &
               dnjs(in,ij,is) * moments_e(ip,is,ikr,ikz,updatetlevel)
            ENDDO
            G_(ikr,ikz) = imagu*kr*G_(ikr,ikz)

          ENDDO kzloope1
        ENDDO krloope1

        CALL convolve_2D_F2F( F_, G_, CONV )
        Sepj(ip,ij,:,:) = Sepj(ip,ij,:,:) + CONV ! Add it to Sepj (Real space here)

        krloope2: DO ikr = 1,Nkr ! Loop over kr
          kzloope2: DO ikz = 1,Nkz ! Loop over kz
            kr     = krarray(ikr)
            kz     = kzarray(ikz)
            be_2   = sigmae2_taue_o2 * (kr**2 + kz**2)
            kerneln = be_2**(in-1)/factn * EXP(-be_2)

            ! First convolution term
            F_(ikr,ikz) = imagu*kr* phi(ikr,ikz) * kerneln
            ! Second convolution term
            G_(ikr,ikz) = 0._dp ! initialization of the sum
            DO is = 1, MIN( in+ij-1, jmaxe+1 ) ! sum truncation on number of moments
              G_(ikr,ikz) = G_(ikr,ikz) + &
               dnjs(in,ij,is) * moments_e(ip,is,ikr,ikz,updatetlevel)
            ENDDO
            G_(ikr,ikz) = imagu*kz*G_(ikr,ikz)

          ENDDO kzloope2
        ENDDO krloope2


        CALL convolve_2D_F2F( F_, G_, CONV )
        Sepj(ip,ij,:,:) = Sepj(ip,ij,:,:) - CONV ! Add it to Sepj (Real space here)

        IF ( in + 1 .LE. jmaxe+1 ) THEN
          factn = real(in,dp) * factn ! factorial(n+1)
        ENDIF
      ENDDO nloope

      !! Antialiasing
      DO ikr = ikrs,ikre
        DO ikz = ikzs, ikze
          Sepj(ip,ij,ikr,ikz) = Sepj(ip,ij,ikr,ikz) * AA_r(ikr) * AA_z(ikz)
        ENDDO
      ENDDO

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

        krloopi1: DO ikr = 1,Nkr ! Loop over kr
          kzloopi1: DO ikz = 1,Nkz ! Loop over kz
            kr      = krarray(ikr)
            kz      = kzarray(ikz)
            bi_2    = sigmai2_taui_o2 * (kr**2 + kz**2)
            kerneln = bi_2**(in-1)/factn * EXP(-bi_2)

            ! First convolution term
            F_(ikr,ikz) = imagu*kz*phi(ikr,ikz) * kerneln
            ! Second convolution term
            G_(ikr,ikz) = 0._dp ! initialization of the sum
            DO is = 1, MIN( in+ij-1, jmaxi+1 )
              G_(ikr,ikz) = G_(ikr,ikz) + &
               dnjs(in,ij,is) * moments_i(ip,is,ikr,ikz,updatetlevel)
            ENDDO
            G_(ikr,ikz) = imagu*kr*G_(ikr,ikz)

          ENDDO kzloopi1
        ENDDO krloopi1

        CALL convolve_2D_F2F( F_, G_, CONV ) ! Convolve and back to Fourier
        Sipj(ip,ij,:,:) = Sipj(ip,ij,:,:) + CONV ! Add it to Sipj (Real space here)

        krloopi2: DO ikr = 1,Nkr ! Loop over kr
          kzloopi2: DO ikz = 1,Nkz ! Loop over kz
            kr      = krarray(ikr)
            kz      = kzarray(ikz)
            bi_2    = sigmai2_taui_o2 * (kr**2 + kz**2)
            kerneln = bi_2**(in-1)/factn * EXP(-bi_2)

            ! First convolution term
            F_(ikr,ikz) = imagu*kr*phi(ikr,ikz) * kerneln
            ! Second convolution term
            G_(ikr,ikz) = 0._dp ! initialization of the sum
            DO is = 1, MIN( in+ij-1, jmaxi+1 )
              G_(ikr,ikz) = G_(ikr,ikz) + &
               dnjs(in,ij,is) * moments_i(ip,is,ikr,ikz,updatetlevel)
            ENDDO
            G_(ikr,ikz) = imagu*kz*G_(ikr,ikz)

          ENDDO kzloopi2
        ENDDO krloopi2

        CALL convolve_2D_F2F( F_, G_, CONV ) ! Convolve and back to Fourier
        Sipj(ip,ij,:,:) = Sipj(ip,ij,:,:) - CONV ! Add it to Sipj (Real space here)

        IF ( in + 1 .LE. jmaxi+1 ) THEN
          factn = real(in,dp) * factn ! factorial(n+1)
        ENDIF

      ENDDO nloopi

      !! Antialiasing
      DO ikr = ikrs,ikre
        DO ikz = ikzs, ikze
          Sipj(ip,ij,ikr,ikz) = Sipj(ip,ij,ikr,ikz) * AA_r(ikr) * AA_z(ikz)
        ENDDO
      ENDDO

    ENDDO jloopi
  ENDDO ploopi
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE compute_Sapj
