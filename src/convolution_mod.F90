MODULE convolution
  USE prec_const
  implicit none

  CONTAINS

  !!! Convolution 2D Fourier to Fourier
  !   - Compute the convolution using the convolution theorem and MKL
  SUBROUTINE convolve_2D_F2F( F_2D, G_2D, C_2D )

      USE prec_const
      USE fourier_grid, ONLY : Nkr, Nkz, Pad
      USE MKL_DFTI

      IMPLICIT NONE

      COMPLEX(dp), DIMENSION(Nkr,Nkz)     :: F_2D, G_2D, C_2D
      COMPLEX(dp), DIMENSION(Pad*Nkr*Pad*Nkz) :: F_1D, G_1D, C_1D

      INTEGER :: ix, iy, Mkr, Mkz
      INTEGER :: Status, L(2), L_pad(2)

      type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle,My_Desc2_Handle,My_Desc3_Handle

      L     = (/ Nkr, Nkz /)
      Mkr = Pad*Nkr; Mkz = Pad*Nkz
      L_pad = (/ Mkr, Mkz /)

      !WRITE(*,*) 'Reshaping..'
      DO ix = 1,Mkr
        DO iy = 1,Mkz
          IF ( (ix .LE. Nkr) .AND. (iy .LE. Nkz) ) THEN
              F_1D(Mkr*(iy-1) + ix) = F_2D(ix,iy)
              G_1D(Mkr*(iy-1) + ix) = G_2D(ix,iy)
          ELSE
              F_1D(Mkr*(iy-1) + ix) = 0._dp
              G_1D(Mkr*(iy-1) + ix) = 0._dp
          ENDIF
        ENDDO
      ENDDO

      Status = DftiCreateDescriptor(My_Desc1_Handle, DFTI_DOUBLE,DFTI_COMPLEX, 2, L_pad)
      Status = DftiCommitDescriptor(My_Desc1_Handle)
      !WRITE(*,*) 'Backward FFT on F..'
      Status = DftiComputeBackward (My_Desc1_Handle, F_1D)
      Status = DftiFreeDescriptor  (My_Desc1_Handle)

      Status = DftiCreateDescriptor(My_Desc2_Handle, DFTI_DOUBLE,DFTI_COMPLEX, 2, L_pad)
      Status = DftiCommitDescriptor(My_Desc2_Handle)
      !WRITE(*,*) 'Backward FFT on G..'
      Status = DftiComputeBackward (My_Desc2_Handle, G_1D)
      Status = DftiFreeDescriptor  (My_Desc2_Handle)

      !CALL backward_FFT( F_2D )
      !CALL backward_FFT( G_2D )
      
      !WRITE(*,*) 'C =  F/(Mkr *Mkz) x G/(Mkr *Mkz)..'
      DO ix=1,Mkr*Mkz
        C_1D(ix) = F_1D(ix)/Mkr/Mkz * G_1D(ix)/Mkr/Mkz
      ENDDO

      !WRITE(*,*) 'Forward FFT on C..'
      Status = DftiCreateDescriptor( My_Desc3_Handle, DFTI_DOUBLE,&
                DFTI_COMPLEX, 2, L_pad)
      Status = DftiCommitDescriptor(My_Desc3_Handle)
      Status = DftiComputeForward  (My_Desc3_Handle, C_1D)
      Status = DftiFreeDescriptor  (My_Desc3_Handle)

      DO ix=1,Nkr
        DO iy=1,Nkz
          C_2D(ix,iy) = C_1D(Mkr*(iy-1) + ix)
        ENDDO
      ENDDO

  END SUBROUTINE convolve_2D_F2F

  !!! Convolution 2D Fourier to Real
  !   - Same as convolve_2D_F2F but does not compute the forward FFT at the end
  !   - Used to same computation in the compute_Sapj loop
  !   - Output : C_2D that is defined in the real space
  SUBROUTINE convolve_2D_F2R( F_2D, G_2D, C_2D )

    USE prec_const
    USE fourier_grid, ONLY : Nkr, Nkz, Pad
    USE MKL_DFTI

    IMPLICIT NONE

    COMPLEX(dp), DIMENSION(Nkr,Nkz)     :: F_2D, G_2D, C_2D
    COMPLEX(dp), DIMENSION(Pad*Nkr*Pad*Nkz) :: F_1D, G_1D, C_1D

    INTEGER :: ix, iy, Mkr, Mkz
    INTEGER :: Status, L(2), L_pad(2)

    type(DFTI_DESCRIPTOR), POINTER :: My_Desc1_Handle,My_Desc2_Handle

    L     = (/ Nkr, Nkz /)
    Mkr = Pad*Nkr; Mkz = Pad*Nkz
    L_pad = (/ Mkr, Mkz /)

    !WRITE(*,*) 'Reshaping..'
    DO ix = 1,Mkr
      DO iy = 1,Mkz
        IF ( (ix .LE. Nkr) .AND. (iy .LE. Nkz) ) THEN
            F_1D(Mkr*(iy-1) + ix) = F_2D(ix,iy)
            G_1D(Mkr*(iy-1) + ix) = G_2D(ix,iy)
        ELSE
            F_1D(Mkr*(iy-1) + ix) = 0._dp
            G_1D(Mkr*(iy-1) + ix) = 0._dp
        ENDIF
      ENDDO
    ENDDO

    CALL backward_FFT( F_2D )
    CALL backward_FFT( G_2D )

    ! We adapt the 1D output vector to a 2D form
    DO ix=1,Nkr
      DO iy=1,Nkz
        !C_2D(ix,iy) = F_1D(Mkr*(iy-1) + ix)/Mkr/Mkz * G_1D(Mkr*(iy-1) + ix)/Mkr/Mkz
        C_2D(ix,iy) = F_2D(ix,iy) * G_2D(ix,iy)
      ENDDO
    ENDDO

  END SUBROUTINE convolve_2D_F2R


  !! Compute a forward FFT to go back to Fourier space 
  !  - used to save computation after the sum of convolution_2D_F2R in compute_Sapj
  SUBROUTINE forward_FFT( C_2D )

    USE prec_const
    USE fourier_grid, ONLY : Nkr, Nkz, Pad
    USE MKL_DFTI

    IMPLICIT NONE

    COMPLEX(dp), DIMENSION(Nkr,Nkz)         :: C_2D
    COMPLEX(dp), DIMENSION(Pad*Nkr*Pad*Nkz) :: C_1D

    INTEGER :: ix, iy, Mkr, Mkz
    INTEGER :: Status, L(2), L_pad(2)

    type(DFTI_DESCRIPTOR), POINTER :: My_Desc_Handle

    L     = (/ Nkr, Nkz /)
    Mkr = Pad*Nkr; Mkz = Pad*Nkz
    L_pad = (/ Mkr, Mkz /)

    !WRITE(*,*) 'Reshaping..'
    DO ix = 1,Mkr
      DO iy = 1,Mkz
        IF ( (ix .LE. Nkr) .AND. (iy .LE. Nkz) ) THEN
            C_1D(Mkr*(iy-1) + ix) = C_2D(ix,iy)
        ELSE
            C_1D(Mkr*(iy-1) + ix) = 0._dp
        ENDIF
      ENDDO
    ENDDO

    Status = DftiCreateDescriptor( My_Desc_Handle, DFTI_DOUBLE,&
              DFTI_COMPLEX, 2, L_pad)
    Status = DftiCommitDescriptor(My_Desc_Handle)
    Status = DftiComputeForward  (My_Desc_Handle, C_1D)
    Status = DftiFreeDescriptor  (My_Desc_Handle)

    DO ix=1,Nkr
      DO iy=1,Nkz
        C_2D(ix,iy) = C_1D(Mkr*(iy-1) + ix)
      ENDDO
    ENDDO

  END SUBROUTINE forward_FFT

    !!! Convolution 2D Fourier to Fourier
  !   - Compute the convolution using the convolution theorem and MKL
  SUBROUTINE backward_FFT( C_2D )

    USE prec_const
    USE fourier_grid, ONLY : Nkr, Nkz, Pad
    USE MKL_DFTI

    IMPLICIT NONE

    COMPLEX(dp), DIMENSION(Nkr,Nkz)     :: C_2D
    COMPLEX(dp), DIMENSION(Pad*Nkr*Pad*Nkz) ::  C_1D

    INTEGER :: ix, iy, Mkr, Mkz
    INTEGER :: Status, L(2), L_pad(2)

    type(DFTI_DESCRIPTOR), POINTER :: My_Desc_Handle

    L     = (/ Nkr, Nkz /)
    Mkr = Pad*Nkr; Mkz = Pad*Nkz
    L_pad = (/ Mkr, Mkz /)

    !WRITE(*,*) 'Reshaping..'
    DO ix = 1,Mkr
      DO iy = 1,Mkz
        IF ( (ix .LE. Nkr) .AND. (iy .LE. Nkz) ) THEN
            C_1D(Mkr*(iy-1) + ix) = C_2D(ix,iy)
        ELSE
            C_1D(Mkr*(iy-1) + ix) = 0._dp
        ENDIF
      ENDDO
    ENDDO

    Status = DftiCreateDescriptor(My_Desc_Handle, DFTI_DOUBLE,DFTI_COMPLEX, 2, L_pad)
    Status = DftiCommitDescriptor(My_Desc_Handle)
    Status = DftiComputeBackward (My_Desc_Handle, C_1D)
    Status = DftiFreeDescriptor  (My_Desc_Handle)

    DO ix=1,Nkr
      DO iy=1,Nkz
        C_2D(ix,iy) = C_1D(Mkr*(iy-1) + ix)/Mkr/Mkz
      ENDDO
    ENDDO

  END SUBROUTINE backward_FFT

END MODULE convolution
