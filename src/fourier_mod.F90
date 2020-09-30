MODULE fourier

  USE prec_const
  USE grid
  use, intrinsic :: iso_c_binding
  implicit none

  INCLUDE 'fftw3.f03'

  PRIVATE
  PUBLIC :: fft_r2cc
  PUBLIC :: ifft_cc2r
  PUBLIC :: fft2_r2cc
  PUBLIC :: ifft2_cc2r
  PUBLIC :: convolve_2D_F2F
  PUBLIC :: set_descriptors, free_descriptors ! Temporary to switch easily with MKL DFTI

  CONTAINS

  SUBROUTINE fft_r2cc(fx_in, Fk_out)

    IMPLICIT NONE
    REAL(dp),    DIMENSION(Nr),  INTENT(IN) :: fx_in
    COMPLEX(dp), DIMENSION(Nkr), INTENT(OUT):: Fk_out
    integer*8 plan

    call dfftw_plan_dft_r2c_1d(plan,Nr,fx_in,Fk_out,FFTW_FORWARD,FFTW_ESTIMATE)
    call dfftw_execute_dft_r2c(plan, fx_in, Fk_out)
    call dfftw_destroy_plan(plan)

  END SUBROUTINE fft_r2cc



  SUBROUTINE ifft_cc2r(Fk_in, fx_out)

    IMPLICIT NONE
    COMPLEX(dp), DIMENSION(Nkr),  INTENT(IN) :: Fk_in
    REAL(dp),    DIMENSION(Nr),   INTENT(OUT):: fx_out
    integer*8 plan

    call dfftw_plan_dft_c2r_1d(plan,Nr,Fk_in,fx_out,FFTW_BACKWARD,FFTW_ESTIMATE)
    call dfftw_execute_dft_c2r(plan, Fk_in, fx_out)
    call dfftw_destroy_plan(plan)
    fx_out = fx_out/Nr

  END SUBROUTINE ifft_cc2r



  SUBROUTINE fft2_r2cc( ffx_in, FFk_out )

      IMPLICIT NONE
      REAL(dp),    DIMENSION(Nr,Nz), INTENT(IN)   :: ffx_in
      COMPLEX(dp), DIMENSION(Nkr,Nkz), INTENT(OUT):: FFk_out
      integer*8 plan

      !!! 2D Forward FFT ________________________!
      call dfftw_plan_dft_r2c_2d(plan,Nr,Nz,ffx_in,FFk_out,FFTW_FORWARD,FFTW_ESTIMATE)
      call dfftw_execute_dft_r2c(plan,ffx_in,FFk_out)
      call dfftw_destroy_plan(plan)

  END SUBROUTINE fft2_r2cc



  SUBROUTINE ifft2_cc2r( FFk_in, ffx_out )

      IMPLICIT NONE
      COMPLEX(dp), DIMENSION(Nkr,Nkz), INTENT(IN) :: FFk_in
      REAL(dp),    DIMENSION(Nr,Nz), INTENT(OUT)  :: ffx_out
      COMPLEX(dp), DIMENSION(Nkr,Nkz) :: tmp_c
      integer*8 plan

      tmp_c = FFk_in
      !!! 2D Backward FFT ________________________!
      call dfftw_plan_dft_c2r_2d(plan,Nr,Nz,tmp_c,ffx_out,FFTW_BACKWARD,FFTW_ESTIMATE)
      call dfftw_execute_dft_c2r(plan,tmp_c,ffx_out)
      call dfftw_destroy_plan(plan)
      ffx_out = ffx_out/Nr/Nz

  END SUBROUTINE ifft2_cc2r



  !!! Convolution 2D Fourier to Fourier
  !   - Compute the convolution using the convolution theorem and MKL
  SUBROUTINE convolve_2D_F2F( F_2D, G_2D, C_2D )

    USE grid, ONLY : AA_r, AA_z
    IMPLICIT NONE
    COMPLEX(dp), DIMENSION(Nkr,Nkz), INTENT(IN)  :: F_2D, G_2D  ! input fields
    COMPLEX(dp), DIMENSION(Nkr,Nkz), INTENT(OUT) :: C_2D  ! output convolutioned field
    REAL(dp), DIMENSION(Nr,Nz) :: ff, gg ! iFFT of input fields
    REAL(dp), DIMENSION(Nr,Nz) :: ffgg  ! will receive the product of f*g in Real
    INTEGER :: ikr, ikz
    REAL    :: a_r

    ! 2D inverse Fourier transform
    CALL ifft2_cc2r(F_2D,ff);
    CALL ifft2_cc2r(G_2D,gg);

    ! Product in physical space
    ffgg = ff * gg;

    ! 2D Fourier tranform
    CALL fft2_r2cc(ffgg,C_2D);

    ! Anti aliasing (2/3 rule)
    DO ikr = 1,Nkr
      a_r = AA_r(ikr)
      DO ikz = 1,Nkz
         C_2D(ikr,ikz) = C_2D(ikr,ikz) * a_r * AA_z(ikz)
      ENDDO
    ENDDO

END SUBROUTINE convolve_2D_F2F


! Empty set/free routines to switch easily with MKL DFTI
SUBROUTINE set_descriptors

END SUBROUTINE set_descriptors

SUBROUTINE free_descriptors

END SUBROUTINE free_descriptors

END MODULE fourier
