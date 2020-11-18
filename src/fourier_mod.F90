MODULE fourier

  USE prec_const
  USE grid
  use, intrinsic :: iso_c_binding
  implicit none

  ! INCLUDE 'fftw3.f03'
  INCLUDE 'fftw3-mpi.f03'

  PRIVATE
  PUBLIC :: fft_r2cc
  PUBLIC :: ifft_cc2r
  PUBLIC :: fft2_r2cc
  PUBLIC :: ifft2_cc2r
  PUBLIC :: convolve_2D_F2F


  CONTAINS

  SUBROUTINE fft_r2cc(fx_in, Fk_out)

    IMPLICIT NONE
    REAL(dp),    DIMENSION(Nr),  INTENT(IN) :: fx_in
    COMPLEX(dp), DIMENSION(Nkr), INTENT(OUT):: Fk_out
    integer*8 plan

    CALL dfftw_plan_dft_r2c_1d(plan,Nr,fx_in,Fk_out,FFTW_FORWARD,FFTW_ESTIMATE)
    CALL dfftw_execute_dft_r2c(plan, fx_in, Fk_out)
    CALL dfftw_destroy_plan(plan)

  END SUBROUTINE fft_r2cc

  SUBROUTINE ifft_cc2r(Fk_in, fx_out)

    IMPLICIT NONE
    COMPLEX(dp), DIMENSION(Nkr),  INTENT(IN) :: Fk_in
    REAL(dp),    DIMENSION(Nr),   INTENT(OUT):: fx_out
    integer*8 plan

    CALL dfftw_plan_dft_c2r_1d(plan,Nr,Fk_in,fx_out,FFTW_BACKWARD,FFTW_ESTIMATE)
    CALL dfftw_execute_dft_c2r(plan, Fk_in, fx_out)
    CALL dfftw_destroy_plan(plan)
    fx_out = fx_out/Nr

  END SUBROUTINE ifft_cc2r



  SUBROUTINE fft2_r2cc( ffx_in, FFk_out )

      IMPLICIT NONE
      REAL(dp),    DIMENSION(Nr,Nz), INTENT(IN)   :: ffx_in
      COMPLEX(dp), DIMENSION(Nkr,Nkz), INTENT(OUT):: FFk_out
      integer*8 plan

      !!! 2D Forward FFT ________________________!
      CALL dfftw_plan_dft_r2c_2d(plan,Nr,Nz,ffx_in,FFk_out,FFTW_FORWARD,FFTW_ESTIMATE)
      CALL dfftw_execute_dft_r2c(plan,ffx_in,FFk_out)
      CALL dfftw_destroy_plan(plan)

  END SUBROUTINE fft2_r2cc



  SUBROUTINE ifft2_cc2r( FFk_in, ffx_out )

      IMPLICIT NONE
      COMPLEX(dp), DIMENSION(Nkr,Nkz), INTENT(IN) :: FFk_in
      REAL(dp),    DIMENSION(Nr,Nz), INTENT(OUT)  :: ffx_out
      COMPLEX(dp), DIMENSION(Nkr,Nkz) :: tmp_c
      integer*8 plan

      tmp_c = FFk_in
      !!! 2D Backward FFT ________________________!
      CALL dfftw_plan_dft_c2r_2d(plan,Nr,Nz,tmp_c,ffx_out,FFTW_BACKWARD,FFTW_ESTIMATE)
      CALL dfftw_execute_dft_c2r(plan,tmp_c,ffx_out)
      CALL dfftw_destroy_plan(plan)
      ffx_out = ffx_out/Nr/Nz

  END SUBROUTINE ifft2_cc2r


  !!! Convolution 2D Fourier to Fourier
  !   - Compute the convolution using the convolution theorem and MKL
  SUBROUTINE convolve_2D_F2F( F_2D, G_2D, C_2D )
    USE basic
    USE grid, ONLY : AA_r, AA_z, Lr, Lz
    IMPLICIT NONE
    COMPLEX(dp), DIMENSION(Nkr,Nkz), INTENT(IN)  :: F_2D, G_2D  ! input fields
    COMPLEX(dp), DIMENSION(Nkr,Nkz), INTENT(OUT) :: C_2D  ! output convolutioned field
    REAL(dp), DIMENSION(Nr,Nz) :: ff, gg ! iFFT of input fields
    REAL(dp), DIMENSION(Nr,Nz) :: ffgg  ! will receive the product of f*g in Real
    INTEGER :: ikr, ikz
    REAL    :: a_r

    ! 2D inverse Fourier transform
    IF ( num_procs .EQ. 1 ) THEN
      CALL ifft2_cc2r(F_2D,ff)
      CALL ifft2_cc2r(G_2D,gg)
    ELSE
      CALL ifft2_cc2r_mpi(F_2D,ff)
      CALL ifft2_cc2r_mpi(G_2D,gg)
    ENDIF

    ! Product in physical space
    ffgg = ff * gg;

    ! 2D Fourier tranform
    IF ( num_procs .EQ. 1 ) THEN
      CALL fft2_r2cc(ffgg,C_2D)
    ELSE
      CALL fft2_r2cc_mpi(ffgg,C_2D)
    ENDIF

    ! Anti aliasing (2/3 rule)
    DO ikr = ikrs,ikre
      a_r = AA_r(ikr)
      DO ikz = ikzs,ikze
         C_2D(ikr,ikz) = C_2D(ikr,ikz) * a_r * AA_z(ikz)
      ENDDO
    ENDDO

END SUBROUTINE convolve_2D_F2F

!! MPI routines
SUBROUTINE fft2_r2cc_mpi( ffx_in, FFk_out )

    IMPLICIT NONE
    REAL(dp),    DIMENSION(Nr,Nz), INTENT(IN)   :: ffx_in
    COMPLEX(dp), DIMENSION(Nkr,Nkz), INTENT(OUT):: FFk_out
    REAL(dp),    DIMENSION(Nkr,Nkz) :: tmp_r
    type(C_PTR) :: plan
    integer(C_INTPTR_T) :: L
    integer(C_INTPTR_T) :: M

    ! L = Nr; M = Nz
    !
    ! tmp_r = ffx_in
    ! !!! 2D Forward FFT ________________________!
    ! plan = fftw_mpi_plan_dft_r2c_2d(L, M, tmp_r, FFk_out, MPI_COMM_WORLD, FFTW_ESTIMATE)
    !
    ! if ((.not. c_associated(plan))) then
    !  write(*,*) "plan creation error!!"
    !  stop
    ! end if
    !
    ! CALL fftw_mpi_execute_dft_r2c(plan,tmp_r,FFk_out)
    ! CALL fftw_destroy_plan(plan)

END SUBROUTINE fft2_r2cc_mpi



SUBROUTINE ifft2_cc2r_mpi( FFk_in, ffx_out )

    IMPLICIT NONE
    COMPLEX(dp), DIMENSION(Nkr,Nkz), INTENT(IN) :: FFk_in
    REAL(dp),    DIMENSION(Nr,Nz), INTENT(OUT)  :: ffx_out
    COMPLEX(dp), DIMENSION(Nkr,Nkz) :: tmp_c
    type(C_PTR) :: plan
    integer(C_INTPTR_T) :: L
    integer(C_INTPTR_T) :: M

    ! L = Nr; M = Nz
    !
    ! tmp_c = FFk_in
    ! !!! 2D Backward FFT ________________________!
    ! plan = fftw_mpi_plan_dft_c2r_2d(L, M, tmp_c, ffx_out, MPI_COMM_WORLD, FFTW_ESTIMATE)
    !
    ! if ((.not. c_associated(plan))) then
    !  write(*,*) "plan creation error!!"
    !  stop
    ! end if
    !
    ! CALL fftw_mpi_execute_dft_c2r(plan,tmp_c,ffx_out)
    ! CALL fftw_destroy_plan(plan)
    !
    ! ffx_out = ffx_out/Nr/Nz

END SUBROUTINE ifft2_cc2r_mpi

! Empty set/free routines to switch easily with MKL DFTI
SUBROUTINE set_descriptors

END SUBROUTINE set_descriptors

SUBROUTINE free_descriptors

END SUBROUTINE free_descriptors

END MODULE fourier
