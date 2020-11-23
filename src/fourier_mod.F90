MODULE fourier

  USE prec_const
  USE grid
  USE basic
  use, intrinsic :: iso_c_binding
  implicit none

  ! INCLUDE 'fftw3.f03'
  INCLUDE 'fftw3-mpi.f03'

  PRIVATE
  type(C_PTR)                        :: planf, planb, real_ptr, cmpx_ptr
  real(C_DOUBLE),            pointer :: real_data_1(:,:), real_data_2(:,:)
  complex(C_DOUBLE_complex), pointer :: cmpx_data(:,:)
  integer (C_INTPTR_T)               :: nx, ny, nh

  PUBLIC :: initialize_FFT, finalize_FFT
  PUBLIC :: convolve_2D_F2F


  CONTAINS

    SUBROUTINE initialize_FFT
      IMPLICIT NONE
      nx = Nr; ny = Nz; nh = ( nx / 2 ) + 1

      IF ( my_id .EQ. 1)write(*,*) 'Initialize FFT..'
      CALL fftw_mpi_init

      IF ( my_id .EQ. 1)write(*,*) 'distribution of the array along kr..'
      alloc_local = fftw_mpi_local_size_2d(nh, ny, MPI_COMM_WORLD, local_nkr, local_kr_start)
      write(*,*) 'local_nkr', local_nkr, 'local_kr_start', local_kr_start

      real_ptr = fftw_alloc_real(2 * alloc_local)
      cmpx_ptr = fftw_alloc_complex(alloc_local)

      call c_f_pointer(real_ptr, real_data_1, [2*local_nkr, nx])
      call c_f_pointer(real_ptr, real_data_2, [2*local_nkr, nx])
      call c_f_pointer(cmpx_ptr, cmpx_data,   [  local_nkr, nx])

      IF ( my_id .EQ. 1)write(*,*) 'setting FFT iFFT 2D plans..'
      ! Backward FFT plan
      planb = fftw_mpi_plan_dft_c2r_2d(ny, nx, cmpx_data, real_data_1, MPI_COMM_WORLD, &
                                                                 FFTW_PATIENT)
      ! Forward FFT plan
      planf = fftw_mpi_plan_dft_r2c_2d(ny, nx, real_data_1, cmpx_data, MPI_COMM_WORLD, &
                                                                 FFTW_PATIENT)
    END SUBROUTINE initialize_FFT

    SUBROUTINE finalize_FFT
      IMPLICIT NONE

      IF ( my_id .EQ. 0)write(*,*) 'finalize FFTW'
      CALL dfftw_destroy_plan(planf)
      CALL dfftw_destroy_plan(planb)
      call fftw_mpi_cleanup
      call fftw_free(real_ptr)
      call fftw_free(cmpx_ptr)

    END SUBROUTINE finalize_FFT


  !!! Convolution 2D Fourier to Fourier
  !   - Compute the convolution using the convolution theorem
  SUBROUTINE convolve_2D_F2F( F_2D, G_2D, C_2D )
    IMPLICIT NONE
    COMPLEX(dp), DIMENSION(ikrs:ikre,ikzs:ikze), INTENT(IN)  :: F_2D, G_2D  ! input fields
    COMPLEX(dp), DIMENSION(ikrs:ikre,ikzs:ikze), INTENT(OUT) :: C_2D  ! output convolutioned field

    ! initialize data to some function my_function(i,j)
    do ikr = ikrs,ikre
      do ikz = ikzs,ikze
      cmpx_data(ikr-local_kr_start, ikz) = F_2D(ikr, ikz)
      end do
    end do
    CALL mpi_barrier(MPI_COMM_WORLD, ierr)

    ! 2D backward Fourier transform
    ! write(*,*) my_id, 'iFFT(F) ..'
    call fftw_mpi_execute_dft_c2r(planb, cmpx_data, real_data_1)

    ! initialize data to some function my_function(i,j)
    do ikr = ikrs,ikre
      do ikz = ikzs,ikze
      cmpx_data(ikr-local_kr_start, ikz) = G_2D(ikr, ikz)
      end do
    end do
    ! 2D inverse Fourier transform
    ! write(*,*) my_id, 'iFFT(G) ..'
    call fftw_mpi_execute_dft_c2r(planb, cmpx_data, real_data_2)

    ! Product in physical space
    ! write(*,*) 'multply..'
    real_data_1 = real_data_1 * real_data_2

    ! 2D Fourier tranform
    ! write(*,*) my_id, 'FFT(f*g) ..'
    call fftw_mpi_execute_dft_r2c(planf, real_data_1, cmpx_data)

    ! COPY TO OUTPUT ARRAY
    ! write(*,*) my_id, 'out ..'
    do ikr = ikrs,ikre
      do ikz = ikzs,ikze
        cmpx_data(ikr-local_kr_start,ikz) = C_2D(ikr,ikz)
      end do
    end do

    ! Anti aliasing (2/3 rule)
    DO ikr = ikrs,ikre
      DO ikz = ikzs,ikze
         C_2D(ikr,ikz) = C_2D(ikr,ikz) * AA_r(ikr) * AA_z(ikz)
      ENDDO
    ENDDO

END SUBROUTINE convolve_2D_F2F

END MODULE fourier
