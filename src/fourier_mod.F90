MODULE fourier

  USE prec_const
  USE grid
  USE basic
  use, intrinsic :: iso_c_binding
  implicit none

  ! INCLUDE 'fftw3.f03'
  INCLUDE 'fftw3-mpi.f03'

  PRIVATE

  PUBLIC :: init_grid_distr_and_plans, convolve_2D_F2F, finalize_plans

  real(C_DOUBLE), pointer, PUBLIC            :: real_data_f(:,:), real_data_g(:,:), real_data_c(:,:)
  complex(C_DOUBLE_complex), pointer, PUBLIC :: cmpx_data_f(:,:), cmpx_data_g(:,:), cmpx_data_c(:,:)
  type(C_PTR)                            :: cdatar_f, cdatar_g, cdatar_c
  type(C_PTR)                            :: cdatac_f, cdatac_g, cdatac_c
  type(C_PTR) , PUBLIC                   :: planf, planb
  integer(C_INTPTR_T)                    :: i, ix, iy
  integer(C_INTPTR_T), PUBLIC            :: alloc_local_1, alloc_local_2
  integer(C_INTPTR_T)                    :: NR_, NZ_, NR_halved

  ! many plan data variables
  integer(C_INTPTR_T) :: howmany=9 ! numer of eleemnt of the tensor
  integer :: rank=3                ! rank of the transform
  integer(C_INTPTR_T), dimension(2) :: fft_dims ! array containing data extent

  CONTAINS

  SUBROUTINE init_grid_distr_and_plans(Nr,Nz)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: Nr,Nz
    NR_ = Nr; NZ_ = Nz
    NR_halved = NR_/2 + 1

    !! Complex arrays F, G
    ! Compute the room to allocate
    alloc_local_1 = fftw_mpi_local_size_2d(NR_halved, NZ_, MPI_COMM_WORLD, local_nkr, local_nkr_offset)
    ! alloc_local_1 = fftw_mpi_local_size_2d_many(2, NR_halved, NZ_, MPI_COMM_WORLD, local_nkr, local_nkr_offset)
    ! Initalize pointers to this room
    cdatac_f = fftw_alloc_complex(alloc_local_1)
    cdatac_g = fftw_alloc_complex(alloc_local_1)
    cdatac_c = fftw_alloc_complex(alloc_local_1)
    ! Initalize the arrays with the rooms pointed
    call c_f_pointer(cdatac_f, cmpx_data_f, [NZ_ ,local_nkr])
    call c_f_pointer(cdatac_g, cmpx_data_g, [NZ_ ,local_nkr])
    call c_f_pointer(cdatac_c, cmpx_data_c, [NZ_ ,local_nkr])

    !! Real arrays iFFT(F), iFFT(G)
    ! Compute the room to allocate
    alloc_local_2 = fftw_mpi_local_size_2d(NZ_, NR_halved, MPI_COMM_WORLD, local_nz, local_nz_offset)
    ! alloc_local_2 = fftw_mpi_local_size_2d_many(2, NZ_, NR_halved, MPI_COMM_WORLD, local_nz, local_nz_offset)
    ! Initalize pointers to this room
    cdatar_f = fftw_alloc_real(2*alloc_local_2)
    cdatar_g = fftw_alloc_real(2*alloc_local_2)
    cdatar_c = fftw_alloc_real(2*alloc_local_2)
    ! Initalize the arrays with the rooms pointed
    call c_f_pointer(cdatar_f, real_data_f, [2*(NR_/2  + 1),local_nz])
    call c_f_pointer(cdatar_g, real_data_g, [2*(NR_/2  + 1),local_nz])
    call c_f_pointer(cdatar_c, real_data_c, [2*(NR_/2  + 1),local_nz])

    ! Plan Creation (out-of-place forward and backward FFT)
    planf = fftw_mpi_plan_dft_r2c_2d(NZ_, NR_, real_data_f, cmpx_data_f, MPI_COMM_WORLD,  ior(FFTW_MEASURE, FFTW_MPI_TRANSPOSED_OUT))
    planb = fftw_mpi_plan_dft_c2r_2d(NZ_, NR_, cmpx_data_f, real_data_f, MPI_COMM_WORLD,  ior(FFTW_MEASURE, FFTW_MPI_TRANSPOSED_IN))

   if ((.not. c_associated(planf)) .OR. (.not. c_associated(planb))) then
      write(*,*) "plan creation error!!"
      stop
   end if

  END SUBROUTINE init_grid_distr_and_plans


  !!! Convolution 2D Fourier to Fourier
  !   - Compute the convolution using the convolution theorem and MKL
  SUBROUTINE convolve_2D_F2F( F_2D, G_2D, C_2D )
    IMPLICIT NONE

    COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(ikrs:ikre,ikzs:ikze), INTENT(IN)  :: F_2D, G_2D  ! input fields
    COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(ikrs:ikre,ikzs:ikze), INTENT(OUT) :: C_2D  ! output convolutioned field

    do ikr = ikrs, ikre
      do ikz = ikzs, ikze
        cmpx_data_f(ikz,ikr-local_nkr_offset) = F_2D(ikr,ikz)
        cmpx_data_g(ikz,ikr-local_nkr_offset) = G_2D(ikr,ikz)
      end do
    end do

    call fftw_mpi_execute_dft_c2r(planb, cmpx_data_f, real_data_f)

    call fftw_mpi_execute_dft_c2r(planb, cmpx_data_g, real_data_g)

    real_data_c = real_data_f/NZ_/NR_  * real_data_g/NZ_/NR_

    call fftw_mpi_execute_dft_r2c(planf, real_data_c, cmpx_data_c)

    ! Retrieve convolution in input format
    do ikr = ikrs, ikre
      do ikz = ikzs, ikze
        C_2D(ikr,ikz) = cmpx_data_c(ikz,ikr-local_nkr_offset)*AA_r(ikr)*AA_z(ikz)
      end do
    end do

END SUBROUTINE convolve_2D_F2F

SUBROUTINE finalize_plans
  IMPLICIT NONE

  IF (my_id .EQ. 0) write(*,*) '..plan Destruction.'
  call fftw_destroy_plan(planb)
  call fftw_destroy_plan(planf)
  call fftw_mpi_cleanup()
  call fftw_free(cdatar_f)
  call fftw_free(cdatar_g)
  call fftw_free(cdatar_c)
  call fftw_free(cdatac_f)
  call fftw_free(cdatac_g)
  call fftw_free(cdatac_c)
END SUBROUTINE finalize_plans

END MODULE fourier
