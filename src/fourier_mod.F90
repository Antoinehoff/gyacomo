MODULE fourier
  USE prec_const
  USE parallel
  use, intrinsic :: iso_c_binding
  implicit none

  ! INCLUDE 'fftw3.f03'
  INCLUDE 'fftw3-mpi.f03'

  PRIVATE

  PUBLIC :: init_grid_distr_and_plans, poisson_bracket_and_sum, finalize_plans

  real(C_DOUBLE), pointer, PUBLIC            :: real_data_f(:,:), real_data_g(:,:), bracket_sum_r(:,:)
  complex(C_DOUBLE_complex), pointer, PUBLIC :: cmpx_data_f(:,:), cmpx_data_g(:,:), bracket_sum_c(:,:)
  type(C_PTR)                 :: cdatar_f, cdatar_g, cdatar_c
  type(C_PTR)                 :: cdatac_f, cdatac_g, cdatac_c
  type(C_PTR) ,        PUBLIC :: planf, planb
  integer(C_INTPTR_T)         :: i, ix, iy
  integer(C_INTPTR_T), PUBLIC :: alloc_local_1, alloc_local_2
  integer(C_INTPTR_T)         :: NX_, NY_, NY_halved
  ! many plan data variables
  integer(C_INTPTR_T) :: howmany=9 ! numer of element of the tensor
  integer :: rank=3                ! rank of the transform
  integer(C_INTPTR_T), dimension(2) :: fft_dims ! array containing data extent

  CONTAINS

  SUBROUTINE init_grid_distr_and_plans(Nx,Ny,communicator,local_nkx_ptr,local_nkx_ptr_offset,local_nky_ptr,local_nky_ptr_offset)
    IMPLICIT NONE
    INTEGER, INTENT(IN)  :: Nx,Ny, communicator
    INTEGER(C_INTPTR_T), INTENT(OUT) :: local_nkx_ptr,local_nkx_ptr_offset,local_nky_ptr,local_nky_ptr_offset
    NX_ = Nx; NY_ = Ny
    NY_halved = NY_/2 + 1
    !! Complex arrays F, G
    ! Compute the room to allocate
    alloc_local_1 = fftw_mpi_local_size_2d(NY_halved, NX_, communicator, local_nky_ptr, local_nky_ptr_offset)
    ! Initalize pointers to this room
    cdatac_f = fftw_alloc_complex(alloc_local_1)
    cdatac_g = fftw_alloc_complex(alloc_local_1)
    cdatac_c = fftw_alloc_complex(alloc_local_1)
    ! Initalize the arrays with the rooms pointed
    call c_f_pointer(cdatac_f,   cmpx_data_f, [NX_ ,local_nky_ptr])
    call c_f_pointer(cdatac_g,   cmpx_data_g, [NX_ ,local_nky_ptr])
    call c_f_pointer(cdatac_c, bracket_sum_c, [NX_ ,local_nky_ptr])

    !! Real arrays iFFT(F), iFFT(G)
    ! Compute the room to allocate
    alloc_local_2 = fftw_mpi_local_size_2d(NX_, NY_halved, communicator, local_nkx_ptr, local_nkx_ptr_offset)
    ! Initalize pointers to this room
    cdatar_f = fftw_alloc_real(2*alloc_local_2)
    cdatar_g = fftw_alloc_real(2*alloc_local_2)
    cdatar_c = fftw_alloc_real(2*alloc_local_2)
    ! Initalize the arrays with the rooms pointed
    call c_f_pointer(cdatar_f,   real_data_f, [2*(NY_/2 + 1),local_nkx_ptr])
    call c_f_pointer(cdatar_g,   real_data_g, [2*(NY_/2 + 1),local_nkx_ptr])
    call c_f_pointer(cdatar_c, bracket_sum_r, [2*(NY_/2 + 1),local_nkx_ptr])

    ! Plan Creation (out-of-place forward and backward FFT)
    planf = fftw_mpi_plan_dft_r2c_2D(NX_, NY_, real_data_f, cmpx_data_f, communicator,  ior(FFTW_MEASURE, FFTW_MPI_TRANSPOSED_OUT))
    planb = fftw_mpi_plan_dft_c2r_2D(NX_, NY_, cmpx_data_f, real_data_f, communicator,  ior(FFTW_MEASURE, FFTW_MPI_TRANSPOSED_IN))

   if ((.not. c_associated(planf)) .OR. (.not. c_associated(planb))) then
      ERROR STOP '>> ERROR << plan creation error!!'
   end if

  END SUBROUTINE init_grid_distr_and_plans

  !!! Compute the poisson bracket of [F,G] to real space
  !   - Compute the convolution using the convolution theorem
  SUBROUTINE poisson_bracket_and_sum(ky_, kx_, inv_Ny, inv_Nx, AA_y, AA_x,&
                                     local_nky_ptr, local_nkx_ptr, F_, G_, sum_real_)
    IMPLICIT NONE
    INTEGER(C_INTPTR_T),                  INTENT(IN) :: local_nkx_ptr,local_nky_ptr
    REAL(dp),                             INTENT(IN) :: inv_Nx, inv_Ny
    REAL(dp), DIMENSION(local_nky_ptr),   INTENT(IN) :: ky_, AA_y
    REAL(dp), DIMENSION(local_nkx_ptr),   INTENT(IN) :: kx_, AA_x
    COMPLEX(C_DOUBLE_COMPLEX), DIMENSION(local_nky_ptr,local_nkx_ptr), &
                                          INTENT(IN) :: F_(:,:), G_(:,:)
    real(C_DOUBLE), pointer,              INTENT(INOUT) :: sum_real_(:,:)
    INTEGER :: ikx,iky
    ! First term df/dx x dg/dy
    DO ikx = 1,local_nkx_ptr
      DO iky = 1,local_nky_ptr
        cmpx_data_f(ikx,iky) = imagu*kx_(ikx)*F_(iky,ikx)*AA_x(ikx)*AA_y(iky) 
        cmpx_data_g(ikx,iky) = imagu*ky_(iky)*G_(iky,ikx)*AA_x(ikx)*AA_y(iky)
      ENDDO
    ENDDO
    call fftw_mpi_execute_dft_c2r(planb, cmpx_data_f, real_data_f)
    call fftw_mpi_execute_dft_c2r(planb, cmpx_data_g, real_data_g)
    sum_real_ = sum_real_ + real_data_f*real_data_g*inv_Ny*inv_Nx
    ! Second term -df/dy x dg/dx
    DO ikx = 1,local_nkx_ptr
      DO iky = 1,local_nky_ptr
        cmpx_data_f(ikx,iky) = &
              imagu*ky_(iky)*F_(iky,ikx)*AA_x(ikx)*AA_y(iky)
        cmpx_data_g(ikx,iky) = &
              imagu*kx_(ikx)*G_(iky,ikx)*AA_x(ikx)*AA_y(iky)
      ENDDO
    ENDDO
    call fftw_mpi_execute_dft_c2r(planb, cmpx_data_f, real_data_f)
    call fftw_mpi_execute_dft_c2r(planb, cmpx_data_g, real_data_g)
    sum_real_ = sum_real_ - real_data_f*real_data_g*inv_Ny*inv_Nx
END SUBROUTINE poisson_bracket_and_sum


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
