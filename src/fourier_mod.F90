MODULE fourier
    USE prec_const, ONLY: xp, c_xp_c, c_xp_r, imagu
    use, intrinsic :: iso_c_binding
    implicit none

    ! INCLUDE 'fftw3.f03'
    INCLUDE 'fftw3-mpi.f03'

    PRIVATE

    !! Module parameter
    ! 2D fft routines or 1D methods (for ExBshear simulations, 1D is required)
    ! The 2D fft is using fftw mpi optimization
    ! The 1D is not using mpi and does a data transfer with redundant computations
    !   (each process computes the convolution)
    LOGICAL, PUBLIC, PROTECTED :: FFT2D = .TRUE.

    !! Module accessible routines
    PUBLIC :: init_grid_distr_and_plans, poisson_bracket_and_sum, finalize_plans

    !! Module variables
    CHARACTER(2)                :: FFT_ALGO ! use of 2D or 1D routines
    type(C_PTR)                 :: cdatar_f, cdatar_g, cdatar_c
    type(C_PTR)                 :: cdatac_f, cdatac_g, cdatac_c
    type(C_PTR) ,        PUBLIC :: planf, planb
    integer(C_INTPTR_T)         :: i, ix, iy
    integer(C_INTPTR_T), PUBLIC :: alloc_local_1, alloc_local_2
    integer(C_INTPTR_T)         :: NX_, NY_, NY_halved 
    real   (c_xp_r), pointer, PUBLIC :: real_data_f(:,:), real_data_g(:,:), bracket_sum_r(:,:)
    complex(c_xp_c), pointer, PUBLIC :: cmpx_data_f(:,:), cmpx_data_g(:,:), bracket_sum_c(:,:)
    !! 1D fft specific variables
    type(C_PTR), PUBLIC :: plan_kx2x_c2c ! transform from (kx,ky) to ( x,ky) (complex to complex)
    type(C_PTR), PUBLIC :: plan_ky2y_c2r ! transform from ( x,ky) to ( x, y) (complex to real)
    type(C_PTR), PUBLIC :: plan_y2ky_r2c ! transform from ( x, y) to ( x,ky) (real to complex)
    type(C_PTR), PUBLIC :: plan_x2kx_c2c ! transform from ( x,ky) to (kx,ky) (complex to complex)
    complex(c_xp_c), pointer, PUBLIC :: ky_x_data(:,:)

    CONTAINS
    !******************************************************************************!
    !------------- Initialize the grid distribution and plans -------------
    ! If we use the 2D fftw routines, the fftw library decides the best data distribution
    SUBROUTINE init_grid_distr_and_plans(FFT2D,Nx,Ny,communicator,&
                    local_nkx_ptr,local_nkx_ptr_offset,local_nky_ptr,local_nky_ptr_offset)
        USE basic, ONLY: speak
        IMPLICIT NONE
        LOGICAL, INTENT(IN)  :: FFT2D
        INTEGER, INTENT(IN)  :: Nx,Ny, communicator
        INTEGER(C_INTPTR_T), INTENT(OUT) :: local_nkx_ptr,local_nkx_ptr_offset
        INTEGER(C_INTPTR_T), INTENT(OUT) :: local_nky_ptr,local_nky_ptr_offset
        NX_ = Nx; NY_ = Ny
        NY_halved = NY_/2 + 1
        IF(FFT2D) THEN
            FFT_ALGO = '2D'
        ELSE
            FFT_ALGO = '1D'
        ENDIF
        CALL speak('FFT algorithm :' // FFT_ALGO)

        SELECT CASE (FFT_ALGO)
        CASE ('2D')
            CALL fft2D_distr_and_plans(Nx,Ny,communicator,&
                    local_nkx_ptr,local_nkx_ptr_offset,local_nky_ptr,local_nky_ptr_offset)
        CASE ('1D')
            CALL fft1D_distr_and_plans(Nx,Ny,communicator,&
                    local_nkx_ptr,local_nkx_ptr_offset,local_nky_ptr,local_nky_ptr_offset)
        END SELECT
    END SUBROUTINE init_grid_distr_and_plans

    !------------- 2D fft initialization and mpi distribution
    SUBROUTINE fft2D_distr_and_plans(Nx,Ny,communicator,&
                local_nkx_ptr,local_nkx_ptr_offset,local_nky_ptr,local_nky_ptr_offset)
        IMPLICIT NONE
        INTEGER, INTENT(IN)  :: Nx,Ny, communicator
        INTEGER(C_INTPTR_T), INTENT(OUT) :: local_nkx_ptr,local_nkx_ptr_offset
        INTEGER(C_INTPTR_T), INTENT(OUT) :: local_nky_ptr,local_nky_ptr_offset
        !! Complex arrays F, G
        ! Compute the room to allocate
#ifdef SINGLE_PRECISION
        alloc_local_1 = fftwf_mpi_local_size_2d(NY_halved, NX_,communicator,&
                local_nky_ptr, local_nky_ptr_offset)
#else
        alloc_local_1 = fftw_mpi_local_size_2d(NY_halved, NX_,communicator,&
                local_nky_ptr, local_nky_ptr_offset)
#endif
    ! Initalize pointers to this room
#ifdef SINGLE_PRECISION
        cdatac_f = fftwf_alloc_complex(alloc_local_1)
        cdatac_g = fftwf_alloc_complex(alloc_local_1)
        cdatac_c = fftwf_alloc_complex(alloc_local_1)
#else
        cdatac_f = fftw_alloc_complex(alloc_local_1)
        cdatac_g = fftw_alloc_complex(alloc_local_1)
        cdatac_c = fftw_alloc_complex(alloc_local_1)
#endif
        ! Initalize the arrays with the rooms pointed
        call c_f_pointer(cdatac_f,   cmpx_data_f, [NX_ ,local_nky_ptr])
        call c_f_pointer(cdatac_g,   cmpx_data_g, [NX_ ,local_nky_ptr])
        call c_f_pointer(cdatac_c, bracket_sum_c, [NX_ ,local_nky_ptr])
        !! Real arrays iFFT(F), iFFT(G)
        ! Compute the room to allocate
        alloc_local_2 = fftw_mpi_local_size_2d(NX_,NY_halved,communicator,&
                local_nkx_ptr, local_nkx_ptr_offset)
        ! Initalize pointers to this room
#ifdef SINGLE_PRECISION
        cdatar_f = fftwf_alloc_real(2*alloc_local_2)
        cdatar_g = fftwf_alloc_real(2*alloc_local_2)
        cdatar_c = fftwf_alloc_real(2*alloc_local_2)
#else
        cdatar_f = fftw_alloc_real(2*alloc_local_2)
        cdatar_g = fftw_alloc_real(2*alloc_local_2)
        cdatar_c = fftw_alloc_real(2*alloc_local_2)
#endif
        ! Initalize the arrays with the rooms pointed
        call c_f_pointer(cdatar_f,   real_data_f, [2*(NY_/2 + 1),local_nkx_ptr])
        call c_f_pointer(cdatar_g,   real_data_g, [2*(NY_/2 + 1),local_nkx_ptr])
        call c_f_pointer(cdatar_c, bracket_sum_r, [2*(NY_/2 + 1),local_nkx_ptr])
        ! Plan Creation (out-of-place forward and backward FFT)
#ifdef SINGLE_PRECISION
        planf = fftwf_mpi_plan_dft_r2c_2D(NX_,NY_,real_data_f,cmpx_data_f,communicator, &
                                            ior(FFTW_MEASURE, FFTW_MPI_TRANSPOSED_OUT))
        planb = fftwf_mpi_plan_dft_c2r_2D(NX_,NY_,cmpx_data_f,real_data_f,communicator, &
                                            ior(FFTW_MEASURE, FFTW_MPI_TRANSPOSED_IN))
#else
        planf = fftw_mpi_plan_dft_r2c_2D(NX_,NY_,real_data_f,cmpx_data_f,communicator, &
                                            ior(FFTW_MEASURE, FFTW_MPI_TRANSPOSED_OUT))
        planb = fftw_mpi_plan_dft_c2r_2D(NX_,NY_,cmpx_data_f,real_data_f,communicator, &
                                            ior(FFTW_MEASURE, FFTW_MPI_TRANSPOSED_IN))
#endif
        if ((.not. c_associated(planf)) .OR. (.not. c_associated(planb))) then
            ERROR STOP '>> ERROR << plan creation error!!'
        end if
    END SUBROUTINE fft2D_distr_and_plans
    !******************************************************************************!

    !******************************************************************************!
    !------------- 1D initialization with balanced data distribution
    SUBROUTINE fft1D_distr_and_plans(Nx,Ny,communicator,&
                local_nkx_ptr,local_nkx_ptr_offset,local_nky_ptr,local_nky_ptr_offset)
        USE utility,  ONLY: decomp1D
        USE parallel, ONLY: num_procs_ky, rank_ky
        IMPLICIT NONE
        INTEGER, INTENT(IN)  :: Nx,Ny, communicator
        INTEGER(C_INTPTR_T), INTENT(OUT) :: local_nkx_ptr,local_nkx_ptr_offset
        INTEGER(C_INTPTR_T), INTENT(OUT) :: local_nky_ptr,local_nky_ptr_offset
        ! local var
        INTEGER :: is,ie    !start and end indices
        INTEGER :: rank     ! rank of each 1D fourier transforms
        INTEGER :: n        ! size of the data to fft
        INTEGER :: howmany  ! howmany 1D fourier transforms
        COMPLEX, DIMENSION(:,:), ALLOCATABLE:: in, out
        INTEGER :: inembed, onembed
        INTEGER :: istride, ostride
        INTEGER :: idist, odist
        INTEGER :: sign
        INTEGER :: flags
        COMPLEX(xp), DIMENSION(:,:), ALLOCATABLE:: fkxky
        COMPLEX(xp), DIMENSION(:,:), ALLOCATABLE:: fxky_l, fxky_g
        REAL(xp),    DIMENSION(:,:), ALLOCATABLE:: fxy, fkxy

        ! number of kx (no data distr.)
        local_nkx_ptr        = Nx
        local_nkx_ptr_offset = 0

        !! Distributinon of data and definition of the size of the arrays that will be worked on
        ! balanced distribution among the processes for ky
        CALL  decomp1D( (Ny/2+1), num_procs_ky, rank_ky, is, ie )
        local_nky_ptr        = ie - is + 1
        local_nky_ptr_offset = is - 1
        ! give the rest of the points to the last process
        if (rank_ky .EQ. num_procs_ky-1) local_nky_ptr = (Ny/2+1)-local_nky_ptr_offset

        !! Allocate temporary array for plan making
        ALLOCATE(  fkxky(Nx,local_nky_ptr))
        ALLOCATE( fxky_l(Nx,local_nky_ptr))
        ALLOCATE( fxky_g(Nx,Ny/2+1))
        ALLOCATE(    fxy(Nx,Ny))
        !! Plan of the 4 many transforms required
        ! 1. (kx,ky) -> (x,ky), C -> C, transforms
        rank    = 1              ! 1D transform
        n       = Nx             ! all kx modes
        howmany = local_nky_ptr  ! all local ky
        inembed = Nx             ! all data must be transformed
        onembed = Nx
        idist   = Nx             ! distance between data to transforms (x columns)
        odist   = Nx
        istride = 1              ! contiguous data
        ostride = 1
#ifdef SINGLE_PRECISION
        CALL sfftw_plan_many_dft(plan_kx2x_c2c, rank, n, howmany,&
                                 fkxky, inembed, istride, idist,&
                                fxky_l, onembed, ostride, odist,& 
                                 FFTW_BACKWARD, FFTW_PATIENT)                
#else
        CALL dfftw_plan_many_dft(plan_kx2x_c2c, rank, n, howmany,&
                                 fkxky, inembed, istride, idist,&
                                fxky_l, onembed, ostride, odist,& 
                                 FFTW_BACKWARD, FFTW_PATIENT)    
#endif
        ! 1.5 MPI communication along ky (from fxky_l to fxky_g)
        ! 2. (x,ky) -> (x,y), C -> R, transforms
        rank    = 1              ! 1D transform
        n       = Ny             ! all ky modes
        howmany = Nx             ! all kx
        inembed = Ny/2+1         ! all ky must be transformed
        onembed = Ny             ! to all y
        idist   = 1              ! distance between two slice to transforms (y row)
        odist   = 1
        istride = Nx             ! non contiguous data
        ostride = Nx
#ifdef SINGLE_PRECISION
        CALL sfftw_plan_many_dft_c2r(plan_ky2y_c2r, rank, n, howmany,&
                                   fxky_g, inembed, istride, idist,&
                                      fxy, onembed, ostride, odist,& 
                                     FFTW_BACKWARD, FFTW_PATIENT)                
#else
        CALL dfftw_plan_many_dft_c2r(plan_ky2y_c2r, rank, n, howmany,&
                                   fxky_g, inembed, istride, idist,&
                                      fxy, onembed, ostride, odist,& 
                                     FFTW_BACKWARD, FFTW_PATIENT)    
#endif
        ! 3. (x,y) -> (x,ky), R -> C, transforms
        rank    = 1              ! 1D transform
        n       = Ny             ! all y
        howmany = Nx             ! all x
        inembed = Ny             ! all y must be used
        onembed = Ny/2+1         ! to all ky
        idist   = 1              ! distance between two slice to transforms (y row)
        odist   = 1
        istride = Nx             ! non contiguous data
        ostride = Nx
#ifdef SINGLE_PRECISION
        CALL sfftw_plan_many_dft_r2c(plan_y2ky_r2c, rank, n, howmany,&
                                     fxy, inembed, istride, idist,&
                                  fxky_g, onembed, ostride, odist,& 
                                    FFTW_FORWARD, FFTW_PATIENT)                
#else
        CALL dfftw_plan_many_dft_r2c(plan_y2ky_r2c, rank, n, howmany,&
                                     fxy, inembed, istride, idist,&
                                  fxky_g, onembed, ostride, odist,& 
                                    FFTW_FORWARD, FFTW_PATIENT)    
#endif
        ! 3.5 MPI splitting along ky (from fxky_g to fxky_l)
        ! 4. (x,ky) -> (kx,ky), C -> C, transforms
        rank    = 1              ! 1D transform
        n       = Nx             ! all x
        howmany = local_nky_ptr  ! only local ky
        inembed = Nx             ! all x must be used
        onembed = local_nky_ptr  ! to the local ky
        idist   = 1              ! distance between two slice to transforms (x row)
        odist   = 1
        istride = Nx             ! non contiguous data
        ostride = Nx
#ifdef SINGLE_PRECISION
        CALL sfftw_plan_many_dft(plan_y2ky_r2c, rank, n, howmany,&
                                fxky_l, inembed, istride, idist,&
                                 fkxky, onembed, ostride, odist,& 
                                FFTW_FORWARD, FFTW_PATIENT)                
#else
        CALL dfftw_plan_many_dft(plan_y2ky_r2c, rank, n, howmany,&
                                fxky_l, inembed, istride, idist,&
                                 fkxky, onembed, ostride, odist,& 
                                FFTW_FORWARD, FFTW_PATIENT)    
#endif
    END SUBROUTINE fft1D_distr_and_plans
    !******************************************************************************!

    !******************************************************************************!
    ! High level routine to ifft a 2D comple array into a real one
    ! It uses variables from the module as the plans
    SUBROUTINE iFFT_2D_c2r
        IMPLICIT NONE
        SELECT CASE (FFT_ALGO)
        CASE ('2D')
        CASE ('1D')
        END SELECT
    END SUBROUTINE iFFT_2D_c2r
    !******************************************************************************!

    !******************************************************************************!
    SUBROUTINE FFT_2D_r2c
        IMPLICIT NONE
        SELECT CASE (FFT_ALGO)
        CASE ('2D')
        CASE ('1D')
        END SELECT
    END SUBROUTINE FFT_2D_r2c  
    !******************************************************************************!

    !******************************************************************************!
    !!! Compute the poisson bracket to real space and sum it to the bracket_sum_r
    !   module variable (convolution theorem)
    SUBROUTINE poisson_bracket_and_sum(ky_, kx_, inv_Ny, inv_Nx, AA_y, AA_x,&
                                        local_nky_ptr, local_nkx_ptr, F_, G_, sum_real_)
        IMPLICIT NONE
        INTEGER(C_INTPTR_T),                  INTENT(IN) :: local_nkx_ptr,local_nky_ptr
        REAL(xp),                             INTENT(IN) :: inv_Nx, inv_Ny
        REAL(xp), DIMENSION(local_nky_ptr),   INTENT(IN) :: ky_, AA_y, AA_x
        REAL(xp), DIMENSION(local_nky_ptr,local_nkx_ptr), INTENT(IN) :: kx_
        COMPLEX(c_xp_c), DIMENSION(local_nky_ptr,local_nkx_ptr) &
                                                         :: F_(:,:), G_(:,:)
        real(c_xp_r), pointer,             INTENT(INOUT) :: sum_real_(:,:)
        INTEGER :: ikx,iky
        !! Anti aliasing
        DO ikx = 1,local_nkx_ptr
            F_(:,ikx) = F_(:,ikx)*AA_y(:)*AA_x(ikx)
            G_(:,ikx) = G_(:,ikx)*AA_y(:)*AA_x(ikx)
        ENDDO
        !------------------------------------------------------------------

        !-------------------- First term df/dx x dg/dy --------------------
        DO ikx = 1,local_nkx_ptr
        DO iky = 1,local_nky_ptr
            cmpx_data_f(ikx,iky) = imagu*kx_(iky,ikx)*F_(iky,ikx)
            cmpx_data_g(ikx,iky) = imagu*ky_(iky)    *G_(iky,ikx)
        ENDDO
        ENDDO

        !CALL iFFT_2D_c2r(cmpx_data_f,real_data_f)

#ifdef SINGLE_PRECISION
        call fftwf_mpi_execute_dft_c2r(planb, cmpx_data_f, real_data_f)
        call fftwf_mpi_execute_dft_c2r(planb, cmpx_data_g, real_data_g)
#else
        call fftw_mpi_execute_dft_c2r(planb, cmpx_data_f, real_data_f)
        call fftw_mpi_execute_dft_c2r(planb, cmpx_data_g, real_data_g)
#endif
        sum_real_ = sum_real_ + real_data_f*real_data_g*inv_Ny*inv_Nx
        !--------------------------------------------------------------------

        !-------------------- Second term -df/dy x dg/dx --------------------
        DO ikx = 1,local_nkx_ptr
        DO iky = 1,local_nky_ptr
            cmpx_data_f(ikx,iky) = imagu*ky_(iky)    *F_(iky,ikx)
            cmpx_data_g(ikx,iky) = imagu*kx_(iky,ikx)*G_(iky,ikx)
        ENDDO
        ENDDO
#ifdef SINGLE_PRECISION
        call fftwf_mpi_execute_dft_c2r(planb, cmpx_data_f, real_data_f)
        call fftwf_mpi_execute_dft_c2r(planb, cmpx_data_g, real_data_g)
#else
        call fftw_mpi_execute_dft_c2r(planb, cmpx_data_f, real_data_f)
        call fftw_mpi_execute_dft_c2r(planb, cmpx_data_g, real_data_g)
#endif
        sum_real_ = sum_real_ - real_data_f*real_data_g*inv_Ny*inv_Nx
    END SUBROUTINE poisson_bracket_and_sum


    SUBROUTINE finalize_plans
        USE basic, ONLY: speak
        IMPLICIT NONE
        CALL speak('..plan Destruction.')

        SELECT CASE (FFT_ALGO)
        CASE ('2D')
            call fftw_destroy_plan(planb)
            call fftw_destroy_plan(planf)
            call fftw_mpi_cleanup()
            call fftw_free(cdatar_f)
            call fftw_free(cdatar_g)
            call fftw_free(cdatar_c)
            call fftw_free(cdatac_f)
            call fftw_free(cdatac_g)
            call fftw_free(cdatac_c)
        CASE ('1D')
        END SELECT
    END SUBROUTINE finalize_plans

END MODULE fourier
