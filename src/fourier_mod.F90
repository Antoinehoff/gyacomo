MODULE fourier
    USE prec_const, ONLY: xp, c_xp_c, c_xp_r, imagu, mpi_xp_c
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
    PUBLIC :: init_grid_distr_and_plans, poisson_bracket_and_sum, finalize_plans, apply_inv_ExB_NL_factor

    !! Module variables
    !! 2D fft specific variables (C interface)
    type(C_PTR)                 :: cdatar_f, cdatar_g, cdatar_c
    type(C_PTR)                 :: cdatac_f, cdatac_g, cdatac_c
    type(C_PTR) ,        PUBLIC :: planf, planb
    integer(C_INTPTR_T), PUBLIC :: alloc_local_1, alloc_local_2
    integer(C_INTPTR_T)         :: NX_, NY_, NY_halved, local_nky_, local_nx_ 
    real   (c_xp_r), pointer, PUBLIC :: real_data_f(:,:), real_data_g(:,:), bracket_sum_r(:,:)
    complex(c_xp_c), pointer, PUBLIC :: cmpx_data_f(:,:), cmpx_data_g(:,:), bracket_sum_c(:,:)
    REAL(xp),                 PUBLIC :: inv_Nx_, inv_Ny_
    !! 1D fft specific variables (full fortran interface)
    type(C_PTR), PUBLIC :: plan_ky2y_c2r ! transform from ( x,ky) to ( x, y) (complex to real)
    type(C_PTR), PUBLIC :: plan_y2ky_r2c ! transform from ( x, y) to ( x,ky) (real to complex)
    type(C_PTR), PUBLIC :: plan_kx2x_c2c ! transform from (kx,ky) to ( x,ky) (complex to complex)
    type(C_PTR), PUBLIC :: plan_x2kx_c2c ! transform from ( x,ky) to (kx,ky) (complex to complex)
    COMPLEX(xp), DIMENSION(:,:), ALLOCATABLE :: f_kxky_l ! working arrays
    COMPLEX(xp), DIMENSION(:,:), ALLOCATABLE :: f_xky_l  ! temp. 1D ifft algo stage (local mpi)
    REAL(xp),    DIMENSION(:,:), ALLOCATABLE :: bracket_sum_yx_g   ! poisson bracket sum in real space
    COMPLEX(xp), DIMENSION(:,:), ALLOCATABLE :: bracket_sum_kyx_g  ! final poisson bracket in complex space

    CONTAINS
    !******************************************************************************!
    !------------- Initialize the grid distribution and plans -------------
    ! If we use the 2D fftw routines, the fftw library decides the best data distribution
    SUBROUTINE init_grid_distr_and_plans(Nx,Ny,communicator,&
            local_nx_ptr,local_nx_ptr_offset,local_nky_ptr,local_nky_ptr_offset)
        USE basic, ONLY: speak
        IMPLICIT NONE
        INTEGER, INTENT(IN)  :: Nx,Ny, communicator
        INTEGER(C_INTPTR_T), INTENT(OUT) :: local_nx_ptr,local_nx_ptr_offset
        INTEGER(C_INTPTR_T), INTENT(OUT) :: local_nky_ptr,local_nky_ptr_offset
        NX_        = Nx; NY_ = Ny
        inv_Nx_    = 1._xp/REAL(NX_,xp)
        NY_halved  = NY_/2 + 1
        inv_Ny_    = 1._xp/REAL(2*NY_halved,xp)     
        ! Call FFTW 2D mpi routines to distribute the data and init 2D MPI FFTW plans
        CALL fft2D_distr_and_plans(communicator,&
                local_nx_ptr,local_nx_ptr_offset,local_nky_ptr,local_nky_ptr_offset)
        local_nx_  = local_nx_ptr  ! store number of local x  in the module
        local_nky_ = local_nky_ptr ! store number of local ky in the module
                ! Init 1D MPI FFTW plans for ExB rate correction factor
        CALL fft1D_plans
        ! store data distr. in the module for the poisson_bracket function
    END SUBROUTINE init_grid_distr_and_plans

    !------------- 2D fft initialization and mpi distribution
    SUBROUTINE fft2D_distr_and_plans(communicator,&
                local_nx_ptr,local_nx_ptr_offset,local_nky_ptr,local_nky_ptr_offset)
        IMPLICIT NONE
        INTEGER, INTENT(IN)  :: communicator
        ! Distribution in the real space along x
        INTEGER(C_INTPTR_T), INTENT(OUT) :: local_nx_ptr, local_nx_ptr_offset
        ! Distribution in the fourier space along ky
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
                local_nx_ptr, local_nx_ptr_offset)
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
        call c_f_pointer(cdatar_f,   real_data_f, [2*(NY_/2 + 1),local_nx_ptr])
        call c_f_pointer(cdatar_g,   real_data_g, [2*(NY_/2 + 1),local_nx_ptr])
        call c_f_pointer(cdatar_c, bracket_sum_r, [2*(NY_/2 + 1),local_nx_ptr])
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
            ERROR STOP '>> ERROR << 2D MPI plan creation error!!'
        end if
    END SUBROUTINE fft2D_distr_and_plans
    !******************************************************************************!

    !******************************************************************************!
    !------------- 1D initialization with balanced data distribution
    SUBROUTINE fft1D_plans
        USE utility,  ONLY: decomp1D
        IMPLICIT NONE
        ! local var
        integer(C_INTPTR_T) :: rank     ! rank of each 1D fourier transforms
        integer(C_INTPTR_T) :: n        ! size of the data to fft
        integer(C_INTPTR_T) :: howmany  ! howmany 1D fourier transforms
        integer(C_INTPTR_T) :: inembed, onembed
        integer(C_INTPTR_T) :: istride, ostride
        integer(C_INTPTR_T) :: idist, odist

        !! Plan of the 4 1D many transforms required
        !----------- 1: FFTx and inv through local ky data
        !1.1 (kx,ky) -> (x,ky), C -> C, transforms
        ! in:
        ALLOCATE( f_kxky_l(NX_,local_nky_))
        ! out:
        ALLOCATE(  f_xky_l(NX_,local_nky_))
        ! transform parameters
        rank    = 1              ! 1D transform
        n       = NX_            ! all kx modes
        howmany = local_nky_     ! all local ky
        inembed = NX_            ! all data must be transformed
        onembed = NX_
        idist   = NX_            ! distance between data to transforms (x columns)
        odist   = NX_
        istride = 1              ! contiguous data
        ostride = 1
#ifdef SINGLE_PRECISION
        CALL sfftw_plan_many_dft(plan_kx2x_c2c, rank, n, howmany,&
                                 f_kxky_l, inembed, istride, idist,&
                                  f_xky_l, onembed, ostride, odist,& 
                                 FFTW_BACKWARD, FFTW_PATIENT)                
#else
        CALL dfftw_plan_many_dft(plan_kx2x_c2c, rank, n, howmany,&
                                 f_kxky_l, inembed, istride, idist,&
                                  f_xky_l, onembed, ostride, odist,& 
                                 FFTW_BACKWARD, FFTW_PATIENT)    
#endif
        ! 1.2 (x,ky) -> (kx,ky), C -> C, transforms
        ! in:  f_xky_l
        ! out: f_kxky_l
        ! transform parameters
        rank    = 1              ! 1D transform
        n       = NX_            ! all kx modes
        howmany = local_nky_     ! all local ky
        inembed = NX_            ! all data must be transformed
        onembed = NX_
        idist   = NX_            ! distance between data to transforms (x columns)
        odist   = NX_
        istride = 1              ! contiguous data
        ostride = 1
#ifdef SINGLE_PRECISION
        CALL sfftw_plan_many_dft(plan_x2kx_c2c, rank, n, howmany,&
                                 f_xky_l, inembed, istride, idist,&
                                f_kxky_l, onembed, ostride, odist,& 
                                FFTW_FORWARD, FFTW_PATIENT)                
#else
        CALL dfftw_plan_many_dft(plan_x2kx_c2c, rank, n, howmany,&
                                 f_xky_l, inembed, istride, idist,&
                                f_kxky_l, onembed, ostride, odist,& 
                                FFTW_FORWARD, FFTW_PATIENT)    
#endif

        !----------- 2: FFTy and inv through global ky data
        ! 2.1 (y,x) -> (ky,x), R -> C, transforms (bplan_y in GENE)
        ! in:
        ALLOCATE(bracket_sum_yx_g(2*NY_halved,local_nx_))
        ! out:
        ALLOCATE(bracket_sum_kyx_g(NY_halved+1,local_nx_))
        ! transform parameters
        rank    = 1              ! 1D transform
        n       = 2*NY_halved    ! all y
        howmany = local_nx_      ! all x
        inembed = 2*NY_halved    ! all y must be used
        istride = 1              ! contiguous data
        idist   = 2*NY_halved    ! distance between two slice to transforms (y row)
        onembed = NY_halved+1
        ostride = 1
        odist   = NY_halved+1
#ifdef SINGLE_PRECISION
        CALL sfftw_plan_many_dft_r2c(plan_y2ky_r2c, rank, n, howmany,&
                             bracket_sum_yx_g, inembed, istride, idist,&
                            bracket_sum_kyx_g, onembed, ostride, odist,& 
                                    FFTW_PATIENT)                
#else
        CALL dfftw_plan_many_dft_r2c(plan_y2ky_r2c, rank, n, howmany,&
                             bracket_sum_yx_g, inembed, istride, idist,&
                            bracket_sum_kyx_g, onembed, ostride, odist,& 
                                    FFTW_PATIENT)    
#endif
        !-----------
        ! 2.2 (x,ky) -> (x,y), C -> R, transforms (fplan_y in GENE)
        ! in:  bracket_sum_kyx_g
        ! out: bracket_sum_yx_g
        ! transform parameters
        rank    = 1               ! 1D transform
        n       = 2*NY_halved     ! all ky
        howmany = local_nx_       ! all x
        inembed = NY_halved+1     ! all ky must be used
        istride = 1               ! contiguous data
        idist   = NY_halved+1     ! distance between two slice to transforms (y row)
        onembed = 2*NY_halved     ! to all y
        ostride = 1
        odist   = 2*NY_halved
#ifdef SINGLE_PRECISION
        CALL sfftw_plan_many_dft_c2r(plan_ky2y_c2r, rank, n, howmany,&
                            bracket_sum_kyx_g, inembed, istride, idist,&
                             bracket_sum_yx_g, onembed, ostride, odist,& 
                                    FFTW_PATIENT)                
#else
        CALL dfftw_plan_many_dft_c2r(plan_ky2y_c2r, rank, n, howmany,&
                            bracket_sum_kyx_g, inembed, istride, idist,&
                             bracket_sum_yx_g, onembed, ostride, odist,& 
                                    FFTW_PATIENT)    
#endif
    if ((.not. c_associated(plan_ky2y_c2r)) .OR. &
        (.not. c_associated(plan_y2ky_r2c)) .OR. &
        (.not. c_associated(plan_kx2x_c2c)) .OR. &
        (.not. c_associated(plan_x2kx_c2c))) then
        ERROR STOP '>> ERROR << 1D plan creation error!!'
    end if
    ! Free mem (optional)
    DEALLOCATE(f_xky_l,f_kxky_l)
    DEALLOCATE(bracket_sum_kyx_g,bracket_sum_yx_g)
END SUBROUTINE fft1D_plans
    !******************************************************************************!

    !******************************************************************************!
    !!! Compute the poisson bracket to real space and sum it to the bracket_sum_r
    !   module variable (convolution theorem)
    SUBROUTINE poisson_bracket_and_sum( ky_array, kx_array, inv_Ny, inv_Nx, AA_y, AA_x,&
                                        local_nky, total_nkx, F_, G_,&
                                        ExB, ExB_NL_factor, sum_real_)
        IMPLICIT NONE
        INTEGER,                                  INTENT(IN) :: local_nky,total_nkx
        REAL(xp),                                 INTENT(IN) :: inv_Nx, inv_Ny
        REAL(xp), DIMENSION(local_nky),           INTENT(IN) :: ky_array, AA_y
        REAL(xp), DIMENSION(total_nkx),           INTENT(IN) :: AA_x
        REAL(xp), DIMENSION(local_nky,total_nkx), INTENT(IN) :: kx_array
        COMPLEX(c_xp_c), DIMENSION(local_nky,total_nkx), &
                                            INTENT(IN)      :: F_, G_
        COMPLEX(xp), DIMENSION(total_nkx,local_nky), &
                                            INTENT(IN)      :: ExB_NL_factor
        LOGICAL, INTENT(IN) :: ExB
        real(c_xp_r), pointer,              INTENT(INOUT)   :: sum_real_(:,:)
        ! local variables
        INTEGER :: ikx,iky
        COMPLEX(xp), DIMENSION(total_nkx,local_nky) :: ikxF, ikyG, ikyF, ikxG
        REAL(xp):: kxs, ky
        
        ! Build the fields to convolve
        ! Store df/dx, dg/dy and df/dy, dg/dx
        DO iky = 1,local_nky
                DO ikx = 1,total_nkx
                        ky  = ky_array(iky)
                        kxs = kx_array(iky,ikx)
                        ikxF(ikx,iky) = imagu*kxs*F_(iky,ikx)
                        ikyG(ikx,iky) = imagu*ky *G_(iky,ikx)
                        ikyF(ikx,iky) = imagu*ky *F_(iky,ikx)
                        ikxG(ikx,iky) = imagu*kxs*G_(iky,ikx)
                ENDDO
        ENDDO
        IF(ExB) THEN 
            ! Apply the ExB shear correction factor exp(ixkySJdT)
            CALL apply_ExB_NL_factor(ikxF,ExB_NL_factor)
            CALL apply_ExB_NL_factor(ikyG,ExB_NL_factor)
            CALL apply_ExB_NL_factor(ikyF,ExB_NL_factor)
            CALL apply_ExB_NL_factor(ikxG,ExB_NL_factor)
        ENDIF
        ! Anti Aliasing
        DO iky = 1,local_nky
                ikxF(:,iky) = ikxF(:,iky)*AA_y(iky)*AA_x(:)
                ikyG(:,iky) = ikyG(:,iky)*AA_y(iky)*AA_x(:)
                ikyF(:,iky) = ikyF(:,iky)*AA_y(iky)*AA_x(:)
                ikxG(:,iky) = ikxG(:,iky)*AA_y(iky)*AA_x(:)
        ENDDO
        !-------------------- First term df/dx x dg/dy --------------------
#ifdef SINGLE_PRECISION
        call fftwf_mpi_execute_dft_c2r(planb, ikxF, real_data_f)
        call fftwf_mpi_execute_dft_c2r(planb, ikyG, real_data_g)
#else
        call fftw_mpi_execute_dft_c2r(planb, ikxF, real_data_f)
        call fftw_mpi_execute_dft_c2r(planb, ikyG, real_data_g)
#endif
        sum_real_ = sum_real_ + real_data_f*real_data_g*inv_Ny*inv_Nx
        !--------------------------------------------------------------------

        !-------------------- Second term -df/dy x dg/dx --------------------
#ifdef SINGLE_PRECISION
        call fftwf_mpi_execute_dft_c2r(planb, ikyF, real_data_f)
        call fftwf_mpi_execute_dft_c2r(planb, ikxG, real_data_g)
#else
        call fftw_mpi_execute_dft_c2r(planb, ikyF, real_data_f)
        call fftw_mpi_execute_dft_c2r(planb, ikxG, real_data_g)
#endif
        sum_real_ = sum_real_ - real_data_f*real_data_g*inv_Ny*inv_Nx
    END SUBROUTINE poisson_bracket_and_sum
    !******************************************************************************!

    !******************************************************************************!
    ! Apply the exp(-i*x*ky*S*J*dt) factor to the Poisson bracket fields (ikxF, ikyG, ...) 
    ! (see Mcmillan et al. 2019)
    SUBROUTINE apply_ExB_NL_factor(fkxky,ExB_NL_factor)
        IMPLICIT NONE
        COMPLEX(xp), DIMENSION(NX_,local_nky_), INTENT(INOUT)  :: fkxky
        COMPLEX(xp), DIMENSION(NX_,local_nky_), INTENT(IN)     :: ExB_NL_factor
        ! local variables
        COMPLEX(xp), DIMENSION(NX_,local_nky_) :: tmp_kxky, tmp_xky
        integer(C_INTPTR_T) :: ix,ikx,iky
        ! Fill the buffer
        DO iky = 1,local_nky_
                DO ikx = 1,NX_
                        tmp_kxky(ikx,iky) = fkxky(ikx,iky)
                ENDDO
        ENDDO
        ! go to x,ky
        CALL iFFT_kxky_to_xky(tmp_kxky,tmp_xky)
        ! multiply result by NL factor
        DO iky = 1,local_nky_
            DO ix = 1,NX_
                tmp_xky(ix,iky) = inv_Nx_*tmp_xky(ix,iky)*ExB_NL_factor(ix,iky)
            ENDDO
        ENDDO
        ! back to kx,ky
        CALL FFT_xky_to_kxky(tmp_xky,tmp_kxky)
        ! fill output array
        DO iky = 1,local_nky_
                DO ikx = 1,NX_
                        fkxky(ikx,iky) = tmp_kxky(ikx,iky)
                ENDDO
        ENDDO
    END SUBROUTINE apply_ExB_NL_factor

    ! Apply the exp(i*x*ky*S*J*dt) factor to the real Poisson Bracket sum [f,g]
    SUBROUTINE apply_inv_ExB_NL_factor(fyx,inv_ExB_NL_factor)
        IMPLICIT NONE
        real(c_xp_r), pointer,                        INTENT(INOUT) :: fyx(:,:)
        COMPLEX(xp), DIMENSION(NY_halved,local_nx_),  INTENT(IN)    :: inv_ExB_NL_factor
        ! local variables
        REAL(xp),    DIMENSION(2*NY_halved,local_nx_) :: tmp_yx_1, tmp_yx_2
        COMPLEX(xp), DIMENSION(NY_halved+1,local_nx_) :: tmp_kyx
        integer(C_INTPTR_T) :: ix, iy, iky
        ! Fill buffer
        DO ix = 1,local_nx_
            DO iy = 1,2*NY_halved
                tmp_yx_1(iy,ix) = fyx(iy,ix)
            ENDDO
        ENDDO
        ! Fourier real to complex on the second buffer (first buffer is now unusable)
        CALL FFT_yx_to_kyx(tmp_yx_1,tmp_kyx)
        ! Treat the result with the ExB NL factor
        DO iky = 1,NY_halved
                DO ix = 1,local_nx_
                        tmp_kyx(iky,ix) = tmp_kyx(iky,ix)*inv_ExB_NL_factor(iky,ix)
                ENDDO
        ENDDO
        ! Back to Fourier space in the third buffer
        CALL iFFT_kyx_to_yx(tmp_kyx,tmp_yx_2)
        ! Fill the result array with normalization of iFFT
        DO ix = 1,local_nx_
            DO iy = 1,2*NY_halved
                fyx(iy,ix) = tmp_yx_2(iy,ix)*inv_Ny_
            ENDDO
        ENDDO
    END SUBROUTINE apply_inv_ExB_NL_factor

    !******************************************************************************!
    ! High level FFT routines
    ! 1D inverse Fourier transform from (kx,ky) to (x,ky), complex to complex
    SUBROUTINE iFFT_kxky_to_xky(in_kxky,out_xky)
        IMPLICIT NONE
        COMPLEX(xp), DIMENSION(NX_,local_nky_), INTENT(IN)  :: in_kxky
        COMPLEX(xp), DIMENSION(NX_,local_nky_), INTENT(OUT) :: out_xky
#ifdef SINGLE_PRECISION
        CALL sfftw_execute_dft(plan_kx2x_c2c, in_kxky, out_xky)
#else 
        CALL dfftw_execute_dft(plan_kx2x_c2c, in_kxky, out_xky)
#endif
    END SUBROUTINE iFFT_kxky_to_xky
    ! 1D Fourier transform from kx,ky) to (kx,ky), complex to complex
    SUBROUTINE FFT_xky_to_kxky(in_xky,out_kxky)
        IMPLICIT NONE
        COMPLEX(xp), DIMENSION(NX_,local_nky_), INTENT(IN)  :: in_xky
        COMPLEX(xp), DIMENSION(NX_,local_nky_), INTENT(OUT) :: out_kxky
#ifdef SINGLE_PRECISION
        CALL sfftw_execute_dft(plan_x2kx_c2c, in_xky, out_kxky)
#else 
        CALL dfftw_execute_dft(plan_x2kx_c2c, in_xky, out_kxky)
#endif
    END SUBROUTINE FFT_xky_to_kxky
    ! 1D Fourier transform from (y,x) to (ky,x), real to complex
    SUBROUTINE FFT_yx_to_kyx(in_yx,out_kyx)
        IMPLICIT NONE
        REAL(xp),    DIMENSION(2*NY_halved,local_nx_), INTENT(IN)  :: in_yx
        COMPLEX(xp), DIMENSION(NY_halved+1,local_nx_), INTENT(OUT) :: out_kyx
#ifdef SINGLE_PRECISION
        CALL sfftw_execute_dft_r2c(plan_y2ky_r2c, in_yx, out_kyx)
#else 
        CALL dfftw_execute_dft_r2c(plan_y2ky_r2c, in_yx, out_kyx)
#endif
    END SUBROUTINE FFT_yx_to_kyx
    ! 1D Fourier transform from (ky,x) to (y,x), complex to real
    SUBROUTINE iFFT_kyx_to_yx(in_kyx,out_yx)
        IMPLICIT NONE
        COMPLEX(xp), DIMENSION(NY_halved+1,local_nx_), INTENT(IN)  :: in_kyx
        REAL(xp),    DIMENSION(2*NY_halved,local_nx_), INTENT(OUT) :: out_yx
#ifdef SINGLE_PRECISION
        CALL sfftw_execute_dft_c2r(plan_ky2y_c2r, in_kyx, out_yx)
#else 
        CALL dfftw_execute_dft_c2r(plan_ky2y_c2r, in_kyx, out_yx)
#endif
    END SUBROUTINE iFFT_kyx_to_yx

    !******************************************************************************!

    SUBROUTINE finalize_plans
        USE basic, ONLY: speak
        IMPLICIT NONE
        CALL speak('..plan Destruction.')
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
