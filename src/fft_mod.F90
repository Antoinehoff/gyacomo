module fft
  ! Fast fourier transform module
  ! Used to supress high frequency instabilities

  use prec_const
  use array
  implicit none
  private


  ! Public Functions
  PUBLIC :: fft_recursive, initialize_filter, suppress_high_freq


  contains
 
    ! In place Cooley-Tukey FFT, tsource : rosettacode.org/wiki/Fast_Fourier_transform
    recursive subroutine fft_recursive(x)
      complex(kind=dp), dimension(:), intent(inout)  :: x
      complex(kind=dp)                               :: t
      integer                                        :: N
      integer                                        :: i
      complex(kind=dp), dimension(:), allocatable    :: even, odd
   
      N=size(x)
   
      if(N .le. 1) return
   
      allocate(odd((N+1)/2))
      allocate(even(N/2))
   
      ! divide
      odd =x(1:N:2)
      even=x(2:N:2)
   
      ! conquer
      call fft_recursive(odd)
      call fft_recursive(even)
   
      ! combine
      do i=1,N/2
         t=exp(cmplx(0.0_dp,-2.0_dp*pi*real(i-1,dp)/real(N,dp),kind=dp))*even(i)
         x(i)     = odd(i) + t
         x(i+N/2) = odd(i) - t
      end do
   
      deallocate(odd)
      deallocate(even)
    end subroutine fft_recursive


    subroutine initialize_filter
        use space_grid
        use array
        use model
        implicit none
        integer :: i
        integer :: N
        integer :: center

        N = (ize-izs+1)
        center = int( ceiling( (N+1)*0.5 ) )

        do i=izs,ize
           fft_filter(i) = 1._dp- exp( -1._dp*(i-center)*(i-center)/ (2*fft_sigma*fft_sigma) ) ! or just a step function?
           ! write(*,'("(", F20.15, ",", F20.15, "i )")') fft_filter(i)
        end do
    end subroutine initialize_filter


    subroutine suppress_high_freq(x)
        use space_grid
        use array
        implicit none
        integer :: i
        real(dp), dimension(izs:ize), intent(inout) :: x
        real(dp) :: invN
        
        do i=izs,ize
          fft_input(i) = x(i)
        end do

        call fft_recursive(fft_input) ! Fast Fourier Transform

        do i=izs,ize ! Apply filter
          fft_input(i) = fft_input(i)*fft_filter(i)
        end do

        call fft_recursive(fft_input) ! Invert the FFT

        invN = 1._dp/(ize-izs+1._dp)
        do i=izs,ize
          x(i) = real(fft_input(ize-i+1))*invN 
        end do
    end subroutine suppress_high_freq


end module fft
 
