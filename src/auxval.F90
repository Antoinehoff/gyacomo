subroutine auxval
  !   Set auxiliary values, at beginning of simulation

  USE basic
  USE space_grid
  USE array
  USE model, ONLY: gradient_scheme, diff_theta, diff_temp, diff_vpar, diff_moments
  USE gradients
  ! use mumps_bsplines, only: putrow_mumps => putrow, updtmat_mumps => updtmat, factor_mumps => factor, bsolve_mumps => bsolve ! from BSPLINES
  use fft
  use prec_const
  IMPLICIT NONE
  
  INTEGER :: irows,irowe, irow, icol
  real(dp), DIMENSION(-2:2) :: stencil ! For laplace operator (needed for poisson solver)

  WRITE(*,'(a/)') '=== Set auxiliary values ==='

  call set_pgrid
  call set_zgrid


  neq = (ize-izs+1) + 1 ! Number of equations to solve for poisson = Nz + 1 (additional equation to fix the constant)

  CALL memory ! Allocate memory for global arrays

  CALL set_gradient_scheme(gradient_scheme) ! Point to correct gradient functions
  if ((gradient_scheme .eq. 'up2') .and. ((diff_theta .ne. 0._dp) .or. (diff_temp .ne. 0._dp) .or. (diff_vpar .ne. 0._dp) .or. (diff_moments .ne. 0._dp))) then
    write(*,*) 'ERROR: the scheme is upwind but the numerical diffusion coefficients are non-zero. Program ends.'
    nlend = .true.
  endif


! Setting up laplace matrix, eg (-2,1,0,0; 1,-2,1,0; 0,1,-2,1; 0,0,1,-2), to solve poisson equation

!  call set_laplace_stencil(gradient_scheme, stencil) ! for fa2 it is (0,1,-2,1,0)/dx**2 
  ! call set_gradient_stencil(gradient_scheme, stencil) ! for fa2 it is (0,-1/2,0,1/2,0)/dx 

!  irows = 1
!  irowe = neq-1

!  do irow=irows,irowe
!    do icol=irows,irowe
      ! do isten=max(-2,irows-irow),min(2,irowe-irow) ! loop over stencil
!      if ( (icol-irow .le. 2) .and. (icol-irow .ge. -2) ) then
!        call putele(laplace_smumps, irow, icol, stencil(icol-irow))
!      else
!        call putele(laplace_smumps, irow, icol, 0._dp)
!      endif
!    enddo
!  enddo

  ! Periodic boundary conditions. Eg transform last row from (0,0,0,1,-2) to (1,0,0,1,-2)
!  call putele(laplace_smumps, irows, irowe-1, stencil(-2))
!  call putele(laplace_smumps, irows, irowe, stencil(-1))

!  call putele(laplace_smumps, irows+1, irowe, stencil(-2))

!  call putele(laplace_smumps, irowe-1, irows, stencil(2))

!  call putele(laplace_smumps, irowe, irows, stencil(1))
!  call putele(laplace_smumps, irowe, irows+1, stencil(2))


  ! ! Boundary condition to fix the value of elec at position iz=21
  ! do icol=irows,irowe
  !   call putele(laplace_smumps, irows+20, icol, 0._dp)
  ! enddo
  ! call putele(laplace_smumps, irows+20, irows+20, 1._dp)



  ! ! Boundary condition to fix the mean value of electric field (replaces equation for iz=21)
  ! do icol=irows,irowe
  !   call putele(laplace_smumps, irows+20, icol, 1._dp/real(ize-izs+1,dp) )
  ! enddo


  ! Lagrange multiplier, adds an equation at neq^th row to fix the constant of integration by fixing the average value
!  do icol=irows,irowe
!    call putele(laplace_smumps, icol, neq, 1._dp)
!    call putele(laplace_smumps, neq, icol, deltaz)
!  enddo
!  call putele(laplace_smumps, neq, neq, 0._dp)




!  call factor(laplace_smumps) ! Factor



  CALL initialize_filter 


 

CONTAINS

  subroutine set_laplace_stencil(gradient_scheme, laplace_stencil)

    character(3), intent(in) :: gradient_scheme
    real(dp), dimension(-2:2), intent(out) :: laplace_stencil
    
    laplace_stencil = 0.0_dp

    select case (gradient_scheme)
    ! case('fa2','up2')
    !    laplace_stencil(-1) = 1._dp * deltazisq
    !    laplace_stencil(0) = -2._dp * deltazisq
    !    laplace_stencil(1) =  1._dp * deltazisq
       
    case('fd4','we4')
       laplace_stencil(-2) = -deltazisq/12._dp 
       laplace_stencil(-1) = deltazisq*16._dp/12._dp
       laplace_stencil(0) = -deltazisq*30._dp/12._dp
       laplace_stencil(1) = deltazisq*16._dp/12._dp
       laplace_stencil(2) = -deltazisq/12._dp 

    end select

  end subroutine set_laplace_stencil

END SUBROUTINE auxval
