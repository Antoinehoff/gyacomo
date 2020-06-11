! see "Efficient Implementation of Weighted ENO Schemes", GUANG-SHAN JIANG et all, 1996

SUBROUTINE gradz_we4( f , f_z ) ! implemented for periodic boundary conditions

  use prec_const
  use space_grid
  implicit none
  real(dp), dimension(izs:ize),intent(in) :: f
  real(dp), dimension(izs:ize),intent(out) :: f_z
  real(dp) :: tmp
  integer :: i

  call compute_f_hat_plus_one_half_we4( f, f_z )
  ! f_z contains f_hat_plus_one_half

  tmp = f_z(ize)

  do i=1,ize-izs
    f_z(ize+1-i) = ( f_z(ize+1-i) - f_z(ize-i) )*deltazi
  enddo
  f_z(izs) =  ( f_z(izs)- tmp       )*deltazi

!!  write(*,*) "gradz_we4"

END SUBROUTINE gradz_we4



SUBROUTINE compute_IS_we4( f, IS ) 

  use prec_const
  use space_grid
  use array
  implicit none
  real(dp), dimension(izs:ize),intent(in) :: f
  real(dp), dimension(izs:ize,0:2),intent(out) :: IS
  integer :: iz


  ! SI_0
  tmp1(izs)   = f(ize-1) - 2.0_dp*f(ize)  + f(izs)
  tmp1(izs+1) = f(ize)   - 2.0_dp*f(izs)  + f(izs+1)
  do iz=izs+2,ize
    tmp1(iz)  = f(iz-2)  - 2.0_dp*f(iz-1) + f(iz)
  enddo

  tmp2(izs)   = f(ize-1) - 4.0_dp*f(ize)  + 3.0_dp*f(izs)
  tmp2(izs+1) = f(ize)   - 4.0_dp*f(izs)  + 3.0_dp*f(izs+1)
  do iz=izs+2,ize
    tmp2(iz)  = f(iz-2)  - 4.0_dp*f(iz-1) + 3.0_dp*f(iz)
  enddo

  do iz=izs,ize
    IS(iz,0)  = thirteentwelfths*tmp1(iz)*tmp1(iz) + 0.25_dp*tmp2(iz)*tmp2(iz)
  enddo


  ! SI_1
  tmp1(izs)   = f(ize)   - 2.0_dp*f(izs)  + f(izs+1)
  do iz=izs+1,ize-1
    tmp1(iz)  = f(iz-1)  - 2.0_dp*f(iz)   + f(iz+1)
  enddo
  tmp1(ize)   = f(ize-1) - 2.0_dp*f(ize)  + f(izs)

  tmp2(izs)   = f(ize)  - f(izs+1)
  do iz=izs+1,ize-1
    tmp2(iz)  = f(iz-1) - f(iz+1)
  enddo
  tmp2(ize)   = f(ize-1)- f(izs)

  do iz=izs,ize
    IS(iz,1)  = thirteentwelfths*tmp1(iz)*tmp1(iz) + 0.25_dp*tmp2(iz)*tmp2(iz)
  enddo


  ! SI_2
  do iz=izs,ize-2
    tmp1(iz)  = f(iz)    - 2.0_dp*f(iz+1) + f(iz+2)
  enddo
  tmp1(ize-1) = f(ize-1) - 2.0_dp*f(ize)  + f(izs)
  tmp1(ize)   = f(ize)   - 2.0_dp*f(izs)  + f(izs+1)

  do iz=izs,ize-2
    tmp2(iz)  = 3.0_dp*f(iz)   - 4.0_dp*f(iz+1) + f(iz+2)
  enddo
  tmp2(ize-1) = 3.0_dp*f(ize-1) - 4.0_dp*f(ize)  + f(izs)
  tmp2(ize)   = 3.0_dp*f(ize)   - 4.0_dp*f(izs)  + f(izs+1)

  do iz=izs,ize
    IS(iz,2)  = thirteentwelfths*tmp1(iz)*tmp1(iz) + 0.25_dp*tmp2(iz)*tmp2(iz)
  enddo

END SUBROUTINE compute_IS_we4


SUBROUTINE compute_omegas_we4( f, omegas ) 

  use prec_const
  use space_grid
  implicit none
  real(dp), dimension(izs:ize),intent(in) :: f
  real(dp), dimension(izs:ize,0:2),intent(out) :: omegas
  real(dp) :: epsi = 0.000001_dp ! 1e-6
  ! real(dp) :: p = 2.0_dp ! this is hardcoded to go faster
  real(dp), dimension(0:2) :: C3k = (/ 0.1, 0.6, 0.3 /)
  integer :: iz
  integer :: k
  real(dp) :: sumalphas  ! tmp variable
  

  call compute_IS_we4( f, omegas ) 
  ! omegas now contains IS


  do k=0,2
    do iz=izs,ize
      omegas(iz,k) = C3k(k) / ( (epsi + omegas(iz,k))*(epsi + omegas(iz,k)) )
!      omegas(iz,k) = C3k(k) / power( epsi + IS, p)
    enddo
  enddo
  ! omegas now contains alphas


  do iz=izs,ize
    sumalphas = omegas(iz,0)+omegas(iz,1)+omegas(iz,2)

    do k=0,2
      omegas(iz,k) = omegas(iz,k) / sumalphas
    enddo
  enddo
  ! omegas now contains omegas

END SUBROUTINE compute_omegas_we4




SUBROUTINE compute_f_hat_plus_one_half_we4( f, f_hat ) 

  use prec_const
  use space_grid
  use array
  implicit none
  real(dp), dimension(izs:ize),intent(in) :: f
  real(dp), dimension(izs:ize),intent(out) :: f_hat
  real(dp) :: epsi = 0.000001_dp ! 1e-6
  ! real(dp) :: p = 2.0_dp ! this is hardcoded to go faster
  real(dp), dimension(0:2,0:2) :: a3kl = reshape((/ onethird, -sevensixths, elevensixths, &
                                                   -onesixth, fivesixths, onethird, &
                                                    onethird,  fivesixths,  -onesixth /) &
                                                  , shape(a3kl), order=(/2,1/) )
  integer :: iz


  call compute_omegas_we4( f, omegak_times_q3k ) 
  ! omegak_times_q3k contains omegas


  ! k=0
  omegak_times_q3k(izs,0)  = omegak_times_q3k(izs,0)  *( a3kl(0,0)*f(ize-1)+ a3kl(0,1)*f(ize)  + a3kl(0,2)*f(izs) )
  omegak_times_q3k(izs+1,0)= omegak_times_q3k(izs+1,0)*( a3kl(0,0)*f(ize)  + a3kl(0,1)*f(izs)  + a3kl(0,2)*f(izs+1) )
  do iz=izs+2,ize
    omegak_times_q3k(iz,0) = omegak_times_q3k(iz,0)   *( a3kl(0,0)*f(iz-2) + a3kl(0,1)*f(iz-1) + a3kl(0,2)*f(iz) )
  enddo


  ! k=1
  omegak_times_q3k(izs,1)  = omegak_times_q3k(izs,1)  *( a3kl(1,0)*f(ize)  + a3kl(1,1)*f(izs)  + a3kl(1,2)*f(izs+1) )
  do iz=izs+1,ize-1
    omegak_times_q3k(iz,1) = omegak_times_q3k(iz,1)   *( a3kl(1,0)*f(iz-1) + a3kl(1,1)*f(iz)   + a3kl(1,2)*f(iz+1) )
  enddo
  omegak_times_q3k(ize,1)  = omegak_times_q3k(ize,1)  *( a3kl(1,0)*f(ize-1)+ a3kl(1,1)*f(ize)  + a3kl(1,2)*f(izs) )

  ! k=2
  do iz=izs,ize-2
    omegak_times_q3k(iz,2) = omegak_times_q3k(iz,2)   *( a3kl(2,0)*f(iz)   + a3kl(2,1)*f(iz+1) + a3kl(2,2)*f(iz+2) )
  enddo
  omegak_times_q3k(ize-1,2)= omegak_times_q3k(ize-1,2)*( a3kl(2,0)*f(ize-1)+ a3kl(2,1)*f(ize)  + a3kl(2,2)*f(izs) )
  omegak_times_q3k(ize,2)  = omegak_times_q3k(ize,2)  *( a3kl(2,0)*f(ize)  + a3kl(2,1)*f(izs)  + a3kl(2,2)*f(izs+1) )



  do iz=izs,ize
    f_hat(iz) = omegak_times_q3k(iz,0) + omegak_times_q3k(iz,1) + omegak_times_q3k(iz,2)
  enddo    

END SUBROUTINE compute_f_hat_plus_one_half_we4




