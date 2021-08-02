MODULE numerics
    USE basic
    USE prec_const
    USE grid
    USE utility
    USE coeff
    implicit none

    PUBLIC :: compute_derivatives, build_dnjs_table, evaluate_kernels, compute_lin_coeff

CONTAINS

! Compute the 2D particle temperature for electron and ions (sum over Laguerre)
SUBROUTINE compute_derivatives

END SUBROUTINE compute_derivatives

!******************************************************************************!
!!!!!!! Build the Laguerre-Laguerre coupling coefficient table for nonlin
!******************************************************************************!
SUBROUTINE build_dnjs_table
  USE array, Only : dnjs
  USE coeff
  IMPLICIT NONE

  INTEGER :: in, ij, is, J
  INTEGER :: n_, j_, s_

  J = max(jmaxe,jmaxi)

  DO in = 1,J+1 ! Nested dependent loops to make benefit from dnjs symmetry
    n_ = in - 1
    DO ij = in,J+1
      j_ = ij - 1
      DO is = ij,J+1
        s_ = is - 1

        dnjs(in,ij,is) = TO_DP(ALL2L(n_,j_,s_,0))
        ! By symmetry
        dnjs(in,is,ij) = dnjs(in,ij,is)
        dnjs(ij,in,is) = dnjs(in,ij,is)
        dnjs(ij,is,in) = dnjs(in,ij,is)
        dnjs(is,ij,in) = dnjs(in,ij,is)
        dnjs(is,in,ij) = dnjs(in,ij,is)
      ENDDO
    ENDDO
  ENDDO
END SUBROUTINE build_dnjs_table
!******************************************************************************!

!******************************************************************************!
!!!!!!! Evaluate the kernels once for all
!******************************************************************************!
SUBROUTINE evaluate_kernels
  USE basic
  USE array, Only : kernel_e, kernel_i
  USE grid
  use model, ONLY : tau_e, tau_i, sigma_e, sigma_i, q_e, q_i, lambdaD, CLOS, sigmae2_taue_o2, sigmai2_taui_o2
  IMPLICIT NONE

  REAL(dp)    :: factj, j_dp, j_int
  REAL(dp)    :: be_2, bi_2, alphaD
  REAL(dp)    :: kx, ky, kperp2

  !!!!! Electron kernels !!!!!
  !Precompute species dependant factors
  factj = 1.0 ! Start of the recursive factorial
  DO ij = 1, jmaxe+1
    j_int = jarray_e(ij)
    j_dp = REAL(j_int,dp) ! REAL of degree

    ! Recursive factorial
    IF (j_dp .GT. 0) THEN
      factj = factj * j_dp
    ELSE
      factj = 1._dp
    ENDIF

    DO ikx = ikxs,ikxe
      kx     = kxarray(ikx)
      DO iky = ikys,ikye
        ky    = kyarray(iky)

        be_2  =  (kx**2 + ky**2) * sigmae2_taue_o2

        kernel_e(ij, ikx, iky) = be_2**j_int * exp(-be_2)/factj

      ENDDO
    ENDDO
  ENDDO
  ! Kernels closure
  DO ikx = ikxs,ikxe
    kx     = kxarray(ikx)
    DO iky = ikys,ikye
      ky    = kyarray(iky)
      be_2  =  (kx**2 + ky**2) * sigmae2_taue_o2
      ! Kernel ghost + 1 with Kj+1 = y/(j+1) Kj (/!\ ij = j+1)
      kernel_e(ijeg_e,ikx,iky) = be_2/(real(ijeg_e-1,dp))*kernel_e(ije_e,ikx,iky)
      ! Kernel ghost - 1 with Kj-1 = j/y Kj(careful about the kperp=0)
      IF ( be_2 .NE. 0 ) THEN
        kernel_e(ijsg_e,ikx,iky) = (real(ijsg_e,dp))/be_2*kernel_e(ijs_e,ikx,iky)
      ELSE
        kernel_e(ijsg_e,ikx,iky) = 0._dp
      ENDIF
    ENDDO
  ENDDO

  !!!!! Ion kernels !!!!!
  factj = 1.0 ! Start of the recursive factorial
  DO ij = 1, jmaxi+1
    j_int = jarray_e(ij)
    j_dp = REAL(j_int,dp) ! REAL of degree

    ! Recursive factorial
    IF (j_dp .GT. 0) THEN
      factj = factj * j_dp
    ELSE
      factj = 1._dp
    ENDIF

    DO ikx = ikxs,ikxe
      kx     = kxarray(ikx)
      DO iky = ikys,ikye
        ky    = kyarray(iky)
        bi_2  =  (kx**2 + ky**2) * sigmai2_taui_o2
        kernel_i(ij, ikx, iky) = bi_2**j_int * exp(-bi_2)/factj

      ENDDO
    ENDDO
  ENDDO
  ! Kernels closure
  DO ikx = ikxs,ikxe
    kx     = kxarray(ikx)
    DO iky = ikys,ikye
      ky    = kyarray(iky)
      bi_2  =  (kx**2 + ky**2) * sigmai2_taui_o2
      ! Kernel ghost + 1 with Kj+1 = y/(j+1) Kj
      kernel_i(ijeg_i,ikx,iky) = bi_2/(real(ijeg_i-1,dp))*kernel_i(ije_i,ikx,iky)
      ! Kernel ghost - 1 with Kj-1 = j/y Kj(careful about the kperp=0)
      IF ( bi_2 .NE. 0 ) THEN
        kernel_i(ijsg_i,ikx,iky) = (real(ijsg_i,dp))/bi_2*kernel_i(ijs_i,ikx,iky)
      ELSE
        kernel_i(ijsg_i,ikx,iky) = 0._dp
      ENDIF
    ENDDO
  ENDDO
END SUBROUTINE evaluate_kernels
!******************************************************************************!

SUBROUTINE compute_lin_coeff
  USE array, ONLY: xnepj, xnepp1j, xnepm1j, xnepp2j, xnepm2j, xnepjp1, xnepjm1, &
                   xnipj, xnipp1j, xnipm1j, xnipp2j, xnipm2j, xnipjp1, xnipjm1, &
                   xphij, xphijp1, xphijm1
  USE model, ONLY: taue_qe, taui_qi, sqrtTaue_qe, sqrtTaui_qi, eta_T, eta_n
  USE prec_const
  USE grid,  ONLY: parray_e, parray_i, jarray_e, jarray_i, &
                   ip,ij, ips_e,ipe_e, ips_i,ipe_i, ijs_e,ije_e, ijs_i,ije_i
  IMPLICIT NONE
  INTEGER     :: p_int, j_int ! polynom. degrees
  REAL(dp)    :: p_dp, j_dp
  REAL(dp)    :: kx, ky, z

  !! Electrons linear coefficients for moment RHS !!!!!!!!!!
  DO ip = ips_e, ipe_e
    p_int= parray_e(ip)   ! Hermite degree
    p_dp = REAL(p_int,dp) ! REAL of Hermite degree
    DO ij = ijs_e, ije_e
      j_int= jarray_e(ij)   ! Laguerre degree
      j_dp = REAL(j_int,dp) ! REAL of Laguerre degree
      xnepj(ip,ij) = taue_qe * 2._dp * (p_dp + j_dp + 1._dp)
    ENDDO
  ENDDO
  DO ip = ips_e, ipe_e
    p_int= parray_e(ip)   ! Hermite degree
    p_dp = REAL(p_int,dp) ! REAL of Hermite degree
    xnepp1j(ip) = sqrtTaue_qe * SQRT(p_dp + 1_dp)
    xnepm1j(ip) = sqrtTaue_qe * SQRT(p_dp)
    xnepp2j(ip) = taue_qe * SQRT((p_dp + 1._dp) * (p_dp + 2._dp))
    xnepm2j(ip) = taue_qe * SQRT(p_dp * (p_dp - 1._dp))
  ENDDO
  DO ij = ijs_e, ije_e
    j_int= jarray_e(ij)   ! Laguerre degree
    j_dp = REAL(j_int,dp) ! REAL of Laguerre degree
    xnepjp1(ij) = -taue_qe * (j_dp + 1._dp)
    xnepjm1(ij) = -taue_qe * j_dp
  ENDDO
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Ions linear coefficients for moment RHS !!!!!!!!!!
  DO ip = ips_i, ipe_i
    p_int= parray_i(ip)   ! Hermite degree
    p_dp = REAL(p_int,dp) ! REAL of Hermite degree
    DO ij = ijs_i, ije_i
      j_int= jarray_i(ij)   ! Laguerre degree
      j_dp = REAL(j_int,dp) ! REAL of Laguerre degree
      xnipj(ip,ij) = taui_qi * 2._dp * (p_dp + j_dp + 1._dp)
    ENDDO
  ENDDO
  DO ip = ips_i, ipe_i
    p_int= parray_i(ip)   ! Hermite degree
    p_dp = REAL(p_int,dp) ! REAL of Hermite degree
    xnipp1j(ip) = sqrtTaui_qi * SQRT(p_dp + 1._dp)
    xnipm1j(ip) = sqrtTaui_qi * SQRT(p_dp)
    xnipp2j(ip) = taui_qi * SQRT((p_dp + 1._dp) * (p_dp + 2._dp))
    xnipm2j(ip) = taui_qi * SQRT(p_dp * (p_dp - 1._dp))
  ENDDO
  DO ij = ijs_i, ije_i
    j_int= jarray_i(ij)   ! Laguerre degree
    j_dp = REAL(j_int,dp) ! REAL of Laguerre degree
    xnipjp1(ij) = -taui_qi * (j_dp + 1._dp)
    xnipjm1(ij) = -taui_qi * j_dp
  ENDDO
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! ES linear coefficients for moment RHS !!!!!!!!!!
  DO ip = ips_i, ipe_i
    p_int= parray_i(ip)   ! Hermite degree
    DO ij = ijs_i, ije_i
      j_int= jarray_i(ij)   ! REALof Laguerre degree
      j_dp = REAL(j_int,dp) ! REALof Laguerre degree
      !! Electrostatic potential pj terms
      IF (p_int .EQ. 0) THEN ! kronecker p0
        xphij(ip,ij)    = eta_n + 2.*j_dp*eta_T
        xphijp1(ip,ij)  =-eta_T*(j_dp+1._dp)
        xphijm1(ip,ij)  =-eta_T* j_dp
      ELSE IF (p_int .EQ. 2) THEN ! kronecker p2
        xphij(ip,ij)    = eta_T/SQRT2
        xphijp1(ip,ij)  = 0._dp; xphijm1(ip,ij)  = 0._dp;
      ELSE
        xphij(ip,ij)    = 0._dp; xphijp1(ip,ij)  = 0._dp
        xphijm1(ip,ij)  = 0._dp;
      ENDIF
    ENDDO
  ENDDO

END SUBROUTINE compute_lin_coeff


END MODULE numerics
