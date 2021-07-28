SUBROUTINE lin_coeff_and_geometry
  USE array, ONLY: xnepj, xnepp1j, xnepm1j, xnepp2j, xnepm2j, xnepjp1, xnepjm1, &
                   xnipj, xnipp1j, xnipm1j, xnipp2j, xnipm2j, xnipjp1, xnipjm1, &
                   xphij, xphijp1, xphijm1, Ckxky
  USE model, ONLY: taue_qe, taui_qi, sqrtTaue_qe, sqrtTaui_qi, eta_T, eta_n
  USE prec_const
  USE grid,  ONLY: parray_e, parray_i, jarray_e, jarray_i, &
                   ip,ij, ips_e,ip_e, ips_i,ipe_i, ijs_e,ije_e, ijs_i,ije_i,&
                   kxarray, kyarray, zarray, &
                   ikx,iky,iz, ikxs,ikxe, ikys,ikye, izs,ize
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
    xnepjp1(ij) = -taui_qi * (j_dp + 1._dp)
    xnepjm1(ij) = -taui_qi * j_dp
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
        xphijp1(ip,ij)  = eta_T*(j_dp+1._dp)
        xphijm1(ip,ij)  = eta_T* j_dp
      ELSE IF (p_int .EQ. 2) THEN ! kronecker p2
        xphij(ip,ij)    =-eta_T/SQRT2
        xphijp1(ip,ij)  = 0._dp; xphijm1(ip,ij)  = 0._dp;
      ELSE
        xphij(ip,ij)    = 0._dp; xphijp1(ip,ij)  = 0._dp
        xphijm1(ip,ij)  = 0._dp;
      ENDIF
    ENDDO
  ENDDO
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Curvature and geometric coefficients !!!!!!!!!!
  DO  iz = izs,ize
    z      = zarray(iz)
    DO ikx = ikxs,ikxe
      kx     = kxarray(ikx)   ! Poloidal wavevector
      DO iky = ikys,ikye
        ky     = kyarray(iky)   ! Toroidal wavevector
        Ckxky(ikx,iky,iz) = SIN(z)*imagu*kx +COS(z)*imagu*ky
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE lin_coeff_and_geometry
