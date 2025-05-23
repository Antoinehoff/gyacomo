MODULE nonlinear
  USE array,       ONLY : dnjs, Sapj, kernel
  USE fourier,     ONLY : bracket_sum_r, bracket_sum_c, planf, planb, poisson_bracket_and_sum,&
                          apply_inv_ExB_NL_factor
  USE fields,      ONLY : phi, psi, moments
  USE grid,        ONLY : local_na, &
                         local_np,ngp_o2,parray,pmax,&
                         local_nj,ngj_o2,jarray,jmax, local_nj_offset,&
                         kyarray, AA_y, local_nky, inv_Ny,&
                         total_nkx,kxarray, AA_x, inv_Nx,local_nx, Ny, &
                         local_nz,ngz_o2,zarray,nzgrid, deltakx, deltaky, iky0, contains_kx0, contains_ky0
  USE model,       ONLY : LINEARITY, EM, ikxZF, ZFrate, ZF_ONLY, ExB_NL_CORRECTION
  USE closure,     ONLY : evolve_mom, nmaxarray
  USE prec_const,  ONLY : xp, imagu
  USE species,     ONLY : sqrt_tau_o_sigma
  USE time_integration, ONLY : updatetlevel
  USE ExB_shear_flow,   ONLY : ExB_NL_factor, inv_ExB_NL_factor
  use, intrinsic :: iso_c_binding

  IMPLICIT NONE

  INCLUDE 'fftw3-mpi.f03'

  COMPLEX(xp), DIMENSION(:,:), ALLOCATABLE :: F_cmpx, G_cmpx
  INTEGER :: in, is, p_int, j_int, n_int
  INTEGER :: smax
  REAL(xp):: sqrt_p, sqrt_pp1, inv_kx_ZF
  PUBLIC  :: compute_Sapj, nonlinear_init

CONTAINS

SUBROUTINE nonlinear_init
  IMPLICIT NONE
  ALLOCATE( F_cmpx(local_nky,total_nkx))
  ALLOCATE( G_cmpx(local_nky,total_nkx))
  inv_kx_ZF = ikxZF/deltakx
END SUBROUTINE nonlinear_init

SUBROUTINE compute_Sapj
  IMPLICIT NONE
  ! This routine is meant to compute the non linear term for each specie and degree
  !! In real space Sapj ~ b*(grad(phi) x grad(g)) which in moments in fourier becomes
  !! Sapj = Sum_n (ikx Kn phi)#(iky Sum_s d_njs Naps) - (iky Kn phi)#(ikx Sum_s d_njs Naps)
  !! where # denotes the convolution.
  SELECT CASE(LINEARITY)
    CASE ('nonlinear')
      CALL compute_nonlinear
    CASE ('linear')
      Sapj = 0._xp
    CASE DEFAULT
      ERROR STOP '>> ERROR << Linearity not recognized '
  END SELECT
END SUBROUTINE compute_Sapj

! Compute the poisson bracket {F,G}
SUBROUTINE compute_nonlinear
  IMPLICIT NONE
  INTEGER :: iz,ij,ip,eo,ia,ikx,iky,izi,ipi,iji,ini,isi
  INTEGER :: ikxExBp, ikxExBn ! Negative and positive ExB flow indices
  ! COMPLEX(xp), DIMENSION(Ny/2+1,local_nx) :: invinvfactor ! TEST
  z:DO iz = 1,local_nz
    izi = iz + ngz_o2
    j:DO ij = 1,local_nj ! Loop over Laguerre moments
      iji = ij + ngj_o2
      j_int=jarray(iji)
      p:DO ip = 1,local_np ! Loop over Hermite moments
        ipi = ip + ngp_o2
        IF(evolve_mom(ipi,iji)) THEN !compute for every moments except for closure 1
        p_int    = parray(ipi)
        sqrt_p   = SQRT(REAL(p_int,xp))
        sqrt_pp1 = SQRT(REAL(p_int,xp) + 1._xp)
        eo       = min(nzgrid,MODULO(parray(ip),2)+1)
        a:DO ia = 1,local_na
            ! Set non linear sum truncation
            bracket_sum_r = 0._xp ! initialize sum over real nonlinear term
            n:DO in = 1,nmaxarray(ij)+1 ! Loop over laguerre for the sum
              ini = in+ngj_o2
  !-----------!! ELECTROSTATIC CONTRIBUTION
              ! First convolution terms
              DO ikx = 1,total_nkx
                DO iky = 1,local_nky
                  F_cmpx(iky,ikx) = phi(iky,ikx,izi) * kernel(ia,ini,iky,ikx,izi,eo)
                ENDDO
              ENDDO
              ! ExB shearing as a additional zonal mode in the ES potential
              IF(abs(ZFrate) .GT. 0) THEN
                ikxExBp = ikxZF + 1
                ikxExBn = total_nkx - (ikxExBp - 2)
                IF(ZF_ONLY) F_cmpx = 0._xp
                IF(contains_kx0 .AND. contains_ky0) THEN
                  F_cmpx(iky0,ikxExBp) = F_cmpx(iky0,ikxExBp) + 8._xp*ZFrate*inv_kx_ZF * kernel(ia,ini,iky0,ikxExBp,izi,eo)
                  F_cmpx(iky0,ikxExBn) = F_cmpx(iky0,ikxExBn) + 8._xp*ZFrate*inv_kx_ZF * kernel(ia,ini,iky0,ikxExBn,izi,eo)
                ENDIF
              ENDIF
              ! Second convolution terms
              G_cmpx = 0._xp ! initialization of the sum
              ! We set smax as min(n+j,Jmax) for safety reason but nmaxarray may prevent it as well (see anti-Laguerre-aliasing truncation)
              smax   = MIN( jarray(ini)+jarray(iji), jmax );
              s1:DO is = 1, smax+1 ! sum truncation on number of moments
                isi = is + ngj_o2
                G_cmpx(:,:) = G_cmpx(:,:) + &
                  dnjs(in,ij,is) * moments(ia,ipi,isi,:,:,izi,updatetlevel)
              ENDDO s1
              ! this function adds its result to bracket_sum_r
                CALL poisson_bracket_and_sum( kyarray,kxarray,inv_Ny,inv_Nx,AA_y,AA_x,&
                                              local_nky,total_nkx,F_cmpx,G_cmpx,&
                                              ExB_NL_CORRECTION, ExB_NL_factor, bracket_sum_r)
  !-----------!! ELECTROMAGNETIC CONTRIBUTION -sqrt(tau)/sigma*{Sum_s dnjs [sqrt(p+1)Nap+1s + sqrt(p)Nap-1s], Kernel psi}
              IF(EM) THEN
                ! First convolution terms
                F_cmpx(:,:) = -sqrt_tau_o_sigma(ia) * psi(:,:,izi) * kernel(ia,ini,:,:,izi,eo)
                ! Second convolution terms
                G_cmpx = 0._xp ! initialization of the sum
                s2:DO is = 1, smax+1 ! sum truncation on number of moments
                  isi = is + ngj_o2
                  G_cmpx(:,:)  = G_cmpx(:,:) + &
                    dnjs(in,ij,is) * (sqrt_pp1*moments(ia,ipi+1,isi,:,:,izi,updatetlevel)&
                                    +sqrt_p  *moments(ia,ipi-1,isi,:,:,izi,updatetlevel))
                ENDDO s2
                ! this function adds its result to bracket_sum_r
                CALL poisson_bracket_and_sum( kyarray,kxarray,inv_Ny,inv_Nx,AA_y,AA_x,&
                                              local_nky,total_nkx,F_cmpx,G_cmpx,&
                                              ExB_NL_CORRECTION, ExB_NL_factor, bracket_sum_r)
              ENDIF
            ENDDO n
            ! Apply the ExB shearing rate factor before going back to k-space
            IF (ExB_NL_CORRECTION) THEN
              CALL apply_inv_ExB_NL_factor(bracket_sum_r,inv_ExB_NL_factor)
            ENDIF
            ! Put the real nonlinear product back into k-space
#ifdef SINGLE_PRECISION
            call fftwf_mpi_execute_dft_r2c(planf, bracket_sum_r, bracket_sum_c)
#else
            call  fftw_mpi_execute_dft_r2c(planf, bracket_sum_r, bracket_sum_c)
#endif
            ! Retrieve convolution in input format and apply anti aliasing
            DO ikx = 1,total_nkx
              DO iky = 1,local_nky
                Sapj(ia,ip,ij,iky,ikx,iz) = bracket_sum_c(ikx,iky)*AA_x(ikx)*AA_y(iky)
              ENDDO
            ENDDO
          ENDDO a
          ELSE
            Sapj(:,ip,ij,:,:,iz) = 0._xp
          ENDIF
      ENDDO p
    ENDDO j
  ENDDO z
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE compute_nonlinear
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE nonlinear
