MODULE nonlinear
  USE array,       ONLY : dnjs, Sapj, kernel
  USE initial_par, ONLY : ACT_ON_MODES
  USE fourier,     ONLY : bracket_sum_r, bracket_sum_c, planf, planb, poisson_bracket_and_sum
  USE fields,      ONLY : phi, psi, moments
  USE grid,        ONLY: local_na, &
                         local_np,ngp,parray,pmax,&
                         local_nj,ngj,jarray,jmax, local_nj_offset, dmax,&
                         kyarray, AA_y, local_nky_ptr, local_nky_ptr_offset,inv_Ny,&
                         local_nkx_ptr,kxarray, AA_x, inv_Nx,&
                         local_nz,ngz,zarray,nzgrid
  USE model,       ONLY : LINEARITY, CLOS, NL_CLOS, EM
  USE prec_const,  ONLY : dp
  USE species,     ONLY : sqrt_tau_o_sigma
  USE time_integration, ONLY : updatetlevel
  use, intrinsic :: iso_c_binding

  IMPLICIT NONE

  INCLUDE 'fftw3-mpi.f03'

  COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE :: F_cmpx, G_cmpx
  COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE :: Fx_cmpx, Gy_cmpx
  COMPLEX(dp), DIMENSION(:,:), ALLOCATABLE :: Fy_cmpx, Gx_cmpx, F_conv_G
  INTEGER :: in, is, p_int, j_int, n_int
  INTEGER :: nmax, smax
  REAL(dp):: sqrt_p, sqrt_pp1
  PUBLIC  :: compute_Sapj, nonlinear_init

CONTAINS

SUBROUTINE nonlinear_init
  IMPLICIT NONE
  ALLOCATE( F_cmpx(local_nky_ptr,local_nkx_ptr))
  ALLOCATE( G_cmpx(local_nky_ptr,local_nkx_ptr))

  ALLOCATE(Fx_cmpx(local_nky_ptr,local_nkx_ptr))
  ALLOCATE(Gy_cmpx(local_nky_ptr,local_nkx_ptr))
  ALLOCATE(Fy_cmpx(local_nky_ptr,local_nkx_ptr))
  ALLOCATE(Gx_cmpx(local_nky_ptr,local_nkx_ptr))

  ALLOCATE(F_conv_G(local_nky_ptr,local_nkx_ptr))
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
      Sapj = 0._dp
    CASE DEFAULT
      ERROR STOP '>> ERROR << Linearity not recognized '
  END SELECT
END SUBROUTINE compute_Sapj

SUBROUTINE compute_nonlinear
  IMPLICIT NONE
  INTEGER :: iz,ij,ip,eo,ia,ikx,iky,izi,ipi,iji,ini,isi
  DO iz = 1,local_nz
    izi = iz + ngz/2
    DO ij = 1,local_nj ! Loop over Laguerre moments
      iji = ij + ngj/2
      j_int=jarray(iji)
      DO ip = 1,local_np ! Loop over Hermite moments
        ipi = ip + ngp/2
        p_int    = parray(ipi)
        sqrt_p   = SQRT(REAL(p_int,dp))
        sqrt_pp1 = SQRT(REAL(p_int,dp) + 1._dp)
        eo       = min(nzgrid,MODULO(parray(ip),2)+1)
        DO ia = 1,local_na
          IF((CLOS .NE. 1) .OR. (p_int+2*j_int .LE. dmax)) THEN !compute for every moments except for closure 1
            ! Set non linear sum truncation
            IF (NL_CLOS .EQ. -2) THEN
              nmax = Jmax
            ELSEIF (NL_CLOS .EQ. -1) THEN
              nmax = Jmax-j_int
            ELSE
              nmax = min(NL_CLOS,Jmax-j_int)
            ENDIF
            bracket_sum_r = 0._dp ! initialize sum over real nonlinear term
            DO in = 1,nmax+1 ! Loop over laguerre for the sum
              ini = in+ngj/2
  !-----------!! ELECTROSTATIC CONTRIBUTION
              ! First convolution terms
              F_cmpx(:,:) = phi(:,:,izi) * kernel(ia,ini,:,:,izi,eo)
              ! Second convolution terms
              G_cmpx = 0._dp ! initialization of the sum
              smax   = MIN( (in-1)+(ij-1), Jmax );
              DO is = 1, smax+1 ! sum truncation on number of moments
                isi = is + ngj/2
                G_cmpx(:,:) = G_cmpx(:,:) + &
                  dnjs(in,ij,is) * moments(ia,ipi,isi,:,:,izi,updatetlevel)
              ENDDO
              ! this function add its result to bracket_sum_r
              CALL poisson_bracket_and_sum(kyarray,kxarray,inv_Ny,inv_Nx,AA_y,AA_x,local_nky_ptr,local_nkx_ptr,F_cmpx,G_cmpx,bracket_sum_r)
  !-----------!! ELECTROMAGNETIC CONTRIBUTION -sqrt(tau)/sigma*{Sum_s dnjs [sqrt(p+1)Nap+1s + sqrt(p)Nap-1s], Kernel psi}
              IF(EM) THEN
              ! First convolution terms
              F_cmpx(:,:) = -sqrt_tau_o_sigma(ia) * psi(:,:,izi) * kernel(ia,ini,:,:,izi,eo)
              ! Second convolution terms
              G_cmpx = 0._dp ! initialization of the sum
              smax   = MIN( (in-1)+(ij-1), Jmax );
              DO is = 1, smax+1 ! sum truncation on number of moments
                isi = is + ngj/2
                G_cmpx(:,:)  = G_cmpx(:,:) + &
                  dnjs(in,ij,is) * (sqrt_pp1*moments(ia,ipi+1,isi,:,:,izi,updatetlevel)&
                                   +sqrt_p  *moments(ia,ipi-1,isi,:,:,izi,updatetlevel))
              ENDDO
              ! this function add its result to bracket_sum_r
              CALL poisson_bracket_and_sum(kyarray,kxarray,inv_Ny,inv_Nx,AA_y,AA_x,local_nky_ptr,local_nkx_ptr,F_cmpx,G_cmpx,bracket_sum_r)
              ENDIF
            ENDDO
            ! Put the real nonlinear product into k-space
            call fftw_mpi_execute_dft_r2c(planf, bracket_sum_r, bracket_sum_c)
            ! Retrieve convolution in input format and apply anti aliasing
            DO ikx = 1,local_nkx_ptr
              DO iky = 1,local_nky_ptr
                Sapj(ia,ip,ij,iky,ikx,iz) = bracket_sum_c(ikx,iky)*AA_x(ikx)*AA_y(iky)
              ENDDO
            ENDDO
          ELSE
            Sapj(ia,ip,ij,:,:,iz) = 0._dp
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END SUBROUTINE compute_nonlinear
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE nonlinear
