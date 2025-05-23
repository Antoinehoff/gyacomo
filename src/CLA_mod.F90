module CLA
   ! This si the Computational Linear Algebra module.
   ! It contains tools as routines to DO singular value decomposition and filtering
   ! and matrices inversion to compute coefficient for the monomial
   ! closure
   USE prec_const
   USE parallel
   implicit none
   ! local variables for DLRA and SVD
   COMPLEX(xp), DIMENSION(:,:), ALLOCATABLE :: A_buff,Bf,Br ! buffer and full/reduced rebuilt matrices
   COMPLEX(xp), DIMENSION(:,:), ALLOCATABLE :: U            ! basis
   REAL(xp),    DIMENSION(:),   ALLOCATABLE :: Sf, Sr       ! full and reduced singular values
   COMPLEX(xp), DIMENSION(:,:), ALLOCATABLE :: VT           ! basis
   INTEGER :: lda, ldu, ldvt, info, lwork, i, m,n, nsv_filter
   COMPLEX(xp), DIMENSION(:), ALLOCATABLE :: work
   REAL(xp),    DIMENSION(:), ALLOCATABLE :: rwork

   ! local variables for monomial truncation (Hermite and Laguerre coeff)
    ! Coeff for the linear combination of the monomial truncation scheme
   REAL(xp), PUBLIC, PROTECTED, DIMENSION(:), ALLOCATABLE :: c_Pp1, c_Pp2, c_Jp1
   
#ifdef TEST_SVD
   ! These routines are meant for testing SVD filtering
   PUBLIC :: init_CLA, filter_sv_moments_ky_pj, test_SVD
#endif

   PUBLIC :: set_monomial_trunc_coeff

CONTAINS
   !--------------Routines to compute coeff for monomial truncation----------
   SUBROUTINE set_monomial_trunc_coeff(Pmax,Jmax)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: Pmax,Jmax ! truncation degree (p,j>d are not evolved)
      REAL(xp), DIMENSION(Pmax+3,Pmax+3) :: HH, iHH   ! Hermite matrix to invert  (+1 dim for p=0 and +2 dim for coeff P+1,P+2)
      REAL(xp), DIMENSION(Jmax+2,Jmax+2) :: LL, iLL   ! Laguerre matrix to invert (+1 dim for j=0 and +1 dim for coeff J+1)
      INTEGER :: i,n
      ! Hermite
      HH = 0._xp
      DO i = 1,Pmax+3
         n = i-1 !poly. degree
         HH(1:i,i) = get_hermite_coeffs(n)
      ENDDO
      CALL inverse_triangular(Pmax+3,HH,iHH)

      ! Laguerre
      LL = 0._xp
      DO i = 1,Jmax+2
         n = i-1 !poly. degree
         LL(1:i,i) = get_laguerre_coeffs(n)
      ENDDO
      CALL inverse_triangular(Jmax+2,LL,iLL)

      ! Allocate and fill the monomial truncation coefficients
      ALLOCATE(c_Pp1(Pmax+1))
      DO i = 1,Pmax+1
         c_Pp1(i) = -iHH(i,Pmax+2)/iHH(Pmax+2,Pmax+2)
      ENDDO
      ALLOCATE(c_Pp2(Pmax+1))
      DO i = 1,Pmax+1
         c_Pp2(i) = -iHH(i,Pmax+3)/iHH(Pmax+3,Pmax+3)
      ENDDO
      ALLOCATE(c_Jp1(Jmax+1))
      DO i = 1,Jmax+1
         c_Jp1(i) = 0*iLL(i,Jmax+2)/iLL(Jmax+2,Jmax+2)
      ENDDO
   END SUBROUTINE

   ! Invert an upper triangular matrix
   SUBROUTINE inverse_triangular(N,U,invU)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: N ! size of the matrix
      REAL(xp), DIMENSION(N,N), INTENT(IN)  :: U   ! Matrix to invert
      REAL(xp), DIMENSION(N,N), INTENT(OUT) :: invU! Result
      ! local variables
      INTEGER :: info
      invU = U
#ifdef LAPACKDIR
#ifdef SINGLE_PRECISION
      CALL STRTRI('U','N',N,invU,N,info)
#else
      CALL DTRTRI('U','N',N,invU,N,info)
#endif
#else
      ERROR STOP "Cannot use monomial truncation without LAPACK"
#endif
      IF (info .LT. 0) THEN
         print*, info
         ERROR STOP "STOP, LAPACK TRTRI: illegal value at index"
      ELSEIF(info .GT. 0) THEN
         ERROR STOP "STOP, LAPACK TRTRI: singular matrix"
      ENDIF
   END SUBROUTINE

   ! Compute the kth coefficient of the nth degree Hermite 
   ! adapted from LaguerrePoly.m by David Terr, Raytheon, 5-10-04
   FUNCTION get_hermite_coeffs(n) RESULT(hn)
      IMPLICIT NONE
      INTEGER,  INTENT(IN)     :: n ! Polyn. Deg. and index
      REAL(xp), DIMENSION(n+1) :: hn ! Hermite coefficients
      REAL(xp), DIMENSION(n+1)  :: hnm2, hnm1, tmp
      INTEGER :: k, e, factn
      IF (n .EQ. 0) THEN 
         hn = 1._xp
      ELSEIF (n .EQ. 1) THEN
         hn = [2._xp, 0._xp]
      ELSE
         hnm2      = 0._xp
         hnm2(n+1) = 1._xp
         hnm1      = 0._xp
         hnm1(n)   = 2._xp
         DO k = 2,n
             hn = 0._xp
             DO e = n-k+1,n,2
                 hn(e) = REAL(2*(hnm1(e+1) - (k-1)*hnm2(e)),xp)
               END DO
             hn(n+1) = REAL(-2*(k-1)*hnm2(n+1),xp)
             IF (k .LT. n) THEN
                 hnm2 = hnm1
                 hnm1 = hn
             END IF
            END DO
      END IF
      ! Normalize the polynomials 
      factn = 1
      DO k = 1,n
         factn = k*factn
      ENDDO
      hn = hn/SQRT(REAL(2**n*factn,xp))
      ! reverse order
      tmp = hn
      DO k = 1,n+1
         hn(n+1-(k-1)) = tmp(k)
      ENDDO
   END FUNCTION

   ! Compute the kth coefficient of the nth degree Laguerre 
   ! adapted from LaguerrePoly.m by David Terr, Raytheon, 5-11-04
   ! lnk = (-1)^k/k! binom(n,k), s.t. Ln(x) = sum_k lnk x^k
   FUNCTION get_laguerre_coeffs(n) RESULT(ln)
      IMPLICIT NONE
      INTEGER,  INTENT(IN)     :: n ! Polyn. Deg. and index
      REAL(xp), DIMENSION(n+1) :: ln ! Laguerre coefficients
      REAL(xp), DIMENSION(n+1) :: lnm2, lnm1, tmp
      INTEGER :: k, e
      IF (n .EQ. 0) THEN
          ln = 1._xp
      ELSEIF (n .EQ. 1) THEN
          ln = [-1._xp, 1._xp]
      ELSE
          lnm2      = 0._xp
          lnm2(n+1) = 1._xp
          lnm1      = 0._xp
          lnm1(n)   =-1._xp
          lnm1(n+1) = 1._xp
          DO k = 2, n
              ln = 0
              DO e = n-k+1, n
                  ln(e) = REAL((2*k-1)*lnm1(e) - lnm1(e+1) + (1-k)*lnm2(e),xp)
              END DO
              ln(n+1) = REAL((2*k-1)*lnm1(n+1) + (1-k)*lnm2(n+1),xp)
              ln = ln/REAL(k,xp)
              IF (k < n) THEN
                  lnm2 = lnm1
                  lnm1 = ln
              END IF
          END DO
      END IF
      ! reverse order
      tmp = ln
      DO k = 1,n+1
         ln(n+1-(k-1)) = tmp(k)
      ENDDO
   END FUNCTION

   !--------------Routines to test singular value filtering -----------------
#ifdef TEST_SVD
   SUBROUTINE init_CLA(m_,n_)
      USE basic
      IMPLICIT NONE
      ! ARGUMENTS
      INTEGER, INTENT(IN) :: m_,n_                      ! dimensions of the input array
      ! read the input
      INTEGER :: lun   = 90              ! File duplicated from STDIN
      NAMELIST /CLA/ nsv_filter
      READ(lun,CLA)
      m   = m_
      n   = n_
      info = 1
      ! Allocate the matrices
      ALLOCATE(A_buff(m,n),Bf(m,n),Br(m,n))
      ALLOCATE(U(m,m),Sf(MIN(m,n)),Sr(MIN(m,n)),VT(n,n))
      ! Set the leading dimensions for the input and output arrays
      lda  = MAX(1, m)
      ldu  = MAX(1, m)
      ldvt = MAX(1, n)
      ALLOCATE(work(5*n), rwork(5*MIN(m,n)))
      ! Compute the optimal workspace size
      lwork = -1
#ifdef SINGLE_PRECISION
      CALL CGESVD('A', 'A', m, n, A_buff, lda, Sf, U, ldu, VT, ldvt, work, lwork, rwork, info)
#else
      CALL ZGESVD('A', 'A', m, n, A_buff, lda, Sf, U, ldu, VT, ldvt, work, lwork, rwork, info)
#endif
      ! Allocate memory for the workspace arrays
      lwork = 2*CEILING(REAL(work(1)))
      DEALLOCATE(work)
      ALLOCATE(work(lwork))
   END SUBROUTINE init_CLA

   SUBROUTINE filter_sv_moments_ky_pj
      USE fields,           ONLY: moments
      USE grid,             ONLY: total_nky, total_np, total_nj,ngz_o2,&
         local_np, local_nj, local_nz, local_nkx, local_nky, local_na
      USE time_integration, ONLY: updatetlevel
      USE basic,            ONLY: start_chrono, stop_chrono, chrono_CLA
      IMPLICIT NONE

      ! Arguments
      INTEGER :: nsv_filter           ! number of singular values to keep
      ! Local variables
      COMPLEX(xp), DIMENSION(:,:), ALLOCATABLE :: moments_lky_lpj ! local ky and local pj data
      COMPLEX(xp), DIMENSION(:,:), ALLOCATABLE :: moments_gky_lpj ! global ky, local pj data
      COMPLEX(xp), DIMENSION(:,:), ALLOCATABLE :: moments_gky_gpj ! full gathered data for SVD (input of SVD)
      INTEGER :: ia,ix,iz, m,n, ip, ij, iy
      CALL start_chrono(chrono_CLA)
      ALLOCATE(moments_lky_lpj(local_nky,local_np*local_nj))
      ALLOCATE(moments_gky_lpj(total_nky,local_np*local_nj))
      ALLOCATE(moments_gky_gpj(total_nky,total_np*total_nj))
      DO iz = 1+ngz_o2,local_nz+ngz_o2
         DO ix = 1,local_nkx
            DO ia = 1,local_na
               ! Build the local slice explicitely
               DO ip = 1,local_np
                  DO ij = 1,local_nj
                     DO iy = 1,local_nky
                        moments_lky_lpj(iy,local_np*(ij-1)+ip) = moments(ia,ip+ngp_o2,ij+ngj_o2,iy,ix,iz,updatetlevel)
                     ENDDO
                  ENDDO
               ENDDO
               ! Gather ky data
               IF(num_procs_ky .GT. 1) THEN
                  ! MPI communication
               ELSE
                  moments_gky_lpj = moments_lky_lpj
               ENDIF
               ! Gather p data
               IF(num_procs_p .GT. 1) THEN
                  ! MPI communication
               ELSE
                  moments_gky_gpj = moments_gky_lpj
               ENDIF
               ! The process 0 performs the SVD
               IF(my_id .EQ. 0) THEN
                  m = total_nky
                  n = total_np*total_nj
                  nsv_filter = -1
                  CALL filter_singular_value(moments_gky_gpj)
               ENDIF
               ! Distribute ky data
               IF(num_procs_ky .GT. 1) THEN
                  ! MPI communication
               ELSE
                  moments_gky_lpj = moments_gky_gpj
               ENDIF
               ! Distribute p data
               IF(num_procs_p .GT. 1) THEN
                  ! MPI communication
               ELSE
                  moments_lky_lpj = moments_gky_lpj
               ENDIF
               ! Put back the data into the moments array
               DO ip = 1,local_np
                  DO ij = 1,local_nj
                     DO iy = 1,local_nky
                        moments(ia,ip+ngp_o2,ij+ngj_o2,iy,ix,iz,updatetlevel) = moments_lky_lpj(iy,local_np*(ij-1)+ip)
                     ENDDO
                  ENDDO
               ENDDO

            ENDDO
         ENDDO
      ENDDO
      CALL stop_chrono(chrono_CLA)
   END SUBROUTINE filter_sv_moments_ky_pj

   SUBROUTINE filter_singular_value(A)
      IMPLICIT NONE
      ! ARGUMENTS
      COMPLEX(xp), DIMENSION(:,:), INTENT(INOUT) :: A ! Array to filter
      ! copy of the input since it is changed in the svd procedure
      A_buff = A;
#ifdef SINGLE_PRECISION
      CALL CGESVD('A', 'A', m, n, A_buff, lda, Sf, U, ldu, VT, ldvt, work, lwork, rwork, info)
#else
      CALL ZGESVD('A', 'A', m, n, A_buff, lda, Sf, U, ldu, VT, ldvt, work, lwork, rwork, info)
#endif
      IF(info.GT.0) print*,"SVD did not converge, info=",info
      ! Filter the values
      Sr = 0._xp
      IF (nsv_filter .GT. 0) THEN
         DO i=1,MIN(m,n)-nsv_filter
            Sr(i) = Sf(i)
         ENDDO
      ELSE ! DO not filter IF nsv_filter<0
         Sr = Sf
      ENDIF
      ! Reconstruct A from its reduced SVD
      A = MATMUL(U, MATMUL(diagmat(Sr,m,n), VT)) ! reduced
   END SUBROUTINE filter_singular_value

   FUNCTION diagmat(S, m, n) RESULT(D)
      IMPLICIT NONE
      ! INPUT
      REAL(xp), DIMENSION(:), INTENT(IN) :: S
      INTEGER, INTENT(IN) :: m, n
      ! OUTPUT
      COMPLEX(xp), DIMENSION(m,n) :: D
      ! Local variables
      INTEGER :: i
      ! Initialize the output array to zero
      D = 0.0
      ! Fill the diagonal elements of the output array with the input vector S
      DO i = 1, MIN(m,n)
         D(i,i) = S(i)
      END DO
   END FUNCTION diagmat

   SUBROUTINE test_svd
      ! SpecIFy the dimensions of the input matrix A
      INTEGER :: m,n

      ! Declare the input matrix A and reconstructed matrix B
      COMPLEX(xp), DIMENSION(:,:), ALLOCATABLE :: A,B,A_buff

      ! OUTPUT
      COMPLEX(xp), DIMENSION(:,:), ALLOCATABLE :: U
      REAL(xp),    DIMENSION(:),   ALLOCATABLE :: S
      COMPLEX(xp), DIMENSION(:,:), ALLOCATABLE :: VT
      ! local variables
      INTEGER :: lda, ldu, ldvt, info, lwork, i, j
      COMPLEX(xp), DIMENSION(:), ALLOCATABLE :: work
      REAL(xp),    DIMENSION(:), ALLOCATABLE :: rwork

      m = 3
      n = 2
      ALLOCATE(A_buff(m,n),A(m,n),B(m,n),U(m,m),S(MIN(m,n)),VT(n,n))
      ! Set the leading dimensions for the input and output arrays
      lda  = MAX(1, m)
      ldu  = MAX(1, m)
      ldvt = MAX(1, n)

      ! Define the input matrix A
      A = RESHAPE((/ (1._xp,0.1_xp), (2._xp,0.2_xp), (3._xp,0.3_xp), (4._xp,0.4_xp), (5._xp,0.5_xp), (6._xp,0.6_xp) /), SHAPE(A))
      ! copy of the input since it is changed in the svd procedure
      A_buff = A;
      ! Print the input matrix A
      WRITE(*,*) 'Input matrix A = '
      DO i = 1, m
         WRITE(*,*) ('(',REAL(A(i,j)), AIMAG(A(i,j)),')', j=1,n)
      END DO

      ! CALL svd(A, U, S, VT, m, n)
      ALLOCATE(work(5*n), rwork(5*n))
      ! Compute the optimal workspace size
      lwork = -1
#ifdef SINGLE_PRECISION
      CALL CGESVD('A', 'A', m, n, A_buff, lda, S, U, ldu, VT, ldvt, work, lwork, rwork, info)
#else
      CALL ZGESVD('A', 'A', m, n, A_buff, lda, S, U, ldu, VT, ldvt, work, lwork, rwork, info)
#endif

      ! Allocate memory for the workspace arrays
      lwork = CEILING(REAL(work(1)), KIND=SELECTED_REAL_KIND(1, 6))

      ! Compute the SVD of A using the LAPACK subroutine CGESVD
#ifdef SINGLE_PRECISION
      CALL CGESVD('A', 'A', m, n, A_buff, lda, S, U, ldu, VT, ldvt, work, lwork, rwork, info)
#else
      CALL ZGESVD('A', 'A', m, n, A_buff, lda, S, U, ldu, VT, ldvt, work, lwork, rwork, info)
#endif
      ! Print the results
      ! WRITE(*,*) 'U = '
      ! DO i = 1, m
      !     WRITE(*,*) ('(',REAL(U(i,j)), AIMAG(U(i,j)),')', j=1,m)
      ! END DO
      ! WRITE(*,*)
      WRITE(*,*) 'S = '
      WRITE(*,'(2F8.3)') (S(i), i=1,n)
      WRITE(*,*)
      ! WRITE(*,*) 'VT = '
      ! DO i = 1, n
      !     WRITE(*,*) ('(',REAL(VT(i,j)), AIMAG(VT(i,j)),')', j=1,n)
      ! END DO

      ! Reconstruct A from its SVD
      B = MATMUL(U, MATMUL(diagmat(S,m,n), VT))
      ! Print the error with the intput matrix
      WRITE(*,*) '||A-USVT||=', sum(abs(A-B))
      stop
   END SUBROUTINE test_svd
#endif
END module CLA
