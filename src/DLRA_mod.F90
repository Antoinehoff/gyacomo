module DLRA
    USE prec_const
    implicit none

    PUBLIC :: filter_singular_value_ky_pj, test_SVD

    CONTAINS

    SUBROUTINE filter_singular_value_ky_pj(nsv,array_ky_pj)
        IMPLICIT NONE
        ! ARGUMENTS
        INTEGER, INTENT(IN) :: nsv                                ! number of singular values to keep
        COMPLEX(xp), DIMENSION(:,:), INTENT(INOUT) :: array_ky_pj ! Array to filter
        !

        ! Singular value decomposition
        ! CALL SVD(array_ky_pj,singular_values)
    END SUBROUTINE

    SUBROUTINE test_svd
#ifdef TEST_SVD
        ! Program to perform Singular Value Decomposition (SVD)
        ! using LAPACK library
        
        ! Specify the dimensions of the input matrix A
        INTEGER, PARAMETER :: m = 3, n = 2
        
        ! Declare the input matrix A
        COMPLEX, DIMENSION(m,n) :: A
        
        ! OUTPUT
        COMPLEX, DIMENSION(m,m) :: U
        REAL,    DIMENSION(MIN(m,n)) :: S
        COMPLEX, DIMENSION(n,n) :: VT

        ! local variables
        INTEGER :: lda, ldu, ldvt, info, lwork, i, j
        COMPLEX, DIMENSION(:), ALLOCATABLE :: work
        REAL,    DIMENSION(:), ALLOCATABLE :: rwork
    
        ! Set the leading dimensions for the input and output arrays
        lda = MAX(1, m)
        ldu = MAX(1, m)
        ldvt = MAX(1, n)
        
        ! Define the input matrix A
        A = RESHAPE((/ (1.0,0.1), (2.0,0.2), (3.0,0.3), (4.0,0.4), (5.0,0.5), (6.0,0.6) /), SHAPE(A))

        ! Print the input matrix A
        WRITE(*,*) 'Input matrix A = '
        DO i = 1, m
          WRITE(*,*) ('(',REAL(A(i,j)), AIMAG(A(i,j)),')', j=1,n)
        END DO

        ALLOCATE(work(5*n), rwork(5*n))
        ! Compute the optimal workspace size
        lwork = -1
        CALL CGESVD('A', 'A', m, n, A, lda, S, U, ldu, VT, ldvt, work, lwork, rwork, info)

        ! Allocate memory for the workspace arrays
        lwork = CEILING(REAL(work(1)), KIND=SELECTED_REAL_KIND(1, 6))

        ! Compute the SVD of A using the LAPACK subroutine CGESVD
        CALL CGESVD('A', 'A', m, n, A, lda, S, U, ldu, VT, ldvt, work, lwork, rwork, info)
        
        ! Print the results
        WRITE(*,*) 'U = '
        DO i = 1, m
            WRITE(*,*) ('(',REAL(U(i,j)), AIMAG(U(i,j)),')', j=1,m)
        END DO
        WRITE(*,*)
        WRITE(*,*) 'S = '
        WRITE(*,'(2F8.3)') (S(i), i=1,n)
        WRITE(*,*)
        WRITE(*,*) 'VT = '
        DO i = 1, n
            WRITE(*,*) ('(',REAL(VT(i,j)), AIMAG(VT(i,j)),')', j=1,n)
        END DO
      
        ! Reconstruct A from its SVD
        A = MATMUL(U, MATMUL(diagmat(S,m,n), VT))
      
        ! Print the reconstructed matrix A
        WRITE(*,*) 'Reconstructed matrix A = '
        DO i = 1, m
          WRITE(*,*) ('(',REAL(A(i,j)), AIMAG(A(i,j)),')', j=1,n)
        END DO
        stop
#endif
      END SUBROUTINE test_svd
      
    ! SUBROUTINE test_svd
    !     ! Program to perform Singular Value Decomposition (SVD)
    !     ! using LAPACK library
        
    !     ! Specify the dimensions of the input matrix A
    !       INTEGER, PARAMETER :: m = 3, n = 2
    !       INTEGER, PARAMETER :: lda = m
        
    !     ! Declare the input matrix A
    !       REAL, DIMENSION(lda,n) :: A
        
    !     ! Specify the dimensions of the output matrices
    !       INTEGER, PARAMETER :: ldu = m, ldvt = n
    !       INTEGER, PARAMETER :: lwork = 5*n
    !       REAL, DIMENSION(ldu,m) :: U
    !       REAL, DIMENSION(n) :: S
    !       REAL, DIMENSION(ldvt,n) :: VT
    !       REAL, DIMENSION(lwork) :: work
    !       INTEGER :: info,i,j
        
    !     ! Define the input matrix A
    !       A = RESHAPE((/ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 /), SHAPE(A))
        
    !     ! Compute the SVD of A using the LAPACK subroutine SGESVD
    !       CALL SGESVD('A', 'A', m, n, A, lda, S, U, ldu, VT, ldvt, work, lwork, info)
        
    !     ! Print the results
    !       WRITE(*,*) 'U = '
    !       DO i = 1, m
    !          WRITE(*,'(6F8.3)') (U(i,j), j=1,m)
    !       END DO
    !       WRITE(*,*)
    !       WRITE(*,*) 'S = '
    !       WRITE(*,'(2F8.3)') (S(i), i=1,n)
    !       WRITE(*,*)
    !       WRITE(*,*) 'VT = '
    !       DO i = 1, n
    !          WRITE(*,'(6F8.3)') (VT(i,j), j=1,n)
    !       END DO

    !     ! Reconstruct A from its SVD
    !     A = MATMUL(U, MATMUL(diagmat(S,m,n), TRANSPOSE(VT)))

    !     ! Print the reconstructed matrix A
    !     WRITE(*,*) 'Reconstructed matrix A = '
    !     DO i = 1, m
    !     WRITE(*,'(2X, 3F8.3)') (A(i,j), j=1,n)
    !     END DO
    !     stop
    !     END SUBROUTINE test_svd

    ! SUBROUTINE test_svd
    !     IMPLICIT NONE
    !     INTEGER, PARAMETER :: m = 3, n = 2
    !     COMPLEX(xp), DIMENSION(m, n) :: A
    !     COMPLEX(xp), DIMENSION(m, m) :: U
    !     COMPLEX(xp), DIMENSION(n, n) :: VT
    !     REAL(xp), DIMENSION(MIN(m,n)) :: S
    !     INTEGER :: i, j
        
    !     ! Initialize A
    !     A = RESHAPE((/ (1.0_xp, 2.0_xp), (3.0_xp, 4.0_xp), (5.0_xp, 6.0_xp),&
    !                    (1.5_xp, 2.5_xp), (3.5_xp, 4.5_xp), (5.5_xp, 6.5_xp) /), [m, n])

    !     ! Print input
    !     WRITE(*,*) "A = "
    !     DO i = 1, m
    !         WRITE(*,"(3F8.3)") (REAL(A(i,j)), AIMAG(A(i,j)), j=1,n)
    !     END DO
    !     WRITE(*,*)
        
    !     ! Call the SVD subroutine
    !     CALL svd(A, U, S, VT, m, n)
        
    !     ! Print the resultas
    !     WRITE(*,*) "U = "
    !     DO i = 1, m
    !         WRITE(*,"(3F8.3)") (REAL(U(i,j)), AIMAG(U(i,j)), j=1,m)
    !     END DO
    !     WRITE(*,*)
        
    !     WRITE(*,*) "S = ", S
    !     WRITE(*,*)
        
    !     WRITE(*,*) "VT = "
    !     DO i = 1, n
    !         WRITE(*,"(3F8.3)") (REAL(VT(i,j)), AIMAG(VT(i,j)), j=1,n)
    !     END DO
    ! END SUBROUTINE test_svd
    

    SUBROUTINE svd(A, U, S, VT, m, n)
      IMPLICIT NONE
      ! INPUT
      COMPLEX(xp), DIMENSION(m,n), INTENT(IN)  :: A
      INTEGER, INTENT(IN) :: m, n
      ! OUTPUT
      COMPLEX(xp), DIMENSION(m,m), INTENT(OUT) :: U
      REAL(xp),    DIMENSION(MIN(m,n)), INTENT(OUT) :: S
      COMPLEX(xp), DIMENSION(n,n), INTENT(OUT) :: VT
      ! local variables
      INTEGER :: lda, ldu, ldvt, info, lwork
      COMPLEX(xp), DIMENSION(:), ALLOCATABLE :: work
      REAL(xp),    DIMENSION(:), ALLOCATABLE :: rwork
#ifdef LAPACKDIR      
        ! Set the leading dimensions for the input and output arrays
        lda = MAX(1, m)
        ldu = MAX(1, m)
        ldvt = MAX(1, n)
    
        ! Compute the optimal workspace size
        lwork = -1
        CALL ZGESVD('A', 'A', m, n, A, lda, S, U, ldu, VT, ldvt, work, lwork, rwork, info)
    
        ! Allocate memory for the workspace arrays
        lwork = CEILING(REAL(work(1)), KIND=SELECTED_REAL_KIND(1, 6))
        ALLOCATE(work(lwork), rwork(5*MIN(m,n)))
    
        ! Compute the SVD
        CALL ZGESVD('A', 'A', m, n, A, lda, S, U, ldu, VT, ldvt, work, lwork, rwork, info)
    
        ! Free the workspace arrays
        DEALLOCATE(work, rwork)
#endif
    END SUBROUTINE svd
    
      
    FUNCTION diagmat(S, m, n) RESULT(D)
      IMPLICIT NONE
      ! INPUT
      REAL, DIMENSION(:), INTENT(IN) :: S
      INTEGER, INTENT(IN) :: m, n
      ! OUTPUT
      COMPLEX, DIMENSION(m,n) :: D
      ! Local variables
      INTEGER :: i, j
  
      ! Initialize the output array to zero
      D = 0.0
  
      ! Fill the diagonal elements of the output array with the input vector S
      DO i = 1, MIN(m,n)
          D(i,i) = S(i)
      END DO
  
  END FUNCTION diagmat

end module DLRA