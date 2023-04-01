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

    SUBROUTINE test_SVD
        ! Program to perform Singular Value Decomposition (SVD)
        ! using LAPACK library

        ! Specify the dimensions of the input matrix A
        INTEGER, PARAMETER :: m = 3, n = 2
        INTEGER, PARAMETER :: lda = m

        ! Declare the input matrix A
        REAL, DIMENSION(lda,n) :: A

        ! Specify the dimensions of the output matrices
        INTEGER, PARAMETER :: ldu = m, ldvt = n
        INTEGER, PARAMETER :: lwork = 5*n
        REAL, DIMENSION(ldu,m) :: U
        REAL, DIMENSION(n) :: S
        REAL, DIMENSION(ldvt,n) :: VT
        REAL, DIMENSION(lwork) :: work
        INTEGER :: info, i,j

        ! Define the input matrix A
        A = RESHAPE((/ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 /), SHAPE(A))

        WRITE(*,*) 'Input matrix A = '

        DO i = 1, m
            WRITE(*,'(2X, 3F8.3)') (A(i,j), j=1,n)
        END DO

        ! Compute the SVD of A using the LAPACK subroutine SGESVD
        CALL SGESVD('A', 'A', m, n, A, lda, S, U, ldu, VT, ldvt, work, lwork, info)

        ! Print the results
        WRITE(*,*) 'U = '
        DO i = 1, m
            WRITE(*,'(6F8.3)') (U(i,j), j=1,m)
        END DO
        WRITE(*,*)
        WRITE(*,*) 'S = '
        WRITE(*,'(2F8.3)') (S(i), i=1,n)
        WRITE(*,*)
        WRITE(*,*) 'VT = '
        DO i = 1, n
            WRITE(*,'(6F8.3)') (VT(i,j), j=1,n)
        END DO

        ! Reconstruct A from its SVD
        A = MATMUL(U, MATMUL(diagmat(S,m,n), TRANSPOSE(VT)))
        ! Print the reconstructed matrix A

        WRITE(*,*) 'Reconstructed matrix A = '

        DO i = 1, m
            WRITE(*,'(2X, 3F8.3)') (A(i,j), j=1,n)
        END DO

        print*, "this was a test of the SVD using LAPACK. End run."
        stop
    END SUBROUTINE test_SVD

    FUNCTION diagmat(v, m, n) RESULT(A)
        REAL, DIMENSION(:), INTENT(IN) :: v
        INTEGER, INTENT(IN) :: m, n
        REAL, DIMENSION(m,n) :: A
        INTEGER :: i, j
      
        A = 0.0   ! Initialize A to a zero matrix
      
        DO i = 1, MIN(m,n)
          A(i,i) = v(i)   ! Set the diagonal elements of A to the values in v
        END DO
      
      END FUNCTION diagmat

end module DLRA