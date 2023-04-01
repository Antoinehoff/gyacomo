module DLRA
    USE prec_const
    USE lapack
    implicit none

    PUBLIC :: filter_singular_value_ky_pj

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
end module DLRA