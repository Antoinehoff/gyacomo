MODULE coeff

! this module contains routines to compute the Laguerre-Laguerre products coefficients
! dnjs for the nonlinear term. ALL2L(n,j,s,0) = dnjs
! Lpjl(p,j,l) returns the l-th coefficient of the j-th order associated Laguerre_{p-1/2} poynomial
  ! the canonical basis (p=0) writes L_j(x) = sum_{l=0}^j c_l x^l

USE PREC_CONST
use BASIC
USE FMZM

PUBLIC

CONTAINS

  !-----------------------------------------------------------
  !> @brief 
  !> Computes the associated-Laguerre/Laguerre to Laguerre basis transformation coefficient \f$ \overline{d}^{|m|}_{nkf} \f$,
  !>
  !> \f[ \overline{d}^{|m|}_{nkf} = \sum_{n_1 = 0}^n \sum_{k_1 =0}^k   \sum_{f_1=0}^f L_{kk_1}^{-1/2}  L_{nn_1}^{|m|-1/2}  L_{f f_1}^{-1/2} (n_1 + k_1 + f_1)!. \f]
  !> @todo Checked against Mathematica
  !> @param[in] n
  !> @param[in] k 
  !> @param[in] s 
  !> @param[in] m 
  !> @param[out] XOUT 
  !-----------------------------------------------------------
    FUNCTION ALL2L(n,k,s,m) RESULT(XOUT)
      !
      ! Computes the associated-Laguerre/Laguerre to Laguerre coefficient, i.e. 
      !
      !    L_n^m L_k = sum_{s=0}^{n+k}\overline{d}_{nks}^m L_s
      !
      ! 
      ! 
      IMPLICIT NONE    
      !
      INTEGER, intent(in) :: n,k,s,m  ! ... degree and orders of Laguerre polynomials
      TYPE(FM)  :: XOUT ! ... coefficient values
      TYPE(FM), SAVE :: AUXN,AUXK,AUXS
      INTEGER :: n1,k1,s1
      !
      CALL FM_ENTER_USER_FUNCTION(XOUT)
      !              Compute coefficients
      !
      XOUT = TO_FM('0.0')
      !
      DO n1=0,n 
         AUXN =Lpjl(REAL(m,xp)-0.5_xp,REAL(n,xp),REAL(n1,xp))
         DO k1=0,k
            AUXK =Lpjl(-0.5_xp,REAL(k,xp),REAL(k1,xp))
            DO s1=0,s
               AUXS = Lpjl(-0.5_xp,REAL(s,xp),REAL(s1,xp))
               XOUT = XOUT + FACTORIAL(TO_FM(n1 + k1 + s1 ))*AUXN*AUXK*AUXS 
            ENDDO
         ENDDO
      ENDDO
      !
      CALL FM_EXIT_USER_FUNCTION(XOUT)
      !
    END FUNCTION ALL2L

    !-----------------------------------------------------------
    !> @brief
    !> Computes the associated Laguerre serie coefficient  \f$ L_{jl}^p \f$
    !>
    !> \f[ L_{jl}^{p}  = \frac{(-1)^l (p + j + 1/2)!}{(j - l)!(l +p + 1/2 )! l!}. \f]
    !>  such that
    !> \f[ L^{p+1/2}_j(x) = \sum_{l=0}^j L^p_{jl} x^l \f]
    !> @param[in] XOUT
    !-----------------------------------------------------------
    FUNCTION Lpjl(p,j,l) RESULT(XOUT)
      !
      ! Computes the associated Laguerre serie coefficient,
      !

      IMPLICIT NONE

      !
      REAL(xp), intent(in) :: p,j,l
      TYPE(FM) :: XOUT
      !
      CALL FM_ENTER_USER_FUNCTION(XOUT)

      !            compute coeff

      XOUT = TO_FM('0.0')

      XOUT = TO_FM(((-1)**l))*FACTORIAL(TO_FM(2*p + 2*j + 1)/2)/&
            (FACTORIAL(TO_FM( j -l ))*FACTORIAL(TO_FM(2*l + 2*p  + 1)/2)*FACTORIAL(TO_FM(l)))

      CALL FM_EXIT_USER_FUNCTION(XOUT)

    END FUNCTION Lpjl
    
END MODULE coeff
