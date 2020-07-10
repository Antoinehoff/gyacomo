MODULE coeff

! this module contains routines to compute normalization coefficients and velocity integrals 
! 

USE PREC_CONST 
use BASIC
USE MODEL
USE FMZM
USE basis_transformation

PUBLIC

CONTAINS 

  !-----------------------------------------------------------
  !> @brief 
  !> Computes the associated-Laguerre/Laguerre/x to Laguerre basis transformation coefficient \f$ d^{|m|}_{nks} \f$,
  !>
  !> \f[  L_n^{|m|}(x) L_k(x) x^{|m|} =\sum_{s=0}^{n+|m| +k} d^{|m|}_{nks} L_s(x),\f]
  !> where 
  !> \f[d^{|m|}_{nks} = \int_0^\infty d x e^{-x} L_n^{|m|}(x) L_s(x) L_k(x) x^{|m|}. \f]
  !> yielding the closed analyticalform,
  !> \f[  d_{nks}^{|m|} =\sum_{n_1=0}^{n} \sum_{k_1=0}^k \sum_{s_1=0}^s L_{kk_1}^{-1/2} L_{nn_1}^{|m|-1/2} L_{ss_1}^{-1/2} (n_1 + k_1 + s_1 + |m|)!, \f]
  !
  !> @param[in] n
  !> @param[in] k 
  !> @param[in] s 
  !> @param[in] m 
  !> 
  !> @param[out] XOUT Output result
  !-----------------------------------------------------------
  FUNCTION ALLX2L(n,k,s,m) RESULT(XOUT)
    !
    ! Computes the associated-Laguerre/Laguerre/x to Laguerre coefficient, i.e. 
    !
    !    L_n^m L_k x^m = \sum_{s=0}^{n+m+k}d_{nks}^m L_s
    !
    ! 
    ! 
    IMPLICIT NONE
    
    !
    INTEGER, intent(in) :: n,k,s,m  ! ... degree and orders of Laguerre polynomials
    TYPE(FM) :: XOUT ! ... coefficient values
    TYPE(FM),SAVE:: AUXN,AUXK,AUXS
    INTEGER :: n1,k1,s1
    !
    CALL FM_ENTER_USER_FUNCTION(XOUT)
    !
    !              Compute coefficients
    !
    XOUT = TO_FM('0.0') 
    !
    DO n1=0,n 
       AUXN =Lpjl(REAL(m,dp)-0.5_dp,REAL(n,dp),REAL(n1,dp))
       DO k1=0,k
          AUXK =Lpjl(-0.5_dp,REAL(k,dp),REAL(k1,dp))
          DO s1=0,s
             AUXS = Lpjl(-0.5_dp,REAL(s,dp),REAL(s1,dp))
             XOUT = XOUT + FACTORIAL(TO_FM(n1 + k1 + s1 + m))*AUXN*AUXK*AUXS 
          ENDDO
       ENDDO
    ENDDO
    !
    CALL FM_EXIT_USER_FUNCTION(XOUT)
    !
  END FUNCTION ALLX2L


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
       AUXN =Lpjl(REAL(m,dp)-0.5_dp,REAL(n,dp),REAL(n1,dp))
       DO k1=0,k
          AUXK =Lpjl(-0.5_dp,REAL(k,dp),REAL(k1,dp))
          DO s1=0,s
             AUXS = Lpjl(-0.5_dp,REAL(s,dp),REAL(s1,dp))
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
    REAL(dp), intent(in) :: p,j,l
    TYPE(FM) :: XOUT 
    !
    CALL FM_ENTER_USER_FUNCTION(XOUT)

    !            compute coeff
    
    XOUT = TO_FM('0.0')

    XOUT = TO_FM(((-1)**l))*FACTORIAL(TO_FM(2*p + 2*j + 1)/2)/&
          (FACTORIAL(TO_FM( j -l ))*FACTORIAL(TO_FM(2*l + 2*p  + 1)/2)*FACTORIAL(TO_FM(l)))

    CALL FM_EXIT_USER_FUNCTION(XOUT)

  END FUNCTION Lpjl

  !-----------------------------------------------------------
  !> @brief 
  !> Computes the particle species 'a' FLR coefficients
  !>
  !> \f[ \mathcal{A}_{an}^m = \frac{ \sigma_m n!}{(n+|m|)!}\left(\frac{b_a}{2} \right)^{|m|} \mathcal{K}_{n}(b_a) \f]
  !>
  !> @param[in] m Spherical Harmonics order
  !> @param[in] n FLR order
  !> @param[in] ba Thermal Normalized Perpendicular Wavenumber
  !> @param[out] Amn
  !-----------------------------------------------------------
  FUNCTION curlyAmn(m,n,ba) RESULT(XOUT)
    ! compute FLR coeff. for particle a
    !                  CHECKED
    IMPLICIT NONE
    !
    INTEGER,INTENT(in) :: m,n
    TYPE(FM),INTENT(in) :: ba
    ! 
    TYPE(FM):: XOUT
    TYPE(FM), SAVE :: ker
    INTEGER :: sigmam
    !
    !
    CALL FM_ENTER_USER_FUNCTION(XOUT)
    !
    XOUT = TO_FM('0.0')
    !
    sigmam = sign(1,m)**m
    !
    ker = kerneln(ba,n)
    !
    IF( m .eq. 0) THEN ! NO FLR EFFECTS 
       XOUT = ker
    ELSE
       XOUT = sigmam*(FACTORIAL(TO_FM(n))/FACTORIAL(TO_FM(n +ABS(m))))*((ba/2)**ABS(m))*ker
    ENDIF
    !
    !
    CALL FM_EXIT_USER_FUNCTION(XOUT)
    !
  END FUNCTION curlyAmn


  !-----------------------------------------------------------
  !> @brief 
  !> Computes the particle species 'b' FLR coefficients
  !> @param[in] m Spherical Harmonics order
  !> @param[in] n FLR order
  !> @param[in] bb Thermal Normalized Perpendicular Wavenumber
  !> @param[out] Bmn
  !-----------------------------------------------------------
  FUNCTION curlyBmn(m,n,b_) RESULT(XOUT)
    !
    ! compute FLR coeff. for particle b
    !                               CHECKED
    IMPLICIT NONE
    !
    INTEGER,INTENT(in) :: m,n
    TYPE(FM),INTENT(in) :: b_
    ! 
    TYPE(FM):: XOUT
    TYPE(FM),SAVE :: ker
    INTEGER :: sigmam
    !
    CALL FM_ENTER_USER_FUNCTION(XOUT)
    !
    XOUT = TO_FM('0.0')
    !
    sigmam = sign(1,m)**m
    !
    ker = kerneln(b_,n)
    !
    IF(m .eq. 0) THEN ! 
       XOUT = ker 
    ELSE
       XOUT = sigmam*(FACTORIAL(TO_FM(n))/FACTORIAL(TO_FM(n +ABS(m))))*((b_/2)**ABS(m))*ker
    ENDIF
    !
    CALL FM_EXIT_USER_FUNCTION(XOUT)

  END FUNCTION curlyBmn


  !-----------------------------------------------------------
  !> @brief 
  !> Computes the nth-order kernel function
  !> 
  !> \f[ \mathcal{K}_n(x) = \frac{1}{n!} \left(\frac{x}{2} \right)^{2n } e^{- x^2/4} \f]
  !> @param[in] n kernel order
  !> @param[in] b normalized perpendicular wavevector
  !> @param[out] kerneln 
  !-----------------------------------------------------------
  FUNCTION kerneln(b,n) RESULT(XOUT)
    ! compute the nth-order kernel function
    !                       CHECKED
    IMPLICIT NONE 
    !
    INTEGER,INTENT(in):: n
    TYPE(FM), INTENT(in) :: b
    !
    TYPE(FM) :: XOUT
    !
    CALL FM_ENTER_USER_FUNCTION(XOUT)

    XOUT = TO_FM(0)
    !              CHECKED
    
    IF (n .LT. 0) THEN ! ... not defined return 0 
           XOUT = TO_FM('0.0')
    ELSE
       IF( b .eq. TO_FM('0.0') .and. n .eq. 0) THEN ! NO FLR effetcs
          XOUT = TO_FM('1.0')
       ELSE 
          XOUT = ((b/2)**n)*EXP(- (b/2)*(b/2))/FACTORIAL(TO_FM(n))
       ENDIF
    ENDIF
    !
    CALL FM_EXIT_USER_FUNCTION(XOUT)

  END FUNCTION kerneln



  FUNCTION normepm(p,m) RESULT(XOUT)
    ! compute the norm of e^{pm} \cdot e^{pm}
    IMPLICIT NONE
    !
    INTEGER,intent(in) :: p,m
    !
    TYPE(FM) :: XOUT

    CALL FM_ENTER_USER_FUNCTION(XOUT)
    !
    !
    XOUT = TO_FM(0)
    !
    XOUT = 2**(2*p)*FACTORIAL(TO_FM(p))*FACTORIAL( TO_FM(2*p -1)/2)/FACTORIAL(TO_FM(2*p))/SQRT(PIFM)

    !
    CALL FM_EXIT_USER_FUNCTION(XOUT)
    !

  END FUNCTION normepm

  FUNCTION Injlkmp(n,j,l,k,m,p) RESULT(XOUT)
    ! Compute speed integrals in Lorentz operator
    ! Note: Only the zeroth-order harmonics, i.e. consider  m=0 !!
    ! Not defined if p ==0
  IMPLICIT NONE
  !
  INTEGER,INTENT(in) :: n,j,l,k,m,p
  TYPE(FM) :: XOUT
  TYPE(FM),SAVE :: NX1,NX2,T4inv,d0nkf_,L1,L2
  CHARACTER(len=T4len) :: T4str_
  !
  INTEGER :: f,g,h,h1,j1
  CALL FM_ENTER_USER_FUNCTION(XOUT)
  !
  XOUT = TO_FM(0)
  
  floop:DO f=0,n+k
     !                  Associated Laguerre basis transformation
     d0nkf_ = ALL2L(n,k,f,0)
     gloop:DO g=0,l+2*f
        IF( g .eq. p) THEN
           hloop:DO h=0,f+l/2
              !            Inverse Basis transformation
              ! _______________________________________________________________________________________
              ! Previous method: very slow at the execution time
              !  T4inv =  READ_T4_FROM_CSV(l,f,p,h,1) ! ... (T^-1)^gh_lf = ... T^lf_gh
              !
              !_______________________________________________________________________________________
              ! New method: Much faster than previous method 
              !                                     write(*,*) l,f,g,h
              T4str_ =  GET_T4(l,f,g,h) ! ... (T^-1)^gh_lf = ... T^lf_gh

              T4inv = SQRT(PIFM)*(2**l)*FACTORIAL(TO_FM(l))*FACTORIAL(TO_FM(h))*(2*g + 1)/2/FACTORIAL(TO_FM(2*g+2*h +1)/2)*TO_FM(T4str_)
              !________________________________________________________________________________________
              IF( T4inv .ne. TO_FM(0)) THEN 

                 h1loop:DO h1=0,h
                    L1= Lpjl(real(p,dp),real(h,dp),real(h1,dp))
                    j1loop: DO j1=0,j
                       L2 = Lpjl(real(p,dp),real(j,dp),real(j1,dp))
                       NX1 = FACTORIAL(TO_FM(p + h1 + j1 -1))
                       XOUT = XOUT + L2*L1*T4inv*d0nkf_*NX1*2/TO_FM(2*p + 1)/SQRT(PIFM)
                    ENDDO j1loop
                 ENDDO h1loop
              ENDIF
           ENDDO hloop
        ENDIF
     ENDDO gloop
  ENDDO floop
  !
  CALL FM_EXIT_USER_FUNCTION(XOUT)
  !
  !
END FUNCTION Injlkmp 



  !-----------------------------------------------------------
  !> @brief 
  !> Computes the spherical harmonics normalization coefficient
  !> 
  !> \f[ \f]
  !> @param[in] p
  !> @param[in] j
  !> @param[out] Results
  !-----------------------------------------------------------
  FUNCTION sigmapj(p,j) RESULT(XOUT)
    ! compute spherical harmonics normalization coefficient
    !
    IMPLICIT NONE 
    !
    INTEGER,INTENT(in):: p
    INTEGER, INTENT(in) :: j
    !
    TYPE(FM) :: XOUT
    !
    CALL FM_ENTER_USER_FUNCTION(XOUT)

    XOUT = TO_FM(0)
    !
    XOUT = FACTORIAL(TO_FM(p))*FACTORIAL(TO_FM(2*p+2*j+1)/2)/(TO_FM(2)**p)/FACTORIAL(TO_FM(2*p+1)/2)/FACTORIAL(TO_FM(j))
    !
    CALL FM_EXIT_USER_FUNCTION(XOUT)

  END FUNCTION sigmapj
  !
  ! _________________ LORENTZ COEFFICIENTS_________________________

  !-----------------------------------------------------------
  !> @brief 
  !> Computes the drift kinetic \mathcal{A}_{\parallel ei}^{lk0}
  !> 
  !> \f[ \f]
  !> @param[in] p
  !> @param[in] j
  !> @param[out] Results
  !-----------------------------------------------------------
  FUNCTION curlyAlk0pareiDK(l,k) RESULT(XOUT)
    IMPLICIT NONE 
    !
    INTEGER,INTENT(in):: l,k
    !
    TYPE(FM) :: XOUT
    !
    ! LOCAL VARIABLES
    INTEGER :: h
    TYPE(FM),SAVE ::  NX0,NX1,NX2,T4inv
    CHARACTER(len=T4len) :: T4invstr_

    
    CALL FM_ENTER_USER_FUNCTION(XOUT)
    !
    XOUT = TO_FM(0)
    !
    NX0 = 1/SQRT(FACTORIAL(TO_FM(l))*TO_FM(2)**l)
    !
    IF(l .ge. 1 .or. k .ge. 1) THEN 
    DO h=0, k + l/2
       !          INVERSE BASIS TRANSFORMATION
       T4invstr_ =  GET_T4(1,h,l,k) ! ... (T^-1)^pt_lk = ... T^lk_pt
       T4inv =TO_FM(T4invstr_)*SQRT(PIFM)*(TO_FM(2)**l)*FACTORIAL(TO_FM(l))*FACTORIAL(TO_FM(h))*TO_FM(2 + 1)/2/FACTORIAL(TO_FM(2+2*h +1)/2)

    XOUT = XOUT + T4inv*NX0*FACTORIAL(TO_FM(2*h + 1)/2)/FACTORIAL(TO_FM(h))
    !
    ENDDO
    ENDIF
    !
    CALL FM_EXIT_USER_FUNCTION(XOUT)
    !
  END FUNCTION curlyAlk0pareiDK



  !-----------------------------------------------------------
  !> @brief 
  !> Computes the gyrokinetic \mathcal{A}_{\parallel ei}^{lk0}
  !> 
  !> \f[ \f]
  !> @param[in] p
  !> @param[in] j
  !> @param[out] Results
  !-----------------------------------------------------------
  FUNCTION curlyAlk0pareiGK(l,k) RESULT(XOUT)
    IMPLICIT NONE 
    !
    INTEGER,INTENT(in):: l,k
    !
    TYPE(FM) :: XOUT
    !
    ! LOCAL VARIABLES
    INTEGER :: h,g,f,n
    TYPE(FM),SAVE ::  NX0,T4inv,ker_,d0nkf_
    CHARACTER(len=T4len) :: T4invstr_

    
    CALL FM_ENTER_USER_FUNCTION(XOUT)
    !
    XOUT = TO_FM(0)
    !
    NX0 = 1/SQRT(FACTORIAL(TO_FM(l))*TO_FM(2)**l)
    !
    nloop: DO n=0, nimaxxFLR
       ker_ = kerneln(bbfm_,n)
       floop: DO f=0, n + k
          d0nkf_ = ALLX2L(n,k,f,0)
          gloop: DO g=0, l+2*f
             IF(g .eq. 1) THEN 
              hloop:  DO h=0, f + l/2
                   !          INVERSE BASIS TRANSFORMATION
                   T4invstr_ =  GET_T4(1,h,l,f) ! ... (T^-1)^pt_lk = ... T^lk_pt
                   T4inv =TO_FM(T4invstr_)*SQRT(PIFM)*(TO_FM(2)**l)*FACTORIAL(TO_FM(l))*FACTORIAL(TO_FM(h))*TO_FM(2*1 + 1)/2/FACTORIAL(TO_FM(2*1+2*h +1)/2)
                   
                   XOUT = XOUT + T4inv*NX0*d0nkf_*ker_*FACTORIAL(TO_FM(2*h + 1)/2)/FACTORIAL(TO_FM(h))
                   !
                ENDDO hloop
             ENDIF
          ENDDO gloop
       ENDDO floop
    ENDDO nloop
    !
    CALL FM_EXIT_USER_FUNCTION(XOUT)
    !
  END FUNCTION curlyAlk0pareiGK

  !-----------------------------------------------------------
  !> @brief 
  !> Computes the gyrokinetic \mathcal{A}_{\perp ei}^{lk0}
  !> 
  !> \f[ \f]
  !> @param[in] p
  !> @param[in] j
  !> @param[out] Results
  !-----------------------------------------------------------
  FUNCTION curlyAlk0perpeiGK(l,k) RESULT(XOUT)
    IMPLICIT NONE 
    !
    INTEGER,INTENT(in):: l,k
    !
    TYPE(FM) :: XOUT
    !
    ! LOCAL VARIABLES
    INTEGER :: h,g,f,n,h1
    TYPE(FM),SAVE ::  NX0,T4inv,ker_,L0hh1,d1nkf_
    CHARACTER(len=T4len) :: T4invstr_

    
    CALL FM_ENTER_USER_FUNCTION(XOUT)
    !
    XOUT = TO_FM(0)
    !
    NX0 = 1/SQRT(FACTORIAL(TO_FM(l))*TO_FM(2)**l)
    !
    nloop: DO n=0, nimaxxFLR
       ker_ = kerneln(bbfm_,n)
       floop: DO f=0, n + k +1
          d1nkf_ = ALLX2L(n,k,f,1)
              hloop:  DO h=0, f + l/2
                   !          INVERSE BASIS TRANSFORMATION
                   T4invstr_ =  GET_T4(0,h,l,f) ! ... (T^-1)^pt_lk = ... T^lk_pt
                   T4inv =TO_FM(T4invstr_)*SQRT(PIFM)*(TO_FM(2)**l)*FACTORIAL(TO_FM(l))*FACTORIAL(TO_FM(h))*TO_FM( 1)/2/FACTORIAL(TO_FM(2*h +1)/2)
                 h1loop: DO h1=1,h
                   L0hh1  =Lpjl(real(0,dp),real(h,dp),real(h1,dp))                  
                   XOUT = XOUT + T4inv*NX0*d1nkf_*ker_*bbfm_*L0hh1/2/(n+1)*FACTORIAL(TO_FM(h1-1))
                   !
                   ENDDO h1loop
                ENDDO hloop
       ENDDO floop
    ENDDO nloop
    !
    CALL FM_EXIT_USER_FUNCTION(XOUT)
    !
  END FUNCTION curlyAlk0perpeiGK
  !
  !__________________________________________________________________
  !
  FUNCTION X2prodLegendre(g,h) RESULT(XOUT)
  !
  ! compute \int d x P_g(x) P_h(x) x^2 = d_h^g
  !
  ! Ref : Arfkten G. B. and Weber H. J., Mathematical methods for physicist, p. 806
  !
  IMPLICIT NONE
  !
  INTEGER,INTENT(in) :: g,h
  TYPE(FM) :: XOUT
  !
  CALL FM_ENTER_USER_FUNCTION(XOUT)
  !
  XOUT = TO_FM(0)
  !
  IF(g .eq. h) THEN 
     XOUT = 2*TO_FM((2*TO_FM(h)**2 + 2*h -1))/(2*h-1)/(2*h +1 )/(2*h +3)
  ELSEIF (g .eq. h -2) THEN 
     XOUT  = TO_FM(2*h*(h-1))/(2*h -3)/(2*h-1)/(2*h+1)
  ELSEIF( g .eq. h +2) THEN
     XOUT = TO_FM(2*(h+1)*(h+2))/(2*h +1)/(2*h + 3)/(2*h +5)
  ENDIF
  !
  CALL FM_EXIT_USER_FUNCTION(XOUT)
  !
  END FUNCTION X2prodLegendre
  !
  !__________________________________________________________________
  !
  FUNCTION CJ(j) RESULT(XOUT)
  ! Associated Laguerre Normalization factor 
  IMPLICIT NONE 
  !
  INTEGER,INTENT(in) :: j
  TYPE(FM) :: XOUT
  !
  CALL FM_ENTER_USER_FUNCTION(XOUT)
  !
  XOUT = FACTORIAL( TO_FM(2*j+3))*3*TO_FM(2)**(3*j +3)*FACTORIAL(TO_FM(j))/FACTORIAL(TO_FM(4*j +6))
  !
  CALL FM_EXIT_USER_FUNCTION(XOUT)
  !  
  END FUNCTION
  !____________________________________________________
  !
  FUNCTION BinomialFM( n, k ) result(XOUT )
 ! Compute the binomial coefficient for positive and negative integers
            implicit none
            integer, intent(in) :: n, k
            TYPE(FM) :: XOUT
            type(FM), save :: f1, f2, f3
            ! See http://mathworld.wolfram.com/BinomialCoefficient.html

            call FM_ENTER_USER_FUNCTION(XOUT)

            XOUT = TO_FM(0)
            if( k>=0 .and. n>=k ) then
              f1 = FACTORIAL(TO_FM(n))
              f2 = FACTORIAL(TO_FM(k))
              f3 = FACTORIAL(TO_FM(n-k))
              XOUT = f1/(f2*f3)
            else if( n<0 ) then
              if( k>=0 ) then
                f1 = FACTORIAL(TO_FM( -n+k-1 ))
                f2 = FACTORIAL( TO_FM(k ))
                f3 = FACTORIAL( TO_FM(-n-1 ))
                XOUT = TO_FM(-1)**k * f1 / ( f2*f3 )
              else if( k<=n ) then
                f1 = FACTORIAL( TO_FM(-k-1 ))
                f2 = FACTORIAL( TO_FM(n-k ))
                f3 = FACTORIAL( -n-1 )
                XOUT = TO_FM(-1)**(n-k)*f1/(f2*f3)
              end if
            end if
            call FM_EXIT_USER_FUNCTION(XOUT)
          END FUNCTION BinomialFM
!___________________________________________________
END MODULE coeff























