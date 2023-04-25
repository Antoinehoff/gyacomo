MODULE species
  ! Module for diagnostic parameters
  USE prec_const
  IMPLICIT NONE
  PRIVATE
  !! Input parameters
  CHARACTER(len=32) :: name_               ! name of the species
  REAL(xp)          :: tau_                ! Temperature
  REAL(xp)          :: sigma_              ! sqrt mass ratio
  REAL(xp)          :: q_                  ! Charge
  REAL(xp)          :: k_N_                ! density drive (L_ref/L_Ni)
  REAL(xp)          :: k_T_                ! temperature drive (L_ref/L_Ti)
  !! Arrays to store all species features
  CHARACTER(len=32),&
            ALLOCATABLE, DIMENSION(:),  PUBLIC, PROTECTED :: name               ! name of the species
  REAL(xp), ALLOCATABLE, DIMENSION(:),  PUBLIC, PROTECTED :: tau                ! Temperature
  REAL(xp), ALLOCATABLE, DIMENSION(:),  PUBLIC, PROTECTED :: sigma              ! sqrt mass ratio
  REAL(xp), ALLOCATABLE, DIMENSION(:),  PUBLIC, PROTECTED :: q                  ! Charge
  REAL(xp), ALLOCATABLE, DIMENSION(:),  PUBLIC, PROTECTED :: k_N                ! density drive (L_ref/L_Ni)
  REAL(xp), ALLOCATABLE, DIMENSION(:),  PUBLIC, PROTECTED :: k_T                ! temperature drive (L_ref/L_Ti)
  REAL(xp), ALLOCATABLE, DIMENSION(:,:),PUBLIC, PROTECTED :: nu_ab              ! Collision frequency tensor
  !! Auxiliary variables to store precomputation
  REAL(xp), ALLOCATABLE, DIMENSION(:),PUBLIC, PROTECTED :: tau_q              ! factor of the magnetic moment coupling
  REAL(xp), ALLOCATABLE, DIMENSION(:),PUBLIC, PROTECTED :: q_tau              ! charge/temp ratio
  REAL(xp), ALLOCATABLE, DIMENSION(:),PUBLIC, PROTECTED :: sqrtTau_q          ! factor of parallel moment term
  REAL(xp), ALLOCATABLE, DIMENSION(:),PUBLIC, PROTECTED :: q_sigma_sqrtTau    ! factor of parallel phi term
  REAL(xp), ALLOCATABLE, DIMENSION(:),PUBLIC, PROTECTED :: sigma2_tau_o2      ! factor of the Kernel argument
  REAL(xp), ALLOCATABLE, DIMENSION(:),PUBLIC, PROTECTED :: sqrt_sigma2_tau_o2 ! to avoid multiple SQRT eval
  REAL(xp), ALLOCATABLE, DIMENSION(:),PUBLIC, PROTECTED :: q2_tau             ! factor of the gammaD sum
  REAL(xp), ALLOCATABLE, DIMENSION(:),PUBLIC, PROTECTED :: q_o_sqrt_tau_sigma ! For psi field terms
  REAL(xp), ALLOCATABLE, DIMENSION(:),PUBLIC, PROTECTED :: sqrt_tau_o_sigma   ! For Ampere eq
  REAL(xp), ALLOCATABLE, DIMENSION(:),PUBLIC, PROTECTED :: xpdx               ! radial pressure gradient
  !! Accessible routines
  PUBLIC :: species_readinputs, species_outputinputs
CONTAINS

  SUBROUTINE species_readinputs
    !    Read the input parameters
    USE basic, ONLY : lu_in
    USE model, ONLY : Na, nu, ADIAB_E
    USE prec_const
    IMPLICIT NONE
    INTEGER :: ia,ib
    ! expected namelist in the input file
    NAMELIST /SPECIES/ &
      name_, tau_, sigma_, q_, k_N_, k_T_
    ! allocate the arrays of species parameters
    CALL species_allocate
    ! loop over the species namelists in the input file
    DO ia = 1,Na
      ! default parameters
      name_  = 'ions'
      tau_   = 1._xp
      sigma_ = 1._xp
      q_     = 1._xp
      k_N_   = 2.22_xp
      k_T_   = 6.96_xp
      ! read input
      READ(lu_in,species)
      ! place values found in the arrays
      name(ia)               = name_
      tau(ia)                = tau_
      sigma(ia)              = sigma_
      q(ia)                  = q_
      k_N(ia)                = k_N_
      k_T(ia)                = k_T_
      tau_q(ia)              = tau_/q_
      ! precompute factors
      q_tau(ia)              = q_/tau_
      sqrtTau_q(ia)          = sqrt(tau_)/q_
      q_sigma_sqrtTau(ia)    = q_/sigma_/SQRT(tau_)
      sigma2_tau_o2(ia)      = sigma_**2 * tau_/2._xp
      sqrt_sigma2_tau_o2(ia) = SQRT(sigma_**2 * tau_/2._xp)
      q2_tau(ia)             = (q_**2)/tau_
      q_o_sqrt_tau_sigma(ia) = q_/SQRT(tau_)/sigma_
      sqrt_tau_o_sigma(ia)   = SQRT(tau_)/sigma_
      xpdx(ia)               = 0._xp !not implemented yet
      ! We remove the adiabatic electron flag if electrons are included
      SELECT CASE (name_)
      CASE ('electrons','e','electron')
        ADIAB_E = .FALSE.
      END SELECT
    ENDDO
    !! Set collision frequency tensor
    IF (nu .EQ. 0) THEN
      nu_ab = 0
    ELSE
      DO ia = 1,Na
        DO ib = 1,Na
          !! We use the ion-ion collision as normalization with definition
          !   nu_ii = 4 sqrt(pi)/3 T_i^(-3/2) m_i^(-1/2) q^4 n_i0 ln(Lambda)
          SELECT CASE (name(ia))
          CASE ('electrons','e','electron') ! e-e and e-i collision
            nu_ab(ia,ib) = nu/sigma(ia) * (tau(ia))**(3._xp/2._xp)  ! (where already multiplied by 0.532)
          CASE ('ions','ion','i') ! i-e and i-i collision
            nu_ab(ia,ib) = nu
          CASE DEFAULT
            ERROR STOP "!! No collision model for these species interactions"
          END SELECT
          ! I think we can just write
          ! nu_ab(ia,ib) = nu/sigma(ia) * (tau(ia))**(3._xp/2._xp)
        ENDDO
      ENDDO
    ENDIF
    ! nu_e            = nu/sigma_e * (tau_e)**(3._xp/2._xp)  ! electron-ion collision frequency (where already multiplied by 0.532)
    ! nu_i            = nu ! ion-ion collision frequ.
    ! nu_ee           = nu_e ! e-e coll. frequ.
    ! nu_ie           = nu_i ! i-e coll. frequ.

  END SUBROUTINE species_readinputs


  SUBROUTINE species_outputinputs(fid)
    !    Write the input parameters to the results_xx.h5 file
    USE futils, ONLY: attach, creatd
    USE model,  ONLY: Na
    IMPLICIT NONE
    INTEGER, INTENT(in) :: fid
    INTEGER :: ia
    CHARACTER(len=256)  :: str
    DO ia = 1,Na
      WRITE(str,'(a,a)') '/data/input/', name(ia)
      CALL creatd(fid, 0,(/0/),TRIM(str),'Species Input')
      CALL attach(fid, TRIM(str),    "name",  name(ia))
      CALL attach(fid, TRIM(str),     "tau",   tau(ia))
      CALL attach(fid, TRIM(str),   "sigma", sigma(ia))
      CALL attach(fid, TRIM(str),       "q",     q(ia))
      CALL attach(fid, TRIM(str),      "k_N",  k_N(ia))
      CALL attach(fid, TRIM(str),      "k_T",  k_T(ia))
    ENDDO
  END SUBROUTINE species_outputinputs

  SUBROUTINE species_allocate
    USE model, ONLY : Na
    IMPLICIT NONE
    !! Allocate the arrays
      ALLOCATE(              name(1:Na))
      ALLOCATE(        nu_ab(1:Na,1:Na))
      ALLOCATE(               tau(1:Na))
      ALLOCATE(             sigma(1:Na))
      ALLOCATE(                 q(1:Na))
      ALLOCATE(               k_N(1:Na))
      ALLOCATE(               k_T(1:Na))
      ALLOCATE(             tau_q(1:Na))
      ALLOCATE(             q_tau(1:Na))
      ALLOCATE(         sqrtTau_q(1:Na))
      ALLOCATE(   q_sigma_sqrtTau(1:Na))
      ALLOCATE(     sigma2_tau_o2(1:Na))
      ALLOCATE(sqrt_sigma2_tau_o2(1:Na))
      ALLOCATE(            q2_tau(1:Na))
      ALLOCATE(q_o_sqrt_tau_sigma(1:Na))
      ALLOCATE(  sqrt_tau_o_sigma(1:Na))
      ALLOCATE(              xpdx(1:Na))
  END SUBROUTINE species_allocate

END MODULE species
