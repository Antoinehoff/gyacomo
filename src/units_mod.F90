MODULE units

    USE prec_const
  
    IMPLICIT NONE
    PRIVATE
    !! Define the reference values for transforming simulation data into physical units
    !  The default value, commented here, are corresponding to the #186473 DIII-D discharge at rho=0.95.
    REAL(xp), PUBLIC, PROTECTED :: n_ref = -1!25      ! density     (electrons) [1e19m-3]
    REAL(xp), PUBLIC, PROTECTED :: T_ref = -1!0.3     ! temperature (electrons) [keV]
    REAL(xp), PUBLIC, PROTECTED :: R_ref = -1!1.7     ! length   (major radius) [m]     
    REAL(xp), PUBLIC, PROTECTED :: B_ref = -1!2.0     ! magnetic field          [T]
    REAL(xp), PUBLIC, PROTECTED :: m_ref = -1!1       ! mass                    [nucleid mass]
    REAL(xp), PUBLIC, PROTECTED :: q_ref = -1!1       ! charge                  [elementary charge]
    ! Input flag to write heat flux in MW in the std terminal
    LOGICAL,  PUBLIC, PROTECTED :: WRITE_MW = .false.
    ! Physical constants
    REAL(xp), PUBLIC, PROTECTED :: kB    = 1.4e-23 ! Boltzmann constant  [J/K]
    REAL(xp), PUBLIC, PROTECTED :: mp    = 1.7e-27 ! proton mass at rest [kg]
    REAL(xp), PUBLIC, PROTECTED :: ec    = 1.6e-19 ! elementary charge   [C]

    ! Derived normalization quantities
    REAL(xp), PUBLIC, PROTECTED :: c_s0     ! sound speed             [m/s]
    REAL(xp), PUBLIC, PROTECTED :: r_s0     ! sound Larmor radius     [m/s]
    REAL(xp), PUBLIC, PROTECTED :: w_ci     ! ion cyclotron frequency [1/s]
    REAL(xp), PUBLIC, PROTECTED :: A_fs     ! flux surface area (torus approx) [m2]

    ! Physical unit converters
    REAL(xp), PUBLIC, PROTECTED :: hf_ref   ! heat flux     [W/m2]
    REAL(xp), PUBLIC, PROTECTED :: pf_ref   ! particle flux [1/s/m2]
    REAL(xp), PUBLIC, PROTECTED :: mf_ref   ! momentum flux [kg (m/s)2]
    REAL(xp), PUBLIC, PROTECTED :: ph_ref   ! voltage       [V]
    REAL(xp), PUBLIC, PROTECTED :: pow_ref  ! power         [MW]

    PUBLIC :: units_readinputs, units_outputinputs

  CONTAINS
  
    SUBROUTINE units_readinputs
      ! Read the input parameters
      USE basic,      ONLY: lu_in
      USE geometry,   ONLY: eps, geom
      USE prec_const, ONLY: xp, pi
      IMPLICIT NONE
      REAL(xp) :: T0,n0, p0, m0, q0
      INTEGER  :: fu
      ! Define namelist to read
      NAMELIST /UNITS/ n_ref, T_ref, R_ref, B_ref, m_ref, q_ref, WRITE_MW
      ! Read namelist in input file
      READ(unit=lu_in,nml=units,iostat=fu)
      ! Convert in mksa
      T0      = ec*(1000*T_ref)/kB                ! temperature [K]
      n0      = n_ref*1e19                        ! density     [m-3]
      m0      = m_ref*mp                          ! mass        [kg]
      q0      = q_ref*ec                          ! charge      [C]
      ! Setup derived normalization quantities
      c_s0  = sqrt(kB*T0/m0)                      ! Sound speed                       [m/s]
      w_ci  = q0*B_ref/m0                         ! Ion cyclotron frequ.              [Hz]
      r_s0  = c_s0/w_ci                           ! Ion Larmor radius                 [m]
      SELECT CASE(geom)
        CASE('z-pinch','Z-pinch','zpinch','Zpinch')
            A_fs  = 2*pi*R_ref*R_ref              ! Area of a cylinder of height R_ref [m2]
        CASE DEFAULT
            A_fs  = 4*pi**2*eps*R_ref**2          ! Area of the toroidal flux surface [m2]
      END SELECT
      ! Setup physical unit converters
      p0      = n0*kB*T0                          ! pressure [Pa] 
      hf_ref  = p0*c_s0*r_s0**2/R_ref**2          ! Gyrobohm heat flux [W/m2]
      pow_ref = hf_ref*A_fs/1e6                   ! Gyrobohm power [MW]
      pf_ref  = n0*c_s0*r_s0**2/R_ref**2          ! Particle flux [1/s/m2]
      mf_ref  = n0*m0*c_s0**2*r_s0**2/R_ref**2    ! Momentum flux [kg (m/s)2]
    END SUBROUTINE units_readinputs
  
  
    SUBROUTINE units_outputinputs(fid)
      ! Write the input parameters to the results_xx.h5 file
      USE futils, ONLY: attach, creatd
      IMPLICIT NONE
      INTEGER, INTENT(in) :: fid
      CHARACTER(len=256)  :: str
      WRITE(str,'(a)') '/data/input/units'
      CALL creatd(fid, 0,(/0/),TRIM(str),'Units Input')
      CALL attach(fid, TRIM(str), "density          (electrons)",   n_ref)
      CALL attach(fid, TRIM(str), "temperature      (electrons)",   T_ref)
      CALL attach(fid, TRIM(str), "length           (major radius)",R_ref)
      CALL attach(fid, TRIM(str), "magnetic field   (electrons)",   B_ref)
      CALL attach(fid, TRIM(str), "mass             (proton)",      m_ref)
      CALL attach(fid, TRIM(str), "charge           (elementary)",  q_ref)
      ! Write the mksa converter if units input are provided
      if ( (hf_ref > 0) .AND. (pf_ref > 0) .AND. (mf_ref > 0)) THEN
        CALL attach(fid, TRIM(str), "heat flux        (gyroBohm)",    hf_ref)
        CALL attach(fid, TRIM(str), "power            (gyroBohm)",    pow_ref)
        CALL attach(fid, TRIM(str), "particle flux    (gyroBohm)",    pf_ref)
        CALL attach(fid, TRIM(str), "momentum flux    (gyroBohm)",    mf_ref)
      ENDIF
    END SUBROUTINE units_outputinputs

END MODULE units