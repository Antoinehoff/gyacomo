module collision
! contains the Hermite-Laguerre collision operators without the use of COSOLVER
USE prec_const, ONLY : xp
IMPLICIT NONE
PRIVATE
CHARACTER(len=32),    PUBLIC, PROTECTED :: collision_model      ! (Lenard-Bernstein: 'LB', Dougherty: 'DG', Sugama: 'SG', Lorentz: 'LR', Landau: 'LD')
LOGICAL,              PUBLIC, PROTECTED :: GK_CO        =.true. ! activates GK effects on CO
LOGICAL,              PUBLIC, PROTECTED :: INTERSPECIES =.true. ! activates interpecies collision
CHARACTER(len=128),   PUBLIC, PROTECTED :: mat_file    ! COSOlver matrix file names
REAL(xp),             PUBLIC, PROTECTED :: collision_kcut = 100.0
LOGICAL,              PUBLIC, PROTECTED :: cosolver_coll ! check if cosolver matrices are used

PUBLIC :: init_collision
PUBLIC :: collision_readinputs, coll_outputinputs
PUBLIC :: compute_Capj

CONTAINS
  SUBROUTINE init_collision
    USE cosolver_interface, ONLY: load_COSOlver_mat
    IMPLICIT NONE
    ! Load the COSOlver collision operator coefficients
    IF(cosolver_coll) &
    CALL load_COSOlver_mat(GK_CO,INTERSPECIES,mat_file,collision_kcut)
  END SUBROUTINE init_collision

  SUBROUTINE collision_readinputs
    ! Read the input parameters
    USE basic, ONLY: lu_in
    IMPLICIT NONE
    NAMELIST /COLLISION/ collision_model, GK_CO, INTERSPECIES, mat_file, collision_kcut
    READ(lu_in,collision)
    SELECT CASE(collision_model)
      ! Lenhard bernstein
      CASE   ('LB','lenhard-bernstein','Lenhard-Bernstein')
        collision_model = 'LB'
        cosolver_coll   = .false.
        INTERSPECIES    = .false.
      ! Ivanov model (Ivanov et al. 2022)
      CASE ('IV','ivanov','Ivanov')
        collision_model = 'IV'
        cosolver_coll   = .false.
        INTERSPECIES    = .false.
      ! Dougherty
      CASE ('DG','dougherty','Dougherty')
        collision_model = 'DG'
        cosolver_coll   = .false.
        INTERSPECIES    = .false.
      ! Sugama
      CASE ('SG','sugama','Sugama')
        collision_model = 'SG'
        cosolver_coll   = .true.
      ! Lorentz (Pitch angle)
      CASE ('LR','lorentz','Lorentz','PA','pitch-angle')
        collision_model = 'LR'
        cosolver_coll   = .true.
        INTERSPECIES    = .false.
      ! Landau named also Coulomb or Fokker-Planck
      CASE ('LD','landau','Landau','CL','coulomb','Coulomb','FP','fokker-planck','Fokker-Planck')
        collision_model = 'LD'
        cosolver_coll   = .true.
      CASE ('none','hypcoll','dvpar4')
        collision_model = 'NO'
        cosolver_coll   = .false.
        INTERSPECIES    = .false.
      CASE DEFAULT
        ERROR STOP '>> ERROR << collision model not recognized!!'
    END SELECT

  END SUBROUTINE collision_readinputs

  SUBROUTINE coll_outputinputs(fidres)
    !    Write the input parameters to the results_xx.h5 file
    USE futils, ONLY: attach, creatd
    IMPLICIT NONE
    INTEGER, INTENT(in) :: fidres
    CHARACTER(len=256) :: str
    CHARACTER(len=2)   :: gkco = 'DK'
    CHARACTER(len=2)   :: abco = 'aa'
    CHARACTER(len=6)   :: coname
    IF (GK_CO)   gkco = 'GK'
    IF (INTERSPECIES) abco = 'ab'
    WRITE(coname,'(a2,a2,a2)') collision_model,gkco,abco
    WRITE(str,'(a)') '/data/input/coll'
    CALL creatd(fidres, 0,(/0/),TRIM(str),'Collision Input')
    CALL attach(fidres, TRIM(str),         "CO", coname)
    CALL attach(fidres, TRIM(str), "matfilename",mat_file)
  END SUBROUTINE coll_outputinputs

  SUBROUTINE compute_Capj
    USE array, ONLY: Capj
    USE model, ONLY: nu
    USE cosolver_interface, ONLY: compute_cosolver_coll
    IMPLICIT NONE
    IF (nu .NE. 0) THEN
      SELECT CASE(collision_model)
        CASE ('LB')
          CALL compute_lenard_bernstein
        CASE('IV')
          CALL compute_ivanov
        CASE ('DG')
          CALL compute_dougherty
        CASE ('SG','LR','LD')
          CALL compute_cosolver_coll(GK_CO)
        CASE ('none','hypcoll','dvpar4')
          Capj = 0._xp
        CASE DEFAULT
          ERROR STOP '>> ERROR << collision operator not recognized!!'
      END SELECT
    ELSE
      Capj = 0._xp
    ENDIF
  END SUBROUTINE compute_Capj

  !******************************************************************************!
  !! Lenard Bernstein collision operator
  !******************************************************************************!
  SUBROUTINE compute_lenard_bernstein
    USE grid, ONLY: ias,iae, ips,ipe, ijs,ije, parray, jarray, &
                    ikys,ikye, ikxs,ikxe, izs,ize, kp2array, ngp, ngj, ngz
    USE species,          ONLY: sigma2_tau_o2, nu_ab
    USE time_integration, ONLY: updatetlevel
    USE array,            ONLY: Capj
    USE fields,           ONLY: moments
    IMPLICIT NONE
    COMPLEX(xp) :: TColl_
    REAL(xp)    :: j_xp, p_xp, ba_2, kp2
    INTEGER     :: iz,ikx,iky,ij,ip,ia,eo,ipi,iji,izi
    DO iz = izs,ize; izi = iz + ngz/2
    DO ikx = ikxs, ikxe;
    DO iky = ikys, ikye;
    DO ij = ijs,ije; iji = ij + ngj/2
    DO ip = ips,ipe; ipi = ip + ngp/2
    DO ia = ias, iae
      !** Auxiliary variables **
      eo   = MODULO(parray(ipi),2)
      p_xp = REAL(parray(ipi),xp)
      j_xp = REAL(jarray(iji),xp)
      kp2  = kp2array(iky,ikx,izi,eo)
      ba_2 = kp2 * sigma2_tau_o2(ia) ! this is (ba/2)^2

      !** Assembling collison operator **
      ! -nuee (p + 2j) Nepj
      TColl_ = -nu_ab(ia,ia) * (p_xp + 2._xp*j_xp)*moments(ia,ipi,iji,iky,ikx,izi,updatetlevel)
      IF(GK_CO) THEN
        TColl_ = TColl_ - nu_ab(ia,ia) *2._xp*ba_2*moments(ia,ipi,iji,iky,ikx,izi,updatetlevel)
      ENDIF
      Capj(ia,ip,ij,iky,ikx,iz) = TColl_
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    END SUBROUTINE compute_lenard_bernstein

      !******************************************************************************!
  !! Doughtery drift-kinetic collision operator (species like)
  !******************************************************************************!
  SUBROUTINE compute_ivanov
    USE grid, ONLY: local_na, local_np, local_nj, parray, jarray, ngp, ngj, &
                    ip0, ip2, ij0, ij1, ieven, &
                    local_nky, local_nkx, local_nz, ngz,&
                    kp2array
    USE species,          ONLY: nu_ab, tau, q_tau
    USE time_integration, ONLY: updatetlevel
    USE array,            ONLY: Capj
    USE fields,           ONLY: moments, phi
    USE prec_const,       ONLY: xp, SQRT2, twothird
    IMPLICIT NONE
    COMPLEX(xp) :: Tmp
    INTEGER     :: iz,ikx,iky,ij,ip,ia, ipi,iji,izi
    REAL(xp)    :: j_xp, p_xp, kp2_xp
    DO iz = 1,local_nz
      izi = iz + ngz/2
      DO ikx = 1,local_nkx
        DO iky = 1,local_nky 
          kp2_xp = kp2array(iky,ikx,izi,ieven)
          DO ij = 1,local_nj
            iji = ij + ngj/2 
            DO ip = 1,local_np
              ipi = ip + ngp/2
              DO ia = 1,local_na
                !** Auxiliary variables **
                p_xp      = REAL(parray(ipi),xp)
                j_xp      = REAL(jarray(iji),xp)
                !** Assembling collison operator **
                IF( (p_xp .EQ. 0._xp) .AND. (j_xp .EQ. 0._xp)) THEN !Ca00
                  Tmp = tau(ia)**2 * kp2_xp**2*(&
                    67._xp/120._xp      *moments(ia,ip0,ij0,iky,ikx,izi,updatetlevel)&
                   +67._xp*SQRT2/240._xp*moments(ia,ip2,ij0,iky,ikx,izi,updatetlevel)&
                   -67._xp/240._xp      *moments(ia,ip0,ij1,iky,ikx,izi,updatetlevel)&
                   -3._xp/10._xp        *q_tau(ia)*phi(iky,ikx,izi))
                ELSEIF( (p_xp .EQ. 2._xp) .AND. (j_xp .EQ. 0._xp)) THEN ! Ca20
                  Tmp = tau(ia) * kp2_xp*(&
                   -SQRT2*twothird*moments(ia,ip0,ij0,iky,ikx,izi,updatetlevel)&
                   -2._xp*twothird*moments(ia,ip2,ij0,iky,ikx,izi,updatetlevel))
                ELSEIF( (p_xp .EQ. 0._xp) .AND. (j_xp .EQ. 1._xp)) THEN ! Ca01
                  Tmp = tau(ia) * kp2_xp*(&
                    2._xp*twothird*moments(ia,ip0,ij0,iky,ikx,izi,updatetlevel)&
                   -2._xp*twothird*moments(ia,ip0,ij1,iky,ikx,izi,updatetlevel))
                ELSE
                  Tmp = 0._xp
                ENDIF
                Capj(ia,ip,ij,iky,ikx,iz) = nu_ab(ia,ia) * Tmp
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE compute_ivanov

  !******************************************************************************!
  !! Doughtery collision operator
  !******************************************************************************!
  SUBROUTINE compute_dougherty
    IMPLICIT NONE
    IF(GK_CO) THEN
      CALL Dougherty_GK
    ELSE
      CALL Dougherty_DK
    ENDIF
  END SUBROUTINE compute_dougherty

  !******************************************************************************!
  !! Doughtery drift-kinetic collision operator (species like)
  !******************************************************************************!
  SUBROUTINE Dougherty_DK
    USE grid, ONLY: local_na, local_np, local_nj, parray, jarray, ngp, ngj, &
                    ip0, ip1, ip2, ij0, ij1, &
                    local_nky, local_nkx, local_nz, ngz
    USE species,          ONLY: nu_ab
    USE time_integration, ONLY: updatetlevel
    USE array,            ONLY: Capj
    USE fields,           ONLY: moments
    USE prec_const,       ONLY: xp, SQRT2, twothird
    IMPLICIT NONE
    COMPLEX(xp) :: Tmp
    INTEGER     :: iz,ikx,iky,ij,ip,ia, ipi,iji,izi
    REAL(xp)    :: j_xp, p_xp
    DO iz = 1,local_nz
      izi = iz + ngz/2
      DO ikx = 1,local_nkx
        DO iky = 1,local_nky 
          DO ij = 1,local_nj
            iji = ij + ngj/2 
            DO ip = 1,local_np
              ipi = ip + ngp/2
              DO ia = 1,local_na
                !** Auxiliary variables **
                p_xp      = REAL(parray(ipi),xp)
                j_xp      = REAL(jarray(iji),xp)
                !** Assembling collison operator **
                Tmp = -(p_xp + 2._xp*j_xp)*moments(ia,ipi,iji,iky,ikx,izi,updatetlevel)
                IF( (p_xp .EQ. 1._xp) .AND. (j_xp .EQ. 0._xp)) THEN !Ce10
                  Tmp = Tmp +  moments(ia,ip1,ij1,iky,ikx,iz,updatetlevel)
                ELSEIF( (p_xp .EQ. 2._xp) .AND. (j_xp .EQ. 0._xp)) THEN ! Ca20
                  Tmp = Tmp +       twothird*moments(ia,ip2,ij0,iky,ikx,izi,updatetlevel) &
                            - SQRT2*twothird*moments(ia,ip0,ij1,iky,ikx,izi,updatetlevel)
                ELSEIF( (p_xp .EQ. 0._xp) .AND. (j_xp .EQ. 1._xp)) THEN ! Ca01
                  Tmp = Tmp +  2._xp*twothird*moments(ia,ip0,ij1,iky,ikx,izi,updatetlevel) &
                              -SQRT2*twothird*moments(ia,ip2,ij0,iky,ikx,izi,updatetlevel)
                ENDIF
                Capj(ia,ip,ij,iky,ikx,iz) = nu_ab(ia,ia) * Tmp
              ENDDO
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  END SUBROUTINE Dougherty_DK
  !******************************************************************************!
  !! Doughtery gyrokinetic collision operator (species like)
  !******************************************************************************!
  SUBROUTINE Dougherty_GK
    USE grid, ONLY: local_na, local_np, local_nj, parray, jarray, ngp, ngj, &
                    ip0, ip1, ip2, &
                    local_nky, local_nkx, local_nz, ngz, kparray, nzgrid
    USE species,          ONLY: sigma2_tau_o2, sqrt_sigma2_tau_o2, nu_ab
    USE array,            ONLY: Capj, nadiab_moments, kernel
    USE prec_const,       ONLY: xp, SQRT2, twothird
    IMPLICIT NONE
    !! Local variables
    COMPLEX(xp) :: dens_,upar_,uperp_,Tpar_,Tperp_,Temp_
    COMPLEX(xp) :: nadiab_moment_0j, Tmp
    REAL(xp)    :: Knp0, Knp1, Knm1, kp
    REAL(xp)    :: n_xp, j_xp, p_xp, ba, ba_2
    INTEGER     :: iz,ikx,iky,ij,ip,ia,eo,in, ipi,iji,izi,ini
    DO iz = 1,local_nz
      izi = iz + ngz/2
      DO ikx = 1,local_nkx
        DO iky = 1,local_nky 
          DO ij = 1,local_nj
            iji = ij + ngj/2 
            DO ip = 1,local_np
              ipi = ip + ngp/2
              DO ia = 1,local_na
    !** Auxiliary variables **
    p_xp      = REAL(parray(ipi),xp)
    j_xp      = REAL(jarray(iji),xp)
    eo        = MIN(nzgrid,MODULO(parray(ipi),2)+1)
    kp        = kparray(iky,ikx,izi,eo)
    ba_2      = kp**2 * sigma2_tau_o2(ia) ! this is (l_a/2)^2
    ba        = 2_xp*kp * sqrt_sigma2_tau_o2(ia)  ! this is l_a
    !** Assembling collison operator **
    ! Velocity-space diffusion (similar to Lenard Bernstein)
    ! -nu (p + 2j + b^2/2) Napj
    Tmp = -(p_xp + 2._xp*j_xp + 2._xp*ba_2)*nadiab_moments(ia,ipi,iji,iky,ikx,izi)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Non zero term for p = 0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    IF( p_xp .EQ. 0 ) THEN ! Kronecker p0
        !** build required fluid moments **
        dens_  = 0._xp
        upar_  = 0._xp; uperp_ = 0._xp
        Tpar_  = 0._xp; Tperp_ = 0._xp
        DO in = 1,local_nj
          ini = in + ngj/2
          n_xp = REAL(jarray(ini),xp)
          ! Store the kernels for sparing readings
          Knp0 =  kernel(ia,ini  ,iky,ikx,izi,eo)
          Knp1 =  kernel(ia,ini+1,iky,ikx,izi,eo)
          Knm1 =  kernel(ia,ini-1,iky,ikx,izi,eo)
          ! Nonadiabatic moments (only different from moments when p=0)
          nadiab_moment_0j   = nadiab_moments(ia,ip0,ini,iky,ikx,izi)
          ! Density
          dens_  = dens_  + Knp0 * nadiab_moment_0j
          ! Perpendicular velocity
          uperp_ = uperp_ + ba*0.5_xp*(Knp0 - Knm1) * nadiab_moment_0j
          ! Parallel temperature
          Tpar_  = Tpar_  + Knp0 * (SQRT2*nadiab_moments(ia,ip2,ini,iky,ikx,izi) + nadiab_moment_0j)
          ! Perpendicular temperature
          Tperp_ = Tperp_ + ((2._xp*n_xp+1._xp)*Knp0 - (n_xp+1._xp)*Knp1 - n_xp*Knm1)*nadiab_moment_0j
        ENDDO
      Temp_  = (Tpar_ + 2._xp*Tperp_)/3._xp - dens_
      ! Add energy restoring term
      Tmp = Tmp + Temp_* 4._xp *  j_xp          * kernel(ia,iji  ,iky,ikx,izi,eo)
      Tmp = Tmp - Temp_* 2._xp * (j_xp + 1._xp) * kernel(ia,iji+1,iky,ikx,izi,eo)
      Tmp = Tmp - Temp_* 2._xp *  j_xp          * kernel(ia,iji-1,iky,ikx,izi,eo)
      Tmp = Tmp + uperp_*ba* (kernel(ia,iji,iky,ikx,izi,eo) - kernel(ia,iji-1,iky,ikx,izi,eo))
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Non zero term for p = 1 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSEIF( p_xp .eq. 1 ) THEN ! kronecker p1
      !** build required fluid moments **
      upar_  = 0._xp
      DO in = 1,local_nj
        ini = in + ngj/2
        n_xp = REAL(jarray(ini),xp)
        ! Parallel velocity
        upar_  = upar_  + kernel(ia,ini,iky,ikx,izi,eo) * nadiab_moments(ia,ip1,ini,iky,ikx,izi)
      ENDDO
      Tmp = Tmp + upar_*kernel(ia,iji,iky,ikx,izi,eo)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Non zero term for p = 2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ELSEIF( p_xp .eq. 2 ) THEN ! kronecker p2
      !** build required fluid moments **
      dens_  = 0._xp
      upar_  = 0._xp; uperp_ = 0._xp
      Tpar_  = 0._xp; Tperp_ = 0._xp
      DO in = 1,local_nj
        ini = in + ngj/2
        n_xp = REAL(jarray(ini),xp)
        ! Store the kernels for sparing readings
        Knp0 =  kernel(ia,ini  ,iky,ikx,izi,eo)
        Knp1 =  kernel(ia,ini+1,iky,ikx,izi,eo)
        Knm1 =  kernel(ia,ini-1,iky,ikx,izi,eo)
        ! Nonadiabatic moments (only different from moments when p=0)
        nadiab_moment_0j = nadiab_moments(ia,ip0,ini,iky,ikx,izi)
        ! Density
        dens_     = dens_     + Knp0 * nadiab_moment_0j
        ! Parallel temperature
        Tpar_  = Tpar_  + Knp0 * (SQRT2*nadiab_moments(ia,ip2,ini,iky,ikx,izi) + nadiab_moment_0j)
        ! Perpendicular temperature
        Tperp_ = Tperp_ + ((2._xp*n_xp+1._xp)*Knp0 - (n_xp+1._xp)*Knp1 - n_xp*Knm1)*nadiab_moment_0j
      ENDDO
      Temp_  = (Tpar_ + 2._xp*Tperp_)/3._xp - dens_
      Tmp = Tmp + Temp_*SQRT2*kernel(ia,iji,iky,ikx,izi,eo)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ENDIF
    ! Multiply by collision parameter
    Capj(ia,ip,ij,iky,ikx,iz) = nu_ab(ia,ia) * Tmp
    ENDDO;ENDDO;ENDDO
    ENDDO;ENDDO
    ENDDO
  END SUBROUTINE Dougherty_GK

end module collision
