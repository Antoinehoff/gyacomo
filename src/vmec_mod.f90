MODULE vmec
  ! Flux-tube geometry from a VMEC equilibrium.
  !
  ! Ported from stella's geometry_vmec_read_netCDF_file.f90 (M. Landreman 2017,
  ! M. Barnes 2018-2019, H. Thienpondt 2024-2025) and expressed in GYACOMO's
  ! GENE-like metric (gxx..gzz, hatB, Jacobian, dBdx, dBdy, dBdz).
  !
  ! The equilibrium is read from an HDF5 file written by
  ! scripts/python_utilities/wout2h5.py (VMEC variable names, Fourier arrays
  ! stored as (mnmax, ns) like LIBSTELL's read_wout_mod).
  !
  ! Normalisation (stella/GIST): a = Aminor_p, psi_edge = phi(ns)*signgs/(2 pi),
  ! B_ref = 2|psi_edge|/a^2, L_ref = Rmajor_p (Lref_opt='R', default) or a.
  !
  ! Field-aligned coordinates on the surface s0 = torflux (rho0 = sqrt(s0)):
  !   x = a (rho - rho0)                       dx/dpsit = c_x = sgn/(rho0 a B_ref)
  !   y = sigma_y a rho0 (alpha - alpha0)      dy/dalpha = c_y = sigma_y a rho0
  !   z = sigma_z (zeta - zeta_center) nfp/P   dz/dzeta  = c_z = sigma_z nfp/P
  ! with alpha = theta_vmec + Lambda - iota*zeta, sgn = sign(psi_edge),
  ! sigma_y = sgn so that B/B_ref = grad x X grad y exactly (C_xy = 1) and
  ! sigma_z = sign(B^zeta) so that b.grad z > 0. P = nfield_periods, so that
  ! z in [-pi,pi] covers P field periods of the device (nfp of them per turn).
  USE prec_const, ONLY: xp, dp
  IMPLICIT NONE
  PRIVATE
  REAL(dp), PARAMETER :: PI_DP = 3.141592653589793238462643383279502884197_dp
  !----------------------------------------------------------------------------
  ! Equilibrium data (names as in the wout file)
  INTEGER,  PUBLIC, PROTECTED :: vmec_nfp = 0, vmec_ns = 0, vmec_mnmax = 0, vmec_mnmax_nyq = 0, vmec_signgs = 0
  LOGICAL,  PUBLIC, PROTECTED :: vmec_lasym = .FALSE.
  REAL(dp), PUBLIC, PROTECTED :: vmec_a_minor = 0._dp, vmec_R_major = 0._dp
  REAL(dp), PUBLIC, PROTECTED :: vmec_psi_edge = 0._dp ! toroidal flux at the LCFS / 2 pi, signed
  REAL(dp), PUBLIC, PROTECTED :: vmec_B_ref    = 0._dp ! 2|psi_edge|/a^2 [T]
  REAL(dp), PUBLIC, PROTECTED, DIMENSION(:), ALLOCATABLE :: vmec_xm, vmec_xn, vmec_xm_nyq, vmec_xn_nyq
  REAL(dp), PUBLIC, PROTECTED, DIMENSION(:), ALLOCATABLE :: vmec_phi, vmec_phips, vmec_iotas, vmec_iotaf, vmec_presf
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: rmnc, zmns, lmns, bmnc, gmnc, bsupumnc, bsupvmnc, bsubumnc, bsubvmnc, bsubsmns
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: rmns, zmnc, lmnc, bmns, gmns, bsupumns, bsupvmns, bsubumns, bsubvmns, bsubsmnc
  LOGICAL :: file_ready = .FALSE.
  !----------------------------------------------------------------------------
  ! Chosen flux surface
  REAL(dp), PUBLIC, PROTECTED :: vmec_s0 = 0._dp, vmec_rho0 = 0._dp, vmec_ds = 0._dp
  REAL(dp), PUBLIC, PROTECTED :: vmec_iota = 0._dp, vmec_diota_ds = 0._dp, vmec_dpres_ds = 0._dp
  REAL(dp), PUBLIC, PROTECTED :: vmec_q0 = 0._dp, vmec_shat = 0._dp
  INTEGER,  PUBLIC, PROTECTED :: vmec_sign_torflux = 0
  INTEGER  :: index_full(2) = 1, index_half(2) = 1
  REAL(dp) :: weight_full(2) = 0._dp, weight_half(2) = 0._dp
  ! Fourier coefficients interpolated to the surface (and their s-derivatives)
  REAL(dp), DIMENSION(:), ALLOCATABLE :: rmnc_s, zmns_s, lmns_s, drmnc_s, dzmns_s, dlmns_s
  REAL(dp), DIMENSION(:), ALLOCATABLE :: rmns_s, zmnc_s, lmnc_s, drmns_s, dzmnc_s, dlmnc_s
  REAL(dp), DIMENSION(:), ALLOCATABLE :: bmnc_s, gmnc_s, bsupumnc_s, bsupvmnc_s, bsubumnc_s, bsubvmnc_s, bsubsmns_s, dbmnc_s
  REAL(dp), DIMENSION(:), ALLOCATABLE :: bmns_s, gmns_s, bsupumns_s, bsupvmns_s, bsubumns_s, bsubvmns_s, bsubsmnc_s, dbmns_s
  LOGICAL :: surface_ready = .FALSE.
  !----------------------------------------------------------------------------
  ! Flux tube
  REAL(dp), PUBLIC, PROTECTED :: vmec_alpha0 = 0._dp, vmec_zeta_center = 0._dp
  REAL(dp), PUBLIC, PROTECTED :: vmec_nfield_periods = 0._dp ! effective value (nfp if input <= 0)
  REAL(dp), PUBLIC, PROTECTED :: vmec_L_ref = 0._dp, vmec_Lref_over_a = 0._dp ! L_ref [m], L_ref/a
  REAL(dp), PUBLIC, PROTECTED :: vmec_cx = 0._dp, vmec_cy = 0._dp, vmec_cz = 0._dp
  INTEGER,  PUBLIC, PROTECTED :: vmec_sigma_y = 1, vmec_sigma_z = 1
  LOGICAL :: tube_ready = .FALSE.
  !----------------------------------------------------------------------------
  ! VMEC quantities at one point (theta, zeta) of the surface, SI units
  TYPE, PUBLIC :: vmec_point_t
    REAL(dp) :: theta = 0._dp, zeta = 0._dp
    REAL(dp) :: R = 0._dp, dRdt = 0._dp, dRdz = 0._dp, dRds = 0._dp ! R and d/dtheta, d/dzeta, d/ds
    REAL(dp) :: Z = 0._dp, dZdt = 0._dp, dZdz = 0._dp, dZds = 0._dp
    REAL(dp) :: Lt = 0._dp, Lz = 0._dp, Ls = 0._dp                 ! dLambda/dtheta, dzeta, ds
    REAL(dp) :: B = 0._dp, dBdt = 0._dp, dBdz = 0._dp, dBds = 0._dp ! |B| and derivatives at fixed VMEC coordinates
    REAL(dp) :: sqrtg = 0._dp                                       ! VMEC Jacobian (signed)
    REAL(dp) :: Bsup_t = 0._dp, Bsup_z = 0._dp                      ! B^theta, B^zeta
    REAL(dp) :: Bsub_t = 0._dp, Bsub_z = 0._dp, Bsub_s = 0._dp      ! B_theta, B_zeta, B_s
  END TYPE vmec_point_t
  ! GYACOMO geometry at one point z of the flux tube (normalised)
  TYPE, PUBLIC :: vmec_geo_t
    REAL(dp) :: zeta = 0._dp, theta = 0._dp, R = 0._dp, Z = 0._dp
    REAL(dp) :: gxx = 0._dp, gxy = 0._dp, gxz = 0._dp, gyy = 0._dp, gyz = 0._dp, gzz = 0._dp
    REAL(dp) :: hatB = 0._dp, Jacobian = 0._dp, dBdx = 0._dp, dBdy = 0._dp, dBdz = 0._dp
    REAL(dp) :: dxdR = 0._dp, dxdZ = 0._dp
    REAL(dp) :: gradx(3) = 0._dp, grady(3) = 0._dp, gradz(3) = 0._dp ! Cartesian (X,Y,Z): grad x, grad y, L_ref grad z
    REAL(dp) :: Bvec(3) = 0._dp                                     ! B/B_ref from B = grad psit x grad alpha
    REAL(dp) :: Bsup_z = 0._dp                                      ! B^zeta [T/m]
  END TYPE vmec_geo_t
  !----------------------------------------------------------------------------
  PUBLIC :: vmec_read_h5, vmec_select_surface, vmec_setup_tube
  PUBLIC :: vmec_eval_point, vmec_lambda, vmec_theta_of, vmec_geometry_at_z
  PUBLIC :: get_vmec_geometry

CONTAINS

  !============================================================================
  ! Read the equilibrium from the HDF5 file produced by wout2h5.py
  !============================================================================
  SUBROUTINE vmec_read_h5(filename)
    USE futils,   ONLY: openf, closef, getarr, isdataset
    USE parallel, ONLY: comm0, my_id
    USE basic,    ONLY: speak, str
    IMPLICIT NONE
    CHARACTER(len=*), INTENT(IN) :: filename
    INTEGER  :: fid, j, itmp(1)
    REAL(dp) :: rtmp(1), dphi
    LOGICAL  :: file_exists
    CHARACTER(len=32), PARAMETER :: needed(27) = [CHARACTER(len=32) :: &
      'nfp','ns','mnmax','mnmax_nyq','signgs','lasym','Aminor_p','Rmajor_p',&
      'xm','xn','xm_nyq','xn_nyq','phi','phips','iotas','iotaf','presf',&
      'rmnc','zmns','lmns','bmnc','gmnc','bsupumnc','bsupvmnc','bsubumnc','bsubvmnc','bsubsmns']
    INQUIRE(file=TRIM(filename), exist=file_exists)
    IF (.NOT. file_exists) THEN
      IF (my_id .EQ. 0) WRITE(*,*) '>> ERROR << VMEC file not found: '//TRIM(filename)
      ERROR STOP '>> ERROR << VMEC file not found (vmec_file in the GEOMETRY namelist, HDF5 from wout2h5.py)'
    ENDIF
    CALL speak('Reading VMEC equilibrium '//TRIM(filename),2)
    CALL openf(TRIM(filename), fid, 'r', 'D', mpicomm=comm0)
    DO j = 1,SIZE(needed)
      IF (.NOT. isdataset(fid, '/'//TRIM(needed(j)))) THEN
        IF (my_id .EQ. 0) WRITE(*,*) '>> ERROR << dataset /'//TRIM(needed(j))//' missing in '//TRIM(filename)
        ERROR STOP '>> ERROR << incomplete VMEC HDF5 file (regenerate it with wout2h5.py)'
      ENDIF
    ENDDO
    CALL getarr(fid, '/nfp',       itmp); vmec_nfp       = itmp(1)
    CALL getarr(fid, '/ns',        itmp); vmec_ns        = itmp(1)
    CALL getarr(fid, '/mnmax',     itmp); vmec_mnmax     = itmp(1)
    CALL getarr(fid, '/mnmax_nyq', itmp); vmec_mnmax_nyq = itmp(1)
    CALL getarr(fid, '/signgs',    itmp); vmec_signgs    = itmp(1)
    CALL getarr(fid, '/lasym',     itmp); vmec_lasym     = (itmp(1) .NE. 0)
    CALL getarr(fid, '/Aminor_p',  rtmp); vmec_a_minor   = rtmp(1)
    CALL getarr(fid, '/Rmajor_p',  rtmp); vmec_R_major   = rtmp(1)
    IF (ALLOCATED(vmec_xm)) CALL vmec_free_file
    ALLOCATE(vmec_xm(vmec_mnmax), vmec_xn(vmec_mnmax), vmec_xm_nyq(vmec_mnmax_nyq), vmec_xn_nyq(vmec_mnmax_nyq))
    ALLOCATE(vmec_phi(vmec_ns), vmec_phips(vmec_ns), vmec_iotas(vmec_ns), vmec_iotaf(vmec_ns), vmec_presf(vmec_ns))
    CALL getarr(fid, '/xm',     vmec_xm);     CALL getarr(fid, '/xn',     vmec_xn)
    CALL getarr(fid, '/xm_nyq', vmec_xm_nyq); CALL getarr(fid, '/xn_nyq', vmec_xn_nyq)
    CALL getarr(fid, '/phi',    vmec_phi);    CALL getarr(fid, '/phips',  vmec_phips)
    CALL getarr(fid, '/iotas',  vmec_iotas);  CALL getarr(fid, '/iotaf',  vmec_iotaf)
    CALL getarr(fid, '/presf',  vmec_presf)
    ALLOCATE(rmnc(vmec_mnmax,vmec_ns), zmns(vmec_mnmax,vmec_ns), lmns(vmec_mnmax,vmec_ns))
    ALLOCATE(bmnc(vmec_mnmax_nyq,vmec_ns), gmnc(vmec_mnmax_nyq,vmec_ns))
    ALLOCATE(bsupumnc(vmec_mnmax_nyq,vmec_ns), bsupvmnc(vmec_mnmax_nyq,vmec_ns))
    ALLOCATE(bsubumnc(vmec_mnmax_nyq,vmec_ns), bsubvmnc(vmec_mnmax_nyq,vmec_ns), bsubsmns(vmec_mnmax_nyq,vmec_ns))
    CALL getarr(fid, '/rmnc', rmnc); CALL getarr(fid, '/zmns', zmns); CALL getarr(fid, '/lmns', lmns)
    CALL getarr(fid, '/bmnc', bmnc); CALL getarr(fid, '/gmnc', gmnc)
    CALL getarr(fid, '/bsupumnc', bsupumnc); CALL getarr(fid, '/bsupvmnc', bsupvmnc)
    CALL getarr(fid, '/bsubumnc', bsubumnc); CALL getarr(fid, '/bsubvmnc', bsubvmnc)
    CALL getarr(fid, '/bsubsmns', bsubsmns)
    IF (vmec_lasym) THEN
      CALL speak('VMEC equilibrium is not stellarator symmetric (lasym): asymmetric terms included, untested',1)
      ALLOCATE(rmns(vmec_mnmax,vmec_ns), zmnc(vmec_mnmax,vmec_ns), lmnc(vmec_mnmax,vmec_ns))
      ALLOCATE(bmns(vmec_mnmax_nyq,vmec_ns), gmns(vmec_mnmax_nyq,vmec_ns))
      ALLOCATE(bsupumns(vmec_mnmax_nyq,vmec_ns), bsupvmns(vmec_mnmax_nyq,vmec_ns))
      ALLOCATE(bsubumns(vmec_mnmax_nyq,vmec_ns), bsubvmns(vmec_mnmax_nyq,vmec_ns), bsubsmnc(vmec_mnmax_nyq,vmec_ns))
      CALL getarr(fid, '/rmns', rmns); CALL getarr(fid, '/zmnc', zmnc); CALL getarr(fid, '/lmnc', lmnc)
      CALL getarr(fid, '/bmns', bmns); CALL getarr(fid, '/gmns', gmns)
      CALL getarr(fid, '/bsupumns', bsupumns); CALL getarr(fid, '/bsupvmns', bsupvmns)
      CALL getarr(fid, '/bsubumns', bsubumns); CALL getarr(fid, '/bsubvmns', bsubvmns)
      CALL getarr(fid, '/bsubsmnc', bsubsmnc)
    ELSE ! zero-size arrays so that the evaluation routines stay generic
      ALLOCATE(rmns(0,0), zmnc(0,0), lmnc(0,0), bmns(0,0), gmns(0,0))
      ALLOCATE(bsupumns(0,0), bsupvmns(0,0), bsubumns(0,0), bsubvmns(0,0), bsubsmnc(0,0))
    ENDIF
    CALL closef(fid)
    ! Sanity checks (as stella's sanity_checks_vmec)
    IF (MAXVAL(ABS(vmec_xm_nyq)) .LT. 1._dp .AND. MAXVAL(ABS(vmec_xn_nyq)) .LT. 1._dp) THEN
      IF (vmec_mnmax_nyq .EQ. vmec_mnmax) THEN
        CALL speak('xm_nyq/xn_nyq not populated in the VMEC file, using xm/xn',1)
        vmec_xm_nyq = vmec_xm; vmec_xn_nyq = vmec_xn
      ELSE
        ERROR STOP '>> ERROR << VMEC xm_nyq/xn_nyq not populated and mnmax_nyq /= mnmax'
      ENDIF
    ENDIF
    IF (ABS(vmec_phi(1)) .GT. 1e-14_dp) ERROR STOP '>> ERROR << VMEC phi array does not begin with 0'
    dphi = vmec_phi(2) - vmec_phi(1)
    DO j = 3,vmec_ns
      IF (ABS(vmec_phi(j)-vmec_phi(j-1)-dphi) .GT. 1e-11_dp*MAX(1._dp,ABS(dphi))) &
        ERROR STOP '>> ERROR << VMEC phi array is not uniformly spaced'
    ENDDO
    IF (vmec_xm(1) .NE. 0._dp .OR. vmec_xn(1) .NE. 0._dp .OR. vmec_xm_nyq(1) .NE. 0._dp .OR. vmec_xn_nyq(1) .NE. 0._dp) &
      ERROR STOP '>> ERROR << first VMEC mode must be m = n = 0'
    ! Reference quantities (stella/GIST conventions)
    vmec_psi_edge     = vmec_phi(vmec_ns)/(2._dp*PI_DP)*REAL(vmec_signgs,dp)
    vmec_sign_torflux = INT(SIGN(1.1_dp, vmec_psi_edge))
    vmec_B_ref        = 2._dp*ABS(vmec_psi_edge)/(vmec_a_minor*vmec_a_minor)
    vmec_ds           = 1._dp/REAL(vmec_ns-1,dp)
    file_ready    = .TRUE.
    surface_ready = .FALSE.
    tube_ready    = .FALSE.
    CALL speak('VMEC: nfp = '//str(vmec_nfp)//', ns = '//str(vmec_ns)//', mnmax = '//str(vmec_mnmax)//&
               ', a = '//str(REAL(vmec_a_minor,xp))//' m, R = '//str(REAL(vmec_R_major,xp))//&
               ' m, B_ref = '//str(REAL(vmec_B_ref,xp))//' T',2)
  END SUBROUTINE vmec_read_h5

  SUBROUTINE vmec_free_file
    DEALLOCATE(vmec_xm, vmec_xn, vmec_xm_nyq, vmec_xn_nyq, vmec_phi, vmec_phips, vmec_iotas, vmec_iotaf, vmec_presf)
    DEALLOCATE(rmnc, zmns, lmns, bmnc, gmnc, bsupumnc, bsupvmnc, bsubumnc, bsubvmnc, bsubsmns)
    DEALLOCATE(rmns, zmnc, lmnc, bmns, gmns, bsupumns, bsupvmns, bsubumns, bsubvmns, bsubsmnc)
  END SUBROUTINE vmec_free_file

  !============================================================================
  ! Choose the flux surface and interpolate every Fourier coefficient to it
  !============================================================================
  SUBROUTINE vmec_select_surface(torflux, surface_option)
    USE basic, ONLY: speak, str
    IMPLICIT NONE
    REAL(xp), INTENT(IN) :: torflux
    INTEGER,  INTENT(IN) :: surface_option
    REAL(dp), DIMENSION(:), ALLOCATABLE :: s_full, s_half, tmp
    REAL(dp) :: s, dr2min
    INTEGER  :: j, jmin, ns, imn
    IF (.NOT. file_ready) ERROR STOP '>> ERROR << vmec_select_surface called before vmec_read_h5'
    ns = vmec_ns
    ALLOCATE(s_full(ns), s_half(ns-1), tmp(ns))
    s_full = [(REAL(j-1,dp)/REAL(ns-1,dp), j=1,ns)]
    s_half = 0.5_dp*(s_full(1:ns-1) + s_full(2:ns))
    SELECT CASE(surface_option)
    CASE(0) ! exact value requested
      s = REAL(torflux,dp)
    CASE(1) ! nearest half-grid surface
      jmin = 1; dr2min = (s_half(1)-REAL(torflux,dp))**2
      DO j = 2,ns-1
        IF ((s_half(j)-REAL(torflux,dp))**2 .LT. dr2min) THEN
          jmin = j; dr2min = (s_half(j)-REAL(torflux,dp))**2
        ENDIF
      ENDDO
      s = s_half(jmin)
    CASE(2) ! nearest full-grid surface
      jmin = 1; dr2min = (s_full(1)-REAL(torflux,dp))**2
      DO j = 2,ns
        IF ((s_full(j)-REAL(torflux,dp))**2 .LT. dr2min) THEN
          jmin = j; dr2min = (s_full(j)-REAL(torflux,dp))**2
        ENDIF
      ENDDO
      s = s_full(jmin)
    CASE DEFAULT
      ERROR STOP '>> ERROR << surface_option must be 0 (exact), 1 (nearest half grid) or 2 (nearest full grid)'
    END SELECT
    IF (s .GT. 1._dp .OR. s .LT. 0._dp) ERROR STOP '>> ERROR << torflux = psi_t/psi_t,edge must lie in [0,1]'
    vmec_s0   = s
    vmec_rho0 = SQRT(s)
    ! Interpolation weights on the full grid: Q = Q(i1) w1 + Q(i2) w2
    IF (s .EQ. 1._dp) THEN
      index_full = [ns-1, ns]; weight_full(1) = 0._dp
    ELSE
      index_full(1) = FLOOR(s*REAL(ns-1,dp)) + 1
      index_full(2) = index_full(1) + 1
      weight_full(1) = REAL(index_full(1),dp) - s*REAL(ns-1,dp)
    ENDIF
    ! Half grid: element j of a half-mesh array sits at s_half(j-1); element 1 is a dummy
    IF (s .LT. s_half(1)) THEN
      CALL speak('WARNING: torflux is inside the first VMEC half-grid surface, extrapolating towards the axis',1)
      index_half = [2, 3]; weight_half(1) = (s_half(2)-s)/(s_half(2)-s_half(1))
    ELSEIF (s .GT. s_half(ns-1)) THEN
      CALL speak('WARNING: torflux is outside the last VMEC half-grid surface, extrapolating towards the LCFS',1)
      index_half = [ns-1, ns]; weight_half(1) = (s_half(ns-1)-s)/(s_half(ns-1)-s_half(ns-2))
    ELSEIF (s .EQ. s_half(ns-1)) THEN
      index_half = [ns-1, ns]; weight_half(1) = 0._dp
    ELSE
      index_half(1) = FLOOR(s*REAL(ns-1,dp) + 0.5_dp) + 1
      IF (index_half(1) .LT. 2) index_half(1) = 2
      index_half(2) = index_half(1) + 1
      weight_half(1) = REAL(index_half(1),dp) - s*REAL(ns-1,dp) - 0.5_dp
    ENDIF
    weight_full(2) = 1._dp - weight_full(1)
    weight_half(2) = 1._dp - weight_half(1)
    ! Profiles on the surface: iota from the half grid, derivatives of the full-grid
    ! profiles land on the half grid
    vmec_iota = interp_half(vmec_iotas)
    tmp = 0._dp; tmp(2:ns) = (vmec_iotaf(2:ns) - vmec_iotaf(1:ns-1))/vmec_ds
    vmec_diota_ds = interp_half(tmp)
    tmp = 0._dp; tmp(2:ns) = (vmec_presf(2:ns) - vmec_presf(1:ns-1))/vmec_ds
    vmec_dpres_ds = interp_half(tmp)
    vmec_q0   = 1._dp/vmec_iota
    vmec_shat = -2._dp*s/vmec_iota*vmec_diota_ds ! shat = -(r/iota) diota/dr with r = a sqrt(s)
    ! Fourier coefficients on the surface. R, Z and B_s are on the full mesh,
    ! Lambda, B, sqrt g, B^u, B^v, B_u, B_v on the half mesh; the radial
    ! derivative of a full-mesh quantity lives on the half mesh and vice versa.
    IF (ALLOCATED(rmnc_s)) CALL vmec_free_surface
    ALLOCATE(rmnc_s(vmec_mnmax), zmns_s(vmec_mnmax), lmns_s(vmec_mnmax))
    ALLOCATE(drmnc_s(vmec_mnmax), dzmns_s(vmec_mnmax), dlmns_s(vmec_mnmax))
    DO imn = 1,vmec_mnmax
      rmnc_s(imn) = interp_full(rmnc(imn,:))
      zmns_s(imn) = interp_full(zmns(imn,:))
      lmns_s(imn) = interp_half(lmns(imn,:))
      tmp(1) = 0._dp; tmp(2:ns) = (rmnc(imn,2:ns) - rmnc(imn,1:ns-1))/vmec_ds
      drmnc_s(imn) = interp_half(tmp)
      tmp(1) = 0._dp; tmp(2:ns) = (zmns(imn,2:ns) - zmns(imn,1:ns-1))/vmec_ds
      dzmns_s(imn) = interp_half(tmp)
      tmp(2:ns-1) = (lmns(imn,3:ns) - lmns(imn,2:ns-1))/vmec_ds; tmp(1) = tmp(2); tmp(ns) = tmp(ns-1)
      dlmns_s(imn) = interp_full(tmp)
    ENDDO
    ALLOCATE(bmnc_s(vmec_mnmax_nyq), gmnc_s(vmec_mnmax_nyq), bsupumnc_s(vmec_mnmax_nyq), bsupvmnc_s(vmec_mnmax_nyq))
    ALLOCATE(bsubumnc_s(vmec_mnmax_nyq), bsubvmnc_s(vmec_mnmax_nyq), bsubsmns_s(vmec_mnmax_nyq), dbmnc_s(vmec_mnmax_nyq))
    DO imn = 1,vmec_mnmax_nyq
      bmnc_s(imn)     = interp_half(bmnc(imn,:))
      gmnc_s(imn)     = interp_half(gmnc(imn,:))
      bsupumnc_s(imn) = interp_half(bsupumnc(imn,:))
      bsupvmnc_s(imn) = interp_half(bsupvmnc(imn,:))
      bsubumnc_s(imn) = interp_half(bsubumnc(imn,:))
      bsubvmnc_s(imn) = interp_half(bsubvmnc(imn,:))
      bsubsmns_s(imn) = interp_full(bsubsmns(imn,:))
      tmp(2:ns-1) = (bmnc(imn,3:ns) - bmnc(imn,2:ns-1))/vmec_ds; tmp(1) = tmp(2); tmp(ns) = tmp(ns-1)
      dbmnc_s(imn) = interp_full(tmp)
    ENDDO
    IF (vmec_lasym) THEN
      ALLOCATE(rmns_s(vmec_mnmax), zmnc_s(vmec_mnmax), lmnc_s(vmec_mnmax))
      ALLOCATE(drmns_s(vmec_mnmax), dzmnc_s(vmec_mnmax), dlmnc_s(vmec_mnmax))
      DO imn = 1,vmec_mnmax
        rmns_s(imn) = interp_full(rmns(imn,:))
        zmnc_s(imn) = interp_full(zmnc(imn,:))
        lmnc_s(imn) = interp_half(lmnc(imn,:))
        tmp(1) = 0._dp; tmp(2:ns) = (rmns(imn,2:ns) - rmns(imn,1:ns-1))/vmec_ds
        drmns_s(imn) = interp_half(tmp)
        tmp(1) = 0._dp; tmp(2:ns) = (zmnc(imn,2:ns) - zmnc(imn,1:ns-1))/vmec_ds
        dzmnc_s(imn) = interp_half(tmp)
        tmp(2:ns-1) = (lmnc(imn,3:ns) - lmnc(imn,2:ns-1))/vmec_ds; tmp(1) = tmp(2); tmp(ns) = tmp(ns-1)
        dlmnc_s(imn) = interp_full(tmp)
      ENDDO
      ALLOCATE(bmns_s(vmec_mnmax_nyq), gmns_s(vmec_mnmax_nyq), bsupumns_s(vmec_mnmax_nyq), bsupvmns_s(vmec_mnmax_nyq))
      ALLOCATE(bsubumns_s(vmec_mnmax_nyq), bsubvmns_s(vmec_mnmax_nyq), bsubsmnc_s(vmec_mnmax_nyq), dbmns_s(vmec_mnmax_nyq))
      DO imn = 1,vmec_mnmax_nyq
        bmns_s(imn)     = interp_half(bmns(imn,:))
        gmns_s(imn)     = interp_half(gmns(imn,:))
        bsupumns_s(imn) = interp_half(bsupumns(imn,:))
        bsupvmns_s(imn) = interp_half(bsupvmns(imn,:))
        bsubumns_s(imn) = interp_half(bsubumns(imn,:))
        bsubvmns_s(imn) = interp_half(bsubvmns(imn,:))
        bsubsmnc_s(imn) = interp_full(bsubsmnc(imn,:))
        tmp(2:ns-1) = (bmns(imn,3:ns) - bmns(imn,2:ns-1))/vmec_ds; tmp(1) = tmp(2); tmp(ns) = tmp(ns-1)
        dbmns_s(imn) = interp_full(tmp)
      ENDDO
    ELSE
      ALLOCATE(rmns_s(0), zmnc_s(0), lmnc_s(0), drmns_s(0), dzmnc_s(0), dlmnc_s(0))
      ALLOCATE(bmns_s(0), gmns_s(0), bsupumns_s(0), bsupvmns_s(0), bsubumns_s(0), bsubvmns_s(0), bsubsmnc_s(0), dbmns_s(0))
    ENDIF
    DEALLOCATE(s_full, s_half, tmp)
    surface_ready = .TRUE.
    tube_ready    = .FALSE.
    CALL speak('VMEC surface s = '//str(REAL(vmec_s0,xp))//' (rho = '//str(REAL(vmec_rho0,xp))//&
               '): iota = '//str(REAL(vmec_iota,xp))//', q = '//str(REAL(vmec_q0,xp))//&
               ', shat = '//str(REAL(vmec_shat,xp)),2)
  CONTAINS
    REAL(dp) FUNCTION interp_full(q)
      REAL(dp), DIMENSION(:), INTENT(IN) :: q
      interp_full = q(index_full(1))*weight_full(1) + q(index_full(2))*weight_full(2)
    END FUNCTION interp_full
    REAL(dp) FUNCTION interp_half(q)
      REAL(dp), DIMENSION(:), INTENT(IN) :: q
      interp_half = q(index_half(1))*weight_half(1) + q(index_half(2))*weight_half(2)
    END FUNCTION interp_half
  END SUBROUTINE vmec_select_surface

  SUBROUTINE vmec_free_surface
    DEALLOCATE(rmnc_s, zmns_s, lmns_s, drmnc_s, dzmns_s, dlmns_s)
    DEALLOCATE(rmns_s, zmnc_s, lmnc_s, drmns_s, dzmnc_s, dlmnc_s)
    DEALLOCATE(bmnc_s, gmnc_s, bsupumnc_s, bsupvmnc_s, bsubumnc_s, bsubvmnc_s, bsubsmns_s, dbmnc_s)
    DEALLOCATE(bmns_s, gmns_s, bsupumns_s, bsupvmns_s, bsubumns_s, bsubvmns_s, bsubsmnc_s, dbmns_s)
  END SUBROUTINE vmec_free_surface

  !============================================================================
  ! Lambda(theta, zeta) on the surface (for the field-line root solve)
  !============================================================================
  REAL(dp) FUNCTION vmec_lambda(theta, zeta)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: theta, zeta
    REAL(dp) :: angle
    INTEGER  :: imn
    vmec_lambda = 0._dp
    DO imn = 1,vmec_mnmax
      angle = vmec_xm(imn)*theta - vmec_xn(imn)*zeta
      vmec_lambda = vmec_lambda + lmns_s(imn)*SIN(angle)
      IF (vmec_lasym) vmec_lambda = vmec_lambda + lmnc_s(imn)*COS(angle)
    ENDDO
  END FUNCTION vmec_lambda

  !============================================================================
  ! Fourier sums at one point (theta, zeta) of the surface
  !============================================================================
  SUBROUTINE vmec_eval_point(theta, zeta, pt)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: theta, zeta
    TYPE(vmec_point_t), INTENT(OUT) :: pt
    REAL(dp) :: m, n, angle, c, s
    INTEGER  :: imn
    IF (.NOT. surface_ready) ERROR STOP '>> ERROR << vmec_eval_point called before vmec_select_surface'
    pt = vmec_point_t()
    pt%theta = theta; pt%zeta = zeta
    DO imn = 1,vmec_mnmax
      m = vmec_xm(imn); n = vmec_xn(imn)
      angle = m*theta - n*zeta
      c = COS(angle); s = SIN(angle)
      pt%R    = pt%R    + rmnc_s(imn)*c
      pt%dRdt = pt%dRdt - rmnc_s(imn)*m*s
      pt%dRdz = pt%dRdz + rmnc_s(imn)*n*s
      pt%Z    = pt%Z    + zmns_s(imn)*s
      pt%dZdt = pt%dZdt + zmns_s(imn)*m*c
      pt%dZdz = pt%dZdz - zmns_s(imn)*n*c
      pt%Lt   = pt%Lt   + lmns_s(imn)*m*c
      pt%Lz   = pt%Lz   - lmns_s(imn)*n*c
      pt%dRds = pt%dRds + drmnc_s(imn)*c
      pt%dZds = pt%dZds + dzmns_s(imn)*s
      pt%Ls   = pt%Ls   + dlmns_s(imn)*s
      IF (vmec_lasym) THEN
        pt%R    = pt%R    + rmns_s(imn)*s
        pt%dRdt = pt%dRdt + rmns_s(imn)*m*c
        pt%dRdz = pt%dRdz - rmns_s(imn)*n*c
        pt%Z    = pt%Z    + zmnc_s(imn)*c
        pt%dZdt = pt%dZdt - zmnc_s(imn)*m*s
        pt%dZdz = pt%dZdz + zmnc_s(imn)*n*s
        pt%Lt   = pt%Lt   - lmnc_s(imn)*m*s
        pt%Lz   = pt%Lz   + lmnc_s(imn)*n*s
        pt%dRds = pt%dRds + drmns_s(imn)*s
        pt%dZds = pt%dZds + dzmnc_s(imn)*c
        pt%Ls   = pt%Ls   + dlmnc_s(imn)*c
      ENDIF
    ENDDO
    DO imn = 1,vmec_mnmax_nyq
      m = vmec_xm_nyq(imn); n = vmec_xn_nyq(imn)
      angle = m*theta - n*zeta
      c = COS(angle); s = SIN(angle)
      pt%B      = pt%B      + bmnc_s(imn)*c
      pt%dBdt   = pt%dBdt   - bmnc_s(imn)*m*s
      pt%dBdz   = pt%dBdz   + bmnc_s(imn)*n*s
      pt%sqrtg  = pt%sqrtg  + gmnc_s(imn)*c
      pt%Bsup_t = pt%Bsup_t + bsupumnc_s(imn)*c
      pt%Bsup_z = pt%Bsup_z + bsupvmnc_s(imn)*c
      pt%Bsub_t = pt%Bsub_t + bsubumnc_s(imn)*c
      pt%Bsub_z = pt%Bsub_z + bsubvmnc_s(imn)*c
      pt%Bsub_s = pt%Bsub_s + bsubsmns_s(imn)*s
      pt%dBds   = pt%dBds   + dbmnc_s(imn)*c
      IF (vmec_lasym) THEN
        pt%B      = pt%B      + bmns_s(imn)*s
        pt%dBdt   = pt%dBdt   + bmns_s(imn)*m*c
        pt%dBdz   = pt%dBdz   - bmns_s(imn)*n*c
        pt%sqrtg  = pt%sqrtg  + gmns_s(imn)*s
        pt%Bsup_t = pt%Bsup_t + bsupumns_s(imn)*s
        pt%Bsup_z = pt%Bsup_z + bsupvmns_s(imn)*s
        pt%Bsub_t = pt%Bsub_t + bsubumns_s(imn)*s
        pt%Bsub_z = pt%Bsub_z + bsubvmns_s(imn)*s
        pt%Bsub_s = pt%Bsub_s + bsubsmnc_s(imn)*c
        pt%dBds   = pt%dBds   + dbmns_s(imn)*s
      ENDIF
    ENDDO
  END SUBROUTINE vmec_eval_point

  !============================================================================
  ! VMEC poloidal angle theta on the field line alpha at the toroidal angle zeta:
  ! root of theta + Lambda(theta,zeta) - (alpha + iota zeta) = 0 (Brent's method)
  !============================================================================
  REAL(dp) FUNCTION vmec_theta_of(alpha, zeta) RESULT(theta)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: alpha, zeta
    INTEGER,  PARAMETER  :: itmax_bracket = 10, itmax_root = 200
    REAL(dp), PARAMETER  :: tol = 1e-13_dp
    REAL(dp) :: theta_p, a, b, c, d, e, fa, fb, fc, p, q, r, s, tol1, xm
    INTEGER  :: it
    IF (.NOT. surface_ready) ERROR STOP '>> ERROR << vmec_theta_of called before vmec_select_surface'
    theta_p = alpha + vmec_iota*zeta
    a = theta_p - 0.3_dp; b = theta_p + 0.3_dp
    fa = residual(a); fb = residual(b)
    DO it = 1,itmax_bracket
      IF (fa*fb .GT. 0._dp) THEN
        a = a - 0.3_dp; b = b + 0.3_dp
        fa = residual(a); fb = residual(b)
      ELSE
        EXIT
      ENDIF
    ENDDO
    IF (fa*fb .GT. 0._dp) ERROR STOP '>> ERROR << vmec_theta_of: could not bracket the root for theta'
    ! Brent (Numerical Recipes zbrent)
    c = b; fc = fb; d = b - a; e = d
    DO it = 1,itmax_root
      IF (fb*fc .GT. 0._dp) THEN
        c = a; fc = fa; d = b - a; e = d
      ENDIF
      IF (ABS(fc) .LT. ABS(fb)) THEN
        a = b; b = c; c = a
        fa = fb; fb = fc; fc = fa
      ENDIF
      tol1 = 2._dp*EPSILON(1._dp)*ABS(b) + 0.5_dp*tol
      xm   = 0.5_dp*(c - b)
      IF (ABS(xm) .LE. tol1 .OR. fb .EQ. 0._dp) THEN
        theta = b
        RETURN
      ENDIF
      IF (ABS(e) .GE. tol1 .AND. ABS(fa) .GT. ABS(fb)) THEN
        s = fb/fa
        IF (a .EQ. c) THEN
          p = 2._dp*xm*s; q = 1._dp - s
        ELSE
          q = fa/fc; r = fb/fc
          p = s*(2._dp*xm*q*(q - r) - (b - a)*(r - 1._dp))
          q = (q - 1._dp)*(r - 1._dp)*(s - 1._dp)
        ENDIF
        IF (p .GT. 0._dp) q = -q
        p = ABS(p)
        IF (2._dp*p .LT. MIN(3._dp*xm*q - ABS(tol1*q), ABS(e*q))) THEN
          e = d; d = p/q
        ELSE
          d = xm; e = d
        ENDIF
      ELSE
        d = xm; e = d
      ENDIF
      a = b; fa = fb
      IF (ABS(d) .GT. tol1) THEN
        b = b + d
      ELSE
        b = b + SIGN(tol1, xm)
      ENDIF
      fb = residual(b)
    ENDDO
    ERROR STOP '>> ERROR << vmec_theta_of: Brent iteration did not converge'
  CONTAINS
    REAL(dp) FUNCTION residual(t)
      REAL(dp), INTENT(IN) :: t
      residual = t + vmec_lambda(t, zeta) - theta_p
    END FUNCTION residual
  END FUNCTION vmec_theta_of

  !============================================================================
  ! Define the flux tube (field line, toroidal extent, reference length, signs)
  !============================================================================
  SUBROUTINE vmec_setup_tube(alpha0, zeta_center, nfield_periods, Lref_opt)
    USE basic, ONLY: speak, str
    IMPLICIT NONE
    REAL(xp), INTENT(IN) :: alpha0, zeta_center, nfield_periods
    CHARACTER(len=*), INTENT(IN) :: Lref_opt
    TYPE(vmec_point_t) :: pt
    REAL(dp) :: theta
    IF (.NOT. surface_ready) ERROR STOP '>> ERROR << vmec_setup_tube called before vmec_select_surface'
    vmec_alpha0      = REAL(alpha0,dp)
    vmec_zeta_center = REAL(zeta_center,dp)
    IF (nfield_periods .GT. 0._xp) THEN
      vmec_nfield_periods = REAL(nfield_periods,dp)
    ELSE
      vmec_nfield_periods = REAL(vmec_nfp,dp) ! one full toroidal turn
    ENDIF
    SELECT CASE(TRIM(Lref_opt))
    CASE('R','r','major','Rmajor','R_major')
      vmec_L_ref = vmec_R_major
    CASE('a','A','minor','Aminor','a_minor')
      vmec_L_ref = vmec_a_minor
    CASE DEFAULT
      ERROR STOP '>> ERROR << Lref_opt must be ''R'' (major radius) or ''a'' (minor radius)'
    END SELECT
    vmec_Lref_over_a = vmec_L_ref/vmec_a_minor
    vmec_sigma_y = vmec_sign_torflux
    vmec_cx = REAL(vmec_sign_torflux,dp)/(vmec_rho0*vmec_a_minor*vmec_B_ref)
    vmec_cy = REAL(vmec_sigma_y,dp)*vmec_a_minor*vmec_rho0
    ! sigma_z from the sign of B^zeta at the centre of the tube
    theta = vmec_theta_of(vmec_alpha0, vmec_zeta_center)
    CALL vmec_eval_point(theta, vmec_zeta_center, pt)
    vmec_sigma_z = INT(SIGN(1.1_dp, pt%Bsup_z))
    vmec_cz = REAL(vmec_sigma_z,dp)*REAL(vmec_nfp,dp)/vmec_nfield_periods
    tube_ready = .TRUE.
    CALL speak('VMEC flux tube: alpha0 = '//str(REAL(vmec_alpha0,xp))//', zeta_center = '//str(REAL(vmec_zeta_center,xp))//&
               ', nfield_periods = '//str(REAL(vmec_nfield_periods,xp))//', L_ref = '//str(REAL(vmec_L_ref,xp))//&
               ' m, sigma_y = '//str(vmec_sigma_y)//', sigma_z = '//str(vmec_sigma_z),2)
  END SUBROUTINE vmec_setup_tube

  !============================================================================
  ! GYACOMO geometry at the point z of the flux tube
  !============================================================================
  SUBROUTINE vmec_geometry_at_z(z, g)
    IMPLICIT NONE
    REAL(dp), INTENT(IN) :: z
    TYPE(vmec_geo_t), INTENT(OUT) :: g
    TYPE(vmec_point_t) :: pt
    REAL(dp) :: zeta, theta, cosz, sinz
    REAL(dp) :: dXds, dXdz, dXdt, dYds, dYdz, dYdt, dZds, dZdz, dZdt
    REAL(dp) :: grad_s(3), grad_t(3), grad_z(3), grad_psit(3), grad_alpha(3)
    REAL(dp) :: Lref, cx, cy, cz, fac_s
    IF (.NOT. tube_ready) ERROR STOP '>> ERROR << vmec_geometry_at_z called before vmec_setup_tube'
    Lref = vmec_L_ref; cx = vmec_cx; cy = vmec_cy; cz = vmec_cz
    zeta  = vmec_zeta_center + z/cz
    theta = vmec_theta_of(vmec_alpha0, zeta)
    CALL vmec_eval_point(theta, zeta, pt)
    g = vmec_geo_t()
    g%zeta = zeta; g%theta = theta; g%R = pt%R; g%Z = pt%Z; g%Bsup_z = pt%Bsup_z
    ! Cartesian derivatives of the position, X = R cos zeta, Y = R sin zeta
    cosz = COS(zeta); sinz = SIN(zeta)
    dXds = pt%dRds*cosz; dXdz = pt%dRdz*cosz - pt%R*sinz; dXdt = pt%dRdt*cosz
    dYds = pt%dRds*sinz; dYdz = pt%dRdz*sinz + pt%R*cosz; dYdt = pt%dRdt*sinz
    dZds = pt%dZds;      dZdz = pt%dZdz;                   dZdt = pt%dZdt
    ! Dual relations: grad s = (e_theta x e_zeta)/sqrt g, grad theta = (e_zeta x e_s)/sqrt g
    grad_s(1) = (dYdt*dZdz - dZdt*dYdz)/pt%sqrtg
    grad_s(2) = (dZdt*dXdz - dXdt*dZdz)/pt%sqrtg
    grad_s(3) = (dXdt*dYdz - dYdt*dXdz)/pt%sqrtg
    grad_t(1) = (dYdz*dZds - dZdz*dYds)/pt%sqrtg
    grad_t(2) = (dZdz*dXds - dXdz*dZds)/pt%sqrtg
    grad_t(3) = (dXdz*dYds - dYdz*dXds)/pt%sqrtg
    grad_z    = [-sinz/pt%R, cosz/pt%R, 0._dp] ! exact for the cylindrical angle
    ! grad psi_t and grad alpha, alpha = theta + Lambda - iota zeta
    fac_s      = pt%Ls - zeta*vmec_diota_ds
    grad_psit  = vmec_psi_edge*grad_s
    grad_alpha = fac_s*grad_s + (1._dp + pt%Lt)*grad_t + (pt%Lz - vmec_iota)*grad_z
    ! GYACOMO coordinates: grad x, grad y dimensionless, L_ref grad z
    g%gradx = cx*grad_psit
    g%grady = cy*grad_alpha
    g%gradz = Lref*cz*grad_z
    g%gxx = DOT_PRODUCT(g%gradx, g%gradx)
    g%gxy = DOT_PRODUCT(g%gradx, g%grady)
    g%gyy = DOT_PRODUCT(g%grady, g%grady)
    g%gxz = DOT_PRODUCT(g%gradx, g%gradz)
    g%gyz = DOT_PRODUCT(g%grady, g%gradz)
    g%gzz = DOT_PRODUCT(g%gradz, g%gradz)
    g%hatB     = pt%B/vmec_B_ref
    g%Jacobian = 1._dp/(Lref*cx*cy*cz*pt%Bsup_z)
    ! Derivatives of hatB at fixed (x,y), (x,z) and (y,z), in L_ref units:
    !   dtheta/dzeta|_(s,alpha) = (iota - Lz)/(1+Lt), dtheta/dalpha = 1/(1+Lt),
    !   dtheta/ds|_(alpha,zeta) = -(Ls - zeta iota')/(1+Lt), ds/dpsit = 1/psi_edge
    g%dBdz = (pt%dBdz + pt%dBdt*(vmec_iota - pt%Lz)/(1._dp + pt%Lt))/(cz*vmec_B_ref)
    g%dBdy = Lref/cy*pt%dBdt/((1._dp + pt%Lt)*vmec_B_ref)
    g%dBdx = Lref/cx/vmec_psi_edge*(pt%dBds - fac_s*pt%dBdt/(1._dp + pt%Lt))/vmec_B_ref
    ! For output/visualisation: dx/dR, dx/dZ in the poloidal plane
    g%dxdR = g%gradx(1)*cosz + g%gradx(2)*sinz
    g%dxdZ = g%gradx(3)
    ! B/B_ref from B = grad psi_t x grad alpha = psi_edge/sqrt g [(1+Lt) e_zeta - (Lz - iota) e_theta]
    g%Bvec(1) = vmec_psi_edge/pt%sqrtg*((1._dp + pt%Lt)*dXdz - (pt%Lz - vmec_iota)*dXdt)/vmec_B_ref
    g%Bvec(2) = vmec_psi_edge/pt%sqrtg*((1._dp + pt%Lt)*dYdz - (pt%Lz - vmec_iota)*dYdt)/vmec_B_ref
    g%Bvec(3) = vmec_psi_edge/pt%sqrtg*((1._dp + pt%Lt)*dZdz - (pt%Lz - vmec_iota)*dZdt)/vmec_B_ref
  END SUBROUTINE vmec_geometry_at_z

  !============================================================================
  ! Fill GYACOMO's geometry arrays on the local z grid(s), ghosts included, and
  ! compute the twist-and-shift factor Gamma (Delta kx = Gamma ky at z = +-pi)
  !============================================================================
  SUBROUTINE get_vmec_geometry(local_nz, ngz, nzgrid, zarray, twist_shift_opt, &
                               gxx_, gxy_, gxz_, gyy_, gyz_, gzz_, dBdx_, dBdy_, dBdz_, hatB_, Jacobian_, &
                               hatR_, hatZ_, Rc_, phic_, Zc_, dxdR_, dxdZ_, zeta_, theta_, twist_shift_fac)
    USE basic, ONLY: speak, str
    IMPLICIT NONE
    INTEGER,  INTENT(IN) :: local_nz, ngz, nzgrid
    REAL(xp), DIMENSION(local_nz+ngz,nzgrid), INTENT(IN) :: zarray
    CHARACTER(len=*), INTENT(IN) :: twist_shift_opt
    REAL(xp), DIMENSION(local_nz+ngz,nzgrid), INTENT(OUT) :: gxx_, gxy_, gxz_, gyy_, gyz_, gzz_
    REAL(xp), DIMENSION(local_nz+ngz,nzgrid), INTENT(OUT) :: dBdx_, dBdy_, dBdz_, hatB_, Jacobian_
    REAL(xp), DIMENSION(local_nz+ngz,nzgrid), INTENT(OUT) :: hatR_, hatZ_, Rc_, phic_, Zc_, dxdR_, dxdZ_, zeta_, theta_
    REAL(xp), INTENT(OUT) :: twist_shift_fac
    TYPE(vmec_geo_t) :: g, gp, gm
    INTEGER  :: iz, eo
    REAL(dp) :: lam_p, lam_m, gam, tol
    LOGICAL  :: sign_ok
    IF (.NOT. tube_ready) ERROR STOP '>> ERROR << get_vmec_geometry called before vmec_setup_tube'
    sign_ok = .TRUE.
    DO eo = 1,nzgrid
      DO iz = 1,local_nz+ngz
        CALL vmec_geometry_at_z(REAL(zarray(iz,eo),dp), g)
        gxx_(iz,eo) = REAL(g%gxx,xp); gxy_(iz,eo) = REAL(g%gxy,xp); gxz_(iz,eo) = REAL(g%gxz,xp)
        gyy_(iz,eo) = REAL(g%gyy,xp); gyz_(iz,eo) = REAL(g%gyz,xp); gzz_(iz,eo) = REAL(g%gzz,xp)
        dBdx_(iz,eo) = REAL(g%dBdx,xp); dBdy_(iz,eo) = REAL(g%dBdy,xp); dBdz_(iz,eo) = REAL(g%dBdz,xp)
        hatB_(iz,eo) = REAL(g%hatB,xp); Jacobian_(iz,eo) = REAL(g%Jacobian,xp)
        hatR_(iz,eo) = REAL(g%R/vmec_L_ref,xp); hatZ_(iz,eo) = REAL(g%Z/vmec_L_ref,xp)
        Rc_(iz,eo) = REAL(g%R,xp); Zc_(iz,eo) = REAL(g%Z,xp); phic_(iz,eo) = REAL(g%zeta,xp)
        dxdR_(iz,eo) = REAL(g%dxdR,xp); dxdZ_(iz,eo) = REAL(g%dxdZ,xp)
        zeta_(iz,eo) = REAL(g%zeta,xp); theta_(iz,eo) = REAL(g%theta,xp)
        IF (INT(SIGN(1.1_dp, g%Bsup_z)) .NE. vmec_sigma_z) sign_ok = .FALSE.
      ENDDO
    ENDDO
    IF (.NOT. sign_ok) ERROR STOP '>> ERROR << B.grad zeta changes sign along the flux tube'
    ! Twist-and-shift factor from the two ends of the tube (same physical point
    ! z = +-pi on every rank, so Lx and the kx map are consistent across ranks)
    CALL vmec_geometry_at_z( PI_DP, gp)
    CALL vmec_geometry_at_z(-PI_DP, gm)
    lam_p = gp%gxy/gp%gxx
    lam_m = gm%gxy/gm%gxx
    SELECT CASE(TRIM(twist_shift_opt))
    CASE('local','stellarator') ! Martin et al. PPCF 2018: local metric at the ends
      gam = lam_p - lam_m
    CASE('global','linked')     ! global shear over the tube: 2 pi shat iota P/nfp
      gam = REAL(vmec_sigma_y*vmec_sigma_z,dp)*2._dp*PI_DP*vmec_shat*vmec_iota*vmec_nfield_periods/REAL(vmec_nfp,dp)
    CASE DEFAULT
      ERROR STOP '>> ERROR << twist_shift_opt must be ''local'' or ''global'''
    END SELECT
    twist_shift_fac = REAL(gam,xp)
    tol = 1e-3_dp
    IF (ABS(gp%gxx - gm%gxx) .GT. tol*MAX(gp%gxx, gm%gxx) .OR. &
        ABS((gp%gyy - gp%gxy**2/gp%gxx) - (gm%gyy - gm%gxy**2/gm%gxx)) .GT. tol*MAX(gp%gyy, gm%gyy) .OR. &
        ABS(gp%gxy + gm%gxy) .GT. tol*MAX(ABS(gp%gxy), ABS(gm%gxy), 1e-10_dp)) THEN
      CALL speak('WARNING: the flux tube is not stellarator symmetric about its ends,',1)
      CALL speak('         k_perp will be discontinuous across the parallel boundary',1)
    ENDIF
    CALL speak('VMEC twist-and-shift factor ('//TRIM(twist_shift_opt)//'): Gamma = '//str(twist_shift_fac)//&
               '  (local ends: '//str(REAL(lam_p-lam_m,xp))//', global: '//&
               str(REAL(REAL(vmec_sigma_y*vmec_sigma_z,dp)*2._dp*PI_DP*vmec_shat*vmec_iota*vmec_nfield_periods/REAL(vmec_nfp,dp),xp))//')',2)
  END SUBROUTINE get_vmec_geometry

END MODULE vmec
