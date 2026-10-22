PROGRAM test_vmec
  ! Checks the VMEC flux-tube geometry (vmec_mod) on the W7-X "kjm" equilibrium
  ! that stella ships in its test suite (unit/fixtures/wout_w7x_kjm.h5, made by
  ! scripts/python_utilities/wout2h5.py).
  !
  ! 1. the HDF5 reader and the surface selection reproduce the scalars that
  !    stella prints for torflux = 0.49 (rho = 0.7)
  ! 2. every derivative of the Fourier sums matches a finite difference
  ! 3. the geometry on the tube satisfies the identities that GYACOMO relies on:
  !    B = grad x X grad y (C_xy = 1, so hatB^2 = gxx gyy - gxy^2), J > 0,
  !    b.grad z > 0, dBdz is the derivative of hatB along z, dBdx and dBdy are
  !    the derivatives across surfaces and field lines
  ! 4. in stella's units (L_ref = a) hatB, b.grad z, the metric and the grad-B
  !    drift coefficients match stella's EXPECTED_OUTPUT.vmec_geometry.geometry
  !    table (4 digits) once the sign conventions sigma_y, sigma_z are applied
  ! 5. through the full input pipeline (readinputs + auxval) the twist-and-shift
  !    factor, the kx grid and the parallel-boundary kx map come out consistent
  !
  ! Single rank.
  USE prec_const, ONLY: xp, dp, PI
  USE basic,      ONLY: lu_in
  USE parallel,   ONLY: ppinit
  USE vmec
  USE test_harness
  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER :: fixture = __GYACDIR__//'/unit/fixtures/wout_w7x_kjm.h5'
  ! stella reference, test_2_geometry/EXPECTED_OUTPUT.vmec_geometry.geometry
  ! (torflux 0.49, nfield_periods 6.073867, alpha0 0, zeta_center 0, nzed 24):
  ! zed, bmag, bdot_grad_z, gds2, gds21, gds22, gbdrift, gbdrift0
  INTEGER,  PARAMETER :: NREF = 25
  REAL(xp), PARAMETER :: ref(8,NREF) = RESHAPE([ &
      -3.1420E+00_xp, 1.2600E+00_xp, 6.6130E-02_xp, 1.5020E+00_xp, -1.8660E-01_xp, 3.0020E-02_xp, -3.4000E-01_xp, -4.2450E-04_xp,&
      -2.8800E+00_xp, 1.0940E+00_xp, 7.0490E-02_xp, 1.9210E+00_xp, -1.0520E-01_xp, 9.7980E-03_xp, -1.5330E-01_xp, -1.3470E-02_xp,&
      -2.6180E+00_xp, 9.2290E-01_xp, 8.0050E-02_xp, 1.2780E+00_xp, -1.6730E-02_xp, 4.5310E-03_xp, 1.4430E-02_xp, -3.2010E-03_xp,&
      -2.3560E+00_xp, 1.0050E+00_xp, 7.5140E-02_xp, 1.7490E+00_xp, 6.2410E-03_xp, 3.7580E-03_xp, -9.4990E-02_xp, 9.3930E-03_xp,&
      -2.0940E+00_xp, 1.2100E+00_xp, 6.7710E-02_xp, 2.8990E+00_xp, 1.2340E-01_xp, 8.5210E-03_xp, -2.4260E-01_xp, 1.1340E-02_xp,&
      -1.8330E+00_xp, 1.1520E+00_xp, 6.8030E-02_xp, 8.0380E-01_xp, 5.8470E-02_xp, 1.4950E-02_xp, -3.7610E-01_xp, 7.0690E-04_xp,&
      -1.5710E+00_xp, 9.5370E-01_xp, 7.5620E-02_xp, 1.1310E+00_xp, -5.5330E-02_xp, 7.9130E-03_xp, -3.5610E-02_xp, -1.3790E-02_xp,&
      -1.3090E+00_xp, 9.8070E-01_xp, 7.3210E-02_xp, 1.7440E+00_xp, -7.4550E-02_xp, 6.7580E-03_xp, 2.6180E-03_xp, 6.2240E-03_xp,&
      -1.0470E+00_xp, 1.1410E+00_xp, 6.6990E-02_xp, 3.9790E+00_xp, -1.4150E-01_xp, 7.1530E-03_xp, -1.1090E-01_xp, 1.4840E-02_xp,&
      -7.8540E-01_xp, 1.0950E+00_xp, 6.7720E-02_xp, 3.3820E+00_xp, -2.5110E-02_xp, 2.4830E-03_xp, -1.7130E-01_xp, 1.3720E-02_xp,&
      -5.2360E-01_xp, 9.8730E-01_xp, 7.2280E-02_xp, 1.0680E+00_xp, 6.5350E-02_xp, 9.9110E-03_xp, -1.5280E-01_xp, 3.3600E-03_xp,&
      -2.6180E-01_xp, 1.0140E+00_xp, 6.8150E-02_xp, 3.5760E-01_xp, 4.8620E-02_xp, 2.5220E-02_xp, -1.5910E-02_xp, -1.4480E-02_xp,&
      0.0000E+00_xp, 1.0910E+00_xp, 6.4580E-02_xp, 2.1230E-01_xp, -0.0000E+00_xp, 3.6320E-02_xp, 7.0970E-02_xp, 0.0000E+00_xp,&
      2.6180E-01_xp, 1.0140E+00_xp, 6.8150E-02_xp, 3.5760E-01_xp, -4.8620E-02_xp, 2.5220E-02_xp, -1.5910E-02_xp, 1.4480E-02_xp,&
      5.2360E-01_xp, 9.8730E-01_xp, 7.2280E-02_xp, 1.0680E+00_xp, -6.5350E-02_xp, 9.9110E-03_xp, -1.5280E-01_xp, -3.3600E-03_xp,&
      7.8540E-01_xp, 1.0950E+00_xp, 6.7720E-02_xp, 3.3820E+00_xp, 2.5110E-02_xp, 2.4830E-03_xp, -1.7130E-01_xp, -1.3720E-02_xp,&
      1.0470E+00_xp, 1.1410E+00_xp, 6.6990E-02_xp, 3.9790E+00_xp, 1.4150E-01_xp, 7.1530E-03_xp, -1.1090E-01_xp, -1.4840E-02_xp,&
      1.3090E+00_xp, 9.8070E-01_xp, 7.3210E-02_xp, 1.7440E+00_xp, 7.4550E-02_xp, 6.7580E-03_xp, 2.6180E-03_xp, -6.2240E-03_xp,&
      1.5710E+00_xp, 9.5370E-01_xp, 7.5620E-02_xp, 1.1310E+00_xp, 5.5330E-02_xp, 7.9130E-03_xp, -3.5610E-02_xp, 1.3790E-02_xp,&
      1.8330E+00_xp, 1.1520E+00_xp, 6.8030E-02_xp, 8.0380E-01_xp, -5.8470E-02_xp, 1.4950E-02_xp, -3.7610E-01_xp, -7.0690E-04_xp,&
      2.0940E+00_xp, 1.2100E+00_xp, 6.7710E-02_xp, 2.8990E+00_xp, -1.2340E-01_xp, 8.5210E-03_xp, -2.4260E-01_xp, -1.1340E-02_xp,&
      2.3560E+00_xp, 1.0050E+00_xp, 7.5140E-02_xp, 1.7490E+00_xp, -6.2410E-03_xp, 3.7580E-03_xp, -9.4990E-02_xp, -9.3930E-03_xp,&
      2.6180E+00_xp, 9.2290E-01_xp, 8.0050E-02_xp, 1.2780E+00_xp, 1.6730E-02_xp, 4.5310E-03_xp, 1.4430E-02_xp, 3.2010E-03_xp,&
      2.8800E+00_xp, 1.0940E+00_xp, 7.0490E-02_xp, 1.9210E+00_xp, 1.0520E-01_xp, 9.7980E-03_xp, -1.5330E-01_xp, 1.3470E-02_xp,&
      3.1420E+00_xp, 1.2600E+00_xp, 6.6130E-02_xp, 1.5020E+00_xp, 1.8660E-01_xp, 3.0020E-02_xp, -3.4000E-01_xp, 4.2450E-04_xp &
      ], [8,NREF])
  REAL(xp), PARAMETER :: ref_rho = 0.7_xp, ref_q = 1.1123926_xp, ref_shat = -0.080451439_xp
  REAL(xp), PARAMETER :: ref_iota = 0.89896317_xp, ref_psi_edge = -0.31830989_xp
  REAL(xp), PARAMETER :: ref_a = 0.49441695_xp, ref_R = 5.4854007_xp, ref_Bref = 2.6043144_xp
  REAL(xp), PARAMETER :: ref_B00 = 2.8411858_xp, ref_Bsupz00 = 0.45080708_xp
  REAL(xp), PARAMETER :: nfp_in = 6.073867_xp

  TYPE(vmec_point_t) :: pt, ptp, ptm
  TYPE(vmec_geo_t)   :: g, gp, gm
  REAL(dp) :: h, th, ze, fd, Bp, Bm, xp_, xm_, trip, gcov_xz, gcov_yz, gcov_zz, Cx, Cy
  REAL(xp) :: tol4, want
  INTEGER  :: j, i

  CALL ppinit
  CALL start_suite('vmec')

  !--------------------------------------------------------------- 1. reader
  CALL section('HDF5 reader: scalars and array layout of the W7-X kjm file')
  CALL vmec_read_h5(fixture)
  CALL check_int('nfp', vmec_nfp, 5)
  CALL check_int('ns', vmec_ns, 99)
  CALL check_int('mnmax', vmec_mnmax, 288)
  CALL check_int('mnmax_nyq', vmec_mnmax_nyq, 288)
  CALL check_int('signgs', vmec_signgs, -1)
  CALL check('lasym is false', .NOT. vmec_lasym)
  CALL check_real('Aminor_p', REAL(vmec_a_minor,xp), ref_a, 1e-7_xp)
  CALL check_real('Rmajor_p', REAL(vmec_R_major,xp), ref_R, 1e-7_xp)
  CALL check_real('phi(ns) = 2', REAL(vmec_phi(vmec_ns),xp), 2._xp, 1e-12_xp)
  CALL check_real('xm(2) = 0 (array layout)', REAL(vmec_xm(2),xp), 0._xp, 1e-14_xp)
  CALL check_real('xn(2) = nfp (array layout)', REAL(vmec_xn(2),xp), 5._xp, 1e-14_xp)
  CALL check_real('psi_edge = phi(ns) signgs/2pi', REAL(vmec_psi_edge,xp), ref_psi_edge, 1e-7_xp)
  CALL check_int('sign_torflux', vmec_sign_torflux, -1)
  CALL check_real('B_ref = 2|psi_edge|/a^2', REAL(vmec_B_ref,xp), ref_Bref, 1e-7_xp)

  !--------------------------------------------------------------- 2. surface
  CALL section('surface s = 0.49: profiles interpolated as in stella')
  CALL vmec_select_surface(0.49_xp, 0)
  CALL check_real('rho0 = sqrt(s)', REAL(vmec_rho0,xp), ref_rho, 1e-12_xp)
  CALL check_real('iota', REAL(vmec_iota,xp), ref_iota, 1e-7_xp)
  CALL check_real('q = 1/iota', REAL(vmec_q0,xp), ref_q, 1e-7_xp)
  CALL check_real('shat = -2 s iota''/iota', REAL(vmec_shat,xp), ref_shat, 1e-6_xp)

  !--------------------------------------------------------------- 3. point values
  CALL section('Fourier sums at (theta,zeta) = (0,0) and derivative checks')
  CALL vmec_eval_point(0._dp, 0._dp, pt)
  CALL check_real('|B|(0,0) [T]', REAL(pt%B,xp), ref_B00, 1e-7_xp)
  CALL check_real('B^zeta(0,0) [T/m]', REAL(pt%Bsup_z,xp), ref_Bsupz00, 1e-7_xp)
  CALL check_real('Lambda(0,0) = 0 (stellarator symmetry)', REAL(vmec_lambda(0._dp,0._dp),xp), 0._xp, 1e-12_xp)
  CALL check_real('theta solver at alpha = 0, zeta = 0', REAL(vmec_theta_of(0._dp,0._dp),xp), 0._xp, 1e-10_xp)
  ! generic point: analytic derivatives against centred differences
  th = 0.7_dp; ze = 0.4_dp; h = 1e-4_dp
  CALL vmec_eval_point(th, ze, pt)
  CALL vmec_eval_point(th+h, ze, ptp); CALL vmec_eval_point(th-h, ze, ptm)
  CALL check_real('dR/dtheta',  REAL(pt%dRdt,xp), REAL((ptp%R-ptm%R)/(2*h),xp), 1e-6_xp)
  CALL check_real('dZ/dtheta',  REAL(pt%dZdt,xp), REAL((ptp%Z-ptm%Z)/(2*h),xp), 1e-6_xp)
  CALL check_real('dB/dtheta',  REAL(pt%dBdt,xp), REAL((ptp%B-ptm%B)/(2*h),xp), 1e-6_xp)
  CALL check_real('dLambda/dtheta', REAL(pt%Lt,xp), REAL((vmec_lambda(th+h,ze)-vmec_lambda(th-h,ze))/(2*h),xp), 1e-6_xp)
  CALL vmec_eval_point(th, ze+h, ptp); CALL vmec_eval_point(th, ze-h, ptm)
  CALL check_real('dR/dzeta',   REAL(pt%dRdz,xp), REAL((ptp%R-ptm%R)/(2*h),xp), 1e-6_xp)
  CALL check_real('dZ/dzeta',   REAL(pt%dZdz,xp), REAL((ptp%Z-ptm%Z)/(2*h),xp), 1e-6_xp)
  CALL check_real('dB/dzeta',   REAL(pt%dBdz,xp), REAL((ptp%B-ptm%B)/(2*h),xp), 1e-6_xp)
  CALL check_real('dLambda/dzeta', REAL(pt%Lz,xp), REAL((vmec_lambda(th,ze+h)-vmec_lambda(th,ze-h))/(2*h),xp), 1e-6_xp)
  ! theta solver: the root satisfies theta + Lambda = alpha + iota zeta
  th = vmec_theta_of(0.3_dp, 1.1_dp)
  CALL check_real('theta solver residual', REAL(th + vmec_lambda(th,1.1_dp) - (0.3_dp + vmec_iota*1.1_dp),xp), 0._xp, 1e-10_xp)

  !--------------------------------------------------------------- 4. tube in stella units
  CALL section('flux tube alpha0 = 0, nfield_periods = 6.073867, L_ref = a: signs and identities')
  CALL vmec_setup_tube(0._xp, 0._xp, nfp_in, 'a')
  CALL check_int('sigma_y = sign(psi_edge)', vmec_sigma_y, -1)
  CALL check_int('sigma_z = sign(B^zeta)', vmec_sigma_z, 1)
  CALL check_real('dz/dzeta = nfp/P', REAL(vmec_cz,xp), 5._xp/nfp_in, 1e-12_xp)
  CALL check_real('L_ref = a', REAL(vmec_L_ref,xp), ref_a, 1e-7_xp)
  DO j = 1,NREF
    CALL vmec_geometry_at_z(zed_of(j), g)
    ! B = grad x X grad y exactly: hatB^2 = gxx gyy - gxy^2 (VMEC-internal consistency, 1e-2)
    CALL check_real('hatB^2 = gxx gyy - gxy^2 at z'//str_of(j), REAL(g%gxx*g%gyy-g%gxy**2,xp), REAL(g%hatB**2,xp), 1e-2_xp)
    ! B/B_ref rebuilt from psi_edge/sqrt g agrees with bmnc
    CALL check_real('|Bvec| = hatB at z'//str_of(j), REAL(SQRT(DOT_PRODUCT(g%Bvec,g%Bvec)),xp), REAL(g%hatB,xp), 1e-2_xp)
    ! Jacobian: J (grad x X grad y . L_ref grad z) = 1
    trip = g%gradx(1)*(g%grady(2)*g%gradz(3)-g%grady(3)*g%gradz(2)) &
         + g%gradx(2)*(g%grady(3)*g%gradz(1)-g%grady(1)*g%gradz(3)) &
         + g%gradx(3)*(g%grady(1)*g%gradz(2)-g%grady(2)*g%gradz(1))
    CALL check_real('J (grad x X grad y . grad z) = 1 at z'//str_of(j), REAL(g%Jacobian*trip,xp), 1._xp, 1e-2_xp)
    CALL check('J > 0 at z'//str_of(j), g%Jacobian .GT. 0._dp)
    CALL check('b.grad z > 0 at z'//str_of(j), 1._dp/(g%Jacobian*g%hatB) .GT. 0._dp)
    CALL check_real('gzz = (nfp/P)^2 (L_ref/R)^2 at z'//str_of(j), REAL(g%gzz,xp), REAL(vmec_cz**2*(vmec_L_ref/g%R)**2,xp), 1e-12_xp)
    ! dBdz against a centred difference of hatB along the tube
    h = 1e-4_dp
    CALL vmec_geometry_at_z(zed_of(j)+h, gp); CALL vmec_geometry_at_z(zed_of(j)-h, gm)
    CALL check_real('dBdz = d hatB/dz at z'//str_of(j), REAL(g%dBdz,xp), REAL((gp%hatB-gm%hatB)/(2*h),xp), 1e-5_xp)
  ENDDO
  ! dBdy: finite difference across field lines at fixed (x,z): y = sigma_y rho0 (alpha - alpha0) in L_ref = a units
  h = 1e-4_dp
  CALL vmec_geometry_at_z(0.9_dp, g)
  CALL vmec_setup_tube(REAL( h,xp), 0._xp, nfp_in, 'a'); CALL vmec_geometry_at_z(0.9_dp, gp)
  CALL vmec_setup_tube(REAL(-h,xp), 0._xp, nfp_in, 'a'); CALL vmec_geometry_at_z(0.9_dp, gm)
  CALL vmec_setup_tube(0._xp, 0._xp, nfp_in, 'a')
  fd = (gp%hatB-gm%hatB)/(2*h*REAL(vmec_sigma_y,dp)*vmec_rho0)
  CALL check_real('dBdy = d hatB/dy across field lines', REAL(g%dBdy,xp), REAL(fd,xp), 1e-5_xp)
  ! dBdx: finite difference across surfaces at fixed (alpha,zeta): x = rho - rho0 in L_ref = a units.
  ! The s-derivatives of the Fourier coefficients are mesh differences, so only a loose agreement is expected
  h = 2e-3_dp
  Bp = 0._dp; Bm = 0._dp
  CALL vmec_select_surface(REAL(0.49_dp+h,xp), 0); CALL vmec_setup_tube(0._xp, 0._xp, nfp_in, 'a')
  xp_ = vmec_rho0; CALL vmec_geometry_at_z(0.9_dp, gp); Bp = gp%hatB
  CALL vmec_select_surface(REAL(0.49_dp-h,xp), 0); CALL vmec_setup_tube(0._xp, 0._xp, nfp_in, 'a')
  xm_ = vmec_rho0; CALL vmec_geometry_at_z(0.9_dp, gm); Bm = gm%hatB
  CALL vmec_select_surface(0.49_xp, 0); CALL vmec_setup_tube(0._xp, 0._xp, nfp_in, 'a')
  CALL check_real('dBdx = d hatB/dx across surfaces (loose)', REAL(g%dBdx,xp), REAL((Bp-Bm)/(xp_-xm_),xp), 5e-2_xp)

  !--------------------------------------------------------------- 4b. stella table
  CALL section('comparison with stella (4-digit table): hatB, b.grad z, metric, grad-B drifts')
  tol4 = 2e-3_xp
  DO j = 1,NREF
    CALL vmec_geometry_at_z(zed_of(j), g)
    CALL check_real('bmag at z'//str_of(j), REAL(g%hatB,xp), ref(2,j), tol4)
    CALL check_real('b.grad z at z'//str_of(j), REAL(1._dp/(g%Jacobian*g%hatB),xp), ref(3,j), tol4)
    CALL check_real('gyy = gds2 at z'//str_of(j), REAL(g%gyy,xp), ref(4,j), tol4)
    CALL check_real('gxy = sigma_y gds21/shat at z'//str_of(j), REAL(g%gxy,xp), REAL(vmec_sigma_y,xp)*ref(5,j)/ref_shat, tol4)
    CALL check_real('gxx = gds22/shat^2 at z'//str_of(j), REAL(g%gxx,xp), ref(6,j)/ref_shat**2, tol4)
    ! grad-B drift coefficients as evaluate_magn_curv builds them (C_xy = 1)
    gcov_xz = g%Jacobian**2*(g%gxy*g%gyz - g%gyy*g%gxz)
    gcov_yz = g%Jacobian**2*(g%gxy*g%gxz - g%gxx*g%gyz)
    gcov_zz = g%Jacobian**2*(g%gxx*g%gyy - g%gxy**2)
    Cx = -(g%dBdy - gcov_yz/gcov_zz*g%dBdz)/g%hatB
    Cy =  (g%dBdx - gcov_xz/gcov_zz*g%dBdz)/g%hatB
    CALL check_real('Cy = sigma_y gbdrift/2 at z'//str_of(j), REAL(Cy,xp), REAL(vmec_sigma_y,xp)*ref(7,j)/2._xp, tol4)
    CALL check_real('Cx = gbdrift0/(2 shat) at z'//str_of(j), REAL(Cx,xp), ref(8,j)/(2._xp*ref_shat), tol4)
  ENDDO

  !--------------------------------------------------------------- 5. full pipeline
  CALL section('full pipeline: twist-and-shift factor, kx grid and parallel kx map')
  CALL pipeline_checks

  CALL finish_suite

CONTAINS

  ! exact grid point j of stella's table (its zed column is printed with 4 digits only)
  REAL(dp) FUNCTION zed_of(j)
    INTEGER, INTENT(IN) :: j
    zed_of = REAL(j-13,dp)*ACOS(-1._dp)/12._dp
  END FUNCTION zed_of

  FUNCTION str_of(j) RESULT(s)
    INTEGER, INTENT(IN) :: j
    CHARACTER(len=8) :: s
    WRITE(s,'(A,I0)') '(', j
    s = TRIM(s)//')'
  END FUNCTION str_of

  SUBROUTINE pipeline_checks
    USE grid,     ONLY: Lx, Ly, deltakx, deltaky, total_nkx, local_nky, ikx0
    USE geometry, ONLY: twist_shift_fac, SHEARED, ikx_zBC_L, ikx_zBC_R, hatB, gxx, gxy, gradz_coeff, dBdy, C_xy, Cyq0_x0
    USE grid,     ONLY: ngz_o2, local_nz, ieven, zarray
    REAL(xp) :: gam
    INTEGER  :: ikx, iz
    CALL write_input_deck
    CALL readinputs
    CALL auxval
    gam = twist_shift_fac
    CALL check_real('Gamma = Lam(+pi) - Lam(-pi) ~ 1 (stella: |fac| = 1.0002)', gam, 1.0002_xp, 2e-3_xp)
    CALL check('SHEARED set from Gamma', SHEARED)
    CALL check_real('C_xy = 1', C_xy, 1._xp, 1e-14_xp)
    CALL check_real('Lx = Nexc Ly/|Gamma|', Lx, Ly/ABS(gam), 1e-12_xp)
    CALL check_real('dkx = |Gamma| dky', deltakx, ABS(gam)*deltaky, 1e-12_xp)
    ! kx map: with Gamma > 0 and Nexc = 1, (kx, ky1) at z = +pi connects to kx + dkx at z = -pi
    ! FFT ordering 0, 1, 2, 3, 4, -3, -2, -1 (Nx = 8): ikx = 1 -> 2 (right), kx = 4 dkx (ikx 5) has no partner
    CALL check_int('ikx_zBC_R(ky1, kx=0) -> kx=dkx', ikx_zBC_R(2,1), 2)
    CALL check_int('ikx_zBC_R(ky1, kx=3dkx) -> kx=4dkx', ikx_zBC_R(2,4), 5)
    CALL check_int('ikx_zBC_R(ky1, kx=4dkx) -> dirichlet', ikx_zBC_R(2,5), -99)
    CALL check_int('ikx_zBC_L(ky1, kx=dkx) -> kx=0', ikx_zBC_L(2,2), 1)
    CALL check_int('ikx_zBC_L(ky1, kx=-3dkx) -> dirichlet', ikx_zBC_L(2,6), -99)
    CALL check_int('ikx_zBC_L(ky0) identity', ikx_zBC_L(1,3), 3)
    ! geometry arrays on the even grid: hatB(z=0) and the midplane values of stella (L_ref = a)
    iz = local_nz/2 + 1 + ngz_o2
    CALL check_real('zarray at the midplane is 0', zarray(iz,ieven), 0._xp, 1e-12_xp)
    CALL check_real('hatB(z=0)', hatB(iz,ieven), 1.091_xp, 2e-3_xp)
    CALL check_real('gradz_coeff(z=0) = b.grad z', gradz_coeff(iz,ieven), 0.06458_xp, 2e-3_xp)
    CALL check_real('gxy(z=0) = 0', gxy(iz,ieven), 0._xp, 1e-6_xp)
    CALL check_real('gxx(z=0) = gds22/shat^2', gxx(iz,ieven), 0.03632_xp/ref_shat**2, 2e-3_xp)
    CALL check('dBdy is not identically zero', MAXVAL(ABS(dBdy)) .GT. 1e-3_xp)
  END SUBROUTINE pipeline_checks

  SUBROUTINE write_input_deck
    OPEN(UNIT=lu_in, STATUS='SCRATCH', ACTION='READWRITE', FORM='FORMATTED')
    WRITE(lu_in,'(A)')       '&GRID'
    WRITE(lu_in,'(A)')       '  pmax = 2'
    WRITE(lu_in,'(A)')       '  jmax = 1'
    WRITE(lu_in,'(A)')       '  Nx   = 8'
    WRITE(lu_in,'(A)')       '  Lx   = 20.943951'
    WRITE(lu_in,'(A)')       '  Ny   = 2'
    WRITE(lu_in,'(A)')       '  Ly   = 20.943951'
    WRITE(lu_in,'(A)')       '  Nz   = 24'
    WRITE(lu_in,'(A)')       '  SG   = .false.'
    WRITE(lu_in,'(A)')       '  Nexc = 1'
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&GEOMETRY'
    WRITE(lu_in,'(A)')       "  geom            = 'vmec'"
    WRITE(lu_in,'(A)')       "  vmec_file       = '"//fixture//"'"
    WRITE(lu_in,'(A)')       '  torflux         = 0.49'
    WRITE(lu_in,'(A)')       '  alpha0          = 0.0'
    WRITE(lu_in,'(A)')       '  zeta_center     = 0.0'
    WRITE(lu_in,'(A)')       '  nfield_periods  = 6.073867'
    WRITE(lu_in,'(A)')       '  surface_option  = 0'
    WRITE(lu_in,'(A)')       "  Lref_opt        = 'a'"
    WRITE(lu_in,'(A)')       "  twist_shift_opt = 'local'"
    WRITE(lu_in,'(A)')       "  parallel_bc     = 'dirichlet'"
    WRITE(lu_in,'(A)')       '  Npol            = 1'
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&DIAGNOSTICS'
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&MODEL'
    WRITE(lu_in,'(A)')       "  LINEARITY = 'linear'"
    WRITE(lu_in,'(A)')       '  Na      = 1'
    WRITE(lu_in,'(A)')       '  ADIAB_E = .true.'
    WRITE(lu_in,'(A)')       '  beta    = 0.0'
    WRITE(lu_in,'(A)')       '  nu      = 0.0'
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&CLOSURE'
    WRITE(lu_in,'(A)')       "  hierarchy_closure = 'truncation'"
    WRITE(lu_in,'(A)')       "  nonlinear_closure = 'truncation'"
    WRITE(lu_in,'(A)')       '  nmax = 0'
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&SPECIES'
    WRITE(lu_in,'(A)')       "  name_  = 'ions'"
    WRITE(lu_in,'(A)')       '  tau_   = 1.0'
    WRITE(lu_in,'(A)')       '  sigma_ = 1.0'
    WRITE(lu_in,'(A)')       '  q_     = 1.0'
    WRITE(lu_in,'(A)')       '  k_N_   = 1.0'
    WRITE(lu_in,'(A)')       '  k_T_   = 3.0'
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&COLLISION'
    WRITE(lu_in,'(A)')       "  collision_model = 'DG'"
    WRITE(lu_in,'(A)')       '  GK_CO           = .true.'
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&INITIAL'
    WRITE(lu_in,'(A)')       "  INIT_OPT = 'phi'"
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&TIME_INTEGRATION'
    WRITE(lu_in,'(A)')       "  numerical_scheme = 'RK4'"
    WRITE(lu_in,'(A)')       '/'
    REWIND(lu_in)
  END SUBROUTINE write_input_deck

END PROGRAM test_vmec
