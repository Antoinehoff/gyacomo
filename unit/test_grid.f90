PROGRAM test_grid
  ! Builds a small grid the way auxval does and checks it end to end.
  !
  ! grid_mod is where the discretisation is decided: spacings, Fourier mode
  ! ordering, ghost cell layout, the index shortcuts (ip0, ij0, ikx0, iz_obmp)
  ! that the rest of the code indexes with, and the Simpson weights used for
  ! every parallel integral. None of that is checked at run time, and a wrong
  ! offset produces a plausible looking but wrong simulation.
  !
  ! The run is single rank, so local quantities must coincide with the totals.
  USE prec_const, ONLY: xp, PI
  USE basic,      ONLY: lu_in
  USE parallel,   ONLY: ppinit
  USE calculus,   ONLY: simpson_rule_z
  USE grid
  USE test_harness
  IMPLICIT NONE

  ! Small but representative: even Nx, Nz a multiple of 2 and at least 4 so the
  ! four point z stencil fits, more than one Hermite and Laguerre degree.
  INTEGER,  PARAMETER :: T_PMAX = 4, T_JMAX = 2
  INTEGER,  PARAMETER :: T_NX = 8, T_NY = 8, T_NZ = 16
  REAL(xp), PARAMETER :: T_LX = 60._xp, T_LY = 60._xp
  REAL(xp), PARAMETER :: T_NPOL = 1._xp, T_SHEAR = 0._xp
  INTEGER,  PARAMETER :: T_NA = 1, T_NHD = 4, T_NHDZ = 4

  REAL(xp) :: dkx, dky, dz, Lz, tol
  REAL(xp) :: want_r(64)
  COMPLEX(xp), ALLOCATABLE :: fz(:)
  COMPLEX(xp) :: intf
  INTEGER  :: i, ip, ij, iz, nk
  LOGICAL  :: ok

  CALL ppinit
  CALL start_suite('grid')
  tol = tol_of(1.0e4_xp)

  ! Feed the GRID namelist through a scratch file so nothing is left on disk.
  OPEN(UNIT=lu_in, STATUS='SCRATCH', ACTION='READWRITE', FORM='FORMATTED')
  WRITE(lu_in,'(A)')      '&GRID'
  WRITE(lu_in,'(A,I0)')   '  pmax = ', T_PMAX
  WRITE(lu_in,'(A,I0)')   '  jmax = ', T_JMAX
  WRITE(lu_in,'(A,I0)')   '  Nx   = ', T_NX
  WRITE(lu_in,'(A,F8.2)') '  Lx   = ', T_LX
  WRITE(lu_in,'(A,I0)')   '  Ny   = ', T_NY
  WRITE(lu_in,'(A,F8.2)') '  Ly   = ', T_LY
  WRITE(lu_in,'(A,I0)')   '  Nz   = ', T_NZ
  WRITE(lu_in,'(A)')      '  SG   = .false.'
  WRITE(lu_in,'(A)')      '  Nexc = 1'
  WRITE(lu_in,'(A)')      '/'
  REWIND(lu_in)
  CALL grid_readinputs
  CLOSE(lu_in)

  ! 'linear' keeps the ky distribution manual instead of handing it to FFTW,
  ! which keeps this test independent of the FFTW MPI plan machinery.
  CALL init_grids_data(T_NA, .FALSE., 'linear')
  CALL set_grids(T_SHEAR, T_NPOL, 'linear', T_NHD, T_NHDZ, .FALSE., T_NA)

  dky = 2._xp*PI/T_LY
  dkx = 2._xp*PI/T_LX
  Lz  = 2._xp*PI*T_NPOL
  dz  = Lz/REAL(T_NZ,xp)

  !------------------------------------------------------------- overall sizes
  CALL section('grid sizes')
  CALL check_int('total_na',  total_na,  T_NA)
  CALL check_int('total_np',  total_np,  T_PMAX+1)
  CALL check_int('total_nj',  total_nj,  T_JMAX+1)
  CALL check_int('total_nky', total_nky, T_NY/2+1)
  CALL check_int('total_nkx', total_nkx, T_NX)
  CALL check_int('total_nz',  total_nz,  T_NZ)

  CALL section('single rank: local sizes equal the totals')
  CALL check_int('local_np',  local_np,  total_np)
  CALL check_int('local_nj',  local_nj,  total_nj)
  CALL check_int('local_nky', local_nky, total_nky)
  CALL check_int('local_nkx', local_nkx, total_nkx)
  CALL check_int('local_nz',  local_nz,  total_nz)
  CALL check_int('no z offset',  local_nz_offset,  0)
  CALL check_int('no ky offset', local_nky_offset, 0)
  ! local_nx is the real space x extent, not a kx split (local_nkx is always
  ! Nx). FFTW sets it in a nonlinear run; a linear run builds no plan, so
  ! grid_mod has to set it itself. It cannot be left undefined because
  ! Setup_ExB_shear_flow sizes inv_ExB_NL_factor with it whatever the
  ! linearity. This run is linear and single rank, hence the whole range.
  CALL check_int('local_nx',     local_nx,         Nx)
  CALL check_int('no x offset',  local_nx_offset,  0)

  CALL section('ghost cell counts')
  ! Four points in p (p+-1 and p+-2 couplings), two in j (j+-1), four in z.
  CALL check_int('ngp',    ngp,    4)
  CALL check_int('ngj',    ngj,    2)
  CALL check_int('ngz',    ngz,    4)
  CALL check_int('ngp_o2', ngp_o2, ngp/2)
  CALL check_int('ngj_o2', ngj_o2, ngj/2)
  CALL check_int('ngz_o2', ngz_o2, ngz/2)
  CALL check_int('nzgrid is 1 without staggering', nzgrid, 1)

  !------------------------------------------------------------- Hermite grid
  CALL section('Hermite grid and its ghosts')
  CALL check_int('interior starts at p = 0',  parray(1+ngp_o2),         0)
  CALL check_int('interior ends at pmax',     parray(local_np+ngp_o2),  T_PMAX)
  CALL check_int('left ghost continues down', parray(1),               -ngp_o2)
  CALL check_int('right ghost continues up',  parray(local_np+ngp),     T_PMAX+ngp_o2)
  ok = .TRUE.
  DO ip = 2,local_np+ngp
    IF (parray(ip) - parray(ip-1) .NE. 1) ok = .FALSE.
  ENDDO
  CALL check('parray increases by deltap = 1', ok)
  CALL check_int('parray_full spans 0..pmax', parray_full(total_np), T_PMAX)

  CALL section('Hermite index shortcuts')
  CALL check_int('parray(ip0) = 0', parray(ip0), 0)
  CALL check_int('parray(ip1) = 1', parray(ip1), 1)
  CALL check_int('parray(ip2) = 2', parray(ip2), 2)
  CALL check_int('parray(ip3) = 3', parray(ip3), 3)
  CALL check_real('kroneck_p0 picks out p=0', SUM(kroneck_p0), 1._xp)
  CALL check_real('kroneck_p0 is at ip0',     kroneck_p0(ip0), 1._xp)
  CALL check_real('kroneck_p2 picks out p=2', SUM(kroneck_p2), 1._xp)
  CALL check_real('kroneck_p2 is at ip2',     kroneck_p2(ip2), 1._xp)

  !------------------------------------------------------------- Laguerre grid
  CALL section('Laguerre grid and its ghosts')
  CALL check_int('interior starts at j = 0', jarray(1+ngj_o2),        0)
  CALL check_int('interior ends at jmax',    jarray(local_nj+ngj_o2), T_JMAX)
  CALL check_int('left ghost is j = -1',     jarray(1),              -1)
  CALL check_int('right ghost is jmax+1',    jarray(local_nj+ngj),    T_JMAX+1)
  CALL check_int('jarray(ij0) = 0', jarray(ij0), 0)
  CALL check_int('jarray(ij1) = 1', jarray(ij1), 1)
  CALL check_real('kroneck_j0 picks out j=0', SUM(kroneck_j0), 1._xp)
  CALL check_real('kroneck_j0 is at ij0',     kroneck_j0(ij0), 1._xp)

  !------------------------------------------------------------- binormal grid
  CALL section('ky grid')
  CALL check_real('deltaky = 2pi/Ly', deltaky, dky, tol)
  nk = total_nky
  DO i = 1,nk
    want_r(i) = REAL(i-1,xp)*dky
  ENDDO
  CALL check_real_arr('kyarray_full is 0, dky, 2dky, ...', kyarray_full(1:nk), want_r(1:nk), tol)
  CALL check_real_arr('kyarray matches the full grid', kyarray(1:local_nky), want_r(1:nk), tol)
  CALL check_real('ky_min = deltaky',        ky_min, dky,                   tol)
  CALL check_real('ky_max = (Nky-1) deltaky',ky_max, REAL(nk-1,xp)*dky,     tol)
  CALL check_int ('iky0 points at ky = 0',   iky0, 1)
  CALL check_real('two_third_kymax',         two_third_kymax, 2._xp/3._xp*ky_max, tol)

  !------------------------------------------------------------- radial grid
  CALL section('kx grid uses the FFT mode ordering')
  ! Modes run 0, dkx, ... up to the Nyquist mode, then wrap to the negatives.
  ! The Nyquist entry itself is stored with a positive sign.
  CALL check_real('deltakx = 2pi/Lx', deltakx, dkx, tol)
  nk = total_nkx
  DO i = 1,nk
    IF (i .LE. nk/2+1) THEN
      want_r(i) = REAL(i-1,xp)*dkx
    ELSE
      want_r(i) = REAL(i-1-nk,xp)*dkx
    ENDIF
  ENDDO
  CALL check_real_arr('kxarray_full', kxarray_full(1:nk), want_r(1:nk), tol)
  CALL check_real_arr('kxarray0 matches the full grid', kxarray0(1:local_nkx), want_r(1:nk), tol)
  CALL check_real_arr('kxarray row for iky=1', kxarray(1,1:local_nkx), want_r(1:nk), tol)
  CALL check_real('kx_max is the Nyquist mode', kx_max,  REAL(nk/2,xp)*dkx,   tol)
  CALL check_real('kx_min',                     kx_min, -REAL(nk/2-1,xp)*dkx, tol)
  CALL check_int ('ikx0 points at kx = 0',      ikx0, 1)
  CALL check_real('two_third_kxmax', two_third_kxmax, 2._xp/3._xp*ABS(kx_min), tol)
  CALL check_real('inv_dkx', inv_dkx*deltakx, 1._xp, tol)

  CALL section('x grid')
  CALL check_real('deltax = Lx/Nx',       deltax, T_LX/REAL(T_NX,xp),  tol)
  CALL check_real('xarray starts at 0',   xarray(1), 0._xp,            tol)
  CALL check_real('xarray stops one cell short of Lx', &
                  xarray(T_NX), T_LX-T_LX/REAL(T_NX,xp), tol)

  !------------------------------------------------------------- parallel grid
  CALL section('z grid')
  CALL check_real('deltaz = 2 pi Npol / Nz', deltaz, dz, tol)
  CALL check_real('inv_deltaz', inv_deltaz*deltaz, 1._xp, tol)
  CALL check_real('z starts at -pi Npol', zarray_full(1), -Lz/2._xp, tol)
  CALL check_real('z stops one cell short of +pi Npol', &
                  zarray_full(T_NZ), Lz/2._xp-dz, tol)
  CALL check_int ('iz_obmp is the midplane index', iz_obmp, T_NZ/2+1)
  CALL check_real('z vanishes at the outboard midplane', zarray_full(iz_obmp), 0._xp, tol)

  CALL section('z ghost cells continue the grid periodically')
  ok = .TRUE.
  DO iz = 1,local_nz
    IF (ABS(zarray(iz+ngz_o2,1) - zarray_full(iz)) .GT. tol) ok = .FALSE.
  ENDDO
  CALL check('interior of zarray matches zarray_full', ok)
  CALL check_real('ghost below the first point',  zarray(ngz_o2,1), zarray_full(1)-dz, tol)
  CALL check_real('two ghosts below',             zarray(1,1),      zarray_full(1)-2._xp*dz, tol)
  CALL check_real('ghost above the last point',   zarray(local_nz+ngz_o2+1,1), &
                  zarray_full(T_NZ)+dz, tol)
  CALL check_real('two ghosts above',             zarray(local_nz+ngz,1), &
                  zarray_full(T_NZ)+2._xp*dz, tol)

  !------------------------------------------------------------- Simpson weights
  CALL section('composite Simpson weights')
  ! Periodic composite Simpson: the two end weights of the open rule merge, so
  ! the pattern is 2h/3, 4h/3, 2h/3, ... and the weights must sum to the domain
  ! length, otherwise every flux diagnostic is off by a constant factor.
  ok = .TRUE.
  DO iz = 1,local_nz
    IF (MODULO(iz,2) .EQ. 1) THEN
      IF (ABS(zweights_SR(iz) - 2._xp*dz/3._xp) .GT. tol) ok = .FALSE.
    ELSE
      IF (ABS(zweights_SR(iz) - 4._xp*dz/3._xp) .GT. tol) ok = .FALSE.
    ENDIF
  ENDDO
  CALL check('weights alternate 2h/3 and 4h/3', ok)
  CALL check_real('weights sum to the domain length', SUM(zweights_SR), Lz, tol)

  CALL section('simpson_rule_z integrates over z')
  ! On a periodic grid this rule is exact for any mode below the Nyquist one,
  ! so these are equalities and not approximations.
  ALLOCATE(fz(local_nz))
  fz = (1._xp, 0._xp)
  CALL simpson_rule_z(local_nz, zweights_SR, fz, intf)
  CALL check_cmpx('integral of 1 is the domain length', intf, CMPLX(Lz,0._xp,xp), tol)
  DO iz = 1,local_nz
    fz(iz) = CMPLX(COS(zarray_full(iz)), SIN(zarray_full(iz)), xp)
  ENDDO
  CALL simpson_rule_z(local_nz, zweights_SR, fz, intf)
  CALL check_cmpx('integral of exp(i z) vanishes', intf, (0._xp,0._xp), tol)
  DO iz = 1,local_nz
    fz(iz) = CMPLX(COS(zarray_full(iz))**2, 0._xp, xp)
  ENDDO
  CALL simpson_rule_z(local_nz, zweights_SR, fz, intf)
  CALL check_cmpx('integral of cos^2 is half the domain', intf, CMPLX(Lz/2._xp,0._xp,xp), tol)
  DEALLOCATE(fz)

  CALL finish_suite

END PROGRAM test_grid
