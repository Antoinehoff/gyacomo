PROGRAM test_ExB_shear_flow
  ! Exercises the background ExB shear machinery of ExB_shear_flow_mod.
  !
  ! The shear is implemented by letting the radial grid drift,
  !
  !   kxarray(iky,ikx) = kxarray0(ikx) + s(ky),      s(ky,t) = -ky gamma_E t,
  !
  ! and, once |s| grows past half a cell, re-indexing the fields by a whole
  ! number of cells so that s can be folded back into [-dkx/2, dkx/2]. That
  ! re-indexing walks an array held in FFT order
  !
  !   0, dkx, ... , kx_max, kx_min, ... , -dkx
  !
  ! so the wrap around offers a mode at the far end of the spectrum whenever the
  ! shift runs off the edge. Those modes do not exist and must be zeroed rather
  ! than copied, and on top of that the nonlinear run only keeps the band that
  ! survives the 2/3 rule. None of this is checked at run time: getting it wrong
  ! leaves a simulation that runs happily while quietly aliasing garbage across
  ! the spectrum at every shift event. It is therefore pinned down here against
  ! a reference built from the wavenumbers themselves, not from index
  ! arithmetic, so that the test cannot reproduce a sign error in the code.
  !
  ! The linearity is taken from the command line because neither
  ! Setup_ExB_shear_flow nor allocate_array deallocates first, so the setup can
  ! only run once per process. "make check" runs this binary once per mode.
  !
  !   usage: test_ExB_shear_flow [linear|nonlinear|off]
  !
  ! "off" is ExBrate = 0 and pins down that the whole module then stands aside.
  USE prec_const, ONLY: xp
  USE basic,      ONLY: lu_in, dt, increase_time
  USE parallel,   ONLY: ppinit
  USE grid,       ONLY: Nx, Ny, total_nkx, total_nky, local_nkx, local_nky, local_nx, &
                        deltakx, deltaky, inv_dkx, kx_max, kx_min, &
                        kxarray, kxarray0, kyarray
  USE fields,     ONLY: moments, phi, psi
  USE model,      ONLY: LINEARITY, ExBrate
  USE geometry,   ONLY: C_y, Cyq0_x0
  USE ExB_shear_flow
  USE test_harness
  IMPLICIT NONE

  ! Nx = 16 leaves room for the 2/3 band to be a proper subset of the grid, and
  ! Ny = 8 gives 5 ky rows, i.e. 5 different shifts in a single update.
  INTEGER,  PARAMETER :: T_PMAX = 2, T_JMAX = 1
  INTEGER,  PARAMETER :: T_NX = 16, T_NY = 8, T_NZ = 1
  REAL(xp), PARAMETER :: T_LX = 60._xp, T_LY = 60._xp
  INTEGER,  PARAMETER :: T_NA = 1
  ! Shift of the ky = deltaky row per time step, in units of deltakx. 1.2 gives
  ! the rows the jumps 0, -1, -2, -4, -5 and leaves every row with a nonzero
  ! residual, while keeping all of them well clear of the half cell rounding
  ! boundary where NINT could go either way.
  REAL(xp), PARAMETER :: T_SHIFT_PER_STEP = 1.2_xp
  ! The jumps the rows above must produce. Hard coded so that a change in the
  ! box or in the rate is caught here instead of silently weakening the test.
  INTEGER,  PARAMETER :: WANT_JUMP(5) = (/ 0, -1, -2, -4, -5 /)
  INTEGER,  PARAMETER :: N_RESTORE_STEPS = 7

  CHARACTER(len=32) :: mode
  REAL(xp) :: tol, rate, dkx, dky, want, kx_src, drift
  REAL(xp) :: sky_kept(64), kx_kept(64)
  REAL(xp), ALLOCATABLE :: kxarray_kept(:,:)
  INTEGER  :: nargs, iky, ikx, s, j, nkept, nzeroed
  LOGICAL  :: shearing, ok, ok_phi, ok_psi, ok_mom, ok_zero

  CALL ppinit
  nargs = COMMAND_ARGUMENT_COUNT()
  IF (nargs .GE. 1) THEN
    CALL GET_COMMAND_ARGUMENT(1, mode)
  ELSE
    mode = 'linear'
  ENDIF
  SELECT CASE (TRIM(mode))
  CASE ('linear', 'nonlinear'); shearing = .TRUE.
  CASE ('off');                 shearing = .FALSE.
  CASE DEFAULT
    WRITE(*,'(A)') 'unknown mode "'//TRIM(mode)//'"'
    CALL start_suite('ExB_shear_flow ['//TRIM(mode)//']')
    CALL check('mode is known to the test', .FALSE.)
    CALL finish_suite
    STOP
  END SELECT

  CALL start_suite('ExB_shear_flow ['//TRIM(mode)//']')
  tol = tol_of(1.0e4_xp)

  dkx = 2._xp*ACOS(-1._xp)/T_LX
  dky = 2._xp*ACOS(-1._xp)/T_LY
  ! basic_data is not called: it would go looking for an input file on disk and
  ! would swallow our command line argument as a file index. dt therefore keeps
  ! the default of basic_mod, and the rate is chosen against it so that the test
  ! stays valid if that default ever moves.
  IF (shearing) THEN
    rate = T_SHIFT_PER_STEP*dkx/(dky*dt)
  ELSE
    rate = 0._xp
  ENDIF

  CALL write_input_deck
  CALL readinputs
  CALL auxval          ! builds the grids, the operators, and calls Setup_ExB_shear_flow

  !--------------------------------------------------------------- Setup
  CALL section('the run was configured the way the test assumes')
  CALL check('LINEARITY was read back', TRIM(LINEARITY) .EQ. &
             TRIM(MERGE('nonlinear', 'linear   ', TRIM(mode) .EQ. 'nonlinear')))
  CALL check_int('single rank: local_nky is the total', local_nky, total_nky)
  CALL check_int('single rank: local_nkx is the total', local_nkx, total_nkx)
  CALL check_int('5 ky rows', total_nky, T_NY/2+1)
  CALL check_int('16 kx modes', total_nkx, T_NX)
  CALL check_real('z-pinch has C_y = 1',     C_y,     1._xp, tol)
  CALL check_real('z-pinch has Cyq0_x0 = 1', Cyq0_x0, 1._xp, tol)

  CALL section('Setup_ExB_shear_flow')
  CALL check_real('gamma_E = ExBrate C_y |Cyq0_x0/C_y|', &
                  gamma_E, ExBrate*C_y*ABS(Cyq0_x0/C_y), tol)
  CALL check('ExB is on exactly when gamma_E does not vanish', ExB .EQV. shearing)
  IF (shearing) THEN
    ! t0 is the time the ky = deltaky row needs to slide by one deltakx.
    CALL check_real('t0 gamma_E deltaky = deltakx', t0*gamma_E*deltaky, deltakx, tol)
    CALL check_real('inv_t0 inverts t0', inv_t0*t0, 1._xp, tol)
  ELSE
    CALL check_real('t0 vanishes without shear',     t0,     0._xp, tol)
    CALL check_real('inv_t0 vanishes without shear', inv_t0, 0._xp, tol)
  ENDIF

  CALL section('Setup_ExB_shear_flow array shapes')
  CALL check_int('sky_ExB is one per local ky',    SIZE(sky_ExB),      local_nky)
  CALL check_int('jump_ExB is one per local ky',   SIZE(jump_ExB),     local_nky)
  CALL check_int('shiftnow_ExB is one per local ky',SIZE(shiftnow_ExB),local_nky)
  ! One entry beyond the last ky: the nonlinear factor needs the row above.
  CALL check_int('sky_ExB_full spans all ky plus one', SIZE(sky_ExB_full), total_nky+1)
  CALL check_int('ExB_NL_factor is (Nx,local_nky) 1', SIZE(ExB_NL_factor,1), Nx)
  CALL check_int('ExB_NL_factor is (Nx,local_nky) 2', SIZE(ExB_NL_factor,2), local_nky)
  CALL check_int('inv_ExB_NL_factor is (Ny/2+1,local_nx) 1', SIZE(inv_ExB_NL_factor,1), Ny/2+1)
  CALL check_int('inv_ExB_NL_factor is (Ny/2+1,local_nx) 2', SIZE(inv_ExB_NL_factor,2), local_nx)

  CALL section('the initial state carries no shift')
  CALL check_real_arr('sky_ExB starts at zero',      sky_ExB,      SPREAD(0._xp,1,local_nky), tol)
  CALL check_real_arr('sky_ExB_full starts at zero', sky_ExB_full, SPREAD(0._xp,1,total_nky+1), tol)
  CALL check('jump_ExB starts at zero',      ALL(jump_ExB .EQ. 0))
  CALL check('shiftnow_ExB starts unset',    ALL(.NOT. shiftnow_ExB))
  CALL check('ExB_NL_factor starts at one',     ALL(ExB_NL_factor     .EQ. (1._xp,0._xp)))
  CALL check('inv_ExB_NL_factor starts at one', ALL(inv_ExB_NL_factor .EQ. (1._xp,0._xp)))
  CALL check_real_arr('kxarray starts on kxarray0', kxarray(1,:), kxarray0, tol)

  CALL section('the band the shift is allowed to keep')
  IF (TRIM(mode) .EQ. 'nonlinear') THEN
    ! A nonlinear run may only keep what survives the 2/3 dealiasing.
    CALL check_real('maximal_kx is the 2/3 rule', maximal_kx, 2._xp/3._xp*kx_max-deltakx, tol)
    CALL check_real('minimal_kx is the 2/3 rule', minimal_kx, 2._xp/3._xp*kx_min+deltakx, tol)
    CALL check('the band is a proper subset of the grid', &
               (maximal_kx .LT. kx_max) .AND. (minimal_kx .GT. kx_min))
  ELSE
    CALL check_real('maximal_kx is kx_max', maximal_kx, kx_max, tol)
    CALL check_real('minimal_kx is kx_min', minimal_kx, kx_min, tol)
  ENDIF

  IF (.NOT. shearing) THEN
    !------------------------------------------------- ExBrate = 0 stands aside
    CALL section('without shear an update changes nothing')
    CALL tag_fields
    CALL Update_ExB_shear_flow(-1)
    CALL check_real_arr('sky_ExB stays at zero', sky_ExB, SPREAD(0._xp,1,local_nky), tol)
    CALL check('jump_ExB stays at zero', ALL(jump_ExB .EQ. 0))
    CALL check_real_arr('kxarray stays on kxarray0', kxarray(1,:), kxarray0, tol)
    ok = .TRUE.
    DO iky = 1,local_nky
      DO ikx = 1,total_nkx
        IF (phi(iky,ikx,1) .NE. tag(iky,ikx,0)) ok = .FALSE.
        IF (psi(iky,ikx,1) .NE. tag(iky,ikx,1)) ok = .FALSE.
      ENDDO
    ENDDO
    CALL check('the fields are left untouched', ok)
  ELSE

    !--------------------------------------------------------------- One shift
    CALL tag_fields
    CALL Update_ExB_shear_flow(-1)
    CALL increase_time

    CALL section('one update produces the expected jump per ky row')
    ok = .TRUE.
    DO iky = 1,local_nky
      ! What the drift of this row asks for, computed from the wavenumbers.
      j = NINT(-kyarray(iky)*gamma_E*dt*inv_dkx)
      IF (j .NE. WANT_JUMP(iky)) ok = .FALSE.
      CALL check_int('jump_ExB row '//istr(iky), jump_ExB(iky), WANT_JUMP(iky))
    ENDDO
    CALL check('the drift of each row asks for those same jumps', ok)
    CALL check('every row that jumped had its flag cleared again', &
               ALL(.NOT. shiftnow_ExB))

    CALL section('the residual shift is folded into half a cell')
    ok = .TRUE.
    DO iky = 1,local_nky
      ! s is the raw drift minus the whole cells that the re-indexing absorbed.
      want = -kyarray(iky)*gamma_E*dt - REAL(WANT_JUMP(iky),xp)*deltakx
      IF (ABS(sky_ExB(iky)-want) .GT. tol*deltakx) ok = .FALSE.
      IF (ABS(sky_ExB(iky)) .GT. 0.5_xp*deltakx*(1._xp+tol)) ok = .FALSE.
    ENDDO
    CALL check('sky_ExB is the drift less the whole cells, and is under dkx/2', ok)
    CALL check_real_arr('sky_ExB_full agrees with sky_ExB on the shared rows', &
                        sky_ExB_full(1:total_nky), sky_ExB(1:local_nky), tol)

    CALL section('the kx grid follows the residual shift')
    ok = .TRUE.
    DO iky = 1,local_nky
      DO ikx = 1,total_nkx
        IF (ABS(kxarray(iky,ikx)-(kxarray0(ikx)+sky_ExB(iky))) .GT. tol*deltakx) ok = .FALSE.
      ENDDO
    ENDDO
    CALL check('kxarray(iky,ikx) = kxarray0(ikx) + sky_ExB(iky)', ok)

    !---------------------------------------------------------------- The remap
    ! The reference. The shear carries the mode sitting at kx0 = v over to
    ! kx0 = v + J dkx, so the slot d ends up holding the mode whose unshifted
    ! wavenumber is kx_src = kxarray0(d) - J dkx. That mode is real data only if
    ! the grid actually has a mode there, which is asked of the wavenumbers by
    ! slot_of_kx rather than of the index arithmetic, and only if it sits inside
    ! the band the run keeps. Everything else is a wrap around artefact and has
    ! to be zeroed.
    !
    ! A row whose jump is zero is a separate case: it is skipped outright, so it
    ! keeps whatever it held, band or no band. The band is a criterion for what a
    ! shift may carry, not a standing filter, and the dealiasing of the fields
    ! themselves belongs to the nonlinear bracket instead.
    CALL section('the remap relabels the surviving modes and zeroes the rest')
    ok_phi = .TRUE.; ok_psi = .TRUE.; ok_mom = .TRUE.; ok_zero = .TRUE.
    ok = .TRUE.
    nkept = 0; nzeroed = 0
    DO iky = 1,local_nky
      IF (WANT_JUMP(iky) .EQ. 0) THEN
        DO ikx = 1,total_nkx
          IF (phi(iky,ikx,1)             .NE. tag(iky,ikx,0)) ok = .FALSE.
          IF (psi(iky,ikx,1)             .NE. tag(iky,ikx,1)) ok = .FALSE.
          IF (moments(1,1,1,iky,ikx,1,1) .NE. tag(iky,ikx,2)) ok = .FALSE.
        ENDDO
        CYCLE
      ENDIF
      DO ikx = 1,total_nkx
        kx_src = kxarray0(ikx) - REAL(WANT_JUMP(iky),xp)*deltakx
        s = slot_of_kx(kx_src)
        IF ((s .GT. 0) .AND. (kx_src .LE. maximal_kx+tol) .AND. &
                             (kx_src .GE. minimal_kx-tol)) THEN
          nkept = nkept + 1
          IF (phi(iky,ikx,1)            .NE. tag(iky,s,0)) ok_phi = .FALSE.
          IF (psi(iky,ikx,1)            .NE. tag(iky,s,1)) ok_psi = .FALSE.
          IF (moments(1,1,1,iky,ikx,1,1).NE. tag(iky,s,2)) ok_mom = .FALSE.
        ELSE
          nzeroed = nzeroed + 1
          IF (phi(iky,ikx,1)             .NE. (0._xp,0._xp)) ok_zero = .FALSE.
          IF (psi(iky,ikx,1)             .NE. (0._xp,0._xp)) ok_zero = .FALSE.
          IF (moments(1,1,1,iky,ikx,1,1) .NE. (0._xp,0._xp)) ok_zero = .FALSE.
        ENDIF
      ENDDO
    ENDDO
    CALL check('a row that does not jump is left alone',         ok)
    CALL check('phi carries the right mode into every slot',     ok_phi)
    CALL check('psi carries the right mode into every slot',     ok_psi)
    CALL check('moments carry the right mode into every slot',   ok_mom)
    CALL check('the slots with no source are zeroed',            ok_zero)
    ! Guard against a reference that kept or dropped everything and so proved
    ! nothing: the rows above must exercise both outcomes.
    CALL check('the reference exercises surviving modes', nkept   .GT. 0)
    CALL check('the reference exercises dropped modes',   nzeroed .GT. 0)

    CALL section('the remap conserves the wavenumber of the modes it keeps')
    ! Read back where each mode landed from its own tag and check that the shear
    ! moved it by exactly -ky gamma_E dt in kx. This is the statement the whole
    ! algorithm exists to make, and it does not mention an index anywhere.
    ok = .TRUE.
    DO iky = 1,local_nky
      drift = -kyarray(iky)*gamma_E*dt
      DO ikx = 1,total_nkx
        IF (phi(iky,ikx,1) .EQ. (0._xp,0._xp)) CYCLE
        s = NINT(AIMAG(phi(iky,ikx,1)))          ! slot the mode started in
        IF (ABS(kxarray(iky,ikx)-(kxarray0(s)+drift)) .GT. tol*deltakx) ok = .FALSE.
      ENDDO
    ENDDO
    CALL check('every surviving mode moved by -ky gamma_E dt', ok)

    CALL section('no mode is duplicated by the in place copy')
    ! The copy runs in place, so a loop that walks the array the wrong way would
    ! smear one mode over several slots instead of moving each one once.
    ok = .TRUE.
    DO iky = 1,local_nky
      nkept = 0
      DO ikx = 1,total_nkx
        IF (phi(iky,ikx,1) .EQ. (0._xp,0._xp)) CYCLE
        nkept = nkept + 1
        kx_kept(nkept) = AIMAG(phi(iky,ikx,1))
      ENDDO
      DO s = 1,nkept
        DO j = s+1,nkept
          IF (kx_kept(s) .EQ. kx_kept(j)) ok = .FALSE.
        ENDDO
      ENDDO
    ENDDO
    CALL check('each surviving mode appears exactly once per ky row', ok)

    !-------------------------------------------------------------- The restart
    CALL section('Restore_ExB_shear_flow rebuilds the accumulated shift')
    ! A checkpoint stores neither s(ky) nor the whole cells already absorbed, so
    ! the restart rebuilds the residual from the time alone. Step the incremental
    ! path forward, then ask the closed form to reproduce where it landed.
    DO s = 1,N_RESTORE_STEPS
      CALL Update_ExB_shear_flow(-1)
      CALL increase_time
    ENDDO
    sky_kept(1:local_nky) = sky_ExB
    ! Every row, not just the first: ky = 0 never moves and would not notice a
    ! restore that rebuilt nothing at all.
    ALLOCATE(kxarray_kept(local_nky,total_nkx))
    kxarray_kept = kxarray
    CALL Restore_ExB_shear_flow
    CALL check_real_arr('the rebuilt sky_ExB matches the accumulated one', &
                        sky_ExB, sky_kept(1:local_nky), tol)
    CALL check('the rebuilt kx grid matches the accumulated one, on every ky row', &
               ALL(ABS(kxarray-kxarray_kept) .LE. tol*deltakx))
    CALL check('the ky rows really had moved, so the check above bites', &
               ANY(ABS(kxarray_kept-SPREAD(kxarray0,1,local_nky)) .GT. tol*deltakx))
    DEALLOCATE(kxarray_kept)
    CALL check('the restore leaves no shift pending', &
               ALL(jump_ExB .EQ. 0) .AND. ALL(.NOT. shiftnow_ExB))
    ok = .TRUE.
    DO iky = 1,local_nky
      IF (ABS(sky_ExB(iky)) .GT. 0.5_xp*deltakx*(1._xp+tol)) ok = .FALSE.
    ENDDO
    CALL check('the rebuilt shift is folded into half a cell', ok)

  ENDIF

  CALL finish_suite

CONTAINS

  !> Unique marker for the mode held in slot (iky,ikx). The imaginary part is
  !! the slot index, so a field can be asked afterwards which mode it ended up
  !! holding. which_ separates phi, psi and moments so that a copy mixing the
  !! three of them up cannot pass.
  COMPLEX(xp) FUNCTION tag(iky_, ikx_, which_)
    INTEGER, INTENT(IN) :: iky_, ikx_, which_
    tag = CMPLX(REAL(iky_,xp)+100._xp*REAL(which_,xp), REAL(ikx_,xp), xp)
  END FUNCTION tag

  !> Mark every slot of every field so the remap can be followed.
  SUBROUTINE tag_fields
    INTEGER :: iy, ix
    DO iy = 1,local_nky
      DO ix = 1,total_nkx
        phi(iy,ix,:)             = tag(iy,ix,0)
        psi(iy,ix,:)             = tag(iy,ix,1)
        moments(:,:,:,iy,ix,:,:) = tag(iy,ix,2)
      ENDDO
    ENDDO
  END SUBROUTINE tag_fields

  !> Slot holding the mode of wavenumber kx_, or 0 when the grid has no such
  !! mode. Searching the wavenumbers themselves keeps the reference independent
  !! of the index arithmetic and of the FFT wrap around that the code under
  !! test relies on, which is what the wrap around bug would hide behind.
  INTEGER FUNCTION slot_of_kx(kx_)
    REAL(xp), INTENT(IN) :: kx_
    INTEGER :: is
    slot_of_kx = 0
    DO is = 1,total_nkx
      IF (ABS(kxarray0(is)-kx_) .LT. 0.5_xp*deltakx) THEN
        slot_of_kx = is
        RETURN
      ENDIF
    ENDDO
  END FUNCTION slot_of_kx

  !> Small integer to string, for assertion labels.
  CHARACTER(len=8) FUNCTION istr(i_)
    INTEGER, INTENT(IN) :: i_
    WRITE(istr,'(I0)') i_
  END FUNCTION istr

  !> Feed the namelists through a scratch file so nothing is left on disk. Only
  !! the entries this test depends on are set; the rest keep their defaults.
  SUBROUTINE write_input_deck
    OPEN(UNIT=lu_in, STATUS='SCRATCH', ACTION='READWRITE', FORM='FORMATTED')
    WRITE(lu_in,'(A)')       '&GRID'
    WRITE(lu_in,'(A,I0)')    '  pmax = ', T_PMAX
    WRITE(lu_in,'(A,I0)')    '  jmax = ', T_JMAX
    WRITE(lu_in,'(A,I0)')    '  Nx   = ', T_NX
    WRITE(lu_in,'(A,F8.2)')  '  Lx   = ', T_LX
    WRITE(lu_in,'(A,I0)')    '  Ny   = ', T_NY
    WRITE(lu_in,'(A,F8.2)')  '  Ly   = ', T_LY
    WRITE(lu_in,'(A,I0)')    '  Nz   = ', T_NZ
    WRITE(lu_in,'(A)')       '  SG   = .false.'
    WRITE(lu_in,'(A)')       '  Nexc = 1'
    WRITE(lu_in,'(A)')       '/'
    ! z-pinch fixes C_y = Cyq0_x0 = 1, so gamma_E is exactly ExBrate and the
    ! jumps below can be predicted from the box alone.
    WRITE(lu_in,'(A)')       '&GEOMETRY'
    WRITE(lu_in,'(A)')       "  geom = 'z-pinch'"
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&DIAGNOSTICS'
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&MODEL'
    WRITE(lu_in,'(A)')       "  LINEARITY = '"//TRIM(MERGE('nonlinear', 'linear   ', &
                                              TRIM(mode) .EQ. 'nonlinear'))//"'"
    WRITE(lu_in,'(A,I0)')    '  Na      = ', T_NA
    WRITE(lu_in,'(A)')       '  ADIAB_E = .true.'
    WRITE(lu_in,'(A)')       '  beta    = 0.0'
    WRITE(lu_in,'(A)')       '  nu      = 0.0'
    WRITE(lu_in,'(A,ES23.15)')'  ExBrate = ', rate
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&CLOSURE'
    WRITE(lu_in,'(A)')       "  hierarchy_closure = 'truncation'"
    WRITE(lu_in,'(A)')       "  nonlinear_closure = 'anti_laguerre_aliasing'"
    WRITE(lu_in,'(A)')       '  nmax = -1'
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&SPECIES'
    WRITE(lu_in,'(A)')       "  name_  = 'ions'"
    WRITE(lu_in,'(A)')       '  tau_   = 1.0'
    WRITE(lu_in,'(A)')       '  sigma_ = 1.0'
    WRITE(lu_in,'(A)')       '  q_     = 1.0'
    WRITE(lu_in,'(A)')       '  k_N_   = 2.0'
    WRITE(lu_in,'(A)')       '  k_T_   = 0.4'
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

END PROGRAM test_ExB_shear_flow
