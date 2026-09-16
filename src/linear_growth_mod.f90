MODULE linear_growth
  ! Measures the complex eigenvalue of a linear run while it runs, and keeps the
  ! mode amplitude inside a fixed band so that a fast growing mode cannot
  ! overflow.
  !
  ! A linear run has no saturation mechanism: the fastest growing mode simply
  ! grows as exp(gamma t) until it reaches the largest representable number and
  ! the run dies on a NaN. The only way to avoid that used to be to guess a tmax
  ! small enough for the growth rate one is trying to measure, which is circular
  ! and has to be redone for every ky. Renormalising the state leaves the linear
  ! solution unchanged (the equations are homogeneous in the fields) and removes
  ! the ceiling altogether, so the run can instead be stopped on the thing one
  ! actually wants, a converged eigenvalue.
  !
  ! The eigenvalue is measured per ky from the projection of the field on its
  ! own value one measurement earlier,
  !
  !   ratio  = <phi_now, phi_prev> / <phi_prev, phi_prev>,   <a,b> = sum conjg(b) a
  !   lambda = log(ratio) / dt_measure
  !
  ! summed over kx and z. With the code's exp(-i omega t + gamma t) convention
  ! that gives gamma = Re(lambda) and omega = -Im(lambda). The projection is a
  ! ratio of two sums over the same points, so no Jacobian weight is needed: any
  ! weight that does not depend on time cancels.
  !
  ! Note that phi_prev is rescaled along with the state. That is what makes the
  ! measurement blind to the renormalisation: both sides of the ratio carry the
  ! same factor and it divides out, so no correction term is ever needed.
  USE prec_const, ONLY: xp, mpi_xp_r, mpi_xp_c
  IMPLICIT NONE
  PRIVATE

  ! Inputs (LINEAR_GROWTH namelist)
  LOGICAL,  PUBLIC, PROTECTED :: LIN_GROWTH  = .FALSE.    ! activates the whole module
  REAL(xp), PUBLIC, PROTECTED :: lin_tol     = 1.0e-3_xp  ! rel. spread of gamma and omega to call it converged
  INTEGER,  PUBLIC, PROTECTED :: lin_nconv   = 5          ! the last lin_nconv+1 measurements must agree within lin_tol
  REAL(xp), PUBLIC, PROTECTED :: lin_ceiling = 1.0e3_xp   ! rescale above this amplitude
  REAL(xp), PUBLIC, PROTECTED :: lin_floor   = 1.0e-3_xp  ! rescale below this amplitude

  ! Measured spectrum, held in full (total_nky) on every process
  REAL(xp), PUBLIC, PROTECTED, ALLOCATABLE :: gamma_ky(:), omega_ky(:)
  ! Cumulative renormalisation: the unscaled field is phi*EXP(lin_lnscale)
  REAL(xp), PUBLIC, PROTECTED :: lin_lnscale   = 0._xp
  LOGICAL,  PUBLIC, PROTECTED :: lin_converged = .FALSE.

  ! Internal state
  ! A component of the eigenvalue smaller than this fraction of |lambda| is judged
  ! against that fraction instead of against itself. Otherwise a mode with no
  ! real frequency would compare round off noise to zero and never converge.
  REAL(xp), PARAMETER :: REL_FLOOR = 1.0e-3_xp

  COMPLEX(xp), ALLOCATABLE :: phi_prev(:,:,:)  ! field at the previous measurement
  REAL(xp),    ALLOCATABLE :: g_hist(:,:)      ! last lin_nconv+1 gammas (ky, slot), ring buffer
  REAL(xp),    ALLOCATABLE :: w_hist(:,:)      ! same for omega
  INTEGER,     ALLOCATABLE :: nfill(:)         ! number of measurements taken so far per ky
  REAL(xp),    ALLOCATABLE :: drift(:)         ! relative spread over the window (displayed)
  LOGICAL,     ALLOCATABLE :: settled(:)       ! the window is full and within lin_tol
  LOGICAL,     ALLOCATABLE :: measured(:)      ! ky carrying a usable signal
  REAL(xp) :: t_prev      = 0._xp
  LOGICAL  :: has_prev    = .FALSE. ! a previous field is available to compare against
  LOGICAL  :: initialized = .FALSE.
  LOGICAL  :: warned_wrap = .FALSE. ! the branch cut warning is only worth giving once

  PUBLIC :: lin_readinputs, lin_outputinputs, lin_init, &
            lin_measure, lin_rescale, lin_print

CONTAINS

  SUBROUTINE lin_readinputs
    ! Read the LINEAR_GROWTH namelist.
    !
    ! The group is optional: every input file written before this module existed
    ! must keep working untouched. REWIND makes the read independent of where the
    ! group sits in the file, unlike the other namelists which are positional.
    USE basic, ONLY: lu_in, speak
    USE model, ONLY: LINEARITY
    IMPLICIT NONE
    INTEGER :: iost
    NAMELIST /LINEAR_GROWTH/ LIN_GROWTH, lin_tol, lin_nconv, lin_ceiling, lin_floor

    REWIND(lu_in)
    READ(lu_in, nml=LINEAR_GROWTH, iostat=iost)
    ! A negative status is the end of the file, i.e. no group, which is fine. A
    ! positive one is a malformed group (a misspelt entry, say), and carrying on
    ! would silently run without the feature the user asked for.
    IF (iost .GT. 0) ERROR STOP '>> ERROR << could not parse the LINEAR_GROWTH namelist'

    IF (LIN_GROWTH) THEN
      IF (lin_tol .LE. 0._xp) ERROR STOP '>> ERROR << lin_tol must be positive'
      IF (lin_nconv .LT. 1)   ERROR STOP '>> ERROR << lin_nconv must be at least 1'
      IF ((lin_floor .LE. 0._xp) .OR. (lin_floor .GE. lin_ceiling)) &
        ERROR STOP '>> ERROR << need 0 < lin_floor < lin_ceiling'
    ENDIF

    IF (LIN_GROWTH .AND. (LINEARITY .NE. 'linear')) THEN
      CALL speak('LIN_GROWTH is only meaningful for a linear run -> disabled',0)
      LIN_GROWTH = .FALSE.
    ENDIF
  END SUBROUTINE lin_readinputs

  SUBROUTINE lin_outputinputs(fid)
    ! Write the input parameters to the results_xx.h5 file
    USE futils, ONLY: attach, creatd
    IMPLICIT NONE
    INTEGER, INTENT(in) :: fid
    CHARACTER(len=256)  :: str
    WRITE(str,'(a)') '/data/input/linear_growth'
    CALL creatd(fid, 0,(/0/),TRIM(str),'Linear growth Input')
    CALL attach(fid, TRIM(str),  "LIN_GROWTH",  LIN_GROWTH)
    CALL attach(fid, TRIM(str),     "lin_tol",     lin_tol)
    CALL attach(fid, TRIM(str),   "lin_nconv",   lin_nconv)
    CALL attach(fid, TRIM(str), "lin_ceiling", lin_ceiling)
    CALL attach(fid, TRIM(str),   "lin_floor",   lin_floor)
  END SUBROUTINE lin_outputinputs

  SUBROUTINE lin_init
    ! Allocate the module state. Must run after the grids are set.
    USE grid, ONLY: local_nky, local_nkx, local_nz, ngz, total_nky
    IMPLICIT NONE
    IF (.NOT. LIN_GROWTH) RETURN
    IF (initialized)      RETURN
    ALLOCATE(phi_prev(local_nky,local_nkx,local_nz+ngz))
    ALLOCATE(gamma_ky(total_nky), omega_ky(total_nky))
    ALLOCATE(g_hist(total_nky,0:lin_nconv), w_hist(total_nky,0:lin_nconv))
    ALLOCATE(nfill(total_nky), drift(total_nky))
    ALLOCATE(settled(total_nky), measured(total_nky))
    phi_prev  = (0._xp,0._xp)
    gamma_ky  = 0._xp
    omega_ky  = 0._xp
    g_hist    = 0._xp
    w_hist    = 0._xp
    nfill     = 0
    drift     = 0._xp
    settled   = .FALSE.
    measured  = .FALSE.
    initialized = .TRUE.
  END SUBROUTINE lin_init

  SUBROUTINE lin_measure
    ! Update gamma_ky and omega_ky from the field, and test for convergence.
    !
    ! Collective: every process must call this. The reduction over z is what
    ! completes the parallel integral, and the one over ky is what gives every
    ! process the whole spectrum so that the convergence decision below is
    ! reached identically everywhere.
    USE basic,    ONLY: time, speak
    USE fields,   ONLY: phi
    USE grid,     ONLY: local_nky, local_nkx, local_nz, ngz_o2, &
                        local_nky_offset, total_nky, contains_ky0, iky0
    USE parallel, ONLY: comm_z, comm_ky, num_procs_z, num_procs_ky
    USE prec_const, ONLY: PI
    USE mpi
    IMPLICIT NONE
    COMPLEX(xp) :: num_loc(local_nky), den_loc(local_nky)
    COMPLEX(xp) :: buf(local_nky), ratio, lambda
    REAL(xp)    :: g_tot(total_nky), w_tot(total_nky), buf_r(total_nky)
    REAL(xp)    :: dt_m, gam, omg, lam_abs, rel_g, rel_w
    LOGICAL     :: seen(total_nky), buf_l(total_nky)
    INTEGER     :: iky, ikx, iz, izi, ikyg, ierr_, nwin

    IF (.NOT. LIN_GROWTH) RETURN
    CALL lin_init

    ! The first call has nothing to compare against: just record the field.
    IF (.NOT. has_prev) THEN
      phi_prev = phi
      t_prev   = time
      has_prev = .TRUE.
      RETURN
    ENDIF

    dt_m = time - t_prev
    IF (dt_m .LE. 0._xp) RETURN

    ! Projection of the field on its previous value, over the interior z points.
    num_loc = (0._xp,0._xp)
    den_loc = (0._xp,0._xp)
    DO iz = 1,local_nz
      izi = iz+ngz_o2
      DO ikx = 1,local_nkx
        DO iky = 1,local_nky
          num_loc(iky) = num_loc(iky) + CONJG(phi_prev(iky,ikx,izi))*phi(iky,ikx,izi)
          den_loc(iky) = den_loc(iky) + CONJG(phi_prev(iky,ikx,izi))*phi_prev(iky,ikx,izi)
        ENDDO
      ENDDO
    ENDDO

    ! Complete the z integral. Only over z: phi is replicated along p by
    ! manual_3D_bcast, so reducing over p would multiply the sum by num_procs_p,
    ! and kx is not distributed.
    IF (num_procs_z .GT. 1) THEN
      CALL MPI_ALLREDUCE(num_loc, buf, local_nky, mpi_xp_c, MPI_SUM, comm_z, ierr_)
      num_loc = buf
      CALL MPI_ALLREDUCE(den_loc, buf, local_nky, mpi_xp_c, MPI_SUM, comm_z, ierr_)
      den_loc = buf
    ENDIF

    ! Place the local ky rows into a full spectrum, zero elsewhere, so that a
    ! plain sum over the ky communicator assembles it.
    g_tot = 0._xp
    w_tot = 0._xp
    seen  = .FALSE.
    DO iky = 1,local_nky
      ! ky = 0 is not an eigenmode of anything here and gives 0/0.
      IF (contains_ky0 .AND. (iky .EQ. iky0)) CYCLE
      ! A mode that never got any amplitude carries no eigenvalue.
      IF (ABS(den_loc(iky)) .LE. TINY(1._xp)) CYCLE
      ratio  = num_loc(iky)/den_loc(iky)
      IF (ABS(ratio) .LE. TINY(1._xp)) CYCLE
      lambda = LOG(ratio)/dt_m
      gam    =  REAL(lambda,xp)
      omg    = -AIMAG(lambda)
      ! log() of a complex number is only unique up to 2 i pi, so a mode turning
      ! by more than pi between two measurements comes back with an aliased
      ! frequency (the growth rate is unaffected). There is no way to tell from
      ! the ratio alone how many turns were missed, so the only cure is a shorter
      ! diagnostic interval. Warn at 0.9 pi: close enough to the limit that a
      ! healthy run stays quiet, far enough to be said before the answer is wrong.
      IF ((ABS(AIMAG(LOG(ratio))) .GT. 0.9_xp*PI) .AND. (.NOT. warned_wrap)) THEN
        CALL speak('LIN_GROWTH: the mode turns by nearly pi between two '//&
                   'measurements, omega may alias -> reduce dtsave_0d',0)
        warned_wrap = .TRUE.
      ENDIF
      ikyg        = iky + local_nky_offset
      g_tot(ikyg) = gam
      w_tot(ikyg) = omg
      seen(ikyg)  = .TRUE.
    ENDDO

    IF (num_procs_ky .GT. 1) THEN
      CALL MPI_ALLREDUCE(g_tot, buf_r, total_nky, mpi_xp_r, MPI_SUM, comm_ky, ierr_)
      g_tot = buf_r
      CALL MPI_ALLREDUCE(w_tot, buf_r, total_nky, mpi_xp_r, MPI_SUM, comm_ky, ierr_)
      w_tot = buf_r
      CALL MPI_ALLREDUCE(seen, buf_l, total_nky, MPI_LOGICAL, MPI_LOR, comm_ky, ierr_)
      seen = buf_l
    ENDIF

    ! Convergence: a ky is converged once the last lin_nconv+1 values of both its
    ! growth rate and its frequency agree to within lin_tol, relative to their
    ! own size. Testing the spread over a window rather than the change between
    ! two neighbours matters: a slow, steady drift gives small consecutive
    ! changes that all pass while the value is still walking away, and the run
    ! would stop with an error several times lin_tol.
    DO ikyg = 1,total_nky
      IF (.NOT. seen(ikyg)) CYCLE
      g_hist(ikyg, MOD(nfill(ikyg), lin_nconv+1)) = g_tot(ikyg)
      w_hist(ikyg, MOD(nfill(ikyg), lin_nconv+1)) = w_tot(ikyg)
      nfill(ikyg) = nfill(ikyg)+1
      nwin    = MIN(nfill(ikyg), lin_nconv+1)
      lam_abs = HYPOT(g_tot(ikyg), w_tot(ikyg))
      rel_g   = (MAXVAL(g_hist(ikyg,0:nwin-1)) - MINVAL(g_hist(ikyg,0:nwin-1))) &
              / MAX(ABS(g_tot(ikyg)), REL_FLOOR*lam_abs, TINY(1._xp))
      rel_w   = (MAXVAL(w_hist(ikyg,0:nwin-1)) - MINVAL(w_hist(ikyg,0:nwin-1))) &
              / MAX(ABS(w_tot(ikyg)), REL_FLOOR*lam_abs, TINY(1._xp))
      ! A single value has no spread yet; show it as unconverged, not as zero.
      IF (nwin .LT. 2) THEN
        drift(ikyg) = 1._xp
      ELSE
        drift(ikyg) = MAX(rel_g, rel_w)
      ENDIF
      settled(ikyg)  = (nwin .EQ. lin_nconv+1) .AND. (drift(ikyg) .LT. lin_tol)
      measured(ikyg) = .TRUE.
    ENDDO

    gamma_ky = g_tot
    omega_ky = w_tot

    lin_converged = ANY(measured) .AND. ALL(settled .OR. .NOT. measured)

    phi_prev = phi
    t_prev   = time
  END SUBROUTINE lin_measure

  SUBROUTINE lin_rescale
    ! Renormalise the state if its amplitude has left [lin_floor, lin_ceiling].
    !
    ! The factor is a single global real number so that the relative amplitude of
    ! the different ky is untouched. Rescaling up as well as down matters: a
    ! damped branch would otherwise decay into the denormals and stop being
    ! measurable, which is the same failure as an overflow at the other end.
    USE basic,  ONLY: speak, str
    USE fields, ONLY: moments, phi, psi
    USE grid,   ONLY: local_nky, local_nkx, local_nz, ngz, ngz_o2
    USE time_integration, ONLY: ntimelevel
    USE mpi
    IMPLICIT NONE
    REAL(xp) :: amp_loc, amp, factor
    INTEGER  :: iz, izi, ikx, iky, ierr_

    IF (.NOT. LIN_GROWTH) RETURN
    CALL lin_init

    ! Measure on the interior points only: the ghosts are a copy of a neighbour
    ! and are refreshed from it anyway.
    amp_loc = 0._xp
    DO iz = 1,local_nz
      izi = iz+ngz_o2
      DO ikx = 1,local_nkx
        DO iky = 1,local_nky
          amp_loc = MAX(amp_loc, ABS(phi(iky,ikx,izi)))
        ENDDO
      ENDDO
    ENDDO
    CALL MPI_ALLREDUCE(amp_loc, amp, 1, mpi_xp_r, MPI_MAX, MPI_COMM_WORLD, ierr_)

    IF (amp .LE. TINY(1._xp))                    RETURN ! nothing to rescale
    IF ((amp .LT. lin_ceiling) .AND. (amp .GT. lin_floor)) RETURN

    factor = 1._xp/amp

    ! Every time level, not just the one carrying the state between steps: the
    ! Runge-Kutta stages are all read from within a step and leaving them at the
    ! old scale would mix two normalisations.
    moments(:,:,:,:,:,:,1:ntimelevel) = factor*moments(:,:,:,:,:,:,1:ntimelevel)
    phi = factor*phi
    psi = factor*psi
    ! The reference of the measurement moves with the state, so the projection
    ! ratio, and hence the eigenvalue, is unaffected by this call.
    IF (has_prev) phi_prev = factor*phi_prev

    ! Undo factor for post-processing: the unscaled field is phi*EXP(lin_lnscale).
    lin_lnscale = lin_lnscale - LOG(factor)
    CALL speak('LIN_GROWTH: rescaled by '//str(factor),2)
  END SUBROUTINE lin_rescale

  SUBROUTINE lin_print
    ! Report the spectrum on stdout. Rank 0 only; every rank holds the same
    ! values by the time this is called.
    USE basic,    ONLY: time
    USE grid,     ONLY: total_nky, kyarray_full
    USE parallel, ONLY: my_id
    IMPLICIT NONE
    INTEGER :: ikyg, idom

    IF (.NOT. LIN_GROWTH) RETURN
    IF (my_id .NE. 0)     RETURN
    IF (.NOT. initialized) RETURN
    IF (.NOT. ANY(measured)) RETURN

    ! Dominant mode: the fastest growing of those actually measured.
    idom = 0
    DO ikyg = 1,total_nky
      IF (.NOT. measured(ikyg)) CYCLE
      IF (idom .EQ. 0) THEN
        idom = ikyg
      ELSE IF (gamma_ky(ikyg) .GT. gamma_ky(idom)) THEN
        idom = ikyg
      ENDIF
    ENDDO
    IF (idom .EQ. 0) RETURN

    WRITE(*,'(A,F8.2,A,F6.3,A,SP,F8.4,A,SP,F8.4,A,SS,ES9.2)') &
      '|t = ', time, '| ky = ', kyarray_full(idom), &
      ' | g = ', gamma_ky(idom), ' | w = ', omega_ky(idom), &
      ' | dr = ', drift(idom)

    ! With a single measured ky the line above already says everything.
    IF (COUNT(measured) .GT. 1) THEN
      DO ikyg = 1,total_nky
        IF (.NOT. measured(ikyg)) CYCLE
        WRITE(*,'(A,F6.3,A,SP,F8.4,A,SP,F8.4,A,SS,ES9.2)') &
          '    ky=', kyarray_full(ikyg), '  g=', gamma_ky(ikyg), &
          '  w=', omega_ky(ikyg), '  dr=', drift(ikyg)
      ENDDO
    ENDIF
  END SUBROUTINE lin_print

END MODULE linear_growth
