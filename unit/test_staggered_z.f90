PROGRAM test_staggered_z
  ! Pins the direction in which processing_mod carries the moments from one
  ! staggered z grid to the other.
  !
  ! With SG = .true. the equation of N^p reads its p+/-1 neighbours on the other
  ! grid, through interp_napj (mirror and trapping terms) and ddz_napj (Landau
  ! damping). compute_interp_z and compute_gradients_z pick the stencil of every
  ! p from its parity, eo = MODULO(p,2)+1, which in calculus_mod means
  !
  !   p even : eo = 1, odd to even rule, result half a cell below the input point
  !   p odd  : eo = 2, even to odd rule, result half a cell above the input point
  !
  ! where "the input point" is the sample of the same moment at the same
  ! interior index iz, i.e. storage index iz + ngz/2. The two routines must
  ! agree point for point, otherwise the mirror force and the Landau damping of
  ! one equation would sample the same neighbour at two different z.
  !
  ! Every moment is filled with a cubic in z, on which both 4th order midpoint
  ! rules are exact, so the only thing that can move the answer is which way
  ! the half cell goes. A moment that is not interpolated at all (eo left
  ! undefined falls through to the identity branch of interp_z) is off by the
  ! cubic's variation over half a cell and fails.
  !
  ! Single rank, z-pinch, one species with adiabatic electrons.
  USE prec_const, ONLY: xp
  USE basic,      ONLY: lu_in
  USE parallel,   ONLY: ppinit
  USE grid,       ONLY: local_na, local_np, local_nj, local_nky, local_nkx, local_nz, &
                        ngp, ngj, ngz, ngz_o2, total_nz, nzgrid, parray, zarray, ieven, deltaz
  USE fields,     ONLY: moments
  USE array,      ONLY: nadiab_moments, interp_napj, ddz_napj
  USE time_integration, ONLY: set_updatetlevel
  USE processing, ONLY: compute_interp_z, compute_gradients_z
  USE test_harness
  IMPLICIT NONE

  INTEGER,  PARAMETER :: T_PMAX = 4, T_JMAX = 1
  INTEGER,  PARAMETER :: T_NX = 4, T_NY = 4, T_NZ = 16
  REAL(xp), PARAMETER :: T_LX = 40._xp, T_LY = 40._xp
  ! Reference cubic, real and imaginary part
  REAL(xp), PARAMETER :: ar(0:3) = (/ 0.30_xp,-0.80_xp, 0.45_xp,-0.25_xp /)
  REAL(xp), PARAMETER :: ai(0:3) = (/-0.15_xp, 0.60_xp,-0.35_xp, 0.18_xp /)

  COMPLEX(xp), ALLOCATABLE :: got_i(:), want_i(:), got_d(:), want_d(:)
  REAL(xp) :: tol, zout
  INTEGER  :: ia, ip, ij, iky, ikx, iz, k, n, parity
  CHARACTER(len=4) :: pname

  CALL ppinit
  CALL start_suite('staggered_z')
  tol = tol_of(1.0e4_xp)

  CALL write_input_deck
  CALL readinputs
  CALL auxval
  CALL set_updatetlevel(1)

  !--------------------------------------------------------------- Setup
  CALL section('the run was configured the way the test assumes')
  CALL check_int('two staggered z grids', nzgrid, 2)
  CALL check_int('four z ghost cells', ngz, 4)
  CALL check_int('single rank: local_nz is the total', local_nz, total_nz)
  CALL check('enough z points for the 3D branch of compute_gradients_z', total_nz .GT. 4)
  CALL check('p ghosts present on both sides', ngp .GE. 2)

  !--------------------------------------------------------------- Fill
  ! Every (a,p,j,ky,kx) gets the cubic times its own weight, so a result
  ! written to the wrong slot cannot pass. The cubic is sampled on the storage
  ! lattice of grid 1; both grids share its spacing and index layout, and the
  ! stencils only see indices, so which grid a given p physically sits on does
  ! not enter.
  moments = 0._xp
  DO iz = 1,local_nz+ngz
    DO ikx = 1,local_nkx
    DO iky = 1,local_nky
    DO ij = 1,local_nj+ngj
    DO ip = 1,local_np+ngp
    DO ia = 1,local_na
      nadiab_moments(ia,ip,ij,iky,ikx,iz) = weight(ia,ip,ij,iky,ikx)*poly(zarray(iz,ieven))
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO
  ENDDO

  CALL compute_interp_z
  CALL compute_gradients_z

  !--------------------------------------------------------------- Check
  ! Both routines also run over the p and j ghosts, which the RHS reads as
  ! p+/-1 and j+/-1, so those are checked too.
  n = local_na*(local_nj+ngj)*local_nky*local_nkx*local_nz
  ALLOCATE(got_i(n), want_i(n), got_d(n), want_d(n))
  DO parity = 0,1
    IF (parity .EQ. 0) THEN
      CALL section('p even: carried half a cell down (odd to even rule)')
    ELSE
      CALL section('p odd: carried half a cell up (even to odd rule)')
    ENDIF
    DO ip = 1,local_np+ngp
      IF (MODULO(parray(ip),2) .NE. parity) CYCLE
      k = 0
      DO iz = 1,local_nz
        zout = zarray(iz+ngz_o2,ieven) + REAL(2*parity-1,xp)*0.5_xp*deltaz
        DO ikx = 1,local_nkx
        DO iky = 1,local_nky
        DO ij = 1,local_nj+ngj
        DO ia = 1,local_na
          k = k + 1
          got_i(k)  = interp_napj(ia,ip,ij,iky,ikx,iz)
          want_i(k) = weight(ia,ip,ij,iky,ikx)*poly(zout)
          got_d(k)  = ddz_napj(ia,ip,ij,iky,ikx,iz)
          want_d(k) = weight(ia,ip,ij,iky,ikx)*dpoly(zout)
        ENDDO
        ENDDO
        ENDDO
        ENDDO
      ENDDO
      WRITE(pname,'(I0)') parray(ip)
      CALL check_cmpx_arr('interp_napj, p = '//TRIM(pname), got_i, want_i, tol)
      CALL check_cmpx_arr('ddz_napj,    p = '//TRIM(pname), got_d, want_d, tol)
    ENDDO
  ENDDO

  DEALLOCATE(got_i, want_i, got_d, want_d)
  CALL finish_suite

CONTAINS

  COMPLEX(xp) FUNCTION poly(x)
    REAL(xp), INTENT(IN) :: x
    INTEGER :: k_
    poly = (0._xp,0._xp)
    DO k_ = 0,3
      poly = poly + CMPLX(ar(k_), ai(k_), xp)*x**k_
    ENDDO
  END FUNCTION poly

  COMPLEX(xp) FUNCTION dpoly(x)
    REAL(xp), INTENT(IN) :: x
    INTEGER :: k_
    dpoly = (0._xp,0._xp)
    DO k_ = 1,3
      dpoly = dpoly + REAL(k_,xp)*CMPLX(ar(k_), ai(k_), xp)*x**(k_-1)
    ENDDO
  END FUNCTION dpoly

  !> A distinct O(1) factor per (a,p,j,ky,kx) slot.
  COMPLEX(xp) FUNCTION weight(ia_, ip_, ij_, iky_, ikx_)
    INTEGER, INTENT(IN) :: ia_, ip_, ij_, iky_, ikx_
    weight = CMPLX(1._xp + 0.1_xp*REAL(ip_,xp) + 0.01_xp*REAL(iky_,xp), &
                   0.2_xp*REAL(ij_,xp) - 0.03_xp*REAL(ikx_,xp) + 0.5_xp*REAL(ia_-1,xp), xp)
  END FUNCTION weight

  !> Feed the namelists through a scratch file so nothing is left on disk.
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
    WRITE(lu_in,'(A)')       '  SG   = .true.'
    WRITE(lu_in,'(A)')       '  Nexc = 1'
    WRITE(lu_in,'(A)')       '/'
    WRITE(lu_in,'(A)')       '&GEOMETRY'
    WRITE(lu_in,'(A)')       "  geom = 'z-pinch'"
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

END PROGRAM test_staggered_z
