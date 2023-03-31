module cosolver_interface
! contains the Hermite-Laguerre collision operators solved using COSOlver.
USE prec_const, ONLY: xp
IMPLICIT NONE
PRIVATE
PUBLIC :: load_COSOlver_mat, compute_cosolver_coll

CONTAINS
  !******************************************************************************!
  !! compute the collision terms in a (Np x Nj x Nkx x Nky) matrix all at once
  !******************************************************************************!
  SUBROUTINE compute_cosolver_coll(GK_CO)
    USE parallel,    ONLY: num_procs_p, comm_p,dsp_p,rcv_p
    USE grid,        ONLY: &
      local_na, &
      local_np, ngp, total_np, total_nj, ngj,&
      local_nkx, local_nky, local_nz, bar
    USE array,       ONLY: Capj
    USE MPI
    USE closure,     ONLY: evolve_mom
    IMPLICIT NONE
    LOGICAL, INTENT(IN) :: GK_CO
    COMPLEX(xp), DIMENSION(total_np)    :: local_coll, buffer
    COMPLEX(xp), DIMENSION(local_np)    :: TColl_distr
    COMPLEX(xp) :: Tmp_
    INTEGER :: iz,ikx,iky,ij,ip,ia,ikx_C,iky_C,iz_C
    INTEGER :: ierr
    z:DO iz = 1,local_nz
      x:DO ikx = 1,local_nkx
        y:DO iky = 1,local_nky
          a:DO ia = 1,local_na
            j:DO ij = 1,total_nj
              p:DO ip = 1,total_np
              IF(evolve_mom(ip+ngp/2,ij+ngj/2)) THEN !compute for every moments except for closure 1
                  !! Take GK or DK limit
                  IF (GK_CO) THEN ! GK operator (k-dependant)
                    ikx_C = ikx; iky_C = iky; iz_C = iz;
                  ELSE ! DK operator (only one mat for every k)
                    ikx_C = 1;   iky_C = 1; iz_C = 1;
                  ENDIF
                  !! Apply the cosolver collision matrix
                  CALL apply_cosolver_mat(ia,ip,ij,iky,ikx,iz,ikx_C,iky_C,iz_C,Tmp_)
                  local_coll(ip) = Tmp_
                ELSE
                  local_coll(ip) = 0._xp
                ENDIF
              ENDDO p
              IF (num_procs_p .GT. 1) THEN
                ! Reduce the local_sums to root = 0
                CALL MPI_REDUCE(local_coll, buffer, total_np, MPI_DOUBLE_COMPLEX, MPI_SUM, 0, comm_p, ierr)
                ! buffer contains the entire collision term along p, we scatter it between
                ! the other processes (use of scatterv since Pmax/Np is not an integer)
                CALL MPI_SCATTERV(buffer, rcv_p, dsp_p, MPI_DOUBLE_COMPLEX,&
                                  TColl_distr, local_np, MPI_DOUBLE_COMPLEX, &
                                  0, comm_p, ierr)
              ELSE
                TColl_distr = local_coll
              ENDIF
              ! Write in output variable
              DO ip = 1,local_np
                Capj(ia,ip,ij,iky,ikx,iz) = TColl_distr(ip)
              ENDDO
            ENDDO j
          ENDDO a
        ENDDO y
      ENDDO x
    ENDDO z
  END SUBROUTINE compute_cosolver_coll

    !******************************************************************************!
  !! compute the collision terms in a (Np x Nj x Nkx x Nky) matrix all at once
  !******************************************************************************!
  SUBROUTINE apply_cosolver_mat(ia,ip,ij,iky,ikx,iz,ikx_C,iky_C,iz_C,local_coll)
    USE grid,        ONLY: &
      total_na, &
      local_np, parray,parray_full, ngp,&
      total_nj, jarray,jarray_full, ngj, bar, ngz
    USE array,       ONLY: nuCself, Cab_F, nadiab_moments
    USE species,     ONLY: nu_ab
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ia,ip,ij,iky,ikx,iz,ikx_C,iky_C,iz_C
    COMPLEX(xp), INTENT(OUT) :: local_coll
    INTEGER :: ib,iq,il
    INTEGER :: p_int,q_int,j_int,l_int
    INTEGER :: izi, iqi, ili
    izi = iz+ngz/2
    p_int = parray_full(ip)
    j_int = jarray_full(ij)
    !! Apply the cosolver collision matrix
    local_coll = 0._xp ! Initialization
    q:DO iq = 1,local_np
      iqi   = iq + ngp/2
      q_int = parray(iqi)
      l:DO il = 1,total_nj
        ili   = il + ngj/2
        l_int = jarray(ili)
        ! self interaction + test interaction
        local_coll = local_coll + nadiab_moments(ia,iqi,ili,iky,ikx,izi) &
              * nuCself(ia,bar(p_int,j_int), bar(q_int,l_int), iky_C, ikx_C, iz_C)
        ! sum the contribution over the other species
        b:DO ib = 1,total_na
          IF(ib .NE. ia) THEN
            ! Field contribution
            local_coll = local_coll + nadiab_moments(ib,iqi,ili,iky,ikx,izi) &
                * nu_ab(ia,ib) * Cab_F(ia,ib,bar(p_int,j_int), bar(q_int,l_int), iky_C, ikx_C, iz_C)
          ENDIF
        ENDDO b
      ENDDO l
    ENDDO q
  END SUBROUTINE apply_cosolver_mat

    !******************************************************************************!
    !!!!!!! Load the collision matrix coefficient table from COSOlver results
    !******************************************************************************!
    SUBROUTINE load_COSOlver_mat(GK_CO,INTERSPECIES,matfile,collision_kcut) ! Load a sub matrix from iCa files (works for pmax,jmax<=P_full,J_full)
      USE basic,       ONLY: allocate_array, speak
      USE parallel,    ONLY: comm_p, my_id
      USE grid,        ONLY: &
        local_na, total_na, &
        local_nkx,local_nky, kparray,&
        local_nz, ngz, bar,&
        pmax, jmax, ieven
      USE array,       ONLY: Caa, Cab_T, Cab_F, nuCself
      USE MPI
      USE model,       ONLY: Na, LINEARITY
      USE species,     ONLY: name, nu_ab
      USE futils
      IMPLICIT NONE
      ! Input
      LOGICAL,            INTENT(IN) :: GK_CO, INTERSPECIES
      CHARACTER(len=128), INTENT(IN) :: matfile    ! COSOlver matrix file names
      REAL(xp),           INTENT(IN) :: collision_kcut
      ! Local variables
      REAL(xp), DIMENSION(:,:),       ALLOCATABLE :: Caa_full,CabT_full, CabF_full  ! To load the self entire matrices
      REAL(xp), DIMENSION(:,:,:,:),   ALLOCATABLE :: Caa__kp         ! To store the coeff that will be used along kperp
      REAL(xp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: CabF_kp,CabT_kp ! ''
      REAL(xp), DIMENSION(:),         ALLOCATABLE :: kp_grid_mat     ! kperp grid of the matrices
      REAL(xp), DIMENSION(2) :: dims
      ! Indices for row and columns of the COSOlver matrix (4D compressed 2D matrices)
      INTEGER  :: irow_sub, irow_full, icol_sub, icol_full
      INTEGER  :: fid                                       ! file indexation
      INTEGER  :: ip, ij, il, ikx, iky, iz, ik, ikp, ikps_C,ikpe_C,ia,ib  ! indices for loops
      INTEGER  :: pdim, jdim                                ! dimensions of the COSOlver matrices
      INTEGER  :: ikp_next, ikp_prev, nkp_mat, ikp_mat
      REAL(xp) :: kp_max,kperp_sim, kperp_mat, zerotoone
      CHARACTER(len=128) :: var_name, ikp_string, name_a, name_b
      CHARACTER(len=1)   :: letter_a, letter_b
      ! Opening the compiled cosolver matrices results
      CALL openf(matfile,fid, 'r', 'D', mpicomm=comm_p);
      ! Get matrices dimensions (polynomials degrees and kperp grid)
      CALL getarr(fid, '/dims_i', dims) ! Get the ion polynomial degrees (consider same for electrons)
      pdim = dims(1); jdim = dims(2);
      !! Here we stop if the matrix is too small, we could put zero to these coefficients otherwise?
      IF ( ((pdim .LT. pmax) .OR. (jdim .LT. jmax)) .AND. (my_id .EQ. 0)) ERROR STOP '>> ERROR << P,J Matrix too small'
      ! Get the dimension kperp grid of the matrices for GK operator
      CALL getsize(fid, '/coordkperp', nkp_mat)
      CALL allocate_array(kp_grid_mat, 1,nkp_mat)
      CALL getarr(fid, '/coordkperp', kp_grid_mat)
      kp_max = MAXVAL(kparray)
      ! check that we have enough kperps mat, if not we apply the kpmax matrix to all k>kpmax
      IF (LINEARITY .NE. 'linear') THEN
        IF ( (kp_grid_mat(nkp_mat) .LT. 2./3.*kp_max) .AND. (my_id .EQ. 0)) WRITE(*,*) 'warning: Matrix kperp grid too small'
      ELSE
        IF ( (kp_grid_mat(nkp_mat) .LT. kp_max) .AND. (my_id .EQ. 0)) WRITE(*,*) 'warning: Matrix kperp grid too small !!'
      ENDIF
      IF (GK_CO) THEN ! GK operator (k-dependant)
        ikps_C = 1; ikpe_C = nkp_mat
      ELSE ! DK operator, only the k=0 mat applied to all k
        ikps_C = 1; ikpe_C = 1
      ENDIF
      ! allocate the temporary matrices
      CALL allocate_array(  Caa__kp, 1,Na, 1,(pmax+1)*(jmax+1), 1,(pmax+1)*(jmax+1), ikps_C,ikpe_C)
      CALL allocate_array(  CabF_kp, 1,Na, 1,Na, 1,(pmax+1)*(jmax+1), 1,(pmax+1)*(jmax+1), ikps_C,ikpe_C)
      CALL allocate_array(  CabT_kp, 1,Na, 1,Na, 1,(pmax+1)*(jmax+1), 1,(pmax+1)*(jmax+1), ikps_C,ikpe_C)
      CALL allocate_array(  Caa_full,1,(pdim+1)*(jdim+1), 1,(pdim+1)*(jdim+1))
      CALL allocate_array( CabT_full,1,(pdim+1)*(jdim+1), 1,(pdim+1)*(jdim+1))
      CALL allocate_array( CabF_full,1,(pdim+1)*(jdim+1), 1,(pdim+1)*(jdim+1))
      ! Loop over every kperp values that will get the collision operator
      kp:DO ikp = ikps_C,ikpe_C
        ! Kperp value in string format to select in cosolver hdf5 file
        IF (GK_CO) THEN
          write(ikp_string,'(i5.5)') ikp-1
        ELSE
          write(ikp_string,'(i5.5)') 0
        ENDIF
        a:DO ia = 1,local_na
          name_a = name(ia); letter_a = name_a(1:1)
          ! get the self colision matrix
          ! Naming of the array to load (kperp dependant)
          ! we look for the data stored at e.g. '00001/Caapj/Ceepj'
          WRITE(var_name,'(a,a,a,a,a)') TRIM(ADJUSTL(ikp_string)),'/Caapj/C',letter_a,letter_a,'pj'
          CALL getarr(fid, var_name, Caa_full) ! get array
          ! Fill sub array with the usefull polynmial degrees only
          DO ip = 0,pmax ! Loop over rows
            DO ij = 0,jmax
              irow_sub  = (jmax +1)*ip + ij +1
              irow_full = (jdim +1)*ip + ij +1
              DO il = 0,pmax ! Loop over columns
                DO ik = 0,jmax
                  icol_sub  = (jmax +1)*il + ik +1
                  icol_full = (jdim +1)*il + ik +1
                  Caa__kp (ia,irow_sub,icol_sub,ikp) = Caa_full(irow_full,icol_full)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          b: DO ib = 1,local_na
          name_b = name(ib); letter_b = name_b(1:1)
          IF(INTERSPECIES .AND. (ib .NE. ia)) THEN ! Pitch angle is only applied on like-species
              !!!!!!!!!!!!!!! Test and field matrices !!!!!!!!!!!!!!
              ! we look for the data stored at e.g. '00001/Ceipj/CeipjT'
              WRITE(var_name,'(a,a,a,a,a,a,a,a)') TRIM(ADJUSTL(ikp_string)),'/C',letter_a,letter_b,'pj/C',letter_a,letter_b,'pjT'
              CALL getarr(fid, var_name, CabT_full)
              ! we look for the data stored at e.g. '00001/Ceipj/CeipjF'
              WRITE(var_name,'(a,a,a,a,a,a,a,a)') TRIM(ADJUSTL(ikp_string)),'/C',letter_a,letter_b,'pj/C',letter_a,letter_b,'pjF'
              CALL getarr(fid, var_name, CabF_full)
              ! Fill sub array with only usefull polynmials degree
              DO ip = 0,pmax ! Loop over rows
              DO ij = 0,jmax
                    irow_sub  = (jmax +1)*ip + ij +1
                    irow_full = (jdim +1)*ip + ij +1
                    DO il = 0,pmax ! Loop over columns
                    DO ik = 0,jmax
                          icol_sub  = (jmax +1)*il + ik +1
                          icol_full = (jdim +1)*il + ik +1
                          CabT_kp(ia,ib,irow_sub,icol_sub,ikp) = CabT_full(irow_full,icol_full)
                          CabF_kp(ia,ib,irow_sub,icol_sub,ikp) = CabF_full(irow_full,icol_full)
                    ENDDO
                    ENDDO
              ENDDO
              ENDDO
            ELSE
              CabT_kp(ia,ib,:,:,:) = 0._xp; CabF_kp(ia,ib,:,:,:) = 0._xp
            ENDIF
          ENDDO b
        ENDDO a
      ENDDO kp
      DEALLOCATE(Caa_full)
      DEALLOCATE(CabT_full)
      DEALLOCATE(CabF_full)
      CALL closef(fid)

      IF (GK_CO) THEN ! Interpolation of the kperp matrix values on kx ky grid
        IF (my_id .EQ. 0 ) WRITE(*,*) '...Interpolation from matrices kperp to simulation kx,ky...'
        DO ikx = 1,local_nkx
          DO iky = 1,local_nky
            DO iz = 1,local_nz
              ! Check for nonlinear case if we are in the anti aliased domain or the filtered one
              kperp_sim = MIN(kparray(iky,ikx,iz+ngz/2,ieven),collision_kcut) ! current simulation kperp
              ! Find the interval in kp grid mat where kperp_sim is contained
              ! Loop over the whole kp mat grid to find the smallest kperp that is
              ! larger than the current kperp_sim (brute force...)
              DO ikp=1,nkp_mat
                ikp_mat   = ikp ! the first indice of the interval (k0)
                kperp_mat = kp_grid_mat(ikp)
                IF(kperp_mat .GT. kperp_sim) EXIT ! a matrix with kperq > current kx2 + ky2 has been found
              ENDDO
              ! Interpolation
              ! interval boundaries
              ikp_next  = ikp_mat     !index of k1 (larger than kperp_sim thanks to previous loop)
              ikp_prev  = ikp_mat - 1 !index of k0 (smaller neighbour to interpolate inbetween)
              if ( (kp_grid_mat(ikp_prev) .GT. kperp_sim) .OR. (kp_grid_mat(ikp_next) .LT. kperp_sim) ) THEN
                ! write(*,*) 'Warning, linear interp of collision matrix failed!! '
                ! write(*,*) kp_grid_mat(ikp_prev), '<', kperp_sim, '<', kp_grid_mat(ikp_next)
              ENDIF
              ! 0->1 variable for linear interp, i.e. zero2one = (k-k0)/(k1-k0)
              zerotoone = MIN(1._xp,(kperp_sim - kp_grid_mat(ikp_prev))/(kp_grid_mat(ikp_next) - kp_grid_mat(ikp_prev)))
              ! Linear interpolation between previous and next kperp matrix values
              Caa (:,:,:,iky,ikx,iz) = (Caa__kp(:,:,:,ikp_next) - Caa__kp(:,:,:,ikp_prev))*zerotoone + Caa__kp(:,:,:,ikp_prev)
              IF(INTERSPECIES) THEN
                Cab_T(:,:,:,:,iky,ikx,iz) = (CabT_kp(:,:,:,:,ikp_next) - CabT_kp(:,:,:,:,ikp_prev))*zerotoone + CabT_kp(:,:,:,:,ikp_prev)
                Cab_F(:,:,:,:,iky,ikx,iz) = (CabF_kp(:,:,:,:,ikp_next) - CabF_kp(:,:,:,:,ikp_prev))*zerotoone + CabF_kp(:,:,:,:,ikp_prev)
              ELSE
                Cab_T(:,:,:,:,iky,ikx,iz) = 0._xp
                Cab_F(:,:,:,:,iky,ikx,iz) = 0._xp
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ELSE ! DK -> No kperp dep, copy simply to final collision matrices
        Caa (:,:,:,1,1,1)    = Caa__kp(:,:,:,1)
        Cab_T(:,:,:,:,1,1,1) = CabT_kp(:,:,:,:,1)
        Cab_F(:,:,:,:,1,1,1) = CabF_kp(:,:,:,:,1)
      ENDIF
      ! Deallocate auxiliary variables
      DEALLOCATE (Caa__kp); DEALLOCATE (CabT_kp); DEALLOCATE (CabF_kp)

      IF( .NOT. INTERSPECIES ) THEN
        CALL speak("--Like Species operator--")
        Cab_F = 0._xp;
        Cab_T = 0._xp;
      ENDIF
      ! Build the self matrix   
      ! nuCself = nuaa*Caa + sum_b_neq_a nu_ab Cab_T
      DO ia = 1,total_na
        nuCself(ia,:,:,:,:,:) = nu_ab(ia,ia)*Caa(ia,:,:,:,:,:)
        DO ib = 1,total_na
          IF(ib .NE. ia) THEN
            nuCself(ia,:,:,:,:,:) = nuCself(ia,:,:,:,:,:)&
                                    +nu_ab(ia,ib)*Cab_T(ia,ib,:,:,:,:,:)
          ENDIF
        ENDDO
      ENDDO
      CALL speak('============DONE===========')

    END SUBROUTINE load_COSOlver_mat
    !******************************************************************************!
END MODULE cosolver_interface
