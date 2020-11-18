SUBROUTINE PPSETUP
	  !
	  ! sets up the MPI topology and distribute the work betwen the nprocs
	  !
	  use prec_const
	  use basic
	  USE basic_mpi
	  use model
	  !
	  IMPLICIT NONE
	  !
	  INTEGER :: njobsie,njobsei,njobsaa,i
	  INTEGER  :: iproc =0
	  INTEGER, DIMENSION(nprocs) :: procslabels
	  LOGICAL, DIMENSION(:),allocatable :: logical_mask
	  INTEGER, DIMENSION(:), allocatable :: buff_
	  !
	  ! loc. vars.
	  INTEGER :: lbaremax,lbarimax
	  INTEGER, DIMENSION(2) :: pjs,pje
	  INTEGER :: ii,grp_rank
	  !
	  INTEGER :: dims(2)
	  LOGICAL :: is_periodic(2)
	  INTEGER :: Ndim_energy
	  INTEGER :: Nprocs_lbar
	  INTEGER :: shift_
	  LOGICAL :: reorder
	  !
	  IF(nprocs > 1 ) THEN
	     !
	     ! Generate an 1D topology: the collisional rows are shared with procesors
	     !
	     lbaremax = numbe(Pmaxe,Jmaxe)
	     lbarimax = numbi(Pmaxi,Jmaxi)
	     !
	     ! Create MPI group E/I/SELF
	     Ngroups = 0
	     IF(eicolls) Ngroups = Ngroups +1
	     IF(iecolls) Ngroups = Ngroups +1
	     IF(eecolls .or. iicolls) Ngroups = Ngroups +1
	     !
	     ! Assign group ranks
	     !
	     iproc = 0
	     IF(eicolls) THEN
	        ALLOCATE(ranks_e(1:Cenproc))
	        DO ii=1,Cenproc
	           ranks_e(ii) = ii-1
	        ENDDO
	        iproc = Cenproc
	        !!  if(me .eq. 0) write(*,*) ranks_e(:)
	     ENDIF
	     !
	     IF(iecolls) THEN
	        ALLOCATE(ranks_i(1:Cinproc))
	        DO ii=1,Cinproc
	           ranks_i(ii) = iproc + ii -1
	        ENDDO
	        iproc = iproc + Cinproc
	        !!    if(me .eq. 0) write(*,*) ranks_i(:)
	     ENDIF
	     !
	     IF(eecolls .or. iicolls) THEN
	        ALLOCATE(ranks_self(1:Caanproc))
	        DO ii=1,Caanproc
	           ranks_self(ii) = iproc + ii -1
	        ENDDO
	        !!      if(me .eq. 0) write(*,*) ranks_self(:)
	     ENDIF
	     !
	     ! CREATE Groups
	     CALL MPI_COMM_GROUP(MPI_COMM_WORLD,grp_world,ierr)
	     !
	     IF(eicolls) THEN
	        CALL MPI_GROUP_INCL(grp_world,Cenproc,ranks_e,grp_e,ierr)
	        ! Get my rank in MPI group
	        CALL MPI_GROUP_RANK(grp_e,me_e,ierr)
	        CALL MPI_COMM_CREATE_GROUP(MPI_COMM_WORLD,grp_e,0,comm_e,ierr)
	        IF(me_e .ne. MPI_UNDEFINED ) is_e = .true.
	     ENDIF
	     !
	     IF(iecolls) THEN
	        CALL MPI_GROUP_INCL(grp_world,Cinproc,ranks_i,grp_i,ierr)
	        ! Get my rank in MPI group
	        CALL MPI_GROUP_RANK(grp_i,me_i,ierr)
	        CALL MPI_COMM_CREATE_GROUP(MPI_COMM_WORLD,grp_i,0,comm_i,ierr)
	        IF(me_i.ne. MPI_UNDEFINED) is_i = .true.
	     ENDIF
	     !
	     IF(eecolls .or. iicolls) THEN
	        ! Create group
	        CALL MPI_GROUP_INCL(grp_world,Caanproc,ranks_self,grp_self,ierr)
	        ! Get my rank in MPI group
	        CALL MPI_GROUP_RANK(grp_self,me_self,ierr)
	        CALL MPI_COMM_CREATE_GROUP(MPI_COMM_WORLD,grp_self,0,comm_self,ierr)
	        IF(me_self .ne. MPI_UNDEFINED) is_self = .true.
	        !
	        IF( IFCOULOMB .and. IFGK) THEN
	           ! Create Cartesian Topology
	           ! Naanprocs_j_ii sould be a mutliple of Cannproc
	           Nprocs_rows_ii = floor(real(Caanproc)/Nprocs_j_ii)
	           dims = (/Nprocs_rows_ii, Nprocs_j_ii/) !
	           is_periodic = (/.FALSE.,.FALSE./)
	           reorder = .true.
	           !
	           CALL MPI_CART_CREATE(comm_self,2,dims,is_periodic,reorder,comm_cart_self,ierr)
	           CALL MPI_CART_COORDS(comm_cart_self,me_self,2,mycoords_ii,ierr)
	           !
	           !! debug
	           !! print*,Nprocs_rows_ii, dims
	        ENDIF
	        !
	     ENDIF
	     !
	     ! local 1D indices: start form 1
	     ! electrons
	     IF(eicolls .and. is_e ) THEN
	        IF(MPI_chksz .eq. 0 ) THEN
	           lbare_s_l = me_e*floor(real(lbaremax/Cenproc)) + 1
	           IF( me_e .eq. (Cenproc - 1) ) THEN
	              lbare_e_l = numbe(Pmaxe,Jmaxe)
	           ELSE
	              lbare_e_l = lbare_s_l + floor(real(lbaremax/Cenproc))-1
	           ENDIF
	        ELSE
	           ! fill up from the bottom with chunk size =  MPI_chksz
	           shift_ = numbe(Pmaxe,Jmaxe) - MPI_chksz*(Cenproc -1)
	           IF( me_e .eq. 0 ) THEN
	              lbare_s_l = 1
	              lbare_e_l = shift_
	           ELSE
	              lbare_s_l = shift_ + (me_e -1)*MPI_chksz +1
	              lbare_e_l = lbare_s_l + MPI_chksz  -1
	           ENDIF
	           ! Debug
	           !! write(*,*) me_e, lbare_s_l,lbare_e_l,numbe(Pmaxe,Jmaxe)
	        ENDIF
	     ENDIF

	     ! ions
	     IF(iecolls .and. is_i ) THEN
	        IF(MPI_chksz .eq. 0 ) THEN
	           lbari_s_l = me_i*floor(real(lbarimax/Cinproc)) + 1
	           IF( me_i .eq. (Cinproc - 1) ) THEN
	              lbari_e_l = numbi(Pmaxi,Jmaxi)
	           ELSE
	              lbari_e_l = lbari_s_l + floor(real(lbarimax/Cinproc))-1
	           ENDIF
	        ELSE
	           ! fill up from the bottom with chunk size =  MPI_chksz
	           shift_ = numbi(Pmaxi,Jmaxi) - MPI_chksz*(Cinproc -1)
	           IF( me_i .eq. 0 ) THEN
	              lbari_s_l = 1
	              lbari_e_l = shift_
	           ELSE
	              lbari_s_l = shift_ + (me_i -1)*MPI_chksz +1
	              lbari_e_l = lbari_s_l + MPI_chksz  -1
	           ENDIF
	           ! Debug
	           !! write(*,*) me_e, lbare_s_l,lbare_e_l,numbe(Pmaxe,Jmaxe)
	        ENDIF
	     ENDIF

	     ! self electron collisions
	     IF( eecolls .and. is_self ) THEN
	        IF(MPI_chksz .eq. 0) THEN
	           lbaree_s_l = me_self*floor(real(lbaremax/Caanproc)) + 1
	           IF( me_self .eq. (Caanproc - 1) ) THEN
	              lbaree_e_l = numbe(Pmaxe,Jmaxe)
	           ELSE
	              lbaree_e_l = lbaree_s_l + floor(real(lbaremax/Caanproc))-1
	           ENDIF
	        ELSE
	           ! fill up from the bottom with chunk size =  MPI_chksz
	           shift_ = numbe(Pmaxe,Jmaxe) - MPI_chksz*(Caanproc -1)
	           IF( me_self .eq. 0 ) THEN
	              lbaree_s_l = 1
	              lbaree_e_l = shift_
	           ELSE
	              lbaree_s_l = shift_ + (me_self -1)*MPI_chksz +1
	              lbaree_e_l = lbaree_s_l + MPI_chksz  -1
	           ENDIF
	           ! Debug
	           !! write(*,*) me_e, lbare_s_l,lbare_e_l,numbe(Pmaxe,Jmaxe)
	        ENDIF
	     ENDIF
	     !
	     ! self ion collisions
	      IF( iicolls .and. is_self ) THEN
	         IF(MPI_chksz .eq. 0 ) THEN
	            IF(IFCOULOMB .and. IFGK) THEN
	               ! MPI implementation for GK Coulomb (start at 1)
	              lbarii_s_l = mycoords_ii(1)*floor(real(lbarimax/Nprocs_rows_ii)) +1
	              IF(mycoords_ii(1) .eq. Nprocs_rows_ii-1) THEN
	                 lbarii_e_l = numbi(Pmaxi,Jmaxi)
	              ELSE
	                 lbarii_e_l = lbarii_s_l + floor(real(lbarimax/Nprocs_rows_ii)) -1
	              ENDIF
	              !
	              ! local energy component (start at 0)
	              jimaxx_s_l = mycoords_ii(2)*floor(real(JEmaxx/Nprocs_j_ii))
	              IF(mycoords_ii(2) .eq. Nprocs_j_ii-1) THEN
	                 jimaxx_e_l = JEmaxx
	              ELSE
	                 jimaxx_e_l = jimaxx_s_l  + floor(real(JEmaxx/Nprocs_j_ii)) -1
	              ENDIF

	              !! Debug
	              !! print*,mycoords_ii(:), lbarii_s_l,lbarii_e_l,jimaxx_s_l,jimaxx_e_l

	              !
	              ! Local FLR indices: (not parallelized in the FLR dimensions)
	              nimaxxFLR_s_l = 0
	              nimaxxFLR_e_l = nimaxxFLR
	              npimaxxFLR_s_l = 0
	              npimaxxFLR_e_l = npimaxxFLR

	            ELSE
	              lbarii_s_l = me_self*floor(real(lbarimax/Caanproc)) + 1
	              IF( me_self .eq. (Caanproc - 1) ) THEN
	                 lbarii_e_l = numbi(Pmaxi,Jmaxi)
	              ELSE
	                 lbarii_e_l = lbarii_s_l + floor(real(lbarimax/Caanproc))-1
	              ENDIF
	              ! Local FLR indices
	              nimaxxFLR_s_l = 0
	              nimaxxFLR_e_l = nimaxxFLR
	              npimaxxFLR_s_l = 0
	              npimaxxFLR_e_l = npimaxxFLR
	              jimaxx_s_l = 0
	              jimaxx_e_l = JEmaxx
	            ENDIF
	         ELSE
	           ! fill up from the bottom with chunk size =  MPI_chksz
	           ! MPI implementation for GK Coulomb (start at 1)
	           !
	            IF(IFCOULOMB .and. IFGK) THEN
	              shift_ = numbi(Pmaxi,Jmaxi) - MPI_chksz*(Nprocs_rows_ii -1)
	              IF( mycoords_ii(1) .eq. 0 ) THEN
	                 lbarii_s_l = 1
	                 lbarii_e_l = shift_
	              ELSE
	                 lbarii_s_l = shift_ + (mycoords_ii(1) -1)*MPI_chksz +1
	                 lbarii_e_l = lbarii_s_l + MPI_chksz  -1
	              ENDIF
	              !
	              ! local energy component (start at 0)
	              jimaxx_s_l = mycoords_ii(2)*floor(real(JEmaxx/Nprocs_j_ii))
	              IF(mycoords_ii(2) .eq. Nprocs_j_ii-1) THEN
	                 jimaxx_e_l = JEmaxx
	              ELSE
	                 jimaxx_e_l = jimaxx_s_l  + floor(real(JEmaxx/Nprocs_j_ii)) -1
	              ENDIF
	              ! debug
	!!              print*, lbarii_s_l,lbarii_e_l,jimaxx_s_l,jimaxx_e_l
	              !
	              ! Local FLR indices
	              nimaxxFLR_s_l = 0
	              nimaxxFLR_e_l = nimaxxFLR
	              npimaxxFLR_s_l = 0
	              npimaxxFLR_e_l = npimaxxFLR

	            ELSE
	              shift_ = numbi(Pmaxi,Jmaxi) - MPI_chksz*(Caanproc -1)
	              IF( me_self .eq. 0 ) THEN
	                 lbarii_s_l = 1
	                 lbarii_e_l = shift_
	              ELSE
	                 lbarii_s_l = shift_ + (me_self -1)*MPI_chksz +1
	                 lbarii_e_l = lbarii_s_l + MPI_chksz  -1
	              ENDIF
	            ENDIF
	         ENDIF
	      ENDIF
	      !
	     ! Get the local Hermite-Laguerre indices
	     !
	     !                 ELECTRON ELECTRON
	     IF(eicolls .and. (is_e .or. me  .eq. 0)) THEN
	        !
	        pjs = invnumbe(lbare_s_l)
	        pje = invnumbe(lbare_e_l)
	        !
	        Pmaxe_s_l = pjs(1)
	        Jmaxe_s_l = pjs(2)
	        Pmaxe_e_l = pje(1)
	        Jmaxe_e_l = pje(2)

	        CALL ALLOCATE_ARRAY(Je_s_l,Pmaxe_s_l,Pmaxe_e_l)
	        CALL ALLOCATE_ARRAY(Je_e_l,Pmaxe_s_l,Pmaxe_e_l)

	        ! local Laguere indices
	        Je_s_l(Pmaxe_s_l) = Jmaxe_s_l
	        Je_e_l(:) = Jmaxe
	        Je_e_l(Pmaxe_e_l) = Jmaxe_e_l
	        !
	     ENDIF

	     !        ION - ELECTRON
	     !
	     IF(iecolls .and. is_i  ) THEN
	        pjs = invnumbi(lbari_s_l)
	        Pmaxi_s_l = pjs(1)
	        Jmaxi_s_l = pjs(2)
	        pje = invnumbi(lbari_e_l)
	        Pmaxi_e_l = pje(1)
	        Jmaxi_e_l = pje(2)

	        CALL ALLOCATE_ARRAY(Ji_s_l,Pmaxi_s_l,Pmaxi_e_l)
	        CALL ALLOCATE_ARRAY(Ji_e_l,Pmaxi_s_l,Pmaxi_e_l)

	        ! local Laguere indices
	        Ji_s_l(Pmaxi_s_l) = Jmaxi_s_l
	        Ji_e_l(:) = Jmaxi
	        Ji_e_l(Pmaxi_e_l) = Jmaxi_e_l
	        !
	     ENDIF

	     !    Electron-Electron
	     !
	     IF(eecolls .and. is_self) THEN

	        pjs = invnumbe(lbaree_s_l)
	        Pmaxee_s_l = pjs(1)
	        Jmaxee_s_l = pjs(2)

	        pje = invnumbe(lbaree_e_l)
	        Pmaxee_e_l = pje(1)
	        Jmaxee_e_l = pje(2)
	        !
	        CALL ALLOCATE_ARRAY(Jee_s_l,Pmaxee_s_l,Pmaxee_e_l)
	        CALL ALLOCATE_ARRAY(Jee_e_l,Pmaxee_s_l,Pmaxee_e_l)
	        !
	     ENDIF

	     !
	     !             ION - ION
	     IF(iicolls .and. is_self) THEN

	        pjs = invnumbi(lbarii_s_l)
	        Pmaxii_s_l = pjs(1)
	        Jmaxii_s_l = pjs(2)

	        pje = invnumbi(lbarii_e_l)
	        Pmaxii_e_l = pje(1)
	        Jmaxii_e_l = pje(2)
	        !
	        CALL ALLOCATE_ARRAY(Jii_s_l,Pmaxii_s_l,Pmaxii_e_l)
	        CALL ALLOCATE_ARRAY(Jii_e_l,Pmaxii_s_l,Pmaxii_e_l)

	        ! local Laguere indices
	        Jii_s_l(Pmaxii_s_l) = Jmaxii_s_l
	        Jii_e_l(:) = Jmaxi
	        Jii_e_l(Pmaxii_e_l) = Jmaxii_e_l
	        !
	     ENDIF
	     !
	  ELSE   ! ... if NPROC = 1

	     ! set local indices to global
	     ! electrons
	     Pmaxe_s_l =0
	     Pmaxe_e_l = Pmaxe
	     Jmaxe_s_l =0
	     Jmaxe_e_l  = Jmaxe

	     CALL ALLOCATE_ARRAY(Je_s_l,Pmaxe_s_l,Pmaxe_e_l)
	     CALL ALLOCATE_ARRAY(Je_e_l,Pmaxe_s_l,Pmaxe_e_l)
	     Je_s_l(:) = 0
	     Je_e_l(:) = Jmaxe

	     ! ion
	     Pmaxi_s_l =0
	     Pmaxi_e_l = Pmaxi
	     Jmaxi_s_l =0
	     Jmaxi_e_l  = Jmaxi
	     !
	     CALL ALLOCATE_ARRAY(Ji_s_l,Pmaxi_s_l,Pmaxi_e_l)
	     CALL ALLOCATE_ARRAY(Ji_e_l,Pmaxi_s_l,Pmaxi_e_l)
	     Ji_s_l(:) = 0
	     Ji_e_l(:) = Jmaxi
	     !
	     ! self electrons
	     !
	     Pmaxee_s_l = 0
	     Pmaxee_e_l = Pmaxe
	     Jmaxee_s_l = 0
	     Jmaxee_e_l = Jmaxe
	     !
	     CALL ALLOCATE_ARRAY(Jee_s_l,Pmaxee_s_l,Pmaxee_e_l)
	     CALL ALLOCATE_ARRAY(Jee_e_l,Pmaxee_s_l,Pmaxee_e_l)
	     Jee_s_l(:) = 0
	     Jee_e_l(:) = Jmaxe
	     !
	     ! self ions
	     !
	     Pmaxii_s_l = 0
	     Pmaxii_e_l = Pmaxi
	     Jmaxii_s_l = 0
	     Jmaxii_e_l = Jmaxi
	     !
	     CALL ALLOCATE_ARRAY(Jii_s_l,Pmaxii_s_l,Pmaxii_e_l)
	     CALL ALLOCATE_ARRAY(Jii_e_l,Pmaxii_s_l,Pmaxii_e_l)
	     Jii_s_l(:) = 0
	     Jii_e_l(:) = Jmaxi
	     !
	     ! Set local GK indices to global
	     ! Engery component
	     jimaxx_s_l = 0
	     jimaxx_e_l = JEmaxx
	     !
	     ! FLR indices
	     nimaxxFLR_s_l = 0
	     nimaxxFLR_e_l = nimaxxFLR
	     !
	     npimaxxFLR_s_l = 0
	     npimaxxFLR_e_l = nimaxxFLR
	     !
	     ! Set all com to MPI_WORLD_COMM
	     comm_e = MPI_COMM_WORLD
	     comm_i = MPI_COMM_WORLD
	     comm_self = MPI_COMM_WORLD
	     comm_cart_self = MPI_COMM_WORLD

	     !
	     me_e = 0
	     me_i= 0
	     me_self =0
	     !
	     is_e = .true.
	     is_i = .true.
	     is_self = .true.
	     !
	  ENDIF
	  !
	END SUBROUTINE PPSETUP
