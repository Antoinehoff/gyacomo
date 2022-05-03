MODULE restarts
USE basic
USE futils,          ONLY: openf, closef, getarr, getatt, isgroup,&
                           isdataset, getarrnd, putarrnd
USE grid
USE fields
USE diagnostics_par
USE time_integration
USE model, ONLY: KIN_E
IMPLICIT NONE

INTEGER :: rank, sz_, n_
INTEGER :: dims(1) = (/0/)
CHARACTER(LEN=50) :: dset_name
INTEGER :: pmaxe_cp, jmaxe_cp, pmaxi_cp, jmaxi_cp, n0
COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: moments_e_cp
COMPLEX(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: moments_i_cp

PUBLIC :: load_moments!, write_restart

CONTAINS

    !******************************************************************************!
    !!!!!!! Load moments from a previous output file
    !******************************************************************************!
    SUBROUTINE load_moments
        IMPLICIT NONE

        ! Checkpoint filename
        WRITE(rstfile,'(a,a1,i2.2,a3)') TRIM(resfile0),'_',job2load,'.h5'

        IF (my_id .EQ. 0) WRITE(*,'(3x,a)') "Resume from ", rstfile
        ! Open file
        CALL openf(rstfile, fidrst,mpicomm=comm0)
        ! Get the checkpoint moments degrees to allocate memory
        IF (KIN_E) THEN
        CALL getatt(fidrst,"/data/input/" , "pmaxe", pmaxe_cp)
        CALL getatt(fidrst,"/data/input/" , "jmaxe", jmaxe_cp)
        ENDIF
        CALL getatt(fidrst,"/data/input/" , "pmaxi", pmaxi_cp)
        CALL getatt(fidrst,"/data/input/" , "jmaxi", jmaxi_cp)
        IF (my_id .EQ. 0) WRITE(*,*) "Pi_cp = ", pmaxi_cp
        IF (my_id .EQ. 0) WRITE(*,*) "Ji_cp = ", jmaxi_cp
        CALL getatt(fidrst,"/data/input/" , "start_iframe5d", n0)

        IF ((KIN_E .AND. ((pmaxe_cp .NE. pmaxe) .OR. (jmaxe_cp .NE. jmaxe))) .OR.&
         (pmaxi_cp .NE. pmaxi) .OR. (jmaxi_cp .NE. jmaxi)) THEN
         IF(my_id.EQ.0)WRITE(*,*) '! Extending the polynomials basis !'
         CALL load_output_adapt_pj
        ELSE

          ! Find the last results of the checkpoint file by iteration
          n_ = n0+1
          WRITE(dset_name, "(A, '/', i6.6)") "/data/var5d/moments_i", n_ ! start with moments_e/000001
          DO WHILE (isdataset(fidrst, dset_name)) ! If n_ is not a file we stop the loop
          n_ = n_ + 1
          WRITE(dset_name, "(A, '/', i6.6)") "/data/var5d/moments_i", n_ ! updtate file number
          ENDDO
          n_ = n_ - 1 ! n_ is not a file so take the previous one n_-1

          ! Read time dependent attributes to continue simulation
          WRITE(dset_name, "(A, '/', i6.6)") "/data/var5d/moments_i", n_
          CALL getatt(fidrst, dset_name, 'cstep', cstep)
          CALL getatt(fidrst, dset_name, 'time', time)
          CALL getatt(fidrst, dset_name, 'jobnum', jobnum)
          jobnum = jobnum+1
          CALL getatt(fidrst, dset_name, 'iframe2d',iframe2d)
          CALL getatt(fidrst, dset_name, 'iframe5d',iframe5d)
          iframe2d = iframe2d-1; iframe5d = iframe5d-1
          IF(my_id.EQ.0) WRITE(*,*) '.. restart from t = ', time

          ! Read state of system from checkpoint file
          IF (KIN_E) THEN
          WRITE(dset_name, "(A, '/', i6.6)") "/data/var5d/moments_e", n_
          CALL getarrnd(fidrst, dset_name, moments_e(ips_e:ipe_e, ijs_e:ije_e, ikys:ikye, ikxs:ikxe, izs:ize, 1),(/1,3,5/))
          ENDIF
          WRITE(dset_name, "(A, '/', i6.6)") "/data/var5d/moments_i", n_
          CALL getarrnd(fidrst, dset_name, moments_i(ips_i:ipe_i, ijs_i:ije_i, ikys:ikye, ikxs:ikxe, izs:ize, 1),(/1,3,5/))

          CALL closef(fidrst)

          IF (my_id .EQ. 0) WRITE(*,'(3x,a)') "Reading from restart file "//TRIM(rstfile)//" completed!"
        ENDIF

    END SUBROUTINE load_moments
    !******************************************************************************!


    !******************************************************************************!
    !!!!!!! Load moments from a previous output file with possible different PJ
    !******************************************************************************!
    SUBROUTINE load_output_adapt_pj
        IMPLICIT NONE
        INTEGER :: pmaxloop_e, pmaxloop_i, jmaxloop_e, jmaxloop_i, Nkx_cp, Nky_cp, Nz_cp

        ! Checkpoint filename
        WRITE(rstfile,'(a,a1,i2.2,a3)') TRIM(resfile0),'_',job2load,'.h5'

        IF (my_id .EQ. 0) WRITE(*,'(3x,a)') "Resume from ", rstfile
        ! Open file
        CALL openf(rstfile, fidrst,mpicomm=comm0)
        ! Get grid info
        CALL getatt(fidrst,"/data/input/" , "Nkx", Nkx_cp)
        CALL getatt(fidrst,"/data/input/" , "Nky", Nky_cp)
        CALL getatt(fidrst,"/data/input/" ,  "Nz",  Nz_cp)
        !!!!!!!!! Load electron moments
        IF (KIN_E) THEN
        ! Get the checkpoint moments degrees to allocate memory
        CALL getatt(fidrst,"/data/input/" , "pmaxe", pmaxe_cp)
        CALL getatt(fidrst,"/data/input/" , "jmaxe", jmaxe_cp)
        IF (my_id .EQ. 0) WRITE(*,*) "Pe_cp = ", pmaxe_cp
        IF (my_id .EQ. 0) WRITE(*,*) "Je_cp = ", jmaxe_cp
        CALL getatt(fidrst,"/data/input/" , "start_iframe5d", n0)
        ! Allocate the required size to load checkpoints moments
        ! CALL allocate_array(moments_e_cp, 1,pmaxe_cp+1, 1,jmaxe_cp+1, ikxs,ikxe, ikys,ikye, izs,ize)
        CALL allocate_array(moments_e_cp, 1,pmaxe_cp+1, 1,jmaxe_cp+1, 1,Nky_cp, 1,Nkx_cp, 1,Nz_cp)
        ! Find the last results of the checkpoint file by iteration
        n_ = n0+1
        WRITE(dset_name, "(A, '/', i6.6)") "/data/var5d/moments_e", n_ ! start with moments_e/000001
        DO WHILE (isdataset(fidrst, dset_name)) ! If n_ is not a file we stop the loop
        n_ = n_ + 1
        WRITE(dset_name, "(A, '/', i6.6)") "/data/var5d/moments_e", n_ ! updtate file number
        ENDDO
        n_ = n_ - 1 ! n_ is not a file so take the previous one n_-1
        ! Read state of system from checkpoint file and load every moment to change the distribution
        WRITE(dset_name, "(A, '/', i6.6)") "/data/var5d/moments_e", n_
        ! CALL getarrnd(fidrst, dset_name, moments_e_cp(1:pmaxe_cp+1, 1:jmaxe_cp+1, ikys:ikye, ikxs:ikxe, izs:ize),(/1,3/))
        CALL getarr(fidrst, dset_name, moments_e_cp(1:pmaxe_cp+1, 1:jmaxe_cp+1, 1:Nky_cp, 1:Nkx_cp, 1:Nz_cp))
        ! Initialize simulation moments array with checkpoints ones
        ! (they may have a larger number of polynomials, set to 0 at the begining)
        moments_e = 0._dp;
        pmaxloop_e = min(ipe_e,pmaxe_cp+1)
        jmaxloop_e = min(ije_e,jmaxe_cp+1)
        IF (ips_e .LE. pmaxe_cp+1) THEN
          DO ip=ips_e,pmaxloop_e
          IF (ijs_e .LE. jmaxe_cp+1) THEN
            DO ij=ijs_e,jmaxloop_e
                DO ikx=ikxs,ikxe
                DO iky=ikys,ikye
                DO iz = izs,ize
                    moments_e(ip,ij,iky,ikx,iz,:) = moments_e_cp(ip,ij,iky,ikx,iz)
                ENDDO
                ENDDO
                ENDDO
            ENDDO
          ENDIF
          ENDDO
        ENDIF
        ! Deallocate checkpoint arrays
        DEALLOCATE(moments_e_cp)
        ENDIF
        !!!!!!! Load ion moments
        ! Get the checkpoint moments degrees to allocate memory
        CALL getatt(fidrst,"/data/input/" , "pmaxi", pmaxi_cp)
        CALL getatt(fidrst,"/data/input/" , "jmaxi", jmaxi_cp)
        IF (my_id .EQ. 0) WRITE(*,*) "Pi_cp = ", pmaxi_cp
        IF (my_id .EQ. 0) WRITE(*,*) "Ji_cp = ", jmaxi_cp
        CALL getatt(fidrst,"/data/input/" , "start_iframe5d", n0)
        ! Allocate the required size to load checkpoints moments
        ! CALL allocate_array(moments_i_cp, 1,pmaxi_cp+1, 1,jmaxi_cp+1, ikxs,ikxe, ikys,ikye, izs,ize)
        CALL allocate_array(moments_i_cp, 1,pmaxi_cp+1, 1,jmaxi_cp+1, 1,Nky_cp, 1,Nkx_cp, 1,Nz_cp)
        ! Find the last results of the checkpoint file by iteration
        n_ = n0+1
        WRITE(dset_name, "(A, '/', i6.6)") "/data/var5d/moments_i", n_ ! start with moments_e/000001
        DO WHILE (isdataset(fidrst, dset_name)) ! If n_ is not a file we stop the loop
        n_ = n_ + 1
        WRITE(dset_name, "(A, '/', i6.6)") "/data/var5d/moments_i", n_ ! updtate file number
        ENDDO
        n_ = n_ - 1 ! n_ is not a file so take the previous one n_-1

        ! Read state of system from checkpoint file and load every moment to change the distribution
        WRITE(dset_name, "(A, '/', i6.6)") "/data/var5d/moments_i", n_
        ! CALL getarrnd(fidrst, dset_name, moments_i_cp(1:pmaxi_cp+1, 1:jmaxi_cp+1, ikys:ikye, ikxs:ikxe, izs:ize),(/1,3/))
        CALL getarr(fidrst, dset_name, moments_i_cp(1:pmaxi_cp+1, 1:jmaxi_cp+1, 1:Nky_cp, 1:Nkx_cp, 1:Nz_cp))

        ! Initialize simulation moments array with checkpoints ones
        ! (they may have a larger number of polynomials, set to 0 at the begining)
        moments_i = 0._dp;
        pmaxloop_i = min(ipe_i,pmaxi_cp+1)
        jmaxloop_i = min(ije_i,jmaxi_cp+1)
        IF (ips_i .LE. pmaxi_cp+1) THEN
          DO ip=ips_i,pmaxloop_i
          IF (ijs_i .LE. jmaxi_cp+1) THEN
            DO ij=ijs_i,jmaxloop_i
                DO ikx=ikxs,ikxe
                DO iky=ikys,ikye
                DO iz = izs,ize
                    moments_i(ip,ij,iky,ikx,iz,:) = moments_i_cp(ip,ij,iky,ikx,iz)
                ENDDO
                ENDDO
                ENDDO
            ENDDO
          ENDIF
          ENDDO
        ENDIF
        ! Deallocate checkpoint arrays
        DEALLOCATE(moments_i_cp)

        ! Read time dependent attributes to continue simulation
        CALL getatt(fidrst, dset_name, 'cstep', cstep)
        CALL getatt(fidrst, dset_name, 'time', time)
        CALL getatt(fidrst, dset_name, 'jobnum', jobnum)
        jobnum = jobnum+1
        CALL getatt(fidrst, dset_name, 'iframe2d',iframe2d)
        CALL getatt(fidrst, dset_name, 'iframe5d',iframe5d)
        iframe2d = iframe2d-1; iframe5d = iframe5d-1

        CALL closef(fidrst)

        IF (my_id .EQ. 0) WRITE(*,'(3x,a)') "Reading from restart file "//TRIM(rstfile)//" completed!"

    END SUBROUTINE load_output_adapt_pj
    !******************************************************************************!

END MODULE restarts
