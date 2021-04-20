MODULE restarts
USE basic
USE futils,          ONLY: openf, closef, getarr, getatt, isgroup, isdataset,getarrnd
USE grid
USE fields
USE diagnostics_par
USE time_integration

IMPLICIT NONE

INTEGER :: rank, sz_, n_
INTEGER :: dims(1) = (/0/)
CHARACTER(LEN=50) :: dset_name
INTEGER :: pmaxe_cp, jmaxe_cp, pmaxi_cp, jmaxi_cp, n0
COMPLEX(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: moments_e_cp
COMPLEX(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: moments_i_cp

PUBLIC :: load_moments

CONTAINS


    SUBROUTINE load_moments
        CALL load_output_same_dims ! load same dimensions older moments from output file
        ! CALL load_output_adapt_pj  ! load moments with possibly different PJ from output file
        ! CALL load_cp             ! load from checkpoint file (meant to be deleted)
    END SUBROUTINE load_moments

        !******************************************************************************!
    !!!!!!! Load moments from a previous output file with same PJ
    !******************************************************************************!
    SUBROUTINE load_output_same_dims
        IMPLICIT NONE

        ! Checkpoint filename
        WRITE(rstfile,'(a,a1,i2.2,a3)') TRIM(resfile0),'_',job2load,'.h5'

        IF (my_id .EQ. 0) WRITE(*,'(3x,a)') "Resume from ", rstfile
        ! Open file
        CALL openf(rstfile, fidrst,mpicomm=comm0)
        ! Get the checkpoint moments degrees to allocate memory
        CALL getatt(fidrst,"/data/input/" , "pmaxe", pmaxe_cp)
        CALL getatt(fidrst,"/data/input/" , "jmaxe", jmaxe_cp)
        CALL getatt(fidrst,"/data/input/" , "pmaxi", pmaxi_cp)
        CALL getatt(fidrst,"/data/input/" , "jmaxi", jmaxi_cp)
        IF (my_id .EQ. 0) WRITE(*,*) "Pe_cp = ", pmaxe_cp
        IF (my_id .EQ. 0) WRITE(*,*) "Je_cp = ", jmaxe_cp
        CALL getatt(fidrst,"/data/input/" , "start_iframe5d", n0)

        IF ((pmaxe_cp .NE. pmaxe) .OR. (jmaxe_cp .NE. jmaxe) .OR.&
         (pmaxi_cp .NE. pmaxi) .OR. (jmaxi_cp .NE. jmaxi)) THEN
         WRITE(*,*) '! Previous simulation has not the same polynomial basis ! -> EXIT'
        ENDIF

        ! Find the last results of the checkpoint file by iteration
        n_ = n0+1
        WRITE(dset_name, "(A, '/', i6.6)") "/data/var5d/moments_e", n_ ! start with moments_e/000001
        DO WHILE (isdataset(fidrst, dset_name)) ! If n_ is not a file we stop the loop
        n_ = n_ + 1
        WRITE(dset_name, "(A, '/', i6.6)") "/data/var5d/moments_e", n_ ! updtate file number
        ENDDO
        n_ = n_ - 1 ! n_ is not a file so take the previous one n_-1

        ! Read time dependent attributes to continue simulation
        WRITE(dset_name, "(A, '/', i6.6)") "/data/var5d/moments_e", n_
        CALL getatt(fidrst, dset_name, 'cstep', cstep)
        CALL getatt(fidrst, dset_name, 'time', time)
        CALL getatt(fidrst, dset_name, 'jobnum', jobnum)
        jobnum = jobnum+1
        CALL getatt(fidrst, dset_name, 'iframe2d',iframe2d)
        CALL getatt(fidrst, dset_name, 'iframe5d',iframe5d)
        iframe2d = iframe2d-1; iframe5d = iframe5d-1
        IF(my_id.EQ.0) WRITE(*,*) '.. restart from t = ', time

        ! Read state of system from checkpoint file
        WRITE(dset_name, "(A, '/', i6.6)") "/data/var5d/moments_e", n_
        CALL getarrnd(fidrst, dset_name, moments_e(ips_e:ipe_e, ijs_e:ije_e, ikrs:ikre, ikzs:ikze, 1),(/1,3/))
        WRITE(dset_name, "(A, '/', i6.6)") "/data/var5d/moments_i", n_
        CALL getarrnd(fidrst, dset_name, moments_i(ips_i:ipe_i, ijs_i:ije_i, ikrs:ikre, ikzs:ikze, 1),(/1,3/))

        CALL closef(fidrst)

        IF (my_id .EQ. 0) WRITE(*,'(3x,a)') "Reading from restart file "//TRIM(rstfile)//" completed!"

    END SUBROUTINE load_output_same_dims
    !******************************************************************************!


    !******************************************************************************!
    !!!!!!! Load moments from a previous output file with possible different PJ
    !******************************************************************************!
    SUBROUTINE load_output_adapt_pj
        IMPLICIT NONE

        ! Checkpoint filename
        WRITE(rstfile,'(a,a1,i2.2,a3)') TRIM(resfile0),'_',job2load,'.h5'

        IF (my_id .EQ. 0) WRITE(*,'(3x,a)') "Resume from ", rstfile
        ! Open file
        CALL openf(rstfile, fidrst,mpicomm=comm0)
        ! Get the checkpoint moments degrees to allocate memory
        CALL getatt(fidrst,"/data/input/" , "pmaxe", pmaxe_cp)
        CALL getatt(fidrst,"/data/input/" , "jmaxe", jmaxe_cp)
        CALL getatt(fidrst,"/data/input/" , "pmaxi", pmaxi_cp)
        CALL getatt(fidrst,"/data/input/" , "jmaxi", jmaxi_cp)
        IF (my_id .EQ. 0) WRITE(*,*) "Pe_cp = ", pmaxe_cp
        IF (my_id .EQ. 0) WRITE(*,*) "Je_cp = ", jmaxe_cp
        CALL getatt(fidrst,"/data/input/" , "start_iframe5d", n0)

        ! Allocate the required size to load checkpoints moments
        CALL allocate_array(moments_e_cp, 1,pmaxe_cp+1, 1,jmaxe_cp+1, ikrs,ikre, ikzs,ikze)
        CALL allocate_array(moments_i_cp, 1,pmaxi_cp+1, 1,jmaxi_cp+1, ikrs,ikre, ikzs,ikze)
        ! Find the last results of the checkpoint file by iteration
        n_ = n0+1
        WRITE(dset_name, "(A, '/', i6.6)") "/data/var5d/moments_e", n_ ! start with moments_e/000001
        DO WHILE (isdataset(fidrst, dset_name)) ! If n_ is not a file we stop the loop
        n_ = n_ + 1
        WRITE(dset_name, "(A, '/', i6.6)") "/data/var5d/moments_e", n_ ! updtate file number
        ENDDO
        n_ = n_ - 1 ! n_ is not a file so take the previous one n_-1

        ! Read state of system from checkpoint file
        WRITE(dset_name, "(A, '/', i6.6)") "/data/var5d/moments_e", n_
        CALL getarrnd(fidrst, dset_name, moments_e_cp(1:pmaxe_cp+1, 1:jmaxe_cp+1, ikrs:ikre, ikzs:ikze),(/1,3/))
        WRITE(dset_name, "(A, '/', i6.6)") "/data/var5d/moments_i", n_
        CALL getarrnd(fidrst, dset_name, moments_i_cp(1:pmaxi_cp+1, 1:jmaxi_cp+1, ikrs:ikre, ikzs:ikze),(/1,3/))

        ! Initialize simulation moments array with checkpoints ones
        ! (they may have a larger number of polynomials, set to 0 at the begining)
        moments_e = 0._dp; moments_i = 0._dp
        DO ip=1,pmaxe_cp+1
        DO ij=1,jmaxe_cp+1
            DO ikr=ikrs,ikre
            DO ikz=ikzs,ikze
                moments_e(ip,ij,ikr,ikz,:) = moments_e_cp(ip,ij,ikr,ikz)
            ENDDO
            ENDDO
        ENDDO
        ENDDO

        DO ip=1,pmaxi_cp+1
        DO ij=1,jmaxi_cp+1
            DO ikr=ikrs,ikre
            DO ikz=ikzs,ikze
                moments_i(ip,ij,ikr,ikz,:) = moments_i_cp(ip,ij,ikr,ikz)
            ENDDO
            ENDDO
        ENDDO
        ENDDO
        ! Deallocate checkpoint arrays
        DEALLOCATE(moments_e_cp)
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

    !******************************************************************************!
    !!!!!!! Load moments from a previous save
    !******************************************************************************!
    SUBROUTINE load_cp
        IMPLICIT NONE

        ! Checkpoint filename
        WRITE(rstfile,'(a,a1,i2.2,a3)') TRIM(rstfile0),'_',job2load,'.h5'

        IF (my_id .EQ. 0) WRITE(*,'(3x,a)') "Resume from previous run"
        ! Open file
        CALL openf(rstfile, fidrst,mpicomm=MPI_COMM_WORLD)
        ! Get the checkpoint moments degrees to allocate memory
        CALL getatt(fidrst,"/Basic/moments_e/" , "pmaxe", pmaxe_cp)
        CALL getatt(fidrst,"/Basic/moments_e/" , "jmaxe", jmaxe_cp)
        CALL getatt(fidrst,"/Basic/moments_i/" , "pmaxi", pmaxi_cp)
        CALL getatt(fidrst,"/Basic/moments_i/" , "jmaxi", jmaxi_cp)
        IF (my_id .EQ. 0) WRITE(*,*) "Pe_cp = ", pmaxe_cp
        IF (my_id .EQ. 0) WRITE(*,*) "Je_cp = ", jmaxe_cp

        ! Allocate the required size to load checkpoints moments
        CALL allocate_array(moments_e_cp, 1,pmaxe_cp+1, 1,jmaxe_cp+1, ikrs,ikre, ikzs,ikze)
        CALL allocate_array(moments_i_cp, 1,pmaxi_cp+1, 1,jmaxi_cp+1, ikrs,ikre, ikzs,ikze)
        ! Find the last results of the checkpoint file by iteration
        n_ = 0
        WRITE(dset_name, "(A, '/', i6.6)") "/Basic/moments_e", n_ ! start with moments_e/000000
        DO WHILE (isdataset(fidrst, dset_name)) ! If n_ is not a file we stop the loop
        n_ = n_ + 1
        WRITE(dset_name, "(A, '/', i6.6)") "/Basic/moments_e", n_ ! updtate file number
        ENDDO
        n_ = n_ - 1 ! n_ is not a file so take the previous one n_-1

        ! Read state of system from checkpoint file
        WRITE(dset_name, "(A, '/', i6.6)") "/Basic/moments_e", n_
        CALL getarr(fidrst, dset_name, moments_e_cp(1:pmaxe_cp+1, 1:jmaxe_cp+1, ikrs:ikre, ikzs:ikze),pardim=3)
        WRITE(dset_name, "(A, '/', i6.6)") "/Basic/moments_i", n_
        CALL getarr(fidrst, dset_name, moments_i_cp(1:pmaxi_cp+1, 1:jmaxi_cp+1, ikrs:ikre, ikzs:ikze),pardim=3)
        WRITE(dset_name, "(A, '/', i6.6)") "/Basic/phi", n_
        CALL getarr(fidrst, dset_name, phi(ikrs:ikre,ikzs:ikze),pardim=1)

        ! Initialize simulation moments array with checkpoints ones
        ! (they may have a larger number of polynomials, set to 0 at the begining)
        moments_e = 0._dp; moments_i = 0._dp
        DO ip=1,pmaxe_cp+1
        DO ij=1,jmaxe_cp+1
            DO ikr=ikrs,ikre
            DO ikz=ikzs,ikze
                moments_e(ip,ij,ikr,ikz,:) = moments_e_cp(ip,ij,ikr,ikz)
            ENDDO
            ENDDO
        ENDDO
        ENDDO

        DO ip=1,pmaxi_cp+1
        DO ij=1,jmaxi_cp+1
            DO ikr=ikrs,ikre
            DO ikz=ikzs,ikze
                moments_i(ip,ij,ikr,ikz,:) = moments_i_cp(ip,ij,ikr,ikz)
            ENDDO
            ENDDO
        ENDDO
        ENDDO
        ! Deallocate checkpoint arrays
        DEALLOCATE(moments_e_cp)
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

    END SUBROUTINE load_cp
    !******************************************************************************!

END MODULE restarts
