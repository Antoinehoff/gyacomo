!******************************************************************************!
!!!!!! initialize the moments and load/build coeff tables
!******************************************************************************!
SUBROUTINE inital

  USE basic
  USE model, ONLY : CO, NON_LIN
  USE prec_const
  USE coeff
  USE time_integration
  USE array, ONLY : Sepj,Sipj
  USE collision
  USE closure

  implicit none

  CALL set_updatetlevel(1)

  IF (my_id .EQ. 0) WRITE(*,*) 'Evaluate kernels'
  CALL evaluate_kernels

  !!!!!! Set the moments arrays Nepj, Nipj !!!!!!
  ! WRITE(*,*) 'Init moments'
  IF ( RESTART ) THEN
    CALL load_output
    ! CALL load_cp
  ELSE
    CALL init_moments
    !!!!!! Set phi !!!!!!
  ENDIF
  IF (my_id .EQ. 0) WRITE(*,*) 'Init phi'
  CALL poisson

  !!!!!! Set Sepj, Sipj and dnjs coeff table !!!!!!
  IF ( NON_LIN ) THEN;
    IF (my_id .EQ. 0) WRITE(*,*) 'Init Sapj'
    CALL compute_Sapj
    ! WRITE(*,*) 'Building Dnjs table'
    CALL build_dnjs_table
  ENDIF

  !!!!!! Load the full coulomb collision operator coefficients !!!!!!
  IF (CO .EQ. -1) THEN
    IF (my_id .EQ. 0) WRITE(*,*) '=== Load Full Coulomb matrix ==='
    CALL load_FC_mat
    IF (my_id .EQ. 0) WRITE(*,*) '..done'
  ENDIF

END SUBROUTINE inital
!******************************************************************************!

!******************************************************************************!
!!!!!!! Initialize the moments randomly
!******************************************************************************!
SUBROUTINE init_moments
  USE basic
  USE grid
  USE fields
  USE prec_const
  USE utility, ONLY: checkfield
  USE initial_par
  USE model, ONLY : NON_LIN
  IMPLICIT NONE

  REAL(dp) :: noise
  REAL(dp) :: kr, kz, sigma, gain, kz_shift
  INTEGER, DIMENSION(12) :: iseedarr

  ! Seed random number generator
  iseedarr(:)=iseed
  CALL RANDOM_SEED(PUT=iseedarr+my_id)

    !**** Broad noise initialization *******************************************
    DO ip=ips_e,ipe_e
      DO ij=ijs_e,ije_e

        DO ikr=ikrs,ikre
          DO ikz=ikzs,ikze
            CALL RANDOM_NUMBER(noise)
            moments_e( ip,ij, ikr,ikz, :) = (initback_moments + initnoise_moments*(noise-0.5_dp))
          END DO
        END DO

        IF ( ikrs .EQ. 1 ) THEN
          DO ikz=2,Nkz/2 !symmetry at kr = 0
            moments_e( ip,ij,1,ikz, :) = moments_e( ip,ij,1,Nkz+2-ikz, :)
          END DO
        ENDIF

      END DO
    END DO

    DO ip=ips_i,ipe_i
      DO ij=ijs_i,ije_i

        DO ikr=ikrs,ikre
          DO ikz=ikzs,ikze
            CALL RANDOM_NUMBER(noise)
            moments_i( ip,ij, ikr,ikz, :) = (initback_moments + initnoise_moments*(noise-0.5_dp))
          END DO
        END DO

        IF ( ikrs .EQ. 1 ) THEN
          DO ikz=2,Nkz/2 !symmetry at kr = 0
            moments_i( ip,ij,1,ikz, :) = moments_i( ip,ij,1,Nkz+2-ikz, :)
          END DO
        ENDIF

      END DO
    END DO

    ! Putting to zero modes that are not in the 2/3 Orszag rule
    IF (NON_LIN) THEN
      DO ip=ips_e,ipe_e
      DO ij=ijs_e,ije_e
      DO ikr=ikrs,ikre
      DO ikz=ikzs,ikze
        moments_e( ip,ij,ikr,ikz, :) = moments_e( ip,ij,ikr,ikz, :)*AA_r(ikr)*AA_z(ikz)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      DO ip=ips_i,ipe_i
      DO ij=ijs_i,ije_i
      DO ikr=ikrs,ikre
      DO ikz=ikzs,ikze
        moments_i( ip,ij,ikr,ikz, :) = moments_i( ip,ij,ikr,ikz, :)*AA_r(ikr)*AA_z(ikz)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDIF

END SUBROUTINE init_moments
!******************************************************************************!

!******************************************************************************!
!!!!!!! Load moments from a previous output file
!******************************************************************************!
SUBROUTINE load_output
  USE basic
  USE futils,          ONLY: openf, closef, getarr, getatt, isgroup, isdataset
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
  ! Checkpoint filename
  WRITE(rstfile,'(a,a1,i2.2,a3)') TRIM(resfile0),'_',job2load,'.h5'

  IF (my_id .EQ. 0) WRITE(*,'(3x,a)') "Resume from ", rstfile
  ! Open file
  CALL openf(rstfile, fidrst,mpicomm=MPI_COMM_WORLD)
  ! Get the checkpoint moments degrees to allocate memory
  CALL getatt(fidrst,"/data/input/" , "pmaxe", pmaxe_cp)
  CALL getatt(fidrst,"/data/input/" , "jmaxe", jmaxe_cp)
  CALL getatt(fidrst,"/data/input/" , "pmaxi", pmaxi_cp)
  CALL getatt(fidrst,"/data/input/" , "jmaxi", jmaxi_cp)
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
  CALL getarr(fidrst, dset_name, moments_e_cp(1:pmaxe_cp+1, 1:jmaxe_cp+1, ikrs:ikre, ikzs:ikze),pardim=3)
  WRITE(dset_name, "(A, '/', i6.6)") "/data/var5d/moments_i", n_
  CALL getarr(fidrst, dset_name, moments_i_cp(1:pmaxi_cp+1, 1:jmaxi_cp+1, ikrs:ikre, ikzs:ikze),pardim=3)

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

END SUBROUTINE load_output
!******************************************************************************!

!******************************************************************************!
!!!!!!! Load moments from a previous save
!******************************************************************************!
SUBROUTINE load_cp
  USE basic
  USE futils,          ONLY: openf, closef, getarr, getatt, isgroup, isdataset
  USE grid
  USE fields
  USE diagnostics_par
  USE time_integration
  IMPLICIT NONE

  INTEGER :: rank, sz_, n_
  INTEGER :: dims(1) = (/0/)
  CHARACTER(LEN=50) :: dset_name
  INTEGER :: pmaxe_cp, jmaxe_cp, pmaxi_cp, jmaxi_cp
  COMPLEX(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: moments_e_cp
  COMPLEX(dp), DIMENSION(:,:,:,:), ALLOCATABLE :: moments_i_cp
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


!******************************************************************************!
!!!!!!! Build the Laguerre-Laguerre coupling coefficient table for nonlin
!******************************************************************************!
SUBROUTINE build_dnjs_table
  USE basic
  USE array, Only : dnjs
  USE grid, Only : jmaxe, jmaxi
  USE coeff
  IMPLICIT NONE

  INTEGER :: in, ij, is, J
  INTEGER :: n_, j_, s_

  J = max(jmaxe,jmaxi)

  DO in = 1,J+1 ! Nested dependent loops to make benefit from dnjs symmetry
    n_ = in - 1
    DO ij = in,J+1
      j_ = ij - 1
      DO is = ij,J+1
        s_ = is - 1

        dnjs(in,ij,is) = TO_DP(ALL2L(n_,j_,s_,0))
        ! By symmetry
        dnjs(in,is,ij) = dnjs(in,ij,is)
        dnjs(ij,in,is) = dnjs(in,ij,is)
        dnjs(ij,is,in) = dnjs(in,ij,is)
        dnjs(is,ij,in) = dnjs(in,ij,is)
        dnjs(is,in,ij) = dnjs(in,ij,is)
      ENDDO
    ENDDO
  ENDDO
END SUBROUTINE build_dnjs_table
!******************************************************************************!

!******************************************************************************!
!!!!!!! Evaluate the kernels once for all
!******************************************************************************!
SUBROUTINE evaluate_kernels
  USE basic
  USE array, Only : kernel_e, kernel_i
  USE grid
  use model, ONLY : tau_e, tau_i, sigma_e, sigma_i, q_e, q_i, lambdaD, CLOS
  IMPLICIT NONE

  REAL(dp)    :: factj, j_dp, j_int
  REAL(dp)    :: sigmae2_taue_o2, sigmai2_taui_o2
  REAL(dp)    :: be_2, bi_2, alphaD
  REAL(dp)    :: kr, kz, kperp2

  !!!!! Electron kernels !!!!!
  !Precompute species dependant factors
  sigmae2_taue_o2 = sigma_e**2 * tau_e/2._dp ! factor of the Kernel argument

  factj = 1.0 ! Start of the recursive factorial
  DO ij = 1, jmaxe+1
    j_int = jarray_e(ij)
    j_dp = REAL(j_int,dp) ! REAL of degree

    ! Recursive factorial
    IF (j_dp .GT. 0) THEN
      factj = factj * j_dp
    ELSE
      factj = 1._dp
    ENDIF

    DO ikr = ikrs,ikre
      kr     = krarray(ikr)
      DO ikz = ikzs,ikze
        kz    = kzarray(ikz)

        be_2  =  (kr**2 + kz**2) * sigmae2_taue_o2

        kernel_e(ij, ikr, ikz) = be_2**j_int * exp(-be_2)/factj

      ENDDO
    ENDDO
  ENDDO
  ! Kernels closure
  DO ikr = ikrs,ikre
    kr     = krarray(ikr)
    DO ikz = ikzs,ikze
      kz    = kzarray(ikz)
      be_2  =  (kr**2 + kz**2) * sigmae2_taue_o2
      kernel_e(ijeg_e,ikr,ikz) = be_2/(real(ijeg_e,dp))*kernel_e(ije_e,ikr,ikz)
    ENDDO
  ENDDO

  !!!!! Ion kernels !!!!!
  sigmai2_taui_o2 = sigma_i**2 * tau_i/2._dp ! (b_a/2)^2 = (kperp sqrt(2 tau_a) sigma_a/2)^2

  factj = 1.0 ! Start of the recursive factorial
  DO ij = 1, jmaxi+1
    j_int = jarray_e(ij)
    j_dp = REAL(j_int,dp) ! REAL of degree

    ! Recursive factorial
    IF (j_dp .GT. 0) THEN
      factj = factj * j_dp
    ELSE
      factj = 1._dp
    ENDIF

    DO ikr = ikrs,ikre
      kr     = krarray(ikr)
      DO ikz = ikzs,ikze
        kz    = kzarray(ikz)

        bi_2  =  (kr**2 + kz**2) * sigmai2_taui_o2

        kernel_i(ij, ikr, ikz) = bi_2**j_int * exp(-bi_2)/factj

      ENDDO
    ENDDO
  ENDDO
  ! Kernels closure
  DO ikr = ikrs,ikre
    kr     = krarray(ikr)
    DO ikz = ikzs,ikze
      kz    = kzarray(ikz)
      bi_2  =  (kr**2 + kz**2) * sigmai2_taui_o2
      kernel_i(ijeg_i,ikr,ikz) = bi_2/(real(ijeg_i,dp))*kernel_e(ije_i,ikr,ikz)
    ENDDO
  ENDDO
END SUBROUTINE evaluate_kernels
!******************************************************************************!
