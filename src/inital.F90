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

  implicit none

  real :: start_init, end_init, time_estimation

  CALL set_updatetlevel(1)

  CALL evaluate_kernels

  !!!!!! Set the moments arrays Nepj, Nipj !!!!!!
  ! WRITE(*,*) 'Init moments'
  IF ( RESTART ) THEN
    CALL load_cp
  ELSE
    CALL init_moments
    !!!!!! Set phi !!!!!!
    IF (my_id .EQ. 0) WRITE(*,*) 'Init phi'
    CALL poisson
  ENDIF


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
            moments_e( ip,ij,     ikr,    ikz, :) = (initback_moments + initnoise_moments*(noise-0.5_dp))*AA_r(ikr)*AA_z(ikz)
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
            moments_i( ip,ij,ikr,ikz, :) = (initback_moments + initnoise_moments*(noise-0.5_dp))*AA_r(ikr)*AA_z(ikz)
          END DO
        END DO

        IF ( ikrs .EQ. 1 ) THEN
          DO ikz=2,Nkz/2 !symmetry at kr = 0
            moments_i( ip,ij,1,ikz, :) = moments_i( ip,ij,1,Nkz+2-ikz, :)
          END DO
        ENDIF

      END DO
    END DO

  ! ENDIF
END SUBROUTINE init_moments
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
  INTEGER ::  dims(1) = (/0/)
  CHARACTER(LEN=50) :: dset_name

  WRITE(rstfile,'(a,a1,i2.2,a3)') TRIM(rstfile0),'_',job2load,'.h5'

  IF (my_id .EQ. 0) WRITE(*,'(3x,a)') "Resume from previous run"

  CALL openf(rstfile, fidrst,mpicomm=MPI_COMM_WORLD)

    n_ = 0
    WRITE(dset_name, "(A, '/', i6.6)") "/Basic/moments_e", n_
    DO WHILE (isdataset(fidrst, dset_name))
      n_ = n_ + 1
      WRITE(dset_name, "(A, '/', i6.6)") "/Basic/moments_e", n_
    ENDDO
    n_ = n_ - 1
    WRITE(dset_name, "(A, '/', i6.6)") "/Basic/moments_e", n_
    WRITE(dset_name, "(A, '/', i6.6)") "/Basic/moments_i", n_

    ! Read state of system from restart file
    WRITE(dset_name, "(A, '/', i6.6)") "/Basic/moments_i", n_
    CALL getarr(fidrst, dset_name, moments_i(ips_i:ipe_i,ijs_i:ije_i,ikrs:ikre,ikzs:ikze,1),pardim=3)
    WRITE(dset_name, "(A, '/', i6.6)") "/Basic/moments_e", n_
    CALL getarr(fidrst, dset_name, moments_e(ips_e:ipe_e,ijs_e:ije_e,ikrs:ikre,ikzs:ikze,1),pardim=3)
    WRITE(dset_name, "(A, '/', i6.6)") "/Basic/phi", n_
    CALL getarr(fidrst, dset_name, phi(ikrs:ikre,ikzs:ikze),pardim=1)

    ! Read time dependent attributes
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
  use model, ONLY : tau_e, tau_i, sigma_e, sigma_i, q_e, q_i, lambdaD, DK
  IMPLICIT NONE

  REAL(dp)    :: factj, j_dp, j_int
  REAL(dp)    :: sigmae2_taue_o2, sigmai2_taui_o2
  REAL(dp)    :: be_2, bi_2, alphaD
  REAL(dp)    :: kr, kz, kperp2

  !!!!! Electron kernels !!!!!
  !Precompute species dependant factors
  sigmae2_taue_o2 = sigma_e**2 * tau_e/2._dp ! factor of the Kernel argument

  factj = 1.0 ! Start of the recursive factorial
  DO ij = ijs_e, ije_e
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

        kernel_e(ij, ikr, ikz) = be_2**(ij)/factj * exp(-be_2)

      ENDDO
    ENDDO
  ENDDO

  !!!!! Ion kernels !!!!!
  sigmai2_taui_o2 = sigma_i**2 * tau_i/2._dp ! (b_a/2)^2 = (kperp sqrt(2 tau_a) sigma_a/2)^2
  factj = 1.0 ! Start of the recursive factorial
  DO ij = ijs_i, ije_i
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

        kernel_i(ij, ikr, ikz) = bi_2**(ij)/factj * exp(-bi_2)

      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE evaluate_kernels
!******************************************************************************!
