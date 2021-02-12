module closure
! Contains the routines to define closures
USE basic
USE model,  ONLY: CLOS
USE grid
USE array,  ONLY: kernel_e,  kernel_i
USE fields, ONLY: moments_e, moments_i

IMPLICIT NONE

PUBLIC :: apply_closure_model

CONTAINS

  SUBROUTINE apply_closure_model
    IMPLICIT NONE

    ! Negative out of bounds indices are put to zero (analytically correct)
    DO ikr = ikrs,ikre
      DO ikz = ikzs,ikze

        DO ip = ipsg_e,ipeg_e
          moments_e(ip,ijs_e-1,ikr,ikz,:) = 0._dp
        ENDDO
        DO ij = ijsg_e,ijeg_e
          moments_e(ips_e-1,ij,ikr,ikz,:) = 0._dp
          moments_e(ips_e-2,ij,ikr,ikz,:) = 0._dp
        ENDDO
        kernel_e(ijs_e-1,ikr,ikz)      = 0._dp

        DO ip = ipsg_i,ipeg_i
          moments_i(ip,ijs_i-1,ikr,ikz,:) = 0._dp
        ENDDO
        DO ij = ijsg_i,ijeg_i
          moments_i(ips_i-1,ij,ikr,ikz,:) = 0._dp
          moments_i(ips_i-2,ij,ikr,ikz,:) = 0._dp
        ENDDO
        kernel_i(ijs_i-1,ikr,ikz)      = 0._dp

      ENDDO
    ENDDO
    ! Positive Oob indices are approximated with a model
    IF (CLOS .EQ. 0) THEN
      ! zero truncation, An+1=0 for n+1>nmax
      DO ikr = ikrs,ikre
        DO ikz = ikzs,ikze

          DO ip = ipsg_e,ipeg_e
            moments_e(ip,ije_e+1,ikr,ikz,:) = 0._dp
          ENDDO
          DO ij = ijsg_e,ijeg_e
            moments_e(ipe_e+1,ij,ikr,ikz,:) = 0._dp
            moments_e(ipe_e+2,ij,ikr,ikz,:) = 0._dp
          ENDDO
          kernel_e(ije_e+1,ikr,ikz)      = 0._dp

          DO ip = ipsg_i,ipeg_i
            moments_i(ip,ije_i+1,ikr,ikz,:) = 0._dp
          ENDDO
          DO ij = ijsg_i,ijeg_i
            moments_i(ipe_i+1,ij,ikr,ikz,:) = 0._dp
            moments_i(ipe_i+2,ij,ikr,ikz,:) = 0._dp
          ENDDO
          kernel_i(ije_i+1,ikr,ikz)      = 0._dp
          
        ENDDO
      ENDDO

    ELSEIF (CLOS .EQ. 1) THEN
      ! Copy truncation with n+1 = min(nmax,n+1)
      ! here pmax+1 and pmax_2 are mapped to pmax
      moments_e(ipe_e+1,:,:,:,:) = moments_e(ipe_e,:,:,:,:)
      moments_e(ipe_e+2,:,:,:,:) = moments_e(ipe_e,:,:,:,:)
      moments_e(:,ije_e+1,:,:,:) = moments_e(:,ije_e,:,:,:)
      kernel_e(ije_e+1,:,:)      = kernel_e(ije_e,:,:)

      moments_i(ipe_i+1,:,:,:,:) = moments_i(ipe_i,:,:,:,:)
      moments_i(ipe_i+2,:,:,:,:) = moments_i(ipe_i,:,:,:,:)
      moments_i(:,ije_i+1,:,:,:) = moments_i(:,ije_i,:,:,:)
      kernel_i(ije_i+1,:,:)      = kernel_i(ije_i,:,:)


    ELSEIF (CLOS .EQ. 2) THEN
      ! Copy truncation with special treatment for Hermite
      ! here pmax+1 is mapped to pmax-1 and pmax+2 to pmax
      moments_e(ipe_e+1,:,:,:,:) = moments_e(ipe_e-1,:,:,:,:)
      moments_e(ipe_e+2,:,:,:,:) = moments_e(ipe_e,:,:,:,:)
      moments_e(:,ije_e+1,:,:,:) = moments_e(:,ije_e,:,:,:)
      kernel_e(ije_e+1,:,:)      = kernel_e(ije_e,:,:)

      moments_i(ipe_i+1,:,:,:,:) = moments_i(ipe_i-1,:,:,:,:)
      moments_i(ipe_i+2,:,:,:,:) = moments_i(ipe_i,:,:,:,:)
      moments_i(:,ije_i+1,:,:,:) = moments_i(:,ije_i,:,:,:)
      kernel_i(ije_i+1,:,:)      = kernel_i(ije_i,:,:)
    ENDIF

  END SUBROUTINE apply_closure_model
END module closure
