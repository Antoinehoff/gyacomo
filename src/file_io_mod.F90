MODULE file_io
  USE h5fortran
  USE mpi
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: attach, append, creatf, creatg, creatd, closef, \
            getatt, init_hdf, openf, putarr, putfile,

CONTAINS

  SUBROUTINE init_hdf()
    INTEGER(HID_T) :: error
    CALL h5open_f(error)
    if (error /= 0) then
      PRINT *, "Error opening HDF5 library"
      STOP
    end if
  END SUBROUTINE init_hdf

  SUBROUTINE attach(fid, path, name, value)
    INTEGER, INTENT(IN) :: fid
    CHARACTER(len=*), INTENT(IN) :: path, name
    CHARACTER(len=*), INTENT(IN) :: value
    INTEGER(HID_T) :: attr_id, space_id, error

    CALL h5screate_f(H5S_SCALAR_F, space_id, error)
    CALL h5acreate_f(fid, TRIM(path)//"/"//TRIM(name), H5T_NATIVE_CHARACTER, space_id, attr_id, error)
    CALL h5awrite_f(attr_id, H5T_NATIVE_CHARACTER, value, error)
    CALL h5aclose_f(attr_id, error)
    CALL h5sclose_f(space_id, error)
  END SUBROUTINE attach

  SUBROUTINE append(fid, path, value, ionode)
    INTEGER, INTENT(IN) :: fid, ionode
    CHARACTER(len=*), INTENT(IN) :: path
    REAL, INTENT(IN) :: value
    INTEGER(HID_T) :: dset_id, space_id, plist_id, error
    INTEGER(HID_T) :: file_space, mem_space
    INTEGER(HSSIZE_T) :: dims(1), maxdims(1), newdims(1)

    CALL h5dopen_f(fid, TRIM(path), dset_id, error)
    CALL h5dget_space_f(dset_id, space_id, error)
    CALL h5sget_simple_extent_dims_f(space_id, dims, maxdims, error)
    newdims = dims + 1
    CALL h5dset_extent_f(dset_id, newdims, error)
    CALL h5sclose_f(space_id, error)

    CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, value, newdims, error)
    CALL h5dclose_f(dset_id, error)
  END SUBROUTINE append

  SUBROUTINE creatf(filename, fid, real_prec, mpicomm)
    CHARACTER(len=*), INTENT(IN) :: filename
    INTEGER, INTENT(OUT) :: fid
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: real_prec
    INTEGER, INTENT(IN), OPTIONAL :: mpicomm
    INTEGER(HID_T) :: error, plist_id

    CALL h5fcreate_f(TRIM(filename), H5F_ACC_TRUNC_F, fid, error)
  END SUBROUTINE creatf

  SUBROUTINE creatg(fid, path, name)
    INTEGER, INTENT(IN) :: fid
    CHARACTER(len=*), INTENT(IN) :: path, name
    INTEGER(HID_T) :: group_id, error

    CALL h5gcreate_f(fid, TRIM(path)//"/"//TRIM(name), group_id, error)
    CALL h5gclose_f(group_id, error)
  END SUBROUTINE creatg

  SUBROUTINE creatd(fid, rank, dims, path, name)
    INTEGER, INTENT(IN) :: fid, rank
    INTEGER, DIMENSION(*), INTENT(IN) :: dims
    CHARACTER(len=*), INTENT(IN) :: path, name
    INTEGER(HID_T) :: dspace_id, dset_id, error

    CALL h5screate_simple_f(rank, dims, dspace_id, error)
    CALL h5dcreate_f(fid, TRIM(path)//"/"//TRIM(name), H5T_NATIVE_REAL, dspace_id, dset_id, error)
    CALL h5sclose_f(dspace_id, error)
    CALL h5dclose_f(dset_id, error)
  END SUBROUTINE creatd

  SUBROUTINE closef(fid)
    INTEGER, INTENT(IN) :: fid
    INTEGER(HID_T) :: error
    CALL h5fclose_f(fid, error)
  END SUBROUTINE closef

  SUBROUTINE getatt(fid, path, name, value)
    INTEGER, INTENT(IN) :: fid
    CHARACTER(len=*), INTENT(IN) :: path, name
    INTEGER, INTENT(OUT) :: value
    INTEGER(HID_T) :: attr_id, error

    CALL h5aopen_f(fid, TRIM(path)//"/"//TRIM(name), attr_id, error)
    CALL h5aread_f(attr_id, H5T_NATIVE_INTEGER, value, error)
    CALL h5aclose_f(attr_id, error)
  END SUBROUTINE getatt

  SUBROUTINE openf(filename, fid)
    CHARACTER(len=*), INTENT(IN) :: filename
    INTEGER, INTENT(OUT) :: fid
    INTEGER(HID_T) :: error
    CALL h5fopen_f(TRIM(filename), H5F_ACC_RDWR_F, fid, error)
  END SUBROUTINE openf

  SUBROUTINE putarr(fid, path, array, ionode)
    INTEGER, INTENT(IN) :: fid, ionode
    CHARACTER(len=*), INTENT(IN) :: path
    REAL, DIMENSION(:), INTENT(IN) :: array
    INTEGER(HID_T) :: dset_id, error

    CALL h5dopen_f(fid, TRIM(path), dset_id, error)
    CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, array, SHAPE(array), error)
    CALL h5dclose_f(dset_id, error)
  END SUBROUTINE putarr

  SUBROUTINE putfile(fid, path, filename, ionode)
    INTEGER, INTENT(IN) :: fid, ionode
    CHARACTER(len=*), INTENT(IN) :: path, filename
    INTEGER(HID_T) :: file_id, error

    CALL h5fopen_f(TRIM(filename), H5F_ACC_RDONLY_F, file_id, error)
    CALL h5fclose_f(file_id, error)
  END SUBROUTINE putfile

END MODULE file_io
