
module paraview
!
   use hdf5_wrapper
   use HDF5
!
   implicit none
!
   save
!
   character(len=128) :: FILE_VIS = '../../files/vis'
   character(len=2  ) :: VLEVEL = '2'
   character(len=128) :: PARAVIEW_DIR = 'vtk/'
   integer,parameter  :: PARAVIEW_IO = 22
   logical            :: PARAVIEW_DUMP_GEOM = .FALSE.
   logical            :: PARAVIEW_DUMP_ATTR = .FALSE.
   logical            :: SECOND_ORDER_VIS   = .FALSE.
!
!..Flag for switching on VTU format for output meshes and solution Fields. FALSE switches to XDMF.
   logical            :: VIS_VTU = .FALSE.
!
!  this is matching EXGEOM flag in "control" module
   integer,parameter  :: PARAVIEW_ISOGEOM = 0
   integer,parameter  :: PARAVIEW_EXGEOM  = 1
!
!  this value is set to EXGEOM flag in "read_control" routine
   integer            :: PARAVIEW_GEOM
!
!  set the domain to be output (PARAVIEW_DOMAIN = 0 means output all)
   integer            :: PARAVIEW_DOMAIN = 0
!
!  Geometry objects (elements) and points (vertices)
   integer, allocatable :: GEOM_OBJ(:)
   real(8), allocatable :: GEOM_PTS(:,:)
   real(8), allocatable :: ATTR_VAL(:,:)
!
!..Variables for VTU Format
   integer, allocatable :: IPARATTR_VTU(:)
   integer, allocatable :: ELEM_TYPES(:)
!
!
   contains
!
!-----------------------------------------------------------------------------
!> Purpose : initialize paraview output (open hdf5 file)
!!
!> @date Mar 2023
!-----------------------------------------------------------------------------
   subroutine paraview_init()
!  ...only call hdf5_w_finalize when xdmf output is chosen
      if (.not. VIS_VTU) call hdf5_w_init()
   end subroutine paraview_init



!-----------------------------------------------------------------------------
!> Purpose : finalize paraview output (close hdf5; deallocate data arrays)
!!
!> @date Mar 2023
!-----------------------------------------------------------------------------
   subroutine paraview_finalize()
!  ...only call hdf5_w_finalize when xdmf output is chosen
      if (.not. VIS_VTU) call hdf5_w_finalize()
      if (allocated(GEOM_OBJ)) deallocate(GEOM_OBJ)
      if (allocated(GEOM_PTS)) deallocate(GEOM_PTS)
      if (allocated(GEOM_PTS)) deallocate(ATTR_VAL)
   end subroutine paraview_finalize



!-----------------------------------------------------------------------------
!> Purpose : initialize data arrays for outputting paraview mesh
!!
!> @date Feb 2023
!-----------------------------------------------------------------------------
   subroutine geometry_init(Nobj,Npts)
!
      integer, intent(in) :: Nobj,Npts
!
      if (allocated(GEOM_OBJ)) deallocate(GEOM_OBJ)
      if (allocated(GEOM_PTS)) deallocate(GEOM_PTS)
      allocate(GEOM_OBJ(Nobj))
      allocate(GEOM_PTS(3,Npts))
!
   end subroutine geometry_init



!-----------------------------------------------------------------------------
!> Purpose : finalize geometry data arrays
!!
!> @date Feb 2023
!-----------------------------------------------------------------------------
   subroutine geometry_close()
      if (allocated(GEOM_OBJ)) deallocate(GEOM_OBJ)
      if (allocated(GEOM_PTS)) deallocate(GEOM_PTS)
   end subroutine geometry_close



!-----------------------------------------------------------------------------
!> Purpose : write mesh to paraview (xdmf)
!!
!> @date Feb 2023
!-----------------------------------------------------------------------------
   subroutine geometry_write(Sname,LSname,Sfile,LSfile)
!
      integer, intent(in)   :: LSname,LSfile
!
      character(len=LSname) :: Sname
      character(len=LSfile) :: Sfile
!
      integer(HID_T) :: file_id   ! File identifier
      integer(HID_T) :: dset_id   ! Dataset identifier
      integer(HID_T) :: dspace_id ! Dataspace identifier
!
      integer        :: ierr
!
!  ...Dataset dimensions and rank
      integer(HSIZE_T), dimension(2) :: pts_dims
      integer(HSIZE_T), dimension(1) :: obj_dims
      integer :: pts_rank = 2
      integer :: obj_rank = 1
!
      pts_dims(1) = size(GEOM_PTS,dim=1)
      pts_dims(2) = size(GEOM_PTS,dim=2)
      obj_dims = size(GEOM_OBJ)
!
!  -----------------------------------------------------------------------------
!
      if (.not. HDF5_IS_INIT) then
         write(*,*) 'geometry_write: HDF5_IS_INIT = .false.; returning...'
         return
      endif
!
!  ...Create file
      call h5fcreate_f(Sfile,H5F_ACC_TRUNC_F, file_id,ierr);
      call hdf5_w_handle_err(ierr,'geometry_write: h5fcreate_f')
!
!      write(*,*) 'Sname  = ', Sname
!      write(*,*) 'Sfile  = ', Sfile
!
!  -----------------------------------------------------------------------------
!  1. Write geometry file (coordinates)
!  ...Create dataspace
      call h5screate_simple_f(pts_rank,pts_dims, dspace_id,ierr)
      call hdf5_w_handle_err(ierr,'geometry_write: h5screate_simple_f')
!
!  ...Create dataset
      call h5dcreate_f(file_id,'/Coords',H5T_NATIVE_DOUBLE,dspace_id, &
                       dset_id,ierr,H5P_DEFAULT_F,H5P_DEFAULT_F,H5P_DEFAULT_F)
      call hdf5_w_handle_err(ierr,'geometry_write: h5dcreate_f')
!
!  ...Write dataset
      call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,GEOM_PTS,pts_dims, ierr)
      call hdf5_w_handle_err(ierr,'geometry_write: h5dwrite_f')
!
!  ...Close dataset
      call h5dclose_f(dset_id, ierr)
      call hdf5_w_handle_err(ierr,'geometry_write: h5dclose_f')
!
!  ...Close dataspace
      call h5sclose_f(dspace_id, ierr)
      call hdf5_w_handle_err(ierr,'geometry_write: h5sclose_f')
!
!  -----------------------------------------------------------------------------
!  2. Write topology file (objects)
!  ...Create dataspace
      call h5screate_simple_f(obj_rank,obj_dims, dspace_id,ierr)
      call hdf5_w_handle_err(ierr,'geometry_write: h5screate_simple_f')
!
!  ...Create dataset
      call h5dcreate_f(file_id,'/Objects',H5T_NATIVE_INTEGER,dspace_id, &
                       dset_id,ierr,H5P_DEFAULT_F,H5P_DEFAULT_F,H5P_DEFAULT_F)
      call hdf5_w_handle_err(ierr,'geometry_write: h5dcreate_f')
!
!  ...Write dataset
      call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,GEOM_OBJ,obj_dims, ierr)
      call hdf5_w_handle_err(ierr,'geometry_write: h5dwrite_f')
!
!  ...Close dataset
      call h5dclose_f(dset_id, ierr)
      call hdf5_w_handle_err(ierr,'geometry_write: h5dclose_f')
!
!  ...Close dataspace
      call h5sclose_f(dspace_id, ierr)
      call hdf5_w_handle_err(ierr,'geometry_write: h5sclose_f')
!
!  -----------------------------------------------------------------------------
!  ...Close file
      call h5fclose_f(file_id, ierr)
      call hdf5_w_handle_err(ierr,'geometry_write: h5fclose_f')
!
   end subroutine geometry_write



!-----------------------------------------------------------------------------
!> Purpose : initialize data arrays for outputting attribute fields
!!
!> @date Feb 2023
!-----------------------------------------------------------------------------
   subroutine attr_init(Ndim,Nval)
!
      integer, intent(in) :: Ndim,Nval
!
      if (allocated(ATTR_VAL)) deallocate(ATTR_VAL)
      allocate(ATTR_VAL(Ndim,Nval))
!
   end subroutine attr_init



!-----------------------------------------------------------------------------
!> Purpose : deallocate attribute arrays
!!
!> @date Feb 2023
!-----------------------------------------------------------------------------
   subroutine attr_close()
      if (allocated(ATTR_VAL)) deallocate(ATTR_VAL)
   end subroutine attr_close



!-----------------------------------------------------------------------------
!> Purpose : write attribute field data to paraview (xdmf)
!!
!> @date Feb 2023
!-----------------------------------------------------------------------------
   subroutine attr_write(Sname,LSname,Sfile,LSfile,Snick,LSnick)
!
      integer, intent(in)   :: LSname,LSfile,LSnick
      character(len=LSname) :: Sname
      character(len=LSfile) :: Sfile
      character(len=LSnick) :: Snick
!
      integer(HID_T) :: file_id   ! File identifier
      integer(HID_T) :: dset_id   ! Dataset identifier
      integer(HID_T) :: dspace_id ! Dataspace identifier
!
      integer        :: ierr
!
!  ...Dataset dimensions and rank
      integer(HSIZE_T), dimension(2) :: attr_dims
      integer :: attr_rank = 2
!
      attr_dims(1) = size(ATTR_VAL,dim=1)
      attr_dims(2) = size(ATTR_VAL,dim=2)
!
!  -----------------------------------------------------------------------------
!
      if (.not. HDF5_IS_INIT) then
         write(*,*) 'attr_write: HDF5_IS_INIT = .false.; returning...'
         return
      endif
!
!  ...Create file
      call h5fcreate_f(Sfile,H5F_ACC_TRUNC_F, file_id,ierr);
      call hdf5_w_handle_err(ierr,'attr_write: h5fcreate_f')
!
!      write(*,*) 'Sname  = ', Sname
!      write(*,*) 'Sfile  = ', Sfile
!      write(*,*) 'Snick  = ', Snick
!
!  -----------------------------------------------------------------------------
!  Write attribute file (values)
!  ...Create dataspace
      call h5screate_simple_f(attr_rank,attr_dims, dspace_id,ierr)
      call hdf5_w_handle_err(ierr,'attr_write: h5screate_simple_f')
!
!  ...Create dataset
      call h5dcreate_f(file_id,Snick,H5T_NATIVE_DOUBLE,dspace_id, &
                       dset_id,ierr,H5P_DEFAULT_F,H5P_DEFAULT_F,H5P_DEFAULT_F)
      call hdf5_w_handle_err(ierr,'attr_write: h5dcreate_f')
!
!  ...Write dataset
      call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,ATTR_VAL,attr_dims, ierr)
      call hdf5_w_handle_err(ierr,'attr_write: h5dwrite_f')
!
!  ...Close dataset
      call h5dclose_f(dset_id, ierr)
      call hdf5_w_handle_err(ierr,'attr_write: h5dclose_f')
!
!  ...Close dataspace
      call h5sclose_f(dspace_id, ierr)
      call hdf5_w_handle_err(ierr,'attr_write: h5sclose_f')
!
!  -----------------------------------------------------------------------------
!  ...Close file
      call h5fclose_f(file_id, ierr)
      call hdf5_w_handle_err(ierr,'attr_write: h5fclose_f')
!
   end subroutine attr_write
!
end module paraview
