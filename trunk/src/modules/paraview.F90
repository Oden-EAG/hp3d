!-------------------------------------------------------------------------------------------
!> @brief   paraview module for exporting mesh and solutions for visualization
!> @date    Sep 2023
!-------------------------------------------------------------------------------------------
module paraview
!
   use hdf5_wrapper
   use HDF5
!
   implicit none
!
   save
!
!..FLAGS:
!  - VLEVEL
!     linear upscaling of geometry output via "VLEVEL" uniform refinements
!     only for linear geometry output (i.e., when SECOND_ORDER_VIS = .false.)
!  - PARAVIEW_DUMP_GEOM
!     if enabled : paraview_geometry writes geometry on every call
!     if disabled: paraview_geometry writes geometry on first call only (supported with XDMF)
!  - PARAVIEW_DUMP_ATTR
!     if enabled : paraview_driver writes all physical attributes (solutions)
!     if disabled: paraview_driver does not write any solution fields
   character(len=128) :: FILE_VIS = '../../files/vis'
   character(len=2  ) :: VLEVEL = '2'
   character(len=128) :: PARAVIEW_DIR = 'vtk/'
   integer,parameter  :: PARAVIEW_IO = 22
   logical            :: PARAVIEW_DUMP_GEOM = .true.
   logical            :: PARAVIEW_DUMP_ATTR = .true.
!
!..Flag for enabling second-order geometry output (both XDMF and VTU)
!  NOTE: Only VTU supports higher-order geometry for hybrid meshes
   logical            :: SECOND_ORDER_VIS   = .false.
!
!..Flag for switching between different formats for output meshes and solution fields
!  ...VIS_VTU = .false. : uses XDMF format (default)
!  ...VIS_VTU = .true.  : uses VTU format
   logical            :: VIS_VTU = .false.
!
!..PVD_IO: File identifier for writting pvd file.
   integer,parameter  :: PVD_IO  = 23
!
!..PARAVIEW_TIME: time tag set by the user
   real(8)            :: PARAVIEW_TIME = -1.d0
!
!..maximum number of time steps anticipated
!  (used to pre-allocate array for saving time tags)
   integer            :: PARAVIEW_MAX_TIMESTEPS = 10000
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
   integer, allocatable :: VTU_ELEM_TYPES(:)
!
!..Advanced control over exporting variables
!     PARAVIEW_LOAD enables/disables writing specific solution vectors
!     PARAVIEW_ATTR enables/disables writing specific physics variables
!     PARAVIEW_COMP_REAL enables/disables writing real part of specific components
!     PARAVIEW_COMP_IMAG enables/disables writing imag part of specific components
   logical, allocatable :: PARAVIEW_LOAD(:)      ! 1...NRRHS
   logical, allocatable :: PARAVIEW_ATTR(:)      ! 1...NR_PHYSA
   logical, allocatable :: PARAVIEW_COMP_REAL(:) ! 1...NRINDEX
   logical, allocatable :: PARAVIEW_COMP_IMAG(:) ! 1...NRINDEX
!
   logical              :: PARAVIEW_IS_INIT = .false.
!
   contains



!-----------------------------------------------------------------------------
!> @brief     paraview compatibility checks
!> @date      Sep 2023
!-----------------------------------------------------------------------------
   subroutine paraview_check
      use GMP, only: NRHEXAS,NRPRISM,NRTETRA,NRPYRAM
      integer :: num_types
!
!  ...check for compatibility of writing the mesh only at first call with output format
!     (VTU format requires the mesh to be written each time the solution is exported)
      if (.not.PARAVIEW_DUMP_GEOM .and. VIS_VTU) then
         if (RANK .eq. ROOT) then
            write(*,*) 'paraview_check: PARAVIEW_DUMP_GEOM disabled is not supported for VTU.'
            write(*,*) '                Setting PARAVIEW_DUMP_GEOM=.true.'
            write(*,*) '                (writing mesh each time the solution is exported).'
         endif
         PARAVIEW_DUMP_GEOM = .true.
      endif
!
!  ...check for compatibility of upscaling with output order
!     (VLEVEL upscaling is only supported for linear element output)
      if (SECOND_ORDER_VIS .and. VLEVEL.ne.'0') then
         if (RANK .eq. ROOT) then
            write(*,*) 'paraview_check: VLEVEL upscaling is only supported for linear output.'
            write(*,*) '                Setting VLEVEL to "0" for SECOND_ORDER_VIS.'
         endif
         VLEVEL = '0'
      endif
!
!  ...check for compatibility of mixed mesh with output format
!     (XDMF does not support mixed meshes with higher-order)
      num_types = 0
      if (NRHEXAS .gt. 0) num_types = num_types + 1
      if (NRPRISM .gt. 0) num_types = num_types + 1
      if (NRTETRA .gt. 0) num_types = num_types + 1
      if (NRPYRAM .gt. 0) num_types = num_types + 1
!
      if (SECOND_ORDER_VIS .and. .not.VIS_VTU) then
         if (num_types .gt. 1) then
            if (RANK .eq. ROOT) then
               write(*,*) 'paraview_check: XDMF does not support hybrid mesh with SECOND_ORDER_VIS.'
               write(*,*) '                Changing output format from XDMF to VTU.'
            endif
            VIS_VTU = .true.
         endif
      endif
!
   end subroutine paraview_check



!-----------------------------------------------------------------------------
!> @brief     select solution vectors for paraview export
!> @date      Sep 2023
!-----------------------------------------------------------------------------
   subroutine paraview_select_load(EnableLoad)
      use parameters, only: NRRHS
      logical, intent(in) :: EnableLoad(NRRHS)
      integer :: iload
      if (.not.PARAVIEW_IS_INIT) call paraview_initialize
!
      do iload = 1,NRRHS
         PARAVIEW_LOAD(iload) = EnableLoad(iload)
      enddo
   end subroutine paraview_select_load


!-----------------------------------------------------------------------------
!> @brief     select physics attributes for paraview export
!> @date      Sep 2023
!-----------------------------------------------------------------------------
   subroutine paraview_select_attr(EnableAttr)
      use physics, only: NR_PHYSA
      logical, intent(in) :: EnableAttr(NR_PHYSA)
      integer :: iattr
      if (.not.PARAVIEW_IS_INIT) call paraview_initialize
!
      do iattr = 1,NR_PHYSA
         PARAVIEW_ATTR(iattr) = EnableAttr(iattr)
      enddo
   end subroutine paraview_select_attr


!-----------------------------------------------------------------------------
!> @brief     select components (real part) for paraview export
!> @date      Sep 2023
!-----------------------------------------------------------------------------
   subroutine paraview_select_comp_real(EnableCompReal)
      use physics, only: NR_PHYSA,NR_COMP,NRINDEX
      logical, intent(in) :: EnableCompReal(NRINDEX)
      integer :: iattr,icomp,jcomp
      if (.not.PARAVIEW_IS_INIT) call paraview_initialize
!
      jcomp = 0
      do iattr = 1,NR_PHYSA
         do icomp = 1,NR_COMP(iattr)
            jcomp = jcomp+1
            PARAVIEW_COMP_REAL(jcomp) = EnableCompReal(jcomp)
         enddo
      enddo
   end subroutine paraview_select_comp_real


!-----------------------------------------------------------------------------
!> @brief     select components (imaginary part) for paraview export
!> @date      Sep 2023
!-----------------------------------------------------------------------------
   subroutine paraview_select_comp_imag(EnableCompImag)
      use physics, only: NR_PHYSA,NR_COMP,NRINDEX
      logical, intent(in) :: EnableCompImag(NRINDEX)
      integer :: iattr,icomp,jcomp
      if (.not.PARAVIEW_IS_INIT) call paraview_initialize
!
      jcomp = 0
      do iattr = 1,NR_PHYSA
         do icomp = 1,NR_COMP(iattr)
            jcomp = jcomp+1
            PARAVIEW_COMP_IMAG(jcomp) = EnableCompImag(jcomp)
         enddo
      enddo
   end subroutine paraview_select_comp_imag



!-----------------------------------------------------------------------------
!> Purpose : initialize paraview output (open hdf5 file)
!!
!> @date Mar 2023
!-----------------------------------------------------------------------------
   subroutine paraview_data_init()
!  ...only call hdf5_w_finalize when xdmf output is chosen
      if (.not. VIS_VTU) call hdf5_w_init()
   end subroutine paraview_data_init



!-----------------------------------------------------------------------------
!> Purpose : finalize paraview output (close hdf5; deallocate data arrays)
!!
!> @date Mar 2023
!-----------------------------------------------------------------------------
   subroutine paraview_data_finalize()
!  ...only call hdf5_w_finalize when xdmf output is chosen
      if (.not. VIS_VTU) call hdf5_w_finalize()
      if (allocated(GEOM_OBJ)) deallocate(GEOM_OBJ)
      if (allocated(GEOM_PTS)) deallocate(GEOM_PTS)
      if (allocated(GEOM_PTS)) deallocate(ATTR_VAL)
   end subroutine paraview_data_finalize



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
