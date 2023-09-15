!-------------------------------------------------------------------------------------------
!> @brief   initialize paraview vis files and auxiliary arrays
!> @date    Sep 2023
!-------------------------------------------------------------------------------------------
   subroutine paraview_initialize()
      use control         , only: EXGEOM
      use data_structure3D, only: NRCOMS
      use node_types      , only: TETR,PRIS,BRIC
      use physics         , only: NR_PHYSA,NRINDEX
      use upscale         , only: TETR_VIS,PRIS_VIS,HEXA_VIS,load_vis
      use paraview
!
!  ...check compatibility of paraview input flags
      call paraview_check
!
!  ...paraview ex./iso. geometry flag
      PARAVIEW_GEOM = EXGEOM
!
!  ...load visualization files
      if (SECOND_ORDER_VIS) then
         call load_vis(TETR_VIS,trim(FILE_VIS)//'/tetra10',TETR)
         call load_vis(PRIS_VIS,trim(FILE_VIS)//'/prism18',PRIS)
         call load_vis(HEXA_VIS,trim(FILE_VIS)//'/hexa27' ,BRIC)
      else
         call load_vis(TETR_VIS,trim(FILE_VIS)//'/tetra_'//trim(VLEVEL),TETR)
         call load_vis(PRIS_VIS,trim(FILE_VIS)//'/prism_'//trim(VLEVEL),PRIS)
         call load_vis(HEXA_VIS,trim(FILE_VIS)//'/hexa_'//trim(VLEVEL),BRIC)
      endif
!
!  ...check if already initialized
      if (PARAVIEW_IS_INIT) goto 90
!
!  ...allocate array for selectively exporting variables
      allocate(PARAVIEW_LOAD(NRCOMS))      ; PARAVIEW_LOAD = .true.
      allocate(PARAVIEW_ATTR(NR_PHYSA))    ; PARAVIEW_ATTR = .true.
      allocate(PARAVIEW_COMP_REAL(NRINDEX)); PARAVIEW_COMP_REAL = .true.
      allocate(PARAVIEW_COMP_IMAG(NRINDEX)); PARAVIEW_COMP_IMAG = .true.
!
      PARAVIEW_IS_INIT = .true.
   90 continue
!
   end subroutine paraview_initialize
!
!
!-------------------------------------------------------------------------------------------
!> @brief     open file, print header
!!
!> @param[in] Id    - integer id to be appended to Fname (postfix)
!> @param[in] Time  - time value for this timestep
!!
!> @date      Sep 2023
!-------------------------------------------------------------------------------------------
subroutine paraview_begin(Id,Time)
!
   use environment, only : PREFIX
   use paraview
!
   implicit none
!
   integer, intent(in) :: Id
   real(8), intent(in) :: Time
!
   character(len=5)     :: postfix
   character(len=8)     :: timestamp
   integer              :: i
!
!..time_tags array saves all previously written time values
   real(8), allocatable, save :: time_tags(:)
!
   integer              :: max_new
   real(8), allocatable :: temp(:)
!
!-------------------------------------------------------------------------------------------
!
   if (.not. allocated(time_tags)) then
      allocate(time_tags(PARAVIEW_MAX_TIMESTEPS))
   endif
!
   if (Id .lt. 0) then
      write(*,*) 'paraview_begin: Id = ',Id
      stop
   endif
!
!..check if time_tags array size needs to be increased
   if (Id .ge. PARAVIEW_MAX_TIMESTEPS) then
!  ...determine size of new time_tags array
      max_new = 2*PARAVIEW_MAX_TIMESTEPS
!  ...allocate new time_tags array
      allocate(temp(max_new))
!  ...copy data from old time_tags array into the new one
      temp(1:PARAVIEW_MAX_TIMESTEPS) = time_tags(1:PARAVIEW_MAX_TIMESTEPS)
!  ...move time_tags pointer to the new array
!     and deallocate the old array
      call move_alloc(temp, time_tags)
      PARAVIEW_MAX_TIMESTEPS = max_new
   endif
!
   time_tags(Id+1) = Time
!
   call paraview_data_init
!
!..XDMF format
   if (.not. VIS_VTU) then
!
!  ...convert integer to string
      write(postfix,"(I5.5)") Id
!
!  ...open .xmf file
      open(unit=PARAVIEW_IO, &
           file=trim(PARAVIEW_DIR)//trim(PREFIX)//"_"//trim(postfix)//'.xmf')
!
      write(PARAVIEW_IO, 1001)
      write(PARAVIEW_IO, 1002)
      write(PARAVIEW_IO, 1003)
!
!  ...non-negative time is provided
      if (Time >= 0.d0) then
         write(PARAVIEW_IO, 1004) Time
      end if
!
!  ...HEADER of .xmf file
      1001 format("<Xdmf xmlns:xi='http://www.w3.org/2003/XInclude' Version='2.1'>")
      1002 format("  <Domain>")
      1003 format("    <Grid Name='Geometry' GridType='Uniform'>")
      1004 format("    <Time Value='", f8.4, "' />")
!
!..VTU format
   else
!  ...opening VTU file in binary format
      write(postfix,"(I5.5)") Id
      open(unit=PARAVIEW_IO,                                                  &
           file=trim(PARAVIEW_DIR)//trim(PREFIX)//"_"//trim(postfix)//'.vtu', &
           status = 'replace', form = 'unformatted', access = 'stream')
!
!  ...opening PVD file
      open(unit=PVD_IO,                                                       &
           file=trim(PARAVIEW_DIR)//trim(PREFIX)//'.pvd',                     &
           status = 'replace', form = 'unformatted', access = 'stream')
!
      write(PVD_IO) '<VTKFile type="Collection" version="1.0" '//             &
                             'byte_order="LittleEndian" '//                   &
                             'header_type="UInt64">'//char(10)
      write(PVD_IO) '  <Collection>'//char(10)
!
!  ...for every time step, PVD file is written from scratch
      do i = 0,Id
         if(time_tags(i+1) >= 0.d0) then
            write(timestamp,"(f8.4)") time_tags(i+1)
         else
            write(timestamp,"(I5.5)") i
         endif
         write(postfix,"(I5.5)") i
         write(PVD_IO) '    <DataSet timestep='//'"'//trim(timestamp)//'"'    &
                                 //' part="0" file="'//trim(PREFIX)//"_"      &
                                 //trim(postfix)//'.vtu" />'//char(10)
      enddo
!
!  ...closing PVD file
      write(PVD_IO) '  </Collection>'//char(10)
      write(PVD_IO) '</VTKFile>'//char(10)
      close(PVD_IO)
   endif
!
end subroutine paraview_begin
!
!
!-------------------------------------------------------------------------------------------
!> @brief   print footer, close file
!> @date    Sep 2023
!-------------------------------------------------------------------------------------------
subroutine paraview_end
!
   use paraview
!
   implicit none
!
!..XDMF format
   if (.not. VIS_VTU) then
!  ...FOOTER of .xmf file
      write(PARAVIEW_IO, 1004)
      write(PARAVIEW_IO, 1005)
      write(PARAVIEW_IO, 1006)
      1004 format("    </Grid>")
      1005 format("  </Domain>")
      1006 format("</Xdmf>")
!
!..VTU format
   else
!  ...closing strings for the VTU file
      write(PARAVIEW_IO) char(10) // '' // '</AppendedData>' // char(10)
      write(PARAVIEW_IO)       '' // '</VTKFile>'      // char(10)
   endif
!
   close(PARAVIEW_IO)
!
   call paraview_data_finalize
!
end subroutine paraview_end
