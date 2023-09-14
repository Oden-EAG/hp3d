!-----------------------------------------------------------------------------
!> @brief   initialize paraview vis files and auxiliary arrays
!> @date    Sep 2023
!-----------------------------------------------------------------------------
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
!> @param[in] Time  - currently only supported with XDMF output format
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
   character(len=5) :: postfix, timestamp
   integer          :: i
!..time_tags saves all time values until the current step.
   real(8),save     :: time_tags(PARAVIEW_TIME_STEPS_NUM)
!
!-------------------------------------------------------------------------------------------
!
   time_tags(Id+1) = Time
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
      1004 format("    <Time Value='", f14.10, "' />")
!
!..VTU format
   else
!  ...opening the VTU file in binary format.
      write(postfix,"(I5.5)") Id
      open(unit=PARAVIEW_IO,                                                  &
           file=trim(PARAVIEW_DIR)//trim(PREFIX)//"_"//trim(postfix)//'.vtu', &
           status = 'replace', form = 'unformatted', access = 'stream')
      
      open(unit = PVD_IO,file=trim(PARAVIEW_DIR)//trim(PREFIX)//'.pvd',status = 'replace', form = 'unformatted', access = 'stream')
      write(PVD_IO) '<VTKFile type="Collection" version="1.0" byte_order="LittleEndian" header_type="UInt64">' //char(10)
      write(PVD_IO) '  <Collection>'//char(10)
!..For every time step, pvd file is written from scratch.
      do i = 0,Id
         if(time_tags(i+1) .gt. -1.d0) then
            write(timestamp,"(f5.3)") time_tags(i+1)
         else
            write(timestamp,"(I5.5)") i
         endif
         write(postfix,"(I5.5)") i
         write(PVD_IO) '     <DataSet timestep=' // '"'//trim(timestamp)//'"'//'  part="0"   file="'//trim(PREFIX)//"_"//trim(postfix)//'.vtu" />'//char(10)
      enddo
!..closing pvd file
      write(PVD_IO) '  </Collection>'//char(10)
      write(PVD_IO) '</VTKFile>'//char(10)
      close(PVD_IO)
      write(postfix,"(I5.5)") Id

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

