!-------------------------------------------------------------------------------------------
!> Purpose : main driver
!-------------------------------------------------------------------------------------------
!
subroutine paraview_driver(IParAttr)
!
   use upscale
   use physics
   use data_structure3D , only : NRCOMS
   use environment,        only: QUIET_MODE
   use paraview,           only: PARAVIEW_DUMP_ATTR,FILE_VIS,     &
                                 VLEVEL,PARAVIEW_DUMP_GEOM,       &
                                 IPARATTR_VTU,SECOND_ORDER_VIS,   &
                                 VIS_VTU,ELEM_TYPES,PARAVIEW_DOMAIN
   use mpi_param,          only: RANK,ROOT
!
   implicit none
!
   integer, intent(in) :: IParAttr(NR_PHYSA)
!
   real(8) :: time
   integer :: idx,iattr,iload,icomp
!
   integer, save :: id = -1
   logical, save :: initialized = .false.
!
!-------------------------------------------------------------------------------------------
!
   PARAVIEW_DUMP_GEOM = .true.
!
!..load files for visualization upscale
   if (.not. initialized) then
      if (SECOND_ORDER_VIS) then
         call load_vis(TETR_VIS,trim(FILE_VIS)//'/tetra10',TETR)
         call load_vis(PRIS_VIS,trim(FILE_VIS)//'/prism18',PRIS)
         call load_vis(HEXA_VIS,trim(FILE_VIS)//'/hexa27' ,BRIC)
      else
         call load_vis(TETR_VIS,trim(FILE_VIS)//'/tetra_'//trim(VLEVEL),TETR)
         call load_vis(PRIS_VIS,trim(FILE_VIS)//'/prism_'//trim(VLEVEL),PRIS)
         call load_vis(HEXA_VIS,trim(FILE_VIS)//'/hexa_'//trim(VLEVEL),BRIC)
      endif
      initialized = .true.
   endif
!
   time=-1.d0
!
!..integer id to append to Fname
   id=id+1
!
   if (RANK .eq. ROOT) then
      call paraview_begin(id,time) ! [OPENS THE XMF FILE, WRITES HEADER]
   endif
!
!  -- GEOMETRY --
!..allocation only if using VTU format
   if (VIS_VTU) then
      allocate(IPARATTR_VTU(NR_PHYSA))
      IPARATTR_VTU = iParAttr(1:NR_PHYSA)
   endif
!
   call paraview_geom
!
!  -- PHYSICAL ATTRIBUTES --
   if (.not. PARAVIEW_DUMP_ATTR) goto 80
!
   if (.not. QUIET_MODE .and. RANK .eq. ROOT) then
      if (PARAVIEW_DOMAIN.eq.0) then
         write(*,*)'Dumping to Paraview...'
      else
         write(*,*)'Dumping to Paraview from domain ',PARAVIEW_DOMAIN,'...'
      endif
      write(*,*)'--------------------------------------'
      write(*,*)'ATTRIBUTE | DISC. SPACE | COMP. | LOAD'
   endif
!
!..loop over rhs's
   do iload=1,NRCOMS
!
!  ...loop over attributes
      do iattr=1,NR_PHYSA
!
!     ...loop over components
         do icomp=1,NR_COMP(iattr)
!
!        ...encode iload, iattr, icomp into a single attribute's index
            idx = iload*100 + iattr*10 + icomp*1
!
            select case(D_TYPE(iattr))
!
!              -- H1 --
               case(CONTIN) ; call paraview_attr_scalar(id,idx)
!
!              -- H(curl) --
               case(TANGEN) ; call paraview_attr_vector(id,idx)
!
!              -- H(div) --
               case(NORMAL) ; call paraview_attr_vector(id,idx)
!
!              -- L2 --
               case(DISCON) ; call paraview_attr_scalar(id,idx)
!
            end select
!
            if (.not. QUIET_MODE .and. RANK .eq. ROOT) then
               write(*,8000) PHYSA(iattr),S_DType(D_TYPE(iattr)),icomp,iload
         8000  format(1x,a5,7x,a6,12x,i1,5x,i2)
            endif
!
!     ...end loop over components
         enddo
!  ...end loop over attributes
      enddo
!..loop over rhs's
   enddo
!
   if (.not. QUIET_MODE .and. RANK .eq. ROOT) then
      write(*,*)'--------------------------------------'
      write(*,*)''
   endif
!
   80 continue
!..deallocation only if using VTU output format
   if (VIS_VTU) then
      deallocate(IPARATTR_VTU)
   endif
!
   if (RANK .eq. ROOT) then
      call paraview_end ! [CLOSES THE XMF FILE, WRITES FOOTER]
   endif
!
end subroutine paraview_driver
!
!
!-------------------------------------------------------------------------------------------
!> Purpose : open file, print header
!!
!> @param[in] Id    - integer id to be appended to Fname (postfix)
!> @param[in] Time  - currently non supported, set to -1.d0
!!
!> @date Mar 2023
!-------------------------------------------------------------------------------------------
subroutine paraview_begin(Id,Time)
!
   use paraview    , only : PARAVIEW_IO,PARAVIEW_DIR,VIS_VTU,paraview_init
   use environment , only : PREFIX
!
   implicit none
!
   integer, intent(in) :: Id
   real(8), intent(in) :: Time
!
   character(len=5) :: postfix
!
!-------------------------------------------------------------------------------------------
!
!..XMDF format
   if (.not. VIS_VTU) then
      call paraview_init
!
!  ...convert integer to string
      write(postfix,"(I5.5)") Id
!
!  ...open .xmf file
      open(unit=PARAVIEW_IO , file=trim(PARAVIEW_DIR)//trim(PREFIX)//"_"//trim(postfix)//'.xmf')
!
      write(PARAVIEW_IO, 1001)
      write(PARAVIEW_IO, 1002)
      write(PARAVIEW_IO, 1003)
!
!  ...non negative time is provided
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
      open(unit=PARAVIEW_IO , file=trim(PARAVIEW_DIR)//trim(PREFIX)//"_"//trim(postfix)//'.vtu', status  = 'replace',        &
            form    = 'unformatted',    &
            access  = 'stream')
   endif
!
end subroutine paraview_begin
!
!
!-------------------------------------------------------------------------------------------
!> Purpose : print footer, close file
!!
!> @date Mar 2023
!-------------------------------------------------------------------------------------------
subroutine paraview_end
!
   use paraview , only : PARAVIEW_IO, paraview_finalize,VIS_VTU
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
   call paraview_finalize
!
end subroutine paraview_end
