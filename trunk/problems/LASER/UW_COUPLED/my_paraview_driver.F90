!-------------------------------------------------------------------------------------------
!> @brief   Custom interface for paraview export
!> @date    Sep 2023
!-------------------------------------------------------------------------------------------
subroutine my_paraview_driver
!
   use upscale
   use physics
   use environment, only: QUIET_MODE
   use commonParam, only: TIMESTEP
   use laserParam,  only: DELTA_T
   use mpi_param,   only: RANK,ROOT
   use parameters,  only: NRRHS
   use paraview
!
   implicit none
!
   integer :: idx,iattr,iload,icomp,jcomp
!
   integer, save :: id = -1
   real(8), save :: time = 0.d0
!
!-------------------------------------------------------------------------------------------
!
!..check compatibility of paraview input flags, and
!  initialize visualization files and auxiliary arrays
   call paraview_initialize
!
!.."time" value is only written to file if "time" is non-negative
!  (currently, only supported with XDMF output)
!  note: DELTA_T is non-dimensional time increment;
!        DELTA_T * TIME_0 is time increment in seconds.
   time = TIMESTEP * DELTA_T
!
!..integer id to append to Fname
   id=id+1
!
!..TODO no need for extra routines for paraview_begin/end
!  (simply modify prefix) -> does not work for geom b/c of id counter
!  WRITE GEOMETRY WITH VIS 0 TO FILE (VISUALIZE ADAPTIVE REFS)
!   if (RANK .eq. ROOT) then
!      call paraview_begin_vis0(id,time)
!   endif
!   vis_level = VLEVEL; VLEVEL = '0'
!   call load_vis(TETR_VIS,trim(FILE_VIS)//'/tetra_'//trim(VLEVEL),TETR)
!   call load_vis(PRIS_VIS,trim(FILE_VIS)//'/prism_'//trim(VLEVEL),PRIS)
!   call load_vis(HEXA_VIS,trim(FILE_VIS)//'/hexa_'//trim(VLEVEL),BRIC)
!   call paraview_geom_vis0
!   VLEVEL = vis_level
!   call load_vis(TETR_VIS,trim(FILE_VIS)//'/tetra_'//trim(VLEVEL),TETR)
!   call load_vis(PRIS_VIS,trim(FILE_VIS)//'/prism_'//trim(VLEVEL),PRIS)
!   call load_vis(HEXA_VIS,trim(FILE_VIS)//'/hexa_'//trim(VLEVEL),BRIC)
!   if (RANK .eq. ROOT) then
!      call paraview_end_vis0
!   endif
!..END ADAPTIVE REFS GEOMETRY WRITING
!
   if (RANK .eq. ROOT) then
      call paraview_begin(id,time) ! [OPENS THE XMF FILE, WRITES HEADER]
   endif
!
!  -- GEOMETRY --
   call paraview_geometry
!
!  -- PHYSICAL ATTRIBUTES --
   if (.not. PARAVIEW_DUMP_ATTR) goto 80
!
   if (.not.QUIET_MODE .and. RANK.eq.ROOT) then
      if (PARAVIEW_DOMAIN.eq.0) then
         write(*,*)'Dumping to Paraview...'
      else
         write(*,*)'Dumping to Paraview from domain ',PARAVIEW_DOMAIN,'...'
      endif
      write(*,*)'--------------------------------------'
      write(*,*)'ATTRIBUTE | DISC. SPACE | COMP. | LOAD'
   endif
!
!..loop over multiple loads
   do iload=1,NRRHS
!
!  ...skip selected loads
      if (.not. PARAVIEW_LOAD(iload)) cycle
!
      jcomp = 0
!
!  ...loop over attributes
      do iattr=1,NR_PHYSA
!
!     ...skip selected physics variables
         if (.not. PARAVIEW_ATTR(iattr)) then
            jcomp = jcomp + NR_COMP(iattr)
            cycle
         endif
!
!     ...loop over components
         do icomp=1,NR_COMP(iattr)
            jcomp = jcomp+1
!
!        ...skip selected components
#if HP3D_COMPLEX
            if (.not. PARAVIEW_COMP_REAL(jcomp) .and. &
                .not. PARAVIEW_COMP_IMAG(jcomp)) cycle
#else
            if (.not. PARAVIEW_COMP_REAL(jcomp)) cycle
#endif
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
!..end loop over multiple loads
   enddo
!
   if (.not.QUIET_MODE .and. RANK.eq.ROOT) then
      write(*,*)'--------------------------------------'
      write(*,*)''
   endif
!
   80 continue
!
   if (RANK .eq. ROOT) then
      call paraview_end ! [CLOSES THE XMF FILE, WRITES FOOTER]
   endif
!
end subroutine my_paraview_driver



!-------------------------------------------------------------------------------------------
! TESTING
!-------------------------------------------------------------------------------------------
subroutine paraview_geom_vis0
!
   use environment , only : PREFIX
   use paraview    , only : PARAVIEW_IO,PARAVIEW_DUMP_GEOM,PARAVIEW_DIR
   use mpi_param   , only : RANK,ROOT
!
   implicit none
   integer,         save :: id=-1
   character(len=5),save :: postfix
   integer,         save :: ice, icn, icp
!
!-------------------------------------------------------------------------------------------
!
!...h5 file is produced on 1st visit, OR as required by the user
   if (PARAVIEW_DUMP_GEOM .or. (id == -1)) then
!
!  ...increment visitation flag
      id=id+1
!
!  ...convert integer to string
      write(postfix,"(I5.5)") id
!
!  ...produce .h5 file
      call geom2vtk("Geometry",trim(PARAVIEW_DIR)//trim(PREFIX)//"geomV0_"//trim(postfix)//".h5", ice,icn,icp)
!
   endif
!
   if (RANK .ne. ROOT) goto 90
!
!..write to .xmf file
   write(PARAVIEW_IO,1012) ice
   write(PARAVIEW_IO,1013) (ice+icn)
   write(PARAVIEW_IO,1014) trim(PREFIX), trim(postfix)
   write(PARAVIEW_IO,1015)
   write(PARAVIEW_IO,1016)
   write(PARAVIEW_IO,1017)
   write(PARAVIEW_IO,1018) icp
   write(PARAVIEW_IO,1019) trim(PREFIX), trim(postfix)
   write(PARAVIEW_IO,1020)
   write(PARAVIEW_IO,1021)
!
 1012 format("      <Topology TopologyType='Mixed' NumberOfElements='",i8,"'>")
 1013 format("        <DataItem Dimensions='",i12,"' NumberType='Int' Precision='4' Format='HDF'>")
 1014 format("        ",a,"geomV0_",a,".h5:/Objects")
 1015 format("        </DataItem>")
 1016 format("      </Topology>")
 1017 format("      <Geometry GeometryType='XYZ'>")
 1018 format("        <DataItem Dimensions='",i10, " 3' NumberType='Float' Precision='4' Format='HDF'>")
 1019 format("        ",a,"geomV0_",a,".h5:/Coords")
 1020 format("        </DataItem>")
 1021 format("      </Geometry>")
!
   90 continue
!
end subroutine paraview_geom_vis0
!
!
!-------------------------------------------------------------------------------------------
subroutine paraview_begin_vis0(Id,Time)
!
   use paraview    , only : PARAVIEW_IO,PARAVIEW_DIR,paraview_data_init
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
   call paraview_data_init
!
!..convert integer to string
   write(postfix,"(I5.5)") Id
!
!..open .xmf file
   open(unit=PARAVIEW_IO , file=trim(PARAVIEW_DIR)//trim(PREFIX)//"V0_"//trim(postfix)//'.xmf')
!
   write(PARAVIEW_IO, 1001)
   write(PARAVIEW_IO, 1002)
   write(PARAVIEW_IO, 1003)
!
!..non negative time is provided
   if (Time >= 0.d0) then
      write(PARAVIEW_IO, 1004) Time
   end if
!
!..HEADER of .xmf file
 1001 format("<Xdmf xmlns:xi='http://www.w3.org/2003/XInclude' Version='2.1'>")
 1002 format("  <Domain>")
 1003 format("    <Grid Name='Geometry' GridType='Uniform'>")
 1004 format("    <Time Value='", f14.10, "' />")
!
end subroutine paraview_begin_vis0
!
!
!-------------------------------------------------------------------------------------------
subroutine paraview_end_vis0
!
   use paraview , only : PARAVIEW_IO, paraview_data_finalize
!
   implicit none
!
!..FOOTER of .xmf file
   write(PARAVIEW_IO, 1004)
   write(PARAVIEW_IO, 1005)
   write(PARAVIEW_IO, 1006)
 1004 format("    </Grid>")
 1005 format("  </Domain>")
 1006 format("</Xdmf>")
!
   close(PARAVIEW_IO)
!
   call paraview_data_finalize
!
end subroutine paraview_end_vis0
