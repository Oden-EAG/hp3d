!
#include "typedefs.h"
!
!-------------------------------------------------------------------------------------------
!> Purpose : paraview driver
!-------------------------------------------------------------------------------------------
subroutine my_paraview_driver(IParAttr)
!
   use upscale
   use physics,            only: D_TYPE,NR_COMP,NR_PHYSA,PHYSA
   use data_structure3D
   use environment,        only: PREFIX,QUIET_MODE
   use paraview,           only: PARAVIEW_DUMP_ATTR,FILE_VIS, &
                                 PARAVIEW_DOMAIN,VLEVEL, &
                                 PARAVIEW_DUMP_GEOM
   use mpi_param,          only: RANK,ROOT
   use par_mesh,           only: DISTRIBUTED,HOST_MESH
!
   implicit none
!
   integer, intent(in) :: IParAttr(NR_PHYSA)
!
   real(8) :: time
   integer :: idx,iphys,iload,icomp
!
   integer, save :: id = -1
   logical, save :: initialized = .false.
!
   character(len=2) :: vis_level
!
!-------------------------------------------------------------------------------------------
!
!..Set true to write out geometry file on every call to this routine
!  Set false to write out geometry file only on first call to this routine
   PARAVIEW_DUMP_GEOM = .true.
!
!..load files for visualization upscale
   if (SECOND_ORDER_VIS) then
      if (.not. initialized) then
         call load_vis(TETR_VIS,trim(FILE_VIS)//'/tetra10',TETR)
         call load_vis(PRIS_VIS,trim(FILE_VIS)//'/prism18',PRIS)
         call load_vis(HEXA_VIS,trim(FILE_VIS)//'/hexa27' ,BRIC)
         initialized = .true.
      endif
   else
      if (.not. initialized) then
         call load_vis(TETR_VIS,trim(FILE_VIS)//'/tetra_'//trim(VLEVEL),TETR)
         call load_vis(PRIS_VIS,trim(FILE_VIS)//'/prism_'//trim(VLEVEL),PRIS)
         call load_vis(HEXA_VIS,trim(FILE_VIS)//'/hexa_'//trim(VLEVEL),BRIC)
         initialized = .true.
      endif
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
   call paraview_geometry
!
!  -- PHYSICAL ATTRIBUTES --
   if (.not. PARAVIEW_DUMP_ATTR) goto 90
!
!..loop over rhs's
   do iload=1,1 !NR_RHS
!
!  ...loop over physics variables
      do iphys=1,NR_PHYSA

         if (IParAttr(iphys) .eq. 0) cycle
!
!     ...loop over components
         do icomp=1,NR_COMP(iphys)
            if (IParAttr(iphys) .ge. icomp) then
!
!           ...encode iload, iphys, icomp into a single attribute's index
               idx = iload*100 + iphys*10 + icomp*1
!
               select case(D_TYPE(iphys))
!
!                 -- H1 --
                  case(CONTIN)
                     call paraview_attr_scalar(id,idx)
!
!                 -- H(curl) --
                  case(TANGEN)
                     call paraview_attr_vector(id,idx)
!
!                 -- H(div) --
                  case(NORMAL)
                     call paraview_attr_vector(id,idx)
!
!                 -- L2 --
                  case(DISCON)
                     call paraview_attr_scalar(id,idx)
!
               end select
!
               if (.not. QUIET_MODE .and. RANK .eq. ROOT) then
                  write(*,8000) PHYSA(iphys),S_dtype(D_TYPE(iphys)),icomp,iload
 8000             format(/,1x,a5,7x,a6,12x,i1,5x,i2)
               endif
!
            endif
!     ...end loop over components of physics variable
         enddo
!  ...end loop over physics variables
      enddo
!..end loop over rhs's
   enddo
!
  90 continue
!
   if (RANK .eq. ROOT) then
      call paraview_end ! [CLOSES THE XMF FILE, WRITES FOOTER]
   endif
!
end subroutine my_paraview_driver
