!-------------------------------------------------------------------------------------------
!> Purpose : paraview driver
!-------------------------------------------------------------------------------------------
!
subroutine my_paraview_driver(IParAttr)
!
   use upscale
   use physics,            only: DTYPE,NR_COMP,NR_PHYSA,PHYSA
   !use data_structure3D,   only: NRCOMS
   use environment,        only: QUIET_MODE
   use paraview,           only: PARAVIEW_DUMP_ATTR,FILE_VIS, &
                                 VLEVEL,PARAVIEW_DUMP_GEOM, VIS_FORMAT,SECOND_ORDER_VIS
   use mpi_param,          only: RANK,ROOT
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
!-------------------------------------------------------------------------------------------
!
   PARAVIEW_DUMP_GEOM = .true.
!
!..load files for visualization upscale
   if (.not. initialized) then
      if(SECOND_ORDER_VIS .eq. 0) then
         call load_vis(TETR_VIS,trim(FILE_VIS)//'/tetra_'//trim(VLEVEL),'tetr')
         call load_vis(PRIS_VIS,trim(FILE_VIS)//'/prism_'//trim(VLEVEL),'pris')
         call load_vis(HEXA_VIS,trim(FILE_VIS)//'/hexa_'//trim(VLEVEL),'hexa')
      else
         call load_vis(TETR_VIS,trim(FILE_VIS)//'/tetra_10','tetr')
         call load_vis(PRIS_VIS,trim(FILE_VIS)//'/prism_18','pris')
         call load_vis(HEXA_VIS,trim(FILE_VIS)//'/hexa_27','hexa')
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
   call paraview_geom

!
!  -- PHYSICAL ATTRIBUTES --
   if (.not. PARAVIEW_DUMP_ATTR) goto 90
!
!..loop over rhs's
   do iload=1,1 !NRCOMS

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
               select case(DTYPE(iphys))
!
               !  -- H1 --
                  case('contin')
                     call paraview_attr_scalar(id,idx)
!
!                 -- H(curl) --
                  case('tangen')
                     call paraview_attr_vector(id,idx)
! !
!                 -- H(div) --
                  case('normal')
                     call paraview_attr_vector(id,idx)
! !
!                 -- L2 --
                  case('discon')
                     call paraview_attr_scalar(id,idx)
!
               end select
!
               if ((.not. QUIET_MODE) .and. (RANK.eq.ROOT)) then
                  write(*,8000) PHYSA(iphys),DTYPE(iphys),icomp,iload
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
! it's a hack for poly order (only works for VTU)
 call paraview_attr_Order(id,111)
! !
  90 continue
! !
   if (RANK .eq. ROOT) then
      call paraview_end ! [CLOSES THE XMF FILE, WRITES FOOTER]
   endif
!
!
end subroutine my_paraview_driver
