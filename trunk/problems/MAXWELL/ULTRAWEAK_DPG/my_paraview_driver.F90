!-------------------------------------------------------------------------------------------
!> @brief   Interface for paraview export
!> @date    Sep 2023
!-------------------------------------------------------------------------------------------
subroutine my_paraview_driver(IParAttr)
!
   use upscale
   use physics
   use data_structure3D, only: NRCOMS
   use environment,      only: QUIET_MODE
   use mpi_param,        only: RANK,ROOT
   use paraview
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
!..check compatibility of paraview input flags
   call paraview_check
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
!.."time" value is only written to file if "time" is non-negative
!  (currently, only supported with XDMF output)
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
   if (VIS_VTU) then
      allocate(IPARATTR_VTU(NR_PHYSA))
      IPARATTR_VTU = iParAttr(1:NR_PHYSA)
   endif
!
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
   if (.not.QUIET_MODE .and. RANK.eq.ROOT) then
      write(*,*)'--------------------------------------'
      write(*,*)''
   endif
!
   80 continue
!
   if (VIS_VTU) then
      deallocate(IPARATTR_VTU)
   endif
!
   if (RANK .eq. ROOT) then
      call paraview_end ! [CLOSES THE XMF FILE, WRITES FOOTER]
   endif
!
end subroutine my_paraview_driver
