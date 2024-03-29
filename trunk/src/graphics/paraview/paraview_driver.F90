!-------------------------------------------------------------------------------------------
!> @brief   Default interface for paraview export
!> @date    Sep 2023
!-------------------------------------------------------------------------------------------
subroutine paraview_driver
!
   use upscale
   use physics
   use environment, only: QUIET_MODE
   use mpi_param,   only: RANK,ROOT
   use parameters,  only: NRRHS
   use paraview
!
   implicit none
!
   integer :: idx,iattr,iload,icomp,jcomp
!
   integer, save :: id = -1
!
!-------------------------------------------------------------------------------------------
!
!..check compatibility of paraview input flags, and
!  initialize visualization files and auxiliary arrays
   call paraview_initialize
!
!..integer id to append to Fname
   id=id+1
!
   if (RANK .eq. ROOT) then
!  ...OPENS XMF OR VTU/PVD FILE, WRITES HEADER
      call paraview_begin(id,PARAVIEW_TIME)
   endif
!
!  -- GEOMETRY --
   call paraview_geometry
!
!  -- PHYSICAL ATTRIBUTES --
   if (.not. PARAVIEW_DUMP_ATTR) goto 80
!
   if (.not.QUIET_MODE .and. RANK.eq.ROOT) then
      write(*,*)'Dumping ParaView files to ',PARAVIEW_DIR
      if (PARAVIEW_DOMAIN.ne.0) then
         write(*,*)'Dumping ParaView files for domain ',PARAVIEW_DOMAIN
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
end subroutine paraview_driver
