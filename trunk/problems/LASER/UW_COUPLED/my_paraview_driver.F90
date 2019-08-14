!-------------------------------------------------------------------------------------------
!> Purpose : paraview driver
!-------------------------------------------------------------------------------------------
!
subroutine my_paraview_driver(IParAttr)
!
   use upscale
   use physics,            only: DTYPE,NR_COMP,NR_PHYSA,PHYSA
   use data_structure3D,   only: NRCOMS
   use environment,        only: PREFIX,QUIET_MODE
   use paraview,           only: PARAVIEW_DUMP_ATTR,VLEVEL,FILE_VIS,PARAVIEW_DOMAIN
!
   implicit none
!
   integer, intent(in) :: IParAttr(NR_PHYSA)
!
   real*8  :: time
   integer :: idx,iphys,iload,icomp      
!
   integer, save :: id = -1
   logical, save :: initialized = .false.
!
!..timer variables
   integer(kind=8) :: t1,t2,clock_rate,clock_max
!
!-------------------------------------------------------------------------------------------
!
!..load files for visualization upscale
   if (.not. initialized) then
      call load_vis(TETR_VIS,trim(FILE_VIS)//'/tetra_'//trim(VLEVEL),'tetr')
      call load_vis(PRIS_VIS,trim(FILE_VIS)//'/prism_'//trim(VLEVEL),'pris')
      call load_vis(HEXA_VIS,trim(FILE_VIS)//'/hexa_'//trim(VLEVEL),'hexa')
      call vis_initialize
      initialized = .true.
   endif
!
   time=-1.d0
!
!..integer id to append to Fname
   id=id+1
   call paraview_begin(id,time)
!
!  -- GEOMETRY --
   call system_clock( t1, clock_rate, clock_max )
   call paraview_geom
   call system_clock( t2, clock_rate, clock_max )
   if (.not. QUIET_MODE) then
      write(*,1010) real(t2 - t1,8)/real(clock_rate,8)
 1010 format(' Paraview geometry: ',f12.5,'  seconds',/)
   endif
!
!  -- PHYSICAL ATTRIBUTES --
   if (.not. PARAVIEW_DUMP_ATTR) goto 90
!
   if (.not. QUIET_MODE) then
      if (PARAVIEW_DOMAIN.eq.0) then
         write(*,*)'Dumping to Paraview...'
      else 
         write(*,*)'Dumping to Paraview from domain',PARAVIEW_DOMAIN,'...'
      endif
      write(*,*)'--------------------------------------'
      write(*,*)'ATTRIBUTE | DISC. SPACE | COMP. | LOAD'
   endif
!
   call system_clock( t1, clock_rate, clock_max )
!..loop over rhs's
   do iload=1,1 !NRCOMS
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
               select case(DTYPE(iphys))
!
!                 -- H1 --
                  case('contin')
                     call paraview_attr_scalar(id,idx)
!
!                 -- H(curl) --
                  case('tangen')
                     call paraview_attr_vector(id,idx)
!
!                 -- H(div) --
                  case('normal')
                     call paraview_attr_vector(id,idx)
!
!                 -- L2 --
                  case('discon')
                     call paraview_attr_scalar(id,idx)
!
               end select
!
               if (.not.QUIET_MODE) then
                  write(*,8000) PHYSA(iphys),DTYPE(iphys),icomp,iload
 8000             format(1x,a5,7x,a6,12x,i1,5x,i2)
               endif
!
            endif
!     ...end loop over components of physics variable
         enddo
!  ...end loop over physics variables
      enddo
!..end loop over rhs's
   enddo
   call system_clock( t2, clock_rate, clock_max )
!
   if (.not.QUIET_MODE) then
      write(*,*)'--------------------------------------'
      write(*,1020) real(t2 - t1,8)/real(clock_rate,8)
 1020 format(' Paraview fields  : ',f12.5,'  seconds',/)
      write(*,*)''
   endif
!
  90 continue
!
   call paraview_end
!
!
end subroutine my_paraview_driver
