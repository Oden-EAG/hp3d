!-------------------------------------------------------------------------------------------
!> @brief   custom paraview driver
!> @date    Sep 2023
!-------------------------------------------------------------------------------------------
subroutine paraview_custom_driver
!
      use upscale
      use environment, only : QUIET_MODE
      use mpi_param  , only : RANK,ROOT
      !
      implicit none
      character*8, parameter :: sAttr(3) = (/'sigma_rr','sigma_rt','sigma_tt'/)
      real(8) :: time
      integer :: scomp
!
      integer,save :: id = -1
!
!-------------------------------------------------------------------------------------------
!
!..check compatibility of paraview input flags, and
!  initialize visualization files and auxiliary arrays
   call paraview_initialize
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
   call paraview_geometry
!
!  loop over custom components (that we care about)
   do scomp=1,3
      call paraview_attr_custom(id,scomp)
!
      if ((.not. QUIET_MODE) .and. (RANK.eq.ROOT)) then
                write(*,8000) sAttr(scomp)
       8000     format(1x,a8)
      endif
!
   enddo
!
   if (RANK .eq. ROOT) then
      call paraview_end ! [CLOSES THE XMF FILE, WRITES FOOTER]
   endif
!
end subroutine paraview_custom_driver
