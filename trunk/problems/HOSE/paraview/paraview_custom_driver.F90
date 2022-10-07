!-------------------------------------------------------------------------------------------
!> Purpose : main driver
!-------------------------------------------------------------------------------------------
!
subroutine paraview_custom_driver
!
      use upscale
      use environment      , only : QUIET_MODE
      use paraview         , only : VLEVEL,FILE_VIS
      use mpi_param        , only : RANK,ROOT
      !
      implicit none
      character*8, parameter, dimension(3) :: sAttr = (/'sigma_rr','sigma_rt','sigma_tt'/)
      real*8 :: time
      integer :: scomp
!
      integer,save :: id = -1
      logical,save :: initialized = .FALSE.
!
!-------------------------------------------------------------------------------------------
!
!..load files for visualization upscale
      if (.not. initialized) then
         call load_vis(TETR_VIS,trim(FILE_VIS)//'/tetra_'//trim(VLEVEL),'tetr')
         call load_vis(PRIS_VIS,trim(FILE_VIS)//'/prism_'//trim(VLEVEL),'pris')
         call load_vis(HEXA_VIS,trim(FILE_VIS)//'/hexa_'//trim(VLEVEL),'hexa')
         initialized = .true.
      endif
!
      time=-1.d0 ! set to negative for now
!
!     integer id to append to Fname
      id=id+1
!
      if (RANK .eq. ROOT) then
         call paraview_begin(id,time) ! [OPENS THE XMF FILE, WRITES HEADER]
      endif
!
!     -- GEOMETRY --
      call paraview_geom
!
!     loop over custom components (that we care about)
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
endsubroutine paraview_custom_driver
