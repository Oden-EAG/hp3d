!-------------------------------------------------------------------------------------------
!> Purpose : main driver
!-------------------------------------------------------------------------------------------
!
subroutine paraview_custom_driver
!
      use upscale
      use physics          , only : DTYPE,NR_COMP,NR_PHYSA,PHYSA
      use data_structure3D , only : NRCOMS
      use environment      , only : PREFIX,QUIET_MODE
      use paraview         , only : PARAVIEW_DUMP_ATTR,VLEVEL,FILE_VIS,PARAVIEW_DOMAIN
!
      implicit none
      character*8, parameter, dimension(4) :: sAttr = (/'sigma_rr','sigma_rt','sigma_tt','displ'/)
      real*8 :: time
      integer :: idx,iattr,iload,icomp,scomp
      integer,save :: id = -1
!
      logical,save :: initialized = .FALSE.
!
!-------------------------------------------------------------------------------------------
!
!     load files for visualization upscale
      if (.NOT. initialized) then
        call load_vis(TETR_VIS,trim(FILE_VIS)//'/tetra_'//trim(VLEVEL),'tetr')
        call load_vis(PRIS_VIS,trim(FILE_VIS)//'/prism_'//trim(VLEVEL),'pris')
        call load_vis(HEXA_VIS,trim(FILE_VIS)//'/hexa_' //trim(VLEVEL),'hexa')
        call vis_initialize
!
        initialized = .TRUE.
      endif
!
      time=-1.d0 ! set to negative for now
!
!     integer id to append to Fname
      id=id+1
      call paraview_custom_begin(id,time)
!
!     -- GEOMETRY --
      call paraview_geom
!
!     -- PHYSICAL ATTRIBUTES --
      if (PARAVIEW_DUMP_ATTR) then
!
IF (.NOT.QUIET_MODE) THEN
        ! write(*,*)'Dumping to Paraview...'
        write(*,*)'-----------------------'
        if (PARAVIEW_DOMAIN.eq.0) then
          write(*,*)'       ATTRIBUTE       '
        elseif (PARAVIEW_DOMAIN.eq.1) then
          write(*,*)'    STEEL ATTRIBUTE    '
        elseif (PARAVIEW_DOMAIN.eq.2) then
          write(*,*)'   RUBBER ATTRIBUTE    '
        endif
ENDIF
!
!       loop over custom components (that we care about)
        do scomp=1,4
          call paraview_attr_custom(id,scomp)
!
IF (.NOT.QUIET_MODE) THEN
          write(*,8000) sAttr(scomp)
 8000     format(1x,a8)
ENDIF
!
        enddo
!
IF (.NOT.QUIET_MODE) THEN
        write(*,*)'-----------------------'
        write(*,*)''
ENDIF
!
      endif
!
      call paraview_custom_end
!
endsubroutine paraview_custom_driver
!
!
!
!-------------------------------------------------------------------------------------------
!> Purpose : open file, print header
!!
!> @param[in] Id    - integer id to be appended to Fname (postfix)
!> @param[in] Time  - currently non supported, set to -1.d0
!> @param[in] Fname - file name (i.e., physical attribute name)
!!
!> @date Nov 14
!-------------------------------------------------------------------------------------------
!
subroutine paraview_custom_begin(Id, Time)
!
      use paraview    , only : PARAVIEW_IO,PARAVIEW_DIR
      use environment , only : PREFIX
!
      implicit none
      integer         ,intent(in) :: Id
      real*8          ,intent(in) :: Time
!
      character(len=5) :: postfix
!-------------------------------------------------------------------------------------------
!
!     convert integer to string
      write(postfix,"(I5.5)") Id
!
!     open .xmf file
      open(unit=PARAVIEW_IO , file = trim(PARAVIEW_DIR) // trim(PREFIX) &
                                                        // trim(postfix) // 'custom.xmf')
!
      write(PARAVIEW_IO, 1001)
      write(PARAVIEW_IO, 1002)
      write(PARAVIEW_IO, 1003)
!
!     non negative time is provided
      if (Time >= 0.d0) then
         write(PARAVIEW_IO, 1004) Time
      end if
!
!     HEADER of .xmf file
 1001 format("<Xdmf xmlns:xi='http://www.w3.org/2003/XInclude' Version='2.1'>")
 1002 format("  <Domain>")
 1003 format("    <Grid Name='Geometry' GridType='Uniform'>")
 1004 format("    <Time Value='", f14.10, "' />")
!
!
endsubroutine paraview_custom_begin
!
!
!
!-------------------------------------------------------------------------------------------
!> Purpose : print footer, close file
!-------------------------------------------------------------------------------------------
!
subroutine paraview_custom_end
!
      use paraview , only : PARAVIEW_IO
!
      implicit none
!-------------------------------------------------------------------------------------------
!
!     FOOTER of .xmf file
      write(PARAVIEW_IO, 1004)
      write(PARAVIEW_IO, 1005)
      write(PARAVIEW_IO, 1006)
 1004 format("    </Grid>")
 1005 format("  </Domain>")
 1006 format("</Xdmf>")
!
      close(PARAVIEW_IO)
!
!
endsubroutine paraview_custom_end
!
!-------------------------------------------------------------------------------------------
!> Purpose : Dump data for steel and rubber separately 
!!
!> @date June 16
!-------------------------------------------------------------------------------------------
!
subroutine paraview_custom_dump
!
      use environment, only : PREFIX
      use paraview   , only : PARAVIEW_DOMAIN
!
      implicit none
      character(len=128) :: PREFIX_tmp
!
!     save PREFIX
      PREFIX_tmp=PREFIX
!
!     dump paraview data
      PARAVIEW_DOMAIN=1
      PREFIX="Steel_"
      call paraview_custom_driver
      PARAVIEW_DOMAIN=2
      PREFIX="Rubber_"
      call paraview_custom_driver
!
!     reset
      PREFIX=PREFIX_tmp
      PARAVIEW_DOMAIN=0
!
end subroutine paraview_custom_dump