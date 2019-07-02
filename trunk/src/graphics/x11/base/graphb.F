c----------------------------------------------------------------------
c> @Purpose  : routine controls interactively calls to various hp3d
c               graphics routines
c
c> @Revision : May 12
c----------------------------------------------------------------------
c
      subroutine graphb
c
      use graphmod
ccc      use data_structure3D
c
#include "syscom.blk"
c
c  ...initialized dimensions for workspaces, if not set by the user
      if (.NOT. INITIALIZED) then
c                                MXIGTR // MXIGSTR // MXRGTRZ 
        call set_x11_workspace(10000000,10000000,10000000)
      endif
c
c  ...allocate workspaces for graphics      
      allocate(IGTRCU(MXIGTR ), stat=istat)
      if (istat.ne.0) then
        write(*,*)'graphb: IGTRCU not allocated!' ; stop
      endif
      allocate(IGTRNO(MXIGTR ), stat=istat)
      if (istat.ne.0) then
        write(*,*)'graphb: IGTRNO not allocated!' ; stop
      endif
      allocate(IGTR(  MXIGTR ), stat=istat)
      if (istat.ne.0) then
        write(*,*)'graphb: IGTR not allocated!'   ; stop
      endif
      allocate(IGSTR( MXIGSTR), stat=istat)
       if (istat.ne.0) then
        write(*,*)'graphb: IGSTR not allocated!'  ; stop
      endif
      allocate(RGTRZ( MXIGTR ), stat=istat)
        if (istat.ne.0) then
        write(*,*)'graphb: RGTRZ not allocated!'  ; stop
      endif
      allocate(RGTR(  MXRGTRZ), stat=istat)
      if (istat.ne.0) then
        write(*,*)'graphb: RGTR not allocated!'   ; stop
      endif
c
c  ...select window
ccc      write(*,*) 'PLEASE SELECT KIND OF WINDOW FOR GRAPHICS'
ccc      write(*,*) '1-B&W, 2-Gray16, 3-Gray256, 4-Color16, 5-Color256'
ccc      read(*,*) iwin
      iwin = 5 ! who uses anything but 5 these days? - jz
      call initwin(iwin)
c
c  ...display in infinite loop      
   10 continue
      write(*,*)
      write(*,*)'graphb: SELECT:'
      write(*,*)'        RETURN..............................0'
      write(*,*)'        BODY GRAPHICS.......................1'
ccc      write(*,*)'        INTERACTIVE MESH MODIFICATION.......4'

      read(*,*) idec
c
      select case(idec)
      case(0)
        deallocate(IGTRCU,IGTRNO,IGTR,IGSTR,RGTRZ,RGTR, stat=istat)
        if (istat.ne.0) then
          write(*,*)'graphb: workspaces not deallocated!' ; stop
        endif
        return
c
      case(1)
        call grbody(iwin)
      endselect

      goto 10
c
c
      end
