c----------------------------------------------------------------------
c
c   routine name       - graphb
c
c----------------------------------------------------------------------
c
c   latest revision    - Jan 08
c
c   purpose            - routine controls interactively calls to
c                        various hp3D graphics routines
c
c   usage              - call graphb
c
c   arguments          - none
c
c----------------------------------------------------------------------
c
      subroutine graphbEM
c
      use graphmod
      use data_structure3D
c
#include "syscom.blk"

      save irout5, iwin, idec, ifirst
      data irout5/0/
c  ...initial lengths of the arrays
c  ...default values:
      if(irout5.ne.100)then
        irout5 = 100
        iwin = 5
        idec = 1
        ifirst = 1
c  ...initial lengths of the arrays
        MXIGTR   = 80000000
        MXIGSTR  = 160000000
        MXRGTRZ  = 80000000
      endif





c  ...select window
   12 write(*,*)'*****************************************************'
      write(*,*)' '
      write(*,*)'           R O U T I N E   G R A P H B E M'
      write(*,*)' '
      write(*,*)'*****************************************************'
      write(*,*)' '
      write(*,*)'WINDOW=?...................KIND OF WINDOW FOR GRAPHICS'
      write(*,*)'BODY.....................................BODY GRAPHICS'
      write(*,*)'MXIGTR   = ?..........................default=',MXIGTR
      write(*,*)'MXIGSTR  = ?..........................default=',MXIGSTR
      write(*,*)'MXRGTRZ  = ?..........................default=',MXRGTRZ
      write(*,*)'END...............................................EXIT'




      call ffldin(5,6,25,ier)
      if (ier.eq.1) go to 12
c---------------------------------------------------------------------
      call findop('WINDOW' ,noc,loc)
      if (noc.gt.0)
     .  then
            iwin = int(rfldata(loc))
        endif
c---------------------------------------------------------------------
      call findop('BODY' ,noc,loc)
      if (noc.gt.0)
     .  then
            idec = 1
        endif
c---------------------------------------------------------------------
      call findop('MXRGTRZ' ,noc,loc)
      if (noc.gt.0)
     .  then
            MXRGTRZ  = int(rfldata(loc))
        endif
c---------------------------------------------------------------------
      call findop('MXIGSTR' ,noc,loc)
      if (noc.gt.0)
     .  then
            MXIGSTR  = int(rfldata(loc))
        endif
c---------------------------------------------------------------------
      call findop('MXIGTR' ,noc,loc)
      if (noc.gt.0)
     .  then
            MXIGTR  = int(rfldata(loc))
        endif
c---------------------------------------------------------------------
      call findop('END' ,noc,loc)
      if (noc.gt.0)
     .  then

c      ...clear space if it was allocated:
        if(ifirst.eq.0)then
          deallocate( IGTRCU,  STAT = i1)
          deallocate( IGTRNO,  STAT = i2)
          deallocate( IGTR,    STAT = i3)
          deallocate( IGSTR,   STAT = i4)
          deallocate( RGTRZ,   STAT = i5)
          deallocate( RGTR,    STAT = i6)
          if( (i1+i2+i3+i4+i5+i6)/=0 )then
             print*,'GRAPH:'
             print*,'sth NOT deallocated i1-6 = ',i1,i2,i3,i4,i5,i6
             stop
          endif
c      ...indicate that space was deallocated:
          ifirst = 1
        endif
        return

        endif
c---------------------------------------------------------------------


c
c  ...select window
c      write(*,*) 'PLEASE SELECT KIND OF WINDOW FOR GRAPHICS'
c      write(*,*) '1-B&W, 2-Gray16, 3-Gray256, 4-Color16, 5-Color256'
c      read(*,*) iwin
      call initwin(iwin)
c
c  ...allocate ALL arrays used in graphics package
      if(ifirst.eq.1)then
       allocate(IGTRCU(MXIGTR))
       allocate(IGTRNO(MXIGTR))
       allocate(IGTR(MXIGTR))
       allocate(IGSTR(MXIGSTR))
       allocate(RGTRZ(MXIGTR))
       allocate(RGTR(MXRGTRZ))
c    ...indicate fact of allocation:
        ifirst = 0
      endif
c
c
c   10 continue
c      write(*,*)
c      write(*,*)'graphb: SELECT:'
c      write(*,*)'        RETURN..............................0'
c      write(*,*)'        BODY GRAPHICS.......................1'
c      read(*,*) idec



c
      select case(idec)
      case(0)
        deallocate(IGTRCU,IGTRNO,IGTR,IGSTR,RGTRZ,RGTR)
        return
c
      case(1)
        call grbodyEM(iwin)
      end select
      go to 12
c
c
      end

