c----------------------------------------------------------------------
c
c   routine name       - nodcor_vert
c
c----------------------------------------------------------------------
c
c   latest revision    - Oct 09
c
c   purpose            - routine returns vertex coordinates for
c                        a 3D element; this routine is for debugging
c                        mesh refinements only
c
c   arguments :
c     in:
c           Mdle       - middle node of an element
c     out:
c           Xnod       - the element local geometry dof
c
c   required  routines - 
c
c----------------------------------------------------------------------
c
      subroutine nodcor_vert(Mdle, Xnod)
c
      use data_structure3D 
      use element_data
#include "syscom.blk"
c
      dimension Xnod(NDIMEN,MAXbrickH)
c
c  ...element nodes and their orientation
      dimension nodesl(27),norientl(27)
c
      iprint=0
c      
ccc      if (Mdle.eq.147) iprint=1
ccc      else iprint=0
      if (iprint.eq.1) write(*,*) 'nodcor_vert: ENTERING'
c
      call elem_nodes(Mdle, nodesl,norientl)
c
      do iv=1,nvert(NODES(Mdle)%type)
        nod = nodesl(iv)
c  .....following initialization is necessary owing to compiler bug
c       don't know why it is, but without touching Xnod, it keep 
c       generating segmentation error...
        Xnod(1:NDIMEN,iv) = 0.d0
        Xnod(1:NDIMEN,iv) = NODES(nod)%coord(1:NDIMEN,1)
      enddo
c
      if (iprint.eq.1) then
        write(*,7001) Mdle
 7001   format('nodcor_vert: VERTEX COORDINATES FOR Mdle = ',i6)
        do iv=1,nvert(NODES(Mdle)%type)
          write(*,7002) iv,Xnod(1:NDIMEN,iv)
 7002     format('iv = ',i2,2x,3f8.3)
        enddo
      endif
      if (iprint.eq.1) write(*,*) 'nodcor_vert: EXITING'
c
c
      end subroutine
