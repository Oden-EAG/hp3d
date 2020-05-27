c----------------------------------------------------------------------
c
c   routine name       - find_nodes
c
c----------------------------------------------------------------------
c
c   latest revision    - Dec 07
c
c   purpose            - routine returns nodes for an element
c
c   arguments
c     in:
c            Mdle      - (middle node of) element
c     out:
c
c            Nel_nodes - element nodes
c
c----------------------------------------------------------------------
c
      subroutine find_nodes(Mdle, Nel_nodes)
c
      use data_structure3D
c
#include "syscom.blk"
c
      dimension Nel_nodes(27)
      dimension norientl(27)
c
c
      iprint=0
      call elem_nodes(Mdle, Nel_nodes,norientl)
      if (iprint.eq.1) then
        write(*,7001) Mdle
 7001   format('find_nodes: NODES FOR Mdle = ',i6)
        write(*,7002) Nel_nodes(1:15)
 7002   format(4i6,2x,6i6,2x,4i6,2x,i6)
      endif
c
      end


