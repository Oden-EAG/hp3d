c----------------------------------------------------------------------
c
c   routine name       - find_orient
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 10
c
c   purpose            - routine returns orientation for an element
c                        edge and face nodes
c
c   arguments 
c     in:
c         Mdle         - (middle node of) element
c     out:
c                
c         Nedge_orient - egde orientations (0,1)
c         Nface_orient - face orientations (0,1,...,7)
c
c----------------------------------------------------------------------
c
      subroutine find_orient(Mdle, Nedge_orient,Nface_orient)
c
      use data_structure3D
c
#include "syscom.blk"
c
      dimension Nedge_orient(12), Nface_orient(6)
c
c  ...element nodes and their orientation
      dimension nodesl(27),norientl(27)
c
      call elem_nodes(Mdle, nodesl,norientl)
      nrv = NVERT(NODES(Mdle)%type)
      nre = NEDGE(NODES(Mdle)%type)
      nrf = NFACE(NODES(Mdle)%type)
      Nedge_orient(1:nre) = norientl(nrv+1:nrv+nre)
      Nface_orient(1:nrf) = norientl(nrv+nre+1:nrv+nre+nrf)
c
c
      end


      
