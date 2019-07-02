c----------------------------------------------------------------------
c
c   routine name       - verify_orient
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 10
c
c   purpose            - routine verifies consistency of node
c                        orientations as seen by elements
c
c   arguments          - none
c
c----------------------------------------------------------------------
c
      subroutine verify_orient
c
      use data_structure3D
      use element_data
c
#include "syscom.blk"
c
c  ...work space for elem_nodes
      dimension nodesl(27),norientl(27)
c
c  ...vertex nodes on an edge or face
      dimension nv(4)
c
      integer, allocatable ::  node_orient(:,:,:)
c     for a mid-edge or mid0face node 'nod' and i-th element
c     connected to the node,
c        node_orient(1,i,nod) = element (middle node) number
c        node_orient(2:3(4,5),i,nod) = vertices on the edge
c                               (face) listed in the order
c                               dictated by the global orientation
c                               (as seen by the element)
c
      allocate(node_orient(5,10,NRNODS))
      node_orient = 0
c
c  ...loop through elements and record mid-edge node and mid-face
c     node orientations as seen by the elements
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        nrv = nvert(NODES(mdle)%type)
        nre = nedge(NODES(mdle)%type)
        nrf = nface(NODES(mdle)%type)
        call elem_nodes(mdle, nodesl,norientl)
c
c  .....loop through the element edges
        do ie=1,nre
          in = nrv+ie
          nod = nodesl(in)
          ic=1
          do while(node_orient(1,ic,nod).ne.0)
            ic=ic+1
            if (ic.gt.10) then
              write(*,*) 'verify_orient: DIMENSION EXCEEDED'
              stop 1
            endif
          enddo
          node_orient(1,ic,nod) = mdle
          call edge_to_vert(Nodes(mdle)%type,ie, nv(1),nv(2))
          nv(1:2) = nodesl(nv(1:2))
          call apply_orient(norientl(in), NODES(nod)%type,
     .                      nv, node_orient(2:3,ic,nod))
        enddo
c
c  .....loop through the element faces
        do ifc=1,nrf
          in = nrv+nre+ifc
          nod = nodesl(in)
          ic = 1
          do while (node_orient(1,ic,nod).ne.0)
            ic=ic+1
            if (ic.gt.10) then
              write(*,*) 'verify_orient: DIMENSION EXCEEDED'
              stop 1
            endif
          enddo
          node_orient(1,ic,nod) = mdle
          call face_to_vert(Nodes(mdle)%type,ifc,
     .                      nv(1),nv(2),nv(3),nv(4))
          if (nv(1).eq.nv(4)) then
            nv(1:3) = nodesl(nv(1:3))
            call apply_orient(norientl(in), NODES(nod)%type,
     .                        nv, node_orient(2:4, ic,nod))

          else
            nv(1:4) = nodesl(nv(1:4))
            call apply_orient(norientl(in), NODES(nod)%type,
     .                        nv, node_orient(2:5, ic,nod))
          endif
        enddo
      enddo
c
c  ...check consistency of orientations
      do nod=1,NRNODS
        if (node_orient(1,1,nod).eq.0) cycle
c
        nflag=0
        ic=1
        do while(node_orient(1,ic,nod).ne.0)
          nsum = 0
          do iv=2,5
            nsum = nsum +
     .             abs(node_orient(iv,ic,nod)-node_orient(iv,1,nod))
          enddo
          if (nsum.gt.0) nflag = 1
          ic=ic+1
        enddo
c
        if (nflag.eq.1) then
          write(*,7011) nod,NODES(nod)%type
 7011     format(' verify_orient: nod,type = ',i6,2x,a4)
          ic=1
          do while(node_orient(1,ic,nod).ne.0)
            iaux=node_orient(1,ic,nod)
            write(*,7012) iaux,NODES(iaux)%type,node_orient(2:5,ic,nod)
 7012 format(' mdle = ',i6,' , type = ',a4,' ; VERTICES = ',4(i6,2x))
            ic=ic+1
          enddo
          call pause
        endif
      enddo
      deallocate(node_orient)
c
c
      end
c
c
c----------------------------------------------------------------------
c
c   routine name       - apply_orient
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 10
c
c   purpose            - apply orientation
c
c----------------------------------------------------------------------
c
      subroutine apply_orient(Norient, Type, Nv1, Nv2)
c
      use element_data , only : EDGE_L2G , TRIAN_L2G , QUADR_L2G
c
#include "syscom.blk"
c  ...workspace
      character(len=4) :: Type
      dimension Nv1(4), Nv2(4)
c
      iprint = 1
c
      select case (Type)
      case('medg')
        do iv=1,2
          Nv2(iv) = Nv1(EDGE_L2G(iv, Norient))
        enddo
      case('mdlt','tria')
        do iv=1,3
          Nv2(iv) = Nv1(TRIAN_L2G(iv, Norient))
        enddo
      case('mdlq','quad')
        do iv=1,4
          Nv2(iv) = Nv1(QUADR_L2G(iv, Norient))
        enddo
      end select
c
c
      end subroutine
