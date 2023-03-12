


#include "typedefs.h"

   subroutine print_orient


   use data_structure3D
   use element_data

   IMPLICIT NONE

!..element order and orientation
   integer :: nedge_orient(12), nface_orient(6), nv1, nv2, nv(2)
   integer :: nodesl(27), norientl(27), nrfn, nfver(4)
   integer :: mdle, iel,i ,nrv, nre, nrf,nface_or, nface_nodes(9),j,nod
   integer :: ntype

   mdle = 0
   do iel = 1, NRELES
      call nelcon(mdle,mdle)
      ntype=NODES(mdle)%ntype
      nrv = nvert(ntype)
      nre = nedge(ntype)
      nrf = nface(ntype)
      call elem_nodes(mdle, nodesl,norientl)
      call find_orient(mdle, nedge_orient,nface_orient)

      write(*,*) 'mdle = ', mdle
      write(*,*) 'edge orient = ',  nedge_orient
      write(*,*) 'face orient = ',  nface_orient

      do i = 1,nre

         if (nedge_orient(i) .ne. 0 ) then 
            call edge_to_vert(ntype,i, nv(1),nv(2))
            write(*,*) 'edge orient = ',  nedge_orient(i)
            write(*,*) 'edge vert   = ',  nodesl(nv(1:2))-1
            call pause
         endif   
      enddo   


      do i = 1, nrf
         nface_or = nface_orient(i)
         if (nface_or .ne. 0) then 
            call face_nodes(ntype,i, nface_nodes,nrfn)
            call face_to_vert_nos(MDLQ,nface_or, Nfver)
            write(*,*) 'nface_or   = ', nface_or
            write(*,*) 'nfver      = ', nfver
            write(*,*) 'face_nodes = ', nface_nodes(Nfver(1:4))
            write(*,*) 'face_vert  = ', nodesl(nface_nodes(Nfver(1:4)))-1
            call pause
         endif   
    
      enddo   
   enddo   


   end subroutine print_orient
