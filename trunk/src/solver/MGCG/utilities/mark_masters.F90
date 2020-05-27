!-----------------------------------------------------------------------
!
!    routine name       - mark_masters
!
!-----------------------------------------------------------------------
!
!    latest revision    - Jan 2018
!
!    purpose            - for a given mesh, routine marks all ACTIVE
!                         nodes in the mesh as `master nodes'
!
!    arguments          - none
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
#include "implicit_none.h"
!
   subroutine mark_masters(igrid)
!

   use data_structure3D,  only: NODES, MAXNODM, MAXNODS, NR_PHYSA
   use mg_data_structure, only: GRID

   use assembly,          only: assembly_begin, assembly_end
   use mg_data_structure, only: NODES_MG
!
   implicit none
!
!-----------------------------------------------------------------------
!
   integer, intent(in) :: igrid

!..work space for celem
   integer    :: nrdofs(NR_PHYSA)
   integer    :: nodm(MAXNODM),   ndofmH(MAXNODM)
   integer    :: ndofmE(MAXNODM), ndofmV(MAXNODM), ndofmQ(MAXNODM)
   integer    :: nrdofm, nrdofc, nrnodm
   integer    :: nod, mdle, iel, i, idec
   VTYPE      :: zvoid
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
   call assembly_begin
!
!..clean up and copy order
   do nod=1,MAXNODS
      NODES_MG(nod)%master(igrid)=0
!  ...save the current order
      NODES_MG(nod)%orderC(igrid) = NODES(nod)%order
   enddo
!
!..loop through active elements
   mdle=0
   do iel=1,GRID(igrid)%nreles
      mdle = GRID(igrid)%mdlel(iel)
!
!  ...determine nodes of the modified element
      idec=1
      call celem_mg(iel,-1,mdle,idec, nrdofs,nrdofm,nrdofc,    &
                    nodm,ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,   &
                    zvoid,zvoid)
!
!  ...loop through the nodes of the element and mark them
      do i=1,nrnodm
         nod = nodm(i)
         NODES_MG(nod)%master(igrid)=1
      enddo
      NODES_MG(mdle)%master(igrid)=1
   enddo
!
   call assembly_end
!
!
   end subroutine mark_masters
