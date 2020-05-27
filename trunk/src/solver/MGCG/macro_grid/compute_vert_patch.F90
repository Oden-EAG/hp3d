! -----------------------------------------------------------------------
!
!    routine name       - compute_vert_patch
!
! -----------------------------------------------------------------------
!
!    latest revision    - MAY 18
!
!    purpose            - computes the list of the active vertex NODES
!                         for the current mesh. Each active vertex NODES
!                         defines one patch of the additive Swartz smoother
!    in:                -
!            Igrid      - grid index
!
!-----------------------------------------------------------------------
!
#include "implicit_none.h"
!
   subroutine compute_vert_patch(Igrid)
!
   use data_structure3D,   ONLY: NODES
   use mg_data_structure,  ONLY: GRID, NODES_MG, mg_reset_visit
   use patch_info,         ONLY: CGRID_VERTICES, compute_patch_mdle
   use assembly,           ONLY: MAXNODM,NR_PHYSA, assembly_begin, assembly_end
!
   implicit none
!
!-----------------------------------------------------------------------
!
   integer, intent(in) :: Igrid
!
!..work space for celem
!..number of variables for each physics attribute for an element
   integer :: nrdofs(NR_PHYSA)
   integer :: nodm(MAXNODM),  ndofmH(MAXNODM), &
              ndofmE(MAXNODM),ndofmV(MAXNODM),ndofmQ(MAXNODM)

   integer :: nrpatch, iel,mdle,nrdofm,nrdofc,nrnodm
   integer :: i, nod
   VTYPE   :: zvoid(1)

!
!-----------------------------------------------------------------------
!

!
!..create new list of vertices for the refined mesh
   call mg_reset_visit
!
   nrpatch = 0
!
   call assembly_begin
!
!..loop through coarse grid elements
   do iel=1,GRID(Igrid+1)%nreles
!
!  ...pick up mdle node number
      mdle = GRID(Igrid+1)%mdlel(iel)
!
!  ...get information from celem
      call celem_mg(-1,-1,mdle,1,nrdofs,nrdofm,nrdofc,nodm,ndofmH, &
                    ndofmE,ndofmV,ndofmQ,nrnodm,zvoid,zvoid)
!
!  ...loop through active nodes
      do i = 1,nrnodm
!
!     ...pick up the node number
         nod = nodm(i)
!
!     ...avoid repetition
         if (NODES_MG(nod)%visit.gt.0) cycle
!
!     ...check if the node is a vertex
         if (NODES(nod)%type .eq. 'vert') then
            nrpatch = nrpatch + 1
         endif
!
!     ...raise visitation flag
         NODES_MG(nod)%visit = 1
      enddo
!
!..end of loop through coarse grid elements
   enddo

   call mg_reset_visit
!
   if (allocated(CGRID_VERTICES)) deallocate(CGRID_VERTICES)
!
   allocate(CGRID_VERTICES(nrpatch))
   GRID(Igrid+1)%nrpatch = nrpatch
   nrpatch = 0
!
   do iel=1,GRID(Igrid+1)%nreles
      mdle = GRID(Igrid+1)%mdlel(iel)
!  ...get information from celem
      call celem_mg(-1,-1,mdle,1,nrdofs,nrdofm,nrdofc,nodm,ndofmH, &
                    ndofmE,ndofmV,ndofmQ,nrnodm,zvoid,zvoid)
!
!  ...loop through active nodes
      do i = 1,nrnodm
!
!     ...pick up the node number
         nod = nodm(i)
!
!     ...avoid repetition
         if (NODES_MG(nod)%visit.ne.0) cycle
!
!     ...check if the node is a vertex
         if (NODES(nod)%type .eq. 'vert') then
            nrpatch = nrpatch + 1
            CGRID_VERTICES(nrpatch) = nod
            NODES_MG(nod)%visit = nrpatch
         else
!        ...raise visitation flag
            NODES_MG(nod)%visit = -1
         endif
      enddo
   enddo
!
   call assembly_end
!
   call compute_patch_mdle(igrid+1)
!
!
   end subroutine compute_vert_patch
