!
#include "typedefs.h"
!
!----------------------------------------------------------------------------
!
!   routine name       - propagate_flag
!
!----------------------------------------------------------------------------
!
!   latest revision    - June 2021
!
!   purpose            - Propagates Nflag (customized BC flag 2-9) from
!                        element faces to element edges and vertices.
!                        The flag is passed provided ALL adjacent faces
!                        share the flag (or have a Dirichlet flag).
!                        Since a node stores only a 0/1 BC flag per component,
!                        propagating Nflag corresponds to setting BC flag = 1
!                        on the corresponding component (Icomp).
!                        Background:
!                        Impedance boundary condition (BC) should not be
!                        inherited by edges and vertices from a face
!                        during refinement, unless the edge/vertex is only
!                        adjacent to impedance and dirichlet faces.
!
!   arguments:
!     in:
!              Icomp   - Physics attribute component number (1,..,NRINDEX)
!              Nflag   - A custom BC flag (2-9); e.g., impedance BC flag
!
!----------------------------------------------------------------------------
!
subroutine propagate_flag(Icomp,Nflag)
!
   use data_structure3D
   use commonParam, only: IBCFLAG
!
   implicit none
!
   integer, intent(in) :: Icomp,Nflag
!
   character(len=4) :: etype
   integer :: iel,mdle,ifc,nrfn,i,j,nod
!
!..element nodes and orientations, face nodes
   integer :: nodesl(27),norientl(27),nface_nodes(9)
!
!..element face BC flags, decoded BC flag for a node
   integer :: ibc(6,NRINDEX),nodflag(NRINDEX)
!
!----------------------------------------------------------------------------
!
   if (IBCFLAG .ne. 3) then
      write(*,*) 'propagate_flag called for IBCFLAG.ne.3, returning...'
      return
   endif
!
!..loop through active elements
!$OMP PARALLEL                                     &
!$OMP PRIVATE(etype,mdle,ifc,nrfn,i,j,nod,nodesl,  &
!$OMP         norientl,nface_nodes,ibc,nodflag)
!
!$OMP DO
   do iel=1,NRELES
      mdle = ELEM_ORDER(iel)
      etype = NODES(mdle)%type
!
!  ...determine element nodes
      call elem_nodes(mdle, nodesl,norientl)
!
!  ...get the element boundary conditions flags
      call find_bc(mdle, ibc)
!
!  ...loop through element faces
      do ifc=1,nface(etype)
!
!     ...if face has a Dirichlet BC flag on this component,
!        then neither propagate Nflag from this face to its edges/vertices,
!        nor prohibit another face from passing Nflag to the edges/vertices.
         if (ibc(ifc,Icomp).eq.1) cycle
!
!     ...determine face node numbers
         call face_nodes(etype,ifc, nface_nodes,nrfn)
!
!     ...loop through the face nodes
!$OMP CRITICAL
         do i=1,nrfn !-1
            j = nface_nodes(i)
            nod = nodesl(j)
!
!        ...if node belongs to a face that has impedance BC (Nflag),
!           then propagate the flag unless prohibited by another adjacent face
            if (ibc(ifc,Icomp).eq.Nflag) then
               if (NODES(nod)%visit.ne.-Nflag) then
                  NODES(nod)%visit = Nflag
               endif
!        ...prohibit the flag to be passed to the node
!           (if node belongs to a face that has no impedance or Dirichlet BC)
            else
               NODES(nod)%visit = -Nflag
            endif
         enddo
!$OMP END CRITICAL
      enddo
   enddo
!$OMP END DO
!
!..change -Nflag to zero
!$OMP DO
   do nod=1,NRNODS
      if (NODES(nod)%visit.eq.0) cycle
      call decod(NODES(nod)%bcond,2,NRINDEX, nodflag)
      if (NODES(nod)%visit.eq.-Nflag) then
         nodflag(Icomp) = 0
      elseif (NODES(nod)%visit.eq.Nflag) then
         nodflag(Icomp) = 1
      endif
      call encod(nodflag,2,NRINDEX, NODES(nod)%bcond)
!
   enddo
!$OMP END DO
!$OMP END PARALLEL
!
   call reset_visit
!
end subroutine propagate_flag
