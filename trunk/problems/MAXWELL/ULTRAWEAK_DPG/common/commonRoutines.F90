!
#include "typedefs.h"
!
!------------------------------------------------------------------------------
!> @brief      Propagate flag from father to son nodes; used to correctly
!!             inherit impedance BCs
!!
!> @param[in]  Icomp  - Physics component on which to propagate BC flag
!> @param[in]  Nflag  - flag to propagate
!!
!> @date       July 2023
!------------------------------------------------------------------------------
   subroutine propagate_flag(Icomp,Nflag)
!
      use data_structure3D
      use commonParam, only: IBCFLAG
!
      implicit none
!
      integer, intent(in) :: Icomp,Nflag
!
      integer :: ntype
      integer :: iel,mdle,ifc,nrfn,i,j,nod
!
!  ...element nodes and orientations, face nodes
      integer :: nodesl(27),norientl(27),nface_nodes(9)
!
!  ...element face BC flags, decoded BC flag for a node
      integer :: ibc(6,NRINDEX),nodflag(NRINDEX_HEV)
!
!----------------------------------------------------------------------------
!
      if (IBCFLAG .ne. 3) then
         write(*,*) 'propagate_flag called for IBCFLAG.ne.3, returning...'
         return
      endif
!
      if ((Icomp.lt.1) .or. (Icomp.gt.NRINDEX_HEV))
         write(*,*) 'propagate_flag: invalid Icomp = ', Icomp
         return
      endif
!
!..loop through active elements
!$OMP PARALLEL                                     &
!$OMP PRIVATE(ntype,mdle,ifc,nrfn,i,j,nod,nodesl,  &
!$OMP         norientl,nface_nodes,ibc,nodflag)
!
!$OMP DO
      do iel=1,NRELES
         mdle = ELEM_ORDER(iel)
         ntype = NODES(mdle)%ntype
!
!     ...determine element nodes
         call elem_nodes(mdle, nodesl,norientl)
!
!     ...get the element boundary conditions flags
         call find_bc(mdle, ibc)
!
!     ...loop through element faces
         do ifc=1,nface(ntype)
!
!        ...if face has a Dirichlet BC flag on this component,
!           then neither propagate Nflag from this face to its edges/vertices,
!           nor prohibit another face from passing Nflag to the edges/vertices.
            if (ibc(ifc,Icomp).eq.1) cycle
!
!        ...determine face node numbers
            call face_nodes(ntype,ifc, nface_nodes,nrfn)
!
!        ...loop through the face nodes
!$OMP CRITICAL
          do i=1,nrfn !-1
             j = nface_nodes(i)
             nod = nodesl(j)
!
!         ...if node belongs to a face that has impedance BC (Nflag),
!            then propagate the flag unless prohibited by another adjacent face
             if (ibc(ifc,Icomp).eq.Nflag) then
                if (NODES(nod)%visit.ne.-Nflag) then
                   NODES(nod)%visit = Nflag
                endif
!         ...prohibit the flag to be passed to the node
!            (if node belongs to a face that has no impedance or Dirichlet BC)
             else
                NODES(nod)%visit = -Nflag
             endif
          enddo
!$OMP END CRITICAL
         enddo
      enddo
!$OMP END DO
!
!  ...change -Nflag to zero
!$OMP DO
      do nod=1,NRNODS
         if (NODES(nod)%visit.eq.0) cycle
         call decod(NODES(nod)%bcond,2,NRINDEX_HEV, nodflag)
         if (NODES(nod)%visit.eq.-Nflag) then
            nodflag(Icomp) = 0
         elseif (NODES(nod)%visit.eq.Nflag) then
            nodflag(Icomp) = 1
         endif
         call encod(nodflag,2,NRINDEX_HEV, NODES(nod)%bcond)
!
      enddo
!$OMP END DO
!$OMP END PARALLEL
!
      call reset_visit
!
   end subroutine propagate_flag




!------------------------------------------------------------------------------
!> @brief      Evaluates unconstrained stiffness matrix and load vector
!!
!> @param[in]  Mdle     - element middle node number
!> @param[in]  x        - physical point to evaluate permittivity
!!
!> @param[in]  eps      - permittivity tensor at point
!!
!> @date       July 2023
!------------------------------------------------------------------------------
   subroutine get_permittivity(mdle,x, eps)
!
      use data_structure3D
      use parameters,          only: ZERO, ZONE
!
      implicit none
!
      integer,    intent(in)  :: mdle
      real(8),    intent(in)  :: x(3)
      complex(8), intent(out) :: eps(3,3)
!
      integer :: i
!
!------------------------------------------------------------------------------
!
! TODO: Implement your own custom permittivity here
!
!  ...set permittivity to identity for now.
      eps = ZERO
      do i=1,3
         eps(i,i) = ZONE
      enddo
!
   end subroutine get_permittivity
