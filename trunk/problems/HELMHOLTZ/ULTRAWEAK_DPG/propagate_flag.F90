!----------------------------------------------------------------------
!> @brief       Routine propagates impedance BC from parent to son nodes.
!!              This is an advanced routine for enforcing impedance
!!              boundary conditions using the elimination strategy; see
!!              notes in parallel hp-book.
!!
!> @param[in]   icomp   - component to propagate flag
!> @param[in]   Nflag   - flag to propagate
!!
!> @date        July 2023
!----------------------------------------------------------------------
   subroutine propagate_flag(Icomp,Nflag)
!
      use data_structure3D
      use common_prob_data_UW, only: IBC_PROB
      use par_mesh,            only: DISTRIBUTED
      use MPI,                 only: MPI_IN_PLACE, MPI_INTEGER, MPI_SUM,   &
                                     MPI_COMM_WORLD
!
      implicit none
!
      integer, intent(in) :: Icomp,Nflag
!
      integer :: ntype
      integer :: iel,mdle,ifc,nrfn,i,j,nod,ierr
!
!  ...element nodes and orientations, face nodes
      integer :: nodesl(27),norientl(27),nface_nodes(9)
!
!  ...element face BC flags, decoded BC flag for a node
      integer :: ibc(6,NRINDEX),nodflag(NRINDEX)
!
!-------------------------------------------------------------------------------
!
      select case(IBC_PROB)
      case(3,4,5,6)
         continue
      case default
         write(*,*) 'propagate_flag called for non-impedance BC, returning...'
         return
      end select
!
      call reset_visit
!
!..loop through active elements
!$OMP PARALLEL                                     &
!$OMP PRIVATE(ntype,mdle,ifc,nrfn,i,j,nod,nodesl,  &
!$OMP         norientl,nface_nodes,ibc,nodflag)
!$OMP DO
      do iel=1,NRELES_SUBD
         mdle = ELEM_SUBD(iel)
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
            do i=1,nrfn-1
               j = nface_nodes(i)
               nod = nodesl(j)
   !
               NODES(nod)%visit = -Nflag
            enddo
   !
            if (ibc(ifc,Icomp).eq.Nflag) then
               nod = nodesl(nface_nodes(nrfn))
               NODES(nod)%visit = Nflag
            endif
!$OMP END CRITICAL
         enddo
      enddo
!$OMP END DO
!$OMP END PARALLEL
!
      if (DISTRIBUTED) then
!     ...Sum used since values are positive and negative;
!        value isn't important, only sign and whether non-zero
         call MPI_Allreduce(MPI_IN_PLACE,NODES(:)%visit,NRNODS,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD, ierr)
      endif
!
!  ...change -Nflag to zero
      do nod=1,NRNODS
         if (NODES(nod)%visit.eq.0) cycle
         call decod(NODES(nod)%bcond,2,NRINDEX, nodflag)
         if (NODES(nod)%visit.le.-Nflag) then
            nodflag(Icomp) = 0
         elseif (NODES(nod)%visit.ge.Nflag) then
            nodflag(Icomp) = 1
         endif
         call encod(nodflag,2,NRINDEX, NODES(nod)%bcond)
!
      enddo
!
      call reset_visit
!
   end subroutine propagate_flag


