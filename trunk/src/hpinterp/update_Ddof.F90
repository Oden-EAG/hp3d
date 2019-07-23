!-----------------------------------------------------------------------
!
!    routine name       - update_Ddof
!
!-----------------------------------------------------------------------
!
!    latest revision    - July 2019
!
!    purpose            - routine updates values of solution degrees
!                         of freedom for Dirichlet nodes
!
!    arguments          - none
!
!-----------------------------------------------------------------------
!
subroutine update_Ddof()
!
   use data_structure3D
   use environment, only: QUIET_MODE
   use par_mesh   , only: DISTRIBUTED
   use MPI        , only: MPI_COMM_WORLD
   use mpi_param  , only: RANK,ROOT
!
   implicit none
!
#include "implicit_none.h"
!
!..orientation of element nodes
   integer, dimension(12) :: nedge_orient
   integer, dimension(6)  :: nface_orient
!
!..element nodes
   integer, dimension(27) :: nodesl, norientl
!
!..order of approximation for an element
   integer, dimension(19) :: norder
!
!..reference coordinates for an element
   real*8, dimension(3,8) :: xsub
!
!..solution dof for an element
   VTYPE, dimension(MAXEQNH,MAXbrickH) :: zdofH
   VTYPE, dimension(MAXEQNE,MAXbrickE) :: zdofE
   VTYPE, dimension(MAXEQNV,MAXbrickV) :: zdofV
   VTYPE, dimension(MAXEQNQ,MAXbrickQ) :: zdofQ
!
!..auxiliary variables for timing
   real(8) :: MPI_Wtime,start_time,end_time
!
!..auxiliary array with active mdle nodes
   integer :: mdlel(NRELES)
!
!..auxiliary variables
   integer :: iel, iv, ie, ifc, ind, iflag, ierr
   integer :: k, mdle, nf, no, nod
!
!..number of elements in subdomain
   integer :: nreles_subd, subd
!
!-----------------------------------------------------------------------
!
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
   start_time = MPI_Wtime()
!
!..fetch active elements
   mdle = 0;
   if (DISTRIBUTED) then
      nreles_subd = 0
      do iel=1,NRELES
         call nelcon(mdle, mdle)
         call get_subd(mdle, subd)
         if (subd .eq. RANK) then
            nreles_subd = nreles_subd + 1
            mdlel(nreles_subd) = mdle
         endif
      enddo
   else
      do iel=1,NRELES
         call nelcon(mdle, mdle)
         mdlel(iel) = mdle
      enddo
      nreles_subd = NRELES
   endif
!
!-----------------------------------------------------------------------
!..Begin parallel OpenMP environment
!..Note that if Dirichlet vertices or edges are shared between elements,
!  then concurrent writes to NODES(nod)&zdofH,E may happen,
!  but they will write the exact same values.
! TODO this appears to make problems (observed for dirichlet data on Hcurl
!      flux with multiple non-zero components Ex,Ey)
!      ..added omp critical sections
!      ..remove critical sections if possible
!      ..concurrent writes in dhpedge, etc., seem to be the issue
!      ..could place OMP CRTICAL there
!
!$OMP PARALLEL                                                 &
!$OMP PRIVATE(mdle,iflag,no,xsub,nodesl,norientl,nedge_orient, &
!$OMP         nface_orient,norder,zdofH,zdofE,zdofV,zdofQ,     &
!$OMP         iel,iv,ie,ifc,k,ind,nod)
!
!  Step 1: Update   V E R T   dof for Dirichlet nodes
!
!$OMP DO SCHEDULE(DYNAMIC)
   do iel=1,nreles_subd
      mdle = mdlel(iel)
      call refel(mdle, iflag,no,xsub)
      call elem_nodes(mdle, nodesl,norientl)
      do iv=1,nvert(NODES(mdle)%type)
         nod = nodesl(iv)
         if (nod.lt.0) cycle
         if (is_dirichlet(nod)) then
!        ...cycle if H1 dofs are not supported by the node
            if (.not.associated(NODES(nod)%zdofH)) cycle
!
            if (is_dirichlet_homogeneous(nod)) then
               NODES(nod)%zdofH = ZERO
            else
               !$OMP CRITICAL
               call dhpvert(mdle,iflag,no,xsub(1:3,iv),NODES(nod)%case, &
                  NODES(nod)%zdofH)
               !$OMP END CRITICAL
            endif
         endif
      enddo
   enddo
!$OMP END DO
!
!-----------------------------------------------------------------------
!
!  Step 2: Update   E D G E   dof for Dirichlet nodes
!
!$OMP DO SCHEDULE(DYNAMIC)
   do iel=1,nreles_subd
      mdle = mdlel(iel)
      call refel(mdle, iflag,no,xsub)
      call elem_nodes(mdle, nodesl,norientl)
      call find_orient(mdle, nedge_orient,nface_orient)
      call find_order(mdle, norder)
!
!  ...compute solution dofs (need for H1 update)
      call solelm(mdle, zdofH,zdofE,zdofV,zdofQ)
      do ie=1,nedge(NODES(mdle)%type)
         ind = nvert(NODES(mdle)%type)+ie
         nod = nodesl(ind)
!
         if (nod.lt.0) cycle
         if (is_dirichlet(nod)) then
!
!        ...update H1 Dirichlet dofs
            if (associated(NODES(nod)%zdofH)) then
               if (is_dirichlet_homogeneous(nod)) then
                  NODES(nod)%zdofH = ZERO
               else
                  !$OMP CRITICAL
                  call dhpedgeH(mdle,iflag,no,xsub,                  &
                                NODES(mdle)%type,NODES(nod)%case,    &
                                nedge_orient,nface_orient,norder,ie, &
                                zdofH, NODES(nod)%zdofH)
                  !$OMP END CRITICAL
               endif
            endif
!
!        ...update H(curl) Dirichlet dofs
            if (associated(NODES(nod)%zdofE)) then
               if (is_dirichlet_homogeneous(nod)) then
                  NODES(nod)%zdofE = ZERO
               else
                  !$OMP CRITICAL
                  call dhpedgeE(mdle,iflag,no,xsub,                  &
                                NODES(mdle)%type,NODES(nod)%case,    &
                                nedge_orient,nface_orient,norder,ie, &
                                NODES(nod)%zdofE)
                  !$OMP END CRITICAL
               endif
            endif
         endif
      enddo
   enddo
!$OMP END DO
!
!-----------------------------------------------------------------------
!
!  Step 3: Update   F A C E   dof for Dirichlet nodes
!
!$OMP DO SCHEDULE(DYNAMIC)
   do iel=1,nreles_subd
      mdle = mdlel(iel)
      call refel(mdle, iflag,no,xsub)
      call elem_nodes(mdle, nodesl,norientl)
      call find_orient(mdle, nedge_orient,nface_orient)
      call find_order(mdle, norder)
!
!  ...compute solution dofs (needed for H1 and Hcurl update)
      call solelm(mdle, zdofH,zdofE,zdofV,zdofQ)
!
!  ...loop over faces
      do ifc=1,nface(NODES(mdle)%type)
!
!     ...get local node number
         ind = nvert(NODES(mdle)%type)+nedge(NODES(mdle)%type)+ifc
!
!     ...get global node number
         nod = nodesl(ind)
         if (nod.lt.0) cycle
         if (is_dirichlet(nod)) then
!
!           -- H1 --
            if (associated(NODES(nod)%zdofH)) then
               if (is_dirichlet_homogeneous(nod)) then
                  NODES(nod)%zdofH = ZERO
               else
                  !$OMP CRITICAL
                  call dhpfaceH(mdle,iflag,no,xsub,                   &
                                NODES(mdle)%type,NODES(nod)%case,     &
                                nedge_orient,nface_orient,norder,ifc, &
                                zdofH, NODES(nod)%zdofH)
                  !$OMP END CRITICAL
               endif
            endif
!
!           -- H(curl) --
            if (associated(NODES(nod)%zdofE)) then
               if (is_dirichlet_homogeneous(nod)) then
                  NODES(nod)%zdofE = ZERO
               else
                  !$OMP CRITICAL
                  call dhpfaceE(mdle,iflag,no,xsub,                   &
                                NODES(mdle)%type,NODES(nod)%case,     &
                                nedge_orient,nface_orient,norder,ifc, &
                                zdofE, NODES(nod)%zdofE)
                  !$OMP END CRITICAL
               endif
            endif
!
!           -- H(div) --
            if (associated(NODES(nod)%zdofV)) then
               if (is_dirichlet_homogeneous(nod)) then
                  NODES(nod)%zdofV = ZERO
               else
                  !$OMP CRITICAL
                  call dhpfaceV(mdle,iflag,no,xsub,                   &
                                NODES(mdle)%type,NODES(nod)%case,     &
                                nedge_orient,nface_orient,norder,ifc, &
                                NODES(nod)%zdofV)
                  !$OMP END CRITICAL
               endif
            endif
         endif
!  ...end of loop over faces
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL

!
!-----------------------------------------------------------------------
!
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
   if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) then
      end_time = MPI_Wtime()
      write(*,8010) end_time-start_time
 8010 format(' update_Ddof: ',f12.5,'  seconds',/)
   endif
!
end subroutine update_Ddof
