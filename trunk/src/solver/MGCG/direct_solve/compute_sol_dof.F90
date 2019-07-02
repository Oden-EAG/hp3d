! ----------------------------------------------------------------------
!
!    routine name       - compute_sol_dof
!
! ----------------------------------------------------------------------
#include 'implicit_none.h'
subroutine compute_sol_dof(Igrid)
!
   use mg_data_structure, only: GRID
   use assembly,          only: assembly_begin, assembly_end
   use macro_grid_info,   only: ZSOL_C
!
   implicit none
!   
!-----------------------------------------------------------------------
!
   integer, intent(in) :: Igrid
!
   integer :: iel, mdle
! 
!----------------------------------------------------------------------
!
!..Step 1: determine the first dof offsets for active nodes
   allocate(ZSOL_C(GRID(Igrid)%nreles))
   call assembly_begin
!
   do iel = 1, GRID(Igrid)%nreles
      mdle = GRID(Igrid)%mdlel(iel)
!  ...get information from celem
      call compute_local_sol_dof(iel,mdle)
   enddo
!
   call assembly_end
!
end subroutine compute_sol_dof
!
!
! ----------------------------------------------------------------------
!
!    routine name       - compute_local_sol_dof
!
! ----------------------------------------------------------------------
subroutine compute_local_sol_dof(Iel,Mdle)
!
   use data_structure3D,  only: NODES,MAXNODM,ZERO,get_index
   use macro_grid_info,   only: ZSOL_C
   use assembly,          only: MAXNODM,NR_RHS
   use physics,           only: NRHVAR,NREVAR,NRVVAR,NRINDEX,NR_PHYSA, &
                                NR_COMP,DTYPE,PHYSAm
!
   implicit none
!
!-----------------------------------------------------------------------
!
   integer, intent(in) :: Iel,Mdle
!
!..locals
!..work space for celem
!..number of variables for each physics attribute for an element
   integer, dimension(MAXNODM) :: nodm,ndofmH,ndofmE,ndofmV,ndofmQ
   integer                     :: nrdofm, nrdofc, nrnodm
   integer                     :: nrdofs(NR_PHYSA)
!
!..decoded index for a node
   integer :: index(NRINDEX)
!
!..active component counter
   integer :: nvarH,nvarE,nvarV
!
!..auxiliary variables
   integer :: nrPhysH,nrPhysE,nrPhysV,nrPhysHE,nrPhysHEV
   integer :: nrVarHE,nrVarHEV
   integer :: nn, load, ivar, i, j, k, nod, iphys, icomp, il
   integer :: ndof, ndofH, ndofE, ndofV
   VTYPE   :: zvoid(1)
! 
!..component counters for the nodes (used in case of multiple loads)
   integer, dimension(MAXNODM) :: mvarH,mvarE,mvarV
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!..initiate component counters
   mvarH(1:MAXNODM)=0
   mvarE(1:MAXNODM)=0
   mvarV(1:MAXNODM)=0
!
!..determine nodes of the modified element
   call celem_mg(Iel,-1,Mdle,1, nrdofs,nrdofm,nrdofc,nodm,ndofmH, &
                 ndofmE,ndofmV,ndofmQ,nrnodm,zvoid,zvoid)
!
!..count number of variables for each physics type
   nrPhysH=0; nrPhysE=0; nrPhysV=0
   do iphys=1,NR_PHYSA
      select case(DTYPE(iphys))
         case('contin')
            nrPhysH=nrPhysH+1
         case('tangen')
            nrPhysE=nrPhysE+1
         case('normal')
            nrPhysV=nrPhysV+1
         case default
      end select
   enddo
   nrPhysHE  = nrPhysH+nrPhysE
   nrPhysHEV = nrPhysH+nrPhysE+nrPhysV
!
   nrVarHE  = NRHVAR+NREVAR
   nrVarHEV = NRHVAR+NREVAR+NRVVAR
!
   ndofH = sum(ndofmH(1:nrnodm))
   ndofE = sum(ndofmE(1:nrnodm))
   ndofV = sum(ndofmV(1:nrnodm))
!
   ndof  = ndofH + ndofE + ndofV
!
   ZSOL_C(Iel)%ndof_coarse = ndof
   allocate(ZSOL_C(Iel)%coarse(ndof))
   ZSOL_C(Iel)%coarse = ZERO
!   
!..initiate dof counter
   nn=0
!
!..initialize the number of H1,H(curl),H(div),L2 stored so far
   ivar=0
!
!..loop through loads
   do load=1,NR_RHS
!
!.....H1 dof .................................
      if (NRHVAR.eq.0) goto 200
!
!  ...loop through nodes of modified element (excluding mdle node)
      do i=1,nrnodm
         nod = nodm(i)
!
!     ...compute the number of active H1 variables for the node
         call get_index(nod, index)
!
         nvarH=0; k=0
         do iphys=1,nrPhysH
            il = NR_COMP(iphys)
            if (.not. PHYSAm(iphys)) then
               k=k+il
               cycle
            endif
            do icomp=1,il
               k=k+1
               if (index(k).eq.2) nvarH=nvarH+1
            enddo
         enddo
         if (nvarH.eq.0) cycle
!
!     ...loop through the nodal dof (potentially NONE)
         do j=1,ndofmH(i)/nvarH
            ivar=mvarH(i); k=0
! 
!        ...loop through physics variables              
            do iphys=1,nrPhysH
               il = NR_COMP(iphys)
               if (index(k+1) .eq. 0) then
                  k=k+il
                  cycle
               endif
               if (.not. PHYSAm(iphys)) then
                  k=k+il; ivar=ivar+il
                  cycle
               endif
!
!           ...loop through components of this physics variable
               do icomp=1,il
                  k=k+1; ivar=ivar+1
                  select case(index(k))
                  case(2)
                     nn=nn+1
!
!                 ...copy the dof
                     ZSOL_C(Iel)%coarse(nn) = NODES(nod)%zdofH(ivar,j)
                  case default
                  end select
!           ...end loop over components
               enddo
!        ...end loop over physics variables
            enddo
!     ...end loop over nodal dof
         enddo
!
!     ...update the number of components stored so far
         mvarH(i) = ivar
!  ...end loop over nodes
      enddo
!
!-----------------------------------------------------------------------
!
!.....H(curl) dof .................................
  200 if (NREVAR.eq.0) goto 300
!
!  ...loop through nodes of modified element (excluding mdle node)
      do i=1,nrnodm
         nod = nodm(i)
!
!     ...compute the number of active H1 variables for the node
         call get_index(nod, index)
         nvarE=0; k=NRHVAR
         do iphys=nrPhysH+1,nrPhysHE
            il = NR_COMP(iphys)
            if (.not. PHYSAm(iphys)) then
               k=k+il
               cycle
            endif
            do icomp=1,il
               k=k+1
               if (index(k).eq.4) nvarE=nvarE+1
            enddo
         enddo
         if (nvarE.eq.0) cycle
!
!     ...loop through the nodal dof
         do j=1,ndofmE(i)/nvarE
            ivar=mvarE(i); k=NRHVAR
! 
!        ...loop through physics variables              
            do iphys=nrPhysH+1,nrPhysHE
               il = NR_COMP(iphys)
               if (index(k+1) .eq. 0) then
                  k=k+il
                  cycle
               endif
               if (.not. PHYSAm(iphys)) then
                  k=k+il; ivar=ivar+il
                  cycle
               endif
!           ...loop through components of this physics variable
               do icomp=1,il
                  k=k+1; ivar=ivar+1
                  select case(index(k))
                  case(4)
                     nn=nn+1
!
!                 ...copy the dof
                     ZSOL_C(Iel)%coarse(nn) = NODES(nod)%zdofE(ivar,j)
                  case default
                  end select
!           ...end loop over components
               enddo
!        ...end loop over physics variables
            enddo
!     ...end loop over nodal dof
         enddo
!     ...update the number of components stored so far
         mvarE(i) = ivar
!  ...end loop over nodes
      enddo
!
!
!-----------------------------------------------------------------------
!
!  .....H(div) dof .................................
  300 if (NRVVAR.eq.0) goto 400
!
!  ...loop through nodes of modified element (excluding mdle node)
      do i=1,nrnodm
         nod = nodm(i)
!
!     ...compute the number of active H1 variables for the node
         call get_index(nod, index)
         nvarV=0; k=nrVarHE
         do iphys=nrPhysHE+1,nrPhysHEV
            il = NR_COMP(iphys)
            if (.not. PHYSAm(iphys)) then
               k=k+il
               cycle
            endif
            do icomp=1,il
               k=k+1
               if (index(k).eq.6) nvarV=nvarV+1
            enddo
         enddo
         if (nvarV.eq.0) cycle
!
!     ...loop through the nodal dof
         do j=1,ndofmV(i)/nvarV
!
!        ...loop through the components
            ivar=mvarV(i); k=nrVarHE
! 
!        ...loop through physics variables              
            do iphys=nrPhysHE+1,nrPhysHEV
               il = NR_COMP(iphys)
               if (index(k+1) .eq. 0) then
                  k=k+il
                  cycle
               endif
               if (.not. PHYSAm(iphys)) then
                  k=k+il; ivar=ivar+il
                  cycle
               endif
!           ...loop through components of this physics variable
               do icomp=1,il
                  k=k+1; ivar=ivar+1
                  select case(index(k))
                  case(6)
                     nn=nn+1
!
!                 ...copy the dof
                     ZSOL_C(Iel)%coarse(nn) = NODES(nod)%zdofV(ivar,j)
case default
                  end select
!           ...end loop over components
               enddo
!        ...end loop over physics variables
            enddo
!     ...end loop over nodal dof
         enddo
!     ...update the number of components stored so far
         mvarV(i) = ivar
!  ...end loop over nodes
      enddo
!
!-----------------------------------------------------------------------
!
  400 continue
!
!..end of loop through loads
   enddo
!
end subroutine compute_local_sol_dof
