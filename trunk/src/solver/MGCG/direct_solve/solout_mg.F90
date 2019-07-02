! -----------------------------------------------------------------------
!
!    routine name       - solout
!
! -----------------------------------------------------------------------
!
!    latest revision    - Sept 2018
!
!    purpose            - an interface routine with multigrid solver,
!                         storing a particular solution in the data
!                         structure arrays
!
!    arguments :
!      in:
!              Ielc     - global coarse grid element number
!              Iel      - if macro solout: local fine grid element number
!                         o/w: -1 
!              Ndof     - number of dof per element (unused)
!              Nrhs     - number of loads
!              Zele     - the element dof
!
! -----------------------------------------------------------------------
#include "implicit_none.h"
subroutine solout_mg(Ielc,Iel,Mdle,Ndof,Nrhs,Zele)
!
   use element_data
   use data_structure3D
   use frsolmod
   use assembly
   use stc, only: stc_bwd_wrapper
!
   implicit none
!
   integer, intent(in) :: Ielc
   integer, intent(in) :: Iel
   integer, intent(in) :: Mdle
   integer, intent(in) :: Ndof
   integer, intent(in) :: Nrhs
   VTYPE  , intent(in) :: Zele(Ndof)
!
!..nodes for a modified element and the corresponding number
!  of H1,H(curl),H(div) and L2 dof
   integer, dimension(MAXNODM) :: nodm,ndofmH,ndofmE,ndofmV,ndofmQ
!
!..number of variables for each physics attribute for an element
   integer :: nrdofs(NR_PHYSA)
!
!..decoded index for a node
   integer :: index(NRINDEX)
!
!..component counters for the nodes (use in case of multiple loads)
   integer, dimension(MAXNODM) :: mvarH,mvarE,mvarV
!
!..dof counters
   integer :: nrdofm,nrdofc,nrnodm
!
!..active component counter
   integer :: nvarH,nvarE,nvarV
!
!..auxiliary variables
   integer :: nrPhysH,nrPhysE,nrPhysV,nrPhysHE,nrPhysHEV
   integer :: nrVarHE,nrVarHEV
   integer :: i,j,k,il,iphys,icomp,ivar,load,nn,nod
   VTYPE   :: zvoid
!
#if DEBUG_MODE
   integer :: iprint=0
#endif
!
!----------------------------------------------------------------------
!
!..initiate component counters
   mvarH(1:MAXNODM)=0
   mvarE(1:MAXNODM)=0
   mvarV(1:MAXNODM)=0
!
!..determine nodes of the modified element
   call celem_mg(Ielc,Iel,Mdle,1,                         &
                 nrdofs,nrdofm,nrdofc,                    &
                 nodm,ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm, &
                 zvoid,zvoid)
!
#if DEBUG_MODE
      if (iprint.ge.1) then
        write(*,7001) Ielc,mdle,nodm(1:nrnodm)
 7001   format(' solout_mg: Ielc,mdle = ',i8,i10,' nodm = ',10(/,10i10))
        write(*,7002) ndofmH(1:nrnodm)
 7002   format(' solout_mg: ndofmH = ',30i3)
        write(*,7003) ndofmE(1:nrnodm)
 7003   format(' solout_mg: ndofmE = ',30i3)
        write(*,7004) ndofmV(1:nrnodm)
 7004   format(' solout_mg: ndofmV = ',30i3)
        write(*,7005) ndofmQ(1:nrnodm)
 7005   format(' solout_mg: ndofmQ = ',30i3)
        call pause
      endif
#endif
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
!-----------------------------------------------------------------------
!
!..store the dof
!
!..initiate dof counter
   nn=0
!
!..initialize the number of H1,H(curl),H(div),L2 stored so far
   ivar=0
!
!..loop through right-hand sides (loads)
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
#if DEBUG_MODE
         if (iprint.eq.1) then
            write(*,7100) nod,index(1:NRINDEX)
 7100       format('solout_mg: nod,index = ',i8,2x,20i2)
         endif
#endif
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
                     NODES(nod)%zdofH(ivar,j) = Zele(nn)
#if DEBUG_MODE
                     if (iprint.eq.1) then
                        write(*,7006) nn,load,Zele(nn)
  7006                  format('solout_mg: nn,load,Zele(nn) = ',i4,i3,x,2e13.5)
                        write(*,7007) nod,j,ivar,NODES(nod)%zdofH(ivar,j)
  7007                  format('solout_mg: nod,j,ivar,NODES(nod)%zdofH(ivar,j)', &
                               ' = ',i5,i3,i3,x,2e13.5)
                     endif
#endif
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
#if DEBUG_MODE
         if (iprint.eq.1) write(*,7100) nod,index
#endif
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
                     NODES(nod)%zdofE(ivar,j) = Zele(nn)
#if DEBUG_MODE
                     if (iprint.eq.1) then
                        write(*,7006) nn,load,Zele(nn)
                        write(*,7009) nod,j,ivar,NODES(nod)%zdofE(ivar,j)
 7009                   format('solout_mg: nod,j,ivar,NODES(nod)%zdofE(ivar,j)', &
                               ' = ',i5,i3,i3,x,2e13.5)
                     endif
#endif
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
#if DEBUG_MODE
         if (iprint.eq.1) write(*,7100) nod,index
#endif
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
                     NODES(nod)%zdofV(ivar,j) = Zele(nn)
#if DEBUG_MODE
                     if (iprint.eq.1) then
                        write(*,7006) nn,load,Zele(nn)
                        write(*,7010) nod,j,ivar,NODES(nod)%zdofV(ivar,j)
 7010                   format('solout_mg: nod,j,ivar,NODES(nod)%zdofV(ivar,j)', &
                               ' = ',i5,i3,i3,x,2e13.5)
                     endif
#endif
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
#if DEBUG_MODE
      if (iprint.eq.1) call pause
#endif
!
!..end of loop through loads
   enddo
!
!..obtain bubble solution
   if (Iel .ge. 0) then
!  ...mg macro grid
      call stc_bwd_wrapper_mg(Ielc,Iel)
   else
!  ...coarse grid
      call stc_bwd_wrapper(Ielc)
   endif
!
end subroutine solout_mg
