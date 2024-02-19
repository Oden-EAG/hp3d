! -----------------------------------------------------------------------
!
!    routine name       - solout
!
! -----------------------------------------------------------------------
!
!    latest revision    - Sept 2018
!
!    purpose            - an interface routine with frontal solver,
!                         storing a particular solution in the data
!                         structure arrays
!
!    arguments :
!      in:
!              Iel      - element number
!              Ndof     - number of dof per element (unused)
!              Nrhs     - number of loads
!              Mdest    - the element destination vectors (unused)
!              Zele     - the element dof
!
! -----------------------------------------------------------------------
#include "typedefs.h"
subroutine solout(Iel,Ndof,Nrhs,Mdest,Zele)
!
   use element_data
   use data_structure3D
   use frsolmod
   use assembly
   use control,  only: ISTC_FLAG
   use par_mesh, only: DISTRIBUTED
   use stc,      only: stc_bwd_wrapper
!
   implicit none
!
   integer, intent(in) :: Iel
   integer, intent(in) :: Ndof
   integer, intent(in) :: Nrhs
   integer, intent(in) :: Mdest
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
   integer                     :: mvarQ
!
!..dof counters
   integer :: nrdofm,nrdofc,nrnodm
!
!..active component counter
   integer :: nvarH,nvarE,nvarV,nvarQ
!
!..auxiliary variables
   integer :: nrPhysH,nrPhysE,nrPhysV,nrPhysHE,nrPhysHEV
   integer :: nrVarHE,nrVarHEV
   integer :: i,j,k,il,iphys,icomp,ivar,load,mdle,nn,nod
   VTYPE   :: zvoid
!
#if DEBUG_MODE
   integer :: iprint
   iprint=0
#endif
!
!----------------------------------------------------------------------
!
!..initiate component counters
   mvarH(1:MAXNODM)=0; mvarE(1:MAXNODM)=0
   mvarV(1:MAXNODM)=0; mvarQ           =0
!
   if (REORDER) then
!  ...find the element number
      mdle = NEW_ELEM_ORDER(Iel)
   else
!  ...find the element number using the standard ordering of elements
      if (DISTRIBUTED) then
         mdle = ELEM_SUBD(Iel)
      else
         mdle = ELEM_ORDER(Iel)
      endif
   endif
!
!..determine nodes of the modified element
   if (ISTC_FLAG) then
      call celem_systemI(Iel,mdle,1,                              &
                         nrdofs,nrdofm,nrdofc,                    &
                         nodm,ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm, &
                         zvoid,zvoid)

   else
      call celem(mdle,1,                                          &
                 nrdofs,nrdofm,nrdofc,                            &
                 nodm,ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,         &
                 zvoid,zvoid)
   endif
!
#if DEBUG_MODE
      if (iprint.ge.1) then
        write(*,7001) Iel,mdle,nodm(1:nrnodm)
 7001   format(' solout: Iel,mdle = ',i8,i10,' nodm = ',10(/,10i10))
        write(*,7002) ndofmH(1:nrnodm)
 7002   format(' solout: ndofmH = ',30i3)
        write(*,7003) ndofmE(1:nrnodm)
 7003   format(' solout: ndofmE = ',30i3)
        write(*,7004) ndofmV(1:nrnodm)
 7004   format(' solout: ndofmV = ',30i3)
        write(*,7005) ndofmQ(1:nrnodm)
 7005   format(' solout: ndofmQ = ',30i3)
        call pause
      endif
#endif
!
!..count number of variables for each physics type
   nrPhysH=0; nrPhysE=0; nrPhysV=0
   do iphys=1,NR_PHYSA
      select case(D_TYPE(iphys))
         case(CONTIN)
            nrPhysH=nrPhysH+1
         case(TANGEN)
            nrPhysE=nrPhysE+1
         case(NORMAL)
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
!..initiate the frontal solver dof counter
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
!  ...loop through nodes of modified element in reversed order
      do i=nrnodm,1,-1
         nod = nodm(i)
!
!     ...compute the number of active H1 variables for the node
         call get_index(nod, index)
#if DEBUG_MODE
         if (iprint.eq.1) then
            write(*,7100) nod,index(1:NRINDEX)
 7100       format('solout: nod,index = ',i8,2x,20i2)
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
                     NODES(nod)%dof%zdofH(ivar,j,N_COMS) = Zele(nn)
#if DEBUG_MODE
                     if (iprint.eq.1) then
                        write(*,7006) nn,load,Zele(nn)
  7006                  format('solout: nn,load,Zele(nn) = ',i4,i3,1x,2e13.5)
                        write(*,7007)   nod,j,ivar,NODES(nod)%dof%zdofH(ivar,j,N_COMS)
  7007                  format('solout: nod,j,ivar,NODES(nod)%dof%zdofH(ivar,j,N_COMS)', &
                               ' = ',i5,i3,i3,1x,2e13.5)
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
!  ...loop through nodes of modified element in reversed order
      do i=nrnodm,1,-1
         nod = nodm(i)
!
!     ...compute the number of active H(curl) variables for the node
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
                     NODES(nod)%dof%zdofE(ivar,j,N_COMS) = Zele(nn)
#if DEBUG_MODE
                     if (iprint.eq.1) then
                        write(*,7006) nn,load,Zele(nn)
                        write(*,7009)   nod,j,ivar,NODES(nod)%dof%zdofE(ivar,j,N_COMS)
 7009                   format('solout: nod,j,ivar,NODES(nod)%dof%zdofE(ivar,j,N_COMS)', &
                               ' = ',i5,i3,i3,1x,2e13.5)
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
!  ...loop through nodes of modified element in reversed order
      do i=nrnodm,1,-1
         nod = nodm(i)
!
!     ...compute the number of active H(div) variables for the node
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
                     NODES(nod)%dof%zdofV(ivar,j,N_COMS) = Zele(nn)
#if DEBUG_MODE
                     if (iprint.eq.1) then
                        write(*,7006) nn,load,Zele(nn)
                        write(*,7010)   nod,j,ivar,NODES(nod)%dof%zdofV(ivar,j,N_COMS)
 7010                   format('solout: nod,j,ivar,NODES(nod)%dof%zdofV(ivar,j,N_COMS)', &
                               ' = ',i5,i3,i3,1x,2e13.5)
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
!  .....L2 dof .................................
  400 if (NRQVAR.eq.0) goto 500
      if (ISTC_FLAG  ) goto 500
!
!  ...middle node only
      i=nrnodm
      nod=nodm(nrnodm)
!
!  ...compute the number of active L2 variables for the node
      call get_index(nod, index)
      nvarQ=0; k=nrVarHEV
      do iphys=nrPhysHEV+1,NR_PHYSA
         il = NR_COMP(iphys)
         if (.not. PHYSAm(iphys)) then
            k=k+il
            cycle
         endif
         do icomp=1,il
            k=k+1
            if (index(k).eq.8) nvarQ=nvarQ+1
         enddo
      enddo
      if (nvarQ.eq.0) goto 500
!
!  ...loop through the nodal dof
      do j=1,ndofmQ(i)/nvarQ
!
!     ...loop through the components
         ivar=mvarQ; k=nrVarHEV
!
!     ...loop through physics variables
         do iphys=nrPhysHEV+1,NR_PHYSA
            il = NR_COMP(iphys)
            if (index(k+1) .eq. 0) then
               k=k+il
               cycle
            endif
            if (.not. PHYSAm(iphys)) then
               k=k+il; ivar=ivar+il
               cycle
            endif
!        ...loop through components of this physics variable
            do icomp=1,il
               k=k+1; ivar=ivar+1
               select case(index(k))
               case(8)
                  nn=nn+1
!
!              ...copy the dof
                  NODES(nod)%dof%zdofQ(ivar,j,N_COMS) = Zele(nn)
#if DEBUG_MODE
                  if (iprint.eq.1) then
                     write(*,7006) nn,load,Zele(nn)
                     write(*,7011)   nod,j,ivar,NODES(nod)%dof%zdofQ(ivar,j,N_COMS)
 7011                format('solout: nod,j,ivar,NODES(nod)%dof%zdofQ(ivar,j,N_COMS)', &
                            ' = ',i5,i3,i3,1x,2e13.5)
                  endif
#endif
               end select
!        ...end loop over components
            enddo
!     ...end loop over physics variables
         enddo
!  ...end loop over nodal dof
      enddo
!  ...update the number of components stored so far
      mvarQ = ivar
 500  continue
!
#if DEBUG_MODE
      if (iprint.eq.1) call pause
#endif
!
!..end of loop through loads
   enddo
!
!..obtain bubble solution
   if (ISTC_FLAG) then
      call stc_bwd_wrapper(Iel)
   endif
!
end subroutine solout
