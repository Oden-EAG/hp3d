!-----------------------------------------------------------------------
!
!   routine name       - solout1
!
!-----------------------------------------------------------------------
!
!   latest revision    - May 10
!
!   purpose            - an interface routine with frontal solver,
!                        storing a particular solution in the data
!                        structure arrays
!
!   arguments :
!     in:
!             Mdle     - element number
!             Ndof     - number of dof per element
!             Nrhs     - number of loads
!             Zele     - the element dof
!
!   required  routines -
!
!-----------------------------------------------------------------------
#include "typedefs.h"
   subroutine solout1(Mdle,Ndof,Nrhs,Mdest,Zele)
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
   integer, intent(in) :: Mdle
   integer, intent(in) :: Ndof
   integer, intent(in) :: Nrhs
   VTYPE  , intent(in) :: Zele(Ndof)
!
!..nodes for a modified element and the corresponding number
!  of H1,H(curl),H(div) and L2 dof
   integer, dimension(MAXNODM) :: nodm,ndofmH,ndofmE,   &
             ndofmV,ndofmQ
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

!..auxiliary variables
   integer :: nrPhysH,nrPhysE,nrPhysV,nrPhysHE,nrPhysHEV
   integer :: nrVarHE,nrVarHEV
   integer :: i,j,k,il,iphys,icomp,ivar,load,mdle,nn,nod
   VTYPE   :: zvoid


!
!----------------------------------------------------------------------
!
   iprint=0
!
!..initiate component counters
   mvarH = 0; mvarE = 0; mvarV = 0; mvarQ = 0
!
!..determine nodes of the modified element
   call celem1(Mdle,1, nrdofs,nrdofm,nrdofc,     &
               nodm,ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm, &
               zvoid,zvoid)
      if (iprint.ge.1) then
        write(*,7001) mdle,nodm(1:nrnodm)
 7001   format(' solout: mdle = ',i8,i10,' nodm = ',10(/,10i10))
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
!
!-----------------------------------------------------------------------
!
!..store the dof
!
!..loop through right-hand sides (loads).............................
!
!..initiate the frontal solver dof counter
   nn=0
!
!..initialize the number of H1,H(curl),H(div),L2 stored so far
   ivar=0
!
!..loop through loads
   do load=1,NR_RHS
!
!  ...H1 dof .................................
      if (NRHVAR.eq.0) goto 200
!
!  ...loop through nodes of modified element
      do i=1,nrnodm
         nod = nodm(i)
!
!     ...compute the number of active H1 variables for the node
         call get_index(nod, index)
         if (iprint.eq.1) then
            write(*,7020) nod,index(1:NRINDEX)
 7020       format('solout1: nod,index = ',i8,2x,10i2)
         endif
         nvarH=0
         do k=1,NRINDEX
            if (index(k).eq.2) nvarH=nvarH+1
         enddo
         if (nvarH.eq.0) cycle
!
!     ...loop through the nodal dof (potentially NONE)
         do j=1,ndofmH(i)/nvarH
!
!        ...loop through the components
            ivar=mvarH(i)
            do k=1,NRINDEX
               select case(index(k))
            case(1)
               ivar=ivar+1
            case(2)
               ivar=ivar+1
               nn=nn+1
!
!           ...copy the dof
               if (iprint.eq.1) then
                  write(*,7021) load,i,ivar,j,nn
 7021             format('solout1: load,i,ivar,j,nn = ',5i5)
               endif
               NODES(nod)%zdofH(ivar,j) = Zele(nn)
               if (iprint.eq.1) then
                  write(*,7006) nn,Zele(nn)
 7006             format('solout1: nn, Zele(nn) = ',i4,2e12.5)
                  write(*,7007) nod,j,ivar,NODES(nod)%zdofH(ivar,j)
 7007             format('solout1: nod,j,ivar,NODES(nod)%zdofH(ivar,1)',   &
                         ' = ',i5,i3,i3,2x,2e12.5)
               endif
            end select
         enddo
      enddo
!
!  ...update the number of components stored so far
      mvarH(i) = ivar
   enddo
   if (iprint.eq.1) write(*,*) 'solout1: H1 solout DONE'
!
!-----------------------------------------------------------------------
!
!..H(curl) dof .................................
 200 if (NREVAR.eq.0) go to 300
!
!..loop through nodes
   do i=1,nrnodm
      nod = nodm(i)
!
!  ...compute the number of active H(curl) variables for the node
      call get_index(nod, index)
      if (iprint.eq.1) then
         write(*,7100) nod,index
 7100    format('solout: nod = ',i5,' index = ',10i2)
      endif
      nvarE=0
      do k=1,NRINDEX
         if (index(k).eq.4) nvarE=nvarE+1
      enddo
      if (nvarE.eq.0) cycle
!
!  ...loop through the nodal dof
      do j=1,ndofmE(i)/nvarE
!
!     ...loop through the components
         ivar=mvarE(i)
         do k=1,NRINDEX
            select case(index(k))
            case(3)
               ivar=ivar+1
            case(4)
               ivar=ivar+1
               nn=nn+1
!
!           ...copy the dof
               NODES(nod)%zdofE(ivar,j) = Zele(nn)
               if (iprint.eq.1) then
                  write(*,7006) nn,Zele(nn)
                  write(*,7009) nod,j,ivar,NODES(nod)%zdofE(ivar,j)
 7009             format('solout: nod,j,ivar,NODES(nod)%zdofE(ivar,1)',   &
                         ' = ',i5,i3,i3,2x,2e12.5)
               endif
            end select
         enddo
      enddo
!
!  ...update the number of components stored so far
      mvarE(i) = ivar
   enddo
!
!-----------------------------------------------------------------------
!
!..H(div) dof .................................
 300 if (NRVVAR.eq.0) go to 400
!
!..loop through nodes
   do i=1,nrnodm
      nod = nodm(i)
!
!  ...compute the number of active H(div) variables for the node
      call get_index(nod, index)
      nvarV=0
      do k=1,NRINDEX
         if (index(k).eq.6) nvarV=nvarV+1
      enddo
      if (nvarV.eq.0) cycle
!
!     ...loop through the nodal dof
         do j=1,ndofmV(i)/nvarV
!
!        ...loop through the components
            ivar=mvarV(i)
            do k=1,NRINDEX
               select case(index(k))
               case(5)
                  ivar=ivar+1
               case(6)
                  ivar=ivar+1
                  nn=nn+1
!
!              ...copy the dof
                  NODES(nod)%zdofV(ivar,j) = Zele(nn)
                  if (iprint.eq.1) then
                     write(*,7006) nn,Zele(nn)
                     write(*,7010) nod,j,ivar,NODES(nod)%zdofV(ivar,j)
 7010                format('solout: nod,j,ivar,NODES(nod)%zdofV(ivar,1)',    &
                            ' = ',i5,i3,i3,2x,2e12.5)
                  endif
               end select
            enddo
         enddo
!
!  ...update the number of components stored so far
      mvarV(i) = ivar
   enddo
!
!-----------------------------------------------------------------------
!
!..L2 dof
  400 if (NRQVAR.eq.0) go to 500
!
!..middle node only
   i=nrnodm
   nod = nodm(i)
!
!..compute the number of active L2 variables for the node
   call get_index(nod, index)
   nvarQ=0
   do k=1,NRINDEX
      if (index(k).eq.8) nvarQ=nvarQ+1
   enddo
   if (nvarQ.eq.0) go to 500
!
!  ...loop through the nodal dof
      do j=1,ndofmQ(i)/nvarQ
!
!     ...loop through the components
         ivar=mvarQ(i)
         do k=1,NRINDEX
            select case(index(k))
            case(7)
               ivar=ivar+1
            case(8)
               ivar=ivar+1
               nn=nn+1
!
!           ...copy the dof
               NODES(nod)%zdofQ(ivar,j) = Zele(nn)
               if (iprint.eq.1) then
                  write(*,7006) nn,Zele(nn)
                  write(*,7011) nod,j,ivar,NODES(nod)%zdofQ(ivar,j)
 7011          format('solout: nod,j,ivar,NODES(nod)%zdofQ(ivar,1)',     &
                      ' = ',i5,i3,i3,2x,2e12.5)
            endif
         end select
      enddo
   enddo
!
!..update the number of components stored so far
   mvarQ(i) = ivar
 500 continue
   if (iprint.eq.1) call pause
!
!..end of loop through loads
   enddo
!
!
   end subroutine solout1







