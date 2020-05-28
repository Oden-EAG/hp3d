c-----------------------------------------------------------------------
c
c   routine name       - solout
c
c-----------------------------------------------------------------------
c
c   latest revision    - Sept 2018
c
c   purpose            - an interface routine with frontal solver,
c                        storing a particular solution in the data
c                        structure arrays
c
c   arguments :
c     in:
c             Iel      - element number
c             Ndof     - number of dof per element (unused)
c             Nrhs     - number of loads
c             Mdest    - the element destination vectors (unused)
c             Zele     - the element dof
c
c   required  routines -
c
c-----------------------------------------------------------------------
c
c
      subroutine solout(Iel,Ndof,Nrhs,Mdest,Zele)
c
      use element_data
      use data_structure3D
      use frsolmod
      use assembly
      use control, only: ISTC_FLAG
      use stc    , only: stc_bwd_wrapper
#include "syscom.blk"
c
      dimension Zele(*)
c
c  ...nodes for a modified element and the corresponding number
c     of H1,H(curl),H(div) and L2 dof
      dimension nodm(MAXNODM),ndofmH(MAXNODM),ndofmE(MAXNODM),
     .                        ndofmV(MAXNODM),ndofmQ(MAXNODM)
c
c  ...number of variables for each physics attribute for an element
      dimension nrdofs(NR_PHYSA)
c
c  ...decoded index for a node
      dimension index(NRINDEX)
c
c  ...component counters for the nodes (use in case of multiple loads)
      dimension mvarH(MAXNODM),mvarE(MAXNODM),
     .          mvarV(MAXNODM),mvarQ(MAXNODM)
c
c----------------------------------------------------------------------
c
      iprint=0
c
c  ...initiate component counters
      mvarH = 0; mvarE = 0; mvarV = 0; mvarQ = 0
c
      if (REORDER) then
c
c  .....find the element number
        mdle = NEW_ELEM_ORDER(Iel)
      else
c
c  .....find the element number using the standard ordering of elements
        mdle=0
        do i=1,Iel
          call nelcon(mdle, mdle)
        enddo
      endif
c
c  ...determine nodes of the modified element
      if (ISTC_FLAG) then
        call celem_systemI(Iel,mdle,1,
     .                     nrdofs,nrdofm,nrdofc,
     .                     nodm,ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,
     .                     zvoid,zvoid)

      else
        call celem(mdle,1,
     .             nrdofs,nrdofm,nrdofc,
     .             nodm,ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,
     .             zvoid,zvoid)
      endif

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
c
c-----------------------------------------------------------------------
c
c  ...store the dof
c
c  ...loop through right-hand sides (loads).............................
c
c  ...initiate the frontal solver dof counter
      nn=0
c
c  ...initialize the number of H1,H(curl),H(div),L2 stored so far
      ivar=0
c
c  ...loop through loads
      do load=1,NR_RHS
c
c  .....H1 dof .................................
        if (NRHVAR.eq.0) goto 200
c
c  .....loop through nodes of modified element in reversed order
        do i=nrnodm,1,-1
          nod = nodm(i)
c
c  .......compute the number of active H1 variables for the node
          call get_index(nod, index)
          if (iprint.eq.1) then
            write(*,7020) nod,index(1:NRINDEX)
 7020       format('solout: nod,index = ',i8,2x,10i2)
          endif
          nvarH=0
          do k=1,NRINDEX
            if (index(k).eq.2) nvarH=nvarH+1
          enddo
          if (nvarH.eq.0) cycle
c
c  .......loop through the nodal dof (potentially NONE)
          do j=1,ndofmH(i)/nvarH
c
c  .........loop through the components
            ivar=mvarH(i)
            do k=1,NRINDEX
              select case(index(k))
              case(1)
                ivar=ivar+1
              case(2)
                ivar=ivar+1
                nn=nn+1
c
c  .............copy the dof
                if (iprint.eq.1) then
                  write(*,7021) load,i,ivar,j,nn
 7021             format('solout: load,i,ivar,j,nn = ',5i5)
                endif
                NODES(nod)%zdofH(ivar,j) = Zele(nn)
                if (iprint.eq.1) then
                  write(*,7006) nn,Zele(nn)
 7006             format('solout: nn, Zele(nn) = ',i4,2e12.5)
                  write(*,7007) nod,j,ivar,NODES(nod)%zdofH(ivar,j)
 7007             format('solout: nod,j,ivar,NODES(nod)%zdofH(ivar,j)',
     .                   ' = ',i5,i3,i3,2x,2e12.5)
                endif
              end select
            enddo
          enddo
c
c  .......update the number of components stored so far
          mvarH(i) = ivar
        enddo
        if (iprint.eq.1) write(*,*) 'solout: H1 solout DONE'
c
c-----------------------------------------------------------------------
c
c  .....H(curl) dof .................................
  200   if (NREVAR.eq.0) go to 300
c
c  .....loop through nodes in the reversed order
        do i=nrnodm,1,-1
          nod = nodm(i)
c
c  .......compute the number of active H(curl) variables for the node
          call get_index(nod, index)
          if (iprint.eq.1) then
            write(*,7100) nod,index
 7100       format('solout: nod = ',i5,' index = ',10i2)
          endif
          nvarE=0
          do k=1,NRINDEX
            if (index(k).eq.4) nvarE=nvarE+1
          enddo
          if (nvarE.eq.0) cycle
c
c  .......loop through the nodal dof
          do j=1,ndofmE(i)/nvarE
c
c  .........loop through the components
            ivar=mvarE(i)
            do k=1,NRINDEX
              select case(index(k))
              case(3)
                ivar=ivar+1
              case(4)
                ivar=ivar+1
                nn=nn+1
c
c  .............copy the dof
                NODES(nod)%zdofE(ivar,j) = Zele(nn)
                if (iprint.eq.1) then
                  write(*,7006) nn,Zele(nn)
                  write(*,7009) nod,j,ivar,NODES(nod)%zdofE(ivar,j)
 7009             format('solout: nod,j,ivar,NODES(nod)%zdofE(ivar,j)',
     .                   ' = ',i5,i3,i3,2x,2e12.5)
                endif
              end select
            enddo
          enddo
c
c  .......update the number of components stored so far
          mvarE(i) = ivar
        enddo
c
c-----------------------------------------------------------------------
c
c  .....H(div) dof .................................
  300   if (NRVVAR.eq.0) go to 400
c
c  .....loop through nodes in the reversed order
        do i=nrnodm,1,-1
          nod = nodm(i)
c
c  .......compute the number of active H(div) variables for the node
          call get_index(nod, index)
          nvarV=0
          do k=1,NRINDEX
            if (index(k).eq.6) nvarV=nvarV+1
          enddo
          if (nvarV.eq.0) cycle
c
c  .......loop through the nodal dof
          do j=1,ndofmV(i)/nvarV
c
c  .........loop through the components
            ivar=mvarV(i)
            do k=1,NRINDEX
              select case(index(k))
              case(5)
                ivar=ivar+1
              case(6)
                ivar=ivar+1
                nn=nn+1
c
c  .............copy the dof
                NODES(nod)%zdofV(ivar,j) = Zele(nn)
                if (iprint.eq.1) then
                  write(*,7006) nn,Zele(nn)
                  write(*,7010) nod,j,ivar,NODES(nod)%zdofV(ivar,j)
 7010             format('solout: nod,j,ivar,NODES(nod)%zdofV(ivar,j)',
     .                   ' = ',i5,i3,i3,2x,2e12.5)
                endif
              end select
            enddo
          enddo
c
c  .......update the number of components stored so far
          mvarV(i) = ivar
        enddo
c
c-----------------------------------------------------------------------
c
c  .....L2 dof .................................
  400   if (NRQVAR.eq.0) go to 500
        if (ISTC_FLAG) go to 500
c
c  .....middle node only
        i=nrnodm
        nod = nodm(i)
c
c  .....compute the number of active L2 variables for the node
        call get_index(nod, index)
        nvarQ=0
        do k=1,NRINDEX
          if (index(k).eq.8) nvarQ=nvarQ+1
        enddo
        if (nvarQ.eq.0) go to 500
c
c  .....loop through the nodal dof
        do j=1,ndofmQ(i)/nvarQ
c
c  .......loop through the components
          ivar=mvarQ(i)
          do k=1,NRINDEX
            select case(index(k))
            case(7)
              ivar=ivar+1
            case(8)
              ivar=ivar+1
              nn=nn+1
c
c  ...........copy the dof
              NODES(nod)%zdofQ(ivar,j) = Zele(nn)
              if (iprint.eq.1) then
                write(*,7006) nn,Zele(nn)
                write(*,7011) nod,j,ivar,NODES(nod)%zdofQ(ivar,j)
 7011           format('solout: nod,j,ivar,NODES(nod)%zdofQ(ivar,j)',
     .                 ' = ',i5,i3,i3,2x,2e12.5)
              endif
            end select
          enddo
        enddo
c
c  .....update the number of components stored so far
        mvarQ(i) = ivar
 500    continue
        if (iprint.eq.1) call pause
c
c  ...end of loop through loads
      enddo
c
c  ...obtain bubble solution
      if (ISTC_FLAG) then
         call stc_bwd_wrapper(Iel)
      endif
c
      end subroutine solout


