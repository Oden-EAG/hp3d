c----------------------------------------------------------------------
c
c   routine name       - nodcor
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 10
c
c   purpose            - routine calculates unconstrained geometry dof
c                        for a 3D element - need revision for TETRA
c
c   arguments :
c     in:
c           Mdle       - middle node of an element
c     out:
c           Xnod       - the element unconstrained geometry dof
c
c----------------------------------------------------------------------
c
      subroutine nodcor(Mdle, Xnod)
c
      use data_structure3D
      use element_data
#include "syscom.blk"
c
ccc      common / comm_nodcor / iprint_nodcor
c
      dimension Xnod(NDIMEN,MAXbrickH)
c
c  ...element order of approximation
      dimension norder(19)
c
c  ...modified element nodes and corresponding number of dof
      dimension nodm(MAXNODM),ndofmH(MAXNODM),ndofmE(MAXNODM),
     .                        ndofmV(MAXNODM),ndofmQ(MAXNODM)
c
      dimension nrconH(MAXbrickH),nrconE(MAXbrickE),nrconV(MAXbrickV),
     .          nacH(NACDIM,MAXbrickH),constrH(NACDIM,MAXbrickH),
     .          nacE(NACDIM,MAXbrickE),constrE(NACDIM,MAXbrickE),
     .          nacV(NACDIM,MAXbrickV),constrV(NACDIM,MAXbrickV)
c
c  ...modified element dof
      dimension val(NDIMEN,2*MAXbrickH)
c
c  ...element nodes
      dimension nodesl(27), norientl(27)
c
c---------------------------------------------------------------------
c
ccc      if (Mdle.eq.1194) then
ccc        iprint=iprint_nodcor
ccc      else
ccc        iprint=iprint_nodcor
ccc      endif
       if (Mdle.eq.1) then
         iprint=0
       else
         iprint=0
       endif
c
c  ...determine order of approximation
      call find_order(Mdle, norder)
c
      if (iprint.eq.1) then
        write(*,7006) Mdle
 7006   format('nodcor: Mdle   = ',i6)
        call print_order(NODES(Mdle)%type,norder)
      endif
c
c  ...determine number of local dof
      call celndof(NODES(Mdle)%type,norder,
     .             nrdoflH,nrdoflE,nrdoflV,nrdoflQ)
c
c  ...determine constraints' coefficients
      call logic(Mdle,2,
     .           nodm,ndofmH,ndofmE,ndofmV,nrnodm,
     .           nrconH,nacH,constrH,
     .           nrconE,nacE,constrE,
     .           nrconV,nacV,constrV)
c
c
      if (iprint.eq.1) then
        write(*,*) 'nodcor: nodm,ndofmH = '
        write(*,7004) (nodm(i),i=1,nrnodm)
        write(*,7004) (ndofmH(i),i=1,nrnodm)
 7004   format(20i5)
        call pause
      endif
c
c---------------------------------------------------------------------
c
c  ...copy the global dof into the local array........
c
c  ...initiate the local dof index
      k=0
c
c  ...loop through nodes
      do j=1,nrnodm
        do i=1,ndofmH(j)
          k=k+1
          if (k.gt.2*MAXbrickH) then
            write(*,*) 'nodcor: DIMENSIONS EXCEEDED !!'
            stop 1
          endif
          val(1:NDIMEN,k) = NODES(nodm(j))%coord(1:NDIMEN,i)
        enddo
      enddo
c
c---------------------------------------------------------------------
c
c  ...calculate the local dof
c
c  ...initiate the dof
      Xnod(1:NDIMEN,1:nrdoflH)=0.d0
c
c  ...loop through the local dof
      do k=1,nrdoflH
        if (iprint.eq.1) then
          write(*,*) 'k = ',k
          write(*,7001) (nacH(l,k),l=1,nrconH(k))
 7001     format(' nacH   (:,k) = ',10i5)
          write(*,7002) (constrH(l,k),l=1,nrconH(k))
 7002     format(' constrH(:,k) = ',10f5.2)
          do ivar=1,NDIMEN
            write(*,7005)ivar,(val(ivar,nacH(l,k)),l=1,nrconH(k))
 7005       format(' ivar,val(ivar,:) = ',i1,3x,10(e12.5,2x))
          enddo
        endif
c
c  .....accumulate for the values
        do kp=1,nrconH(k)
          l = nacH(kp,k)
          do ivar=1,NDIMEN
            Xnod(ivar,k) = Xnod(ivar,k) + constrH(kp,k)*val(ivar,l)
          enddo
        enddo
        if (iprint.eq.1) then
           write(*,7003) (Xnod(ivar,k),ivar=1,NDIMEN)
 7003      format('Xnod = ',3f10.4)
           call pause
        endif
c
c  ...end of loop through local dof
      enddo
c
c  ...printing
      if (iprint.eq.2) then
        call elem_nodes(Mdle, nodesl, norientl)
        write(*,*) 'nodcor: VERTEX dof'
        inod=0; k=0
        do i=1,nvert(NODES(Mdle)%Type)
          inod=inod+1; k=k+1
          write(*,7003) (Xnod(ivar,k),ivar=1,NDIMEN)
        enddo
        write(*,*) 'nodcor: EDGE dof'
        do i=1,nedge(NODES(Mdle)%Type)
          inod=inod+1; nod = nodesl(inod)
          call ndof_nod(NODES(nod)%Type,NODES(nod)%order,
     .                  ndofH,ndofE,ndofV,ndofQ)
          write(*,*) 'EDGE = ',i
          do j=1,ndofH
            k=k+1
            write(*,7003) (Xnod(ivar,k),ivar=1,NDIMEN)
          enddo
        enddo
        write(*,*) 'nodcor: FACE dof'
        do i=1,nface(NODES(Mdle)%Type)
          inod=inod+1; nod = nodesl(inod)
          call ndof_nod(NODES(nod)%Type,NODES(nod)%order,
     .                  ndofH,ndofE,ndofV,ndofQ)
          write(*,*) 'FACE = ',i
          do j=1,ndofH
            k=k+1
            write(*,7003) (Xnod(ivar,k),ivar=1,NDIMEN)
          enddo
        enddo
        write(*,*) 'nodcor: MIDDLE NODE dof'
        call ndof_nod(NODES(Mdle)%Type,NODES(nod)%order,
     .                ndofH,ndofE,ndofV,ndofQ)
        do j=1,ndofH
          k=k+1
          write(*,7003) (Xnod(ivar,k),ivar=1,NDIMEN)
        enddo
        write(*,*) 'nodcor: TOTAL NUMBER OF DOF = ',k
        call pause
      endif
c
c
      end
c
c
      subroutine test_nodcor(Mdle)

      use data_structure3D
      use control , only : EXGEOM
      implicit none
!-------------------------------------------------------------------
      ! order of approximation
      integer, dimension(19) :: norder
      integer, dimension(5)  :: nordf
      ! reference and physical coordinates
      real*8, dimension(3)   :: xi,x
      real*8, dimension(3,3) :: dxdxi,dxidx
      ! Gauss points and weights
      real*8, dimension(3,MAX_NINT3) :: xiloc
      real*8, dimension(  MAX_NINT3) :: wxi
      ! miscellanea
      real*8 :: rjac
      integer :: mdle, i, nint, l, iflag, ndom, nrdofH, k, j, ifig,
     .           iprint,nv,icheck
      character(len=4) :: ftype

      ! shape function
      real*8, dimension(  MAXbrickH) :: vshapH
      real*8, dimension(3,MAXbrickH) :: dvshapH
      integer, dimension(12)         :: nedge_orient
      integer, dimension(6)          :: nface_orient
      real*8, dimension(3,MAXbrickH) :: xnod
      !-------------------------------------------------------------------

      if (EXGEOM.ne.0) then
        write(*,5000)Mdle
 5000   format(' test_nodcor: EXACT geometry element! Mdle = ',i10)
        return
      endif

      call find_domain(mdle, ndom)
      write(*,7000) mdle,NODES(mdle)%type,ndom
 7000 format(' test_nodcor: Mdle,type,ndom = ',i8,2x,a4,2x,i2)

      !  ...order, orientations, gdof's
      call find_order( Mdle, norder)
      call find_orient(Mdle, nedge_orient,nface_orient)
      call nodcor(     Mdle, xnod)

      !  ...integration points
      call set_3Dint(NODES(mdle)%type,norder, nint,xiloc,wxi)

      !  ...number of element vertices, need to perform incremental check
      nv=nvert(NODES(mdle)%type)

      !  ...loop over integration points
      do l=1,nint

        !  ...integration point
        xi(1:3) = xiloc(1:3,l)
        write(*,7001)xi(1:3)
 7001   format(' xi = ',3(e12.5,2x))

        !  ...shape functions
        call shape3H(NODES(Mdle)%type,xi,norder,nedge_orient,
     .                nface_orient, nrdofH,vshapH,dvshapH)

        !  ...accumulate
        x(1:3)=0.d0 ; dxdxi(1:3,1:3)=0.d0
        do k=1,nrdofH
          x(1:3) = x(1:3) + xnod(1:3,k)*vshapH(k)
          do j=1,3
            dxdxi(1:3,j) = dxdxi(1:3,j) + xnod(1:3,k)*dvshapH(j,k)
          enddo

          !  ...incremental check of godf's and jacobian
          rjac=0.d0 ; if (k.ge.nv) call geom(dxdxi, dxidx,rjac,iflag)
          write(*,6999)k,xnod(1:3,k),rjac
 6999     format(' k,xnod(1:3,k),rjac = ',i3,2x,3(e12.5,2x),2x,e12.5)
        enddo

      !  ...loop over integration points
      enddo


      end subroutine test_nodcor
