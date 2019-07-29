c----------------------------------------------------------------------
c
c     routine name       - hpmdle
c
c----------------------------------------------------------------------
c
c     latest revision    - Oct 11
c
c     purpose            - routine determines geometry dof for a middle
c     node of an element by performing H_0^1
c     projection in the GMP REFERENCE space
c
c     arguments :
c     in:
c     Iflag        - a flag specifying which of the objects the
c     mid-face node is generated on
c     = 5  for a prism
c     = 6  for a hexahedron
c     = 7  for a tetrahedron
c     = 8  for a pyramid
c     No           - number of a specific object (prism,hexa...)
c     Etav         - reference coordinates of the element vertices
c     Type         - element (middle node) type
c     Nedge_orient - edge orientation
c     Nface_orient - face orientation
c     Norder       - element order
c     Xnod         - element geometry dof so far
c     out:
c     Xdof         - geometry dof for the middle node
c
c----------------------------------------------------------------------
c
      subroutine hpmdle(Iflag,No,Etav,
     .  Type,Nedge_orient,Nface_orient,Norder,
     .  Xnod,Xdof)
c
      use element_data
      use parameters
c
#include "syscom.blk"
      common /cshape3H/ iprint_shape3H
c
      dimension Etav(3,8),Xnod(3,MAXbrickH),Xdof(3,*)
      character(len=4) :: Type
      dimension Nedge_orient(12),Nface_orient(6),Norder(19)
c
c     ...local element order
      dimension nordloc(19)
c
c     ...element shape functions
      dimension shapH(MAXbrickH),dshapH(3,MAXbrickH)
c
c     ...geometry
      dimension xi(3),eta(3),x(3),
     .  dxdt(3),dxidt(3,2),
     .  dxdxi(3,3),dxdeta(3,3),
     .  detadxi(3,3),dxideta(3,3),gg(3,3),svec(3)
c
c     ...quadrature data
      dimension xiloc(3,MAX_NINT3),waloc(MAX_NINT3)
c
c     ...linear system matrices
      dimension aa(MAXmdlbH,MAXmdlbH),bb(MAXmdlbH,3)
      dimension ipiv(MAXmdlbH)
c
c----------------------------------------------------------------------
c
c      if (Type.eq.'mdld') then
c        iprint=0
c      else
c        iprint=0
c      endif
      iprint = 0
ccc   iprint_shape3H=iprint
      if (iprint.eq.1) then
        write(*,7001) Norder(1:19)
 7001   format('hpmdle: Norder = ',19i4)
      endif
c
c     ...set up the quadrature (dependent upon element type)
      call set_3Dint(Type,Norder, nint,xiloc,waloc)
c
      select case(Type)
      case('mdlp'); iv=15
      case('mdlb'); iv=19
      case('mdln'); iv=11
      case('mdld'); iv=14
      end select
c
c     ...number of dof
      nord = Norder(iv)
      call ndof_nod(Type,nord, ndofH,ndofE,ndofV,ndofQ)
      if (ndofH.eq.0) return
c
c     ...initiate matrices
      aa(1:ndofH,1:ndofH) = 0.d0
      bb(1:ndofH,1:3) = 0.d0
c
c     ...nordloc is a copy of Norder, except for the middle node's order
c     which is set to 1
      call initiate_order(Type,nordloc)
      nef = nedge(Type)+nface(Type)
      nordloc(1:nef) = Norder(1:nef)
c
c     ...compute the number of vertex,edge and face dof
      call celndof(Type,nordloc, nrdofH,nrdofE,nrdofV,nrdofQ)
      if (iprint.ge.1) then
        write(*,7003) nordloc(1:iv)
 7003   format('hpmdle: nordloc = ',19i4)
        write(*,7004) nrdofH
 7004   format('hpmdle: nrdofH = ',i4)
      endif
c
c     ...loop through integration points
      do l=1,nint
c
c     .....compute ALL 3D element shape functions
        call shape3H(Type,xiloc(1:3,l),Norder,
     .    Nedge_orient,Nface_orient,
     .    nvoid,shapH,dshapH)
c
c     .....evaluate reference coordinates of the point
        eta(1:3) = 0.d0; detadxi(1:3,1:3)=0.d0
        do k=1,nvert(Type)
          eta(1:3) = eta(1:3) + etav(1:3,k)*shapH(k)
          do i=1,3
            detadxi(1:3,i) = detadxi(1:3,i) + etav(1:3,k)*dshapH(i,k)
          enddo
        enddo
        if (iprint.eq.1) then
          do ivar=1,3
            write(*,7012) ivar,detadxi(ivar,1:3)
 7012       format(' ivar,detadxi(ivar,1:3) = ',i2,3e12.5)
          enddo
          call pause
        endif
c
c     .....compute metric
        call geom(detadxi, dxideta,rjac,iflag1)
        if (iflag1.ne.0) then
          write(*,*) 'hpmdle: iflag1 = ',iflag1; stop 1
        endif
        do i=1,3
          do j=1,3
            s = 0.d0
            do k=1,3
              s = s + dxideta(i,k)*dxideta(j,k)
            enddo
            gg(i,j) = s
          enddo
        enddo
        weight = waloc(l)*rjac
c
c     .....evaluate -dxdxi obtained with vertex, edge and face dof
        dxdxi(1:3,1:3)=0.d0
        do k=1,nrdofH
          do i=1,3
            dxdxi(1:3,i) = dxdxi(1:3,i) - Xnod(1:3,k)*dshapH(i,k)
          enddo
        enddo
        if (iprint.eq.1) then
          do ivar=1,3
            write(*,7011) ivar,dxdxi(ivar,1:3)
          enddo
          call pause
        endif
c
c     .....evaluate derivatives of the physical coordinates wrt
c     the reference coordinate
        select case(Iflag)
      case(5);          call prism(No,eta, x,dxdeta)
      case(6);          call  hexa(No,eta, x,dxdeta)
      case(7);          call tetra(No,eta, x,dxdeta)
      case(8);          call pyram(No,eta, x,dxdeta)
      case default
        write(*,7002) Type,Iflag
 7002   format('hpmdle: Type,Iflag = ',a5,i4)
        stop 1
      end select
c
      if (iprint.eq.1) then
        do ivar=1,3
          write(*,7010) ivar,dxdeta(ivar,1:3)
 7010     format('ivar,dxdeta(ivar,1:3) = ',i2,3e12.5)
        enddo
        call pause
      endif
c
c     .....evaluate derivatives wrt master coordinates
      do i=1,3
        do j=1,3
          dxdxi(1:3,i) = dxdxi(1:3,i) + dxdeta(1:3,j)*detadxi(j,i)
        enddo
      enddo
      if (iprint.eq.1) then
        do ivar=1,3
          write(*,7011) ivar,dxdxi(ivar,1:3)
 7011     format('ivar,dxdxi(ivar,1:3) = ',i2,3e12.5)
        enddo
        call pause
      endif
c
c
c     .....first loop through dof
      do j=1,ndofH
c
c     .......accumulate for the load vector
        svec = 0.d0
        do k1=1,3; do k2=1,3
          svec(1:3) = svec(1:3)
     .      + dxdxi(1:3,k1)*dshapH(k2,nrdofH+j)*gg(k1,k2)
          enddo; enddo
          bb(j,1:3) = bb(j,1:3) + svec(1:3)*weight
c
c     .......accumulate for the stiffness matrix
          do i=1,ndofH
            s = 0.d0
            do k1=1,3; do k2=1,3
              s = s + dshapH(k1,nrdofH+i)*dshapH(k2,nrdofH+j)*gg(k1,k2)
              enddo; enddo
              aa(i,j) = aa(i,j) + s*weight
            enddo
          enddo
c
c     ...end of the loop through integration points
        enddo
c
c
        if (iprint.eq.2) then
          write(*,*) 'hpmdle: LOAD VECTOR = '
          do ii=1,ndofH
            write(*,7005) bb(ii,1:3)
 7005       format(1x,10e12.5)
          enddo
          write(*,*) 'hpmdle: STIFFNESS MATRIX = '
          do ii=1,ndofH
            write(*,7005) aa(ii,1:ndofH)
          enddo
          call pause
        endif
c
c     ...solve the linear system...
        call DGESV( ndofH, 3, aa, MAXmdlbH, ipiv, bb, MAXmdlbH, info )
        do i=1,3
          Xdof(i,1:ndofH) = bb(1:ndofH,i)
        enddo
c
        if (iprint.eq.2) then
          write(*,*) 'hpmdle: Xdof = '
          do i=1,ndofH
            write(*,7005) Xdof(1:3,i)
          enddo
          call pause
        endif
c
c
        endsubroutine hpmdle
