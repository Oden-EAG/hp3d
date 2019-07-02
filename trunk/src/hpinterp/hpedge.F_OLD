c----------------------------------------------------------------------
c   latest revision    - Jan 07
c
c   purpose            - routine determines geometry dof for a mid-edge
c                        node of an element
c
c   arguments :
c     in:
c         Iflag        - a flag specifying which of the objects the
c                        mid-edge node is generated on
c                        = 5  for a prism
c                        = 6  for a hexahedron
c                        = 7  for a tetrahedron
c                        = 8  for a pyramid
c         No           - number of a specific object (curve,triangle,
c                        rectangle...)
c         Etav         - reference coordinates of the element vertices
c         Type         - element (middle node) type
c         Nedge_orient - edge orientation (not needed really)
c         Nface_orient - face orientation (not needed really)
c         Norder       - element order
c         Iedg         - edge number
c         Xnod         - element gdofs so far (vertices)
c     out:
c         Xdof         - geometry dofs for the mid-edge node
c----------------------------------------------------------------------
c
      subroutine hpedge(Iflag,No,Etav,
     .                  Type,Nedge_orient,Nface_orient,Norder,Iedg,
     .                  Xnod,Xdof)
c
      use error
      use element_data
      use parameters , only : MAXbrickH,MAXP
c
#include "syscom.blk"
      common /chpedge/ iprint_hpedge
      common /cshape3H/ iprint_shape3H
      common /cshapHbe/ iprint_shapHbe
      common /cprism/ iprint_prism
      common /cprism_TI/ iprint_prism_TI
c
      dimension Etav(3,8),Xnod(3,MAXbrickH),Xdof(3,*)
      character(len=4) :: Type
      dimension Nedge_orient(12),Nface_orient(6),Norder(19)
c
c  ...local element order
      dimension nordloc(19)
c
c  ...element shape functions
      dimension vshapH(MAXbrickH),dvshapH(3,MAXbrickH)
c
c  ...edge shape functions
      dimension vshap(MAXP-1),dvshap(MAXP-1)
c
c  ...geometry
      dimension eta(3),x(3),dxdeta(3,3)
      dimension xi(3),
     .          dxdt(3),dxidt(3),
     .          dxdxi(3,3),detadxi(3,3)
c
c  ...Gauss points and weights
      dimension tloc(MAXP+1), waloc(MAXP+1)
c
c  ...linear system matrices
      dimension aa(MAXP-1,MAXP-1),bb(MAXP-1,3),uu(MAXP-1)
      dimension ipiv(MAXP-1)
c
      save nord_old,aa,ipiv
      data nord_old /1/
c
c----------------------------------------------------------------------
c
      iprint=iprint_hpedge
ccc      iprint_shape3H=iprint
ccc      iprint_shapHbe=iprint
      if (iprint.eq.1) then
        write(*,7001) Iflag,No
 7001   format('hpedge: Iflag,No = ',i2,i6)
        if (No.eq.1) then
ccc          iprint_prism_TI=1
ccc          iprint_prism=1
        else
ccc          iprint_prism_TI=0
ccc          iprint_prism=0
        endif
      endif
c
      nord = Norder(Iedg)
      if (nord.lt.2) return
c
c  ...Gauss integration points and weights
      call set_1Dint(nord, nint,tloc,waloc)
c
c  ...number of dofs
      nrdof=nord-1
c
c  ...initiate matrices
      if (nord.ne.nord_old) aa(1:nrdof,1:nrdof)=0.d0
      bb(1:nrdof,1:3)=0.d0
c
c  ...loop through integration points
      do l=1,nint
c
c  .....Gauss point and weight
        t=tloc(l) ; wa=waloc(l)
c
c  .....master coordinates of the integration point
        call edge_param(Type,Iedg,t, xi,dxidt)
c
c  .....compute element shape functions
        call initiate_order(Type, nordloc)
        call shape3H(Type,xi,nordloc,Nedge_orient,Nface_orient,
     .               nrdofH,vshapH,dvshapH)
c
c  .....evaluate reference coordinates of the point
        eta(1:3)=0.d0 ; detadxi(1:3,1:3)=0.d0
        do k=1,nrdofH
          eta(1:3) = eta(1:3) + Etav(1:3,k)*vshapH(k)
          do i=1,3
            detadxi(1:3,i) = detadxi(1:3,i) + Etav(1:3,k)*dvshapH(i,k)
          enddo
        enddo
c
c  .....evaluate -dxdxi resulting from vertex dof
        dxdxi(1:3,1:3)=0.d0
        do k=1,nrdofH
          do i=1,3
            dxdxi(1:3,i) = dxdxi(1:3,i) - Xnod(1:3,k)*dvshapH(i,k)
          enddo
        enddo
c
c  .....evaluate derivatives of the physical coordinates wrt
c       the reference coordinate
        select case(Iflag)
        case(5)
          call prism(No,eta, x,dxdeta)
        case(6)
          call  hexa(No,eta, x,dxdeta)
        case(7)
          call tetra(No,eta, x,dxdeta)
        case(8)
          call pyram(No,eta, x,dxdeta)
        case default
          write(*,7002) Type
 7002     format('hpedge: Type = ',a5)
          stop 1
        end select
        if (iprint.eq.1) then
          write(*,7003)l,eta(1:3)
 7003     format('hpedge: l,eta = ',i3,2x,3e12.5)
          do ivar=1,3
            write(*,7004) ivar,x(ivar),dxdeta(ivar,1:3)
 7004       format('hpedge: ivar,x,dxdeta = ',i2,2x,e12.5,2x,3e12.5)
          enddo
        endif
c
c  .....add exact dxdxi
        do i=1,3
          do j=1,3
            dxdxi(1:3,i) = dxdxi(1:3,i) + dxdeta(1:3,j)*detadxi(j,i)
          enddo
        enddo
c
c  .....evaluate derivative wrt edge coordinate
        dxdt(1:3) = 0.d0
        do i=1,3
          dxdt(1:3) = dxdt(1:3) + dxdxi(1:3,i)*dxidt(i)
        enddo
        if (iprint.eq.1) then
          write(*,7005) x(1:3),dxdt(1:3)
 7005     format(' hpedge: x,dxdt = ',2(3e12.5,2x))
        endif
c
c  .....evaluate derivatives of the 1D shape functions at the point
        call shapHbe(t,nord,Nedge_orient(Iedg), nrdof,vshap,dvshap)
c
c  .....first loop through dof
        do j=1,nrdof
c
c  .......accumulate for the load vector
          bb(j,1:3) = bb(j,1:3) + dxdt(1:3)*dvshap(j)*wa
c
c  .......accumulate for the stiffness matrix
          if (nord.ne.nord_old) then
            do i=1,nrdof
              aa(i,j) = aa(i,j) + dvshap(i)*dvshap(j)*wa
            enddo
          endif
        enddo
c
c  ...end of the loop through integration points
      enddo
c
c
      if (iprint.eq.1) then
        write(*,*) 'hpedge: LOAD VECTOR AND STIFFNESS MATRIX = '
        do ii=1,nrdof
          write(*,7006) bb(ii,1:3),aa(ii,1:nrdof)
 7006     format(3(e12.5,2x),';',2x,10(e12.5,2x))
        enddo
      endif
c
c  ...solve the linear system...
c     pg, Apr 12: replaced with calls to Lapack
      na = MAXP-1
      if (nord.ne.nord_old) then
ccc        call tri(aa,na,nrdof)
c
        call dgetrf(nrdof,nrdof,aa,na,ipiv,info)
        if (info.ne.0) then
          write(*,*)'hpedge: H1 DGETRF RETURNED INFO = ',info
          call logic_error(FAILURE,__FILE__,__LINE__)
        endif
c
        nord_old = nord
      endif
c
      do ivar=1,3
ccc        call rhsub(aa,uu,bb(1,ivar),na,nrdof)
c
        uu(1:nrdof)=bb(1:nrdof,ivar)
        call dlaswp(1,uu,na,1,nrdof,ipiv,1)
        call dtrsm('L','L','N','U',nrdof,1,1.d0,aa,na,uu,na)
        call dtrsm('L','U','N','N',nrdof,1,1.d0,aa,na,uu,na)
c
        do i=1,nrdof
          Xdof(ivar,i) = uu(i)
        enddo
      enddo
      if (iprint.eq.1) then
        do ii=1,nrdof
          write(*,7007) ii,(Xdof(jj,ii),jj=1,3)
 7007     format(' i = ',i2,' ; Xnod = ',3(e12.5,2x))
        enddo
        call pause
      endif
c
c
      endsubroutine hpedge
