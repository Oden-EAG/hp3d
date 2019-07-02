c----------------------------------------------------------------------
c
c   routine name       - dhpedge
c
c----------------------------------------------------------------------
c
c   latest revision    - Dec 10
c
c   purpose            - routine determines solution dof for a mid-edge
c                        Dirichlet node of an element
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
c         Icase        - the mid-edge node case
c         Nedge_orient - edge orientation (not needed really)
c         Nface_orient - face orientation (not needed really)
c         Norder       - element order
c         Iedg         - edge number
c         ZdofH        - H1 dof for the element (vertex values really)
c     out:
c         ZnodH        - H1 dof for the edge
c         ZnodE        - H(curl) dof for the edge
c     H(curl) case MISSING !!!!
c
c----------------------------------------------------------------------
c
      subroutine dhpedge(Iflag,No,Etav,
     .                   Type,Icase,
     .                   Nedge_orient,Nface_orient,Norder,Iedg,
     .                   ZdofH,ZnodH,
     .                   ZnodE)
c
      use element_data
      use data_structure3D
#include "syscom.blk"
#include "cint.blk"
c
      dimension Etav(3,8)
      character(len=4) :: Type
      dimension Nedge_orient(12),Nface_orient(6),Norder(19)
      dimension ZdofH(MAXEQNH,MAXbrickH),
     .          ZnodH(NRCOMS*NREQNH(Icase),*),
     .          ZnodE(NRCOMS*NREQNE(Icase),*)

c
c  ...local element order
      dimension nordloc(19)
c
c  ...element shape functions (vertex shape functions only)
      dimension vshapH(8),dvshapH(3,8)
c
c  ...edge shape functions
      dimension vshap(MAXP-1),dvshap(MAXP-1)
c
c  ...geometry
      dimension xi(3),eta(3),x(3),
     .          dxdt(3),dxidt(3),
     .          dxdxi(3,3),dxdeta(3,3),detadxi(3,3)
c
c  ...Dirichlet BC data at a point
      dimension zvalH(MAXEQNH),zdvalH(MAXEQNH,3),
     .          zvalE(3,MAXEQNE),
     .          zvalV(3,MAXEQNV),
     .          zdvalHdxi(MAXEQNH,3),zdvalHdt(MAXEQNH)
c
c  ...physics flags decoded
      dimension ncase(NR_PHYSA)
c
c  ...linear system matrices
      dimension aa(MAXP-1,MAXP-1),zb(MAXP-1,MAXEQNH),zu(MAXP-1,MAXEQNH)
c
      save nord_old,aa
      data nord_old /1/
c
c----------------------------------------------------------------------
c
      iprint=0
c
      nord = Norder(Iedg)
      if (nord.lt.2) return
c
c  ...number of integration points
      nint = nord+1
c
c  ...initiate matrices
      if (nord.ne.nord_old) aa = 0.d0
      zb = ZERO
c
c  ...loop through integration points
      do l=1,nint
c
c  .....determine master coordinates of the integration point
        t = XIGAUS1(l,nint)
        call edge_param(Type,Iedg,t, xi,dxidt)
        wa = WAGAUS1(l,nint)
        if (iprint.eq.1) then
          write(*,3067) l,t
 3067     format('dphmedge: INTEGRATION POINT l,t = ',i3,e12.5,
     .           '  ************************************')
          write(*,7028) xi(1:3),dxidt(1:3)
 7028     format('dhpmedge: xi,dxidt = ',3f8.3,2x,3f8.3)
        endif
c
c  .....compute element shape functions
        call initiate_order(Type, nordloc)
        call shape3H(Type,xi,nordloc,Nedge_orient,Nface_orient,
     .               nrdofH,vshapH,dvshapH)
c
c  .....evaluate reference coordinates of the point
        eta(1:3) = 0.d0; detadxi(1:3,1:3)=0.d0
        do k=1,nrdofH
          eta(1:3) = eta(1:3) + etav(1:3,k)*vshapH(k)
          do i=1,3
            detadxi(1:3,i) = detadxi(1:3,i) + etav(1:3,k)*dvshapH(i,k)
          enddo
        enddo
c
c  .....evaluate -zdvalHdxi resulting from vertex dof
        zdvalHdxi(1:MAXEQNH,1:3) = ZERO
        do k=1,nrdofH
          do i=1,3
            zdvalHdxi(1:MAXEQNH,i) = zdvalHdxi(1:MAXEQNH,i)
     .                             - ZdofH(1:MAXEQNH,k)*dvshapH(i,k)
          enddo
          if (iprint.eq.1) then
            write(*,7039) k,ZdofH(1:MAXEQNH,k)
 7039       format('dhpedge: k,ZdofH(1:MAXEQNH,k) = ',
     .             /,i2,2x,4(2e12.5,2x))
          endif
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
          write(*,7001) Type
 7001     format('dhpedge: Type = ',a5)
          stop 1
        end select
c
c  .....evaluate the Dirichlet data
        call dirichlet(x,Icase, zvalH,zdvalH,zvalE,zvalV)
        if (iprint.eq.1) then
          write(*,7893) x(1:3)
 7893     format('dhpedge: x = ',3f8.3)
          do ivar=1,MAXEQNH
            write(*,7006) ivar,zvalH(ivar),zdvalH(ivar,1:3)
 7006       format('dhpedge: ivar,zvalH(ivar),zdvalH(ivar,1:3) = ',
     .            /,i2,2x,2e12.5,3x,3(2e12.5,2x))
          enddo
        endif
c
c  .....evaluate dxdxi (using the exact geometry)
        dxdxi(1:3,1:3)=0.d0
        do i=1,3
          do j=1,3
            dxdxi(1:3,i) = dxdxi(1:3,i) + dxdeta(1:3,j)*detadxi(j,i)
          enddo
        enddo
c
c  .....add derivatives of Dirichlet data wrt master coordinates
        do i=1,3
          do j=1,3
            zdvalHdxi(1:MAXEQNH,i) = zdvalHdxi(1:MAXEQNH,i)
     .                             + zdvalH(1:MAXEQNH,j)*dxdxi(j,i)
          enddo
        enddo
c
c  .....evalute derivative wrt edge coordinate
ccc        with interpolation done over the master element, dxdt is not needed
ccc        dxdt(1:3) = 0.d0
        zdvalHdt(1:MAXEQNH) = ZERO
        do i=1,3
          zdvalHdt(1:MAXEQNH) = zdvalHdt(1:MAXEQNH)
     .                        + zdvalHdxi(1:MAXEQNH,i)*dxidt(i)
ccc          dxdt(1:3) = dxdt(1:3) + dxdxi(1:3,i)*dxidt(i)
        enddo
        if (iprint.eq.1) then
          write(*,7804) zdvalHdt(1:MAXEQNH)
 7804     format('dhpmedge: zdvalHdt = ',/,4(2e12.5,2x))
        endif

c
c  .....evaluate derivatives of the 1D shape functions at the point
        call shapHbe(t,nord,Nedge_orient(Iedg), nrdof,vshap,dvshap)
c
c  .....first loop through dof
        do j=1,nrdof
c
c  .......accumulate for the load vector
          zb(j,1:MAXEQNH) = zb(j,1:MAXEQNH)
     .                    + zdvalHdt(1:MAXEQNH)*dvshap(j)*wa
ccc          write(*,*) 'j,vshap(j),dvshap(j) = ',j,vshap(j),dvshap(j)
c
c  .......accumulate for the stiffness matrix
          if (nord.ne.nord_old) then
            do i=1,nrdof
              aa(i,j) = aa(i,j) + dvshap(i)*dvshap(j)*wa
            enddo
          endif
        enddo
        if (iprint.eq.1) then
          do j=1,nrdof
            write(*,7378) j,zb(j,1:MAXEQNH)
 7378       format('dhpmedge: j,zb(j,1:MAXEQNH) = ',
     .             /,i3,2x,4(2e12.5,2x))
          enddo
          call pause
        endif
c
c  ...end of the loop through integration points
      enddo
c
c
      if (iprint.eq.1) then
        write(*,*) 'dhpedge: LOAD VECTOR AND STIFFNESS MATRIX = '
        do ii=1,nrdof
          write(*,7036) zb(ii,1:MAXEQNH)
 7036     format(4(2e12.5,2x))
          write(*,7005) aa(ii,1:nrdof)
 7005     format(1x,10f8.3)
        enddo
      endif
c
c  ...solve the linear system...
      na = MAXP-1
      if (nord.ne.nord_old) then
        call tri(aa,na,nrdof)
        nord_old = nord
      endif
      do ivar=1,MAXEQNH
        call zrhsub(aa,zu(1:nrdof,ivar),zb(1:nrdof,ivar),na,nrdof)
      enddo
      if (iprint.eq.1) then
        write(*,*) 'dhpedge: k,zu(k) = '
        do k=1,nrdof
          write(*,*) k,zu(k,1:MAXEQNH)
        enddo
        call pause
      endif
c
c  ...save the dof, skipping irrelevant entries
      call decod(Icase,2,NR_PHYSA, ncase)
      if (iprint.eq.1) then
        write(*,7002) ncase
 7002   format('dhpedge: ncase = ',5i1)
      endif
c
      ivarH=0; nvarH=0
c
c  ...a shortcut !!!
      ivarE=0; nvarE=0; nrdofE = nord
c
c  ...loop through multiple copies of variables
      do j=1,NRCOMS
c
c  .....loop through physical attributes
        do i=1,NR_PHYSA
c
c  .......loop through components
          do k=1,NR_COMP(i)
            select case(DTYPE(i))
            case('contin')
              ivarH=ivarH+1
              if (ncase(i).eq.1) then
                nvarH=nvarH+1
                ZnodH(nvarH,1:nrdof) = zu(1:nrdof,ivarH)
              endif
            case('tangen')
              ivarE=ivarE+1
              if (ncase(i).eq.1) then
                nvarE=nvarE+1
                ZnodE(nvarE,1:nrdofE) = ZERO
              endif
            end select
          enddo
        enddo
      enddo
c
c
      end

