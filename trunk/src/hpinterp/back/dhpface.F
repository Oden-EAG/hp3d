c----------------------------------------------------------------------
c
c   routine name       - dhpface
c
c----------------------------------------------------------------------
c
c   latest revision    - Dec 10
c
c   purpose            - routine determines solution dof for a mid-face
c                        Dirichlet node of an element
c
c   arguments :
c     in:
c         Iflag        - a flag specifying which of the objects the
c                        mid-face node is generated on
c                        = 5  for a prism
c                        = 6  for a hexahedron
c                        = 7  for a tetrahedron
c                        = 8  for a pyramid
c         No           - number of a specific object (curve,triangle,
c                        rectangle...)
c         Etav         - reference coordinates of the element vertices
c         Type         - element (middle node) type
c         Icase        - face-node case
c         Nedge_orient - edge orientation
c         Nface_orient - face orientation (not needed really)
c         Norder       - element order
c         If           - face number
c         ZdofH        - H1 dof for the element (vertex and edge values
c                        really)
c         ZdofE        - H(curl) dof for the element (edge values really)
c     out:
c         ZnodH        - H1 dof for the edge
c         ZnodE        - H(curl) dof for the edge
c         ZnodV        - H(div) dof for the edge
c
c     H(curl) and H(div) cases MISSING !!!!
c
c----------------------------------------------------------------------
c
      subroutine dhpface(Iflag,No,Etav,
     .                   Type,Icase,
     .                   Nedge_orient,Nface_orient,Norder,If,
     .                   ZdofH,ZnodH, ZdofE,ZnodE, ZnodV)
c
      use element_data
      use data_structure3D
c
#include "syscom.blk"
#include "cint.blk"
c
      dimension Etav(3,8)
      character(len=4) :: Type
      dimension Nedge_orient(12),Nface_orient(6),Norder(19)
      dimension ZdofH(MAXEQNH,MAXbrickH),
     .          ZnodH(NRCOMS*NREQNH(Icase),*),
     .          ZdofE(MAXEQNE,MAXbrickE),
     .          ZnodE(NRCOMS*NREQNE(Icase),*),
     .          ZnodV(NRCOMS*NREQNV(Icase),*)
c
c  ...local element order
      dimension nordloc(19)
c
c  ...element shape functions
      parameter(maxvshap=8+12*(MAXP-1)+6*MAXmdlqH)
      dimension vshapH(maxvshap),dvshapH(3,maxvshap)
c
c  ...face shape functions
      dimension vshap(MAXmdlqH),dvshap(2,MAXmdlqH)
c
c  ...geometry
      dimension xi(3),eta(3),x(3),
     .          dxdt(3,2),dxidt(3,2),
     .          dxdxi(3,3),dxdeta(3,3),detadxi(3,3)
c
c  ...Dirichlet BC data at a point
      dimension zvalH(MAXEQNH),zdvalH(MAXEQNH,3),
     .          zvalE(3,MAXEQNE),
     .          zvalV(3,MAXEQNV),
     .          zdvalHdxi(MAXEQNH,3),zdvalHdt(MAXEQNH,2)
c
c  ...physics flags decoded
      dimension ncase(NR_PHYSA)
c
c  ...quadrature data
      dimension tloc(2,MAXquadH),waloc(MAXquadH)
c
c  ...linear system matrices
      dimension aa(MAXmdlqH,MAXmdlqH),zb(MAXmdlqH,MAXEQNH),
     .          zu(MAXmdlqH,MAXEQNH)
c
      integer, parameter :: na = MAXmdlqH
c
      save nord_old,aa
      data nord_old /0/
c
c----------------------------------------------------------------------
c
      iprint=0
c
c  ...printing
      if (iprint.eq.1) then
        write(*,7001) Iflag,No
 7001   format('dhpface: Iflag,No = ',i1,2x,i8,' Etav = ')
        do iv=1,nvert(Type)
          write(*,7002) iv,Etav(1:3,iv)
 7002     format('         iv = ',i1,2x,3f8.3)
        enddo
        write(*,7003) Type,Norder,If
 7003   format('         Type,Norder,If = ',a5,2x,19i4,5x,i1)
ccc        call pause
      endif
c
c  ...get face order of approximation
      in = nedge(Type)+If
      nord = Norder(in)
      select case(face_type(Type,If))
      case('mdlt','tria')
        if (nord.lt.2) return
ccc        if (nord.lt.3) return
      case('mdlq','rect')
        call decode(Nord, nordh,nordv)
        if ((nordh.lt.2).or.(nordv.lt.2)) return
      end select
c
c  ...set quadrature data
      call face_order(Type,If,Norder, nordloc)
      call set_2Dint(face_type(Type,If),nordloc, nint,tloc,waloc)
c
c  ...initiate matrices
      if (nord.ne.nord_old) aa = 0.d0
      zb = ZERO
      call initiate_order(Type, nordloc)
      ne = nedge(Type)
      nordloc(1:ne) = Norder(1:ne)
c
c  ...loop through integration points
      do l=1,nint
c
c  .....determine element master coordinates of the integration point
        call face_param(Type,If,tloc(1:2,l), xi,dxidt)
c
c  .....compute element shape functions
        call shape3H(Type,xi,nordloc,Nedge_orient,Nface_orient,
     .               nrdofH,vshapH,dvshapH)
c
c  .....evaluate reference coordinates of the point
        eta(1:3) = 0.d0; detadxi(1:3,1:3)=0.d0
        do k=1,nvert(Type)
          eta(1:3) = eta(1:3) + etav(1:3,k)*vshapH(k)
          do i=1,3
            detadxi(1:3,i) = detadxi(1:3,i) + etav(1:3,k)*dvshapH(i,k)
          enddo
        enddo
c
c  .....evaluate -zdvalHdxi resulting from vertex and edge dof
        zdvalHdxi(1:MAXEQNH,1:3) = ZERO
        do k=1,nrdofH
          do i=1,3
            zdvalHdxi(1:MAXEQNH,i) = zdvalHdxi(1:MAXEQNH,i)
     .                             - ZdofH(1:MAXEQNH,k)*dvshapH(i,k)
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
          write(*,7010) Type
 7010     format('dhpface: Type = ',a5)
          stop 1
        end select
c
c  .....evaluate derivatives wrt master coordinates
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
c  .....evalute derivatives wrt face coordinates
ccc        with interpolation done over the master element, dxdt is not needed
ccc        dxdt(1:3,1:2) = 0.d0
        zdvalHdt(1:MAXEQNH,1:2) = ZERO
        do i=1,2
          do j=1,3
            zdvalHdt(1:MAXEQNH,i) = zdvalHdt(1:MAXEQNH,i)
     .                            + zdvalHdxi(1:MAXEQNH,j)*dxidt(j,i)
ccc            dxdt(1:3,i) = dxdt(1:3,i) + dxdxi(1:3,j)*dxidt(j,i)
          enddo
        enddo
c
c  .....evaluate the 2D middle node shape functions at the point
        call shapHb2(face_type(Type,If),tloc(1:2,l),nord,
     .               Nface_orient(If), nrdof,vshap,dvshap)
c
c  .....first loop through dof
        do j=1,nrdof
c
c  .......accumulate for the load vector
          zb(j,1:MAXEQNH) = zb(j,1:MAXEQNH)
     .                    + (zdvalHdt(1:MAXEQNH,1)*dvshap(1,j)
     .                     + zdvalHdt(1:MAXEQNH,2)*dvshap(2,j))*waloc(l)
c
c  .......accumulate for the stiffness matrix
          if (nord.ne.nord_old) then
            do i=1,nrdof
              aa(i,j) = aa(i,j) + (dvshap(1,i)*dvshap(1,j)
     .                           + dvshap(2,i)*dvshap(2,j))*waloc(l)
            enddo
          endif
        enddo
c
c  ...end of loop through integration points
      enddo
c
c  ...printing
      if (iprint.eq.1) then
        write(*,*) 'dhpface: LOAD VECTOR AND STIFFNESS MATRIX = '
        do ii=1,nrdof
          write(*,7015) zb(ii,1:MAXEQNH),aa(ii,1:nrdof)
ccc 7015     format(1x,4(2e12.5,2x),3x,10f8.3)
 7015     format(1x,e12.5,4x,10(e12.5,2x))
        enddo
      endif
c
c  ...solve the linear system...
ccc      na=MAXmdlqH
      if (nord.ne.nord_old) then
        call tri(aa,na,nrdof)
ccc        write(*,*)'aa = '
ccc        do ii=1,nrdof
ccc          write(*,8000)aa(ii,1:nrdof)
ccc8000      format(1x,10(e12.5,2x))
ccc        enddo
        nord_old = nord
      endif
      do ivar=1,MAXEQNH
       call zrhsub(aa,zu(1:nrdof,ivar),zb(1:nrdof,ivar),na,nrdof)
      enddo
c
c  ...save the dof skipping irrelevant entries
      call decod(Icase,2,NR_PHYSA, ncase)
      if (iprint.eq.1) then
        write(*,7020) ncase
 7020   format('dhpface: ncase = ',5i1)
      endif
c
      ivarH=0; nvarH=0
c
c  ...a shortcut !!!
      ivarE=0; nvarE=0; nrdofE=nord*(nord+2)
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

ccc                write(*,*)'nvarH = ',nvarH
ccc                write(*,*)'nrdof = ',nrdof

                ZnodH(nvarH,1:nrdof) = zu(1:nrdof,ivarH)

ccc                write(*,*)'ivarH = ',ivarH
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
      end subroutine
