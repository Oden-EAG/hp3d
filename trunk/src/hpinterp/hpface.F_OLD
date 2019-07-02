c----------------------------------------------------------------------
c   latest revision    - Mar 12
c
c   purpose            - routine determines geometry dof for a mid-face
c                        node of an element by performing H^1_0
c                        projection on the REFERENCE (NOT master)
c                        element
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
c         Nedge_orient - edge orientation
c         Nface_orient - face orientation (not needed really)
c         Norder       - element order
c         If           - face number
c         Xnod         - element gdofs so far (vertices and edges)
c     out:
c         Xdof         - geometry dofs for the mid-face node
c----------------------------------------------------------------------
c
      subroutine hpface(Iflag,No,Etav,Type,
     .                  Nedge_orient,Nface_orient,Norder,If,
     .                  Xnod,Xdof)
c
      use element_data
      use parameters
      use cross_product_module
#include "syscom.blk"
      common /cshapHb2/ iprint_shapHb2
      common /chpface/ iprint_hpface
c
      dimension Etav(3,8),Xnod(3,MAXbrickH),Xdof(3,*)
      character(len=4) :: Type
      dimension Nedge_orient(12),Nface_orient(6),Norder(19)
c
c  ...local element order
      dimension nordloc(19)
c
c  ...3D element shape functions
      dimension shapH(maxbrickH),dshapH(3,maxbrickH)
c
c  ...mid-face shape functions
      dimension shap(MAXmdlqH),dshap(2,MAXmdlqH)
c
c  ...geometry
      dimension xi(3),eta(3),x(3),
     .          dxdt(3,2),dxidt(3,2),
     .          dxdxi(3,3),dxdeta(3,3),detadxi(3,3),
     .          detadt(3,2),rnjac(3)
c
c  ...cobasis vectors, metric for surface jacobians
      dimension dual(3,2),gg(2,2)
c
c  ...quadrature data
      dimension tloc(2,MAXquadH),waloc(MAXquadH)
c
c  ...linear system matrices
      dimension aa(MAXmdlqH,MAXmdlqH),bb(MAXmdlqH,3),uu(MAXmdlqH)
      dimension ipiv(MAXmdlqH)
c
c----------------------------------------------------------------------
c
      iprint=iprint_hpface
      iprint_shapHb2=0
      if (iprint.eq.1) then
        write(*,7001) Iflag,No
 7001   format('hpface: Iflag,No = ',i1,2x,i8,' Etav = ')
        do iv=1,nvert(Type)
          write(*,7002) iv,Etav(1:3,iv)
 7002     format('        iv = ',i1,2x,3f8.3)
        enddo
        write(*,7003) Type,Norder,If
 7003   format('        Type,Norder,If = ',a4,2x,19i4,5x,i1)
        call pause
      endif
c
      in = nedge(Type)+If
      nord = Norder(in)
      select case(face_type(Type,If))
      case('mdlt','tria')
        if (nord.lt.3) return
      case('mdlq','rect')
        call decode(nord, nordh,nordv)
        if ((nordh.lt.2).or.(nordv.lt.2)) return
      end select
c
c  ...set quadrature data
      call face_order(Type,If,Norder, nordloc)
      if (iprint.eq.1) then
        write(*,6999)nordloc(1:5)
 6999   format('        nordloc = ',5i4)
      endif
      call set_2Dint(face_type(Type,If),nordloc, nint,tloc,waloc)
c
c  ...initiate matrices
      aa = 0.d0
      bb = 0.d0
c
c  ...set the order for the element
      call initiate_order(Type, nordloc)
      ne = nedge(Type)
      nordloc(1:ne) = Norder(1:ne)
c
c  ...loop through integration points
      do l=1,nint
c
c  .....determine element master coordinates of the integration point
        call face_param(Type,If,tloc(1:2,l), xi,dxidt)
        if (iprint.eq.1) then
          write(*,7045) l,xi(1:3)
 7045     format('hpface: l, xi = ',i4,2x,3f8.3)
          do ivar=1,3
            write(*,7046) ivar,dxidt(ivar,1:2)
 7046       format('hpface: dxidt(',i1,',:) = ',2f8.3)
          enddo
          call pause
        endif
c
c
c  .....compute element shape functions
        call shape3H(Type,xi,nordloc,Nedge_orient,Nface_orient,
     .               nrdofH,shapH,dshapH)
c
c  .....evaluate reference coordinates of the point
        eta(1:3) = 0.d0; detadxi(1:3,1:3)=0.d0
        do k=1,nvert(Type)
          eta(1:3) = eta(1:3) + Etav(1:3,k)*shapH(k)
          do i=1,3
            detadxi(1:3,i) = detadxi(1:3,i) + Etav(1:3,k)*dshapH(i,k)
          enddo
        enddo
c
c  .....evaluate surface jacobian
        do i=1,3
          do j=1,2
            s = 0.d0
            do k=1,3
              s = s + detadxi(i,k)*dxidt(k,j)
            enddo
            detadt(i,j) = s
          enddo
        enddo
        call cross_product3D(dxidt(1:3,1),dxidt(1:3,2), rnjac)
        call norm(rnjac, Sjac)
        rnjac(1:3) = rnjac(1:3)/Sjac
c
c  .....evaluate dual basis vectors
        call cross_product3D(dxidt(1:3,2),rnjac, dual(1:3,1))
        dual(1:3,1) = dual(1:3,1)/Sjac
        call cross_product3D(rnjac,dxidt(1:3,1), dual(1:3,2))
        dual(1:3,2) = dual(1:3,2)/Sjac
c
c  .....evaluate the dot product of the dual basis vectors
        do i=1,2
          do j=1,2
            call scalar_product(dual(1:3,i),dual(1:3,j), gg(i,j))
          enddo
        enddo
        weight = waloc(l)*Sjac
c
c  .....evaluate -dxdxi obtained with vertex and edge dof
        dxdxi(1:3,1:3)=0.d0
        do k=1,nrdofH
          if (iprint.eq.1) then
            write(*,7053) k,Xnod(1:3,k)
 7053       format('k = ',i2,2x,' Xnod = ',3e12.5)
          endif
          do i=1,3
            dxdxi(1:3,i) = dxdxi(1:3,i) - Xnod(1:3,k)*dshapH(i,k)
          enddo
        enddo
        if (iprint.eq.1) then
          write(*,*) 'hpface: dxdxi = '
          do i=1,3
            write(*,7054) i,dxdxi(i,1:3)
 7054       format('ivar = ',i2,3x,3e12.5)
          enddo
          call pause
        endif
c
c  .....evaluate derivatives of the physical coordinates wrt
c       the reference coordinates
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
 7010     format('hpface: Type = ',a5)
          stop 1
        end select
        if (iprint.eq.1) then
          write(*,*) 'eta = ',eta
          do ivar=1,3
            write(*,7011) ivar,x(ivar),dxdeta(ivar,1:3)
 7011       format('hpface: ivar,x(ivar),dxdeta(ivar,1:3) = ',
     .                      i2,2x,f10.5,3x,3f10.5)
          enddo
        endif
c
c  .....add to the -vertex/edge derivatives the exact (GMP) dxdxi
        do i=1,3
          do j=1,3
            dxdxi(1:3,i) = dxdxi(1:3,i) + dxdeta(1:3,j)*detadxi(j,i)
          enddo
        enddo
c
c  .....evalute derivatives wrt face coordinates
        dxdt(1:3,1:2) = 0.d0
        do i=1,2
          do j=1,3
            dxdt(1:3,i) = dxdt(1:3,i) + dxdxi(1:3,j)*dxidt(j,i)
          enddo
        enddo
c
c  .....evaluate the 2D middle node shape functions at the point
        call shapHb2(face_type(Type,If),tloc(1:2,l),nord,
     .               Nface_orient(If), nrdof,shap,dshap)
c
c  .....first loop through dof
        do j=1,nrdof
c
c  .......accumulate for the load vector
          bb(j,1:3) = bb(j,1:3)
     .              + (dxdt(1:3,1)*dshap(1,j)*gg(1,1)
     .               + dxdt(1:3,1)*dshap(2,j)*gg(1,2)
     .               + dxdt(1:3,2)*dshap(1,j)*gg(2,1)
     .               + dxdt(1:3,2)*dshap(2,j)*gg(2,2))*weight
c
c  .......accumulate for the stiffness matrix
          do i=1,nrdof
            aa(i,j) = aa(i,j)
     .              + (dshap(1,i)*dshap(1,j)*gg(1,1)
     .               + dshap(1,i)*dshap(2,j)*gg(1,2)
     .               + dshap(2,i)*dshap(1,j)*gg(2,1)
     .               + dshap(2,i)*dshap(2,j)*gg(2,2))*weight
          enddo
        enddo
c
c  ...end of the loop through integration points
      enddo
c
c
      if (iprint.eq.1) then
        write(*,*) 'hpface: LOAD VECTOR  = '
        do ii=1,nrdof
          write(*,7015) bb(ii,1:3)
 7015     format(1x,10e12.5)
        enddo
        call pause
        write(*,*) 'hpface: STIFFNESS MATRIX = '
        do ii=1,nrdof
          write(*,7016) ii,aa(ii,1:nrdof)
 7016     format(i2,2x,10(10e12.5,/,4x))
          call pause
        enddo
        call pause
      endif
c
c  ...solve the linear system...
c     pg, Apr 12: replaced with calls to Lapack
      na = MAXmdlqH
ccc      call tri(aa,na,nrdof)
c
      call dgetrf(nrdof,nrdof,aa,na,ipiv,info)
      if (info.ne.0) then
        write(*,*)'dhpedge: H1 DGETRF RETURNED INFO = ',info
        call logic_error(FAILURE,__FILE__,__LINE__)
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
      if (iprint.ge.1) then
        do k=1,nrdof
          write(*,7020) k,Xdof(1:3,k)
 7020     format('hpface: k,Xdof(1:3,k) = ',i3,2x,3e12.5)
        enddo
        call pause
      endif
c
c
      end subroutine hpface
