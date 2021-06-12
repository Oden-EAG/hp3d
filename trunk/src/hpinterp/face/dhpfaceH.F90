!
#include "typedefs.h"
!
!-----------------------------------------------------------------------
!> Purpose : determine H1 face dof interpolating H1 Dirichlet data using
!            PB interpolation
!  NOTE:     the interpolation (projection) is done in the reference space
!!
!! @param[in]  Iflag        - a flag specifying which of the objects the
!!                            face is on: 5 pris, 6 hexa, 7 tetr, 8 pyra
!! @param[in]  No           - number of a specific object
!! @param[in]  Etav         - reference coordinates of the element vertices
!! @param[in]  Type         - element (middle node) type
!! @param[in]  Icase        - the face node case
!! @param[in]  Bcond        - the face node BC flag
!! @param[in]  Nedge_orient - edge orientation
!! @param[in]  Nface_orient - face orientation
!! @param[in]  Norder       - element order
!! @param[in]  Iface        - face number
!! @param[in]  ZdofH        - H1 dof for the element (vertex and edge values)
!!
!! @param[out] ZnodH        - H1 dof for the face
!-----------------------------------------------------------------------
  subroutine dhpfaceH(Mdle,Iflag,No,Etav,Type,Icase,Bcond, &
                      Nedge_orient,Nface_orient,Norder,Iface, &
                      ZdofH, ZnodH)
  use control
  use parameters
  use physics
  use element_data
  implicit none
!
! ** Arguments
!-----------------------------------------------------------------------
  integer,                                    intent(in)  :: Iflag,No,Mdle
  integer,                                    intent(in)  :: Icase,Bcond,Iface
  real(8), dimension(3,8),                    intent(in)  :: Etav
  character(len=4),                           intent(in)  :: Type
  integer, dimension(12),                     intent(in)  :: Nedge_orient
  integer, dimension(6),                      intent(in)  :: Nface_orient
  integer, dimension(19),                     intent(in)  :: Norder
!
  VTYPE,   dimension(MAXEQNH,MAXbrickH),      intent(in)    :: ZdofH
  VTYPE,   dimension(NRCOMS*NREQNH(Icase),*), intent(inout) :: ZnodH
!
! ** Locals
!-----------------------------------------------------------------------
!
! face nodes order (to set the quadrature)
  integer, dimension(5)                 :: norder_face
!
! quadrature
  integer                               :: l,nint
  real(8), dimension(2, MAXquadH)       :: xi_list
  real(8), dimension(   MAXquadH)       :: wa_list
  real(8)                               :: wa, weight
!
! work space for shape3DH
  integer                               :: nrdofH
  integer, dimension(19)                :: norder_1
  real(8), dimension(MAXbrickH)         :: shapH
  real(8), dimension(3,MAXbrickH)       :: gradH
!
! derivatives of a shape function wrt reference coordinates
  real(8), dimension(3)                 :: duHdeta,dvHdeta
!
! dot product
  real(8)                               :: prod
!
! geometry
  real(8)                               :: rjac,bjac
  real(8), dimension(2)                 :: t
  real(8), dimension(3)                 :: xi,eta,rn,x
  real(8), dimension(3,2)               :: dxidt,detadt
  real(8), dimension(3,3)               :: detadxi,dxideta,dxdeta
!
! Dirichlet BC data at a point
  VTYPE :: zvalH(  MAXEQNH), zdvalH(  MAXEQNH,3), &
           zvalE(3,MAXEQNE), zdvalE(3,MAXEQNE,3), &
           zvalV(3,MAXEQNV), zdvalV(3,MAXEQNV,3)
!
! derivatives of Dirichlet date wrt reference coordinates
  VTYPE :: zdvalHdeta(MAXEQNH,3)
!
! work space for linear solvers
  integer                               :: naH,info
  real(8), dimension(MAXMdlqH,MAXMdlqH) :: aaH
  integer, dimension(MAXMdlqH)          :: ipivH
!
! load vector and solution
  VTYPE,   dimension(MAXMdlqH,MAXEQNH)  :: zbH,zuH
#if C_MODE
  real(8), dimension(MAXMdlqH,MAXEQNH)  :: duH_real,duH_imag
#endif
!
! decoded case and BC flags for the edge node
  integer, dimension(NR_PHYSA)          :: ncase
  integer, dimension(NRINDEX)           :: ibcnd
!
! misc work space
  integer :: nrv,nre,nrf,i,j,k,ie,ivarH,nvarH,kj,ki,&
             ndofH_face,ndofE_face,ndofV_face,ndofQ_Face,nsign,ic
!
  logical :: is_homD
!
#if DEBUG_MODE
  integer :: iprint = 0
#endif
!
!-----------------------------------------------------------------------
!
  nrv = nvert(Type); nre = nedge(Type); nrf = nface(Type)
!
#if DEBUG_MODE
  if (iprint.eq.1) then
     write(*,7010) Mdle,Iflag,No,Icase,Iface,Type
7010 format('dhpfaceH: Mdle,Iflag,No,Icase,Iface,Type = ',5i4,a4)
     write(*,7020) Etav(1:3,1:nrv)
7020 format('          Etav = ',8(3f6.2,1x))
     write(*,7030) Nedge_orient(1:nre)
7030 format('          Nedge_orient = ',12i2)
     write(*,7040) Nface_orient(1:nrf)
7040 format('          Nface_orient = ',6i2)
     write(*,7050) Norder(1:nre+nrf+1)
7050 format('          Norder = ',19i4)
     call pause
  endif
#endif
!
! check if a homogeneous Dirichlet node
  call homogenD('contin',Icase,Bcond, is_homD,ncase,ibcnd)
  if (is_homD) then
    zuH = ZERO
    go to 100
  endif
!
! determine # of dof for the face node
  call ndof_nod(face_type(Type,Iface),Norder(nre+Iface), &
                ndofH_face,ndofE_face,ndofV_face,ndofQ_face)
!
! if # of dof is zero, return, nothing to do
  if (ndofH_face.eq.0) return
!
! set order and orientation for all element edge nodes and the face node
  call initiate_order(Type, norder_1)
  do ie=1,nre
    norder_1(ie) = Norder(ie)
  enddo
  norder_1(nre+Iface) = Norder(nre+Iface)
!
#if DEBUG_MODE
  if (iprint.eq.1) then
     write(*,7060) norder_1; call pause
7060 format('dhpfaceH: norder_1 = ',20i4)
  endif
#endif
!
! get face order to find out quadrature information
  call face_order(Type,Iface,Norder, norder_face)
  INTEGRATION=1   ! overintegrate
  call set_2Dint(face_type(Type,Iface),norder_face, &
                 nint,xi_list,wa_list)
  INTEGRATION=0   ! reset
!
! initiate stiffness matrix and load vector
  zbH = ZERO; aaH = 0.d0
!
! loop through integration points
  do l=1,nint
!
!   integration point local face coordinates
    t(1:2) = xi_list(1:2,l)
    wa     = wa_list(l)
!
!   get the corresponding master element coordinates and Jacobian
    call face_param(Type,Iface,t, xi,dxidt)
!
!   compute element H1 shape functions
    call shape3DH(Type,xi,norder_1,Nedge_orient,Nface_orient, &
                  nrdofH,shapH,gradH)
!
!   evaluate reference coordinates of the point as needed by GMP
    nsign = nsign_param(Type,Iface)
    call brefgeom3D(Mdle,xi,Etav,shapH,gradH,nrv,dxidt,nsign, &
                    eta,detadxi,dxideta,rjac,detadt,rn,bjac)
    weight = wa*bjac
!
!   call GMP routines to evaluate physical coordinates and their
!   derivatives wrt reference coordinates
    select case(Iflag)
    case(5);        call prism(No,eta, x,dxdeta)
    case(6);        call  hexa(No,eta, x,dxdeta)
    case(7);        call tetra(No,eta, x,dxdeta)
    case(8);        call pyram(No,eta, x,dxdeta)
    case default
      write(*,*) 'dhpfaceH: Type = ', Type
      call logic_error(ERR_INVALID_VALUE,__FILE__,__LINE__)
    end select
!
!   get Dirichlet data at the point
    call dirichlet(Mdle,x,Icase, &
                   zvalH,zdvalH, zvalE,zdvalE, zvalV,zdvalV)
!
!   evaluate derivatives of Dirichlet data wrt reference coordinates
    zdvalHdeta(1:MAXEQNH,1:3) = ZERO
    do i=1,3
      do j=1,3
        zdvalHdeta(1:MAXEQNH,i) = zdvalHdeta(1:MAXEQNH,i) &
                                + zdvalH(1:MAXEQNH,j)*dxdeta(j,i)
      enddo
    enddo
!
!   remove the contributions from vertices and edges...
!
!   loop through vertex and edge shape functions
    nrdofH = nrdofH - ndofH_face
    do k=1,nrdofH
!
!     compute derivatives of the shape function wrt reference coordinates
      duHdeta(1:3) = gradH(1,k)*dxideta(1,1:3) &
                   + gradH(2,k)*dxideta(2,1:3) &
                   + gradH(3,k)*dxideta(3,1:3)
!
!     subtract...
      do i=1,3
        zdvalHdeta(1:MAXEQNH,i) = zdvalHdeta(1:MAXEQNH,i) &
                                - ZdofH(1:MAXEQNH,k)*duHdeta(i)
      enddo
    enddo
!
!   loop through element face test functions
    do j=1,ndofH_face
      kj = nrdofH + j
!
!     compute gradient wrt reference coordinates
      dvHdeta(1:3) = gradH(1,kj)*dxideta(1,1:3) &
                   + gradH(2,kj)*dxideta(2,1:3) &
                   + gradH(3,kj)*dxideta(3,1:3)
!
!     subtract normal component
      call dot_product(dvHdeta,rn, prod)
      dvHdeta(1:3) = dvHdeta(1:3) - prod*rn(1:3)
!
!     accumulate for the load vector
      zbH(j,1:MAXEQNH) = zbH(j,1:MAXEQNH) &
                       + (zdvalHdeta(1:MAXEQNH,1)*dvHdeta(1) &
                       +  zdvalHdeta(1:MAXEQNH,2)*dvHdeta(2) &
                       +  zdvalHdeta(1:MAXEQNH,3)*dvHdeta(3))*weight
!
!     loop through element face trial functions
      do i=1,ndofH_face
        ki = nrdofH + i
!
!       compute gradient wrt reference coordinates
        duHdeta(1:3) = gradH(1,ki)*dxideta(1,1:3) &
                     + gradH(2,ki)*dxideta(2,1:3) &
                     + gradH(3,ki)*dxideta(3,1:3)
!
!       no need to subtract normal component...
!
!       accumulate for the stiffness matrix
        aaH(j,i) = aaH(j,i) &
                 + (dvHdeta(1)*duHdeta(1) &
                 +  dvHdeta(2)*duHdeta(2) &
                 +  dvHdeta(3)*duHdeta(3))*weight
!
!     end of loop through trial functions
      enddo
!
!   end of loop through test functions
    enddo
!
! end of loop through integration points
  enddo
!
#if DEBUG_MODE
  if (iprint.eq.1) then
    write(*,*) 'dhpfaceH: LOAD VECTOR AND STIFFNESS MATRIX FOR ', &
               'ndofH_face = ',ndofH_face
    do j=1,ndofH_face
      write(*,7015) j, zbH(j,1:MAXEQNH)
      write(*,7016) aaH(j,1:ndofH_face)
    enddo
# if C_MODE
7015    format(i5, 6(2e10.3,2x))
# else
7015    format(i5, 10e12.5)
# endif
7016    format(10e12.5)
  endif
#endif
!
!-----------------------------------------------------------------------
!
! invert the stiffness matrix
  naH = MAXmdlqH
  call dgetrf(ndofH_face,ndofH_face,aaH,naH,ipivH,info)
  if (info.ne.0) then
    write(*,*)'dhpfaceH: H1 DGETRF RETURNED INFO = ',info
    call logic_error(FAILURE,__FILE__,__LINE__)
  endif
!
! back substitute  why double calls ?????????
  zuH(1:ndofH_face,:) = zbH(1:ndofH_face,:)
#if C_MODE
  call zlaswp(MAXEQNH,zuH,naH,1,ndofH_face,ipivH,1)
  duH_real(1:ndofH_face,:) = real(zuH(1:ndofH_face,:))
  duH_imag(1:ndofH_face,:) = aimag(zuH(1:ndofH_face,:))
!
  call dtrsm('L','L','N','U',ndofH_face,MAXEQNH,1.d0,aaH,naH, &
             duH_real,naH)
  call dtrsm('L','U','N','N',ndofH_face,MAXEQNH,1.d0,aaH,naH, &
             duH_real,naH)
!
  call dtrsm('L','L','N','U',ndofH_face,MAXEQNH,1.d0,aaH,naH, &
             duH_imag,naH)
  call dtrsm('L','U','N','N',ndofH_face,MAXEQNH,1.d0,aaH,naH, &
             duH_imag,naH)
!
  zuH(1:ndofH_face,:) &
          = dcmplx(duH_real(1:ndofH_face,:), duH_imag(1:ndofH_face,:))
#else
  call dlaswp(MAXEQNH,zuH,naH,1,ndofH_face,ipivH,1)
  call dtrsm('L','L','N','U',ndofH_face,MAXEQNH,1.d0,aaH,naH, &
             zuH,naH)
  call dtrsm('L','U','N','N',ndofH_face,MAXEQNH,1.d0,aaH,naH, &
             zuH,naH)
#endif
!
#if DEBUG_MODE
  if (iprint.eq.1) then
   write(*,*) 'dhpfaceH: k,zu(k) = '
   do k=1,ndofH_face
     write(*,7015) k,zuH(k,1:MAXEQNH)
   enddo
   call pause
  endif
#endif
!
!  ...save the DOFs, skipping irrelevant entries
  100 continue
!
#if DEBUG_MODE
      if (iprint.eq.1) then
        write(*,*) 'dhpfaceH: ncase = ', ncase
      endif
#endif
!
!  ...initialize global variable counter, and node local variable counter
      ivarH=0 ; nvarH=0
!
!  ...loop through multiple copies of variables
      do j=1,NRCOMS
!
!  .....initiate the BC component counter
        ic=0
!
!  .....loop through physical attributes
        do i=1,NR_PHYSA
!
!  .......loop through components
          do k=1,NR_COMP(i)
!
!  .........if the variable is supported by the node, update the BC component counter
            if (ncase(i).eq.1) ic=ic+1
!
!  .........select the discretization type
            select case(DTYPE(i))
!
!  .........H1 component
            case('contin')
!
!  ...........update global counter
              ivarH = ivarH + 1
!
!  ...........if the variable is supported by the node
              if (ncase(i).eq.1) then
!
!  .............update the node local conter
                nvarH = nvarH + 1
!
!  .............store Dirichlet dof
                if (ibcnd(ic).eq.1) ZnodH(nvarH,1:ndofH_face) = zuH(1:ndofH_face,ivarH)
              endif
            end select
          enddo
        enddo
      enddo
!
#if DEBUG_MODE
      if (iprint.eq.1) call result
#endif
!
      end subroutine dhpfaceH

