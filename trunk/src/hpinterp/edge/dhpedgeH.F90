!-----------------------------------------------------------------------
!
!> @brief    update H1 edge dof interpolating H1 Dirichlet data using
!            PB interpolation
!!
!! @param[in]  Iflag        - a flag specifying which of the objects the
!!                            edge is on: 5 pris, 6 hexa, 7 tetr, 8 pyra
!! @param[in]  No           - number of a specific object
!! @param[in]  Etav         - reference coordinates of the element vertices
!! @param[in]  Ntype        - element (middle node) type
!! @param[in]  Icase        - the edge node case
!! @param[in]  Bcond        - the edge node BC flag
!! @param[in]  Nedge_orient - edge orientation
!! @param[in]  Nface_orient - face orientation (not used)
!! @param[in]  Norder       - element order
!! @param[in]  Iedge        - edge number
!! @param[in]  ZdofH        - H1 dof for the element (vertex values)
!!
!! @param[in,out] ZnodH     - H1 dof for the edge
!!
!> @date Feb 2023
!-----------------------------------------------------------------------
!
#include "typedefs.h"
  subroutine dhpedgeH(Mdle,Iflag,No,Etav,Ntype,Icase,Bcond,&
                      Nedge_orient,Nface_orient,Norder,Iedge,&
                      ZdofH, ZnodH)
!
  use control
  use parameters
  use physics
  use element_data
  implicit none
!
! ** Arguments
!-----------------------------------------------------------------------
!
  integer, intent(in)    :: Iflag,No,Mdle
  integer, intent(in)    :: Icase,Bcond,Iedge
  real(8), intent(in)    :: Etav(3,8)
  integer, intent(in)    :: Ntype
  integer, intent(in)    :: Nedge_orient(12)
  integer, intent(in)    :: Nface_orient(6)
  integer, intent(in)    :: Norder(19)
  VTYPE,   intent(in)    :: ZdofH(MAXEQNH,MAXbrickH)
  VTYPE,   intent(inout) :: ZnodH(NRCOMS*NREQNH(Icase),*)
!
! ** Locals
!-----------------------------------------------------------------------
!
! quadrature
  integer                               :: l,nint
  real(8), dimension(MAX_NINT1)         :: xi_list
  real(8), dimension(MAX_NINT1)         :: wa_list
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
! geometry
  real(8)                               :: t,rjac,bjac,prod
  real(8), dimension(3)                 :: xi,eta,x
  real(8), dimension(3)                 :: dxidt,detadt,rt
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
! decoded case and BC flags for the edge node
  integer, dimension(NR_PHYSA)          :: ncase
  integer, dimension(NRINDEX)           :: ibcnd
!
! work space for linear solvers
  integer                               :: naH,info
  real(8), dimension(MAXP-1,MAXP-1)     :: aaH
  integer, dimension(MAXP-1)            :: ipivH
!
! load vector and solution
  VTYPE,   dimension(MAXP-1,MAXEQNH)    :: zbH,zuH
#if C_MODE
  real(8), dimension(MAXP-1,MAXEQNH)    :: duH_real,duH_imag
#endif
!
! misc work space
  integer :: nrv,nre,nrf,i,j,k,ivarH,nvarH,kj,ki,&
             ndofH_edge,ndofE_edge,ndofV_edge,ndofQ_Edge,iflag1,ic
!
  logical :: is_homD
!
#if DEBUG_MODE
  integer :: iprint
  iprint=0
#endif
!
!----------------------------------------------------------------------
!
  nrv = nvert(Ntype); nre = nedge(Ntype); nrf = nface(Ntype)
!
#if DEBUG_MODE
  if (iprint.eq.1) then
     write(*,7010) Mdle,Iflag,No,Icase,Iedge,S_Type(Ntype)
7010 format('dhpedgeH: Mdle,Iflag,No,Icase,Iedge,Type = ',5i4,2x,a4)
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
! # of edge dof
  call ndof_nod(MEDG,norder(Iedge), &
                ndofH_edge,ndofE_edge,ndofV_edge,ndofQ_edge)
!
! check if a homogeneous Dirichlet node
  call homogenD(CONTIN,Icase,Bcond, is_homD,ncase,ibcnd)
  if (is_homD) then
    zuH = ZERO
    go to 100
  endif
!
! if # of dof is zero, return, nothing to do
  if (ndofH_edge.eq.0) return
!
! set order and orientation for the edge node
  call initiate_order(Ntype, norder_1)
  norder_1(Iedge) = Norder(Iedge)
!
! 1D integration rule
  INTEGRATION=1   ! overintegrate
  call set_1Dint(Norder(Iedge), nint, xi_list, wa_list)
  INTEGRATION=0   ! reset
!
! initialize
  aaH = 0.d0; zbH = ZERO
!
! loop over integration points
  do l=1,nint
!
!   Gauss point and weight
    t  = xi_list(l)
    wa = wa_list(l)
!
!   determine edge parameterization for line integral
    call edge_param(Ntype,Iedge,t, xi,dxidt)
!
!   compute element H1 shape functions
    call shape3DH(Ntype,xi,norder_1,Nedge_orient,Nface_orient, &
                  nrdofH,shapH,gradH)
!
!   compute reference geometry
    call refgeom3D(Mdle,xi,Etav,shapH,gradH,nrv, &
                   eta,detadxi,dxideta,rjac,iflag1)
!
!   compute unit tangent vector
    detadt(1:3) = detadxi(1:3,1)*dxidt(1) &
                + detadxi(1:3,2)*dxidt(2) &
                + detadxi(1:3,3)*dxidt(3)
    call norm(detadt, bjac)
    rt(1:3) = detadt(1:3)/bjac
    weight = wa*bjac
!
!   call GMP routines to evaluate physical coordinates and their
!   derivatives wrt reference coordinates
    select case(Iflag)
    case(5) ; call prism(No,eta, x,dxdeta)
    case(6) ; call  hexa(No,eta, x,dxdeta)
    case(7) ; call tetra(No,eta, x,dxdeta)
    case(8) ; call pyram(No,eta, x,dxdeta)
    case default
      write(*,*) 'dhpedgeH: Type,Iflag = ', S_Type(Ntype),Iflag
      call logic_error(ERR_INVALID_VALUE,__FILE__,__LINE__)
    endselect
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
!   remove the contributions from vertices
!
!   loop through vertex shape functions
    nrdofH = nrdofH - ndofH_edge
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
!   loop through element edge test functions
    do j=1,ndofH_edge
      kj = nrdofH + j
!
!     compute gradient wrt reference coordinates
      dvHdeta(1:3) = gradH(1,kj)*dxideta(1,1:3) &
                   + gradH(2,kj)*dxideta(2,1:3) &
                   + gradH(3,kj)*dxideta(3,1:3)
!
!     leave tangential component only
      call dot_product(dvHdeta,rt, prod)
      dvHdeta(1:3) = prod*rt(1:3)
!
!     accumulate for the load vector
      zbH(j,1:MAXEQNH) = zbH(j,1:MAXEQNH) &
                       + (zdvalHdeta(1:MAXEQNH,1)*dvHdeta(1) &
                       +  zdvalHdeta(1:MAXEQNH,2)*dvHdeta(2) &
                       +  zdvalHdeta(1:MAXEQNH,3)*dvHdeta(3))*weight
!
!     loop through element face trial functions
      do i=1,ndofH_edge
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
    write(*,*) 'dhpedgeH: LOAD VECTOR AND STIFFNESS MATRIX FOR ', &
               'ndofH_edge = ',ndofH_edge
    do j=1,ndofH_edge
      write(*,7015) j, zbH(j,1:MAXEQNH)
      write(*,7016) aaH(j,1:ndofH_edge)
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
! projection matrix leading dimension (maximum number of 1D bubbles)
  naH=MAXP-1
!
! over-write aaH with its LU factorization
  call dgetrf(ndofH_edge,ndofH_edge,aaH,naH,ipivH,info)
!
! check that factorization was successful
  if (info /= 0) then
    write(*,*)'dhpedge: H1 DGETRF RETURNED INFO = ',info
    call logic_error(FAILURE,__FILE__,__LINE__)
  endif
!
! solve linear system
!
! copy load vector
  zuH(1:ndofH_edge,:) = zbH(1:ndofH_edge,:)
!
#if C_MODE
!
! apply pivots to load vector
  call zlaswp(MAXEQNH,zuH,naH,1,ndofH_edge,ipivH,1)
!
! compute real part
  duH_real(1:ndofH_edge,:) = real( zuH(1:ndofH_edge,:))
!
! compute imaginary part
  duH_imag(1:ndofH_edge,:) = aimag(zuH(1:ndofH_edge,:))
!
! triangular solves for real part
  call dtrsm('L','L','N','U',ndofH_edge,MAXEQNH,1.d0,aaH,naH, &
             duH_real,naH)
  call dtrsm('L','U','N','N',ndofH_edge,MAXEQNH,1.d0,aaH,naH, &
             duH_real,naH)
!
! triangular solves for imaginary part
  call dtrsm('L','L','N','U',ndofH_edge,MAXEQNH,1.d0,aaH,naH, &
             duH_imag,naH)
  call dtrsm('L','U','N','N',ndofH_edge,MAXEQNH,1.d0,aaH,naH, &
             duH_imag,naH)
!
! combine real and imaginary parts by forcing type to DOUBLE
! precision complex
  zuH(1:ndofH_edge,:) &
   = dcmplx(duH_real(1:ndofH_edge,:), duH_imag(1:ndofH_edge,:))
!
#else
!
! apply pivots to load vector
  call dlaswp(MAXEQNH,zuH,naH,1,ndofH_edge,ipivH,1)
!
! triangular solves
  call dtrsm('L','L','N','U',ndofH_edge,MAXEQNH,1.d0,aaH,naH, zuH,naH)
  call dtrsm('L','U','N','N',ndofH_edge,MAXEQNH,1.d0,aaH,naH, zuH,naH)
!
#endif
!
#if DEBUG_MODE
  if (iprint.eq.1) then
   write(*,*) 'dhpedgeH: k,zu(k) = '
   do k=1,ndofH_edge
     write(*,7015) k,zuH(k,1:MAXEQNH)
   enddo
   call pause
  endif
#endif
!
!  ...save the DOFs, skipping irrelevant entries
 100  continue
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
!  .......loop through components of physical attribute
          do k=1,NR_COMP(i)
!
!  .........if the variable is supported by the node, update the BC component counter
            if (ncase(i).eq.1) ic=ic+1
!
!  .........select the discretization type
            select case(D_TYPE(i))
!
!  .........H1 component
            case(CONTIN)
!
!  ...........update global counter
              ivarH = ivarH + 1
!
!  ...........if the variable is supported by the node
              if (ncase(i).eq.1) then
!
!  .............update the node local counter
                nvarH = nvarH + 1
!
!  .............do not write dof if physics attribute is deactivated
                if (.not. PHYSAm(i)) exit
!
!  .............store Dirichlet dof
                if (ibcnd(ic).eq.1) ZnodH(nvarH,1:ndofH_edge) = zuH(1:ndofH_edge,ivarH)
!
              endif
            end select
          enddo
        enddo
      enddo
!
!
      end subroutine dhpedgeH
