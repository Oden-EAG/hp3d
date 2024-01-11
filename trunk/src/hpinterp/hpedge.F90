!
#include "typedefs.h"
!
!-----------------------------------------------------------------------
!> @brief      update H1 geometry dof interpolating GMP reference map
!              using PB interpolation
!!
!! @param[in]  Iflag        - a flag specifying which of the objects the
!!                            edge is on: 5 pris, 6 hexa, 7 tetr, 8 pyra
!! @param[in]  No           - number of a specific object
!! @param[in]  Etav         - reference coordinates of the element vertices
!! @param[in]  Ntype        - element (middle node) type
!! @param[in]  Nedge_orient - edge orientation
!! @param[in]  Nface_orient - face orientation (not used)
!! @param[in]  Norder       - element order
!! @param[in]  Iedge        - edge number
!! @param[in]  Xnod         - geometry dof for the element (vertex values)
!!
!! @param[out] Xdof         - geometry dof for the edge
!!
!> @date       Feb 2023
!-----------------------------------------------------------------------
  subroutine hpedge(Mdle,Iflag,No,Etav,Ntype, &
                    Nedge_orient,Nface_orient,Norder,Iedge,&
                    Xnod, Xdof)
!
  use parameters
  use element_data
  implicit none
!
! ** Arguments
!-----------------------------------------------------------------------
!
  integer,                         intent(in)  :: Iflag,No,Mdle
  integer,                         intent(in)  :: Iedge
  real(8), dimension(3,8),         intent(in)  :: Etav
  integer,                         intent(in)  :: Ntype
  integer, dimension(12),          intent(in)  :: Nedge_orient
  integer, dimension(6),           intent(in)  :: Nface_orient
  integer, dimension(19),          intent(in)  :: Norder
  real(8), dimension(3,MAXbrickH), intent(in)  :: Xnod
  real(8), dimension(3,*),         intent(out) :: Xdof
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
! work space for shape3H
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
! work space for linear solvers
  integer                               :: naH,info
  real(8), dimension(MAXP-1,MAXP-1)     :: aaH
  integer, dimension(MAXP-1)            :: ipivH
!
! load vector and solution
  real(8), dimension(MAXP-1,3)          :: bb,uu
!
! misc work space
  integer :: nrv,nre,nrf,i,j,k,kj,ki,&
             ndofH_edge,ndofE_edge,ndofV_edge,ndofQ_Edge,iflag1
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
     write(*,7010) Mdle,Iflag,No,Iedge,S_Type(Ntype)
7010 format('hpedge: Mdle,Iflag,No,Iedge,Type = ',4i4,2x,a4)
     write(*,7020) Etav(1:3,1:nrv)
7020 format('        Etav = ',8(3f6.2,1x))
     write(*,7030) Nedge_orient(1:nre)
7030 format('        Nedge_orient = ',12i2)
     write(*,7040) Nface_orient(1:nrf)
7040 format('        Nface_orient = ',6i2)
     write(*,7050) Norder(1:nre+nrf+1)
7050 format('        Norder = ',19i4)
     call pause
  endif
#endif
!
! # of edge dof
  call ndof_nod(MEDG,norder(Iedge), &
                ndofH_edge,ndofE_edge,ndofV_edge,ndofQ_edge)
!
! if # of dof is zero, return, nothing to do
  if (ndofH_edge.eq.0) return
!
! set order and orientation for the edge node
  call initiate_order(Ntype, norder_1)
  norder_1(Iedge) = Norder(Iedge)
!
! 1D integration rule
  call set_1Dint(Norder(Iedge), nint,xi_list,wa_list)
!
! initialize
  aaH = 0.d0; bb = 0.d0
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
      write(*,*) 'hpedge: Type,Iflag = ', S_Type(Ntype),Iflag
      call logic_error(ERR_INVALID_VALUE,__FILE__,__LINE__)
    endselect
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
        dxdeta(1:3,i) = dxdeta(1:3,i) &
                      - Xnod(1:3,k)*duHdeta(i)
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
      bb(j,1:3) = bb(j,1:3) &
                + (dxdeta(1:3,1)*dvHdeta(1) &
                +  dxdeta(1:3,2)*dvHdeta(2) &
                +  dxdeta(1:3,3)*dvHdeta(3))*weight
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
    write(*,*) 'hpedge: LOAD VECTOR AND STIFFNESS MATRIX FOR ', &
               'ndofH_edge = ',ndofH_edge
    do j=1,ndofH_edge
      write(*,7015) j, bb(j,1:3)
      write(*,7016) aaH(j,1:ndofH_edge)
    enddo
7015    format(i5, 10e12.5)
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
  uu(1:ndofH_edge,:) = bb(1:ndofH_edge,:)
!
! apply pivots to load vector
  call dlaswp(3,uu,naH,1,ndofH_edge,ipivH,1)
!
! triangular solves
  call dtrsm('L','L','N','U',ndofH_edge,3,1.d0,aaH,naH, uu,naH)
  call dtrsm('L','U','N','N',ndofH_edge,3,1.d0,aaH,naH, uu,naH)
!
#if DEBUG_MODE
  if (iprint.eq.1) then
   write(*,*) 'hpedge: k,uu(k) = '
   do k=1,ndofH_edge
     write(*,7015) k,uu(k,1:3)
   enddo
   call pause
  endif
#endif
!
! store geometry dof
  do i=1,3
    Xdof(i,1:ndofH_edge) = uu(1:ndofH_edge,i)
  enddo
!
!
  end subroutine hpedge
