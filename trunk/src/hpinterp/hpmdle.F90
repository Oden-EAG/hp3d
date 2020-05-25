!
#include "implicit_none.h"
!
!-----------------------------------------------------------------------
!> Purpose : determine middle node geometry dof interpolating GMP map using
!            PB interpolation
!  NOTE:     the interpolation (projection) is done in the reference space
!!
!! @param[in]  Iflag        - a flag specifying which of the objects the
!!                            face is on: 5 pris, 6 hexa, 7 tetr, 8 pyra
!! @param[in]  No           - number of a specific object
!! @param[in]  Etav         - reference coordinates of the element vertices
!! @param[in]  Type         - element (middle node) type
!! @param[in]  Nedge_orient - edge orientation
!! @param[in]  Nface_orient - face orientation
!! @param[in]  Norder       - element order
!! @param[in]  Xnod         - geometry dof for the element (vertex,edge
!!                            and face  values)
!!
!! @param[out] Xdof         - geometry dof for the middle node
!-----------------------------------------------------------------------
#include "implicit_none.h"
  subroutine hpmdle(Mdle,Iflag,No,Etav,Type, &
                    Nedge_orient,Nface_orient,Norder, &
                    Xnod, Xdof)
  use parameters
  use physics
  use element_data
  implicit none
!
! ** Arguments
!-----------------------------------------------------------------------
  integer,                         intent(in)  :: Iflag,No,Mdle
  real(8),  dimension(3,8),        intent(in)  :: Etav
  character(len=4),                intent(in)  :: Type
  integer, dimension(12),          intent(in)  :: Nedge_orient
  integer, dimension(6),           intent(in)  :: Nface_orient
  integer, dimension(19),          intent(in)  :: Norder
!
  real(8), dimension(3,MAXbrickH), intent(in)  :: Xnod
  real(8), dimension(3,*),         intent(out) :: Xdof
!
! ** Locals
!-----------------------------------------------------------------------
!
! quadrature
  integer                               :: l,nint
  real(8), dimension(3, MAX_NINT3)      :: xi_list
  real(8), dimension(   MAX_NINT3)      :: wa_list
  real(8)                               :: wa, weight
!
! work space for shape3H
  integer                               :: nrdofH
  real(8), dimension(MAXbrickH)         :: shapH
  real(8), dimension(3,MAXbrickH)       :: gradH
!
! derivatives of a shape function wrt reference coordinates
  real(8), dimension(3)                 :: duHdeta,dvHdeta
!
! geometry
  real(8)                               :: rjac
  real(8), dimension(3)                 :: xi,eta,x
  real(8), dimension(3,3)               :: detadxi,dxideta,dxdeta
!
! work space for linear solvers
  integer                               :: naH,info
  real(8), dimension(MAXMdlbH,MAXMdlbH) :: aaH
  integer, dimension(MAXMdlbH)          :: ipivH
!
! load vector and solution
  real(8), dimension(MAXMdlbH,3)        :: bb,uu
!
! misc work space
  integer :: iprint,nrv,nre,nrf,i,j,k,kj,ki,&
             ndofH_mdle,ndofE_mdle,ndofV_mdle,ndofQ_Mdle,iflag1
!
!-----------------------------------------------------------------------
  iprint=0
!
  nrv = nvert(Type); nre = nedge(Type); nrf = nface(Type)
  if (iprint.eq.1) then
     write(*,7010) Mdle,Iflag,No, Type
7010 format('hpmdle: Mdle,Iflag,No, Type = ',3i4,2x,a4)
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
!
! determine # of dof for the mdle node
  call ndof_nod(Type,Norder(nre+nrf+1), &
                ndofH_mdle,ndofE_mdle,ndofV_mdle,ndofQ_mdle)
!
! if # of dof is zero, return, nothing to do
  if (ndofH_mdle.eq.0) return
!
! get quadrature
  call set_3Dint(Type,Norder, nint,xi_list,wa_list)
!
! initiate stiffness matrix and load vector
  bb = ZERO; aaH = 0.d0
!
! loop through integration points
  do l=1,nint
!
!   integration point local face coordinates
    xi(1:3) = xi_list(1:3,l)
    wa      = wa_list(l)
!
!   compute element H1 shape functions
    call shape3H(Type,xi,Norder,Nedge_orient,Nface_orient, &
                 nrdofH,shapH,gradH)
!
!   compute reference geometry
    call refgeom3D(Mdle,xi,Etav,shapH,gradH,nrv, &
                   eta,detadxi,dxideta,rjac,iflag1)
    weight = wa*rjac
!
!   call GMP routines to evaluate physical coordinates and their
!   derivatives wrt reference coordinates
    select case(Iflag)
    case(5);        call prism(No,eta, x,dxdeta)
    case(6);        call  hexa(No,eta, x,dxdeta)
    case(7);        call tetra(No,eta, x,dxdeta)
    case(8);        call pyram(No,eta, x,dxdeta)
    case default
      write(*,*) 'hpmdle: Type = ', Type
      call logic_error(ERR_INVALID_VALUE,__FILE__,__LINE__)
    end select
!
!   remove the contributions from vertices, edges and faces
!
!   loop through vertex and edge shape functions
    nrdofH = nrdofH - ndofH_mdle
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
!   loop through element face test functions
    do j=1,ndofH_mdle
      kj = nrdofH + j
!
!     compute gradient wrt reference coordinates
      dvHdeta(1:3) = gradH(1,kj)*dxideta(1,1:3) &
                   + gradH(2,kj)*dxideta(2,1:3) &
                   + gradH(3,kj)*dxideta(3,1:3)
!
!     accumulate for the load vector
      bb(j,1:3) = bb(j,1:3) &
                + (dxdeta(1:3,1)*dvHdeta(1) &
                +  dxdeta(1:3,2)*dvHdeta(2) &
                +  dxdeta(1:3,3)*dvHdeta(3))*weight
!
!     loop through element face trial functions
      do i=1,ndofH_mdle
        ki = nrdofH + i
!
!       compute gradient wrt reference coordinates
        duHdeta(1:3) = gradH(1,ki)*dxideta(1,1:3) &
                     + gradH(2,ki)*dxideta(2,1:3) &
                     + gradH(3,ki)*dxideta(3,1:3)
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
  if (iprint.eq.1) then
    write(*,*) 'hpmdle: LOAD VECTOR AND STIFFNESS MATRIX FOR ', &
               'ndofH_mdle = ',ndofH_mdle
    do j=1,ndofH_mdle
      write(*,7015) j, bb(j,1:3)
      write(*,7016) aaH(j,1:ndofH_mdle)
    enddo
7015    format(i5, 10e12.5)
7016    format(10e12.5)
  endif
!
!-----------------------------------------------------------------------
!
! invert the stiffness matrix
  naH = MAXmdlbH
  call dgetrf(ndofH_mdle,ndofH_mdle,aaH,naH,ipivH,info)
  if (info.ne.0) then
    write(*,*)'hpmdle:  DGETRF RETURNED INFO = ',info
    call logic_error(FAILURE,__FILE__,__LINE__)
  endif
!
! back substitute  why double calls ?????????
  uu(1:ndofH_mdle,:) = bb(1:ndofH_mdle,:)
  call dlaswp(3,uu(1:ndofH_mdle,:),naH,1,ndofH_mdle,ipivH,1)
  call dtrsm('L','L','N','U',ndofH_mdle,3,1.d0,aaH,naH, &
             uu,naH)
  call dtrsm('L','U','N','N',ndofH_mdle,3,1.d0,aaH,naH, &
             uu,naH)
!
  if (iprint.eq.1) then
   write(*,*) 'hpmdle: k,uu(k) = '
   do k=1,ndofH_mdle
     write(*,7015) k,uu(k,1:3)
   enddo
   call pause
  endif
!
! save the dof
  do i=1,3
    Xdof(i,1:ndofH_mdle) = uu(1:ndofH_mdle,i)
  enddo
!
end subroutine hpmdle

