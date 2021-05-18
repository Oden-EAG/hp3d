!-----------------------------------------------------------------------
!
!> Purpose : update H(curl) edge dof interpolating H1 Dirichlet data
!            using PB interpolation
!!
!! @param[in]  Iflag        - a flag specifying which of the objects the
!!                            edge is on: 5 pris, 6 hexa, 7 tetr, 8 pyra
!! @param[in]  No           - number of a specific object
!! @param[in]  Etav         - reference coordinates of the element vertices
!! @param[in]  Type         - element (middle node) type
!! @param[in]  Icase        - the edge node case
!! @param[in]  Bcond        - the edge node BC flag
!! @param[in]  Nedge_orient - edge orientation
!! @param[in]  Nface_orient - face orientation (not used)
!! @param[in]  Norder       - element order
!! @param[in]  Iedge        - edge number
!!
!! @param[out] ZnodE        - H(curl) dof for the edge
!!
!> @date Mar 17
!-----------------------------------------------------------------------
!
#include "typedefs.h"
  subroutine dhpedgeE(Mdle,Iflag,No,Etav,Type,Icase,Bcond,&
                      Nedge_orient,Nface_orient,Norder,Iedge,&
                      ZnodE)
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
  integer,          intent(in)    :: Iflag,No,Mdle
  integer,          intent(in)    :: Icase,Bcond,Iedge
  real(8),          intent(in)    :: Etav(3,8)
  character(len=4), intent(in)    :: Type
  integer,          intent(in)    :: Nedge_orient(12)
  integer,          intent(in)    :: Nface_orient(6)
  integer,          intent(in)    :: Norder(19)
  VTYPE,            intent(inout) :: ZnodE(NRCOMS*NREQNE(Icase),*)
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
! work space for shape3DE
  integer                               :: nrdofE
  real(8), dimension(3,MAXbrickE)       :: shapE
  real(8), dimension(3,MAXbrickE)       :: curlE
!
! H(curl) shape functions in reference coordinates
  real(8), dimension(3)                 :: uEeta,vEeta
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
! Dirichlet data in reference coordinates
  VTYPE :: zvalEeta(3,MAXEQNE)
!
! decoded case and BC flag for the face node
  integer, dimension(NR_PHYSA)          :: ncase
  integer, dimension(NRINDEX)           :: ibcnd
!
! work space for linear solvers
  integer                               :: naE,info
  real(8), dimension(MAXP,MAXP)         :: aaE
  integer, dimension(MAXP)              :: ipivE
!
! load vector and solution
  VTYPE,   dimension(MAXP,MAXEQNE)      :: zbE,zuE
#if C_MODE
  real(8), dimension(MAXP,MAXEQNE)      :: uE_real,uE_imag
#endif
!
! misc work space
  integer :: iprint,nrv,nre,nrf,i,j,k,ivarE,nvarE,kj,ki,&
             ndofH_edge,ndofE_edge,ndofV_edge,ndofQ_Edge,iflag1,ic
!
!----------------------------------------------------------------------
!
  nrv = nvert(Type); nre = nedge(Type); nrf = nface(Type)
!
  iprint = 0
#if DEBUG_MODE
  if (iprint.eq.1) then
     write(*,7010) Mdle,Iflag,No,Icase,Iedge,Type
7010 format('dhpedgeE: Mdle,Iflag,No,Icase,Iedge,Type = ',5i4,2x,a4)
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
  call ndof_nod('medg',norder(Iedge), &
                ndofH_edge,ndofE_edge,ndofV_edge,ndofQ_edge)
!
! # of dof cannot be zero
  if (ndofE_edge.le.0) then
    write(*,*) 'dhpedgeE: ndofE_edge = ',ndofE_edge
    stop 1
  endif
!
! set order and orientation for the edge node
  call initiate_order(Type, norder_1)
  norder_1(Iedge) = Norder(Iedge)
!
! 1D integration rule
  INTEGRATION=1   ! overintegrate
  call set_1Dint(Norder(Iedge), nint, xi_list, wa_list)
  INTEGRATION=0   ! reset
!
! initialize
  aaE = 0.d0; zbE = ZERO
!
! loop over integration points
  do l=1,nint
!
!   Gauss point and weight
    t  = xi_list(l)
    wa = wa_list(l)
!
!   determine edge parameterization for line integral
    call edge_param(Type,Iedge,t, xi,dxidt)
!
!   compute element H1 shape functions (for geometry)
    call shape3DH(Type,xi,norder_1,Nedge_orient,Nface_orient, &
                  nrdofH,shapH,gradH)
!
!   compute element Hcurl shape functions
    call shape3DE(Type,xi,norder_1,Nedge_orient,Nface_orient, &
                  nrdofE,shapE,curlE)
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
      write(*,*) 'dhpedgeE: Type,Iflag = ', Type,Iflag
      call logic_error(ERR_INVALID_VALUE,__FILE__,__LINE__)
    endselect
!
!   get Dirichlet data at the point
    call dirichlet(Mdle,x,Icase, &
                   zvalH,zdvalH, zvalE,zdvalE, zvalV,zdvalV)
!
!   transform Dirichlet data to reference coordinates
    zvalEeta(1:3,1:MAXEQNE) = ZERO
    do i=1,3
      do j=1,3
        zvalEeta(i,1:MAXEQNE) = zvalEeta(i,1:MAXEQNE) &
                              + zvalE(j,1:MAXEQNE)*dxdeta(j,i)
      enddo
    enddo
!
!   loop through element edge test functions
    do j=1,ndofE_edge
      kj = Iedge-1 + j
!
!     compute test shape function in reference coordinates
      vEeta(1:3) = shapE(1,kj)*dxideta(1,1:3) &
                 + shapE(2,kj)*dxideta(2,1:3) &
                 + shapE(3,kj)*dxideta(3,1:3)
!
!     leave tangential component only
      call dot_product(vEeta,rt, prod)
      vEeta(1:3) = prod*rt(1:3)
!
!     accumulate for the load vector
      zbE(j,1:MAXEQNE) = zbE(j,1:MAXEQNE) &
                       + (zvalEeta(1,1:MAXEQNE)*vEeta(1) &
                       +  zvalEeta(2,1:MAXEQNE)*vEeta(2) &
                       +  zvalEeta(3,1:MAXEQNE)*vEeta(3))*weight
!
!     loop through element face trial functions
      do i=1,ndofE_edge
        ki = Iedge-1 + i
!
!       compute gradient wrt reference coordinates
        uEeta(1:3) = shapE(1,ki)*dxideta(1,1:3) &
                   + shapE(2,ki)*dxideta(2,1:3) &
                   + shapE(3,ki)*dxideta(3,1:3)
!
!       no need to subtract normal component...
!
!       accumulate for the stiffness matrix
        aaE(j,i) = aaE(j,i) &
                 + (vEeta(1)*uEeta(1) &
                 +  vEeta(2)*uEeta(2) &
                 +  vEeta(3)*uEeta(3))*weight
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
    write(*,*) 'dhpedgeE: LOAD VECTOR AND STIFFNESS MATRIX FOR ', &
               'ndofE_edge = ',ndofE_edge
    do j=1,ndofE_edge
      write(*,7015) j, zbE(j,1:MAXEQNE)
      write(*,7016) aaE(j,1:ndofE_edge)
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
  naE=MAXP
!
! over-write aaE with its LU factorization
  call dgetrf(ndofE_edge,ndofE_edge,aaE,naE,ipivE,info)
!
! check that factorization was successful
  if (info /= 0) then
    write(*,*)'dhpedge: H(curl) DGETRF RETURNED INFO = ',info
    call logic_error(FAILURE,__FILE__,__LINE__)
  endif
!
! solve linear system
!
! copy load vector
  zuE(1:ndofE_edge,:) = zbE(1:ndofE_edge,:)
!
#if C_MODE
!
! apply pivots to load vector
  call zlaswp(MAXEQNE,zuE,naE,1,ndofE_edge,ipivE,1)
!
! compute real part
  uE_real(1:ndofE_edge,:) = real( zuE(1:ndofE_edge,:))
!
! compute imaginary part
  uE_imag(1:ndofE_edge,:) = aimag(zuE(1:ndofE_edge,:))
!
! triangular solves for real part
  call dtrsm('L','L','N','U',ndofE_edge,MAXEQNE,1.d0,aaE,naE,uE_real,naE)
  call dtrsm('L','U','N','N',ndofE_edge,MAXEQNE,1.d0,aaE,naE,uE_real,naE)
!
! triangular solves for imaginary part
  call dtrsm('L','L','N','U',ndofE_edge,MAXEQNE,1.d0,aaE,naE,uE_imag,naE)
  call dtrsm('L','U','N','N',ndofE_edge,MAXEQNE,1.d0,aaE,naE,uE_imag,naE)
!
! combine real and imaginary parts by forcing type to DOUBLE
! precision complex
  zuE(1:ndofE_edge,:) = dcmplx(uE_real(1:ndofE_edge,:), uE_imag(1:ndofE_edge,:))
!
#else
!
! apply pivots to load vector
  call dlaswp(MAXEQNE,zuE,naE,1,ndofE_edge,ipivE,1)
!
! triangular solves
  call dtrsm('L','L','N','U',ndofE_edge,MAXEQNE,1.d0,aaE,naE, zuE,naE)
  call dtrsm('L','U','N','N',ndofE_edge,MAXEQNE,1.d0,aaE,naE, zuE,naE)
!
#endif
!
#if DEBUG_MODE
  if (iprint.eq.1) then
   write(*,*) 'dhpedgeE: k,zuE(k) = '
   do k=1,ndofE_edge
     write(*,7015) k,zuE(k,1:MAXEQNE)
   enddo
   call pause
  endif
#endif
!
!  ...save dof's, skipping irrelevant entries
!
!  ...decode the case and the BC flag
      call decod(Icase,2,NR_PHYSA, ncase)
      call decod(Bcond,2,NRINDEX,  ibcnd)
!
!  ...initialize global variable counter, and node local variable counter
      ivarE=0 ; nvarE=0
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
            select case(DTYPE(i))
!
!  .........H(curl) component
            case('tangen')
!
!  ...........update global counter
              ivarE = ivarE + 1
!
!  ...........if the variable is supported by the node
              if (ncase(i).eq.1) then
!
!  .............update node local conter
                nvarE = nvarE + 1
!
!  .............store Dirichlet dof
                if (ibcnd(ic).eq.1) ZnodE(nvarE,1:ndofE_edge) = zuE(1:ndofE_edge,ivarE)
!
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
      end subroutine dhpedgeE
