!
#include "typedefs.h"
!
!-----------------------------------------------------------------------
!> @brief    determine H(div) face dof interpolating H(div) Dirichlet data
!!           using PB interpolation; NOTE: the interpolation (projection)
!!           is done in the reference space
!!
!! @param[in]  Iflag        - a flag specifying which of the objects the
!!                            face is on: 5 pris, 6 hexa, 7 tetr, 8 pyra
!! @param[in]  No           - number of a specific object
!! @param[in]  Etav         - reference coordinates of the element vertices
!! @param[in]  Ntype        - element (middle node) type
!! @param[in]  Icase        - the face node case
!! @param[in]  Bcond        - the edge node BC flag
!! @param[in]  Nedge_orient - edge orientation, never used
!! @param[in]  Nface_orient - face orientation
!! @param[in]  Norder       - element order
!! @param[in]  Iface        - face number
!!
!! @param[in,out] ZnodV     - H(div) dof for the face
!!
!> @date Feb 2023
!-----------------------------------------------------------------------
  subroutine dhpfaceV(Mdle,Iflag,No,Etav,Ntype,Icase,Bcond, &
                      Nedge_orient,Nface_orient,Norder,Iface, &
                      ZnodV)
  use control
  use parameters
  use physics
  use element_data
  implicit none
!
! ** Arguments
!-----------------------------------------------------------------------
  integer,                 intent(in)  :: Iflag,No,Mdle
  integer,                 intent(in)  :: Icase,Bcond,Iface
  real(8), dimension(3,8), intent(in)  :: Etav
  integer,                 intent(in)  :: Ntype
  integer, dimension(12),  intent(in)  :: Nedge_orient
  integer, dimension(6),   intent(in)  :: Nface_orient
  integer, dimension(19),  intent(in)  :: Norder
!
  VTYPE,   dimension(NRCOMS*NREQNV(Icase),*), intent(inout) :: ZnodV
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
!
  real(8), dimension(MAXbrickH)         :: shapH
  real(8), dimension(3,MAXbrickH)       :: gradH
!
! work space for shape3DV
  integer                               :: nrdofV
  real(8), dimension(3,MAXbrickV)       :: shapV
  real(8), dimension(MAXbrickV)         :: divV
!
! H(div) test and trial shape function in reference coordinates
  real(8), dimension(3)                 :: uVeta,vVeta
!
! dot product
  real(8)                               :: prod
!
! geometry
  real(8)                               :: rjac,bjac,rjacdxdeta
  real(8), dimension(2)                 :: t
  real(8), dimension(3)                 :: xi,eta,rn,x
  real(8), dimension(3,2)               :: dxidt,detadt
  real(8), dimension(3,3)               :: detadxi,dxideta,dxdeta,detadx
!
! Dirichlet BC data at a point
  VTYPE :: zvalH(  MAXEQNH), zdvalH(  MAXEQNH,3), &
           zvalE(3,MAXEQNE), zdvalE(3,MAXEQNE,3), &
           zvalV(3,MAXEQNV), zdvalV(3,MAXEQNV,3)
!
! H(div) Dirichlet data in reference coordinates
  VTYPE :: zvalVeta(3,MAXEQNV)
!
! work space for linear solvers
  integer                               :: naV,info
  real(8), dimension(MAXmdlqV,MAXmdlqV) :: aaV
  integer, dimension(MAXmdlqV)          :: ipivV
!
! load vector and solution
  VTYPE,   dimension(MAXmdlqV,MAXEQNV)  :: zbV,zuV
#if C_MODE
  real(8), dimension(MAXmdlqV,MAXEQNV)  :: uV_real,uV_imag
#endif
!
! decoded case and BC flag for the face node
  integer, dimension(NR_PHYSA)          :: ncase
  integer, dimension(NRINDEX)           :: ibcnd
!
! misc work space
  integer :: nrv,nre,nrf,nsign,nflag,i,j,k,ivarV,nvarV,kj,ki,&
             ndofH_face,ndofE_face,ndofV_face,ndofQ_face,ic
!
  logical :: is_homD
!
#if DEBUG_MODE
  integer :: iprint
  iprint=0
#endif
!
!-----------------------------------------------------------------------
!
  nrv = nvert(Ntype); nre = nedge(Ntype); nrf = nface(Ntype)
!
#if DEBUG_MODE
  if (iprint.eq.1) then
     write(*,7010) Mdle,Iflag,No,Icase,Iface,S_Type(Ntype)
7010 format('dhpfaceV: Mdle,Iflag,No,Icase,Iface,Type = ',5i4,2x,a4)
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
! determine # of dof for the face node
  call ndof_nod(face_type(Ntype,Iface),Norder(nre+Iface), &
                ndofH_face,ndofE_face,ndofV_face,ndofQ_face)
!
! check if a homogeneous Dirichlet node
  call homogenD(NORMAL,Icase,Bcond, is_homD,ncase,ibcnd)
  if (is_homD) then
    zuV = ZERO
    go to 100
  endif
!
! # of dof has to be positive...
  if (ndofV_face.lt.0) then
    write(*,*) 'dhpfaceV: ndofV_face = ',ndofV_face
    stop 1
  endif
!
! set order and orientation for all element edge nodes and the face node
  call initiate_order(Ntype, norder_1)
  norder_1(nre+Iface) = Norder(nre+Iface)
!
#if DEBUG_MODE
  if (iprint.eq.1) then
     write(*,7060) norder_1; call pause
7060 format('dhpfaceV: norder_1 = ',20i4)
  endif
#endif
!
! get face order to find out quadrature information
  call face_order(Ntype,Iface,Norder, norder_face)
  INTEGRATION=1   ! overintegrate
  call set_2Dint(face_type(Ntype,Iface),norder_face, &
                 nint,xi_list,wa_list)
  INTEGRATION=0   ! reset
!
! initiate stiffness matrix and load vector
  zbV = ZERO; aaV = 0.d0
!
! loop through integration points
  do l=1,nint
!
!   integration point local face coordinates
    t(1:2) = xi_list(1:2,l)
    wa     = wa_list(l)
!
!   get the corresponding master element coordinates and Jacobian
    call face_param(Ntype,Iface,t, xi,dxidt)
!
!   compute element H1 shape functions (for geometry)
    call shape3DH(Ntype,xi,norder_1,Nedge_orient,Nface_orient, &
                  nrdofH,shapH,gradH)
!
!   compute element H(div) shape functions
    call shape3DV(Ntype,xi,norder_1,Nface_orient, &
                  nrdofV,shapV,divV)

!   evaluate reference coordinates of the point as needed by GMP
    nsign = nsign_param(Ntype,Iface)
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
      write(*,*) 'dhpfaceV: Type = ', S_Type(Ntype)
      call logic_error(ERR_INVALID_VALUE,__FILE__,__LINE__)
    end select
!
!   compute inverse jacobian (for Piola transform)
    call geom(dxdeta, detadx, rjacdxdeta, nflag)
    if (nflag.ne.0) then
      write(*,*) 'dhpfaceV: rjacdxdeta = ',rjacdxdeta
      stop 1
    endif
!
!   get Dirichlet data at the point
    call dirichlet(Mdle,x,Icase, &
                   zvalH,zdvalH, zvalE,zdvalE, zvalV,zdvalV)
!
!   evaluate H(div) Dirichlet data in reference coordinates
    zvalVeta(1:3,1:MAXEQNV) = ZERO
    do i=1,3
      do j=1,3
        zvalVeta(i,1:MAXEQNV) = zvalVeta(i,1:MAXEQNV) &
                     + detadx(i,j)*zvalV(j,1:MAXEQNV)*rjacdxdeta
      enddo
    enddo
!
!   number of unused shape functions
    nrdofV = Iface-1
!
!   loop through element face test functions
    do j=1,ndofV_face
      kj = nrdofV + j
!
!     compute test function in reference coordinates
      vVeta(1:3) = (detadxi(1:3,1)*shapV(1,kj) &
                 +  detadxi(1:3,2)*shapV(2,kj) &
                 +  detadxi(1:3,3)*shapV(3,kj))/rjac
!
!     leave normal component only
      call dot_product(vVeta,rn, prod)
      vVeta(1:3) = prod*rn(1:3)
!
!     accumulate for the load vector
      zbV(j,1:MAXEQNV) = zbV(j,1:MAXEQNV) &
                       + (zvalVeta(1,1:MAXEQNV)*vVeta(1) &
                       +  zvalVeta(2,1:MAXEQNV)*vVeta(2) &
                       +  zvalVeta(3,1:MAXEQNV)*vVeta(3))*weight
!
!     loop through element face trial functions
      do i=1,ndofV_face
        ki = nrdofV + i
!
!       compute trial function in reference coordinates
        uVeta(1:3) = (detadxi(1:3,1)*shapV(1,ki) &
                   +  detadxi(1:3,2)*shapV(2,ki) &
                   +  detadxi(1:3,3)*shapV(3,ki))/rjac
!
!       no need to evaluate normal component...
!
!       accumulate for the stiffness matrix
        aaV(j,i) = aaV(j,i) &
                 + (vVeta(1)*uVeta(1) &
                 +  vVeta(2)*uVeta(2) &
                 +  vVeta(3)*uVeta(3))*weight
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
    write(*,*) 'dhpfaceV: LOAD VECTOR AND STIFFNESS MATRIX FOR ', &
               'ndofV_face = ',ndofV_face
    do j=1,ndofV_face
      write(*,7015) j, zbV(j,1:MAXEQNV)
      write(*,7016) aaV(j,1:ndofV_face)
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
  naV = MAXmdlqV
  call dgetrf(ndofV_face,ndofV_face,aaV,naV,ipivV,info)
  if (info.ne.0) then
    write(*,*)'dhpfaceV: H(div) DGETRF RETURNED INFO = ',info
    call logic_error(FAILURE,__FILE__,__LINE__)
  endif
!
! copy load vector
  zuV(1:ndofV_face,:) = zbV(1:ndofV_face,:)
!
#if C_MODE
! apply pivots to load vector
  call zlaswp(MAXEQNV,zuV,naV,1,ndofV_face,ipivV,1)
!
  uV_real(1:ndofV_face,:) =  real(zuV(1:ndofV_face,:))
  uV_imag(1:ndofV_face,:) = aimag(zuV(1:ndofV_face,:))
!
! triangular solves
  call dtrsm('L','L','N','U',ndofV_face,MAXEQNV,1.d0,aaV,naV, uV_real,naV)
  call dtrsm('L','U','N','N',ndofV_face,MAXEQNV,1.d0,aaV,naV, uV_real,naV)
!
  call dtrsm('L','L','N','U',ndofV_face,MAXEQNV,1.d0,aaV,naV, uV_imag,naV)
  call dtrsm('L','U','N','N',ndofV_face,MAXEQNV,1.d0,aaV,naV, uV_imag,naV)
!
  zuV(1:ndofV_face,:) = dcmplx(uV_real(1:ndofV_face,:), uV_imag(1:ndofV_face,:))
#else
! apply pivots to load vector
  call dlaswp(MAXEQNV,zuV,naV,1,ndofV_face,ipivV,1)
! triangular solves
  call dtrsm('L','L','N','U',ndofV_face,MAXEQNV,1.d0,aaV,naV, zuV,naV)
  call dtrsm('L','U','N','N',ndofV_face,MAXEQNV,1.d0,aaV,naV, zuV,naV)
#endif
!
#if DEBUG_MODE
  if (iprint.eq.1) then
   write(*,*) 'dhpfaceV: k,zu(k) = '
   do k=1,ndofV_face
     write(*,7015) k,zuV(k,1:MAXEQNV)
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
        write(*,*) 'dhpfaceV: ncase = ', ncase
      endif
#endif
!
!  ...initialize global variable counter, and node local variable counter
      ivarV=0; nvarV=0
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
            select case(D_TYPE(i))
!
!  .........H(div) component
            case(NORMAL)
!
!  ...........update global counter
              ivarV = ivarV + 1
!
!  ...........if the variable is supported by the node
              if (ncase(i).eq.1) then
!
!  .............update node local counter
                nvarV = nvarV + 1
!
!  .............do not write dof if physics attribute is deactivated
                if (.not. PHYSAm(i)) exit
!
!  .............store Dirichlet dof
                if (ibcnd(ic).eq.1) ZnodV(nvarV,1:ndofV_face) = zuV(1:ndofV_face,ivarV)
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
  end subroutine dhpfaceV

