!
#include "typedefs.h"
!
!-----------------------------------------------------------------------
!> @brief    determine H(div) face dof interpolating H(div) Dirichlet data
!!           using PB interpolation; NOTE: the interpolation (projection)
!!           is done in the reference space
!!
!> @param[in]  Mdle         - element (middle node) number
!> @param[in]  Iflag        - a flag specifying which of the objects the
!!                            face is on: 5 pris, 6 hexa, 7 tetr, 8 pyra
!> @param[in]  No           - number of a specific object
!> @param[in]  Etav         - reference coordinates of the element vertices
!> @param[in]  Ntype        - element (middle node) type
!> @param[in]  Icase        - the face node case
!> @param[in]  Bcond        - the edge node BC flag
!> @param[in]  Nedge_orient - edge orientation, never used
!> @param[in]  Nface_orient - face orientation
!> @param[in]  Norder       - element order
!> @param[in]  Iface        - face number
!!
!> @param[in,out] ZnodV     - H(div) dof for the face
!!
!> @date Sep 2023
!-----------------------------------------------------------------------
subroutine dhpfaceV_opt(Mdle,Iflag,No,Etav,Ntype,Icase,Bcond,   &
                        Nedge_orient,Nface_orient,Norder,Iface, &
                        ZnodV)
  use control
  use parameters
  use physics
  use element_data
  use mpi_wrapper
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
  VTYPE,   dimension(NRRHS*NREQNH(Icase),*), intent(inout) :: ZnodV
!
! ** Locals
!-----------------------------------------------------------------------
!
! face nodes order (to set the quadrature)
  integer, dimension(5)                 :: norder_face
!
! quadrature
  integer                               :: l,nint
  real(8), dimension(2, MAX_NINT2)      :: xi_list
  real(8), dimension(   MAX_NINT2)      :: wa_list
  real(8)                               :: wa, weight
!
! work space for shape3DH
  integer                               :: nrdofH
  integer, dimension(19)                :: norder_ifc,norder_1
!
  real(8), dimension(  MAXbrickH)       :: shapH
  real(8), dimension(3,MAXbrickH)       :: gradH
!
! work space for shape3DV
  integer                               :: nrdofV
  real(8), dimension(3,MAXbrickV)       :: shapV
  real(8), dimension(  MAXbrickV)       :: divV
!
! H(div) test and trial shape function in reference coordinates
  real(8), dimension(3)                 :: vVeta
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
  integer                                     :: info
  real(8), dimension(MAXmdlqV*(MAXmdlqV+1)/2) :: aaV_RFP
!
! load vector and solution
  VTYPE,   dimension(MAXmdlqV,MAXEQNV)  :: zbV
#if HP3D_COMPLEX
  real(8), dimension(MAXmdlqV,MAXEQNV)  :: uV_real,uV_imag
#endif
!
!..workspace for auxiliary matrix (storing info at integration points)
   real(8) :: A_TEST(MAXmdlqV,3*MAX_NINT2)
   integer :: lda,nda,noff,nRFP
!
! decoded case and BC flag for the face node
  integer, dimension(NR_PHYSA)          :: ncase
  integer, dimension(NRINDEX_HEV)       :: ibcnd
!
! misc work space
  integer :: nrv,nre,nrf,nsign,nflag,i,j,k,ivarV,nvarV,kj,ic
  integer :: ndofH_face,ndofE_face,ndofV_face,ndofQ_face
!
!..TIMER
!   real(8) :: start_time,end_time
!
  logical :: is_homD
!
#if HP3D_DEBUG
  integer :: iprint
  iprint=0
#endif
!
!-----------------------------------------------------------------------
!..TIMER
!   start_time = MPI_Wtime()
!
  nrv = nvert(Ntype); nre = nedge(Ntype); nrf = nface(Ntype)
!
#if HP3D_DEBUG
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
    zbV = ZERO
    goto 100
  endif
!
! if #dof is zero, return (nothing to do)
  if (ndofV_face.eq.0) return
!
! set order and orientation for all element edge nodes and the face node
  call initiate_order(Ntype, norder_1)
  norder_ifc = norder_1
  norder_ifc(nre+Iface) = Norder(nre+Iface)
!
#if HP3D_DEBUG
  if (iprint.eq.1) then
     write(*,7060) norder_ifc; call pause
7060 format('dhpfaceV: norder_ifc = ',20i4)
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
  zbV = ZERO
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
    call shape3DV(Ntype,xi,norder_ifc,Nface_orient, &
                  nrdofV,shapV,divV)

!   evaluate reference coordinates of the point as needed by GMP
!   brefgeom3D returns the outward normal as seen from the local element
    nsign = nsign_param(Ntype,Iface)
    call brefgeom3D(Mdle,xi,Etav,shapH(1:8),gradH(1:3,1:8),nrv,dxidt, &
                    nsign, eta,detadxi,dxideta,rjac,detadt,rn,bjac)
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
!   compute inverse Jacobian (for Piola transform)
    call geom(dxdeta, detadx,rjacdxdeta,nflag)
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
!   offset in auxiliary matrix
    noff = 3*(l-1)
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
!     fill auxiliary matrix for stiffness
      A_TEST(j,noff+1:noff+3) = vVeta(1:3) * sqrt(weight)
!
!   end of loop through face test functions
    enddo
!
! end of loop through integration points
  enddo
!
!  compute stiffness matrix in RFP format
   nda = 3*nint; lda = MAXmdlqV
   nRFP = ndofV_face*(ndofV_face+1)/2
   call DSFRK('N','U','N',ndofV_face,nda,       &
               1.0d0,A_TEST(1:lda,1:nda),lda,   &
               0.0d0,aaV_RFP(1:nRFP))
!
#if HP3D_DEBUG
  if (iprint.eq.1) then
    write(*,*) 'dhpfaceV: LOAD VECTOR AND STIFFNESS MATRIX FOR ', &
               'ndofV_face = ',ndofV_face
    do j=1,ndofV_face
      write(*,7015) j, zbV(j,1:MAXEQNV)
      !write(*,7016) aaV(j,1:ndofV_face)
    enddo
  endif
#if HP3D_COMPLEX
  7015 format(i5,2x,6(2e10.3,2x))
#else
  7015 format(i5,2x,10e12.5)
#endif
  !7016 format(10e12.5)
#endif
!
!-----------------------------------------------------------------------
!
! compute Cholesky decomposition
  call DPFTRF('N','U',ndofV_face,aaV_RFP(1:nRFP),info)
  if (info.ne.0) then
    write(*,*)'dhpfaceV: DPFTRF RETURNED INFO = ',info
    call logic_error(FAILURE,__FILE__,__LINE__)
  endif
!
#if HP3D_COMPLEX
!
  uV_real(1:ndofV_face,:) =  real(zbV(1:ndofV_face,:))
  uV_imag(1:ndofV_face,:) = aimag(zbV(1:ndofV_face,:))
!
! solve linear system
  call DPFTRS('N','U',ndofV_face,MAXEQNV,aaV_RFP(1:nRFP), &
              uV_real(1:lda,1:MAXEQNV),lda, info)
  call DPFTRS('N','U',ndofV_face,MAXEQNV,aaV_RFP(1:nRFP), &
              uV_imag(1:lda,1:MAXEQNV),lda, info)
  if (info.ne.0) then
    write(*,*)'dhpfaceV: DPFTRS RETURNED INFO = ',info
    call logic_error(FAILURE,__FILE__,__LINE__)
  endif
!
  zbV(1:ndofV_face,:) = dcmplx(uV_real(1:ndofV_face,:), uV_imag(1:ndofV_face,:))
#else
!
! solve linear system
  call DPFTRS('N','U',ndofV_face,MAXEQNV,aaV_RFP(1:nRFP), &
              zbV(1:lda,1:MAXEQNV),lda, info)
  if (info.ne.0) then
    write(*,*)'dhpfaceV: DPFTRS RETURNED INFO = ',info
    call logic_error(FAILURE,__FILE__,__LINE__)
  endif
#endif
!
#if HP3D_DEBUG
  if (iprint.eq.1) then
   write(*,*) 'dhpfaceV: k,zu(k) = '
   do k=1,ndofV_face
     write(*,7015) k,zbV(k,1:MAXEQNV)
   enddo
   call pause
  endif
#endif
!
!  ...save the DOFs, skipping irrelevant entries
  100 continue
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,*) 'dhpfaceV: ncase = ', ncase
      endif
#endif
!
!  ...initialize global variable counter, and node local variable counter
      ivarV=0; nvarV=0
!
!  ...loop through multiple loads
      do j=1,NRRHS
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
                if (ibcnd(ic).eq.1) ZnodV(nvarV,1:ndofV_face) = zbV(1:ndofV_face,ivarV)
!
              endif
            end select
!
!  .......loop through components
          enddo
!  .....loop through physical attributes
        enddo
!  ...loop through multiple loads
      enddo
!
!..TIMER
!   end_time = MPI_Wtime()
!   !$OMP CRITICAL
!   write(*,11) 'dhpfaceV_opt: ', end_time-start_time
!11 format(A,f12.5,' s')
!   !$OMP END CRITICAL
!
#if HP3D_DEBUG
      if (iprint.eq.1) call result
#endif
!
end subroutine dhpfaceV_opt
