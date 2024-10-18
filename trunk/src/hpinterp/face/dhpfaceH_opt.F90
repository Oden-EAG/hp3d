!
#include "typedefs.h"
!
!-----------------------------------------------------------------------
!> @brief    determine H1 face dof interpolating H1 Dirichlet data using
!!           PB interpolation; NOTE: the interpolation (projection) is
!!           done in the reference space
!!
!> @param[in]  Mdle         - element (middle node) number
!> @param[in]  Iflag        - a flag specifying which of the objects the
!!                            face is on: 5 pris, 6 hexa, 7 tetr, 8 pyra
!> @param[in]  No           - number of a specific object
!> @param[in]  Etav         - reference coordinates of the element vertices
!> @param[in]  Ntype        - element (middle node) type
!> @param[in]  Icase        - the face node case
!> @param[in]  Bcond        - the face node BC flag
!> @param[in]  Nedge_orient - edge orientation
!> @param[in]  Nface_orient - face orientation
!> @param[in]  Norder       - element order
!> @param[in]  Iface        - face number
!> @param[in]  ZdofH        - H1 dof for the element (vertex and edge values)
!!
!> @param[in,out] ZnodH     - H1 dof for the face
!!
!> @date Sep 2023
!-----------------------------------------------------------------------
subroutine dhpfaceH_opt(Mdle,Iflag,No,Etav,Ntype,Icase,Bcond,   &
                        Nedge_orient,Nface_orient,Norder,Iface, &
                        ZdofH, ZnodH)
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
  VTYPE,   dimension(MAXEQNH,MAXbrickH),     intent(in)    :: ZdofH
  VTYPE,   dimension(NRRHS*NREQNH(Icase),*), intent(inout) :: ZnodH
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
  integer, dimension(19)                :: norder_ifc
  real(8), dimension(  MAXbrickH)       :: shapH
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
  integer                                     :: info
  real(8), dimension(MAXmdlqH*(MAXmdlqH+1)/2) :: aaH_RFP
!
! load vector and solution
  VTYPE,   dimension(MAXmdlqH,MAXEQNH)  :: zbH
#if HP3D_COMPLEX
  real(8), dimension(MAXmdlqH,MAXEQNH)  :: duH_real,duH_imag
#endif
!
!..workspace for auxiliary matrix (storing info at integration points)
   real(8) :: A_TEST(MAXmdlqH,3*MAX_NINT2)
   integer :: lda,nda,noff,nRFP
!
! decoded case and BC flags for the edge node
  integer, dimension(NR_PHYSA)          :: ncase
  integer, dimension(NRINDEX_HEV)       :: ibcnd
!
! misc work space
  integer :: nrv,nre,nrf,i,j,k,ie,ivarH,nvarH,kj,ic
  integer :: ndofH_face,ndofE_face,ndofV_face,ndofQ_Face,nsign
!
  logical :: is_homD
!
!..TIMER
!   real(8) :: start_time,end_time
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
! determine # of dof for the face node
  call ndof_nod(face_type(Ntype,Iface),Norder(nre+Iface), &
                ndofH_face,ndofE_face,ndofV_face,ndofQ_face)
!
! check if a homogeneous Dirichlet node
  call homogenD(CONTIN,Icase,Bcond, is_homD,ncase,ibcnd)
  if (is_homD) then
    zbH = ZERO
    goto 100
  endif
!
! if #dof is zero, return (nothing to do)
  if (ndofH_face.eq.0) return
!
! set order and orientation for all element edge nodes and the face node
  call initiate_order(Ntype, norder_ifc)
  do ie=1,nre
    norder_ifc(ie) = Norder(ie)
  enddo
  norder_ifc(nre+Iface) = Norder(nre+Iface)
!
#if HP3D_DEBUG
  if (iprint.eq.1) then
     write(*,7060) norder_ifc; call pause
7060 format('dhpfaceH: norder_ifc = ',20i4)
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
  zbH = ZERO
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
!   compute element H1 shape functions
    call shape3DH(Ntype,xi,norder_ifc,Nedge_orient,Nface_orient, &
                  nrdofH,shapH,gradH)
!
!   evaluate reference coordinates of the point as needed by GMP
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
      write(*,*) 'dhpfaceH: Type = ', S_Type(Ntype)
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
!   offset in auxiliary matrix
    noff = 3*(l-1)
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
!     fill auxiliary matrix for stiffness
      A_TEST(j,noff+1:noff+3) = dvHdeta(1:3) * sqrt(weight)
!
!   end of loop through face test functions
    enddo
!
! end of loop through integration points
  enddo
!
!  compute stiffness matrix in RFP format
   nda = 3*nint; lda = MAXmdlqH
   nRFP = ndofH_face*(ndofH_face+1)/2
   call DSFRK('N','U','N',ndofH_face,nda,       &
               1.0d0,A_TEST(1:lda,1:nda),lda,   &
               0.0d0,aaH_RFP(1:nRFP))
!
#if HP3D_DEBUG
  if (iprint.eq.1) then
    write(*,*) 'dhpfaceH: LOAD VECTOR AND STIFFNESS MATRIX FOR ', &
               'ndofH_face = ',ndofH_face
    do j=1,ndofH_face
      write(*,7015) j, zbH(j,1:MAXEQNH)
      !write(*,7016) aaH(j,1:ndofH_face)
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
  call DPFTRF('N','U',ndofH_face,aaH_RFP(1:nRFP),info)
  if (info.ne.0) then
    write(*,*)'dhpfaceH: DPFTRF RETURNED INFO = ',info
    call logic_error(FAILURE,__FILE__,__LINE__)
  endif
!
#if HP3D_COMPLEX
!
  duH_real(1:ndofH_face,:) =  real(zbH(1:ndofH_face,:))
  duH_imag(1:ndofH_face,:) = aimag(zbH(1:ndofH_face,:))
!
! solve linear system
  call DPFTRS('N','U',ndofH_face,MAXEQNH,aaH_RFP(1:nRFP), &
              duH_real(1:lda,1:MAXEQNH),lda, info)
  call DPFTRS('N','U',ndofH_face,MAXEQNH,aaH_RFP(1:nRFP), &
              duH_imag(1:lda,1:MAXEQNH),lda, info)
  if (info.ne.0) then
    write(*,*)'dhpfaceH: DPFTRS RETURNED INFO = ',info
    call logic_error(FAILURE,__FILE__,__LINE__)
  endif
!
  zbH(1:ndofH_face,:) = dcmplx(duH_real(1:ndofH_face,:), duH_imag(1:ndofH_face,:))
#else
!
! solve linear system
  call DPFTRS('N','U',ndofH_face,MAXEQNH,aaH_RFP(1:nRFP), &
              zbH(1:lda,1:MAXEQNH),lda, info)
  if (info.ne.0) then
    write(*,*)'dhpfaceH: DPFTRS RETURNED INFO = ',info
    call logic_error(FAILURE,__FILE__,__LINE__)
  endif
#endif
!
#if HP3D_DEBUG
  if (iprint.eq.1) then
   write(*,*) 'dhpfaceH: k,zu(k) = '
   do k=1,ndofH_face
     write(*,7015) k,zbH(k,1:MAXEQNH)
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
        write(*,*) 'dhpfaceH: ncase = ', ncase
      endif
#endif
!
!  ...initialize global variable counter, and node local variable counter
      ivarH=0 ; nvarH=0
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
                if (ibcnd(ic).eq.1) ZnodH(nvarH,1:ndofH_face) = zbH(1:ndofH_face,ivarH)
!
              endif
!
            end select
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
!   write(*,11) 'dhpfaceH_opt: ', end_time-start_time
!11 format(A,f12.5,' s')
!   !$OMP END CRITICAL
!
#if HP3D_DEBUG
      if (iprint.eq.1) call result
#endif
!
end subroutine dhpfaceH_opt
