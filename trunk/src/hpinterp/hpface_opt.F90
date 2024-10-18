!
#include "typedefs.h"
!
!-----------------------------------------------------------------------
!> @brief      determine face geometry dof interpolating GMP map using
!!             PB interpolation; NOTE: the interpolation (projection)
!!             is done in the reference space
!!
!> @param[in]  Iflag        - a flag specifying which of the objects the
!!                            face is on: 5 pris, 6 hexa, 7 tetr, 8 pyra
!> @param[in]  No           - number of a specific object
!> @param[in]  Etav         - reference coordinates of the element vertices
!> @param[in]  Ntype        - element (middle node) type
!> @param[in]  Nedge_orient - edge orientation
!> @param[in]  Nface_orient - face orientation
!> @param[in]  Norder       - element order
!> @param[in]  Iface        - face number
!> @param[in]  Xnod         - geometry dof for the element (vertex and edge values)
!!
!> @param[out] Xdof         - geometry dof for the face
!!
!> @date       Apr 2024
!-----------------------------------------------------------------------
  subroutine hpface_opt(Mdle,Iflag,No,Etav,Ntype, &
                        Nedge_orient,Nface_orient,Norder,Iface, &
                        Xnod, Xdof)
  use parameters
  use element_data
  use mpi_wrapper
  implicit none
!
! ** Arguments
!-----------------------------------------------------------------------
  integer,                          intent(in)  :: Iflag,No,Mdle
  integer,                          intent(in)  :: Iface
  real(8), dimension(3,8),          intent(in)  :: Etav
  integer,                          intent(in)  :: Ntype
  integer, dimension(12),           intent(in)  :: Nedge_orient
  integer, dimension(6),            intent(in)  :: Nface_orient
  integer, dimension(19),           intent(in)  :: Norder
!
  real(8), dimension(3,MAXbrickH),  intent(in)  :: Xnod
  real(8), dimension(3,*),          intent(out) :: Xdof
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
! work space for shape3H
  integer                               :: nrdofH
  integer, dimension(19)                :: norder_ifc
  real(8), dimension(  MAXbrickH)       :: shapH
  real(8), dimension(3,MAXbrickH)       :: gradH
!
! derivatives of a shape function wrt reference coordinates
  real(8), dimension(3)                 :: dvHdeta
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
! work space for linear solvers
  integer                                     :: info
  real(8), dimension(MAXmdlqH*(MAXmdlqH+1)/2) :: aaH_RFP
!
! load vector
  real(8), dimension(MAXmdlqH,3) :: bb
!
!..workspace for auxiliary matrix (storing info at integration points)
   real(8) :: A_TEST(MAXmdlqH,3*MAX_NINT2)
   integer :: lda,nda,noff,nRFP
!
! misc work space
  integer :: nrv,nre,nrf,i,j,k,ie,kj
  integer :: ndofH_face,ndofE_face,ndofV_face,ndofQ_Face,nsign
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
     write(*,7010) Mdle,Iflag,No,Iface,S_Type(Ntype)
7010 format('hpface: Mdle,Iflag,No,Iface,Type = ',4i4,2x,a4)
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
! determine # of dof for the face node
  call ndof_nod(face_type(Ntype,Iface),Norder(nre+Iface), &
                ndofH_face,ndofE_face,ndofV_face,ndofQ_face)
!
! if # of dof is zero, return, nothing to do
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
7060 format('hpface: norder_ifc = ',20i4)
  endif
#endif
!
! get face order to find out quadrature information
  call face_order(Ntype,Iface,Norder, norder_face)
  call set_2Dint(face_type(Ntype,Iface),norder_face, &
                 nint,xi_list,wa_list)
!
! initiate load vector
  bb = 0.d0
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
      write(*,*) 'hpface: Type = ', S_Type(Ntype)
      call logic_error(ERR_INVALID_VALUE,__FILE__,__LINE__)
    end select
!
!   remove the contributions from vertices and edges...
!
!   loop through vertex and edge shape functions
    nrdofH = nrdofH - ndofH_face
    do k=1,nrdofH
!
!     compute derivatives of the shape function wrt reference coordinates
      dvHdeta(1:3) = gradH(1,k)*dxideta(1,1:3) &
                   + gradH(2,k)*dxideta(2,1:3) &
                   + gradH(3,k)*dxideta(3,1:3)
!
!     subtract...
      do i=1,3
        dxdeta(1:3,i) = dxdeta(1:3,i) &
                      - Xnod(1:3,k)*dvHdeta(i)
      enddo
    enddo
!
!   offset in auxiliary matrix
    noff = 3*(l-1)
!
!   loop through face shape functions
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
      bb(j,1:3) = bb(j,1:3) &
                + (dxdeta(1:3,1)*dvHdeta(1) &
                +  dxdeta(1:3,2)*dvHdeta(2) &
                +  dxdeta(1:3,3)*dvHdeta(3))*weight
!
!     fill auxiliary matrix for stiffness
      A_TEST(j,noff+1:noff+3) = dvHdeta(1:3) * sqrt(weight)
!
!   end of loop through face shape functions
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
    write(*,*) 'hpface: LOAD VECTOR AND STIFFNESS MATRIX FOR ', &
               'ndofH_face = ',ndofH_face
    do j=1,ndofH_face
      write(*,7015) j, bb(j,1:3)
!      write(*,7016) aaH(j,1:ndofH_face)
    enddo
7015    format(i5, 10e12.5)
!7016    format(10e12.5)
  endif
#endif
!
!-----------------------------------------------------------------------
!
! compute Cholesky decomposition
  call DPFTRF('N','U',ndofH_face,aaH_RFP(1:nRFP),info)
  if (info.ne.0) then
    write(*,*)'hpface:  DPFTRF RETURNED INFO = ',info
    call logic_error(FAILURE,__FILE__,__LINE__)
  endif

! solve linear system
  call DPFTRS('N','U',ndofH_face,3,aaH_RFP(1:nRFP),bb(1:lda,1:3),lda, info)
  if (info.ne.0) then
    write(*,*)'hpface:  DPFTRS RETURNED INFO = ',info
    call logic_error(FAILURE,__FILE__,__LINE__)
  endif
!
#if HP3D_DEBUG
  if (iprint.eq.1) then
   write(*,*) 'hpface: k,bb(k) = '
   do k=1,ndofH_face
     write(*,7015) k,bb(k,1:3)
   enddo
   call pause
  endif
#endif
!
! save the dof
  do i=1,3
    Xdof(i,1:ndofH_face) = bb(1:ndofH_face,i)
  enddo
!
!..TIMER
!   end_time = MPI_Wtime()
!   !$OMP CRITICAL
!   write(*,11) 'hpface_opt: ', end_time-start_time
!11 format(A,f12.5,' s')
!   !$OMP END CRITICAL
!
  end subroutine hpface_opt
