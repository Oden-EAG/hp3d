!
#include "typedefs.h"
!
!-----------------------------------------------------------------------
!> @brief      determine middle node geometry dof interpolating GMP map
!!             using PB interpolation; NOTE: the interpolation
!!             (projection) is done in the reference space
!!
!> @param[in]  Iflag        - a flag specifying which of the objects the
!!                            face is on: 5 pris, 6 hexa, 7 tetr, 8 pyra
!> @param[in]  No           - number of a specific object
!> @param[in]  Etav         - reference coordinates of the element vertices
!> @param[in]  Ntype        - element (middle node) type
!> @param[in]  Nedge_orient - edge orientation
!> @param[in]  Nface_orient - face orientation
!> @param[in]  Norder       - element order
!> @param[in]  Xnod         - geometry dof for the element (vertex,edge
!!                            and face  values)
!> @param[out] Xdof         - geometry dof for the middle node
!!
!> @date       Apr 2024
!-----------------------------------------------------------------------
  subroutine hpmdle_opt(Mdle,Iflag,No,Etav,Ntype, &
                        Nedge_orient,Nface_orient,Norder, &
                        Xnod, Xdof)
  use parameters
  use physics
  use element_data
  use mpi_wrapper
  implicit none
!
! ** Arguments
!-----------------------------------------------------------------------
  integer,                         intent(in)  :: Iflag,No,Mdle
  real(8), dimension(3,8),         intent(in)  :: Etav
  integer,                         intent(in)  :: Ntype
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
  real(8), dimension(  MAXbrickH)       :: shapH
  real(8), dimension(3,MAXbrickH)       :: gradH
!
! derivatives of a shape function wrt reference coordinates
  real(8), dimension(3)                 :: dvHdeta
!
! geometry
  real(8)                               :: rjac
  real(8), dimension(3)                 :: xi,eta,x
  real(8), dimension(3,3)               :: detadxi,dxideta,dxdeta
!
! work space for linear solvers
  integer                                     :: info
  real(8), dimension(MAXmdlbH*(MAXmdlbH+1)/2) :: aaH_RFP
!
! load vector
  real(8), dimension(MAXmdlbH,3) :: bb
!
!..workspace for auxiliary matrix (storing info at integration points)
   real(8) :: A_TEST(MAXmdlbH,3*MAX_NINT3)
   integer :: lda,nda,noff,nRFP
!
! misc work space
  integer :: nrv,nre,nrf,i,j,k,kj,iflag1
  integer :: ndofH_mdle,ndofE_mdle,ndofV_mdle,ndofQ_Mdle
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
     write(*,7010) Mdle,Iflag,No,S_Type(Ntype)
7010 format('hpmdle: Mdle,Iflag,No,Type = ',3i4,2x,a4)
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
! determine # of dof for the mdle node
  call ndof_nod(Ntype,Norder(nre+nrf+1), &
                ndofH_mdle,ndofE_mdle,ndofV_mdle,ndofQ_mdle)
!
! if # of dof is zero, return, nothing to do
  if (ndofH_mdle.eq.0) return
!
! get quadrature
  call set_3Dint(Ntype,Norder, nint,xi_list,wa_list)
!
! initiate load vector
  bb = 0.d0
!
! loop through integration points
  do l=1,nint
!
!   integration point local face coordinates
    xi(1:3) = xi_list(1:3,l)
    wa      = wa_list(l)
!
!   compute element H1 shape functions
    call shape3DH(Ntype,xi,Norder,Nedge_orient,Nface_orient, &
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
      write(*,*) 'hpmdle: Type = ', S_Type(Ntype)
      call logic_error(ERR_INVALID_VALUE,__FILE__,__LINE__)
    end select
!
!   remove the contributions from vertices, edges and faces
!
!   loop through vertex, edge, and face shape functions
    nrdofH = nrdofH - ndofH_mdle
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
!   loop through middle node shape functions
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
!     fill auxiliary matrix for stiffness
      A_TEST(j,noff+1:noff+3) = dvHdeta(1:3) * sqrt(weight)
!
!   end of loop through shape functions
    enddo
!
! end of loop through integration points
  enddo
!
!  compute stiffness matrix in RFP format
   nda = 3*nint; lda = MAXmdlbH
   nRFP = ndofH_mdle*(ndofH_mdle+1)/2
   call DSFRK('N','U','N',ndofH_mdle,nda,       &
               1.0d0,A_TEST(1:lda,1:nda),lda,   &
               0.0d0,aaH_RFP(1:nRFP))
!
#if HP3D_DEBUG
  if (iprint.eq.1) then
    write(*,*) 'hpmdle: LOAD VECTOR AND STIFFNESS MATRIX FOR ', &
               'ndofH_mdle = ',ndofH_mdle
    do j=1,ndofH_mdle
      write(*,7015) j, bb(j,1:3)
!      write(*,7016) aaH(j,1:ndofH_mdle)
    enddo
7015    format(i5, 10e12.5)
!7016    format(10e12.5)
  endif
#endif
!
!-----------------------------------------------------------------------
!
! compute Cholesky decomposition
  call DPFTRF('N','U',ndofH_mdle,aaH_RFP(1:nRFP),info)
  if (info.ne.0) then
    write(*,*)'hpmdle:  DPFTRF RETURNED INFO = ',info
    call logic_error(FAILURE,__FILE__,__LINE__)
  endif
!
! solve linear system
  call DPFTRS('N','U',ndofH_mdle,3,aaH_RFP(1:nRFP),bb(1:lda,1:3),lda, info)
  if (info.ne.0) then
    write(*,*)'hpmdle:  DPFTRS RETURNED INFO = ',info
    call logic_error(FAILURE,__FILE__,__LINE__)
  endif
!
#if HP3D_DEBUG
  if (iprint.eq.1) then
   write(*,*) 'hpmdle: k,bb(k) = '
   do k=1,ndofH_mdle
     write(*,7015) k,bb(k,1:3)
   enddo
   call pause
  endif
#endif
!
! save the dof
  do i=1,3
    Xdof(i,1:ndofH_mdle) = bb(1:ndofH_mdle,i)
  enddo
!
!..TIMER
!   end_time = MPI_Wtime()
!   !$OMP CRITICAL
!   write(*,11) 'hpmdle_opt: ', end_time-start_time
!11 format(A,f12.5,' s')
!   !$OMP END CRITICAL
!
  end subroutine hpmdle_opt
