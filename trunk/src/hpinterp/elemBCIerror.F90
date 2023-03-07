!
#include "typedefs.h"
!
!-----------------------------------------------------------------------
!> @brief      compute PB interpolation error for Dirichlet data
!!
!> @param[in]  Mdle   - element (middle node) number
!> @param[out] ErrorH - element interpolation error for H1      BC data
!> @param[out] ErrorE - element interpolation error for H(curl) BC data
!> @param[out] ErrorV - element interpolation error for H(div)  BC data
!!
!> @date       Feb 2023
!-----------------------------------------------------------------------
subroutine elemBCIerror(Mdle, ErrorH,ErrorE,ErrorV)
!
  use control
  use data_structure3D
  use element_data
  implicit none
!
! ** Arguments
!-----------------------------------------------------------------------
!
  integer, intent(in)   :: Mdle
  real(8), intent(out)  :: ErrorH,ErrorE,ErrorV
!
! ** Locals
!-----------------------------------------------------------------------
!
! element type
  integer                               :: ntype
!
! GMP block type and number
  integer                               :: iflag,no
!
! element vertices reference coordinates
  real(8), dimension(3,8)               :: etav
!
! element order
  integer, dimension(19)                :: norder
!
! element orientation for edges
  integer, dimension(12)                :: nedge_orient
!
! element orientation for faces
  integer, dimension(6)                 :: nface_orient
!
! work space for elem nodes
  integer, dimension(27)                :: nodesl,norientl
!
! index for a face node
  integer, dimension(NRINDEX)           :: index
!
! is there any Dirichlet data at the node:
  logical                               :: Dirichlet_node
!
! element solution dof
  VTYPE, dimension(MAXEQNH,MAXbrickH)   :: zdofH
  VTYPE, dimension(MAXEQNE,MAXbrickE)   :: zdofE
  VTYPE, dimension(MAXEQNV,MAXbrickV)   :: zdofV
  VTYPE, dimension(MAXEQNQ,MAXbrickQ)   :: zdofQ
!
! face nodes order (to set the quadrature)
  integer, dimension(5)                 :: norder_face
!
! quadrature
  integer                               :: l,nint
  real(8), dimension(2, MAX_NINT2)      :: xi_list
  real(8), dimension(   MAX_NINT2)      :: wa_list
  real(8)                               :: wa,weight
!
! work space for shape3H
  integer                               :: nrdofH
  real(8), dimension(MAXbrickH)         :: shapH
  real(8), dimension(3,MAXbrickH)       :: gradH
!
! derivatives of a H1 shape function wrt reference coordinates
  real(8), dimension(3)                 :: duHdeta
!
! work space for shape3E
  integer                               :: nrdofE
  real(8), dimension(3,MAXbrickE)       :: shapE
  real(8), dimension(3,MAXbrickE)       :: curlE
!
! H(curl) shape functions in reference coordinates
  real(8), dimension(3)                 :: uEeta,curluEeta
!
! work space for shape3V
  integer                               :: nrdofV
  real(8), dimension(3,MAXbrickH)       :: shapV
  real(8), dimension(MAXbrickH)         :: divV
!
! H(div) test and trial shape function in reference coordinates
  real(8), dimension(3)                 :: uVeta
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
           zvalE(3,MAXEQNE), zdvalE(3,MAXEQNE,3), zcurlE(3,MAXEQNE), &
           zvalV(3,MAXEQNV), zdvalV(3,MAXEQNV,3)
!
! Dirichlet data in reference coordinates
  VTYPE :: zdvalHdeta(MAXEQNH,3)
  VTYPE :: zvalEeta(3,MAXEQNE),zcvalEeta(3,MAXEQNE)
  VTYPE :: zvalVeta(3,MAXEQNV)
!
! approximate solution in reference coordinates
  VTYPE :: zdsolHdeta(MAXEQNH,3)
  VTYPE :: zsolEeta(3,MAXEQNE),zcsolEeta(3,MAXEQNE)
  VTYPE :: zsolVeta(3,MAXEQNV)
!
! difference between the exact and approximate solutions
  VTYPE :: zddifHdeta(MAXEQNH,3)
  VTYPE :: zdifEeta(3,MAXEQNE),zcdifEeta(3,MAXEQNE)
  VTYPE :: zdifVeta(3,MAXEQNV)
!
! error per component
  real(8) :: derrorH(MAXEQNH),derrorE(MAXEQNE),derrorV(MAXEQNV)
!
! misc
  integer :: nrv,nre,nrf,ifc,nod,nflag,kH,kE,kV,i,j,nsign
  integer :: ivarH,ivarE,ivarV
!
#if DEBUG_MODE
  integer :: iprint
#endif
!
!-----------------------------------------------------------------------
!
#if DEBUG_MODE
  select case(Mdle)
    case(355)   ; iprint=0
    case default; iprint=0
  end select
#endif
!
  ntype = NODES(Mdle)%ntype
  nrv = nvert(ntype); nre = nedge(ntype); nrf = nface(ntype)
  call refel(Mdle, iflag,no,etav)
  call elem_nodes(Mdle, nodesl, norientl)
  call find_orient(mdle, nedge_orient,nface_orient)
  call find_order(mdle, norder)
!
#if DEBUG_MODE
  if (iprint.eq.1) then
     write(*,7010) Mdle,S_Type(ntype)
7010 format('elemBCIerror: Mdle,type = ',i4,2x,a4)
     write(*,7020) etav(1:3,1:nrv)
7020 format('          etav = ',8(3f6.2,1x))
     write(*,7030) nedge_orient(1:nre)
7030 format('          nedge_orient = ',12i2)
     write(*,7040) nface_orient(1:nrf)
7040 format('          nface_orient = ',6i2)
     write(*,7050) norder(1:nre+nrf+1)
7050 format('          norder = ',19i4)
     call pause
  endif
#endif
!
! get solution dof
  call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
!
! initialize errors
  derrorH = 0.d0; derrorE = 0.d0; derrorV = 0.d0
  ivarH=0; ivarE=0; ivarV=0
!
! loop through element faces
  do ifc=1,nrf
!
!   face node
    nod = nodesl(nrv+nre+ifc)
!
!   determine if there are any variables are flagged with Dirichlet flags
    Dirichlet_node = .false.
    call get_index(nod, index)
    do i=1,NRINDEX
      select case(index(i))
      case(1,3,5); Dirichlet_node = .true.
      end select
    enddo
    if (.not.Dirichlet_node) cycle
!
!   get face order to find out quadrature information
    call face_order(ntype,ifc,norder, norder_face)
    INTEGRATION=1   ! overintegrate
    call set_2Dint(face_type(ntype,ifc),norder_face, &
                   nint,xi_list,wa_list)
    INTEGRATION=0   ! reset
!
!   loop through integration points
    do l=1,nint
      t(1:2) = xi_list(1:2,l)
      wa     = wa_list(l)
!
!     get the corresponding master element coordinates and Jacobian
      call face_param(ntype,ifc,t, xi,dxidt)
!
!     compute element H1 shape functions
      call shape3DH(ntype,xi,norder,nedge_orient,nface_orient, &
                    nrdofH,shapH,gradH)
!
!     compute element Hcurl shape functions
      call shape3DE(ntype,xi,norder,nedge_orient,nface_orient, &
                    nrdofE,shapE,curlE)
!
!     compute element H(div) shape functions
      call shape3DV(ntype,xi,norder,nface_orient, &
                    nrdofV,shapV,divV)
!
!     evaluate reference coordinates of the point as needed by GMP
      nsign = nsign_param(ntype,ifc)
      call brefgeom3D(Mdle,xi,etav,shapH,gradH,nrv,dxidt,nsign, &
                      eta,detadxi,dxideta,rjac,detadt,rn,bjac)
      weight = wa*bjac
!
#if DEBUG_MODE
      if (iprint.eq.1) then
        write(*,7100) xi(1:2),eta(1:3),detadxi(1:3,1:3),rn(1:3),bjac
7100    format('elemBCIerror: xi,eta  = ',2f8.3,3x,3f8.3,/, &
               '              detadxi = ',3(3f8.3,2x),/, &
               '              rn,bjac = ',3f8.3,3x,f8.3)
      endif
#endif
!
!     call GMP routines to evaluate physical coordinates and their
!     derivatives wrt reference coordinates
      select case(iflag)
      case(5);        call prism(no,eta, x,dxdeta)
      case(6);        call  hexa(no,eta, x,dxdeta)
      case(7);        call tetra(no,eta, x,dxdeta)
      case(8);        call pyram(no,eta, x,dxdeta)
      case default
        write(*,*) 'dhpfaceH: Type = ', S_Type(ntype)
        call logic_error(ERR_INVALID_VALUE,__FILE__,__LINE__)
      end select
!
!     compute inverse jacobian (for transforming the curl)
      call geom(dxdeta, detadx, rjacdxdeta, nflag)
      if (nflag.ne.0) then
        write(*,*) 'elemBCIerror: rjacdxdeta = ',rjacdxdeta
        stop 1
      endif
!
!     get Dirichlet data at the point
      call dirichlet(Mdle,x,NODES(Mdle)%case, &
                     zvalH,zdvalH, zvalE,zdvalE, zvalV,zdvalV)
!
!     compute curl of the Dirichlet data
      zcurlE(1,1:MAXEQNE) = zdvalE(3,1:MAXEQNE,2) - zdvalE(2,1:MAXEQNE,3)
      zcurlE(2,1:MAXEQNE) = zdvalE(1,1:MAXEQNE,3) - zdvalE(3,1:MAXEQNE,1)
      zcurlE(3,1:MAXEQNE) = zdvalE(2,1:MAXEQNE,1) - zdvalE(1,1:MAXEQNE,2)
!
!     use Piola transforms to transform the Dirichlet data to reference
!     coordinates
      zdvalHdeta(1:MAXEQNH,1:3) = ZERO
      zvalEeta (1:3,1:MAXEQNE) = ZERO
      zcvalEeta(1:3,1:MAXEQNE) = ZERO
      zvalVeta(1:3,1:MAXEQNV) = ZERO
      do i=1,3
        do j=1,3
          zdvalHdeta(1:MAXEQNH,i) = zdvalHdeta(1:MAXEQNH,i) &
                                  + zdvalH(1:MAXEQNH,j)*dxdeta(j,i)
          zvalEeta(i,1:MAXEQNE) = zvalEeta(i,1:MAXEQNE) &
                                + zvalE(j,1:MAXEQNE)*dxdeta(j,i)
          zcvalEeta(i,1:MAXEQNE) = zcvalEeta(i,1:MAXEQNE) &
                            + detadx(i,j)*zcurlE(j,1:MAXEQNE)*rjacdxdeta
          zvalVeta(i,1:MAXEQNV) = zvalVeta(i,1:MAXEQNV) &
                     + detadx(i,j)*zvalV(j,1:MAXEQNV)*rjacdxdeta
        enddo
      enddo
!
!     compute the approximate solution at the point
      zdsolHdeta(1:MAXEQNH,1:3) = ZERO
      do kH=1,nrdofH
!
!       compute gradient wrt reference coordinates
        duHdeta(1:3) = gradH(1,kH)*dxideta(1,1:3) &
                     + gradH(2,kH)*dxideta(2,1:3) &
                     + gradH(3,kH)*dxideta(3,1:3)
        do i=1,3
          zdsolHdeta(1:MAXEQNH,i) = zdsolHdeta(1:MAXEQNH,i) &
                  + zdofH(1:MAXEQNH,kH)*duHdeta(i)
        enddo
      enddo
      zsolEeta(1:3,1:MAXEQNE) = ZERO
      zcsolEeta(1:3,1:MAXEQNE) = ZERO
      do kE=1,nrdofE
!
!       transform the shape functions and their curl to reference coordinates
        uEeta(1:3) = shapE(1,kE)*dxideta(1,1:3) &
                   + shapE(2,kE)*dxideta(2,1:3) &
                   + shapE(3,kE)*dxideta(3,1:3)
        curluEeta(1:3) = (detadxi(1:3,1)*curlE(1,kE) &
                       +  detadxi(1:3,2)*curlE(2,kE) &
                       +  detadxi(1:3,3)*curlE(3,kE))/rjac
        do i=1,3
          zsolEeta(i,1:MAXEQNE) = zsolEeta(i,1:MAXEQNE) &
                                + zdofE(1:MAXEQNE,kE)*uEeta(i)
          zcsolEeta(i,1:MAXEQNE) = zcsolEeta(i,1:MAXEQNE) &
                                 + zdofE(1:MAXEQNE,kE)*curluEeta(i)
        enddo
      enddo
      zsolVeta(1:3,1:MAXEQNV) = ZERO
      do kV=1,nrdofV
!
!       compute the shape function in reference coordinates
        uVeta(1:3) = (detadxi(1:3,1)*shapV(1,kV) &
                   +  detadxi(1:3,2)*shapV(2,kV) &
                   +  detadxi(1:3,3)*shapV(3,kV))/rjac
        do i=1,3
          zsolVeta(i,1:MAXEQNV) = zsolVeta(i,1:MAXEQNV) &
                                + zdofV(1:MAXEQNV,kV)*uVeta(i)
        enddo
      enddo
!
!     compute difference of exact and approximate solutions
      zddifHdeta = zdvalHdeta - zdsolHdeta
      zdifEeta = zvalEeta - zsolEeta
      zcdifEeta = zcvalEeta - zcsolEeta
      zdifVeta = zvalVeta - zsolVeta
!
!     eliminate the appropriate components
      do ivarH=1,MAXEQNH
        call dot_product(zddifHdeta(ivarH,1:3),rn, prod)
        zddifHdeta(ivarH,1:3) = zddifHdeta(ivarH,1:3) - prod*rn(1:3)
      enddo
      do ivarE=1,MAXEQNE
        call dot_product(zdifEeta(ivarE,1:3),rn, prod)
        zdifEeta(ivarE,1:3) = zdifEeta(ivarE,1:3) - prod*rn(1:3)
        call dot_product(zcdifEeta(ivarE,1:3),rn, prod)
        zcdifEeta(ivarE,1:3) = prod*rn(1:3)
      enddo
      do ivarV=1,MAXEQNV
        call dot_product(zdifVeta(ivarV,1:3),rn, prod)
        zdifVeta(ivarV,1:3) = prod*rn(1:3)
      enddo
!
!     accumulate for the error
      ivarH=0; ivarE=0; ivarV=0
!
!     loop through multiple copies of variables
      do j=1,NRCOMS
!
!       loop through all variables
        do i=1,NRINDEX
          select case(index(i))
          case(1)
            ivarH=ivarH+1
            derrorH(ivarH) = derrorH(ivarH) &
              + (abs(zddifHdeta(ivarH,1))**2 &
              +  abs(zddifHdeta(ivarH,2))**2 &
              +  abs(zddifHdeta(ivarH,3))**2)*weight
          case(2)
            ivarH=ivarH+1
          case(3)
            ivarE=ivarE+1
            derrorE(ivarE) = derrorE(ivarE) &
              + (abs(zdifEeta(ivarE,1))**2 &
              +  abs(zdifEeta(ivarE,2))**2 &
              +  abs(zdifEeta(ivarE,3))**2)*weight &
              + (abs(zcdifEeta(ivarE,1))**2 &
              +  abs(zcdifEeta(ivarE,2))**2 &
              +  abs(zcdifEeta(ivarE,3))**2)*weight
          case(4)
            ivarE=ivarE+1
          case(5)
            ivarV=ivarV+1
            derrorV(ivarV) = derrorV(ivarV) &
              + (abs(zdifVeta(ivarV,1))**2 &
              +  abs(zdifVeta(ivarV,2))**2 &
              +  abs(zdifVeta(ivarV,3))**2)*weight
          case(6)
            ivarV=ivarV+1
          end select
        enddo
      enddo
!
!   end of loop through integration points
    enddo
!
! end of loop through faces
  enddo
!
#if DEBUG_MODE
  if (iprint.ge.1) then
     write(*,*) 'ivarH,ivarE,ivarV = ',ivarH,ivarE,ivarV
     write(*,7210) derrorH(1:ivarH)
7210 format('elemBCIerror: derrorH = ',10e12.5)
     write(*,7220) derrorE(1:ivarE)
7220 format('elemBCIerror: derrorE = ',10e12.5)
     write(*,7230) derrorV(1:ivarV)
7230 format('elemBCIerror: derrorV = ',10e12.5)
     call pause
  endif
#endif
!
! sum up the errors componentwise
  ErrorH = 0.d0; ErrorE = 0.d0; ErrorV = 0.d0
  do i=1,ivarH
    ErrorH = ErrorH + derrorH(i)
  enddo
  do i=1,ivarE
    ErrorE = ErrorE + derrorE(i)
  enddo
  do i=1,ivarV
    ErrorV = ErrorV + derrorV(i)
  enddo
!
end subroutine elemBCIerror
