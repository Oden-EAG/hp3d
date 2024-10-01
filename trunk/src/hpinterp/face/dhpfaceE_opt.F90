!
#include "typedefs.h"
!
!-----------------------------------------------------------------------
!> @brief    determine H(curl) face dof interpolating H(curl) Dirichlet
!!           data using PB interpolation
!!
!> @param[in]  Mdle         - element (middle node) number
!> @param[in]  Iflag        - a flag specifying the GMP object
!!                            5 pris, 6 hexa, 7 tetr, 8 pyra
!> @param[in]  No           - number of a specific object
!> @param[in]  Etav         - GMP reference coordinates of the element vertices
!> @param[in]  Ntype        - element (middle node) type
!> @param[in]  Icase        - the face node case
!> @param[in]  Bcond        - the edge node BC flag
!> @param[in]  Nedge_orient - edge orientation
!> @param[in]  Nface_orient - face orientation
!> @param[in]  Norder       - element order
!> @param[in]  Iface        - face number
!> @param[in]  ZdofE        - H(curl) dof for the element (edge values only)
!!
!> @param[in,out] ZnodE     - H(curl) dof for the face
!!
!> @date Sep 2023
!-----------------------------------------------------------------------
subroutine dhpfaceE_opt(Mdle,Iflag,No,Etav,Ntype,Icase,Bcond,   &
                        Nedge_orient,Nface_orient,Norder,Iface, &
                        ZdofE, ZnodE)
  use control
  use parameters
  use physics
  use element_data
  use mpi_wrapper
!
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
  VTYPE,   dimension(MAXEQNE,MAXbrickE),     intent(in)    :: ZdofE
  VTYPE,   dimension(NRRHS*NREQNE(Icase),*), intent(inout) :: ZnodE
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
! derivatives of a H1 shape function wrt reference coordinates
  real(8), dimension(3)                 :: dvHdeta
!
! work space for shape3DE
  integer                               :: nrdofE
  real(8), dimension(3,MAXbrickE)       :: shapE
  real(8), dimension(3,MAXbrickE)       :: curlE
!
! H(curl) shape functions in reference coordinates
  real(8), dimension(3)                 :: uE,vE,curluE,curlvE
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
  VTYPE :: zvalEeta(3,MAXEQNE),zcurlEeta(3,MAXEQNE)
!
! linear systems
  real(8), dimension(MAXmdlqH+MAXmdlqE,MAXmdlqH+MAXmdlqE) :: aaE
  integer, dimension(MAXmdlqH+MAXmdlqE) :: ipivE
!
! load vector and solution
  VTYPE,   dimension(MAXmdlqH+MAXmdlqE,MAXEQNE) :: zbE
#if HP3D_COMPLEX
  real(8), dimension(MAXmdlqH+MAXmdlqE,MAXEQNE) :: duE_real, duE_imag
#endif
!
!..workspace for auxiliary matrix (storing info at integration points)
   real(8) :: A_E (MAXmdlqE,3*MAX_NINT2)
   real(8) :: A_CE(MAXmdlqE,3*MAX_NINT2)
   real(8) :: A_GH(MAXmdlqH,3*MAX_NINT2)
   integer :: lda,ldaE,ldaH,nE,nH,nda,noff
!
! decoded case and BC flag for the face node
  integer, dimension(NR_PHYSA)          :: ncase
  integer, dimension(NRINDEX_HEV)       :: ibcnd
!
! misc
  integer :: nrv,nre,nrf,nsign,nflag,i,j,jH,iE,jE,kjE,kjH,k,kE,info,ic
  integer :: ivarE,nvarE,ndofH_face,ndofE_face,ndofV_face,ndofQ_Face,ndofE_tot
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
7010 format('dhpfaceE: Mdle,Iflag,No,Icase,Iface,Type = ',5i4,2x,a4)
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
  call homogenD(TANGEN,Icase,Bcond, is_homD,ncase,ibcnd)
  if (is_homD) then
    zbE = ZERO
    goto 100
  endif
!
! if #dof is zero, return (nothing to do)
  if (ndofE_face.eq.0) return
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
7060 format('dhpfaceE: norder_ifc = ',20i4)
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
! initialize
  aaE = 0.d0; zbE = ZERO
!
! loop through integration points
  do l=1,nint
     t(1:2) = xi_list(1:2,l)
     wa     = wa_list(l)
!
!    get the corresponding master element coordinates and Jacobian
     call face_param(Ntype,Iface,t, xi,dxidt)
!
!    compute element H1 shape functions
     call shape3DH(Ntype,xi, &
                   norder_ifc,Nedge_orient,Nface_orient, &
                   nrdofH,shapH,gradH)
!
!    compute element Hcurl shape functions
     call shape3DE(Ntype,xi, &
                   norder_ifc,Nedge_orient,Nface_orient, &
                   nrdofE,shapE,curlE)
!
!    evaluate reference coordinates of the point as needed by GMP
!    brefgeom3D returns the outward normal as seen from the local element
     nsign = nsign_param(Ntype,Iface)
     call brefgeom3D(Mdle,xi,Etav,shapH(1:8),gradH(1:3,1:8),nrv,dxidt, &
                     nsign, eta,detadxi,dxideta,rjac,detadt,rn,bjac)
     weight = wa*bjac
!
#if HP3D_DEBUG
     if (iprint.eq.1) then
       write(*,7100) xi(1:2),eta(1:3),detadxi(1:3,1:3),rn(1:3),bjac
7100   format('dhpfaceE: xi,eta  = ',2f8.3,3x,3f8.3,/, &
              '          detadxi = ',3(3f8.3,2x),/, &
              '          rn,bjac = ',3f8.3,3x,f8.3)
     endif
#endif
!
!    call GMP routines to evaluate physical coordinates and their
!    derivatives wrt reference coordinates
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
!    compute inverse Jacobian (for transforming the curl)
     call geom(dxdeta, detadx,rjacdxdeta,nflag)
     if (nflag.ne.0) then
       write(*,*) 'dhpfaceE: rjacdxdeta = ',rjacdxdeta
       stop 1
     endif
!
!    get Dirichlet data at the point
     call dirichlet(Mdle,x,Icase, &
                    zvalH,zdvalH, zvalE,zdvalE, zvalV,zdvalV)
!
!    compute curl of the Dirichlet data
     zcurlE(1,1:MAXEQNE) = zdvalE(3,1:MAXEQNE,2) - zdvalE(2,1:MAXEQNE,3)
     zcurlE(2,1:MAXEQNE) = zdvalE(1,1:MAXEQNE,3) - zdvalE(3,1:MAXEQNE,1)
     zcurlE(3,1:MAXEQNE) = zdvalE(2,1:MAXEQNE,1) - zdvalE(1,1:MAXEQNE,2)
!
#if HP3D_DEBUG
     if (iprint.eq.1) then
       write(*,7130) x(1:3)
7130   format('dhpfaceE: x      = ',3f8.3)
       write(*,7140) zvalE(1:3,1:MAXEQNE)
7140   format('          zvalE  = ',2(3(2e10.3,2x),3x))
       write(*,7150) zcurlE(1:3,1:MAXEQNE)
7150   format('          zcurlE = ',2(3(2e10.3,2x),3x))
       write(*,1060) zdvalE(1:3,2,1:3)
1060   format('zdvalE(1:3,2,1:3) = ',3(2e12.5,2x),/,&
              '                    ',3(2e12.5,2x),/,&
              '                    ',3(2e12.5,2x))
     endif
#endif
!
!    Piola transform Dirichlet data to reference coordinates
     zvalEeta (1:3,1:MAXEQNE) = ZERO
     zcurlEeta(1:3,1:MAXEQNE) = ZERO
     do i=1,3
       do j=1,3
         zvalEeta(i,1:MAXEQNE)  = zvalEeta (i,1:MAXEQNE) &
                                + zvalE    (j,1:MAXEQNE)*dxdeta(j,i)
         zcurlEeta(i,1:MAXEQNE) = zcurlEeta(i,1:MAXEQNE) &
                                + zcurlE   (j,1:MAXEQNE)*detadx(i,j) &
                                                        *rjacdxdeta
       enddo
     enddo
!
#if HP3D_DEBUG
     if (iprint.eq.1) then
       write(*,*) 'zvalEeta, zcurlEeta BEFORE EDGE SUBTRACTION = '
       write(*,7110) zvalEeta(1:3,1:MAXEQNE)
7110   format(2(3(2e10.3,2x),3x))
       write(*,7110) zcurlEeta(1:3,1:MAXEQNE)
     endif
#endif
!
!    remove contributions from edges
!
!    loop through edge shape functions
     nrdofH = nrdofH - ndofH_face
     nrdofE = nrdofE - ndofE_face
     do kE=1,nrdofE
!
!      transform the shape functions and their curl to reference coordinates
       uE(1:3) = shapE(1,kE)*dxideta(1,1:3) &
               + shapE(2,kE)*dxideta(2,1:3) &
               + shapE(3,kE)*dxideta(3,1:3)
       curluE(1:3) = (detadxi(1:3,1)*curlE(1,kE) &
                   +  detadxi(1:3,2)*curlE(2,kE) &
                   +  detadxi(1:3,3)*curlE(3,kE))/rjac
       do i=1,3
         zvalEeta(i,1:MAXEQNE) = zvalEeta(i,1:MAXEQNE) &
                               - ZdofE(1:MAXEQNE,kE)*uE(i)
         zcurlEeta(i,1:MAXEQNE) = zcurlEeta(i,1:MAXEQNE) &
                               - ZdofE(1:MAXEQNE,kE)*curluE(i)
       enddo
     enddo
!
#if HP3D_DEBUG
     if (iprint.eq.1) then
       write(*,*) 'zvalEeta, zcurlEeta AFTER EDGE SUBTRACTION = '
       write(*,7110) zvalEeta(1:3,1:MAXEQNE)
       write(*,7110) zcurlEeta(1:3,1:MAXEQNE)
       call pause
     endif
#endif
!
!    offset in auxiliary matrix
     noff = 3*(l-1)
!
!    loop through element face test H(curl) functions
     do jE=1,ndofE_face
       kjE = nrdofE + jE
!
!      transform the shape functions and their curl to reference coordinates
       vE(1:3) = shapE(1,kjE)*dxideta(1,1:3) &
               + shapE(2,kjE)*dxideta(2,1:3) &
               + shapE(3,kjE)*dxideta(3,1:3)
       curlvE(1:3) = (detadxi(1:3,1)*curlE(1,kjE) &
                   +  detadxi(1:3,2)*curlE(2,kjE) &
                   +  detadxi(1:3,3)*curlE(3,kjE))/rjac
!
!      subtract normal component
       call dot_product(vE,rn, prod)
       vE(1:3) = vE(1:3) - prod*rn(1:3)
!
!      leave only the normal component
       call dot_product(curlvE,rn, prod)
       curlvE(1:3) = prod*rn(1:3)
!
!      accumulate for the load vector
       zbE(jE,1:MAXEQNE) = zbE(jE,1:MAXEQNE) &
                         + (curlvE(1)*zcurlEeta(1,1:MAXEQNE) &
                         +  curlvE(2)*zcurlEeta(2,1:MAXEQNE) &
                         +  curlvE(3)*zcurlEeta(3,1:MAXEQNE))*weight
!
!      fill auxiliary matrices for stiffness (E, curl E)
       A_E (jE,noff+1:noff+3) = vE(1:3)     * sqrt(weight)
       A_CE(jE,noff+1:noff+3) = curlvE(1:3) * sqrt(weight)
!
!    end of loop through H(curl) test functions
     enddo
!
!    loop through element face H1 test functions
     do jH=1,ndofH_face
       kjH = nrdofH + jH
!
!      compute gradient wrt reference coordinates
       dvHdeta(1:3) = gradH(1,kjH)*dxideta(1,1:3) &
                    + gradH(2,kjH)*dxideta(2,1:3) &
                    + gradH(3,kjH)*dxideta(3,1:3)
!
!      subtract normal component
       call dot_product(dvHdeta,rn, prod)
       dvHdeta(1:3) = dvHdeta(1:3) - prod*rn(1:3)
!
!      accumulate for the load vector
       zbE(ndofE_face+jH,1:MAXEQNE) = zbE(ndofE_face+jH,1:MAXEQNE) &
                   + (zvalEeta(1,1:MAXEQNE)*dvHdeta(1) &
                   +  zvalEeta(2,1:MAXEQNE)*dvHdeta(2) &
                   +  zvalEeta(3,1:MAXEQNE)*dvHdeta(3))*weight
!
!      fill auxiliary matrix for stiffness (grad H)
       A_GH(jH,noff+1:noff+3) = dvHdeta(1:3) * sqrt(weight)
!
!    end of loop through H1 test functions
     enddo
!
! end of loop through integration points
  enddo
!
! compute stiffness matrix
   nda = 3*nint
   nE = ndofE_face
   nH = ndofH_face
   ldaE = MAXmdlqE
   ldaH = MAXmdlqH
   lda = ldaE + ldaH
!
!  aaE[1,1] = A_CE * A_CE^T (curl E, curl E)
!  - compute upper triangle
   call DSYRK('U','N',nE,nda,                  &
              1.0d0,A_CE(1:ldaE,1:nda), ldaE,  &
              0.0d0,aaE (1:lda ,1:nE ), lda)
!  aaE[1,2] = A_E * A_GH^T  (E, grad u)
   call DGEMM('N','T',nE,nH,nda,                     &
               1.0d0,A_E (1:ldaE,   1:nda  ), ldaE,  &
                     A_GH(1:ldaH,   1:nda  ), ldaH,  &
               0.0d0,aaE (1:lda ,nE+1:nE+nH), lda)
!  use symmetry to complete the stiffness matrix
!  - fill lower triangle of aaE[1,1]
   do iE=1,nE
      aaE(iE+1:nE,iE) = aaE(iE,iE+1:nE)
   enddo
!  - aaE[2,1] = aaE[1,2]^T
   do jH=1,nH
      do iE=1,nE
         aaE(nE+jH,iE) = aaE(iE,nE+jH)
      enddo
   enddo
!
   ndofE_tot = ndofE_face + ndofH_face
!
#if HP3D_DEBUG
  if (iprint.ge.1) then
    write(*,*) 'dhpfaceE: LOAD VECTOR AND STIFFNESS MATRIX FOR ', &
               'ndofE_face,ndofH_face = ',ndofE_face,ndofH_face
    do j=1,ndofE_tot
      write(*,7015) j, zbE(j,1:MAXEQNE)
      write(*,7016) aaE(j,1:ndofE_tot)
    enddo
  endif
#if HP3D_COMPLEX
  7015 format(i5,2x,6(2e10.3,2x))
#else
  7015 format(i5,2x,10e12.5)
#endif
  7016 format(10e12.5)
#endif
!
!------------------------------------------------------
!
! solve the linear system
!
! lu pivot
  call DGETRF(ndofE_tot,ndofE_tot,aaE,lda,ipivE,info)
  if (info.ne.0) then
     write(*,*)'dhpfaceE: DGETRF RETURNED INFO = ',info
     call logic_error(FAILURE,__FILE__,__LINE__)
  endif
!
#if HP3D_COMPLEX
! apply pivots to load vector
  call ZLASWP(MAXEQNE,zbE,lda,1,ndofE_tot,ipivE,1)
!
  duE_real(1:ndofE_tot,:) =  real(zbE(1:ndofE_tot,:))
  duE_imag(1:ndofE_tot,:) = aimag(zbE(1:ndofE_tot,:))
!
! triangular solves
  call DTRSM('L','L','N','U',ndofE_tot,MAXEQNE,1.d0,aaE,lda, duE_real,lda)
  call DTRSM('L','U','N','N',ndofE_tot,MAXEQNE,1.d0,aaE,lda, duE_real,lda)
!
  call DTRSM('L','L','N','U',ndofE_tot,MAXEQNE,1.d0,aaE,lda, duE_imag,lda)
  call DTRSM('L','U','N','N',ndofE_tot,MAXEQNE,1.d0,aaE,lda, duE_imag,lda)
!
  zbE(1:ndofE_tot,:) = dcmplx(duE_real(1:ndofE_tot,:), duE_imag(1:ndofE_tot,:))
#else
! apply pivots to load vector
  call DLASWP(MAXEQNE,zbE,lda,1,ndofE_tot,ipivE,1)
! triangular solves
  call DTRSM('L','L','N','U',ndofE_tot,MAXEQNE,1.d0,aaE,lda, zbE,lda)
  call DTRSM('L','U','N','N',ndofE_tot,MAXEQNE,1.d0,aaE,lda, zbE,lda)
#endif
!
#if HP3D_DEBUG
  if (iprint.ge.1) then
     write(*,*) 'dhpfaceE: k,zbE(k) = '
     do k=1,ndofE_tot
        write(*,7015) k,zbE(k,1:MAXEQNE)
     enddo
     call pause
  endif
#endif
!
!------------------------------------------------------
!
!  ...save the DOFs, skipping irrelevant entries
  100 continue
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,*) 'dhpfaceE: ncase = ', ncase
      endif
#endif
!
!  ...initialize global variable counter, and node local variable counter
      ivarE=0; nvarE=0
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
!  .........H(curl) component
            case(TANGEN)
!
!  ...........update global counter
              ivarE=ivarE+1
!
!  ...........if the variable is supported by the node
              if (ncase(i).eq.1) then
!
!  .............update node local counter
                nvarE = nvarE + 1
!
!  .............do not write dof if physics attribute is deactivated
                if (.not. PHYSAm(i)) exit
!
!  .............store Dirichlet dof
                if (ibcnd(ic).eq.1) ZnodE(nvarE,1:ndofE_face) = zbE(1:ndofE_face,ivarE)
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
!   write(*,11) 'dhpfaceE_opt: ', end_time-start_time
!11 format(A,f12.5,' s')
!   !$OMP END CRITICAL
!
#if HP3D_DEBUG
      if (iprint.eq.1) call result
#endif
!
end subroutine dhpfaceE_opt
