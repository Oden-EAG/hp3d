!--------------------------------------------------------------------
!
!     routine name      - elem_residual_heat
!
!--------------------------------------------------------------------
!
!     latest revision:  - Mar 2019
!
!     purpose:          - routine returns element residual (squared)
!                         for the Primal Heat equation (Poisson)
!
!     arguments:
!        in:
!             Mdle      - element middle node number
!             NrTest    - total number of test dof
!             NrdofHH   - number of H1 test dof
!             NrdofH    - number of H1 trial dof
!             NrdofV    - number of H(div) trial dof
!        out:
!             Resid     - element residual (squared)
!             Nref_flag - suggested h-refinement flag
!
!---------------------------------------------------------------------
!
#include "typedefs.h"
!
subroutine elem_residual_heat(Mdle,                &
                              NrTest,NrdofHH,      &
                              NrdofH,NrdofV,       &
                              Resid,Nref_flag)
!..modules used
   use control
   use parametersDPG
   use element_data
   use data_structure3D
   use commonParam
   use laserParam
!..no implicit statements
   implicit none
!..declare input/output variables
   integer, intent(in)  :: Mdle
   integer, intent(in)  :: NrTest
   integer, intent(in)  :: NrdofHH
   integer, intent(in)  :: NrdofH
   integer, intent(in)  :: NrdofV
   real(8), intent(out) :: Resid
   integer, intent(out) :: Nref_flag
!
!..declare edge/face type variables
   integer :: etype,ftype
!
!..declare element order, orientation for edges and faces
   integer, dimension(19)  :: norder
   integer, dimension(12)  :: norient_edge
   integer, dimension(6)   :: norient_face
!..face order
   integer, dimension(5)   :: norderf
!
!..geometry dof (work space for nodcor)
   real(8) :: xnod(3,MAXbrickH)
!
!..solution dof (work space for solelm)
   VTYPE, dimension(MAXEQNH,MAXbrickH) :: zdofH
   VTYPE, dimension(MAXEQNE,MAXbrickE) :: zdofE
   VTYPE, dimension(MAXEQNV,MAXbrickV) :: zdofV
   VTYPE, dimension(MAXEQNQ,MAXbrickQ) :: zdofQ
!..current H1 solution
   real(8) :: rsolH
!
!..variables for geometry
   real(8), dimension(3)    :: xi,x,rn
   real(8), dimension(3,2)  :: dxidt,dxdt,rt
   real(8), dimension(3,3)  :: dxdxi,dxidx
   real(8), dimension(2)    :: t
!
!..H1 shape functions
   real(8), dimension(MAXbrickH)    :: shapH
   real(8), dimension(3,MAXbrickH)  :: gradH
!
!..enriched H1 shape functions
   real(8), dimension(MAXbrickHH)   :: shapHH
   real(8), dimension(3,MAXbrickHH) :: gradHH
!
!..H(div) shape functions
   real(8), dimension(3,MAXbrickV)  :: shapV
   real(8), dimension(MAXbrickV)    :: divV
!
!..test functions and gradients
   real(8) :: v1,v2,v2n
   real(8), dimension(3) :: dv1,dv2
!
!..Gram matrix in packed format
   !real(8) :: gramP(NrTest*(NrTest+1)/2)
   real(8), allocatable :: gramP(:)
!
!..load vector for the enriched space
   real(8), dimension(NrTest) :: bload_H,bload_Hc
!
!..3D quadrature data
   real(8), dimension(3,MAXNINT3ADD)  :: xiloc
   real(8), dimension(MAXNINT3ADD)    :: waloc
!
!..2D quadrature data
   real(8), dimension(2,MAXNINT2ADD)  :: tloc
   real(8), dimension(MAXNINT2ADD)    :: wtloc
!
!..approximate solution
   real(8), dimension(3) :: rgradHxi,rgradH
   real(8), dimension(3) :: rsolVxi,rsolV
   real(8)               :: rfval,rsolVn
   VTYPE                 :: zfval
!
   VTYPE :: zvoid(3)
!
!..number of faces per element type
   integer :: nrf
!..various variables for the problem
   real(8) :: rjac,weight,wa
   real(8) :: bjac
   integer :: i1,i2,j1,j2,k1,k2,kH,kk,i,j,nint,nint3,iflag,kE,k,l
   integer :: nordP,nsign,ifc,info,icomp,nrdof
!
   integer, external :: ij_upper_to_packed
!
#if DEBUG_MODE
   integer :: iprint
   iprint = 0
#endif
!
!-----------------------------------------------------------------------
!
!..element type
   etype = NODES(Mdle)%ntype
   nrf = nface(etype)
!
!..determine order of approximation
   call find_order(Mdle, norder)
!
!  ...set the enriched order of approximation
   select case(etype)
      case(MDLB)     ; nordP = NODES(Mdle)%order+NORD_ADD*111
      case(MDLP)     ; nordP = NODES(Mdle)%order+NORD_ADD*11
      case(MDLN,MDLD); nordP = NODES(Mdle)%order+NORD_ADD
   end select
!
!..determine edge and face orientations
   call find_orient( Mdle, norient_edge,norient_face)
!
!..determine nodes coordinates
   call nodcor(Mdle, xnod)
!
!..determine solution dof
   call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
!
#if DEBUG_MODE
   if (iprint.eq.1) then
      write(*,7020) xnod(1,1:8),xnod(2,1:8),xnod(3,1:8)
7020  format('elem_residual_heat: xnod  = ',8(f8.3,2x), &
       2(  /,'                            ',8(f8.3,2x)))
      write(*,7030) zdofH(1,1:8),zdofV(1,1:6)
7030  format('elem_residual_heat: zdofH = ',8(e12.5,2x), &
           /,'                    zdofV = ',6(e12.5,2x))
   endif
#endif
!
!..allocate space for auxiliary matrices
   allocate(gramP(NrTest*(NrTest+1)/2))
!
!..clear space for auxiliary matrices
   bload_H = rZERO; gramP = rZERO
!
!-----------------------------------------------------------------------
!.......... element integrals
!-----------------------------------------------------------------------
!
!..use the enriched order to set the quadrature
   INTEGRATION = NORD_ADD
   call set_3D_int_DPG(etype,norder,norient_face, nint3,xiloc,waloc)
   INTEGRATION = 0
!
!..loop over integration points
   do l=1,nint3
      xi(1:3) = xiloc(1:3,l)
      wa = waloc(l)
!
!  ...determine element H1 shape functions
      call shape3DH(etype,xi,norder,norient_edge,norient_face, &
                    nrdof,shapH,gradH)
#if DEBUG_MODE
      if (nrdof .ne. NrdofH) then
         write(*,*) 'elem_residual_heat: INCONSISTENCY NrdofH. stop.'
         stop
      endif
#endif
!
!  ...determine discontinuous H1 shape functions
      call shape3HH(etype,xi,nordP, nrdof,shapHH,gradHH)
#if DEBUG_MODE
      if (nrdof .ne. NrdofHH) then
         write(*,*) 'elem_residual_heat: INCONSISTENCY NrdofHH. stop.'
         stop
      endif
#endif
!
!  ...geometry
      call geom3D(Mdle,xi,xnod,shapH,gradH,NrdofH, &
                  x,dxdxi,dxidx,rjac,iflag)
!
!  ...integration weight
      weight = rjac*wa
!
!  ...compute the approximate solution
      rsolH = 0.d0; rgradHxi(1:3) = 0.d0
      do k=1,NrdofH
         rsolH = rsolH + real(zdofH(1,k))*shapH(k)
         rgradHxi(1:3) = rgradHxi(1:3) + real(zdofH(1,k))*gradH(1:3,k)
      enddo
      rgradH(1:3) = rgradHxi(1)*dxidx(1,1:3) &
                  + rgradHxi(2)*dxidx(2,1:3) &
                  + rgradHxi(3)*dxidx(3,1:3)
!
!  ...get the RHS
      call getf(Mdle,x, zfval,zvoid)
      rfval = real(zfval)
!
!  ...1st loop through enriched H1 test functions
      do k1=1,NrdofHH
!
         v1 = shapHH(k1)
         dv1(1:3) = gradHH(1,k1)*dxidx(1,1:3) &
                  + gradHH(2,k1)*dxidx(2,1:3) &
                  + gradHH(3,k1)*dxidx(3,1:3)
!     ...account for short fiber scaling (anisotropic diffusion operator)
         if (ANISO_HEAT .eq. 1) then
            rgradH(3) = ALPHA_Z*ALPHA_Z*rgradH(3)
         elseif (NO_PROBLEM.eq.2 .and. USE_PML .and. &
                 x(3).ge.(ZL-2.0*PML_FRAC*ZL)) then
!        ...reduce artificial cooling from PML
            rgradH(3) = ALPHA_Z*ALPHA_Z*rgradH(3)
         endif
!
!     ...residual for single step of heat equation
!     ...with thermal load coming from fval (getf)
         bload_H(k1) = bload_H(k1) &
                     + (DELTA_T*ALPHA_0*        &
                        ( rgradH(1)*dv1(1) +    &
                          rgradH(2)*dv1(2) +    &
                          rgradH(3)*dv1(3)      &
                        ) + rsolH*v1-rfval*v1   &
                       )*weight
!
!     ...2nd loop through enriched H1 test functions
         do k2=k1,NrdofHH
            v2 = shapHH(k2)
            dv2(1:3) = gradHH(1,k2)*dxidx(1,1:3) &
                     + gradHH(2,k2)*dxidx(2,1:3) &
                     + gradHH(3,k2)*dxidx(3,1:3)
!
!        ...accumulate for the test stiffness matrix
            k = ij_upper_to_packed(k1,k2)
            select case(INNER_PRODUCT)
               case(1)
                  if (ANISO_HEAT .eq. 1) then
                     dv2(3) = ALPHA_Z*ALPHA_Z*dv2(3)
                  elseif (NO_PROBLEM.eq.2 .and. USE_PML .and. &
                          x(3).ge.(ZL-2.0*PML_FRAC*ZL)) then
!                 ...reduce artificial cooling from PML
                     dv2(3) = ALPHA_Z*ALPHA_Z*dv2(3)
                  endif
                  gramP(k) = gramP(k) &
                           + ( v1*v2 + DELTA_T*ALPHA_0 * &
                                ( dv1(1)*dv2(1) + &
                                  dv1(2)*dv2(2) + &
                                  dv1(3)*dv2(3) ) &
                              ) *weight
            end select
!     ...end 1st loop through enriched H1 test functions
         enddo
!  ...end 2nd loop through enriched H1 test functions
      enddo
!..end loop through integration points
   enddo
!
!
!-----------------------------------------------------------------------
!.......... boundary integrals
!-----------------------------------------------------------------------
!
!..loop through element faces
   do ifc=1,nrf
!
!  ...sign factor to determine the OUTWARD normal unit vector
      nsign = nsign_param(etype,ifc)
!
!  ...face type
      ftype = face_type(etype,ifc)
!
!  ...face order of approximation
      call face_order(etype,ifc,norder, norderf)
!
!  ...set 2D quadrature
      INTEGRATION = NORD_ADD
      call set_2D_int_DPG(ftype,norderf,norient_face(ifc), nint,tloc,wtloc)
      INTEGRATION = 0
!
!  ...loop through integration points
      do l=1,nint
!
!     ...face coordinates
         t(1:2) = tloc(1:2,l)
!
!     ...face parametrization
         call face_param(etype,ifc,t, xi,dxidt)
!
!     ...determine discontinuous H1 shape functions
         call shape3HH(etype,xi,nordP, nrdof,shapHH,gradHH)
#if DEBUG_MODE
         if (nrdof .ne. NrdofHH) then
            write(*,*) 'elem_residual_heat: INCONSISTENCY NrdofHH. stop.'
            stop
         endif
#endif
!
!     ...determine element H1 shape functions (for geometry)
         call shape3DH(etype,xi,norder,norient_edge,norient_face, &
                       nrdof,shapH,gradH)
#if DEBUG_MODE
         if (nrdof .ne. NrdofH) then
            write(*,*) 'elem_residual_heat: INCONSISTENCY NrdofH. stop.'
            stop
         endif
#endif
!
!     ...determine element Hdiv shape functions (for fluxes)
         call shape3DV(etype,xi,norder,norient_face, &
                       nrdof,shapV,divV)
#if DEBUG_MODE
         if (nrdof .ne. NrdofV) then
            write(*,*) 'elem_residual_heat: INCONSISTENCY NrdofV. stop.'
            stop
         endif
#endif
!
!     ...geometry
         call bgeom3D(Mdle,xi,xnod,shapH,gradH,NrdofH,dxidt,nsign, &
                      x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
         weight = bjac*wtloc(l)
!
!     ...compute approximate flux at the point
         rsolVxi(1:3) = 0.d0
         do k=1,nrdofV
            rsolVxi(1:3) = rsolVxi(1:3) + real(zdofV(1,k))*shapV(1:3,k)
         enddo
         rsolV(1:3) = (dxdxi(1:3,1)*rsolVxi(1) &
                      +dxdxi(1:3,2)*rsolVxi(2) &
                      +dxdxi(1:3,3)*rsolVxi(3))/rjac
         rsolVn = rsolV(1)*rn(1)+rsolV(2)*rn(2)+rsolV(3)*rn(3)
!
!     ...loop through enriched test functions
         do k1=1,NrdofHH
            v1 = shapHH(k1)
!        ...accumulate for the load vector
            bload_H(k1) = bload_H(k1) - DELTA_T*ALPHA_0*rsolVn*v1*weight
         enddo
      enddo
   enddo
!
#if DEBUG_MODE
   if (iprint.ge.1) then
      write(*,7015) bload_H(1:NrTest)
7015  format('elem_residual_heat: FINAL bload_H = ',10(/,10(e12.5,2x)))
      call pause
   endif
#endif
!
!-----------------------------------------------------------------------
!
!..factorize the test stiffness matrix
   call DPPTRF('U', NrTest, gramP, info)
   if (info.ne.0) then
      write(*,*) 'elem_residual_heat: info = ',info
      stop
   endif
!
!..save copies of the RHS to compute later the residual
   bload_Hc = bload_H
!
!..compute the product of inverted test Gram matrix with RHS,
!..bload_H is overwritten with the solution
   call DPPTRS('U', NrTest, 1, gramP, bload_H, NrTest, info)
   if (info.ne.0) then
      write(*,*) 'elem_residual_heat: info = ',info
      stop
   endif
!
   deallocate(gramP)
!
!..compute the residual
   Resid = 0.d0
   do k=1,NrdofHH
      Resid = Resid + bload_Hc(k)*bload_H(k)
   enddo
   Nref_flag = 111
!
#if DEBUG_MODE
   if (iprint.ge.1) then
      write(*,7010) Mdle, Resid
 7010 format('elem_residual_heat: Mdle, Resid = ',i5,3x,e12.5)
      call pause
   endif
#endif
!
end subroutine elem_residual_heat

