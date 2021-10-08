!--------------------------------------------------------------------
!
!     routine name      - elem_residual
!
!--------------------------------------------------------------------
!
!     latest revision:  - Oct 2021
!
!     purpose:          - routine returns element residual (squared)
!                         for the Primal Poisson and UW Time Harmonic
!                         Maxwell equation
!
!     arguments:
!        in:
!             Mdle      - element middle node number
!        out:
!             Resid     - element residual (squared)
!             Nref_flag - suggested h-refinement flag
!
!---------------------------------------------------------------------
!
#include "typedefs.h"
!
subroutine elem_residual(Mdle, Resid,Nref_flag)
!..modules used
   use control
   use data_structure3D
   use element_data
   use parametersDPG
!..no implicit statements
   implicit none
!..declare input/output variables
   integer, intent(in)  :: Mdle
   real(8), intent(out) :: Resid
   integer, intent(out) :: Nref_flag
!
!..number of test and trial degrees of freedom
   integer :: nrdofH ,nrdofE ,nrdofV ,nrdofQ
   integer :: nrdofHH,nrdofEE,nrdofVV,nrdofQQ
   integer :: nrTest
!
!..(enriched) order of element nodes
   integer :: norder(19),norderP(19),nordP
!
!..element type
   character(len=4) :: etype
!
!---------------------------------------------------------------------
!
   norder (1:19) = 0
   norderP(1:19) = 0
!
   etype = NODES(Mdle)%type
!..determine order of approximation
   call find_order(Mdle, norder)
!..set the enriched order of appoximation
   select case(etype)
      case('mdlb')
         nordP = NODES(Mdle)%order+NORD_ADD*111
      case('mdlp')
         nordP = NODES(Mdle)%order+NORD_ADD*11
      case('mdln','mdld')
         nordP = NODES(Mdle)%order+NORD_ADD
      case default
         write(*,*) 'elem_residual: invalid etype param. stop.'
         stop
   end select
!..note: compute_enriched_order works is only implemented for hexa currently
   call compute_enriched_order(etype,nordP, norderP)
!..compute nrdof for trial
   call celndof(etype,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
!..compute nrdof for test
   call celndof(etype,norderP, nrdofHH,nrdofEE,nrdofVV,nrdofQQ)
!
   nrTest = nrdofHH
   call elem_residual_poisson(           &
      Mdle,nrTest,nrdofHH,nrdofH,nrdofV, &
      Resid,Nref_flag)
!
end subroutine elem_residual

!--------------------------------------------------------------------
!
!     routine name      - elem_residual_poisson
!
!--------------------------------------------------------------------
!
!     latest revision:  - Oct 2021
!
!     purpose:          - routine returns element residual (squared)
!                         for the Primal DPG Poisson equation
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
subroutine elem_residual_poisson(Mdle,                &
                                 NrTest,NrdofHH,      &
                                 NrdofH,NrdofV,       &
                                 Resid,Nref_flag)
!..modules used
   use control
   use parametersDPG
   use element_data
   use data_structure3D
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
   character(len=4) :: etype,ftype
!
!..element order, orientation for edges and faces
   integer :: norder(19), norient_edge(12), norient_face(6)
!
!..face order
   integer :: norderf(5)
!
!..geometry dof (work space for nodcor)
   real(8) :: xnod(3,MAXbrickH)
!
!..solution dof (work space for solelm)
   real(8) :: dofH(MAXEQNH,MAXbrickH)
   real(8) :: dofE(MAXEQNE,MAXbrickE)
   real(8) :: dofV(MAXEQNV,MAXbrickV)
   real(8) :: dofQ(MAXEQNQ,MAXbrickQ)
!
!..variables for geometry
   real(8) :: xi(3), x(3), dxdxi(3,3), dxidx(3,3), rn(3)
   real(8) :: dxidt(3,2), dxdt(3,2), t(2)
!
!..H1 shape functions
   real(8) :: shapH(MAXbrickH), gradH(3,MAXbrickH)
!
!..H(div) shape functions
   real(8) :: shapV(3,MAXbrickV), divV(MAXbrickV)
!
!..Enriched H1 shape functions
   real(8) :: shapHH(MAXbrickHH), gradHH(3,MAXbrickHH)
!
!..Gram matrix in packed format
   real(8), allocatable :: gramP(:)
!
!..load vector for the enriched space
   real(8) :: bload_H(NrTest),bload_Hc(NrTest)
!
!..3D quadrature data
   real(8) :: xiloc(3,MAXNINT3ADD), waloc(MAXNINT3ADD)
!
!..2D quadrature data
   real(8) :: tloc(2,MAXNINT2ADD), wtloc(MAXNINT2ADD)
!
!..test functions and gradients
   real(8) :: q, v, aux
   real(8) :: dq(3), dv(3)
!
!..approximate solution
   real(8) :: rgradHxi(3), rgradH(3)
   real(8) :: solVxi(3), solV(3)
   real(8) :: solVn
!
!..various variables for the problem
   real(8) :: rjac, bjac, fval, wa, weight
   integer :: iflag, nrf, nint
   integer :: k1, k2, k, l
   integer :: nordP, nrdof, nsign, ifc, info
!
#if DEBUG_MODE
   integer :: iprint = 0
#endif
!
!..for Gram matrix compressed storage format
   integer :: nk
   nk(k1,k2) = (k2-1)*k2/2+k1
!
!-----------------------------------------------------------------------
!
!..element type
   etype = NODES(Mdle)%type
   nrf = nface(etype)
!
!..determine order of approximation
   call find_order(Mdle, norder)
!
!  ...set the enriched order of approximation
   select case(etype)
      case('mdlb'); nordP = NODES(Mdle)%order+NORD_ADD*111
      case('mdln','mdld'); nordP = NODES(Mdle)%order+NORD_ADD
      case('mdlp'); nordP = NODES(Mdle)%order+NORD_ADD*11
   end select
!
!..determine edge and face orientations
   call find_orient( Mdle, norient_edge,norient_face)
!
!..determine nodes coordinates
   call nodcor(Mdle, xnod)
!
!..determine solution dof
   call solelm(Mdle, dofH,dofE,dofV,dofQ)
!
!..allocate space for auxiliary matrices
   allocate(gramP(NrTest*(NrTest+1)/2))
!
!..clear space for auxiliary matrices
   bload_H = ZERO; gramP = ZERO
!
!-----------------------------------------------------------------------
!.......... element integrals
!-----------------------------------------------------------------------
!
!..use the enriched order to set the quadrature
   INTEGRATION = NORD_ADD
   call set_3D_int_DPG(etype,norder,norient_face, nint,xiloc,waloc)
   INTEGRATION = 0
!
!..loop over integration points
   do l=1,nint
      xi(1:3) = xiloc(1:3,l)
      wa = waloc(l)
!
!  ...determine element H1 shape functions
      call shape3DH(etype,xi,norder,norient_edge,norient_face, &
                    nrdof,shapH,gradH)
!
!  ...determine discontinuous H1 shape functions
      call shape3HH(etype,xi,nordP, nrdof,shapHH,gradHH)
!
!  ...geometry
      call geom3D(Mdle,xi,xnod,shapH,gradH,NrdofH, &
                  x,dxdxi,dxidx,rjac,iflag)
!
!  ...integration weight
      weight = rjac*wa
!
!  ...compute the approximate solution
      rgradHxi(1:3) = 0.d0
      do k=1,NrdofH
         rgradHxi(1:3) = rgradHxi(1:3) + dofH(1,k)*gradH(1:3,k)
      enddo
      rgradH(1:3) = rgradHxi(1)*dxidx(1,1:3) &
                  + rgradHxi(2)*dxidx(2,1:3) &
                  + rgradHxi(3)*dxidx(3,1:3)
!
!  ...get the RHS
      call getf(Mdle,x, fval)
!
!  ...1st loop through enriched H1 test functions
      do k1=1,NrdofHH
!     ...Piola transformation
         v = shapHH(k1)
         dv(1:3) = gradHH(1,k1)*dxidx(1,1:3) &
                 + gradHH(2,k1)*dxidx(2,1:3) &
                 + gradHH(3,k1)*dxidx(3,1:3)
!
!     ...residual volume terms: (grad u, grad v) - (f,v)
         aux = (rgradH(1)*dv(1) + rgradH(2)*dv(2) + rgradH(3)*dv(3)) - fval*v
         bload_H(k1) = bload_H(k1) + aux*weight
!
!     ...2nd loop through enriched H1 test functions
         do k2=k1,NrdofHH
!        ...Piola transformation
            q = shapHH(k2)
            dq(1:3) = gradHH(1,k2)*dxidx(1,1:3) &
                    + gradHH(2,k2)*dxidx(2,1:3) &
                    + gradHH(3,k2)*dxidx(3,1:3)
!
!        ...determine index in triangular packed format
            k = nk(k1,k2)
!
!        ...H1 test inner product: (q,v) + (grad_h q, grad_h v)
            aux = q*v + (dq(1)*dv(1) + dq(2)*dv(2) + dq(3)*dv(3))
            gramP(k) = gramP(k) + aux*weight
!
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
!
!     ...determine element H1 shape functions (for geometry)
         call shape3DH(etype,xi,norder,norient_edge,norient_face, &
                       nrdof,shapH,gradH)
!
!     ...determine element Hdiv shape functions (for fluxes)
         call shape3DV(etype,xi,norder,norient_face, &
                       nrdof,shapV,divV)
!
!     ...geometry
         call bgeom3D(Mdle,xi,xnod,shapH,gradH,NrdofH,dxidt,nsign, &
                      x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
         weight = bjac*wtloc(l)
!
!     ...compute approximate flux at the point
         solVxi(1:3) = 0.d0
         do k=1,nrdofV
            solVxi(1:3) = solVxi(1:3) + dofV(1,k)*shapV(1:3,k)
         enddo
!     ...Piola transformation
         solV(1:3) = dxdxi(1:3,1)*solVxi(1) &
                   + dxdxi(1:3,2)*solVxi(2) &
                   + dxdxi(1:3,3)*solVxi(3)
         solV(1:3) = solV(1:3) / rjac
!     ...normal component
         solVn = solV(1)*rn(1) + solV(2)*rn(2) + solV(3)*rn(3)
!
!     ...loop through enriched test functions
         do k1=1,NrdofHH
            v = shapHH(k1)
!        ...accumulate for the load vector
            bload_H(k1) = bload_H(k1) - solVn*v*weight
         enddo
      enddo
   enddo
!
!-----------------------------------------------------------------------
!
!..factorize the test stiffness matrix
   call DPPTRF('U', NrTest, gramP, info)
   if (info.ne.0) then
      write(*,*) 'elem_residual_poisson: info = ',info
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
      write(*,*) 'elem_residual_poisson: info = ',info
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
!
!..set recommended h-refinement flag for this hexa element
   Nref_flag = 111
!
#if DEBUG_MODE
   if (iprint.ge.1) then
      write(*,7010) Mdle, Resid
 7010 format('elem_residual_poisson: Mdle, Resid = ',i5,3x,e12.5)
   endif
#endif
!
end subroutine elem_residual_poisson

