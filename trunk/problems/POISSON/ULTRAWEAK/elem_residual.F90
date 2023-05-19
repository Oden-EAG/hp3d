!--------------------------------------------------------------------
!
!     routine name      - elem_residual
!
!--------------------------------------------------------------------
!
!     latest revision:  - May 2023
!
!     purpose:          - routine returns element residual (squared)
!                         for the Ultraweak Poisson.
!
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
!
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
   integer :: etype
!
!---------------------------------------------------------------------
!
   norder (1:19) = 0
   norderP(1:19) = 0
!
   etype = NODES(Mdle)%ntype
!..determine order of approximation
   call find_order(Mdle, norder)
!..set the enriched order of appoximation
   select case(etype)
      case(MDLB)
         nordP = NODES(Mdle)%order+NORD_ADD*111
      case(MDLP)
         nordP = NODES(Mdle)%order+NORD_ADD*11
      case(MDLN,MDLD)
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
   nrTest = nrdofHH + nrdofVV
   call elem_residual_poisson_UW(        &
      Mdle,nrTest,nrdofHH,nrdofVV,nrdofQ,nrdofH,nrdofV, &
      Resid,Nref_flag)
!
end subroutine elem_residual

!--------------------------------------------------------------------
!
!     routine name      - elem_residual_poisson_UW
!
!--------------------------------------------------------------------
!
!     latest revision:  - May 2023
!
!     purpose:          - routine returns element residual (squared)
!                         for the UltraWeak DPG formulation of the Poisson Problem
!
!     arguments:
!        in:
!             Mdle      - element middle node number
!             NrTest    - total number of test dof
!             NrdofHH   - number of H1 test dof
!             NrdofVV   - number of H(div) test dof
!             NrdofQ    - number of L2 trial dof
!             NrdofH    - number of H1 trial dof
!             NrdofV   - number of H(div) trial dof
!        out:
!             Resid     - element residual (squared)
!             Nref_flag - suggested h-refinement flag
!
!---------------------------------------------------------------------
!
#include "typedefs.h"
!
subroutine elem_residual_poisson_UW(Mdle,                                     &
                                    NrTest,NrdofHH,NrdofVV,                   &
                                    NrdofQ,NrdofH,NrdofV,                     &
                                    Resid,Nref_flag)
!
   use control
   use parametersDPG
   use element_data
   use data_structure3D
!
   implicit none
!..declare input/output variables
   integer, intent(in)  :: Mdle
   integer, intent(in)  :: NrTest
   integer, intent(in)  :: NrdofHH
   integer, intent(in)  :: NrdofVV
   integer, intent(in)  :: NrdofQ
   integer, intent(in)  :: NrdofH
   integer, intent(in)  :: NrdofV
   real(8), intent(out) :: Resid
   integer, intent(out) :: Nref_flag
!
!..declare edge/face type variables
   integer :: etype,ftype
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
!.. L2 shape functions
   real(8) :: shapQ(MAXbrickQ)
!
!..Enriched H1 shape functions
   real(8) :: shapHH(MAXbrickHH), gradHH(3,MAXbrickHH)

!..Enriched H(div) shape functions
   real(8) :: shapVV(3,MAXbrickVV), divVV(MAXbrickVV)
!
!..Gram matrix in packed format
   real(8), allocatable :: gramP(:)
!
!..load vector for the enriched space
   real(8) :: bload_V(NrTest), bload_Vc(NrTest)
!
!..3D quadrature data
   real(8) :: xiloc(3,MAXNINT3ADD), waloc(MAXNINT3ADD)
!
!..2D quadrature data
   real(8) :: tloc(2,MAXNINT2ADD), wtloc(MAXNINT2ADD)
!
!..test functions and gradients
   real(8) :: q, v, aux, divtau_a,divtau_b,tn
   real(8) :: dq(3), dv(3), tau_a(3),tau_b(3)
!
!..approximate solution
   real(8) :: solU
   real(8) :: solVxi(3), solV(3), sigma(3)
   real(8) :: solVn
   real(8) :: solLmb
!
!..various variables for the problem
   real(8) :: rjac, bjac, fval, wa, weight
   integer :: iflag, nrf, nint
   integer :: k1, k2, k, l, i
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
   dofQ = ZERO; dofV = ZERO; dofE = ZERO; dofH = ZERO
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
   call find_orient(Mdle, norient_edge,norient_face)
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
   bload_V = ZERO; gramP = ZERO; bload_Vc = ZERO  
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
!..coordinates and weight of l^{th} integration point
      xi(1:3)=xiloc(1:3,l)
      wa=waloc(l)
!
!..H1 shape functions
      call shape3DH(etype,xi,norder,norient_edge,norient_face, & 
                    nrdof,shapH,gradH)
!
!..L2 shape function
      call shape3DQ(etype,xi,norder, nrdof,shapQ)
!
!..discontinuous H1 shape functions
      call shape3HH(etype,xi,nordP, nrdof,shapHH,gradHH)
!..discontinuous H(div) shape functions
      call shape3VV(etype,xi,nordP, nrdof,shapVV,divVV)
!
!..geometry map
      call geom3D(Mdle,xi,xnod,shapH,gradH,NrdofH, &
                  x,dxdxi,dxidx,rjac,iflag)
!
!..integration weight
      weight = rjac*wa
!
!..compute the approximate solution
      sigma(1:3) = 0.d0
      solU = 0.d0
!
      do k=1,NrdofQ
         solU = solU + dofQ(1,k)*shapQ(k)
      enddo
!..Piola Transformation for L2 variable u
      solU = solU/rjac
!
!..reconstruction of sigma which is a L2 variable  f
      do k = 1,NrdofQ
            sigma(1) = sigma(1) + dofQ(2,k)*shapQ(k)
            sigma(2) = sigma(2) + dofQ(3,k)*shapQ(k)
            sigma(3) = sigma(3) + dofQ(4,k)*shapQ(k)
      enddo
!..Piola transform
      sigma = sigma/rjac
!
!..get the RHS
      call getf(Mdle,x, fval)
!
!..1st loop through enriched H1 test functions
      do k1=1,NrdofHH
!..Piola transformation
         v = shapHH(k1)
         dv(1:3) = gradHH(1,k1)*dxidx(1,1:3) &
                 + gradHH(2,k1)*dxidx(2,1:3) &
                 + gradHH(3,k1)*dxidx(3,1:3)
!
!..residual volume terms: (sigma, grad v) - (f,v)
         aux = (sigma(1)*dv(1) + sigma(2)*dv(2) + sigma(3)*dv(3)) - fval*v
         bload_V(k1) = bload_V(k1) + aux*weight
!
!..2nd loop through enriched H1 test functions
         do k2=k1,NrdofHH
!..Piola transformation
            q = shapHH(k2)
            dq(1:3) = gradHH(1,k2)*dxidx(1,1:3) &
                    + gradHH(2,k2)*dxidx(2,1:3) &
                    + gradHH(3,k2)*dxidx(3,1:3)
!
!..determine index in triangular packed format
            k = nk(k1,k2)
!
!..H1 test inner product: (q,v) + (grad_h q, grad_h v)
            aux = q*v + (dq(1)*dv(1) + dq(2)*dv(2) + dq(3)*dv(3))
            gramP(k) = gramP(k) + aux*weight
!
!     ...end 1st loop through enriched H1 test functions
         enddo
!..cross terms for  graph norm
         do k2 = 1,NrdofVV
            divtau_a = divVV(k2)/rjac
            tau_a(1:3) =     dxdxi(1:3,1) * shapVV(1,k2)      &
                           + dxdxi(1:3,2) * shapVV(2,k2)       &
                           + dxdxi(1:3,3) * shapVV(3,k2)
!            
            tau_a(1:3) = tau_a(1:3)/rjac

            k = nk(k1,NrdofHH+k2)
            aux = dv(1) * tau_a(1) + dv(2) * tau_a(2) + dv(3) * tau_a(3)
            gramP(k) = gramP(k) + aux * weight
         enddo
!..end 2nd loop through enriched H1 test functions
      enddo
!.. starting computations with discontinuous H(div) test functions
      do k1=1,NrdofVV
!..Piola Transform of the divergence
         divtau_a = divVV(k1)/rjac
         tau_a(1:3) =     dxdxi(1:3,1) * shapVV(1,k1)      &
                        + dxdxi(1:3,2) * shapVV(2,k1)       &
                        + dxdxi(1:3,3) * shapVV(3,k1)
!         
         tau_a(1:3) = tau_a(1:3)/rjac
!
         aux = tau_a(1)*sigma(1) + tau_a(2)*sigma(2) + tau_a(3)*sigma(3) + divtau_a * solU
!
         bload_V(NrdofHH + k1) = bload_V(NrdofHH + k1) + aux*weight
!
         do k2 = k1,NrdofVV
            divtau_b = divVV(k2)/rjac
            tau_b(1:3) =     dxdxi(1:3,1) * shapVV(1,k2)       &
                           + dxdxi(1:3,2) * shapVV(2,k2)       &
                           + dxdxi(1:3,3) * shapVV(3,k2)
!   
            tau_b(1:3) = tau_b(1:3)/rjac
!            
            aux = divtau_a * divtau_b + 2.d0* (tau_a(1)*tau_b(1) + tau_a(2)*tau_b(2) + tau_a(3)*tau_b(3))
!
            k = nk(k1 + NrdofHH,k2 + NrdofHH)
            gramP(k) = gramP(k) +  weight * aux
!.. end of first loop through enriched discont H(div) test functions            
         enddo
!.. end of second loop through enriched discont H(div) test functions  
      enddo
!..end loop through integration points
   enddo
!
!-----------------------------------------------------------------------
!.......... boundary integrals
!-----------------------------------------------------------------------
!
!..loop through element faces
   do ifc=1,nrf
!
!..sign factor to determine the OUTWARD normal unit vector
      nsign = nsign_param(etype,ifc)
!
!..face type
      ftype = face_type(etype,ifc)
!
!..face order of approximation
      call face_order(etype,ifc,norder, norderf)
!
!..set 2D quadrature
      INTEGRATION = NORD_ADD
      call set_2D_int_DPG(ftype,norderf,norient_face(ifc), nint,tloc,wtloc)
      INTEGRATION = 0
!
!..loop through integration points
      do l=1,nint
!
!..face coordinates
         t(1:2) = tloc(1:2,l)
!
!...face parametrization
         call face_param(etype,ifc,t, xi,dxidt)
!
!..determine discontinuous H1 shape functions
         call shape3HH(etype,xi,nordP, nrdof,shapHH,gradHH)
!..determine discontinuous H(div) shape functions
         call shape3VV(etype,xi,nordP, nrdof,shapVV,divVV)
!
!..determine element H1 shape functions (for geometry)
         call shape3DH(etype,xi,norder,norient_edge,norient_face, &
                       nrdof,shapH,gradH)
!
!..determine element H(div) shape functions (for fluxes)
         call shape3DV(etype,xi,norder,norient_face, &
                       nrdof,shapV,divV)
!..geometry
         call bgeom3D(Mdle,xi,xnod,shapH,gradH,NrdofH,dxidt,nsign, &
                  x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
         weight = bjac*wtloc(l)
!
!..compute approximate flux at the point
         solVxi(1:3) = 0.d0
         do k=1,NrdofV
            solVxi(1:3) = solVxi(1:3) + dofV(1,k)*shapV(1:3,k)
         enddo
!..Piola transformation
         solV(1:3) = dxdxi(1:3,1)*solVxi(1) &
                   + dxdxi(1:3,2)*solVxi(2) &
                   + dxdxi(1:3,3)*solVxi(3)
         solV(1:3) = solV(1:3) / rjac
!..normal component
         solVn = solV(1)*rn(1) + solV(2)*rn(2) + solV(3)*rn(3)
!
!..loop through enriched test functions
         do k1=1,NrdofHH
            v = shapHH(k1)
!..accumulate for the load vector
            bload_V(k1) = bload_V(k1) - solVn*v*weight
         enddo

!..compute approximate trace at the point
         solLmb = 0.d0
         do k1=1,NrdofH
            solLmb = solLmb + dofH(1,k1) * shapH(k1)
         enddo
!..loop through discont H(div) test functions
         do k1=1,NrdofVV
!..Piola transform for H(div)
            tau_a(1:3) = dxdxi(1:3,1) * shapVV(1,k1)       &
                       + dxdxi(1:3,2) * shapVV(2,k1)       &
                       + dxdxi(1:3,3) * shapVV(3,k1)

            tau_a(1:3) = tau_a(1:3)/rjac
!
            tn = tau_a(1)*rn(1) + tau_a(2)*rn(2) + tau_a(3)*rn(3)
            aux = tn * solLmb
!
            bload_V(NrdofHH + k1) = bload_V(NrdofHH + k1) - weight * aux
!.. end of loop on discont H(div) test functions
         enddo
!.. end of loop over the integeration points
      enddo
!.. end of loop over the faces
   enddo
!
!----------------------------------------------------------------------- 
!..factorize the test stiffness matrix
   call DPPTRF('U', NrTest, gramP, info)
   if (info.ne.0) then
      write(*,*) 'elem_residual_poisson: info = ',info
      stop
   endif
!
!..save copies of the RHS to compute later the residual
   bload_Vc = bload_V
!
!..compute the product of inverted test Gram matrix with RHS,
!..bload_V is overwritten with the solution
   call DPPTRS('U', NrTest, 1, gramP, bload_V, NrTest, info)
   if (info.ne.0) then
      write(*,*) 'elem_residual_poisson: info = ',info
      stop
   endif
!
   deallocate(gramP)
!
!..compute the residual
   Resid = 0.d0
   do k=1,NrTest
      Resid = Resid + bload_Vc(k)*bload_V(k)
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
end subroutine elem_residual_poisson_UW

