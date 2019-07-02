!--------------------------------------------------------------------
!                                                                     
!     routine name      - elem_residual
!                                                                     
!-------------------------------------------------------------------- 
!                                                                     
!     latest revision:  - July 17
!                                                                     
!     purpose:          - routine returns element residual (squared)
!                         for the Ultraweak formulation for acoustics
!                                                                    
!     arguments:                                                     
!                                                                     
!     in:              
!             Mdle      - an element middle node number, identified
!                         with the element
!     out:              
!             Resid     - element residual (squared)
!             Nref_flag - suggested h-refinement flag
!
!---------------------------------------------------------------------
!
   subroutine elem_residual(Mdle, Resid,Nref_flag)

   use control
   use parametersDPG
   use element_data
   use data_structure3D
   use common_prob_data
#include "syscom.blk"
!
   character(len=4) :: etype,ftype
!
!..element order, orientation for edges and faces
   dimension norder(19),norderP(19),norient_edge(12),norient_face(6)
!
!..face order
   dimension norderf(5)
!
!..geometry dof
   dimension xnod(3,MAXbrickH)
!
!..solution dof
   dimension zdofH(MAXEQNH,MAXbrickH), zdofE(MAXEQNE,MAXbrickE),  &
             zdofV(MAXEQNV,MAXbrickV), zdofQ(MAXEQNQ,MAXbrickQ)
!
!..geometry
   dimension xi(3),dxidt(3,2),x(3),dxdxi(3,3),dxidx(3,3),dxdt(3,2),rt(3,2),rn(3),t(2)
!
!..3D quadrature data
   dimension xiloc(3,MAXNINT3ADD),waloc(MAXNINT3ADD)
!
!..2D quadrature data
   dimension tloc(2,MAXNINT2ADD),wtloc(MAXNINT2ADD)
!
!..derivatives wrt physical coordinates, flux
   dimension dv1(3),dv2(3)
!
!..approximate solution
   dimension zgradHxi(3),zgradH(3),zsolVxi(3), zsolV(3)
!
!..H1 shape functions
   real*8  :: shapH (MAXbrickH)  , gradH(3,MAXbrickH)
   real*8  :: shapHH(MAXbrickHH) , gradHH(3,MAXbrickHH)
!..Hdiv shape functions
   real*8  :: shapV (3,MAXbrickV) , divV(MAXbrickV) 
!
!..load vector for the enriched space
   complex*16, allocatable :: BLOADH(:), AP(:)
!   
   character*1 uplo
   nk(k1,k2) = (k2-1)*k2/2+k1
!
!---------------------------------------------------------------------
!
   select case(Mdle)
   case(1)
      iprint=0
   case default
      iprint=0
   end select
!
!..element type
   etype = NODES(Mdle)%type
   nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
!
!..determine order of approximation
   call find_order(Mdle, norder)
!
!..set the enriched order of appoximation
   select case(etype)
   case('mdlb'); nordP = NODES(Mdle)%order+NORD_ADD*111
   case('mdln','mdld'); nordP = NODES(Mdle)%order+NORD_ADD
   case('mdlp'); nordP = NODES(Mdle)%order+NORD_ADD*11
   end select

   norderP(1:nre) = norder(1:nre) + NORD_ADD
   norderP(nre+1:nre+nrf) = norder(nre+1:nre+nrf) + NORD_ADD*11
   norderP(nre+nrf+1) = norder(nre+nrf+1) + NORD_ADD*111   
!
!..determine edge and face orientations
   call find_orient( Mdle, norient_edge,norient_face)
!                                                                     
!..determine nodes coordinates 
   call nodcor(Mdle, xnod)
!
!..determine element size and scale correction
   h_elem = min(abs(xnod(1,2)-xnod(1,1)),abs(xnod(2,4)-xnod(2,1)),   &
                abs(xnod(3,5)-xnod(3,1)))
!
!

   call celndof(NODES(Mdle)%type,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
   call celndof(NODES(Mdle)%type,norderP,nrdofHH,nrdofEE,nrdofVV,nrdofQQ)

   nrTEST  = nrdofHH
   nrTRIAL = nrdofH + nrdofV

!..memory for the matrices
   allocate(BLOADH(nrTEST))               ; BLOADH   = ZERO
   allocate(AP(nrTEST*(nrTEST+1)/2))      ; AP       = ZERO

!..determine solution dof
   call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
!
!-----------------------------------------------------------------------
!
!..element integrals...
!
!..use the enriched order to set the quadrature
   INTEGRATION = NORD_ADD
   call set_3Dint_DPG(etype,norder, nint3,xiloc,waloc)
   INTEGRATION = 0
   do l=1,nint3
      xi(1:3) = xiloc(1:3,l)
      wa = waloc(l)
!
!  ...determine element H1 shape functions
      call shape3H(etype,xi,norder,norient_edge,norient_face, nrdofH,shapH,gradH)
!
!  ...determine discontinuous H1 shape functions
      call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
!
!  ...geometry
      call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,x,dxdxi,dxidx,rjac,iflag)
!
!  ...integration weight 
      weight = rjac*wa
!
!  ...compute the approximate solution
      zsolH = 0.d0; zgradHxi(1:3) = 0.d0
      do k=1,nrdofH
         zsolH = zsolH + zdofH(1,k)*shapH(k)
         zgradHxi(1:3) = zgradHxi(1:3) + zdofH(1,k)*gradH(1:3,k)
      enddo
      zgradH(1:3) = zgradHxi(1)*dxidx(1,1:3)   &
                  + zgradHxi(2)*dxidx(2,1:3)   &
                  + zgradHxi(3)*dxidx(3,1:3) 
!
!  ...get the RHS
      call getf(Mdle,x, zfval)
!
!  ...loop through enriched H1 test functions
      do k1=1,nrdofHH
         v1 = shapHH(k1)
         dv1(1:3) = gradHH(1,k1)*dxidx(1,1:3)    &
                  + gradHH(2,k1)*dxidx(2,1:3)    &
                  + gradHH(3,k1)*dxidx(3,1:3)
!
!     ...compute the RHS
         BLOADH(k1) = BLOADH(k1)                                              &
                    + (zgradH(1)*dv1(1)+zgradH(2)*dv1(2)+zgradH(3)*dv1(3)     &
                    - OMEGA**2*zsolH*v1                                       &            
                    - zfval*v1)*weight
!
!     ...loop through enriched H1 trial functions
         do k2=k1,nrdofHH
            v2 = shapHH(k2)
            dv2(1:3) = gradHH(1,k2)*dxidx(1,1:3)      &
                     + gradHH(2,k2)*dxidx(2,1:3)      &
                     + gradHH(3,k2)*dxidx(3,1:3)
!
!        ...accumulate for the test stiffness matrix
            k = nk(k1,k2)
            AP(k) = AP(k)                                         &
                  + (dv1(1)*dv2(1)+dv1(2)*dv2(2)+dv1(3)*dv2(3)    &
                  +v1*v2)*weight
         enddo
      enddo
   enddo
!
!-----------------------------------------------------------------------
!
!..boundary integrals
!
!..loop through element faces
   do if=1,nrf
!
!  ...sign factor to determine the OUTWARD normal unit vector
      nsign = nsign_param(etype,if)
!
!  ...face type
      ftype = face_type(etype,if)
!
!  ...face order of approximation
      call face_order(etype,if,norder, norderf)
!
!  ...set 2D quadrature
      INTEGRATION = NORD_ADD
      call set_2Dint_DPG(ftype,norderf, nint,tloc,wtloc)
      INTEGRATION = 0
!
!  ...loop through integration points
      do l=1,nint
!
!     ...face coordinates
         t(1:2) = tloc(1:2,l)
!
!     ...face parametrization
         call face_param(etype,if,t, xi,dxidt)
!
!     ...determine discontinuous H1 shape functions
         call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
!
!     ...determine element H1 shape functions (for geometry)
         call shape3H(etype,xi,norder,norient_edge,norient_face, nrdofH,shapH,gradH)
!
!     ...determine element Hdiv shape functions (for fluxes)
!
         call shape3V(etype,xi,norder,norient_face, nrdofV,shapV,divV)
!
!     ...geometry
         call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign,     &
                      x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
         weight = bjac*wtloc(l)
!
!     ...compute approximate flux at the point
         zsolVxi(1:3) = 0.d0
         do k=1,nrdofV
            zsolVxi(1:3) = zsolVxi(1:3) + zdofV(1,k)*shapV(1:3,k)
         enddo
         zsolV(1:3) = (dxdxi(1:3,1)*zsolVxi(1)                 &
                    + dxdxi(1:3,2)*zsolVxi(2)                  &
                    + dxdxi(1:3,3)*zsolVxi(3))/rjac            
!         
         zsolVn = ZsolV(1)*rn(1)+ZsolV(2)*rn(2)+ZsolV(3)*rn(3)
!     ...loop through enriched test functions
         do k1=1,nrdofHH
            v1 = shapHH(k1)
!
!        ...accumulate for the load vector
            BLOADH(k1) = BLOADH(k1) - zsolVn*v1*weight
         enddo
      enddo
   enddo
!
!-----------------------------------------------------------------------
!
!..factorize the test stiffness matrix
   uplo = 'U'

   call ZPPTRF(uplo, nrTEST, AP, info) 
   if (info.ne.0) then
      write(*,*) 'elem_residual: AP info = ',info
      stop 1
   endif
!
   call ZTPTRS('U','C','N',nrTEST,1,AP,BLOADH,nrTEST,info1)
   if (info.ne.0) then
      write(*,*) 'elem: AP ZTPTRS: Mdle,info = ',Mdle,info
      stop 2
   endif
! 
   Resid = zdotc(nrTEST,BLOADH,1,BLOADH,1)
!
   Nref_flag = 111

   deallocate(AP,BLOADH)

!
   end subroutine elem_residual

