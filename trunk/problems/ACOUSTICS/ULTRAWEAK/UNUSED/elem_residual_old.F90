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

   use data_structure3D
   use control
   use parametersDPG
   use uweak_acoustics_module
   use common_prob_data
#include "syscom.blk"
!
  character(len=4) :: etype,ftype
!
!..element order, orientation for edges and faces
   dimension norder(19),norient_edge(12),norient_face(6)
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
!..approximate solution
   dimension zsolH(MAXEQNH),zdsolH(MAXEQNH,3),zsolV(MAXEQNV,3),   &
             zdivV(MAXEQNV),zsolQ(MAXEQNQ)
!
!..geometry
   dimension xi(3),dxidt(3,2),x(3),dxdxi(3,3),dxidx(3,3),dxdt(3,2),rt(3,2),rn(3),t(2)
!..BC's flags
   dimension ibc(6,NR_PHYSA)   
!
!..H1 shape functions
   dimension shapH(MAXbrickH),gradH(3,MAXbrickH)
!
!..Hdiv shape functions
   dimension shapV(3,MAXbrickV),divV(MAXbrickV)
!
!..L2 shape functions 
   dimension shapQ(MAXbrickQ)
!   
!..3D quadrature data
   dimension xiloc(3,MAXNINT3ADD),waloc(MAXNINT3ADD)
!
!..2D quadrature data
   dimension tloc(2,MAXNINT2ADD),wtloc(MAXNINT2ADD)
!
!..workspace for trial and test variables
   dimension dq(3) , u(3), dp(1:3), v(3), vec(3), zu(3), zvec(3)
!
!..source
   dimension zf(4)
!
!..imaginary unit
   complex*16, parameter :: zi = (0.0d0,1.0d0)   
!
   character*1 uplo
   nk(k1,k2) = (k2-1)*k2/2+k1
!
!----------------------------------------------------------------------------------------   
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
  omeg = min(OMEGA,6.d0/h_elem)
   ! omeg = OMEGA

!
!..determine solution dof
   call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
!
!..get the element boundary conditions flags
   call find_bc(Mdle, ibc)

   if (iprint.eq.1) then
      write(*,7020) xnod(1,1:8),xnod(2,1:8),xnod(3,1:8)
 7020   format('elem_residual: xnod  = ',8(f8.3,2x),2(  /,'                ',8(f8.3,2x)))
        write(*,7030) zdofH(1,1:8),zdofV(1,1:6)
 7030   format('elem_residual: zdofH = ',8(e12.5,2x),   /,'        zdofV = ',6(e12.5,2x))
   endif
!
!..clear space for auxiliary matrices
   BLOADHV = ZERO; AP = ZERO
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
!     E L E M E N T    I N T E G R A L S                               |
!-----------------------------------------------------------------------
!
!..use the enriched order to set the quadrature
   INTEGRATION = NORD_ADD
   call set_3Dint_DPG(etype,norder, nint,xiloc,waloc)
   INTEGRATION = 0
!      
!..loop over integration points      
   do l=1,nint
!      
      xi(1:3)=xiloc(1:3,l) ; wa=waloc(l)
!
!  ...H1 shape functions (for geometry)
      call shape3H(etype,xi,norder,norient_edge,norient_face,nrdofH,shapH,gradH)
!
!  ...L2 shape functions for the trial space
      call shape3Q(etype,xi,norder, nrdofQ,shapQ)
!
!  ...discontinuous H1 shape functions for the enriched test space
      call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
!
!  ...discontinuous H(div) shape functions for the enriched test space
      call shape3VV(etype,xi,nordP, nrdofVV,shapVV,divVV)      
!
!  ...geometry map
      call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,x,dxdxi,dxidx,rjac,iflag)
!
!  ...integration weight 
      weight = rjac*wa
!
!  ...get the RHS
      call getf(Mdle,x, zf)
!
!  ...evaluate the L2 solution
      zsolQ(1:4) = ZERO
      do k=1,NrdofQ
         zsolQ(1:4) = zsolQ(1:4) + zdofQ(1:4,k)*shapQ(k)
      enddo
!  ...Piola transformation
      zp = zsolQ(1)/rjac
      zu(1:3) = zsolQ(2:4)/rjac
!      
!  ...loop through enriched H1 test functions in the enriched space
      do k1=1,nrdofHH
!
!     ...Piola transformation
         q = shapHH(k1)
         dq(1:3) = gradHH(1,k1)*dxidx(1,1:3) &
                 + gradHH(2,k1)*dxidx(2,1:3) & 
                 + gradHH(3,k1)*dxidx(3,1:3) 
!
!     ...accumulate for the residual    
!
         BLOADHV(k1) = BLOADHV(k1)                          &
                     + (zi*OMEGA*zp*q                       &
                     -  zu(1)*dq(1)-zu(2)*dq(2)-zu(3)*dq(3) &
                     -  zf(1)*q)                            &
                     *  weight
!          
!     ...second loop through enriched test functions
         select case(TEST_NORM)
!     ...standard norm
         case(MATHEMATICIANS)    
!        ...H1 test functions
            do k2=k1,nrdofHH
!           ...Piola transformation
               p = shapHH(k2)
               dp(1:3) = gradHH(1,k2)*dxidx(1,1:3) &
                       + gradHH(2,k2)*dxidx(2,1:3) & 
                       + gradHH(3,k2)*dxidx(3,1:3) 
! 
!           ...determine index in triangular format
               k = nk(k1,k2)
!           ...accumulate for the gram matrix
!  
               AP(k) = AP(k)     &
                     + (q*p + dq(1)*dp(1)+dq(2)*dp(2)+dq(3)*dp(3))*weight
            enddo         
!                     
         case(ADJOINT_GRAPH)
!        ...H1 test functions
            do k2=k1,nrdofHH
!           ...Piola transformation
               p = shapHH(k2)
               dp(1:3) = gradHH(1,k2)*dxidx(1,1:3) &
                       + gradHH(2,k2)*dxidx(2,1:3) & 
                       + gradHH(3,k2)*dxidx(3,1:3) 

!           ...determine index in triangular format
               k = nk(k1,k2)
!           ...accumulate for the gram matrix
!  
               AP(k) = AP(k)                                     &
                     + (dq(1)*dp(1)+dq(2)*dp(2)+dq(3)*dp(3)      &
                     + (omeg**2+ALPHA)*q*p)*weight
            enddo
!        ...H(div) test functions
            do k2=1,nrdofVV
               n2 = nrdofHH + k2
!           ...Piola transformation
               u(1:3) = (dxdxi(1:3,1)*shapVV(1,k2)                 &
                      +  dxdxi(1:3,2)*shapVV(2,k2)                 &
                      +  dxdxi(1:3,3)*shapVV(3,k2))/rjac
               div_u  =  divVV(k2)/rjac
!
               k = nk(k1,n2)
               AP(k) = AP(k)  &
                     + ZI*omeg*(dq(1)*u(1)+dq(2)*u(2)+dq(3)*u(3) - q*div_u)*weight    
            enddo            
         end select   
!     ...end of loop through H1 test functions         
      enddo
!      
!  ...loop through enriched H(div) test functions in the enriched space
      do k1 = 1,nrdofVV
         n1 = nrdofHH+k1
!     ...Piola transformation
         v(1:3) = (dxdxi(1:3,1)*shapVV(1,k1)                 &
                +  dxdxi(1:3,2)*shapVV(2,k1)                 &
                +  dxdxi(1:3,3)*shapVV(3,k1))/rjac
         div_v  =  divVV(k1)/rjac
!         
!     ...accumulate for the load vector        
!
         BLOADHV(n1) = BLOADHV(n1) +                                 &
                     + (zi*OMEGA*(zu(1)*v(1)+zu(2)*v(2)+zu(3)*v(3))  &
                     -  zp*div_v                                     &
                     -  zf(2)*v(1)- zf(3)*v(2) - zf(4)*v(3))         &
                     * weight
!         
!     ...second loop through enriched test functions
         select case(TEST_NORM)
!     ...standard norm
         case(MATHEMATICIANS)    
            do k2 = k1,nrdofVV
               n2 = nrdofHH+k2
!           ...Piola transformation
               u(1:3) = (dxdxi(1:3,1)*shapVV(1,k2)             &
                      +  dxdxi(1:3,2)*shapVV(2,k2)             &
                      +  dxdxi(1:3,3)*shapVV(3,k2))/rjac
               div_u  =  divVV(k2)/rjac

               k = nk(n1,n2)
               AP(k) = AP(k)        &
                     + (v(1)*u(1)+v(2)*u(2)+v(3)*u(3) +div_v*div_u)*weight
            enddo         
!                     
         case(ADJOINT_GRAPH) 
            do k2 = k1,nrdofVV
               n2 = nrdofHH+k2
!
!           ...Piola transformation
               u(1:3) = (dxdxi(1:3,1)*shapVV(1,k2)             &
                      +  dxdxi(1:3,2)*shapVV(2,k2)             &
                      +  dxdxi(1:3,3)*shapVV(3,k2))/rjac
               div_u  =  divVV(k2)/rjac

               k = nk(n1,n2)
               AP(k) = AP(k)        &
                     + ((omeg**2+ALPHA)*(v(1)*u(1)+v(2)*u(2)+v(3)*u(3))  &
                     + (div_v*div_u))*weight
            enddo
         end select   
!     ...end of loop through H(div) test functions      
      enddo
!  ...end of loop through integration points      
   enddo
!
!-----------------------------------------------------------------------
!     B O U N D A R Y    I N T E G R A L S                             |
!-----------------------------------------------------------------------
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
!     ...determine element H1 shape functions (for geometry)
         call shape3H(etype,xi,norder,norient_edge,norient_face, nrdofH,shapH,gradH)
! 
!     ...determine element Hdiv shape functions (for fluxes)
         call shape3V(etype,xi,norder,norient_face, nrdofV,shapV,divV)
!
!     ...determine discontinuous H1 shape functions
         call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
!         
!     ...determine discontinuous H(div) shape functions
         call shape3VV(etype,xi,nordP, nrdofVV,shapVV,divVV)
!
!     ...geometry
         call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign,x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
         weight = bjac*wtloc(l)
! 
!     ...compute approximate trace at the point
         zsolH(1) = ZERO
         do k=1,nrdofH
            zsolH(1) = zsolH(1) + zdofH(1,k)*shapH(k)
         enddo
         zp = zsolH(1)
!
!     ...compute approximate flux at the point
         zvec(1:3) = ZERO
         do k=1,nrdofV
            zvec(1:3) = zvec(1:3) + zdofV(1,k)*shapV(1:3,k)
         enddo
         zsolV(1,1:3) = (dxdxi(1:3,1)*zvec(1)                  &
                      +  dxdxi(1:3,2)*zvec(2)                  &
                      +  dxdxi(1:3,3)*zvec(3))/rjac            
!         
         zun = ZsolV(1,1)*rn(1)+ZsolV(1,2)*rn(2)+ZsolV(1,3)*rn(3)
!
!     ...impedance BC boundary
         if (ibc(if,2).eq.9) then
!
!        ...get the boundary source
            call getg(Mdle,x,rn,ibc(if,2), zg)
!
!        ...loop through H1 enriched test functions 
            do k1=1,nrdofHH
!
!           ...value of the shape function at the point
               q = shapHH(k1)
!
!           ...accumulate for the load vector
               BLOADHV(k1) = BLOADHV(k1) + q*(zp-zg)*weight
            enddo
!
!        ...regular boundary
         else
!        ...loop through enriched H1 test functions
            do k1 = 1,nrdofHH
!               
!           ...value of the shape function at the point
               q  = shapHH(k1)
!               
!           ...accumulate for the load vector
               BLOADHV(k1) = BLOADHV(k1) + zun*q*weight
!
!           ...end of loop through H1 test functions
            enddo
         endif   
!         
!        ...loop through H(div) enriched test functions
         do k1 = 1,nrdofVV
            n1 = nrdofHH+k1
!
!        ...normal component of the test function (Piola transformation at work!)
            vec(1:3) = dxdxi(1:3,1)*shapVV(1,k1)   &
                     + dxdxi(1:3,2)*shapVV(2,k1)   & 
                     + dxdxi(1:3,3)*shapVV(3,k1)
            vec(1:3) = vec(1:3)/rjac
            vn = vec(1)*rn(1)+vec(2)*rn(2)+vec(3)*rn(3)
!
            BLOADHV(n1) = BLOADHV(n1) + zp*vn*weight
!        ...end of loop through H(div) test functions                 
         enddo
!     ...end of loop through integration points
      enddo
!  ...end of loop through faces  
   enddo
!
!-----------------------------------------------------------------------
!
!..factorize the test stiffness matrix
   uplo = 'U'
   nrTEST = nrdofHH+nrdofVV

   call ZPPTRF(uplo, nrTEST, AP, info) 
   if (info.ne.0) then
      write(*,*) 'elem_residual: AP info = ',info
      stop 1
   endif
!
!
   call ZTPTRS('U','C','N',nrTEST,1,AP,BLOADHV,MAXtest,info)
   if (info.ne.0) then
      write(*,*) 'elem: AP ZTPTRS: Mdle,info = ',Mdle,info
      stop 2
   endif
! 
   Resid = zdotc(nrTEST,BLOADHV,1,BLOADHV,1)
!
   Nref_flag = 111

!
   end subroutine elem_residual

