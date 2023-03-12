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
   use common_prob_data_UW
#include "syscom.blk"
!
   integer :: ntype,ftype
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
!..approximate solution
   dimension zsolH(MAXEQNH),zdsolH(MAXEQNH,3),zsolV(MAXEQNV,3),   &
             zdivV(MAXEQNV),zsolQ(MAXEQNQ)
!
!..geometry
   dimension xi(3),dxidt(3,2),x(3),dxdxi(3,3),dxidx(3,3),dxdt(3,2),rt(3,2),rn(3),t(2)
!..BC's flags
   dimension ibc(6,NRINDEX)   
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

!..H1 shape functions
   real*8  :: shapH (MAXbrickH)  , gradH(3,MAXbrickH)
   real*8  :: shapHH(MAXbrickHH) , gradHH(3,MAXbrickHH)
!..Hdiv shape functions
   real*8  :: shapV (3,MAXbrickV) , divV(MAXbrickV) 
   real*8  :: shapVV(3,MAXbrickVV), divVV(MAXbrickVV)
!..L2 shape functions
   real*8  :: shapQ(MAXbrickQ)!
!..load vector for the enriched space
   complex*16, allocatable :: BLOADHV(:), AP(:), temp(:)
!
   character*1 uplo
   nk(k1,k2) = (k2-1)*k2/2+k1
!
!----------------------------------------------------------------------------------------   
!
!..element type
   ntype = NODES(Mdle)%ntype
   nrv = nvert(ntype); nre = nedge(ntype); nrf = nface(ntype)
!
!..determine order of approximation
   call find_order(Mdle, norder)
!
!..set the enriched order of approximation
   select case(ntype)
      case(MDLB);      nordP = NODES(Mdle)%order+NORD_ADD*111
      case(MDLN,MDLD); nordP = NODES(Mdle)%order+NORD_ADD
      case(MDLP);      nordP = NODES(Mdle)%order+NORD_ADD*11
   end select

   call compute_enriched_order(nordP, norderP)
!
!..determine edge and face orientations
   call find_orient( Mdle, norient_edge,norient_face)
!                                                                     
!..determine nodes coordinates 
   call nodcor(Mdle, xnod)
!
!..determine element size and scale correction
   call find_hmin(Mdle,h_elem)
!
   omeg = min(OMEGA,6.d0/h_elem)
   ! omeg = OMEGA


   call celndof(NODES(Mdle)%ntype,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
   call celndof(NODES(Mdle)%ntype,norderP,nrdofHH,nrdofEE,nrdofVV,nrdofQQ)
   call ndof_nod(ntype,norder(nre+nrf+1), ndofHmdl,ndofEmdl,ndofVmdl,ndofQmdl)
   nrdofHi = nrdofH - ndofHmdl
   nrdofVi = nrdofV - ndofVmdl

   nrTest  = nrdofHH + nrdofVV
   nrTrial = nrdofHi + nrdofVi + 4*nrdofQ

!..memory for the matrices
   allocate(BloadHV(nrTest))              ; BloadHV   = ZERO
   allocate(AP(nrTest*(nrTest+1)/2))      ; AP        = ZERO
!
!..determine solution dof
   call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
!
!..get the element boundary conditions flags
   call find_bc(Mdle, ibc)
!
!-----------------------------------------------------------------------
!     E L E M E N T    I N T E G R A L S                               |
!-----------------------------------------------------------------------
!
!..use the enriched order to set the quadrature
   INTEGRATION = NORD_ADD
   call set_3Dint_DPG(ntype,norder, nint,xiloc,waloc)
   INTEGRATION = 0
!      
!..loop over integration points      
   do l=1,nint
!      
      xi(1:3)=xiloc(1:3,l) ; wa=waloc(l)
!
!  ...H1 shape functions (for geometry)
      call shape3DH(ntype,xi,norder,norient_edge,norient_face,nrdofH,shapH,gradH)
!
!  ...L2 shape functions for the trial space
      call shape3DQ(ntype,xi,norder, nrdofQ,shapQ)
!
!  ...discontinuous H1 shape functions for the enriched test space
      call shape3HH(ntype,xi,nordP, nrdofHH,shapHH,gradHH)
!
!  ...discontinuous H(div) shape functions for the enriched test space
      call shape3VV(ntype,xi,nordP, nrdofVV,shapVV,divVV)
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
         BloadHV(k1) = BloadHV(k1)                          &
                     + (ZIMG*OMEGA*zp*q                       &
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
         case default
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
                     + ZIMG*omeg*(dq(1)*u(1)+dq(2)*u(2)+dq(3)*u(3) - q*div_u)*weight    
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
         BloadHV(n1) = BloadHV(n1) +                                 &
                     + (ZIMG*OMEGA*(zu(1)*v(1)+zu(2)*v(2)+zu(3)*v(3))  &
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
         case default
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
      nsign = nsign_param(ntype,if)
!
!  ...face type
      ftype = face_type(ntype,if)
!
!  ...face order of approximation
      call face_order(ntype,if,norder, norderf)
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
         call face_param(ntype,if,t, xi,dxidt)
! 
!     ...determine element H1 shape functions (for geometry)
         call shape3DH(ntype,xi,norder,norient_edge,norient_face, nrdofH,shapH,gradH)
! 
!     ...determine element Hdiv shape functions (for fluxes)
         call shape3DV(ntype,xi,norder,norient_face, nrdofV,shapV,divV)
!
!     ...determine discontinuous H1 shape functions
         call shape3HH(ntype,xi,nordP, nrdofHH,shapHH,gradHH)
!         
!     ...determine discontinuous H(div) shape functions
         call shape3VV(ntype,xi,nordP, nrdofVV,shapVV,divVV)
!
!     ...geometry
         call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign,x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
         weight = bjac*wtloc(l)
! 
!     ...compute approximate trace at the point
         zsolH(1) = ZERO
         do k=1,nrdofHi
            zsolH(1) = zsolH(1) + zdofH(1,k)*shapH(k)
         enddo
         zp = zsolH(1)
!
!     ...compute approximate flux at the point
         zvec(1:3) = ZERO
         do k=1,nrdofVi
            zvec(1:3) = zvec(1:3) + zdofV(1,k)*shapV(1:3,k)
         enddo
         zsolV(1,1:3) = (dxdxi(1:3,1)*zvec(1)                  &
                      +  dxdxi(1:3,2)*zvec(2)                  &
                      +  dxdxi(1:3,3)*zvec(3))/rjac            
!         
         zun = ZsolV(1,1)*rn(1)+ZsolV(1,2)*rn(2)+ZsolV(1,3)*rn(3)
!
!     ...impedance BC boundary
         if (ibc(if,2).eq.3) then
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
               BloadHV(k1) = BloadHV(k1) + q*(zp-zg)*weight
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
               BloadHV(k1) = BloadHV(k1) + zun*q*weight
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
            BloadHV(n1) = BloadHV(n1) + zp*vn*weight
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
   call ZTPTRS('U','C','N',nrTEST,1,AP,BloadHV,nrTest,info)
   if (info.ne.0) then
      write(*,*) 'elem: AP ZTPTRS: Mdle,info = ',Mdle,info
      stop 2
   endif
!
   Resid = ZERO
   do i=1,nrTest
      Resid = Resid + conjg(BloadHV(i))*BloadHV(i)
   enddo
!
   Nref_flag = 111
   deallocate(AP,BloadHV)
!
   end subroutine elem_residual

