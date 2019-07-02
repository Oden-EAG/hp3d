
!----------------------------------------------------------------------
!                                                                     
!     routine name      - elem
!                                                                     
!----------------------------------------------------------------------
!                                                                     
!     latest revision:  - July 17
!                                                                     
!     purpose:          - compute element stiffness and load
!                                                                    
!     arguments:                                                     
!                                                                     
!     in:              
!             Mdle      - an element middle node number, identified
!                         with the element
!
!----------------------------------------------------------------------
!    
   subroutine elem_DPG_UWEAK_ACOUSTICS_old(Mdle)
!
   use data_structure3D
   use control,     only: INTEGRATION
   use parametersDPG
   use physics,     only: NR_PHYSA
   use assembly,    only: ALOC,BLOC,NR_RHS
   use common_prob_data
   use uweak_acoustics_module
!
#include "syscom.blk"
!
!----------------------------------------------------------------------
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
!..geometry
   dimension xi(3),dxidt(3,2),x(3),dxdxi(3,3),dxidx(3,3), dxdt(3,2)
   dimension rt(3,2),rn(3),t(2)
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
!..BC's flags
   dimension ibc(6,NR_PHYSA)
!..workspace for trial and test variables
   dimension dq(3) , u(3), dp(1:3), v(3), vec(3)
!
!..source
   dimension zf(4)
!
!..imaginary unit
   complex*16, parameter :: zi = (0.0d0,1.0d0)  
!   
!..lapack
   character*1 uplo
   nk(k1,k2) = (k2-1)*k2/2+k1
!
!----------------------------------------------------------------------
!
!..element type
   etype = NODES(Mdle)%type
   nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
!
!..determine order of approximation
   call find_order(Mdle, norder)
!
!..set the enriched order of approximation
   select case(etype)
      case('mdlb')        ; nordP = NODES(Mdle)%order + NORD_ADD*111
      case('mdln','mdld') ; nordP = NODES(Mdle)%order + NORD_ADD*1
      case('mdlp')        ; nordP = NODES(Mdle)%order + NORD_ADD*11
   end select
!
!..determine edge and face orientations
   call find_orient(Mdle, norient_edge,norient_face)
!                                                                     
!..determine nodes coordinates 
   call nodcor(Mdle, xnod)
!
!..determine element size and scale correction
   h_elem = min(abs(xnod(1,2)-xnod(1,1)),abs(xnod(2,4)-xnod(2,1)),  &
                abs(xnod(3,5)-xnod(3,1)))
! 
!..get the element boundary conditions flags
   call find_bc(Mdle, ibc)
!                                                                    
!..extended load vector
   BLOADHV=ZERO 
!..extended stiffness matrices
   STIFFHV_H = ZERO; STIFFHV_V = ZERO; STIFFHV_Q = ZERO
!..Gram matrix
   AP=ZERO
!..adjusted frequency for the test space 
   omeg = min(OMEGA,6.d0/h_elem)
   ! omeg = OMEGA
!
!----------------------------------------------------------------------
!     E L E M E N T    I N T E G R A L S                             
!----------------------------------------------------------------------
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
!  ...loop through enriched H1 test functions in the enriched space
      do k1=1,nrdofHH
!
!     ...Piola transformation
         q = shapHH(k1)
         dq(1:3) = gradHH(1,k1)*dxidx(1,1:3) &
                 + gradHH(2,k1)*dxidx(2,1:3) & 
                 + gradHH(3,k1)*dxidx(3,1:3) 
!
!     ...accumulate for the load vector        
!
         BLOADHV(k1) = BLOADHV(k1) + q*zf(1)*weight
!          
!     ...loop through L2 trial shape functions
         do k2=1,nrdofQ
!          
!        ...Piola transformation
            p = shapQ(k2)/rjac ; u(1:3) = p
!           
            m = (k2-1)*4+1
            STIFFHV_Q(k1,m) = STIFFHV_Q(k1,m) + ZI*OMEGA*q*p*weight
!
            m = (k2-1)*4+2           
            STIFFHV_Q(k1,m) = STIFFHV_Q(k1,m) - dq(1)*u(1)*weight         
!
            m = (k2-1)*4+3
            STIFFHV_Q(k1,m) = STIFFHV_Q(k1,m) - dq(2)*u(2)*weight
!
            m = (k2-1)*4+4
            STIFFHV_Q(k1,m) = STIFFHV_Q(k1,m) - dq(3)*u(3)*weight
         enddo   
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
! 
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
         BLOADHV(n1) = BLOADHV(n1) + (v(1)*zf(2)+v(2)*zf(3)+v(3)*zf(4))*weight

!     ...loop through L2 trial shape functions
         do k2=1,nrdofQ
!
!        ...Piola transformation
            p = shapQ(k2)/rjac; u(1:3)=p
!
            m = (k2-1)*4+1
            STIFFHV_Q(n1,m) = STIFFHV_Q(n1,m) - div_v*p*weight
!            
            m = (k2-1)*4+2
            STIFFHV_Q(n1,m) = STIFFHV_Q(n1,m) + ZI*OMEGA*v(1)*u(1)*weight
!            
            m = (k2-1)*4+3
            STIFFHV_Q(n1,m) = STIFFHV_Q(n1,m) + ZI*OMEGA*v(2)*u(2)*weight
!
            m = (k2-1)*4+4
            STIFFHV_Q(n1,m) = STIFFHV_Q(n1,m) + ZI*OMEGA*v(3)*u(3)*weight
!            
         enddo
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
!----------------------------------------------------------------------
!     B O U N D A R Y    I N T E G R A L S                           
!----------------------------------------------------------------------
!
   call celndof(NODES(Mdle)%type,norder, nrdofH,nrdofE,nrdofV,nrdofQ)

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
!     ...check if on impedance boundary
         if (ibc(if,2) .eq. 9) then
!        ...get boundary data
            call getg(mdle,x,rn,ibc(if,2),zg)   
!
!        ...loop through H1 test functions
            do k1 = 1,nrdofHH
!           ...value of the shape function at the point
               q = shapHH(k1)
!           ...accumulate for the load vector
               BLOADHV(k1) = BLOADHV(k1) + q*zg*weight
!           ...loop through H1 trial functions
               do k2 = 1,nrdofH
!              ...value of the shape function at the point
                  p = shapH(k2)
!              ...accumulate for the stiffness matrix
                  STIFFHV_H(k1,k2) = STIFFHV_H(k1,k2) + q*p*weight
               enddo
            enddo                                                
!        ...regular boundary
         else
! 
!        ...loop through enriched H1 test functions
            do k1 = 1,nrdofHH
               q  = shapHH(k1)
!
!           ...loop through H(div) trial functions
               do k2=1,nrdofV
! 
!              ...normal component (Piola transformation)
                  vec(1:3) = dxdxi(1:3,1)*shapV(1,k2)   &
                           + dxdxi(1:3,2)*shapV(2,k2)   & 
                           + dxdxi(1:3,3)*shapV(3,k2)
                  vec(1:3) = vec(1:3)/rjac
                  un = vec(1)*rn(1)+vec(2)*rn(2)+vec(3)*rn(3)
!
!              ...accumulate for the stiffness matrix
!
                  STIFFHV_V(k1,k2) = STIFFHV_V(k1,k2) + q*un*weight
!
!              ...end of loop through H(div) trial functions
               enddo
!           ...end of loop through H1 test functions
            enddo
         endif   
!         
!     ...loop through H(div) enriched test functions
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
!        ...loop through H1 trial functions
            do k2=1,nrdofH
!
!           ...value of the shape function at the point
               p = shapH(k2)
!
!           ...accumulate for the rectangular stiffness matrix
               STIFFHV_H(n1,k2) = STIFFHV_H(n1,k2) + vn*p*weight
!           ...end of loop through H1 trial functions
            enddo
!        ...end of loop through H(div) test functions                 
         enddo
!     ...end of loop through integration points
      enddo
!  ...end of loop through faces  
   enddo
!
!----------------------------------------------------------------------
!      Alternative construction of normal matrix
!----------------------------------------------------------------------
! 
   call celndof(NODES(Mdle)%type,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
! 
   i1 = MAXtest ; j1 = nrdofH ; j2 = nrdofV ; j3 = 4*nrdofQ

   STIFF_ALL(1:i1,1:j1)             = STIFFHV_H(1:i1,1:j1)
   STIFF_ALL(1:i1,j1+1:j1+j2)       = STIFFHV_V(1:i1,1:j2)
   STIFF_ALL(1:i1,j1+j2+1:j1+j2+j3) = STIFFHV_Q(1:i1,1:j3)
   STIFF_ALL(1:i1,j1+j2+j3+1)       = BLOADHV(1:i1)
!
   nrTEST = nrdofHH+nrdofVV
   N     = nrTEST
   NRHS  = nrdofH + nrdofV + 4*nrdofQ + 1

!..diagonal scaling of the gram matrix
   ! call diag_scaling(N,NRHS,AP(1:N*(N+1)/2), STIFF_ALL(1:N,1:NRHS))
   uplo = 'U'
! 
!..factorize the test stiffness matrix
   call ZPPTRF(uplo, nrTEST, AP(1:N*(N+1)/2), info)
   if (info.ne.0) then
      write(*,*) 'elem: AP ZPPTRF: Mdle,info = ',Mdle,info
      stop 1
   endif
! 
   call ZTPTRS('U','C','N',N,NRHS,AP,STIFF_ALL,MAXtest,info)
   if (info.ne.0) then
      write(*,*) 'elem: AP ZTPTRS: Mdle,info = ',Mdle,info
      stop 2
   endif
!
   call ZHERK('U','C',NRHS,N,1.0d0,STIFF_ALL(1:N,1:NRHS),N,0.0d0,   &
              STIFF_ALL(1:NRHS,1:NRHS),NRHS)
! 
   do i=1,NRHS-1
      STIFF_ALL(i+1:NRHS,i) = conjg(STIFF_ALL(i,i+1:NRHS))
   enddo

! 
   BLOC(1)%array(1:j1,1) = STIFF_ALL(1:j1,j1+j2+j3+1)
   BLOC(2)%array(1:j2,1) = STIFF_ALL(j1+1:j1+j2,j1+j2+j3+1) 
   BLOC(3)%array(1:j3,1) = STIFF_ALL(j1+j2+1:j1+j2+j3,j1+j2+j3+1) 
!
   ALOC(1,1)%array(1:j1,1:j1) = STIFF_ALL(1:j1,1:j1)
   ALOC(1,2)%array(1:j1,1:j2) = STIFF_ALL(1:j1,j1+1:j1+j2)
   ALOC(1,3)%array(1:j1,1:j3) = STIFF_ALL(1:j1,j1+j2+1:j1+j2+j3)
!
   ALOC(2,1)%array(1:j2,1:j1) = STIFF_ALL(j1+1:j1+j2,1:j1)
   ALOC(2,2)%array(1:j2,1:j2) = STIFF_ALL(j1+1:j1+j2,j1+1:j1+j2)
   ALOC(2,3)%array(1:j2,1:j3) = STIFF_ALL(j1+1:j1+j2,j1+j2+1:j1+j2+j3)
!
   ALOC(3,1)%array(1:j3,1:j1) = STIFF_ALL(j1+j2+1:j1+j2+j3,1:j1)
   ALOC(3,2)%array(1:j3,1:j2) = STIFF_ALL(j1+j2+1:j1+j2+j3,j1+1:j1+j2)
   ALOC(3,3)%array(1:j3,1:j3) = STIFF_ALL(j1+j2+1:j1+j2+j3,j1+j2+1:j1+j2+j3)
!
! 
!   
   end subroutine elem_DPG_UWEAK_ACOUSTICS_old