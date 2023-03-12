!----------------------------------------------------------------------
!                                                                     
!     routine name      - elem
!                                                                     
!---------------------------------------------------------------------- 
!                                                                     
!     latest revision:  - Sep 17
!                                                                     
!     purpose:          - driver for the element routine
!                                                                    
!     arguments:                                                     
!                                                                     
!     in:              
!             Mdle      - an element middle node number, identified
!                         with the element
!     out:              
!             Itest     - index for assembly
!             Itrial    - index for assembly
!
!----------------------------------------------------------------------
!    
   subroutine elem(Mdle, Itest,Itrial)
!
   use physics,       only : NR_PHYSA
   use data_structure3D
   use parametersDPG
   use control,       only : INTEGRATION
!
!----------------------------------------------------------------------
!
   implicit none
!
   integer,                    intent(in)  :: Mdle
   integer,dimension(NR_PHYSA),intent(out) :: Itest,Itrial
   integer :: ntype
!
!..element order, orientation for edges and faces
   integer :: norder(19),norderP(19), nordP
   integer :: nordx,nordy,nordz
   integer :: nrdofH,nrdofE,nrdofV,nrdofQ
   integer :: nrdofHi,nrdofVi
   integer :: nrTest, nrTrial
   integer :: nrdofHH,nrdofEE,nrdofVV,nrdofQQ
   integer :: ndofHmdl,ndofEmdl,ndofVmdl,ndofQmdl
   integer :: nord_add_local
   integer :: nrv, nre, nrf
!
!..sum factorization
   integer :: nrdofHx, nrdofHy, nrdofHz   
   integer :: nrdofHx_tr, nrdofHy_tr, nrdofHz_tr
   integer :: nrdofQx_tr, nrdofQy_tr, nrdofQz_tr   
!----------------------------------------------------------------------
!
   ntype = NODES(Mdle)%ntype
   nrv = nvert(ntype); nre = nedge(ntype); nrf = nface(ntype)
! 
!----------------------------------------------------------------------
!
   Itest (1:NR_PHYSA) = 0
   Itrial(1:NR_PHYSA) = 0

   select case(NODES(Mdle)%case)
   case(7)
      Itest(1:NR_PHYSA)=1; Itrial(1:NR_PHYSA)=1

!  ...initialize nord_add_local
      nord_add_local = NORD_ADD
!      
!  ...determine order of approximation
      call find_order(Mdle, norder)
!
 10   continue
!  ...set the enriched order of approximation
      select case(ntype)
         case(MDLB)      ; nordP = NODES(Mdle)%order + nord_add_local*111
         case(MDLN,MDLD) ; nordP = NODES(Mdle)%order + nord_add_local*1
         case(MDLP)      ; nordP = NODES(Mdle)%order + nord_add_local*11
      end select
!
      call compute_enriched_order(nordP, norderP)
!     
      INTEGRATION = nord_add_local
      call compute_1D_ord(ntype,norder(19),nordx,nordy,nordz)
      INTEGRATION=0
!  ...# dof for 1D H1 test functions with order p+dp
      nrdofHx=nordx+1; nrdofHy=nordy+1; nrdofHz=nordz+1
!      
!  ...# dof for 1D H1 trial functions with order p
      nrdofHx_tr=nrdofHx-nord_add_local
      nrdofHy_tr=nrdofHy-nord_add_local
      nrdofHz_tr=nrdofHz-nord_add_local
!      
!  ...# dof 1D L2 trial functions with order p
      nrdofQx_tr=nrdofHx_tr-1; nrdofQy_tr=nrdofHy_tr-1; nrdofQz_tr=nrdofHz_tr-1

!  ...total number of trial dof associated with the element
      call celndof(ntype,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
!  ...total number of test dof associated with the element
      call celndof(ntype,norderP,nrdofHH,nrdofEE,nrdofVV,nrdofQQ)
!
      call ndof_nod(ntype,norder(nre+nrf+1),ndofHmdl,ndofEmdl,ndofVmdl,ndofQmdl)
      nrTest = nrdofHH+nrdofVV
      nrdofHi = nrdofH - ndofHmdl
      nrdofVi = nrdofV - ndofVmdl
      nrTrial = nrdofHi + nrdofVi + 4*nrdofQ 
!      
      if (nrTest .le. nrTrial-ndofHmdl-ndofVmdl) then
         nord_add_local = nord_add_local + 1
!$omp critical         
         write(*,*) 'elem: WARNING mdle =  ', mdle
         write(*,*) 'elem: nrTest, nrTrial = ', nrTest, nrTrial-ndofHmdl-ndofVmdl
         write(*,*) 'elem: nord_add_local  = ', nord_add_local
!$omp end critical                  
         go to 10
      endif  
!
      call elem_DPG_UWEAK_ACOUSTICS_fi(Mdle,nord_add_local,nrTest,nrTrial,       &
                                       nrdofHH,nrdofVV,nrdofH,nrdofV,nrdofQ,     &
                                       nrdofHi,nrdofVi,                          &
                                       nrdofHx,nrdofHy,nrdofHz,                  &
                                       nrdofQx_tr,nrdofQy_tr,nrdofQz_tr)

!      call elem_DPG_UWEAK_ACOUSTICS(mdle,nord_add_local,nrTest,nrTrial,    &
!                                    nrdofHH,nrdofVV,nrdofH,nrdofV,nrdofQ, &
!                                    ndofHmdl,ndofVmdl)

   case default
      write(*,*) 'elem: Mdle,NODES(Mdle)%case = ',  &
      Mdle,NODES(Mdle)%case
      call logic_error(ERR_INVALID_VALUE, __FILE__,__LINE__)
   end select
!
   end subroutine elem
!   
!----------------------------------------------------------------------
!                                                                     
!     routine name      - elem
!                                                                     
!----------------------------------------------------------------------
!                                                                     
!     latest revision:  - July 17
!                                                                     
!     purpose:          - compute element Stiffness and load
!                                                                    
!     arguments:                                                     
!                                                                     
!     in:              
!             Mdle      - an element middle node number, identified
!                         with the element
!
!----------------------------------------------------------------------
!    
   subroutine elem_DPG_UWEAK_ACOUSTICS(Mdle,nord_add_local,nrTest,nrTrial,    &
                                       nrdofHH,nrdofVV,nrdofH,nrdofV,nrdofQ,  &
                                       ndofHmdl,ndofVmdl)
!
   use data_structure3D
   use parametersDPG
   use control,       only: INTEGRATION
   use assembly,      only: ALOC,BLOC,NR_RHS
   use common_prob_data_UW
!
!   implicit none
!
#include "syscom.blk"
!
!----------------------------------------------------------------------
!
   integer :: ntype,ftype
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
!..3D quadrature data
   dimension xiloc(3,MAXNINT3ADD),waloc(MAXNINT3ADD)
!
!..2D quadrature data
   dimension tloc(2,MAXNINT2ADD),wtloc(MAXNINT2ADD)
!
!..BC's flags
   dimension ibc(6,NRINDEX)
!..workspace for trial and test variables
   dimension dq(3) , u(3), dp(1:3), v(3), vec(3)
!
!..source
   dimension zf(4)
!
!..enriched order
!  
   integer :: nord_add_local
!   
!
!..H1 shape functions
   real*8  :: shapH (MAXbrickH)  , gradH(3,MAXbrickH)
   real*8  :: shapHH(MAXbrickHH) , gradHH(3,MAXbrickHH)
!..Hdiv shape functions
   real*8  :: shapV (3,MAXbrickV) , divV(MAXbrickV) 
   real*8  :: shapVV(3,MAXbrickVV), divVV(MAXbrickVV)
!..L2 shape functions
   real*8  :: shapQ(MAXbrickQ)
!
!..load vector for the enriched space
   complex*16 :: BloadHV(nrTest), Gram(nrTest,nrTest)
!
!..Stiffness matrices for the enriched test space
   complex*16 :: StiffHV_H(nrdofH,nrTest)
   complex*16 :: StiffHV_V(nrdofV,nrTest)
   complex*16 :: StiffHV_Q(4*nrdofQ,nrTest)
   complex*16 :: Stiff_ALL(nrTest,nrTRIAL+1)
   complex*16 :: Zaloc(nrTRIAL+1,nrTRIAL+1)
!
!..lapack
   character*1 uplo
!
!----------------------------------------------------------------------
!
   iprint = 0
!..element type
   ntype = NODES(Mdle)%ntype
   nrv = nvert(ntype); nre = nedge(ntype); nrf = nface(ntype)
!
!..determine order of approximation
   call find_order(Mdle, norder)
!
!..set the enriched order of approximation
   select case(ntype)
      case(MDLB);      nordP = NODES(Mdle)%order + nord_add_local*111
      case(MDLN,MDLD); nordP = NODES(Mdle)%order + nord_add_local*1
      case(MDLP);      nordP = NODES(Mdle)%order + nord_add_local*11
   end select

!..determine edge and face orientations
   call find_orient(Mdle, norient_edge,norient_face)
!                                                                     
!..determine nodes coordinates 
   call nodcor(Mdle, xnod)
!
!..determine element size and scale correction
   call find_hmin(Mdle,h_elem)
! 
!..get the element boundary conditions flags
   call find_bc(Mdle, ibc)
!                                                                    
!..adjusted frequency for the test space 
!    omeg = min(OMEGA,6.d0/h_elem)
   omeg = OMEGA

!
   BloadHV=ZERO; Gram=ZERO; Zaloc = ZERO
   StiffHV_H=ZERO; StiffHV_V = ZERO; StiffHV_Q = ZERO ; Stiff_ALL = ZERO
!
!----------------------------------------------------------------------
!     E L E M E N T    I N T E G R A L S                             
!----------------------------------------------------------------------
!
!..use the enriched order to set the quadrature
   INTEGRATION = nord_add_local
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
         BloadHV(k1) = BloadHV(k1) + q*zf(1)*weight
!
!     ...loop through L2 trial shape functions
         do k2=1,nrdofQ
!          
!        ...Piola transformation
            p = shapQ(k2)/rjac ; u(1:3) = p
!           
            m = (k2-1)*4+1
            StiffHV_Q(m,k1) = StiffHV_Q(m,k1) + ZIMG*OMEGA*q*p*weight
!
            m = (k2-1)*4+2           
            StiffHV_Q(m,k1) = StiffHV_Q(m,k1) - dq(1)*u(1)*weight
!
            m = (k2-1)*4+3
            StiffHV_Q(m,k1) = StiffHV_Q(m,k1) - dq(2)*u(2)*weight
!
            m = (k2-1)*4+4
            StiffHV_Q(m,k1) = StiffHV_Q(m,k1) - dq(3)*u(3)*weight

         enddo
!
!     ...second loop through enriched test functions
         select case(TEST_NORM)
!     ...standard norm
         case(MATHEMATICIANS)    
!        ...H1 test functions
            do k2=1,k1
!           ...Piola transformation
               p = shapHH(k2)
               dp(1:3) = gradHH(1,k2)*dxidx(1,1:3) &
                       + gradHH(2,k2)*dxidx(2,1:3) & 
                       + gradHH(3,k2)*dxidx(3,1:3) 
! 
!           ...determine index in triangular format
               Gram(k2,k1) = Gram(k2,k1)     &
                           + (q*p + dq(1)*dp(1)+dq(2)*dp(2)+dq(3)*dp(3))*weight                     
            enddo         
!                     
         case default      
!        ...H1 test functions
            do k2=1,k1
!           ...Piola transformation
               p = shapHH(k2)
               dp(1:3) = gradHH(1,k2)*dxidx(1,1:3) &
                       + gradHH(2,k2)*dxidx(2,1:3) & 
                       + gradHH(3,k2)*dxidx(3,1:3) 
! 
!           ...accumulate for the gram matrix
               Gram(k2,k1) = Gram(k2,k1)                              &
                           + (dq(1)*dp(1)+dq(2)*dp(2)+dq(3)*dp(3)   &
                           + (omeg**2+ALPHA)*q*p)*weight
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
         BloadHV(n1) = BloadHV(n1) + (v(1)*zf(2)+v(2)*zf(3)+v(3)*zf(4))*weight

!     ...loop through L2 trial shape functions
         do k2=1,nrdofQ
!
!        ...Piola transformation
            p = shapQ(k2)/rjac; u(1:3)=p
!
            m = (k2-1)*4+1
            StiffHV_Q(m,n1) = StiffHV_Q(m,n1) - div_v*p*weight
!            
            m = (k2-1)*4+2
            StiffHV_Q(m,n1) = StiffHV_Q(m,n1) + ZIMG*OMEGA*v(1)*u(1)*weight
!            
            m = (k2-1)*4+3
            StiffHV_Q(m,n1) = StiffHV_Q(m,n1) + ZIMG*OMEGA*v(2)*u(2)*weight
!
            m = (k2-1)*4+4
            StiffHV_Q(m,n1) = StiffHV_Q(m,n1) + ZIMG*OMEGA*v(3)*u(3)*weight
!            
         enddo
!         
!     ...second loop through enriched test functions
         select case(TEST_NORM)
!     ...standard norm
         case(MATHEMATICIANS)    
            do k2 = 1,k1
               n2 = nrdofHH+k2
!           ...Piola transformation
               u(1:3) = (dxdxi(1:3,1)*shapVV(1,k2)             &
                      +  dxdxi(1:3,2)*shapVV(2,k2)             &
                      +  dxdxi(1:3,3)*shapVV(3,k2))/rjac
               div_u  =  divVV(k2)/rjac

               Gram(n2,n1) = Gram(n2,n1)        &
                           + (v(1)*u(1)+v(2)*u(2)+v(3)*u(3) +div_v*div_u)*weight      
            enddo         
!                     
         case default 
!        ...H1 test functions
            do k2=1,nrdofHH
!           ...Piola transformation
               p = shapHH(k2)
               dp(1:3) = gradHH(1,k2)*dxidx(1,1:3) &
                       + gradHH(2,k2)*dxidx(2,1:3) & 
                       + gradHH(3,k2)*dxidx(3,1:3) 
! 
!           ...accumulate for the gram matrix
               Gram(k2,n1) = Gram(k2,n1)                            &
                           + ZIMG*omeg*(dp(1)*v(1)+dp(2)*v(2)+dp(3)*v(3) - p*div_v)*weight
            enddo
            do k2 = 1,k1
               n2 = nrdofHH+k2
!
!           ...Piola transformation
               u(1:3) = (dxdxi(1:3,1)*shapVV(1,k2)             &
                      +  dxdxi(1:3,2)*shapVV(2,k2)             &
                      +  dxdxi(1:3,3)*shapVV(3,k2))/rjac
               div_u  =  divVV(k2)/rjac

               Gram(n2,n1) = Gram(n2,n1)        &
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
   nrdofHi = nrdofH - ndofHmdl
   nrdofVi = nrdofV - ndofVmdl
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
      INTEGRATION = nord_add_local
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
!     ...check if on impedance boundary
         if (ibc(if,2) .eq. 3) then
!        ...get boundary data
            call getg(mdle,x,rn,ibc(if,2),zg)
!
!        ...loop through H1 test functions
            do k1 = 1,nrdofHH
!           ...value of the shape function at the point
               q = shapHH(k1)
!           ...accumulate for the load vector
               BloadHV(k1) = BloadHV(k1) + q*zg*weight
!           ...loop through H1 trial functions
               do k2 = 1,nrdofHi
!              ...value of the shape function at the point
                  p = shapH(k2)
!              ...accumulate for the Stiffness matrix
                  StiffHV_H(k2,k1) = StiffHV_H(k2,k1) + q*p*weight
               enddo
!
            enddo
!
!        ...regular boundary
         else
! 
!        ...loop through enriched H1 test functions
            do k1 = 1,nrdofHH
               q  = shapHH(k1)
!
!           ...loop through H(div) trial functions
               do k2=1,nrdofVi
! 
!              ...normal component (Piola transformation)
                  vec(1:3) = dxdxi(1:3,1)*shapV(1,k2)   &
                           + dxdxi(1:3,2)*shapV(2,k2)   & 
                           + dxdxi(1:3,3)*shapV(3,k2)
                  vec(1:3) = vec(1:3)/rjac
                  un = vec(1)*rn(1)+vec(2)*rn(2)+vec(3)*rn(3)
!
!              ...accumulate for the Stiffness matrix
                  StiffHV_V(k2,k1) = StiffHV_V(k2,k1) + q*un*weight
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
            do k2=1,nrdofHi
!
!           ...value of the shape function at the point
               p = shapH(k2)
!
!           ...accumulate for the rectangular Stiffness matrix
               StiffHV_H(k2,n1) = StiffHV_H(k2,n1) + vn*p*weight
!           ...end of loop through H1 trial functions
            enddo
!        ...end of loop through H(div) test functions                 
         enddo
!     ...end of loop through integration points
      enddo
!  ...end of loop through faces  
   enddo
!
!
!----------------------------------------------------------------------
!      Alternative construction of normal matrix
!----------------------------------------------------------------------
! 
   call celndof(NODES(Mdle)%ntype,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
! 
   nrTRIAL = nrdofHi + nrdofVi + 4*nrdofQ
   i1 = nrTEST ; j1 = nrdofHi ; j2 = nrdofVi ; j3 = 4*nrdofQ
!
   Stiff_ALL(1:i1,1:j1)             = transpose(StiffHV_H(1:j1,1:i1))
   Stiff_ALL(1:i1,j1+1:j1+j2)       = transpose(StiffHV_V(1:j2,1:i1))
   Stiff_ALL(1:i1,j1+j2+1:j1+j2+j3) = transpose(StiffHV_Q(1:j3,1:i1))
   Stiff_ALL(1:i1,j1+j2+j3+1)       = BloadHV(1:i1)
!
   N     = nrTEST
   NRHS  = nrdofHi + nrdofVi + 4*nrdofQ + 1
!
!   call diag_scaling(N,NRHS,Gram,Stiff_ALL)
!
!..diagonal scaling of the gram matrix
   uplo = 'U'
! 
!..factorize the test Stiffness matrix
!
   call ZPOTRF(uplo,N,Gram,N,info)
!
   if (info.ne.0) then
      write(*,*) 'elem: Gram ZPOTRF: Mdle,info = ',Mdle,info
      write(*,*) 'elem: nrTest, nrTrial = ',nrTest, nrTrial
      stop 1
   endif
!
   call ZTRSM('L',uplo,'C','N',N,NRHS,ZONE,Gram,N,Stiff_ALL,N)
! 
   call ZHERK('U','C',NRHS,N,ZONE,Stiff_ALL,N,ZERO,Zaloc,NRHS)
!
   do i=1,NRHS-1
      Zaloc(i+1:NRHS,i) = conjg(Zaloc(i,i+1:NRHS))
   enddo
!
   N = NRHS-1
!
   BLOC(1)%array(1:j1,1) = Zaloc(1:j1,j1+j2+j3+1)
   BLOC(2)%array(1:j2,1) = Zaloc(j1+1:j1+j2,j1+j2+j3+1) 
   BLOC(3)%array(1:j3,1) = Zaloc(j1+j2+1:j1+j2+j3,j1+j2+j3+1) 
!
   ALOC(1,1)%array(1:j1,1:j1) = Zaloc(1:j1,1:j1)
   ALOC(1,2)%array(1:j1,1:j2) = Zaloc(1:j1,j1+1:j1+j2)
   ALOC(1,3)%array(1:j1,1:j3) = Zaloc(1:j1,j1+j2+1:j1+j2+j3)
!
   ALOC(2,1)%array(1:j2,1:j1) = Zaloc(j1+1:j1+j2,1:j1)
   ALOC(2,2)%array(1:j2,1:j2) = Zaloc(j1+1:j1+j2,j1+1:j1+j2)
   ALOC(2,3)%array(1:j2,1:j3) = Zaloc(j1+1:j1+j2,j1+j2+1:j1+j2+j3)
!
   ALOC(3,1)%array(1:j3,1:j1) = Zaloc(j1+j2+1:j1+j2+j3,1:j1)
   ALOC(3,2)%array(1:j3,1:j2) = Zaloc(j1+j2+1:j1+j2+j3,j1+1:j1+j2)
   ALOC(3,3)%array(1:j3,1:j3) = Zaloc(j1+j2+1:j1+j2+j3,j1+j2+1:j1+j2+j3)
!   
   end subroutine elem_DPG_UWEAK_ACOUSTICS

 

 
   subroutine compute_enriched_order(Nord,Norder)

   use parameters , only : MODORDER
!
   implicit none
   integer, intent(in)  :: Nord
   integer :: temp(2)
   integer :: Norder(19),nordF(3)
   integer :: nordB(3),ndofH(3)

   call decod(Nord,MODORDER,2, temp)
   nordF(1) = temp(1) ; nordB(3) = temp(2)
   call decod(nordF(1),MODORDER,2, nordB(1:2))
   call encod((/nordB(1),nordB(3)/),MODORDER,2, nordF(2))
   call encod(nordB(2:3),MODORDER,2, nordF(3))
   norder(1:4)=(/nordB(1),nordB(2),nordB(1),nordB(2)/)
   norder(5:8)=(/nordB(1),nordB(2),nordB(1),nordB(2)/)
   norder(9:12)=nordB(3)
   norder(13:14)=nordF(1)
   norder(15:18)=(/nordF(2),nordF(3),nordF(2),nordF(3)/)
   norder(19)=Nord


   end subroutine compute_enriched_order




   subroutine compute_1D_ord(ntype,nord,nordx,nordy,nordz)

   use parametersDPG,     only : MAXPP
   use control,           only : INTEGRATION
!
   implicit none

   integer, intent(in)  :: ntype
   integer, intent(in)  :: Nord
   integer, intent(out) :: nordx, nordy, nordz
   integer              :: nordh, nord1,nord2,nord3
!
   call decode(Nord, nordh,nord3)
   call decode(nordh, nord1,nord2)


   nordx=min(nord1+INTEGRATION,MAXPP)
   nordy=min(nord2+INTEGRATION,MAXPP)
   nordz=min(nord3+INTEGRATION,MAXPP)

   end subroutine compute_1D_ord
     




   subroutine propagate_flag(Icomp,Nflag)
!
   use data_structure3D
   use common_prob_data_UW, only: IBC_PROB
   use par_mesh,            only: DISTRIBUTED
   use MPI,                 only: MPI_IN_PLACE, MPI_INTEGER, MPI_SUM,   &
                                  MPI_COMM_WORLD
!
   implicit none
!
   integer, intent(in) :: Icomp,Nflag
!
   integer :: ntype
   integer :: iel,mdle,ifc,nrfn,i,j,nod,ierr
!
!..element nodes and orientations, face nodes
   integer :: nodesl(27),norientl(27),nface_nodes(9)
!
!..element face BC flags, decoded BC flag for a node
   integer :: ibc(6,NRINDEX),nodflag(NRINDEX)
!
!-------------------------------------------------------------------------------
!
   select case(IBC_PROB)
   case(3,4,5,6)

   case default
      write(*,*) 'propagate_flag called for non-impedence BC, returning...'
      return
   end select
!
   call reset_visit
!
!..loop through active elements
!$OMP PARALLEL                                     &
!$OMP PRIVATE(ntype,mdle,ifc,nrfn,i,j,nod,nodesl,  &
!$OMP         norientl,nface_nodes,ibc,nodflag)
!$OMP DO
   do iel=1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
      ntype = NODES(mdle)%ntype
!
!  ...determine element nodes
      call elem_nodes(mdle, nodesl,norientl)
!
!  ...get the element boundary conditions flags
      call find_bc(mdle, ibc)
!
!  ...loop through element faces
      do ifc=1,nface(ntype)
!
!     ...if face has a Dirichlet BC flag on this component,
!        then neither propagate Nflag from this face to its edges/vertices,
!        nor prohibit another face from passing Nflag to the edges/vertices.
         if (ibc(ifc,Icomp).eq.1) cycle
!
!     ...determine face node numbers
         call face_nodes(ntype,ifc, nface_nodes,nrfn)
!
!     ...loop through the face nodes
!$OMP CRITICAL
         do i=1,nrfn-1
            j = nface_nodes(i)
            nod = nodesl(j)
!
            NODES(nod)%visit = -Nflag
         enddo
!
         if (ibc(ifc,Icomp).eq.Nflag) then
            nod = nodesl(nface_nodes(nrfn))
            NODES(nod)%visit = Nflag
         endif
!$OMP END CRITICAL
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL
!
   if (DISTRIBUTED) then
!  ...Sum used since values are positive and negative;
!     value isn't important, only sign and whether non-zero
      call MPI_Allreduce(MPI_IN_PLACE,NODES(:)%visit,NRNODS, MPI_INTEGER,MPI_SUM, MPI_COMM_WORLD,ierr)
   endif
!
!..change -Nflag to zero
   do nod=1,NRNODS
      if (NODES(nod)%visit.eq.0) cycle
      call decod(NODES(nod)%bcond,2,NRINDEX, nodflag)
      if (NODES(nod)%visit.le.-Nflag) then
         nodflag(Icomp) = 0
      elseif (NODES(nod)%visit.ge.Nflag) then
         nodflag(Icomp) = 1
      endif
      call encod(nodflag,2,NRINDEX, NODES(nod)%bcond)
!
   enddo
!
   call reset_visit
!
   end subroutine propagate_flag

