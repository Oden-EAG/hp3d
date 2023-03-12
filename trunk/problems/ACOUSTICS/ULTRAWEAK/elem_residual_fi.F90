!    
!----------------------------------------------------------------------
!                                                                     
!     routine name      - elem_DPG_UWEAK_ACOUSTICS_fi
!                                                                     
!----------------------------------------------------------------------
!                                                                     
!     latest revision:  - July 17
!                                                                     
!     purpose:          - compute element stiffness and load
!                         using fast quadrature for hexahedra
!                                                                    
!     arguments:                                                     
!                                                                     
!     in:              
!             Mdle      - an element middle node number, identified
!                         with the element
!
!----------------------------------------------------------------------
!
#include "typedefs.h"
!
   subroutine elem_residual_fi(Mdle,Resid, Nref_flag)
!
   use data_structure3D
   use control
   use parametersDPG
   use common_prob_data_UW
!
#include "syscom.blk"
!
!----------------------------------------------------------------------
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
!..geometry
   dimension xi(3),dxidt(3,2),x(3),dxdxi(3,3),dxidx(3,3), dxdt(3,2)
   dimension rt(3,2),rn(3),t(2)
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
!..3D quadrature data
   dimension xiloc(3,MAXNINT3ADD),waloc(MAXNINT3ADD)
!
!..2D quadrature data
   dimension tloc(2,MAXNINT2ADD),wtloc(MAXNINT2ADD)
!
!..BC's flags
   dimension ibc(6,NRINDEX)
!..workspace for trial and test variables
   dimension dq(3) , u(3), dp(1:3), v(3), vec(3), zu(3), zvec(3)
!
   complex*16 :: dfpml(3), Jstretch(3,3), JJStretch(3,3), detJ, zA(3)
!..source
   dimension zf(4)
!
!..added to use fast integration
   VTYPE, allocatable :: AUXHH_A(:,:,:), AUXHH_A0(:)
   VTYPE, allocatable :: AUXVH_A(:,:)  , AUXVH_A0(:,:)
   VTYPE, allocatable :: AUXVV_A(:,:,:), AUXVV_A0(:,:,:)
!
   VTYPE, allocatable :: AUXHH_B(:,:,:,:), AUXHH_B0(:,:)
   VTYPE, allocatable :: AUXVH_B(:,:,:)  , AUXVH_B0(:,:,:)
   VTYPE, allocatable :: AUXVV_B(:,:,:,:), AUXVV_B0(:,:,:,:)
!
   integer :: k1,k2,k3,m1,m2,a,b,sa,sb,sc,m1t,m2t
   integer :: l1,l2,l3,px,py,pz,nintx,ninty,nintz,nord1,nord2,nord3, &
              nrdofH1,nrdofH2,nrdofH3,i1,i2,i3,j1,j2,j3,ltrial,kk
   real*8 :: xi1,xi2,xi3,wt1,wt2,wt3,clock1,clock2,tmp,tmpt
   real*8, dimension(MAXPP+1) :: xilocx,xilocy,xilocz
   real*8, dimension(MAXPP+1) :: wlocx,wlocy,wlocz
   real*8, dimension(3,MAXNINT3ADD) :: wloc3
   real*8, dimension(3) :: xip,dHdx,dHHdx,shapVVt
   real*8, dimension(3,3) :: D,C
   real*8, dimension(MAXPP+1,2) :: shapH1,shapH2,shapH3
   real*8, dimension(MAXPP+1,MAXPP+1) :: sH2p,sH3p,dsH2p,dsH3p
   integer, dimension(3,3) :: deltak
!
!..load vector for the enriched space
   VTYPE, allocatable :: BloadHV(:), Gram(:,:)
!
!..Stiffness matrices for the enriched test space
   integer, allocatable  :: HH_index(:)
! 
   real*8 :: tm,tm2,tm3
   integer(kind=8) :: t1,t2,clock_rate,clock_max
   integer :: nord_add_temp
!..lapack
   character*1 uplo
! 
!..Identity/Kronecker delta tensor
   deltak=0
   do a=1,3
     deltak(a,a)=1
   enddo
!
!----------------------------------------------------------------------
!
!..element type
   ntype = NODES(Mdle)%ntype
   nrv = nvert(ntype); nre = nedge(ntype); nrf = nface(ntype)
!
!..determine order of approximation
   call find_order(Mdle, norder)
!
!   nord_add_temp = max(NORD_ADD,2)
   nord_add_temp = NORD_ADD
!
!..set the enriched order of approximation
   select case(ntype)
      case(MDLB);      nordP = NODES(Mdle)%order + nord_add_temp*111
      case(MDLN,MDLD); nordP = NODES(Mdle)%order + nord_add_temp*1
      case(MDLP);      nordP = NODES(Mdle)%order + nord_add_temp*11
   end select
!
   call compute_enriched_order(nordP, norderP)

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
!..determine solution dof
   call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
!
   call celndof(NODES(Mdle)%ntype,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
   call celndof(NODES(Mdle)%ntype,norderP,nrdofHH,nrdofEE,nrdofVV,nrdofQQ)
!
   nrTest  = nrdofHH + nrdofVV
   nrTrial = nrdofH  + nrdofV  + 4*nrdofQ 
!
!..memory for the matrices
   allocate(BloadHV(nrTest))              ; BloadHV   = ZERO
   allocate(Gram(nrTest,nrTest))          ; Gram      = ZERO
!
!..adjusted frequency for the test space 
!   omeg = min(OMEGA,6.d0/h_elem)
    omeg = OMEGA
! !
!----------------------------------------------------------------------
!     E L E M E N T    I N T E G R A L S                             
!----------------------------------------------------------------------
!
!..here begins the setup for tensorized num quadrature for hexahedra
!..set up the element quadrature
   xiloc  = 0.d0 ;  wloc3  = 0.d0
   xilocx = 0.d0 ;  xilocy = 0.d0 ;  xilocz = 0.d0
   wlocx  = 0.d0 ;  wlocy  = 0.d0 ;  wlocz  = 0.d0        
   sa     = 0    ;  sb     = 0    ;  D      = ZERO
!
!..set up the integration   
   INTEGRATION = nord_add_temp
   call set_3Dint_fi(ntype,norder,nord1,nord2,nord3,nintx,ninty,nintz,xiloc,wloc3)
   INTEGRATION = 0
! 
!..set up # dof for each direction for 1D H1 test functions with order p+dp
   nrdofH1=nord1+1; nrdofH2=nord2+1; nrdofH3=nord3+1
!..set up # dof for each direction for 1D H1 trial functions with order p
   nrdofH1_tr=nrdofH1-nord_add_temp; nrdofH2_tr=nrdofH2-nord_add_temp; nrdofH3_tr=nrdofH3-nord_add_temp
!..set up # dof for each direction for 1D L2 trial functions with order p
   nrdofQ1_tr=nrdofH1_tr-1; nrdofQ2_tr=nrdofH2_tr-1; nrdofQ3_tr=nrdofH3_tr-1
!   
!..Allocate the auxiliary arrays for sum factorization
   allocate(AUXHH_A (3,3,nrdofH3**2))
   allocate(AUXHH_A0(nrdofH3**2))
   allocate(AUXHH_B (3,3,nrdofH2**2,nrdofH3**2))
   allocate(AUXHH_B0(nrdofH2**2,nrdofH3**2))
   allocate(AUXVH_A (3,nrdofH3**2))
   allocate(AUXVH_A0(3,nrdofH3**2))
   allocate(AUXVH_B (3,nrdofH2**2,nrdofH3**2))
   allocate(AUXVH_B0(3,nrdofH2**2,nrdofH3**2))
   allocate(AUXVV_A (3,3,nrdofH3**2))
   allocate(AUXVV_A0(3,3,nrdofH3**2))
   allocate(AUXVV_B (3,3,nrdofH2**2,nrdofH3**2))
   allocate(AUXVV_B0(3,3,nrdofH2**2,nrdofH3**2))

!..total number of test functions 
   nrdofHH= nrdofH1*nrdofH2*nrdofH3
   nrdofVV= nrdofH1*nord2*nord3+nord1*nrdofH2*nord3+nord1*nord2*nrdofH3
   nrdofQ = nrdofQ1_tr*nrdofQ2_tr*nrdofQ3_tr
!
!..Quadrature points and weights are copied into 1-dim vectors
   l=0    
   do l3=1,nintz 
      do l2=1,ninty
         do l1=1,nintx
            l=l+1
            xilocx(l1)=xiloc(1,l)
            xilocy(l2)=xiloc(2,l)
            xilocz(l3)=xiloc(3,l)
            wlocx(l1)=wloc3(1,l)
            wlocy(l2)=wloc3(2,l)
            wlocz(l3)=wloc3(3,l)
         enddo 
      enddo
   enddo
!  
   call system_clock ( t1, clock_rate, clock_max )
!..Loop on quadrature points in direction \xi_1
   do px=1,nintx
      xi1=xilocx(px)
      wt1=wlocx(px)
!      
!  ...call 1D shape functions for coordinate 1
      call shape1HH(xi1,nord1, nrdofH1,shapH1(:,1),shapH1(:,2))
!      
!  ...Initialize auxiliary matrices AUX_B
      AUXHH_B  = ZERO ;  AUXHH_B0  = ZERO ;  AUXVH_B   = ZERO
      AUXVH_B0 = ZERO ;  AUXVV_B   = ZERO ;  AUXVV_B0  = ZERO
!
!  ...Loop through quadrature points in direction \xi_2
      do py=1,ninty
         xi2=xilocy(py)
         wt2=wlocy(py)
!     ...Shape function subroutine is called only once, when 
!        px=1 and stored in sH2p(:,py) and dsH2p(:,py)
         if (px.eq.1) then
            sH2p (:,py) = 0.d0
            dsH2p(:,py) = 0.d0
            call shape1HH(xi2,nord2,nrdofH2,sH2p(:,py),dsH2p(:,py))
         endif
!     ...Copy shape functions in coord. 2 previously evaluated
         shapH2(:,1) = sH2p(:,py)
         shapH2(:,2) = dsH2p(:,py)
!     ...Initialize auxiliary matrices AUX_A
         AUXHH_A  = ZERO ;  AUXHH_A0  = ZERO ; AUXVH_A   = ZERO
         AUXVH_A0 = ZERO ;  AUXVV_A   = ZERO ; AUXVV_A0  = ZERO
!
!     ...Loop on quadrature points in direction \xi_3
         do pz=1,nintz
            xi3=xilocz(pz)
            wt3=wlocz(pz)
            xip(1)=xi1 ;  xip(2)=xi2 ; xip(3)=xi3
!        ...Shape function subroutine is called only once, when 
!           px=py=1 and stored in sH3p(:,pz) and dsH3p(:,pz) 
            if (px*py.eq.1) then
               call shape1HH(xi3,nord3,nrdofH3,sH3p(:,pz),dsH3p(:,pz))
            endif
!        ...Copy shape functions in coord. 3 previously evaluated
            shapH3(:,1)=sH3p(:,pz)
            shapH3(:,2)=dsH3p(:,pz)
!        ...Compute shape functions needed for geometry - 3D H1 shape functions        
            call shape3DH(ntype,xip,norder,norient_edge,norient_face,nrdofH,shapH,gradH)
!        ...Geometry map
            call geom3D(Mdle,xip,xnod,shapH,gradH,nrdofH,x,dxdxi,dxidx,rjac,iflag)
!        ...compute total quadrature weight
            wt123=wt1*wt2*wt3
!        ...compute Jacobian determinant * quadrature weight
            weighthh=wt123*rjac
!        ...Determine D
            D(1,1) = weighthh*(dxidx(1,1)**2+dxidx(1,2)**2+dxidx(1,3)**2)
            D(1,2) = weighthh*(dxidx(1,1)*dxidx(2,1)                     &
                   + dxidx(1,2)*dxidx(2,2)+dxidx(1,3)*dxidx(2,3))
            D(1,3) = weighthh*(dxidx(1,1)*dxidx(3,1)                     &
                   + dxidx(1,2)*dxidx(3,2)+dxidx(1,3)*dxidx(3,3))
            D(2,1) = D(1,2)
            D(2,2) = weighthh*(dxidx(2,1)**2+dxidx(2,2)**2+dxidx(2,3)**2)
            D(2,3) = weighthh*(dxidx(2,1)*dxidx(3,1)                     & 
                   + dxidx(2,2)*dxidx(3,2)+dxidx(2,3)*dxidx(3,3))
            D(3,1) = D(1,3)
            D(3,2) = D(2,3)
            D(3,3) = weighthh*(dxidx(3,1)**2+dxidx(3,2)**2+dxidx(3,3)**2)
!        ...compute inverse Jacobian determinant * quadrature weight
            weightvv = wt123/rjac
!        ...Determine C
            C(1,1) = weightvv*(dxdxi(1,1)**2+dxdxi(2,1)**2+dxdxi(3,1)**2)
            C(1,2) = weightvv*(dxdxi(1,1)*dxdxi(1,2)                     &
                   + dxdxi(2,1)*dxdxi(2,2)+dxdxi(3,1)*dxdxi(3,2))
            C(1,3) = weightvv*(dxdxi(1,1)*dxdxi(1,3)                     &
                   + dxdxi(2,1)*dxdxi(2,3)+dxdxi(3,1)*dxdxi(3,3))
            C(2,1) = C(1,2)
            C(2,2) = weightvv*(dxdxi(1,2)**2+dxdxi(2,2)**2+dxdxi(3,2)**2)
            C(2,3) = weightvv*(dxdxi(1,2)*dxdxi(1,3)                     &
                   + dxdxi(2,2)*dxdxi(2,3)+dxdxi(3,2)*dxdxi(3,3))
            C(3,1) = C(1,3)
            C(3,2) = C(2,3)
            C(3,3) = weightvv*(dxdxi(1,3)**2+dxdxi(2,3)**2+dxdxi(3,3)**2)
!           
!        ...put appropriate quadrature weight on Jacobian and its inverse
            dxdxi = dxdxi * weightvv
            dxidx = dxidx * wt123
!
!        ...HERE STARTS COMPUTATION OF VOLUME INTEGRALS: GRAM, STIFFNESS AND BLOAD MATRICES
!        ...COMPUTATION OF AUXILIARY ARRAYS 'A' FOR GRAM MATRIX - H1,H1 BLOCK
!        ...loop on 1D DOFs, coord. 3, test shape function H1 
            do i3=1,nrdofH3
!           ...loop on 1D DOFs, coord. 3, error repr. trial func H1
               do j3=i3,nrdofH3
!           ...combine indices i3, j3 into k3
                  k3=(i3-1)*nrdofH3+j3
!              ...Only upper-triangular part of matrices is computed
                  do b=1,3
                     sb=1+deltak(b,3)
                     do a=1,3                      
                        sa=1+deltak(a,3)
                        AUXHH_A(a,b,k3) = AUXHH_A(a,b,k3)     &
                                        + shapH3(i3,sa) * shapH3(j3,sb)*D(a,b)
                     enddo
                  enddo
                  AUXHH_A0(k3) = AUXHH_A0(k3)     &
                                + (omeg**2+ALPHA)*shapH3(i3,1)*shapH3(j3,1)*weighthh
!                               + ( omeg**2*detJ*conjg(detJ)*shapH3(i3,1)*shapH3(j3,1) +  &
!                                   ALPHA*shapH3(i3,1)*shapH3(j3,1))*weighthh
!    
!              ...end of j3 loop for H1 trial functions
               enddo
            enddo
! 
!        ...COMPUTATION OF AUXILIARY ARRAYS 'A' FOR GRAM MATRIX - H(div),H1 BLOCK
!
!        ...start new loop on i3 - H1 test functions
            do i3=1,nrdofH3
!           ...start new loop on j3 - H(div) trial functions
               do j3=1,nrdofH3
!              ...combine indices i3, j3 into k3
                  k3=(i3-1)*nrdofH3+j3
!              ...loop over vector components of H(div) test functions
                  do b=1,3                    
!                 ...determine indices sb and idxb for trial shape func (j3,b)
                     sb   = 1+1-deltak(b,3)
                     idxb = j3+1-deltak(b,3)
                     if (idxb.le.nrdofH3) then
                        AUXVH_A0(b,k3) = AUXVH_A0(b,k3)                              &
                                        - shapH3(i3,1 )*shapH3(idxb,2 )*ZIMG*omeg*wt123
!                                       - conjg(detJ)*shapH3(i3,1 )*shapH3(idxb,2 )*ZIMG*omeg*wt123
                        sa=1+deltak(b,3)                                              
                        AUXVH_A (b,k3) = AUXVH_A (b,k3)                              &
                                        + shapH3(i3,sa)*shapH3(idxb,sb)*ZIMG*omeg*wt123
!                                       + conjg(zA(b))*shapH3(i3,sa)*shapH3(idxb,sb)*ZIMG*omeg*wt123
                     endif
                  enddo
!              ...end of j3 loop for H(div) trial functions
               enddo
!           ...end of i3 loop for H1 test functions
            enddo
!
!        ...COMPUTATION OF AUXILIARY ARRAYS 'A' FOR GRAM MATRIX - H(div),H(div) BLOCK
!
!        ...start new loop on i3 - H(div) test functions
            do i3=1,nrdofH3
!           ...start new loop on j3 - H(div) trial functions
               do j3=1,nrdofH3
!           ...combine indices i3, j3 into k3
                  k3=(i3-1)*nrdofH3+j3
!              ...loops over vector components of H(div) trial and test functions
                  do a=1,3 
!                 ...determine indices sa and idxa for test shape func (i3,a)
                     sa = 1+1-deltak(a,3)                      
                     idxa = i3+1-deltak(a,3)  
                     if (idxa.le.nrdofH3) then                                   
!                    ...Only upper-triangular part of matrices is computed. b starts at a
                        do b=a,3
!                       ...determine indices sb and idxb for trial shape func (j3,b)
                           sb   = 1+1-deltak(b,3)
                           idxb = j3+1-deltak(b,3)
                           if (idxb.le.nrdofH3) then
                               AUXVV_A (b,a,k3) = AUXVV_A (b,a,k3)   &
                                                + (omeg**2+ALPHA)*shapH3(idxa,sa)*shapH3(idxb,sb)*C(a,b)
!                              AUXVV_A (b,a,k3) = AUXVV_A (b,a,k3)   &
!                                               + omeg**2   &
!                                               * zA(a)*conjg(zA(b))*shapH3(idxa,sa)*shapH3(idxb,sb)*C(a,b) &
!                                               + ALPHA*shapH3(idxa,sa)*shapH3(idxb,sb)*C(a,b)

                              AUXVV_A0(b,a,k3) = AUXVV_A0(b,a,k3)                   &
                                               + shapH3(idxa,2)*shapH3(idxb,2)*weightvv
                           endif
!                       ...end loop for b
                        enddo
                     endif   
!                 ...end loop for a
                  enddo
!              ...end of j3 loop for H(div) trial functions
               enddo
!           ...end of i3 loop for H(div) test functions
            enddo
!
!        ...pz loop ends
         enddo         
!
! COMPUTATION OF AUXILIARY ARRAYS 'B' FOR GRAM MATRIX - H1,H1 BLOCK
!
!     ...loop on 1D DOFs, coord. 3, test shape function - H1  
         do i3=1,nrdofH3
!        ...loop on 1D DOFs, coord. 3, error repr. shape func - H1
            do j3=i3,nrdofH3
!           ...combine indices i3, j3 into k3
               k3=(i3-1)*nrdofH3+j3
!           ...loop on 1D DOFs, coord. 2, test shape function - H1
               do i2=1,nrdofH2
!              ...loop on 1D DOFs, coord. 2, error repr. shape fn - H1
                  do j2=1,nrdofH2
!                 ...combine indices i2, j2 into k2  
                     k2=(i2-1)*nrdofH2+j2
                     do b=1,3
                        sb=1+deltak(b,2)
                        do a=1,3                          
                           sa=1+deltak(a,2)
                           AUXHH_B(a,b,k2,k3) = AUXHH_B(a,b,k2,k3)              &
                                              + shapH2(i2,sa) * shapH2(j2,sb)   &
                                              * AUXHH_A(a,b,k3)
                        enddo
                     enddo
                     AUXHH_B0(k2,k3) = AUXHH_B0(k2,k3)            &
                                     + shapH2(i2,1)*shapH2(j2,1)*AUXHH_A0(k3)
!                 ...j2 loop ends                      
                  enddo
!              ...i2 loop ends
               enddo
!           ...j3 loop ends
            enddo
!        ...i3 loop ends
         enddo
! 
!     ...COMPUTATION OF AUXILIARY ARRAYS 'B' FOR GRAM MATRIX - H(div),H1 BLOCK
! 
!     ...start new loop on i3 - H1 test functions
         do i3=1,nrdofH3
!        ...start new loop on j3 - H(div) trial functions
            do j3=1,nrdofH3
!           ...combine indices i3, j3 into k3
               k3=(i3-1)*nrdofH3+j3
!           ...start new loop on i2 - H1 test functions
               do i2=1,nrdofH2
!              ...start new loop on j2 - H(div) trial functions
                  do j2=1,nrdofH2
!              ...combine indices i2, j2 into k2
                     k2=(i2-1)*nrdofH2+j2
                     do b=1,3                    
!                    ...determine indices sb and idxb for trial shape func (j2,b)
                        sb=1+1-deltak(b,2)
                        idxb=j2+1-deltak(b,2)
                        if (idxb.le.nrdofH2) then 
                           AUXVH_B0(b,k2,k3) = AUXVH_B0(b,k2,k3)              &
                                             + shapH2(i2,1)*shapH2(idxb,2)    &
                                             * AUXVH_A0(b,k3)
                           sa=1+deltak(b,2) 
                           AUXVH_B (b,k2,k3) = AUXVH_B (b,k2,k3)              &
                                             + shapH2(i2,sa)*shapH2(idxb,sb)  &
                                             * AUXVH_A (b,k3)
                        endif
!                    ...end b loop for H(div) shape function 
                     enddo
!                 ...j2 loop ends for H(div) shape function
                  enddo
!              ...i2 loop ends for H1 shape function
               enddo
!           ...j3 loop ends for H(div) shape function
            enddo
!        ...i3 loop ends for H1 shape function
         enddo
!
! COMPUTATION OF AUXILIARY ARRAYS 'B' FOR GRAM MATRIX - H(div),H(div) BLOCK
!
!     ...loop on i3 test function - H(div)
         do i3=1,nrdofH3
!        ...loop on j3 trial functions - H(div)
            do j3=1,nrdofH3
!           ...combine indices i3, j3 into k3
               k3=(i3-1)*nrdofH3+j3
!           ...loop on 1D DOFs, coord. 2, test shape function - H(div)
               do i2=1,nrdofH2
!              ...loop on 1D DOFs, coord. 2, error repr. shape fn - H(div)
                  do j2=1,nrdofH2
!                 ...combine indices i2, j2 into k2
                     k2=(i2-1)*nrdofH2+j2                     
                     do a=1,3    
!                    ...determine indices sa and idxa for test shape func (i2,a)
                        sa=1+1-deltak(a,2)                      
                        idxa=i2+1-deltak(a,2)                  
                        if (idxa.le.nrdofH2) then
                           do b=a,3 
!                          ...determine indices sb and idxb for trial shape func (j2,b)
                              sb=1+1-deltak(b,2)
                              idxb=j2+1-deltak(b,2)
                              if (idxb.le.nrdofH2) then
                                 AUXVV_B (b,a,k2,k3) = AUXVV_B (b,a,k2,k3)             &
                                                     + shapH2(idxa,sa)*shapH2(idxb,sb) &
                                                     * AUXVV_A (b,a,k3)
                                 AUXVV_B0(b,a,k2,k3) = AUXVV_B0(b,a,k2,k3)             &
                                                     + shapH2(idxa,2 )*shapH2(idxb,2 ) &
                                                     * AUXVV_A0(b,a,k3)
                              endif
!                          ...b loop ends
                           enddo
                        endif   
!                    ...a loop ends
                     enddo
!                 ...j2 loop ends
                  enddo
!              ...i2 loop ends
               enddo
!           ...j3 loop ends
            enddo
!        ...i3 loop ends
         enddo
!     ...py loop ends
      enddo
!  
!  ...FINAL COMPUTATION OF GRAM MATRIX - H1,H1 BLOCK
!
!  ...loop on H1 trial shape function identified by indices j1,j2,j3
      do j3=1,nrdofH3
         do j2=1,nrdofH2
            do j1=1,nrdofH1
               m2=j1+(j2-1)*nrdofH1+(j3-1)*nrdofH2*nrdofH1
!           ...loop on H1 test shape function identified by indices i1,i2,i3
               do i3=1,nrdofH3
!              ...combine indices i3, j3 into k3
                  k3=(i3-1)*nrdofH3+j3  
                  do i2=1,nrdofH2
!                 ...combine indices i2, j2 into k2
                     k2=(i2-1)*nrdofH2+j2
                     do i1=1,nrdofH1                      
!                    ...determine unique indices for H1 function
                        m1=i1+(i2-1)*nrdofH1+(i3-1)*nrdofH2*nrdofH1
!                    ...Only upper-triangular part of Gram matrix computed
                        if (m2.ge.m1) then
!                       ...loops over components of gradient of functions J and I
                           do b=1,3
                              sb=1+deltak(b,1)                            
                              do a=1,3
                                 sa=1+deltak(a,1)
                                 Gram(m1,m2) = Gram(m1,m2)       &
                                             + shapH1(i1,sa)*shapH1(j1,sb) &
                                             * AUXHH_B(a,b,k2,k3)
                              enddo
                           enddo                  
                           Gram(m1,m2) = Gram(m1,m2)       &
                                  + shapH1(i1,1)*shapH1(j1,1)*AUXHH_B0(k2,k3)
                        endif
!                    ...i1 loop ends    
                     enddo
!                 ...i2 loop ends
                  enddo
!              ...i3 loop ends
               enddo
!           ...j1 loop ends
            enddo
!        ...j2 loop ends
         enddo
!     ...j3 loop ends
      enddo
!
!  ...FINAL COMPUTATION OF GRAM MATRIX - H(div),H1 BLOCK
!
!  ...start new loop on H(div) trial functions identified by indices j1,j2,j3,b
      do b=1,3
         sb=1+1-deltak(1,b)
         sa=1+deltak(1,b) 
         do j3=1,nrdofH3
            idxb3=j3+1-deltak(3,b)
            if (idxb3.le.nrdofH3) then
               do j2=1,nrdofH2
                  idxb2=j2+1-deltak(2,b)
                  if (idxb2.le.nrdofH2) then
                     do j1=1,nrdofH1
                        idxb1=j1+1-deltak(1,b)
                        if (idxb1.le.nrdofH1) then
!                       ...determine unique index for H(div) function
                           select case(b)
                           case(1)
                              m2 = j1+nrdofH1*(j2-1)+nrdofH1*nord2*(j3-1)
                           case(2)
                              m2 = nrdofH1*nord2*nord3                        &
                                 + j1+nord1*(j2-1)+nord1*nrdofH2*(j3-1)
                           case(3)
                              m2 = nrdofH1*nord2*nord3+nord1*nrdofH2*nord3    &
                                 + j1+nord1*(j2-1)+nord1*nord2*(j3-1)
                           end select
                           n2 = nrdofHH+m2
!                       ...loop on H1 test shape function identified by indices i1,i2,i3
                           do i3=1,nrdofH3
!                          ...combine indices i3, j3 into k3
                              k3=(i3-1)*nrdofH3+j3                     
                              do i2=1,nrdofH2
!                             ...combine indices i2, j2 into k2
                                 k2=(i2-1)*nrdofH2+j2                        
                                 do i1=1,nrdofH1
!                                ...determine unique index for H1 function
                                    m1=i1+(i2-1)*nrdofH1+(i3-1)*nrdofH2*nrdofH1
!
                                    Gram(m1,n2) = Gram(m1,n2)              &
                                                + shapH1(i1,1 )*shapH1(idxb1,2 )*AUXVH_B0(b,k2,k3)
                                    Gram(m1,n2) = Gram(m1,n2)              &
                                                + shapH1(i1,sa)*shapH1(idxb1,sb)*AUXVH_B (b,k2,k3)
!                                ...i1 loop ends
                                 enddo
!                             ...i2 loop ends
                              enddo
!                          ...i3 loop ends
                           enddo
                        endif   
!                    ...j1 loop ends
                     enddo
                  endif   
!              ...j2 loop ends
               enddo
            endif   
!        ...j3 loop ends
         enddo
!        ...b loop ends
      enddo
!
!  ...FINAL COMPUTATION OF GRAM MATRIX - H(div), H(div) BLOCK
!
!  ...loop on trial shape function H(div) of indices by j1,j2,j3,b
      do b=1,3
         sb=1+1-deltak(b,1)
         do j3=1,nrdofH3
            idxb3=j3+1-deltak(b,3)
            if (idxb3.le.nrdofH3) then
               do j2=1,nrdofH2
                  idxb2=j2+1-deltak(b,2)
                  if (idxb2.le.nrdofH2) then
                     do j1=1,nrdofH1
                        idxb1=j1+1-deltak(b,1)
                        if (idxb1.le.nrdofH1) then
!                       ...loop on test shape function H(div) of indices i1,i2,i3,a
                           do a=1,b
                              sa=1+1-deltak(a,1)                      
                              do i3=1,nrdofH3   
                                 idxa3=i3+1-deltak(a,3)
                                 if (idxa3.le.nrdofH3) then
                                    k3=(i3-1)*nrdofH3+j3
                                    do i2=1,nrdofH2                   
                                       idxa2=i2+1-deltak(a,2)
                                       if (idxa2.le.nrdofH2) then
                                          k2=(i2-1)*nrdofH2+j2
                                          do i1=1,nrdofH1
                                             idxa1=i1+1-deltak(a,1)
                                             if (idxa1.le.nrdofH1) then
!!                                           ...determine unique indices for H(div) functions
                                                select case(a)
                                                case(1)
                                                   m1 = i1+nrdofH1*(i2-1)               &
                                                     + nrdofH1*nord2*(i3-1)
                                                case(2)
                                                   m1 = nrdofH1*nord2*nord3             &
                                                      + i1+nord1*(i2-1)                 &
                                                      + nord1*nrdofH2*(i3-1)
                                                case(3)
                                                   m1 = nrdofH1*nord2*nord3   &
                                                      + nord1*nrdofH2*nord3   &
                                                      + i1+nord1*(i2-1)+nord1*nord2*(i3-1)
                                                end select
                                                n1 = nrdofHH+m1
                                                select case(b)
                                                case(1)
                                                   m2 = j1+nrdofH1*(j2-1)+nrdofH1*nord2*(j3-1)
                                                case(2)
                                                   m2 = nrdofH1*nord2*nord3             &
                                                      + j1+nord1*(j2-1)+nord1*nrdofH2*(j3-1)
                                                case(3)
                                                   m2 = nrdofH1*nord2*nord3+nord1*nrdofH2*nord3    &
                                                      + j1+nord1*(j2-1)+nord1*nord2*(j3-1)
                                                end select
                                                n2 = nrdofHH+m2
!                                            ...fill only the upper part of the Gram matrix block    
                                                if (n2.ge.n1) then
                                                   Gram(n1,n2) = Gram(n1,n2) &
                                                      + shapH1(idxa1,sa)*shapH1(idxb1,sb)  &
                                                      * AUXVV_B (b,a,k2,k3)                &
                                                      + shapH1(idxa1,2 )*shapH1(idxb1,2 )  &
                                                      * AUXVV_B0(b,a,k2,k3)
                                                endif                                
                                             endif
            !                             ...i1 loop ends
                                          enddo
                                       endif 
         !                          ...i2 loop ends
                                    enddo
                                 endif 
!                                ...i3 loop ends
                              enddo
!                          ...a loop ends
                           enddo
                        endif 
!                    ...j1 loop ends
                     enddo
                  endif 
!              ...j2 loop ends
               enddo
            endif 
!        ...j3 loop ends
         enddo
!     ...b loop ends
      enddo
!  ...px loop ends
   enddo
! ..deallocate auxiliary arrays memory
   deallocate(AUXHH_A,AUXHH_A0,AUXHH_B,AUXHH_B0)
   deallocate(AUXVH_A,AUXVH_A0,AUXVH_B,AUXVH_B0)
   deallocate(AUXVV_A,AUXVV_A0,AUXVV_B,AUXVV_B0)
!
   INTEGRATION = nord_add_temp
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
                     + (ZIMG*OMEGA*zp*q                     &
                     -  zu(1)*dq(1)-zu(2)*dq(2)-zu(3)*dq(3) &
                     -  zf(1)*q)                            &
                     *  weight
!          
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
                     + (ZIMG*OMEGA*(zu(1)*v(1)+       &
                                    zu(2)*v(2)+       &
                                    zu(3)*v(3))       &
                     -  zp*div_v                                     &
                     -  zf(2)*v(1)- zf(3)*v(2) - zf(4)*v(3))         &
                     * weight
!         
!     ...end of loop through H(div) test functions      
      enddo
!  ...end of loop through integration points      
   enddo

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
      INTEGRATION = nord_add_temp
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
   uplo = 'U'
!
!..factorize the test Stiffness matrix
   call ZPOTRF(uplo,nrTest,Gram,nrTest,info)
!   
   if (info.ne.0) then
      write(*,*) 'elem_residual_fi: Gram ZPOTRF: Mdle,info = ',Mdle,info
      stop 1
   endif
!
   call ZTRSM('L',uplo,'C','N',nrTest,1,ZONE,Gram,nrTest,BloadHV,nrTest)
!
   Resid = ZERO
   do i=1,nrTest
      Resid = Resid + conjg(BloadHV(i))*BloadHV(i)
   enddo
!
   Nref_flag = 111
   deallocate(Gram,BloadHV)


   end subroutine elem_residual_fi


 
