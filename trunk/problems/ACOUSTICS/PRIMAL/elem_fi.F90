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
   subroutine elem_fi_primal_acoustics(Mdle)
!
   use primal_acoustics_module
   use control, only: INTEGRATION
   use parameters
   use parametersDPG
   use data_structure3D
   use element_data
   use physics   , only : NR_PHYSA
   use assembly  , only : ALOC,BLOC,NR_RHS
   use common_prob_data, only: SYMMETRY_TOL, TEST_NORM, OMEGA,ALPHA
!
#include "syscom.blk"
!
!----------------------------------------------------------------------
!
   complex*16, dimension(3,3)                         :: AUX_A
   complex*16, dimension(3,3,(MAXPP+1)**2)            :: AUX_B
   complex*16                                         :: AUX_A0
   complex*16, dimension((MAXPP+1)**2)                :: AUX_B0
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
   dimension xi(3),dxidt(3,2),x(3),dxdxi(3,3),dxidx(3,3), dxdt(3,2),rt(3,2),rn(3),t(2)
!
!..H1 shape functions
   dimension shapH(MAXbrickH),gradH(3,MAXbrickH)
!
!..Hdiv shape functions
   dimension shapV(3,MAXbrickV),divV(MAXbrickV)
!
!..3D quadrature data
   dimension xiloc(3,MAXNINT3ADD),waloc(MAXNINT3ADD)
!
!..2D quadrature data
   dimension tloc(2,MAXNINT2ADD),wtloc(MAXNINT2ADD)
!
!..BC's flags
   dimension ibc(6,NR_PHYSA)
!
!..derivatives wrt physical coordinates, flux
   dimension dv1(3),dv2(3),vec(3)
!
!..for debug printing
   dimension aux(10)
!

   integer :: k1,k2,k3,m1,m2,a,b,sa,sb,sc,m1t,m2t
!  ...added to use fast integration
   integer :: l1,l2,l3,px,py,pz,nintx,ninty,nintz,nord1,nord2,nord3, &
              nrdofH1,nrdofH2,nrdofH3,i1,i2,i3,j1,j2,j3,ltrial,kk
   real*8 :: xi1,xi2,xi3,wt1,wt2,wt3,clock1,clock2,tmp,tmpt
!       real*8 :: E
   real*8, dimension(MAXPP+1) :: xilocx,xilocy,xilocz
   real*8, dimension(MAXPP+1) :: wlocx,wlocy,wlocz
   real*8, dimension(3,MAXNINT3ADD) :: wloc3
   real*8, dimension(3) :: xip,dHdx,dHHdx
   real*8, dimension(3,3) :: D
   real*8, dimension(MAXPP,2) :: shapH1,shapH2,shapH3
   real*8, dimension(MAXPP,MAXPP+1) :: sH2p,sH3p,dsH2p,dsH3p
   integer, dimension(3,3) :: deltak

   character*1 uplo
   nk(k1,k2) = (k2-1)*k2/2+k1

   real*8 :: Mtime(10)
   integer(kind=8) :: t1,t2,clock_rate,clock_max
!
   deltak=ZERO
   do a=1,3
     deltak(a,a)=1
   enddo
!---------------------------------------------------------------------
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
!..get the element boundary conditions flags
   call find_bc(Mdle, ibc)
!                                                                    
!..clear space for stiffness matrix and rhsv:                     
   ZblocH = ZERO; ZblocV = ZERO
   ZalocHH = ZERO; ZalocHV = ZERO; ZalocVH = ZERO; ZalocVV = ZERO
!
!..extended load vector and extended stiffness matrices
   BLOADH=ZERO ; STIFFHH=ZERO ; STIFFHV=ZERO
!
!..Gram matrix
   AP=ZERO ; 
!
!
!-----------------------------------------------------------------------
!     E L E M E N T    I N T E G R A L S                               |
!-----------------------------------------------------------------------
!
!..here begins the setup for tensorized num quadrature for hexahedra
!..set up the element quadrature
   xiloc(:,:)=0.d0
   wloc3(:,:)=0.d0
   xilocx(:)=0.d0
   xilocy(:)=0.d0
   xilocz(:)=0.d0
   wlocx(:)=0.d0
   wlocy(:)=0.d0
   wlocz(:)=0.d0        
   sa=0
   sb=0
   D=ZERO
   INTEGRATION = NORD_ADD
   call set_3Dint_fi(etype,norder,nord1,nord2,nord3,nintx,ninty,nintz,xiloc,wloc3)
   INTEGRATION = 0
! 
!..set up # dof for each direction
   nrdofH1=nord1+1; nrdofH2=nord2+1; nrdofH3=nord3+1
   nrdofHH= nrdofH1*nrdofH2*nrdofH3
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
   xip(:) = 0.d0

   call system_clock ( t1, clock_rate, clock_max )
!..Loop on quadrature points in direction \xi_1
   do px=1,nintx
      xi1=xilocx(px)
      wt1=wlocx(px)
      call shape1H(xi1,nord1, nrdofH1,shapH1(:,1),shapH1(:,2))
!      
!  ...loop on 1D DOFs, coord. 3, test shape function
!  
      do i3=1,nrdofH3
!     ...loop on 1D DOFs, coord. 3, error repr. shape func
         do j3=i3,nrdofH3
!        ...Initialize auxiliary matrices AUX_B
            AUX_B=ZERO
            AUX_B0=ZERO
!        ...combine indices i3, j3 into k3 and its 'transpose' l3
!        ...Loop on quadrature points in direction \xi_2
            do py=1,ninty
               xi2=xilocy(py)
               wt2=wlocy(py)
!           ...Shape function subroutine is called only once, when 
!              i3=j3=1 and stored in sH2p(:,py) and dsH2p(:,py)
               if (i3*j3.eq.1) then
                  sH2p(:,py)=0.d0
                  dsH2p(:,py)=0.d0
                  call shape1H(xi2,nord2,nrdofH2,sH2p(:,py),dsH2p(:,py))
               endif
!           ...Copy shape functions in coord. 2 previously evaluated
               shapH2(:,1)=sH2p(:,py)
               shapH2(:,2)=dsH2p(:,py)
!           ...Initialize auxiliary matrices AUX_A
               AUX_A=ZERO
               AUX_A0=ZERO
!            ..Loop on quadrature points in direction \xi_3
               do pz=1,nintz
                  xi3=xilocz(pz)
                  wt3=wlocz(pz)
                  xip(1)=xi1
                  xip(2)=xi2
                  xip(3)=xi3
! 
                  if (i3*j3.eq.1) then
!                 ...Shape function subroutine is called only once, when 
!                    i3=j3=1 and stored in sH3p(:,pz) and dsH3p(:,pz)
                     call shape1H(xi3,nord3,nrdofH3,sH3p(:,pz),dsH3p(:,pz))
                  endif
!              ...Compute shape functions needed for trial field vars
!              ...H1 (trial/geometry)
                  call shape3H(etype,xip,norder,norient_edge,norient_face,nrdofH,shapH,gradH)
!              ...Geometry map
                  call geom3D(Mdle,xip,xnod,shapH,gradH,nrdofH,x,dxdxi,dxidx,rjac,iflag)
                  if (iflag.ne.0) then
                    write(*,5999) Mdle,rjac
 5999                 format('Negative JacobiancMdle,rjac=',i8,2x,e12.5)
                    stop
                  endif
!              ...Evaluate forcing function f at physical point x
                  call getf(Mdle,x, zfval)
!              ...Determine D
                  D(1,1)=rjac*(dxidx(1,1)**2+dxidx(1,2)**2+dxidx(1,3)**2)
                  D(1,2)=rjac*(dxidx(1,1)*dxidx(2,1)+ &
                          dxidx(1,2)*dxidx(2,2)+dxidx(1,3)*dxidx(2,3))
                  D(1,3)=rjac*(dxidx(1,1)*dxidx(3,1)+ &
                          dxidx(1,2)*dxidx(3,2)+dxidx(1,3)*dxidx(3,3))
                  D(2,1)=D(1,2)
                  D(2,2)=rjac*(dxidx(2,1)**2+dxidx(2,2)**2+dxidx(2,3)**2)
                  D(2,3)=rjac*(dxidx(2,1)*dxidx(3,1)+ &
                          dxidx(2,2)*dxidx(3,2)+dxidx(2,3)*dxidx(3,3))
                  D(3,1)=D(1,3)
                  D(3,2)=D(2,3)
                  D(3,3)=rjac*(dxidx(3,1)**2+dxidx(3,2)**2+dxidx(3,3)**2)
!              ...Copy shape functions in coord. 3 previously evaluated
                  shapH3(:,1)=sH3p(:,pz)
                  shapH3(:,2)=dsH3p(:,pz)
!              ...Only upper-triangular part of matrices is computed
                  do b=1,3
                     sb=1+deltak(b,3)
                     do a=1,3                      
                        sa=1+deltak(a,3)
                        AUX_A(a,b) = AUX_A(a,b)     &
                                   + shapH3(i3,sa) * shapH3(j3,sb)*D(a,b)*wt3
                     enddo
                  enddo
                  AUX_A0 = AUX_A0 + shapH3(i3,1)*shapH3(j3,1)*rjac*wt3
!    
!              ...Integration of enriched stiffness matrix and load vector.
!              ...This is computed only once, when j3=i3 (1st value of j3)
                  if (j3.eq.i3) then
 !                ...loop on 1D DOFs, coord. 2, test shape function
                     do i2=1,nrdofH2
 !                   ...loop on 1D DOFs, coord. 1, test shape function
                        do i1=1,nrdofH1
                           m1=i1+(i2-1)*nrdofH1+(i3-1)*nrdofH2*nrdofH1
                           v2=shapH1(i1,1)*shapH2(i2,1)*shapH3(i3,1)
                           dHHdx(1:3) = shapH1(i1,2)*shapH2(i2,1)*shapH3(i3,1)*  &
                                        dxidx(1,1:3)+  &
                                        shapH1(i1,1)*shapH2(i2,2)*shapH3(i3,1)*  &
                                        dxidx(2,1:3)+  &
                                        shapH1(i1,1)*shapH2(i2,1)*shapH3(i3,2)*  &
                                        dxidx(3,1:3)
                           weight=wt1*wt2*wt3*rjac
 !                      ...loop on trial shape functions
                           do ltrial=1,nrdofH
                              v1=shapH(ltrial)                           
!
                              dHdx(1:3) =  gradH(1,ltrial)*dxidx(1,1:3) &
                                        +gradH(2,ltrial)*dxidx(2,1:3) &
                                        +gradH(3,ltrial)*dxidx(3,1:3)
                              STIFFHH(m1,ltrial) = STIFFHH(m1,ltrial)   &
                                                 + (dHdx(1)*dHHdx(1)      &
                                                 + dHdx(2)*dHHdx(2)       &
                                                 + dHdx(3)*dHHdx(3)       &
                                                 - OMEGA**2*v1*v2)*weight
                           enddo
                           BLOADH(m1)=BLOADH(m1)+ zfval*v2*weight
                        enddo
                     enddo
                  endif
!              ...pz loop ends
               enddo
!           ...loop on 1D DOFs, coord. 2, test shape function
               do j2=1,nrdofH2
!              ...loop on 1D DOFs, coord. 2, error repr. shape fn
                  do i2=1,nrdofH2
!                 ...combine indices i2, j2 into k2
                     k2=(i2-1)*nrdofH2+j2
                     l2=(j2-1)*nrdofH2+i2
                     do b=1,3
                        sb=1+deltak(2,b)
                        do a=1,3                          
                           sa=1+deltak(2,a)
                           AUX_B(a,b,k2) = AUX_B(a,b,k2)     &
                                         + shapH2(i2,sa)*shapH2(j2,sb)*AUX_A(a,b)*wt2
                        enddo
                     enddo
                     if (i2.ge.j2) then
                        AUX_B0(k2)=AUX_B0(k2) + shapH2(i2,1)*shapH2(j2,1)*AUX_A0*wt2
                     endif
!                 ...j2 loop ends                      
                  enddo
!              ...i2 loop ends
               enddo
!           ...py loop ends
            enddo
!  
!           ...loop on 1D DOFs, coord. 2, test shape function
               do j2=1,nrdofH2
!              ...loop on 1D DOFs, coord. 2, error repr. shape fn
                  do i2=1,nrdofH2
!                 ...combine indices i2, j2 into k2 and its 'transpose' l2
                     k2=(i2-1)*nrdofH2+j2
                     l2=(j2-1)*nrdofH2+i2
                     if (i2.lt.j2) then
                        AUX_B0(k2)=AUX_B0(l2) 
                     endif
!                 ...loop on 1D DOFs, coord. 1, test shape function
                  do j1=1,nrdofH1
!                 ...loop on 1D DOFs, coord. 1, error repr. shape fn
                     do i1=1,nrdofH1                      
                        m1=i1+(i2-1)*nrdofH1+(i3-1)*nrdofH2*nrdofH1
                        m2=j1+(j2-1)*nrdofH1+(j3-1)*nrdofH2*nrdofH1                      
!                    ...Only upper-triangular part of Gram matrix computed
                        if (m2.ge.m1) then
                           kk = nk(m1,m2)
                           tmp=0.d0
                           do b=1,3
                              sb=1+deltak(1,b)                            
                              do a=1,3
                                 sa=1+deltak(1,a)
                                 tmp=tmp+shapH1(i1,sa)*shapH1(j1,sb)*AUX_B(a,b,k2)
                              enddo
                           enddo                  
                           tmp=tmp+shapH1(i1,1)*shapH1(j1,1)*AUX_B0(k2)
                           AP(kk)=AP(kk) + ALPHA*tmp*wt1
                        endif
!                    ...j1 loop ends    
                     enddo
!                 ...i1 loop ends
                  enddo
!              ...j2 loop ends
               enddo
!           ...i2 loop ends
            enddo
!        ...j3 loop ends
         enddo
!     ...i3 loop ends
      enddo
!  ...px loop ends
   enddo

   call system_clock ( t2, clock_rate, clock_max )
   Mtime(1) =  real(t2 - t1,8)/real(clock_rate,8)

!
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
      call set_2Dint(ftype,norderf, nint,tloc,wtloc)
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
         call shape3V(etype,xi,norder,norient_face, nrdofV,shapV,divV)
!
!     ...geometry
         call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign,x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
         weight = bjac*wtloc(l)
! 
!     ...loop through enriched H1 test functions
         do k1=1,nrdofHH
            v1 = shapHH(k1)
!
!        ...loop through H(div) trial functions
            do k2=1,nrdofV
! 
!           ...normal component (Piola transformation at work!)
               vec(1:3) = dxdxi(1:3,1)*shapV(1,k2)   &
                        + dxdxi(1:3,2)*shapV(2,k2)   & 
                        + dxdxi(1:3,3)*shapV(3,k2)
               vec(1:3) = vec(1:3)/rjac
               v2n = vec(1)*rn(1)+vec(2)*rn(2)+vec(3)*rn(3)
!
!           ...accumulate for the stiffness matrix
               STIFFHV(k1,k2) = STIFFHV(k1,k2) - v1*v2n*weight
            enddo
         enddo
      enddo
   enddo
!
!-----------------------------------------------------------------------
!      Alternative construction of normal matrix
!-----------------------------------------------------------------------
!
   call celndof(NODES(Mdle)%type,norder, nrdofH,nrdofE,nrdofV,nrdofQ)

   i1 = nrdofHH ; j1 = nrdofH ; j2 = nrdofV 

   STIFF_ALL(1:i1,1:j1) = STIFFHH(1:i1,1:j1)
   STIFF_ALL(1:i1,j1+1:j1+j2) = STIFFHV(1:i1,1:j2)
   STIFF_ALL(1:i1,j1+j2+1) = BLOADH(1:i1)

   nrTEST = nrdofHH
   N      = nrTEST
   NRHS   = nrdofH + nrdofV + 1

!..diagonal scaling of the gram matrix
   call diag_scaling(N,NRHS,AP(1:N*(N+1)/2), STIFF_ALL(1:N,1:NRHS))
   uplo = 'U'
! 
#if C_MODE
!
!..factorize the test stiffness matrix
   call ZPPTRF(uplo, nrTEST, AP, info)
   if (info.ne.0) then
      write(*,*) 'elem: AP ZPPTRF: Mdle,info = ',Mdle,info
      stop 1
   endif

   call ZTPTRS('U','C','N',N,NRHS,AP,STIFF_ALL,MAXbrickHH,info)
   if (info.ne.0) then
      write(*,*) 'elem: AP ZTPTRS: Mdle,info = ',Mdle,info
      stop 2
   endif
   
   call ZHERK('U','C',NRHS,nrTEST,ZONE,STIFF_ALL,MAXbrickHH,ZERO,   &
              STIFF_ALL(1:NRHS,1:NRHS),NRHS)
! 
   do i=1,NRHS-1
      STIFF_ALL(i+1:NRHS,i) = conjg(STIFF_ALL(i,i+1:NRHS))
   enddo
! 
#else
!
!..factorize the test stiffness matrix
   call DPPTRF(uplo, nrTEST, AP, info)
   if (info.ne.0) then
      write(*,*) 'elem: AP DPPTRF: Mdle,info = ',Mdle,info
      stop 1
   endif

   call DTPTRS('U','C','N',N,NRHS,AP,STIFF_ALL,MAXbrickHH,info)
   if (info.ne.0) then
      write(*,*) 'elem: AP DTPTRS: Mdle,info = ',Mdle,info
      stop 2
   endif
   
   call DSYRK('U','C',NRHS,nrTEST,ZONE,STIFF_ALL,MAXbrickHH,ZERO,   &
              STIFF_ALL(1:NRHS,1:NRHS),NRHS)
! 
   do i=1,NRHS-1
      STIFF_ALL(i+1:NRHS,i) = STIFF_ALL(i,i+1:NRHS)
   enddo
#endif
! 
   BLOC(1)%array(1:j1,1)      = STIFF_ALL(1:j1,j1+j2+1)
   BLOC(2)%array(1:j2,1)      = STIFF_ALL(j1+1:j1+j2,j1+j2+1) 
   ALOC(1,1)%array(1:j1,1:j1) = STIFF_ALL(1:j1,1:j1)
   ALOC(1,2)%array(1:j1,1:j2) = STIFF_ALL(1:j1,j1+1:j1+j2)
   ALOC(2,1)%array(1:j2,1:j1) = STIFF_ALL(j1+1:j1+j2,1:j1)
   ALOC(2,2)%array(1:j2,1:j2) = STIFF_ALL(j1+1:j1+j2,j1+1:j1+j2)
! 
! 
   write(*,1001) Mtime(1)
 1001 format('Fast integration            = ', f13.5)  

   end subroutine elem_fi_primal_acoustics


 