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
   subroutine elem_DPG_UWEAK_ACOUSTICS_fi(Mdle,Nord_add_local,NrTest,NrTrial,      &
                                          NrdofHH,NrdofVV,NrdofH,NrdofV,NrdofQ,    &
                                          NrdofHi, NrdofVi,                        &   
                                          NrdofH1,    NrdofH2,    NrdofH3,         &
                                          NrdofQ1_tr, NrdofQ2_tr, NrdofQ3_tr)
!
   use data_structure3D
   use control,     only: INTEGRATION, ISTC_FLAG
   use parametersDPG
   use physics,     only: NR_PHYSA
   use assembly,    only: ALOC,BLOC,NR_RHS
   use common_prob_data_UW
!
!----------------------------------------------------------------------
!
   implicit none

   integer, intent(in) :: Mdle, Nord_add_local 
   integer, intent(in) :: NrTest, NrTrial
   integer, intent(in) :: NrdofHH,NrdofVV, NrdofH, NrdofV, NrdofQ
   integer, intent(in) :: NrdofHi, NrdofVi
   integer, intent(in) :: NrdofH1, NrdofH2, NrdofH3
   integer, intent(in) :: NrdofQ1_tr, NrdofQ2_tr, NrdofQ3_tr

   integer :: ntype,ftype
!
   integer :: nrv, nre, nrf

!
!!..element order, orientation for edges and faces
   integer :: norder(19),norient_edge(12),norient_face(6),norderi(19)
!
!..enriched order
   integer :: nordP, norderP(19)
!..face order
   integer :: norderf(5)
!   
!..geometry dof
   real*8  :: xnod(3,MAXbrickH)
! 
!..geometry
   real*8  :: xi(3),dxidt(3,2),x(3),dxdxi(3,3),dxidx(3,3), dxdt(3,2)
   real*8  :: rt(3,2),rn(3),t(2)
!
!..H1 shape functions
   real*8  :: shapH (MAXbrickH)  , gradH(3,MAXbrickH)
   real*8  :: shapHH(MAXbrickHH) , gradHH(3,MAXbrickHH)
!..Hdiv shape functions
   real*8  :: shapV (3,MAXbrickV) , divV(MAXbrickV) 
   real*8  :: shapVV(3,MAXbrickVV), divVV(MAXbrickVV)
!      
!..3D quadrature data
   integer :: nint
   real*8  :: xiloc(3,MAXNINT3ADD),waloc(MAXNINT3ADD)
!
!..2D quadrature data
   real*8  :: tloc(2,MAXNINT2ADD),wtloc(MAXNINT2ADD)
!
   real*8  :: bjac, rjac, weight, wa
!
!..BC's flags
   integer :: ibc(6,NRINDEX)
!..workspace for trial and test variables
   real*8  :: dq(3) , u(3), dp(1:3), v(3), vec(3), p, q, un, vn
!
!..source
   VTYPE :: zf(4), zg
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
   VTYPE, allocatable :: STIFQH_A(:),      STIFQH_B(:,:)
   VTYPE, allocatable :: STIFQ3H_A(:,:,:), STIFQ3H_B(:,:,:,:)
   VTYPE, allocatable :: STIFQV_A(:,:),    STIFQV_B(:,:,:)
   VTYPE, allocatable :: STIFQ3V_A(:,:,:), STIFQ3V_B(:,:,:,:)
   VTYPE, allocatable :: LOADH_A(:),       LOADH_B(:,:)
   VTYPE, allocatable :: LOADV_A(:,:),     LOADV_B(:,:,:)
!
!..load vector for the enriched space
   VTYPE, allocatable :: bloadHV(:), gram(:,:)
!
!..Stiffness matrices for the enriched test space
   VTYPE, allocatable :: stiffHV_Ht(:,:)
   VTYPE, allocatable :: stiffHV_Vt(:,:)
   VTYPE, allocatable :: stiffHV_Qt(:,:)
   VTYPE, allocatable :: stiff_ALL(:,:)
   VTYPE, allocatable :: zaloc(:,:)
!
   integer :: k1,k2,k3,m1,m2,a,b,sa,sb,sc,m1t,m2t
   integer :: l1,l2,l3,px,py,pz,nintx,ninty,nintz,nord1,nord2,nord3, &
              i1,i2,i3,j1,j2,j3,ltrial,kk,l,iflag
   integer :: alph,beta,idxa,idxb,idxalph,idxbeta,idxb2,idxb3,idxa2,idxa3
   integer :: ia1, ia2, ia3,ik, if, nsign
   integer :: idxa1, idxb1, n1, n2, ii, info, i
!..element mdle node dof   
   integer :: ndofHmdl,ndofEmdl,ndofVmdl,ndofQmdl, nrdofE
   integer :: ndofHHmdl,ndofEEmdl,ndofVVmdl,ndofQQmdl
   integer :: nrdofHint, nrdofVint, nvoid
!   
   real*8  :: weighthh, weightvv, wt123

   real*8  :: xi1,xi2,xi3,wt1,wt2,wt3,clock1,clock2,tmp,tmpt
   real*8, dimension(MAXPP+1) :: xilocx,xilocy,xilocz
   real*8, dimension(MAXPP+1) :: wlocx,wlocy,wlocz
   real*8, dimension(3,MAXNINT3ADD) :: wloc3
   real*8, dimension(3) :: xip,dHdx,dHHdx,shapVVt
   real*8, dimension(3,3) :: D,C
   real*8, dimension(MAXPP+1,2) :: shapH1,shapH2,shapH3
   real*8, dimension(MAXPP+1,MAXPP+1) :: sH2p,sH3p,dsH2p,dsH3p
   integer, dimension(3,3) :: deltak
!
!..DoF maps
   integer, allocatable  :: idxHH(:), idxVV(:)
! 
!..pml
   VTYPE :: dfpml(3), Jstretch(3,3), JJStretch(3,3), detJ, zA(3)
!..Acoustics specific
   real*8 :: h_elem, omeg
   real*8 :: tm,tm2,tm3
   integer(kind=8) :: t1,t2,clock_rate,clock_max
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
!..set the enriched order of approximation
   select case(ntype)
      case(MDLB)     ; nordP = NODES(Mdle)%order + Nord_add_local*111
      case(MDLN,MDLD); nordP = NODES(Mdle)%order + Nord_add_local*1
      case(MDLP)     ; nordP = NODES(Mdle)%order + Nord_add_local*11
   end select
!
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
!..init ALOC,BLOC
   do i=1,3
      do ii = 1,3
         ALOC(i,ii)%array = ZERO
      enddo
      BLOC(i)%array = ZERO
   enddo
!
!..allocate matrices
   allocate(bloadHV(NrTest))
   allocate(gram(NrTest,NrTest))
   allocate(stiffHV_Ht(nrdofHi,NrTest))
   allocate(stiffHV_Vt(nrdofVi,NrTest))
   allocate(stiffHV_Qt(4*nrdofQ,NrTest))
!                                                                    
   bloadHV=ZERO; gram=ZERO
   stiffHV_Ht=ZERO; stiffHV_Vt = ZERO; stiffHV_Qt = ZERO;
!
!..adjusted frequency for the test space 
   ! omeg = min(OMEGA,6.d0/h_elem)
   omeg = OMEGA
!
!..allocate auxiliary arrays for sum-factorization
   allocate(AUXHH_A(3,3,nrdofH3**2))
   allocate(AUXHH_A0(nrdofH3**2))
   allocate(AUXVH_A(3,nrdofH3**2))
   allocate(AUXVH_A0(3,nrdofH3**2))
   allocate(AUXVV_A(3,3,nrdofH3**2))
   allocate(AUXVV_A0(3,3,nrdofH3**2))
!
   allocate(AUXHH_B(3,3,nrdofH2**2,nrdofH3**2))
   allocate(AUXHH_B0(nrdofH2**2,nrdofH3**2))
   allocate(AUXVH_B (3,nrdofH2**2,nrdofH3**2))
   allocate(AUXVH_B0(3,nrdofH2**2,nrdofH3**2))
   allocate(AUXVV_B (3,3,nrdofH2**2,nrdofH3**2))
   allocate(AUXVV_B0(3,3,nrdofH2**2,nrdofH3**2))
!
   allocate(STIFQH_A(nrdofQ3_tr*nrdofH3))
   allocate(STIFQH_B(nrdofQ3_tr*nrdofH3,nrdofQ3_tr*nrdofH3))
   allocate(STIFQ3H_A(3,3,nrdofQ3_tr*nrdofH3))
   allocate(STIFQ3H_B(3,3,nrdofQ2_tr*nrdofH2,nrdofQ3_tr*nrdofH3))
   allocate(STIFQV_A(3,nrdofQ3_tr*nrdofH3))
   allocate(STIFQV_B(3,nrdofQ2_tr*nrdofH2,nrdofQ3_tr*nrdofH3))
   allocate(STIFQ3V_A(3,3,nrdofQ3_tr*nrdofH3))
   allocate(STIFQ3V_B(3,3,nrdofQ2_tr*nrdofH2,nrdofQ3_tr*nrdofH3))
!
   allocate(LOADH_A(nrdofH3))
   allocate(LOADH_B(nrdofH2,nrdofH3))
   allocate(LOADV_A(3,nrdofH3))
   allocate(LOADV_B(3,nrdofH2,nrdofH3))
!
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
   INTEGRATION = Nord_add_local
   call set_3Dint_fi(ntype,norder,nord1,nord2,nord3,nintx,ninty,nintz,xiloc,wloc3)
   INTEGRATION = 0
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
      STIFQH_B = ZERO ;  STIFQ3H_B = ZERO ;  STIFQ3V_B = ZERO
      STIFQV_B = ZERO ;  LOADH_B   = ZERO ;  LOADV_B   = ZERO
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
         STIFQH_A = ZERO ;  STIFQ3H_A = ZERO ; STIFQ3V_A = ZERO
         STIFQV_A = ZERO ;  LOADH_A   = ZERO ; LOADV_A   = ZERO
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
            call shape3DH(ntype,xip,norder,norient_edge,norient_face, nvoid,shapH,gradH)
!        ...Geometry map
            call geom3D(Mdle,xip,xnod,shapH,gradH,nrdofH,x,dxdxi,dxidx,rjac,iflag)
            if (iflag.ne.0) then
                  write(*,5999) Mdle,rjac
 5999             format('Negative JacobiancMdle,rjac=',i8,2x,e12.5)
               stop
            endif
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
!        ...Evaluate forcing function f at physical point x
            call getf( Mdle , x , zf )
!
!        ...get pml function derivatives at physical point x
            dfpml = ZERO
            Jstretch = ZERO
            JJStretch= ZERO
            detJ = ZERO
            call getpml(x,dfpml)
            Jstretch(1,1) = dfpml(1)
            Jstretch(2,2) = dfpml(2)
            Jstretch(3,3) = dfpml(3)

            call zgemm('N', 'N', 3, 3, 3, ZONE, JStretch, 3, JStretch, 3, ZERO, JJStretch, 3)
            detJ = dfpml(1)*dfpml(2)*dfpml(3)
            
            zA(1) = conjg(JJStretch(1,1)/detJ)
            zA(2) = conjg(JJStretch(2,2)/detJ)
            zA(3) = conjg(JJStretch(3,3)/detJ)
! 
            JJStretch = JJStretch/detJ
!            
!        ...HERE STARTS COMPUTATION OF VOLUME INTEGRALS: gram, STIFFNESS AND LOAD MATRICES
!        ...COMPUTATION OF AUXILIARY ARRAYS 'A' FOR gram MATRIX - H1,H1 BLOCK
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
                               ! + (omeg**2+ALPHA)*shapH3(i3,1)*shapH3(j3,1)*weighthh
                               + ( omeg**2*detJ*conjg(detJ)*shapH3(i3,1)*shapH3(j3,1)  &
                               +   ALPHA*shapH3(i3,1)*shapH3(j3,1))*weighthh                               
!    
!              ...end of j3 loop for H1 trial functions
               enddo
            enddo
! 
!        ...Integration of enriched stiffness matrix and load vector - H1 test func
!        ...start new loop on i3 - H1 test functions
            do i3=1,nrdofH3
!           ...load vector inner H1 integral               
               LOADH_A(i3)=LOADH_A(i3)+zf(1)*shapH3(i3,1)*weighthh
!           ...start new loop on j3 - L2 trial functions
               do j3=1,nrdofQ3_tr
!              ...combine indices i3, j3 into k3
                  k3=(i3-1)*nrdofQ3_tr+j3
                  STIFQH_A(k3) = STIFQH_A(k3) &
                               ! + ZIMG*OMEGA*shapH3(j3+1,2)*shapH3(i3,1)*wt123
                               + ZIMG*OMEGA*detJ*shapH3(j3+1,2)*shapH3(i3,1)*wt123

                  do b=1,3
                     do a=1,3
                        sa=1+deltak(a,3)
                        STIFQ3H_A(a,b,k3) = STIFQ3H_A(a,b,k3) &
                                          - (dxidx(a,b)*shapH3(j3+1,2))*shapH3(i3,sa)
                     enddo
                  enddo
               enddo
            enddo
!
!        ...COMPUTATION OF AUXILIARY ARRAYS 'A' FOR gram MATRIX - H(div),H1 BLOCK
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
                                       ! - shapH3(i3,1 )*shapH3(idxb,2 )*ZIMG*omeg*wt123
                                       - conjg(detJ)*shapH3(i3,1 )*shapH3(idxb,2 )*ZIMG*omeg*wt123
                        sa=1+deltak(b,3)                                              
                        AUXVH_A (b,k3) = AUXVH_A (b,k3)                              &
                                       ! + shapH3(i3,sa)*shapH3(idxb,sb)*ZIMG*omeg*wt123
                                       + conjg(zA(b))*shapH3(i3,sa)*shapH3(idxb,sb)*ZIMG*omeg*wt123
                     endif
                  enddo
!              ...end of j3 loop for H(div) trial functions
               enddo
!           ...end of i3 loop for H1 test functions
            enddo
!
!        ...COMPUTATION OF AUXILIARY ARRAYS 'A' FOR gram MATRIX - H(div),H(div) BLOCK
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
                              ! AUXVV_A (b,a,k3) = AUXVV_A (b,a,k3)+(omeg**2+ALPHA)   &
                                               ! * shapH3(idxa,sa)*shapH3(idxb,sb)*C(a,b)
                              AUXVV_A (b,a,k3) = AUXVV_A (b,a,k3)   &
                                               + omeg**2   &
                                               * zA(a)*conjg(zA(b))*shapH3(idxa,sa)*shapH3(idxb,sb)*C(a,b) &
                                               + ALPHA*shapH3(idxa,sa)*shapH3(idxb,sb)*C(a,b)
!
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
!        ...Integration of enriched load vector - H(div) test f
            do i3=1,nrdofH3
               do a=1,3
                  sa=1+1-deltak(a,3)
                  ia3=i3+1-deltak(a,3)
                  if (ia3.le.nrdofH3) then
                     LOADV_A(a,i3) = LOADV_A(a,i3)      &
                                   + (zf(2)*dxdxi(1,a)  &
                                   +  zf(3)*dxdxi(2,a)  &
                                   +  zf(4)*dxdxi(3,a)) &
                                   * shapH3(ia3,sa)*rjac
                  endif
               enddo                 
            enddo   
!        ...Integration of enriched stiffness matrix - H(div) test f
            do i3=1,nrdofH3
               do j3=1,nrdofQ3_tr
                  k3=(i3-1)*nrdofQ3_tr+j3
                  do a=1,3
                     sa=1+1-deltak(a,3)
                     ia3=i3+1-deltak(a,3)
                     if (ia3.le.nrdofH3) then
                        do b=1,3                        
                           STIFQ3V_A(b,a,k3) = STIFQ3V_A(b,a,k3) &
!                                             + ZIMG*OMEGA*dxdxi(b,a) &
!                                             * shapH3(j3+1,2)*shapH3(ia3,sa)
                                             + ZIMG*OMEGA*JJStretch(a,a)*dxdxi(b,a) &
                                             * shapH3(j3+1,2)*shapH3(ia3,sa)

                        enddo                     
                        STIFQV_A(a,k3) = STIFQV_A(a,k3) &
                                       - shapH3(j3+1,2)*shapH3(ia3,2)*weightvv
                     endif
                  enddo
               enddo
            enddo
!        ...pz loop ends
         enddo         
!
!     ...COMPUTATION OF AUXILIARY ARRAYS 'B' FOR gram MATRIX - H1,H1 BLOCK
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

!     ...load vector H1 second integral test functions
         do i3=1,nrdofH3
            do i2=1,nrdofH2
               LOADH_B(i2,i3)=LOADH_B(i2,i3)+shapH2(i2,1)*LOADH_A(i3)
            enddo
         enddo
         do i3=1,nrdofH3
            do j3=1,nrdofQ3_tr
               k3=(i3-1)*nrdofQ3_tr+j3
               do i2=1,nrdofH2
                  do j2=1,nrdofQ2_tr
                     k2=(i2-1)*nrdofQ2_tr+j2
                     STIFQH_B(k2,k3) = STIFQH_B(k2,k3) &
                                     + shapH2(j2+1,2)*shapH2(i2,1)*STIFQH_A(k3)
                     do b=1,3
                        do a=1,3
                           sa=1+deltak(a,2)
                           STIFQ3H_B(a,b,k2,k3) = STIFQ3H_B(a,b,k2,k3) &
                                                + shapH2(j2+1,2)*shapH2(i2,sa) &
                                                * STIFQ3H_A(a,b,k3)
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo 
!
!     ...COMPUTATION OF AUXILIARY ARRAYS 'B' FOR gram MATRIX - H(div),H1 BLOCK
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
!     ...COMPUTATION OF AUXILIARY ARRAYS 'B' FOR gram MATRIX - H(div),H(div) BLOCK
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
! 
!     ...load vector for H(div) middle integral
         do i3=1,nrdofH3
            do i2=1,nrdofH2
               do a=1,3   
                  ia2=i2+1-deltak(a,2)
                  if (ia2.le.nrdofH2) then
                      LOADV_B(a,i2,i3) = LOADV_B(a,i2,i3) &
                                       + shapH2(ia2,sa)*LOADV_A(a,i3)
                  endif                     
               enddo
            enddo
         enddo                          
! 
         do i3=1,nrdofH3
            do j3=1,nrdofQ3_tr
               k3=(i3-1)*nrdofQ3_tr+j3
               do i2=1,nrdofH2
                  do j2=1,nrdofQ2_tr
                     k2=(i2-1)*nrdofQ2_tr+j2
                     do a=1,3
                        sa=1+1-deltak(a,2)
                        ia2=i2+1-deltak(a,2)
                        if (ia2.le.nrdofH2) then
                           do b=1,3                        
                              STIFQ3V_B(b,a,k2,k3) = STIFQ3V_B(b,a,k2,k3) &
                                                   + shapH2(j2+1,2)*shapH2(ia2,sa) &
                                                   * STIFQ3V_A(b,a,k3)
                           enddo                     
                           STIFQV_B(a,k2,k3) = STIFQV_B(a,k2,k3) &
                                             + shapH2(j2+1,2)*shapH2(ia2,2) & 
                                             * STIFQV_A(a,k3)
                        endif
                     enddo
                  enddo
               enddo
            enddo
         enddo
!
!     ...py loop ends
      enddo
!  
!  ...FINAL COMPUTATION OF gram MATRIX - H1,H1 BLOCK
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
!                    ...Only upper-triangular part of gram matrix computed
                        if (m2.ge.m1) then
!                       ...loops over components of gradient of functions J and I
                           do b=1,3
                              sb=1+deltak(b,1)                            
                              do a=1,3
                                 sa=1+deltak(a,1)
                                 gram(m1,m2) = gram(m1,m2)       &
                                             + shapH1(i1,sa)*shapH1(j1,sb) &
                                             * AUXHH_B(a,b,k2,k3)
                              enddo
                           enddo                  
                           gram(m1,m2) = gram(m1,m2)       &
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

      do i3=1,nrdofH3
         do i2=1,nrdofH2
            do i1=1,nrdofH1
               m1=i1+(i2-1)*nrdofH1+(i3-1)*nrdofH1*nrdofH2
!           ...load vector H1 test functions               
               bloadHV(m1)=bloadHV(m1)+shapH1(i1,1)*LOADH_B(i2,i3)
               do j3=1,nrdofQ3_tr
                  k3=(i3-1)*nrdofQ3_tr+j3
                  do j2=1,nrdofQ2_tr
                     k2=(i2-1)*nrdofQ2_tr+j2
                     do j1=1,nrdofQ1_tr
                        m2=1+4*(j1-1+(j2-1)*nrdofQ1_tr+(j3-1)*nrdofQ1_tr*nrdofQ2_tr)
                        sTIFFHV_Qt(m2,m1) = sTIFFHV_Qt(m2,m1) &
                                         + shapH1(j1+1,2)*shapH1(i1,1) &
                                         * STIFQH_B(k2,k3)
                        do b=1,3
                           do a=1,3
                              sa=1+deltak(a,1)
                              sTIFFHV_Qt(m2+b,m1) = sTIFFHV_Qt(m2+b,m1) &
                                                 + shapH1(j1+1,2)*shapH1(i1,sa) &
                                                 * STIFQ3H_B(a,b,k2,k3)
                           enddo
                        enddo                        
                     enddo
                  enddo
               enddo
            enddo
         enddo 
      enddo
!
!  ...FINAL COMPUTATION OF gram MATRIX - H(div),H1 BLOCK
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
                                    gram(m1,n2) = gram(m1,n2)              &
                                                + shapH1(i1,1 )*shapH1(idxb1,2 )*AUXVH_B0(b,k2,k3)
                                    gram(m1,n2) = gram(m1,n2)              &
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
!  ...FINAL COMPUTATION OF gram MATRIX - H(div), H(div) BLOCK
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
!                                            ...fill only the upper part of the gram matrix block    
                                                if (n2.ge.n1) then
                                                   gram(n1,n2) = gram(n1,n2) &
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
!      
      do a=1,3
         sa=1+1-deltak(a,1)                      
         do i3=1,nrdofH3
            ia3=i3+1-deltak(a,3)
            if (ia3.le.nrdofH3) then
               do i2=1,nrdofH2
                  ia2=i2+1-deltak(a,2)
                  if (ia2.le.nrdofH2) then
                     do i1=1,nrdofH1
                        ia1=i1+1-deltak(a,1)
                        if (ia1.le.nrdofH1) then
                           select case(a)                        
                           case(1)
                              m1 = i1+nrdofH1*(i2-1)+nrdofH1*nord2*(i3-1)
                           case(2)
                              m1 = nrdofH1*nord2*nord3                        &
                                 + i1+nord1*(i2-1)+nord1*nrdofH2*(i3-1)
                           case(3)
                              m1 = nrdofH1*nord2*nord3+nord1*nrdofH2*nord3    &
                                 + i1+nord1*(i2-1)+nord1*nord2*(i3-1)
                           end select
                           n1 = nrdofHH+m1
                           do j3=1,nrdofQ3_tr
                              k3=(i3-1)*nrdofQ3_tr+j3
                              do j2=1,nrdofQ2_tr
                                 k2=(i2-1)*nrdofQ2_tr+j2
                                 do j1=1,nrdofQ1_tr
                                    m2=1+4*(j1-1+(j2-1)*nrdofQ1_tr+(j3-1)*nrdofQ1_tr*nrdofQ2_tr)
                                    sTIFFHV_Qt(m2,n1) = sTIFFHV_Qt(m2,n1) &
                                                     + shapH1(j1+1,2)*shapH1(ia1,2) &
                                                     * STIFQV_B(a,k2,k3)
                                    do b=1,3
                                       sTIFFHV_Qt(m2+b,n1) = sTIFFHV_Qt(m2+b,n1) &
                                                          + shapH1(j1+1,2)*shapH1(ia1,sa) &
                                                          * STIFQ3V_B(b,a,k2,k3)
                                    enddo
                                 enddo
                              enddo
                           enddo
                           bloadHV(n1) = bloadHV(n1) + shapH1(ia1,sa)*LOADV_B(a,i2,i3)
                        endif
                     enddo
                  endif   
               enddo
            endif   
         enddo
      enddo
!  ...px loop ends
   enddo
!
!..find ndof associated with the mdle node of the element

   call compute_enriched_order(nordP, norderP)
   call ndof_nod(ntype,norder(nre+nrf+1) ,ndofHmdl,ndofEmdl,ndofVmdl,ndofQmdl)
   call ndof_nod(ntype,norderP(nre+nrf+1),ndofHHmdl,ndofEEmdl,ndofVVmdl,ndofQQmdl)

   nrdofHint = nrdofH-ndofHmdl
   nrdofVint = nrdofV-ndofVmdl
! 
   allocate(idxHH(nrdofHH-ndofHHmdl))
   ik = 0 
   do i3 = 1,nrdofH3
      do i2 = 1,nrdofH2
         do i1 = 1,nrdofH1
            if (i3 .lt.3 .or. i2 .lt.3 .or. i1 .lt. 3) then
               m1=i1+(i2-1)*nrdofH1+(i3-1)*nrdofH2*nrdofH1
               ik = ik +1
               idxHH(ik) = m1
            endif   
         enddo
      enddo
   enddo
!
   allocate(idxVV(NrdofVV-ndofVVmdl))
!
   ik = 0
   do i3 = 1,nrdofH3-1
      do i2 = 1,nrdofH2-1
         do i1 = 1,2
            ik = ik +1
            m1 = i1+nrdofH1*(i2-1) + nrdofH1*nord2*(i3-1)
            idxVV(ik) = m1
         enddo
      enddo
   enddo
   do i3 = 1,nrdofH3-1
      do i2 = 1,2
         do i1 = 1,nrdofH1-1
            ik = ik +1
            m1 = nrdofH1*nord2*nord3 + i1+nord1*(i2-1) + nord1*nrdofH2*(i3-1)
            idxVV(ik) = m1
         enddo
      enddo
   enddo
   do i3 = 1,2
      do i2 = 1,nrdofH2-1
         do i1 = 1,nrdofH1-1
            ik = ik +1
            m1 = nrdofH1*nord2*nord3 + nord1*nrdofH2*nord3 + i1+nord1*(i2-1)+nord1*nord2*(i3-1)
            idxVV(ik) = m1
         enddo
      enddo
   enddo
!
!..deallocate sum-factorization arrays
   deallocate(AUXHH_A, AUXHH_A0)
   deallocate(AUXVH_A, AUXVH_A0)
   deallocate(AUXVV_A, AUXVV_A0)
!
   deallocate(AUXHH_B, AUXHH_B0)
   deallocate(AUXVH_B, AUXVH_B0)
   deallocate(AUXVV_B, AUXVV_B0)
!
   deallocate(STIFQH_A,STIFQH_B)
   deallocate(STIFQ3H_A,STIFQ3H_B)
   deallocate(STIFQV_A,STIFQV_B)
   deallocate(STIFQ3V_A,STIFQ3V_B)
   deallocate(LOADH_A,LOADH_B)
   deallocate(LOADV_A,LOADV_B)
!
!
!----------------------------------------------------------------------
!     B O U N D A R Y    I N T E G R A L S                           
!----------------------------------------------------------------------
!
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
      INTEGRATION = Nord_add_local
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
         call shape3DH(ntype,xi,norder,norient_edge,norient_face, nvoid,shapH,gradH)
! 
!     ...determine element Hdiv shape functions (for fluxes)
         call shape3DV(ntype,xi,norder,norient_face, nvoid,shapV,divV)
!  
!     ...determine discontinuous H1 shape functions
         call shape3HH(ntype,xi,nordP, nvoid,shapHH,gradHH)
!         
!     ...determine discontinuous H(div) shape functions
         call shape3VV(ntype,xi,nordP, nvoid,shapVV,divVV)
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
            do ii = 1,nrdofHH-ndofHHmdl
               k1 = idxHH(ii)
!           ...value of the shape function at the point
               q = shapHH(k1)
!           ...accumulate for the load vector
               bloadHV(k1) = bloadHV(k1) + q*zg*weight
!           ...loop through H1 trial functions
               do k2 = 1,nrdofH-ndofHmdl
!              ...value of the shape function at the point
                  p = shapH(k2)
!              ...accumulate for the Stiffness matrix
                  stiffHV_Ht(k2,k1) = stiffHV_Ht(k2,k1) + q*p*weight
               enddo
            enddo                                                
!        ...regular boundary
         else
! 
!        ...loop through enriched H1 test functions
            do ii = 1,nrdofHH-ndofHHmdl
               k1 = idxHH(ii)
!               
               q  = shapHH(k1)
!
!           ...loop through H(div) trial functions
               do k2=1,nrdofV-ndofVmdl
! 
!              ...normal component (Piola transformation)
                  vec(1:3) = dxdxi(1:3,1)*shapV(1,k2)   &
                           + dxdxi(1:3,2)*shapV(2,k2)   & 
                           + dxdxi(1:3,3)*shapV(3,k2)
                  vec(1:3) = vec(1:3)/rjac
                  un = vec(1)*rn(1)+vec(2)*rn(2)+vec(3)*rn(3)
!
!              ...accumulate for the Stiffness matrix
                  stiffHV_Vt(k2,k1) = stiffHV_Vt(k2,k1) + q*un*weight
!
!              ...end of loop through H(div) trial functions
               enddo
!           ...end of loop through H1 test functions
            enddo
         endif
!         
!     ...loop through H(div) enriched test functions
         do ii = 1,NrdofVV-ndofVVmdl
            k1 = idxVV(ii)
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
            do k2=1,nrdofH-ndofHmdl
!
!           ...value of the shape function at the point
               p = shapH(k2)
!
!           ...accumulate for the rectangular Stiffness matrix
               stiffHV_Ht(k2,n1) = stiffHV_Ht(k2,n1) + vn*p*weight
!           ...end of loop through H1 trial functions
            enddo
!        ...end of loop through H(div) test functions                 
         enddo
!     ...end of loop through integration points
      enddo
!  ...end of loop through faces  
   enddo
!
   deallocate(idxHH)
   deallocate(idxVV)
!
!----------------------------------------------------------------------
!      Alternative construction of normal matrix
!----------------------------------------------------------------------
! 
   i1 = NrTEST ; j1 = NrdofHi ; j2 = NrdofVi ; j3 = 4*NrdofQ
!
   allocate(stiff_ALL(NrTest,NrTrial+1)); stiff_all = ZERO
!
   stiff_ALL(1:i1,1:j1)             = transpose(stiffHV_Ht(1:j1,1:i1))
   stiff_ALL(1:i1,j1+1:j1+j2)       = transpose(stiffHV_Vt(1:j2,1:i1))
   stiff_ALL(1:i1,j1+j2+1:j1+j2+j3) = transpose(stiffHV_Qt(1:j3,1:i1))
   stiff_ALL(1:i1,j1+j2+j3+1)       = bloadHV(1:i1)
!
   deallocate(stiffHV_Ht,stiffHV_Vt,stiffHV_Qt,bloadHV)
!
!    call diag_scaling(N,NRHS,gram,stiff_ALL)
!
!..diagonal scaling of the gram matrix
   uplo = 'U'
! 
!..factorize the test Stiffness matrix
   call ZPOTRF(uplo,NrTEST,gram,NrTEST,info)
   if (info.ne.0) then
      write(*,*) 'elem_fi: gram ZPOTRF: Mdle,info = ',Mdle,info
      stop 1
   endif
!
   call ZTRSM('L',uplo,'C','N',NrTEST,NrTrial+1,ZONE,gram,NrTEST,stiff_ALL,NrTEST)
!
   deallocate(gram)
   allocate(zaloc(NrTrial+1,NrTrial+1)); ZALOC = ZERO
!
   call ZHERK('U','C',NrTrial+1,NrTEST,ZONE,stiff_ALL,NrTEST,ZERO,zaloc,NrTrial+1)
!
   deallocate(stiff_ALL)
!
   do i=1,NrTrial
      zaloc(i+1:NrTrial+1,i) = conjg(zaloc(i,i+1:NrTrial+1))
   enddo
!
!..E. Fill ALOC and BLOC matrices
   BLOC(1)%array(1:j1,1) = zaloc(1:j1,j1+j2+j3+1)
   BLOC(2)%array(1:j2,1) = zaloc(j1+1:j1+j2,j1+j2+j3+1) 
   BLOC(3)%array(1:j3,1) = zaloc(j1+j2+1:j1+j2+j3,j1+j2+j3+1) 
!
   ALOC(1,1)%array(1:j1,1:j1) = zaloc(1:j1,1:j1)
   ALOC(1,2)%array(1:j1,1:j2) = zaloc(1:j1,j1+1:j1+j2)
   ALOC(1,3)%array(1:j1,1:j3) = zaloc(1:j1,j1+j2+1:j1+j2+j3)
!
   ALOC(2,1)%array(1:j2,1:j1) = zaloc(j1+1:j1+j2,1:j1)
   ALOC(2,2)%array(1:j2,1:j2) = zaloc(j1+1:j1+j2,j1+1:j1+j2)
   ALOC(2,3)%array(1:j2,1:j3) = zaloc(j1+1:j1+j2,j1+j2+1:j1+j2+j3)
!
   ALOC(3,1)%array(1:j3,1:j1) = zaloc(j1+j2+1:j1+j2+j3,1:j1)
   ALOC(3,2)%array(1:j3,1:j2) = zaloc(j1+j2+1:j1+j2+j3,j1+1:j1+j2)
   ALOC(3,3)%array(1:j3,1:j3) = zaloc(j1+j2+1:j1+j2+j3,j1+j2+1:j1+j2+j3)
!
   deallocate(zaloc)
!
!   
   end subroutine elem_DPG_UWEAK_ACOUSTICS_fi
