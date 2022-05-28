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
      character(len=4) :: etype
   !..element order, orientation for edges and faces
      integer :: norder(19),norderP(19), nordP
      integer :: nordx,nordy,nordz
      integer :: nrdofH,nrdofE,nrdofV,nrdofQ
      integer :: nrdofEi
      integer :: nrTest, nrTrial
      integer :: nrdofHH,nrdofEE,nrdofVV,nrdofQQ
      integer :: ndofHmdl,ndofEmdl,ndofVmdl,ndofQmdl
      integer :: nord_add_local
      integer :: nrv, nre, nrf
   !..sum factorization
      integer :: nrdofHx, nrdofHy, nrdofHz   
      integer :: nrdofHx_tr, nrdofHy_tr, nrdofHz_tr
      integer :: nrdofQx_tr, nrdofQy_tr, nrdofQz_tr   
   !----------------------------------------------------------------------
   !
      etype = NODES(Mdle)%type
      nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)   
   ! 
   !----------------------------------------------------------------------
   !
      Itest (1:NR_PHYSA) = 0
      Itrial(1:NR_PHYSA) = 0
   
      select case(NODES(Mdle)%case)
      case(3)
         Itest(1:NR_PHYSA)=1; Itrial(1:NR_PHYSA)=1
   
   !  ...initialize nord_add_local
         nord_add_local = NORD_ADD
   !      
   !  ...determine order of approximation
         call find_order(Mdle, norder)
   !
    10   continue
   !  ...set the enriched order of approximation
         select case(etype)
            case('mdlb')        ; nordP = NODES(Mdle)%order + nord_add_local*111
            case('mdln','mdld') ; nordP = NODES(Mdle)%order + nord_add_local*1
            case('mdlp')        ; nordP = NODES(Mdle)%order + nord_add_local*11
         end select
   !
         call compute_enriched_order(nordP, norderP)
   !
         INTEGRATION = nord_add_local
         call compute_1D_ord(etype,norder(19),nordx,nordy,nordz) 
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
         call celndof(etype,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
   !  ...total number of test dof associated with the element
         call celndof(etype,norderP,nrdofHH,nrdofEE,nrdofVV,nrdofQQ)
   !
         call ndof_nod(etype,norder(nre+nrf+1),ndofHmdl,ndofEmdl,ndofVmdl,ndofQmdl)
         nrTest = 2*nrdofEE
         nrdofEi = nrdofE-ndofEmdl
         nrTrial = 2*nrdofEi + 6*nrdofQ
   !
   
         call elem_maxwell_fi(Mdle,nord_add_local,nrTest,nrTrial,   &
                              nrdofEE,    nrdofE,     nrdofQ,       &
                              nrdofEi,                              &
                              nrdofHx,    nrdofHy,    nrdofHz,      &
                              nrdofQx_tr, nrdofQy_tr, nrdofQz_tr)    
   
   
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
!     routine name      - elem_maxwell_fi
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
   subroutine elem_maxwell_fi(Mdle,Nord_add_local,NrTest,NrTrial,   &
                              NrdofEE,    NrdofE,     NrdofQ,       &
                              NrdofEi,                              &
                              NrdofH1,    NrdofH2,    NrdofH3,      &
                              NrdofQ1_tr, NrdofQ2_tr, NrdofQ3_tr)
!
   use data_structure3D
   use control,     only: INTEGRATION, ISTC_FLAG
   use parametersDPG
   use physics,     only: NR_PHYSA
   use assembly,    only: ALOC,BLOC,NR_RHS
   use common_prob_data
!
!----------------------------------------------------------------------
!
   implicit none

   integer, intent(in) :: Mdle, Nord_add_local 
   integer, intent(in) :: NrTest, NrTrial
   integer, intent(in) :: NrdofEE, NrdofE, NrdofQ, NrdofEi
   integer, intent(in) :: NrdofH1, NrdofH2, NrdofH3
   integer, intent(in) :: NrdofQ1_tr, NrdofQ2_tr, NrdofQ3_tr

   character(len=4) :: etype,ftype
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
!..element interface test dof
   integer :: nrdofEEi
!..element mdle node dof   
   integer :: ndofHmdl,ndofEmdl,ndofVmdl,ndofQmdl
   integer :: ndofHHmdl,ndofEEmdl,ndofVVmdl,ndofQQmdl
!
!..geometry dof
   real*8  :: xnod(3,MAXbrickH)
! 
!..geometry
   real*8  :: xi(3),dxidt(3,2),x(3),dxdxi(3,3),dxidx(3,3), dxdt(3,2)
   real*8  :: rt(3,2),rn(3),t(2)
!
!..H1 shape functions (for Geometry)
   real*8  :: shapH (MAXbrickH)  , gradH(3,MAXbrickH)
!..H(curl) shape functions
   real*8  :: shapE (3,MAXbrickE) , curlE(3,MAXbrickE) 
   real*8  :: shapEE(3,MAXbrickEE), curlEE(3,MAXbrickEE)

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

 !..counters and other integers

   integer :: l, if, nsign, iflag   
   integer :: nrdof,nrdofH
!
!..source
   complex*16 :: zf(6), za1, za2, zb1, zb2, zg(3)
!
!..sum factorization
!
   integer    :: nord1, nord2, nord3, nintx,ninty,nintz
   integer    :: px,py,pz
   real*8     :: wloc3(3,MAXNINT3ADD)
   real*8     :: wlocx(MAXPP+1), xilocx(MAXPP+1)
   real*8     :: wlocy(MAXPP+1), xilocy(MAXPP+1)
   real*8     :: wlocz(MAXPP+1), xilocz(MAXPP+1)
   real*8     :: xi1, xi2, xi3, wt1, wt2, wt3, xip(3)
   real*8     :: weighthh, weightvv, wt123
   real*8     :: sH2p(MAXPP+1,MAXPP+1), dsH2p(MAXPP+1,MAXPP+1)
   real*8     :: sH3p(MAXPP+1,MAXPP+1), dsH3p(MAXPP+1,MAXPP+1)
!
   real*8     :: shapH1(MAXPP+1,2), shapH2(MAXPP+1,2), shapH3(MAXPP+1,2)
!
   complex*16, allocatable :: auxFF_A(:,:,:),     auxFF_B(:,:,:,:)
   complex*16, allocatable :: auxCC_A(:,:,:,:,:), auxCC_B(:,:,:,:,:,:)
   complex*16, allocatable :: auxFC_A(:,:,:,:),   auxFC_B(:,:,:,:,:)
   complex*16, allocatable :: auxCF_A(:,:,:,:),   auxCF_B(:,:,:,:,:)
!
   complex*16, allocatable :: stiffQF_A(:,:,:)      
   complex*16, allocatable :: stiffQF_B(:,:,:,:)
   complex*16, allocatable :: stiffQC_A(:,:,:,:)
   complex*16, allocatable :: stiffQC_B(:,:,:,:,:)
   complex*16, allocatable :: load1_A(:,:), load2_A(:,:)
   complex*16, allocatable :: load1_B(:,:,:), load2_B(:,:,:)

   real*8     :: D(3,3),C(3,3)
   integer    :: deltak(3,3), sa, sb
   integer    :: nvoid
   integer    :: i1,i2,i3,j1,j2,j3,k1,k2,k3,k,kk,l1,l2,l3,m1,m2,nk
   integer    :: a,b,alph,beta,idxa,idxb,idxalph,idxbeta,idxb2,idxb3,idxa2,idxa3
!
!..load vector for the enriched space
   complex*16 :: bloadE(NrTest)
!
!..gram matrix in packed format
   complex*16 :: gramP(NrTest*(NrTest+1)/2)
!
!..stiffness matrices for the enriched test space
   complex*16 :: stiffE_Et(2*NrdofEi,NrTest)
   complex*16 :: stiffE_Qt(6*NrdofQ ,NrTest)
   complex*16 :: stiff_ALL(NrTest   ,NrTrial+1)
   complex*16 :: zaloc    (NrTrial+1,NrTrial+1)
!
!..workspace for trial and test variables
   real*8  :: fldE(3), fldH(3), crlE(3), crlH(3), rnE(3), rnH(3), rn2E(3)
   real*8  :: fldF(3), fldG(3), crlF(3), crlG(3), rnHG, rnEF, rn2EG
   integer :: n,m,NRHS,I,info
!..Maxwell specific
   real*8    :: omeg, h_elem
!
!..H(curl) bubble index
   integer, allocatable  :: idxEE(:)
   integer :: ik
!..timers
   !real*8 :: omp_get_wtime, start, finish
!   
!..lapack
   character*1 uplo
   nk(k1,k2) = (k2-1)*k2/2+k1
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
   etype = NODES(Mdle)%type
   nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
!
!..determine order of approximation
   call find_order(Mdle, norder)

!..copy norder and replace mdle node order with p=1 
!..(Avoids computing with bubbles on the boundary)   
   norderi(1:nre+nrf) = norder(1:nre+nrf)
!
!..set the enriched order of approximation
   select case(etype)
   case('mdlb')
      nordP = NODES(Mdle)%order+Nord_add_local*111
      norderi(nre+nrf+1) = 111
   case('mdln','mdld')
      nordP = NODES(Mdle)%order+Nord_add_local
      norderi(nre+nrf+1) = 1
   case('mdlp')
      nordP = NODES(Mdle)%order+Nord_add_local*11
      norderi(nre+nrf+1) = 11
   end select
!   
!..determine edge and face orientations
   call find_orient(Mdle, norient_edge,norient_face)
! !                                                                     
!..determine nodes coordinates 
   call nodcor(Mdle, xnod)
!
!..determine element size and scale correction
   call find_hmin(Mdle,h_elem)
!..get the element boundary conditions flags
   call find_bc(Mdle, ibc)
!                                                                    
!..adjusted frequency for the test space 
   omeg = min(OMEGA,6.d0/h_elem)
   ! omeg = OMEGA
!
   za1 = ZIMG*OMEGA*EPS
   zb1 = ZIMG*omeg*EPS
   za2 = ZIMG*OMEGA*MU
   zb2 = ZIMG*omeg*MU
!
   allocate(auxFF_A(3,3,nrdofH3**2),auxFF_B(3,3,nrdofH2**2,nrdofH3**2))
   allocate(auxCC_A(2,2,3,3,nrdofH3**2),auxCC_B(2,2,3,3,nrdofH2**2,nrdofH3**2))
   allocate(auxFC_A(2,3,3,nrdofH3**2),auxFC_B(2,3,3,nrdofH2**2,nrdofH3**2))
   allocate(auxCF_A(2,3,3,nrdofH3**2),auxCF_B(2,3,3,nrdofH2**2,nrdofH3**2))
!
   allocate(stiffQF_A(3,3,nrdofQ3_tr*nrdofH3)      )
   allocate(stiffQF_B(3,3,nrdofQ2_tr*nrdofH2,nrdofQ3_tr*nrdofH3))
   allocate(stiffQC_A(2,3,3,nrdofQ3_tr*nrdofH3))
   allocate(stiffQC_B(2,3,3,nrdofQ2_tr*nrdofH2,nrdofQ3_tr*nrdofH3))
   allocate(load1_A(3,nrdofH3), load2_A(3,nrdofH3))
   allocate(load1_B(3,nrdofH2,nrdofH3), load2_B(3,nrdofH2,nrdofH3))
!
   bloadE=ZERO; zaloc = ZERO
   stiffE_Et=ZERO; stiffE_Qt = ZERO; stiff_ALL = ZERO
!
   gramP=ZERO
!
!
!----------------------------------------------------------------------
!     E L E M E N T    I N T E G R A L S                             
!----------------------------------------------------------------------
! 
!   
!..initialize 
   xiloc  = 0.0d0;  wloc3  = 0.0d0
   xilocx = 0.d0 ;  xilocy = 0.d0 ;  xilocz = 0.d0
   wlocx  = 0.d0 ;  wlocy  = 0.d0 ;  wlocz  = 0.d0        
   sa     = 0    ;  sb     = 0    ;  D      = ZERO
!
!..set up integration data
   INTEGRATION = Nord_add_local
   call set_3Dint_fi(etype,norder,nord1,nord2,nord3,nintx,ninty,nintz,xiloc,wloc3)
   INTEGRATION = 0
!
!..check dof numbering
   if (NrdofEE .ne.  &
       nord1*nrdofH2*nrdofH3+nrdofH1*nord2*nrdofH3+nrdofH1*nrdofH2*nord3) then 
      write(*,*) 'elem_fi: NrdofEE inconsistency'
      stop 1
   endif

   if (NrdofQ .ne. nrdofQ1_tr*nrdofQ2_tr*nrdofQ3_tr) then 
      write(*,*) 'elem_fi: NrdofQ inconsistency'
      stop 2
   endif
!


   ! start = omp_get_wtime()

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
!-----------------------------------------------------------------------------
!                3D INTEGRATION LOOPS BEGIN
!-----------------------------------------------------------------------------
!
!..1st integration loop (through x direction)
   do px=1,nintx
      xi1=xilocx(px)
      wt1=wlocx(px)
!  ...call 1D shape functions for coordinate 1
      call shape1HH(xi1,nord1, nrdofH1,shapH1(:,1),shapH1(:,2))

!  ...Initialize auxiliary matrices B: Gram matrix
      auxFF_B = ZERO; auxCC_B = ZERO;  auxFC_B = ZERO;  auxCF_B = ZERO
!      
!  ...Initialize auxiliary matrices B: Stiffness and load
      stiffQF_B = ZERO; stiffQC_B = ZERO; load1_B = ZERO; load2_B = ZERO
!
!  ...2nd integration loop (through y direction)
      do py=1,ninty
!         
         xi2=xilocy(py);  wt2=wlocy(py)
!         
!     ...Shape function subroutine is called only once, when 
!        px=1 and stored in sH2p(:,py) and dsH2p(:,py)
         if (px.eq.1) then
            sH2p (:,py) = 0.d0
            dsH2p(:,py) = 0.d0
            call shape1HH(xi2,nord2,nrdofH2,sH2p(:,py),dsH2p(:,py))
         endif
!     ...Copy shape functions in coord. 2 previously evaluated
         shapH2(:,1)=sH2p(:,py); shapH2(:,2)=dsH2p(:,py)
!         
!     ...Initialize auxiliary matrices A: Gram matrix
         auxFF_A = ZERO; auxCC_A =ZERO; auxFC_A =ZERO; auxCF_A =ZERO
!         
!     ...Initialize auxiliary matrices A: Stiffness matrix and load vector
         stiffQF_A = ZERO; stiffQC_A = ZERO; load1_A = ZERO; load2_A = ZERO
!     ...loop over quadrature points in direction \xi_3
         do pz=1,nintz
!            
!        ...read quadrature point location and weight
            xi3=xilocz(pz); wt3=wlocz(pz)
!            
!        ...store 3D quadrature point
            xip(1)=xi1;  xip(2)=xi2; xip(3)=xi3
!            
!        ...Shape function subroutine is called only once, when
!           px=py=1 and stored in sH3p(:,pz) and dsH3p(:,pz)
            if (px*py.eq.1) then
               call shape1HH(xi3,nord3,nrdofH3,sH3p(:,pz),dsH3p(:,pz))
            endif
!        ...Copy shape functions in coord. 3 previously evaluated
            shapH3(:,1)=sH3p(:,pz);  shapH3(:,2)=dsH3p(:,pz)
!
!        ...Compute shape functions needed for geometry - 3D H1 shape functions
            call shape3DH(etype,xip,norder,norient_edge,norient_face, nrdofH,shapH,gradH)
!            
!        ...Geometry map
            call geom3D(Mdle,xip,xnod,shapH,gradH,nrdofH,x,dxdxi,dxidx,rjac,iflag)
            if (iflag.ne.0) then
                  write(*,1000) Mdle,rjac
 1000             format('elem_fi: Negative JacobiancMdle,rjac=',i8,2x,e12.5)
               stop 1
            endif
!            
!        ...compute total quadrature weight
            wt123=wt1*wt2*wt3
!            
!        ...compute Jacobian determinant * quadrature weight
            weighthh=wt123*rjac
!            
!        ...Determine D = J^-1 * J^-T.  Multiply by appropriate weight
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
!            
!        ...compute inverse Jacobian determinant * quadrature weight
            weightvv=wt123/rjac
!        ...Determine C = J^T * J.  Multiply by appropriate weight
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

!        ...HERE STARTS COMPUTATION OF VOLUME INTEGRALS: GRAM, STIFFNESS AND LOAD MATRICES
!
!        ...loop through 1D dof in z direction
            do i3=1,nrdofH3
!           ...loop through 1D dof in z direction
               do j3=1,nrdofH3
!              ...combine indices i3 and j3 into k3
                  k3=(i3-1)*nrdofH3+j3
!              ...loops on vector components a (for test function), b (for trial)
!              ...The way the 3D H(curl) broken functions are organized, allow for
!              ...computing only the upper-triangular blocks of matrices (b>=a)
                  do a=1,3
                     do b=a,3
!                    ...indices sb and sa are 1 or 2, depending on a and b
                        sb=1+deltak(b,3)
                        sa=1+deltak(a,3)
!                    ...indices idxb, idxa, are equal to j3,i3, shifted by 0 or 1,
!                       depending on whether we need the H1 value (j3+0,i3+0)
!                       or the L2 (j3+1,i3+1), which is determined
!                       by the vector component the shape function lies at
                        idxb=j3+deltak(b,3)
                        idxa=i3+deltak(a,3)
!                    ...we need to check that idxb, idxa don't exceed dimension of
!                       1D H1 enriched test space
                        if ((idxb.le.nrdofH3).and.(idxa.le.nrdofH3)) then
!                       ...accumulate innermost 1D integral for FF term in Gram matrix
                           auxFF_A(b,a,k3) = auxFF_A(b,a,k3)    &
                                           + (shapH3(idxa,sa)    &
                                           * shapH3(idxb,sb)*D(a,b))
!                       ...loop over components a+alph, b+beta, (modulo 3)
!                          where the curl of the shape functions lie
                           do beta=1,2
                              do alph=1,2
!                             ...compute a+alph, b+beta, (modulo 3)
                                 idxbeta=mod(b+beta-1,3)+1
                                 idxalph=mod(a+alph-1,3)+1
!                             ...determine indices sb,sa for the curl components
                                 sb=1+1-deltak(idxbeta,3)
                                 sa=1+1-deltak(idxalph,3)
!                             ...accumulate innermost 1D integral for CC term in Gram matrix
                                 auxCC_A(alph,beta,b,a,k3) = auxCC_A(alph,beta,b,a,k3)  & 
                                                     + shapH3(idxa,sa)*shapH3(idxb,sb)  &
                                                     *(-1)**(alph+beta)*C(idxalph,idxbeta)
                              enddo
                           enddo
!                       ...loop over components a+alph, where the curl of
!                          the TEST shape function, for the CF term of Gram matrix
                           do alph=1,2
                              idxalph=mod(a+alph-1,3)+1
                              sb=1+deltak(b,3)
                              sa=1+1-deltak(idxalph,3)
!                          ...the only nonzero is when  a+alph == b
!                             ...accumulate innermost 1D integral for CF term
                                 auxCF_A(alph,b,a,k3)=auxCF_A(alph,b,a,k3) &
                                      +shapH3(idxa,sa)*shapH3(idxb,sb)     &
                                      *(-1)**(alph-1)*wt123
                           enddo
!                       ...loop over components b+beta, where the curl of
!                          the TRIAL shape function, for the FC term of Gram matrix
                           do beta=1,2
                              idxbeta=mod(b+beta-1,3)+1
                              sb=1+1-deltak(idxbeta,3)
                              sa=1+deltak(a,3)
!                          ...the only nonzero is when  b+beta == a
!                          ...accumulate innermost 1D integral for FC term
                              auxFC_A(beta,b,a,k3) = auxFC_A(beta,b,a,k3)             &
                                                   + shapH3(idxa,sa)*shapH3(idxb,sb)  &
                                                   * (-1)**(beta-1)*wt123
!                       ...loop over beta ends
                           enddo
                        endif
!                 ...end of loop through b            
                     enddo      
!              ...end of loop through a
                  enddo                           
!           ...end of loop through 1D dof in z direction
               enddo

!           ...Integration of enriched stiffness matrix and load vector
!           ...start new loop over j3 - L2 trial functions
               do j3=1,nrdofQ3_tr
!              ...combine indices i3, j3 into k3
                  k3=(i3-1)*nrdofQ3_tr+j3
!              ...loop over vector components
                  do b=1,3
                     do a=1,3
!                    ...determine indices for test function
                        sa=1+deltak(a,3)
                        idxa=i3+deltak(a,3)
!                    ...check that idxa don't exceed dimension
!                       of 1D H1 enriched test space
                        if (idxa.le.nrdofH3) then
!                       ...accumulate innermost 1D integral for QE term
                           stiffQF_A(a,b,k3) = stiffQF_A(a,b,k3) &
                                             + (dxidx(a,b)*shapH3(j3+1,2)) &
                                             * shapH3(idxa,sa)
!
!                       ...loop over components a+alph, required for curl of test f
                           do alph=1,2
                              idxalph=mod(a+alph-1,3)+1
                              sa=1+1-deltak(idxalph,3)
!                          ...accumulate innermost 1D integral for QC term
                              stiffQC_A(alph,a,b,k3) = stiffQC_A(alph,a,b,k3) &
                                                     + shapH3(j3+1,2)*dxdxi(b,idxalph) &
                                                     * (-1)**(alph-1)*shapH3(idxa,sa)
                           enddo
                        endif
!                    ...loop over a ends
                     enddo
!                 ...loop over b ends
                  enddo
!              ...loop over j3 ends
               enddo
!           ...loop over vector component a of Hcurl test functions
               do a=1,3
                  sa=1+deltak(a,3)
                  idxa=i3+deltak(a,3)
                  if (idxa.le.nrdofH3) then
!                 ...accumulate innermost integral of load vector
                     load1_A(a,i3)=load1_A(a,i3) &
                                  +(dxidx(a,1)*zf(1)+dxidx(a,2)*zf(2)+dxidx(a,3)*zf(3)) &
                                  *shapH3(idxa,sa)*rjac
                     load2_A(a,i3)=load2_A(a,i3) &
                                  +(dxidx(a,1)*zf(4)+dxidx(a,2)*zf(5)+dxidx(a,3)*zf(6)) &
                                  *shapH3(idxa,sa)*rjac             
                  endif
!              ...loop over a ends
               enddo
!
!        ...end of loop through 1D dof in z direction
            enddo
!
!     ...end of loop through z direction
         enddo   
!         
!     ...Computation of middle 1D integrals
!
!     ...Compute middle integrals in Gram matrix terms
!
!     ...loop over Hcurl i3 test function
         do i3=1,nrdofH3
!        ...loop over Hcurl j3 trial function
            do j3=1,nrdofH3
!           ...combine i3 and j3 into k3
               k3=(i3-1)*nrdofH3+j3
!           ...loop over Hcurl i3 test function
               do i2=1,nrdofH2
!              ...loop over Hcurl j3 trial function
                  do j2=1,nrdofH2
!                 ...combine i2 and j2 into k2
                     k2=(i2-1)*nrdofH2+j2
!                 ...loop over vector components a (test), b(trial)
!                 ...Only upper blocks of matrix ( b>=a )are computed
                     do a=1,3
                        do b=a,3
!                       ...determine indices for shape functions
                           sb=1+deltak(b,2)
                           sa=1+deltak(a,2)
                           idxb=j2+deltak(b,2)
                           idxa=i2+deltak(a,2)
!                       ...check that dimensions are not surpassed by idxa,idxb
                           if((idxb.le.nrdofH2).and.(idxa.le.nrdofH2)) then
!                          ...accumulate middle 1D integral of FF term in Gram
                              auxFF_B(b,a,k2,k3) = auxFF_B(b,a,k2,k3)      &
                                                 + shapH2(idxa,sa)         &
                                                 * shapH2(idxb,sb)         &
                                                 * auxFF_A(b,a,k3)
!                          ...loop over b+beta, a+alph, for curl of trial and test
                              do beta=1,2 
                                 do alph=1,2
                                    idxbeta=mod(b+beta-1,3)+1
                                    idxalph=mod(a+alph-1,3)+1
                                    sb=1+1-deltak(idxbeta,2)
                                    sa=1+1-deltak(idxalph,2)
!                                ...accumulate middle 1D integral of CC term in Gram
                                    auxCC_B(alph,beta,b,a,k2,k3) =          &
                                                       auxCC_B(alph,beta,b,a,k2,k3)     &
                                                     + shapH2(idxa,sa)*shapH2(idxb,sb)  &
                                                     * auxCC_A(alph,beta,b,a,k3)
                                 enddo
                              enddo
!                          ...accumulate for CF term
                              do alph=1,2
                                 idxalph=mod(a+alph-1,3)+1
                                 sb=1+deltak(b,2)
                                 sa=1+1-deltak(idxalph,2)
                                 if (idxalph.eq.b) then
                                    auxCF_B(alph,b,a,k2,k3)=              &
                                          auxCF_B(alph,b,a,k2,k3)         &
                                         +shapH2(idxa,sa)*shapH2(idxb,sb) &
                                         *auxCF_A(alph,b,a,k3)
                                 endif
                              enddo
!                          ...accumulate for FC term
                              do beta=1,2
                                 idxbeta=mod(b+beta-1,3)+1
                                 sb=1+1-deltak(idxbeta,2)
                                 sa=1+deltak(a,2)
                                 if (idxbeta.eq.a) then
                                    auxFC_B(beta,b,a,k2,k3)=              &
                                          auxFC_B(beta,b,a,k2,k3)         &
                                         +shapH2(idxa,sa)*shapH2(idxb,sb) &
                                         *auxFC_A(beta,b,a,k3)
                                 endif
                              enddo
                           endif
!                       ...loop over b ends
                        enddo
!                    ...loop over a ends
                     enddo
!                 ...loop over j2 ends
                  enddo
!              ...loop over i2 ends
               enddo
!           ...loop over j3 ends
            enddo

!        ...Compute middle integrals for stiffness terms
            do j3=1,nrdofQ3_tr
               k3=(i3-1)*nrdofQ3_tr+j3
               do i2=1,nrdofH2
                  do j2=1,nrdofQ2_tr
                     k2=(i2-1)*nrdofQ2_tr+j2
                     do b=1,3
                        do a=1,3
                           sa=1+deltak(a,2)
                           idxa=i2+deltak(a,2)
                           if (idxa.le.nrdofH2) then
                         ! ...accumulate zc*QE term
                              stiffQF_B(a,b,k2,k3) = stiffQF_B(a,b,k2,k3)             &
                                                   + shapH2(j2+1,2)*shapH2(idxa,sa)   &
                                                   * stiffQF_A(a,b,k3)
!                          ...accumulate QC term
                              do alph=1,2
                                 idxalph=mod(a+alph-1,3)+1
                                 sa=1+1-deltak(idxalph,2)
                                 stiffQC_B(alph,a,b,k2,k3) = stiffQC_B(alph,a,b,k2,k3) &
                                                     + shapH2(j2+1,2)*shapH2(idxa,sa)  &
                                                     * stiffQC_A(alph,a,b,k3)
                              enddo
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo

!        ...Compute middle integral of load vector
            do i2=1,nrdofH2
               do a=1,3
                  sa=1+deltak(a,2)
                  idxa=i2+deltak(a,2)
                  if (idxa.le.nrdofH2) then
                     load1_B(a,i2,i3)=load1_B(a,i2,i3)+shapH2(idxa,sa)*load1_A(a,i3)
                     load2_B(a,i2,i3)=load2_B(a,i2,i3)+shapH2(idxa,sa)*load2_A(a,i3)
                  endif
               enddo
            enddo
!     ...loop over i3 ends
         enddo
!
!     ...end of loop through y direction
      enddo   
!
!  ...FINAL COMPUTATION OF GRAM MATRIX
!
!  ...loop over Hcurl trial shape function identified by indices j1,j2,j3,b
      do b=1,3
         do j3=1,nrdofH3
            do j2=1,nrdofH2
               do j1=1,nrdofH1
                  sb=1+deltak(b,1)
                  idxb=j1+deltak(b,1)
                  idxb2=j2+deltak(b,2)
                  idxb3=j3+deltak(b,3)
!              ...determine index of 3D Hcurl trial shape function using j1,j2,j3,b
                  if ((idxb.le.nrdofH1).and.(idxb2.le.nrdofH2).and.(idxb3.le.nrdofH3)) then
                     select case(b)
                        case(1)
                        m2=j1+nord1*(j2-1)+nord1*nrdofH2*(j3-1)
                        case(2)
                        m2=nord1*nrdofH2*nrdofH3 &
                           +j1+nrdofH1*(j2-1)+nrdofH1*nord2*(j3-1)
                        case(3)
                        m2=nord1*nrdofH2*nrdofH3 &
                           +nrdofH1*nord2*nrdofH3 &
                           +j1+nrdofH1*(j2-1)+nrdofH1*nrdofH2*(j3-1)
                     end select
!                     
!                 ...loop over Hcurl test shape function identified by indices i1,i2,i3,a
                     do a=1,3
                        do i3=1,nrdofH3
                           do i2=1,nrdofH2
                              do i1=1,nrdofH1
                                 sa=1+deltak(a,1)
                                 idxa=i1+deltak(a,1)
                                 idxa2=i2+deltak(a,2)
                                 idxa3=i3+deltak(a,3)
                                 if ((idxa.le.nrdofH1).and.(idxa2.le.nrdofH2)  &
                                                            .and.(idxa3.le.nrdofH3)) then
!                                ...combine indices i3, j3 into k3
                                    k3=(i3-1)*nrdofH3+j3
!                                ...combine indices i2, j2 into k2
                                    k2=(i2-1)*nrdofH2+j2
!                                ...determine index of 3D Hcurl trial shape function
                                    select case(a)
                                       case(1)
                                       m1=i1+nord1*(i2-1)+nord1*nrdofH2*(i3-1)
                                       case(2)
                                       m1=nord1*nrdofH2*nrdofH3 &
                                         +i1+nrdofH1*(i2-1)+nrdofH1*nord2*(i3-1)
                                       case(3)
                                       m1=nord1*nrdofH2*nrdofH3+nrdofH1*nord2*nrdofH3 &
                                         +i1+nrdofH1*(i2-1)+nrdofH1*nrdofH2*(i3-1)
                                    end select
!                                ...accumulate integrals only if m1<=m2
                                    if (m1.le.m2) then
                                       sa=1+deltak(a,1)
                                       sb=1+deltak(b,1)
                                       kk = nk(2*m1-1,2*m2-1)
!                                   ...sum FF terms
                                       gramP(kk) = gramP(kk)      &
                                                 + shapH1(idxa,sa)*    &
                                                   shapH1(idxb,sb)*    &
                                                   auxFF_B(b,a,k2,k3)* &
                                                   (abs(zb1)**2 + ALPHA)

                                       ! gram(2*m1-1,2*m2-1) = gram(2*m1-1,2*m2-1) &
                                       !                     + shapH1(idxa,sa)*    &
                                       !                       shapH1(idxb,sb)*    &
                                       !                       auxFF_B(b,a,k2,k3)* &
                                       !                       (abs(zb1)**2 + ALPHA)
!
!                                   ...sum CC terms
                                       do beta=1,2
                                          do alph=1,2
                                             idxbeta=mod(b+beta-1,3)+1
                                             idxalph=mod(a+alph-1,3)+1
                                             sb=1+1-deltak(idxbeta,1)
                                             sa=1+1-deltak(idxalph,1)
                                             gramP(kk) = gramP(kk) +        &
                                                         shapH1(idxa,sa)*   &
                                                         shapH1(idxb,sb)*   &
                                                         auxCC_B(alph,beta,b,a,k2,k3)
                                          enddo
                                       enddo
!
                                       kk = nk(2*m1-1,2*m2)
!                                   ...sum CF terms
                                       do alph=1,2
                                          idxalph=mod(a+alph-1,3)+1
                                          sb=1+deltak(b,1)
                                          sa=1+1-deltak(idxalph,1)
                                          if (idxalph.eq.b) then
                                             gramP(kk) = gramP(kk)       &
                                                   ! +conjg(zb2)*auxCF_B(alph,b,a,k2,k3) &
                                                   +zb2*auxCF_B(alph,b,a,k2,k3) &
                                                   *shapH1(idxa,sa)*shapH1(idxb,sb)
                                          endif
                                       enddo
!                                   ...sum FC terms
                                       do beta=1,2
                                          idxbeta=mod(b+beta-1,3)+1
                                          sb=1+1-deltak(idxbeta,1)
                                          sa=1+deltak(a,1)
                                          if (idxbeta.eq.a) then
                                             gramP(kk) = gramP(kk)    &
                                                   ! +conjg(zb1)*auxFC_B(beta,b,a,k2,k3)  &
                                                   +zb1*auxFC_B(beta,b,a,k2,k3)  &
                                                   *shapH1(idxa,sa)*shapH1(idxb,sb)
                                          endif
                                       enddo

                                       if (m1.ne.m2) then

                                          kk = nk(2*m1,2*m2-1)
!                                      ...sum CF terms
                                          do alph=1,2
                                             idxalph=mod(a+alph-1,3)+1
                                             sb=1+deltak(b,1)
                                             sa=1+1-deltak(idxalph,1)
                                             if (idxalph.eq.b) then
                                                gramP(kk) = gramP(kk)       &
                                                   ! +zb1*auxCF_B(alph,b,a,k2,k3) &
                                                   +conjg(zb1)*auxCF_B(alph,b,a,k2,k3) &
                                                   *shapH1(idxa,sa)*shapH1(idxb,sb)
                                             endif
                                          enddo
!                                      ...sum EF terms
                                          do beta=1,2
                                             idxbeta=mod(b+beta-1,3)+1
                                             sb=1+1-deltak(idxbeta,1)
                                             sa=1+deltak(a,1)
                                             if (idxbeta.eq.a) then
                                                gramP(kk) = gramP(kk)    &
                                                      ! +zb2*auxFC_B(beta,b,a,k2,k3)  &
                                                      +conjg(zb2)*auxFC_B(beta,b,a,k2,k3)  &
                                                      *shapH1(idxa,sa)*shapH1(idxb,sb)
                                             endif
                                          enddo
                                       endif
!
                                       kk = nk(2*m1  ,2*m2  )
!                                   ...sum FF terms
                                       sb=1+deltak(b,1)
                                       sa=1+deltak(a,1)
                                       gramP(kk)= gramP(kk)      &
                                                + shapH1(idxa,sa)*    &
                                                  shapH1(idxb,sb)*    &
                                                  auxFF_B(b,a,k2,k3)* &
                                                  (abs(zb2)**2 + ALPHA)
!                                   ...sum CC terms
                                       do beta=1,2
                                          do alph=1,2
                                             idxbeta=mod(b+beta-1,3)+1
                                             idxalph=mod(a+alph-1,3)+1
                                             sb=1+1-deltak(idxbeta,1)
                                             sa=1+1-deltak(idxalph,1)
                                             gramP(kk) = gramP(kk)+    &
                                                         shapH1(idxa,sa)*   &
                                                         shapH1(idxb,sb)*   &
                                                         auxCC_B(alph,beta,b,a,k2,k3)
                                             ! gram(2*m1,2*m2) = gram(2*m1,2*m2)+    &
                                             !          shapH1(idxa,sa)*   &
                                             !          shapH1(idxb,sb)*   &
                                             !          auxCC_B(alph,beta,b,a,k2,k3)                                                         
                                          enddo
                                       enddo
                                     endif
                                 endif
!
!                          ...end of loop through i1
                              enddo
!                       ...end of loop through i2
                           enddo
!                    ...end of loop through i3
                        enddo
!
!                 ...end of loop through a
                     enddo
                  endif
!           ...end of loop through j1
               enddo
! 
!        ...end of loop through j2
            enddo
!
!     ...end of loop through j3
         enddo
!  ...end of loop through b         
      enddo
!      
!  ...Compute final integral for Stiffness matrix
!
!  ...loop over trial L2 shape functions
      do j3=1,nrdofQ3_tr
         do j2=1,nrdofQ2_tr
            do j1=1,nrdofQ1_tr
!           ...determine 3D L2 shape function index
               m2=j1+(j2-1)*nrdofQ1_tr+(j3-1)*nrdofQ1_tr*nrdofQ2_tr
!           ...loop over Hcurl test shape function identified by indices i1,i2,i3,a
               do a=1,3
                  do i3=1,nrdofH3-deltak(a,3)
                     do i2=1,nrdofH2-deltak(a,2)
                        do i1=1,nrdofH1-deltak(a,1)
                           sa=1+deltak(a,1)
                           idxa=i1+deltak(a,1)
                           idxa2=i2+deltak(a,2)
                           idxa3=i3+deltak(a,3)
                           if ((idxa.le.nrdofH1).and.(idxa2.le.nrdofH2)   &
                               .and.(idxa3.le.nrdofH3)) then
!                          ...combine indices i3, j3 into k3
                              k3=(i3-1)*nrdofQ3_tr+j3
!                          ...combine indices i2, j2 into k2
                              k2=(i2-1)*nrdofQ2_tr+j2
!                          ...determine index for 3D Hcurl test function
                              select case(a)
                                 case(1)
                                 m1=i1+nord1*(i2-1)+nord1*nrdofH2*(i3-1)
                                 case(2)
                                 m1=nord1*nrdofH2*nrdofH3 &
                                   +i1+nrdofH1*(i2-1)+nrdofH1*nord2*(i3-1)
                                 case(3)
                                 m1=nord1*nrdofH2*nrdofH3+nrdofH1*nord2*nrdofH3 &
                                   +i1+nrdofH1*(i2-1)+nrdofH1*nrdofH2*(i3-1)
                              end select
!                          ...loop over vector components of L2 trial functions
                              do b=1,3

!                             ...first 3 trial L2 components 
!                             ...(E,curl(F)) (testing with 1st H(curl))
                                 k = (m2-1)*6 + b
                                 do alph=1,2
                                    idxalph=mod(a+alph-1,3)+1
                                    sa=1+1-deltak(idxalph,1)


                                    stiffE_Qt(k,2*m1-1) = stiffE_Qt(k,2*m1-1) &
                                                       + shapH1(idxa,sa)*shapH1(j1+1,2) &
                                                       * stiffQC_B(alph,a,b,k2,k3)                                                 
                                 enddo
!
!                            ...-iωε(E,G) (testing with 2nd H(curl))
                                 sa=1+deltak(a,1)

                                 stiffE_Qt(k,2*m1) = stiffE_Qt(k,2*m1) &
                                                  + shapH1(idxa,sa)*shapH1(j1+1,2) &
                                                  * stiffQF_B(a,b,k2,k3)*conjg(za1)                 
!
!                             ...second 3 trial L2 components 
                                 k = (m2-1)*6 + 3+ b   
                                 sa=1+deltak(a,1)
!                             ...iωμ(H,F) (testing with 1st H(curl))

                                 stiffE_Qt(k,2*m1-1)= stiffE_Qt(k,2*m1-1) &
                                                    + shapH1(idxa,sa)*shapH1(j1+1,2) &
                                                    * stiffQF_B(a,b,k2,k3)*za2                   

!                             ...(H,curl(G)) (testing with 2st H(curl))
                                 do alph=1,2
                                    idxalph=mod(a+alph-1,3)+1
                                    sa=1+1-deltak(idxalph,1)
                                    stiffE_Qt(k,2*m1) = stiffE_Qt(k,2*m1) &
                                                     + shapH1(idxa,sa)*shapH1(j1+1,2) &
                                                     * stiffQC_B(alph,a,b,k2,k3)                 
                                 enddo
!
!                          ...end of loop through vector components                                 
                              enddo
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
!     Final computation of load vector
!
!  ...loop over Hcurl test shape function identified by indices i1,i2,i3,a
      do a=1,3
!         
         do i3=1,nrdofH3
!
            do i2=1,nrdofH2
!
               do i1=1,nrdofH1
!
                  sa=1+deltak(a,1)
                  idxa=i1+deltak(a,1)
                  idxa2=i2+deltak(a,2)
                  idxa3=i3+deltak(a,3)
!                  
                  if ((idxa.le.nrdofH1).and.(idxa2.le.nrdofH2).and.(idxa3.le.nrdofH3)) then
!                 ...combine indices i3, j3 into k3
                     k3=(i3-1)*nrdofQ3_tr+j3
!                 ...combine indices i2, j2 into k2
                     k2=(i2-1)*nrdofQ2_tr+j2
!                 ...determine index for 3D Hcurl test function
                     select case(a)
                        case(1)
                        m1=i1+nord1*(i2-1)+nord1*nrdofH2*(i3-1)
                        case(2)
                        m1=nord1*nrdofH2*nrdofH3 &
                          +i1+nrdofH1*(i2-1)+nrdofH1*nord2*(i3-1)
                        case(3)
                        m1=nord1*nrdofH2*nrdofH3+nrdofH1*nord2*nrdofH3 &
                          +i1+nrdofH1*(i2-1)+nrdofH1*nrdofH2*(i3-1)
                     end select
!                 ...accumulate outermost integral
                     bloadE(2*m1-1) = bloadE(2*m1-1)+shapH1(idxa,sa)*load1_B(a,i2,i3)
                     bloadE(2*m1)   = bloadE(2*m1)  +shapH1(idxa,sa)*load2_B(a,i2,i3)
                  endif
!                  
!           ...end of loop through i1
               enddo
!               
!        ...end of loop through i2
            enddo
!                 
!     ...end of loop through i3
         enddo
!         
!  ...end of loop through a         
      enddo
!
!..end of loop through x direction
   enddo   
!,*

   ! finish = omp_get_wtime()

   ! write(*,*) 'sum fact time = ', finish-start

!..find ndof associated with the mdle node of the element
!
   call compute_enriched_order(nordP, norderP)
   call ndof_nod(etype,norder(nre+nrf+1) , ndofHmdl,ndofEmdl,ndofVmdl,ndofQmdl)
   call ndof_nod(etype,norderP(nre+nrf+1), ndofHHmdl,ndofEEmdl,ndofVVmdl,ndofQQmdl)
!
   nrdofEEi = NrdofEE-ndofEEmdl
!
   allocate(idxEE(nrdofEEi))
   ik = 0
   do i3 = 1,nrdofH3
      do i2 = 1,nrdofH2
         do i1 = 1,nord1
            if (i2 .lt.3 .or. i3 .lt.3) then
               ik = ik +1
               m1 = i1+nord1*(i2-1) + nord1*nrdofH2*(i3-1)
               idxEE(ik) = m1
            endif
         enddo
      enddo
   enddo
   do i3 = 1,nrdofH3
      do i2 = 1,nord2
         do i1 = 1,nrdofH1
            if (i1 .lt.3 .or. i3 .lt.3) then
               ik = ik +1
               m1 = nord1*nrdofH2*nrdofH3 + i1 + nrdofH1*(i2-1) + nrdofH1*nord2*(i3-1)
               idxEE(ik) = m1
            endif   
         enddo
      enddo
   enddo
   do i3 = 1,nord3
      do i2 = 1,nrdofH2
         do i1 = 1,nrdofH1
            if (i1 .lt.3 .or. i2 .lt.3) then
               ik = ik +1
               m1 = nrdofH1*nord2*nrdofH3 + nord1*nrdofH2*nrdofH3   &
                  + i1+nrdofH1*(i2-1)+nrdofH1*nrdofH2*(i3-1)
               idxEE(ik) = m1
            endif
         enddo
      enddo
   enddo
!
!----------------------------------------------------------------------
!     B O U N D A R Y    I N T E G R A L S                           
!----------------------------------------------------------------------

   ! start = omp_get_wtime()

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
         call face_param(etype,if,t, xi,dxidt)
! 
!     ...broken H(curl) shape functions for the enriched test space
         call shape3EE(etype,xi,nordP, nrdof,shapEE,curlEE)
         if (nrdof .ne. NrdofEE) then
            write(*,*) 'elem_fi: INCONSISTENCY NrdofEE. stop.'
            stop
         endif
!         
!     ...determine element H1 shape functions (for geometry)
         call shape3DH(etype,xi,norder,norient_edge,norient_face, nrdofH,shapH,gradH)
! 
!     ...determine element H(curl) shape functions (for fluxes)
!
         call shape3DE(etype,xi,norderi,norient_edge,norient_face, nrdof,shapE,curlE)
         if (nrdof .ne. NrdofEi) then
            write(*,*) 'elem_fi: INCONSISTENCY NrdofEi. stop.'
            stop
         endif
!
!     ...geometry
         call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign,x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
         weight = bjac*wtloc(l)
!
!     ...check if on impedance boundary
!        ( < n x H , G > = < n x n x E , G > + < zg , G > )
         if (ibc(if,1) .eq. 9) then
!        ...get boundary data
            call getg(mdle,x,rn,ibc(if,1),zg)   

!        ...loop through enriched H(curl) test functions
            do k = 1,nrdofEEi
               k1 = idxEE(k)
               fldF(1:3) = shapEE(1,k1)*dxidx(1,1:3)   &
                         + shapEE(2,k1)*dxidx(2,1:3)   &
                         + shapEE(3,k1)*dxidx(3,1:3)  
               fldG(1:3) = fldF(1:3)        
!               
!           ...accumulate for the load vector     
               n = 2*k1   
               bloadE(n) = bloadE(n) &
                         - (fldG(1)*zg(1)+fldG(2)*zg(2)+fldG(3)*zg(3))*weight
!                         
!           ...loop through H(curl) trial functions
               do k2=1,nrdofEi
                  fldE(1:3) = shapE(1,k2)*dxidx(1,1:3)  &
                            + shapE(2,k2)*dxidx(2,1:3)  &
                            + shapE(3,k2)*dxidx(3,1:3) 
                  fldH(1:3) = fldE(1:3)

                  call cross_product(rn,fldE, rnE)
                  call cross_product(rn,rnE, rn2E)
                  call dot_product(rnE,fldF,rnEF) 
                  call dot_product(rn2E,fldG,rn2EG) 
!
!              ...accumulate for the skeleton stiffness matrix
                  n = 2*k1-1; m = 2*k2-1
                  stiffE_Et(m,n) = stiffE_Et(m,n) + rnEF*weight 
!
                  n = 2*k1; m = 2*k2-1
                  stiffE_Et(m,n) = stiffE_Et(m,n) + rn2EG*weight
!
!           ...end of loop through H(curl) trial functions
               enddo
!               
!        ...end of loop through H(curl) test functions
            enddo
!     ...regular boundary
         else
!
!        ...loop through enriched H(curl) test functions
            do k = 1,nrdofEEi
               k1 = idxEE(k)
               fldF(1:3) = shapEE(1,k1)*dxidx(1,1:3)   &
                         + shapEE(2,k1)*dxidx(2,1:3)   &
                         + shapEE(3,k1)*dxidx(3,1:3)  
               fldG(1:3) = fldF(1:3)          
!
!           ...loop through H(curl) trial functions
               do k2=1,nrdofEi
                  fldE(1:3) = shapE(1,k2)*dxidx(1,1:3)  &
                            + shapE(2,k2)*dxidx(2,1:3)  &
                            + shapE(3,k2)*dxidx(3,1:3) 
                  fldH(1:3) = fldE(1:3)
!                         
                  call cross_product(rn,fldE, rnE)
                  call cross_product(rn,fldH, rnH)
                  call dot_product(rnE,fldF,rnEF)
                  call dot_product(rnH,fldG,rnHG)               
!
!              ...accumulate for the skeleton stiffness matrix
                  n = 2*k1-1; m = 2*k2-1
                  stiffE_Et(m,n) = stiffE_Et(m,n) + rnEF*weight 
!               
                  n = 2*k1; m = 2*k2
                  stiffE_Et(m,n) = stiffE_Et(m,n) + rnHG*weight
!
!           ...end of loop through H(curl) trial functions
               enddo
!        ...end of loop through H(curl) test functions
            enddo
!
         endif            
! 
!  ...end of loop through integration points
      enddo
!..end of loop through faces  
   enddo

   ! finish = omp_get_wtime()

   ! write(*,*) 'boundary int time = ', finish-start

   deallocate(idxEE)

   deallocate(auxFF_A,auxFF_B)
   deallocate(auxCC_A,auxCC_B)
   deallocate(auxFC_A,auxFC_B)
   deallocate(auxCF_A,auxCF_B)
!
   deallocate(stiffQF_A)
   deallocate(stiffQF_B)
   deallocate(stiffQC_A)
   deallocate(stiffQC_B)
   deallocate(load1_A, load2_A)
   deallocate(load1_B, load2_B)
!
!----------------------------------------------------------------------
!      Construction of the DPG system
!----------------------------------------------------------------------
!
   i1 = NrTest ; j1 = 2*NrdofEi ; j2 = 6*NrdofQ
!
   stiff_ALL(1:i1,1:j1)        = transpose(StiffE_Et(1:j1,1:i1))
   stiff_ALL(1:i1,j1+1:j1+j2)  = transpose(StiffE_Qt(1:j2,1:i1))
   stiff_ALL(1:i1,j1+j2+1)     = BloadE(1:i1)
!
!..A. Compute Cholesky factorization of Gram Matrix
!


   call ZPPTRF(uplo,NrTest,gramP,info)
   if (info.ne.0) then
      write(*,*) 'elem_fi: ZPPTRF: Mdle,info = ',Mdle,info
      stop
   endif
!

!..B. Solve triangular system to obtain B
   call ZTPTRS(uplo,'C','N',NrTest,NrTrial+1,gramP,stiff_ALL,NrTest,info)
   if (info.ne.0) then
      write(*,*) 'elem_fi: ZTPTRS: Mdle,info = ',Mdle,info
      stop
   endif

!..C. Matrix multiply: B^* B
   call ZHERK(uplo,'C',NrTrial+1,NrTest,ZONE,stiff_ALL,NrTest,ZERO,zaloc,NrTrial+1)
!
!..D. Fill lower triangular part of Hermitian matrix
   do i=1,NrTrial
      zaloc(i+1:NrTrial+1,i) = conjg(zaloc(i,i+1:NrTrial+1))
   enddo
!
!..E. Fill ALOC and BLOC matrices
   if (.not. ISTC_FLAG) then
!  ...make sure ALOC,BLOC are zero for mdl dofs of interface variables
      ALOC(1,1)%array = ZERO
      ALOC(1,2)%array = ZERO
      ALOC(2,1)%array = ZERO
      BLOC(1)%array = ZERO
   endif
!


   BLOC(1)%array(1:j1,1) = zaloc(1:j1,j1+j2+1)
   BLOC(2)%array(1:j2,1) = zaloc(j1+1:j1+j2,j1+j2+1)
!
   ALOC(1,1)%array(1:j1,1:j1) = zaloc(1:j1,1:j1)
   ALOC(1,2)%array(1:j1,1:j2) = zaloc(1:j1,j1+1:j1+j2)
!
   ALOC(2,1)%array(1:j2,1:j1) = zaloc(j1+1:j1+j2,1:j1)
   ALOC(2,2)%array(1:j2,1:j2) = zaloc(j1+1:j1+j2,j1+1:j1+j2)
   
   end subroutine elem_maxwell_fi




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
   
   
   
   
subroutine compute_1D_ord(etype,nord,nordx,nordy,nordz)

   use parametersDPG,     only : MAXPP
   use control,           only : INTEGRATION
   !
   !
   implicit none

   character(len=4), intent(in)  :: etype
   integer,          intent(in)  :: Nord
   integer,          intent(out) :: nordx, nordy, nordz
   integer                       :: nordh, nord1,nord2,nord3
   !
   call decode(Nord, nordh,nord3)
   call decode(nordh, nord1,nord2)


   nordx=min(nord1+INTEGRATION,MAXPP)
   nordy=min(nord2+INTEGRATION,MAXPP)
   nordz=min(nord3+INTEGRATION,MAXPP)

end subroutine compute_1D_ord

subroutine find_hmin(Mdle,Hmin)
!
   use data_structure3D

   implicit none
!  
   integer, intent(in)  :: Mdle
   real*8,  intent(out) :: Hmin
!
!..geometry dof
   real*8  :: xnod(3,MAXbrickH), h1, h2, h3
   integer :: i
!
!..determine min element size

   call nodcor(Mdle,xnod)

   h1 = dsqrt((xnod(1,2)-xnod(1,1))**2 +           &
               (xnod(2,2)-xnod(2,1))**2 + (xnod(3,2)-xnod(3,1))**2) 

   h2 = dsqrt((xnod(1,4)-xnod(1,1))**2 +           &
               (xnod(2,4)-xnod(2,1))**2 + (xnod(3,4)-xnod(3,1))**2) 

   h3 = dsqrt((xnod(1,5)-xnod(1,1))**2 +           &
               (xnod(2,5)-xnod(2,1))**2 + (xnod(3,5)-xnod(3,1))**2) 


   Hmin = min(h1,h2,h3)

end subroutine find_hmin
