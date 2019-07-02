!    
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
   subroutine elem_primal_acoustics_old(Mdle)
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

   real*8 :: Mtime(10)
   integer(kind=8) :: t1,t2,clock_rate,clock_max

   real*8, allocatable :: AP_temp(:)
!
   character*1 uplo
   nk(k1,k2) = (k2-1)*k2/2+k1
!
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
   AP=ZERO
!
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
   call system_clock ( t1, clock_rate, clock_max )

!..loop over integration points      
   do l=1,nint
!      
      xi(1:3)=xiloc(1:3,l) ; wa=waloc(l)
!
!  ...H1 shape functions
      call shape3H(etype,xi,norder,norient_edge,norient_face,nrdofH,shapH,gradH)

!
!  ...discontinuous H1 shape functions
      call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
!
!  ...geometry map
      call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,x,dxdxi,dxidx,rjac,iflag)
!
!  ...integration weight 
      weight = rjac*wa
!
!  ...get the RHS
      call getf(Mdle,x, zfval)
!
!  ...loop through enriched H1 test functions
      do k1=1,nrdofHH
!       
         v1 = shapHH(k1)
!          
         dv1(1:3) = gradHH(1,k1)*dxidx(1,1:3) &
                  + gradHH(2,k1)*dxidx(2,1:3) & 
                  + gradHH(3,k1)*dxidx(3,1:3) 
!
!     ...accumulate for the load vector 
         BLOADH(k1) = BLOADH(k1) + zfval*v1*weight
!          
!     ...loop through enriched H1 trial functions
         do k2=k1,nrdofHH
!          
            v2 = shapHH(k2)
!            
            dv2(1:3) = gradHH(1,k2)*dxidx(1,1:3)  &
                     + gradHH(2,k2)*dxidx(2,1:3)  &
                     + gradHH(3,k2)*dxidx(3,1:3)
! 
!           
!        ...determine index in triangular format
            k = nk(k1,k2)
!        ...accumulate for the Gram matrix           
            AP(k) = AP(k) &
                   + (dv1(1)*dv2(1)+dv1(2)*dv2(2)+dv1(3)*dv2(3)+1.0d0*v1*v2)*weight
         enddo
!
!     ...loop through H1 trial functions
         do k2=1,nrdofH
!            
            v2 = shapH(k2)
            dv2(1:3) = gradH(1,k2)*dxidx(1,1:3)   &
                     + gradH(2,k2)*dxidx(2,1:3)   &
                     + gradH(3,k2)*dxidx(3,1:3)
! 
!        ...accumulate for the stiffness matrix
            STIFFHH(k1,k2) = STIFFHH(k1,k2)                               &
                           + (dv1(1)*dv2(1)+dv1(2)*dv2(2)+dv1(3)*dv2(3)   &
                           - OMEGA**2*v1*v2)*weight
         enddo
      enddo
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
   ! write(*,1001) Mtime(1)
 ! 1001 format('Conventional            = ', f13.5)  

! 
   end subroutine elem_primal_acoustics_old


 
