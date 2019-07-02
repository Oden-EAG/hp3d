!----------------------------------------------------------------------
!                                                                     
!     routine name      - elem
!                                                                     
!---------------------------------------------------------------------- 
!                                                                     
!     latest revision:  - July 17
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
    
   use physics   , only : NR_PHYSA
   use data_structure3D
!
!----------------------------------------------------------------------
!
   implicit none
!
   integer,                    intent(in)  :: Mdle
   integer,dimension(NR_PHYSA),intent(out) :: Itest,Itrial
! 
!----------------------------------------------------------------------
!
   Itest (1:NR_PHYSA) = 0
   Itrial(1:NR_PHYSA) = 0

   select case(NODES(Mdle)%case)
!      
   case(3)
      Itest(1:NR_PHYSA)=1; Itrial(1:NR_PHYSA)=1
      call elem_primal_acoustics(Mdle)
!
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
   subroutine elem_primal_acoustics(Mdle)
!
   use data_structure3D
   use parametersDPG
   use control,          only : INTEGRATION
   use assembly,         only : ALOC,BLOC,NR_RHS
   use common_prob_data, only : SYMMETRY_TOL, TEST_NORM, OMEGA,ALPHA
!
#include "syscom.blk"
!
!----------------------------------------------------------------------
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
!..geometry
   dimension xi(3),dxidt(3,2),x(3),dxdxi(3,3),dxidx(3,3), dxdt(3,2),rt(3,2),rn(3),t(2)
!
!..H1 shape functions
   dimension shapH(MAXbrickH),gradH(3,MAXbrickH)
   dimension shapHH(MAXbrickHH),gradHH(3,MAXbrickHH)
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
   complex*16, allocatable :: BLOADH(:), AP(:,:)
   complex*16, allocatable :: STIFFHH(:,:)
   complex*16, allocatable :: STIFFHV(:,:)
   complex*16, allocatable :: STIFF_ALL(:,:)
   complex*16, allocatable :: Zaloc(:,:)
!
   character*1 uplo
!
!---------------------------------------------------------------------
!
   iprint = 0
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
!..enriched order   
   norderP(1:nre) = norder(1:nre) + NORD_ADD
   norderP(nre+1:nre+nrf) = norder(nre+1:nre+nrf) + NORD_ADD*11
   norderP(nre+nrf+1) = norder(nre+nrf+1) + NORD_ADD*111 

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
!..number of trial dof associated with the mdle node
   call ndof_nod(etype,norder(nre+nrf+1),ndofHmdl,ndofEmdl,ndofVmdl,ndofQmdl)
!..total number of trial dof associated with the element
   call celndof(etype,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
!..total number of test dof associated with the element
   call celndof(etype,norderP,nrdofHH,nrdofEE,nrdofVV,nrdofQQ)   
!
!..allocate memory for local matrices
   nrTEST  = nrdofHH
   nrTRIAL = nrdofH + nrdofV  
!
   if (iprint .eq. 1)  then
!..check dimensions
      if (nrTEST .lt. nrTRIAL-ndofVmdl) then
         write(*,*) 'elem: WARNING: nrTEST smaller than nrTRIAL'
         write(*,*) 'mdle, nrTEST, nrTRIAL = ', mdle, nrTEST, nrTRIAL-ndofVmdl
      endif   
   endif   
 
   allocate(BLOADH(nrTEST))          ; BLOADH    = ZERO
   allocate(STIFFHH(nrTEST,nrdofH))  ; STIFFHH   = ZERO
   allocate(STIFFHV(nrTEST,nrdofV))  ; STIFFHV   = ZERO
   allocate(AP(nrTEST,nrTEST))       ; AP        = ZERO

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
!        ...accumulate for the Gram matrix           
            AP(k1,k2) = AP(k1,k2) &
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

   nrTRIAL = nrdofH + nrdofV
   i1 = nrTEST ; j1 = nrdofH ; j2 = nrdofV 

   allocate(STIFF_ALL(nrTEST,nrTRIAL+1))     ; STIFF_ALL = ZERO

   STIFF_ALL(1:i1,1:j1) = STIFFHH(1:i1,1:j1)
   STIFF_ALL(1:i1,j1+1:j1+j2) = STIFFHV(1:i1,1:j2)
   STIFF_ALL(1:i1,j1+j2+1) = BLOADH(1:i1)

   deallocate(STIFFHH,STIFFHV,BLOADH)
!   
   N      = nrTEST
   NRHS   = nrdofH + nrdofV + 1
!
!..diagonal scaling of the gram matrix
   ! call diag_scaling(N,NRHS,AP(1:N*(N+1)/2), STIFF_ALL(1:N,1:NRHS))
   uplo = 'U'
! 
#if C_MODE
!
!..factorize the test stiffness matrix
   call ZPOTRF(uplo,N,AP,N,info)
   if (info.ne.0) then
      write(*,*) 'elem: AP ZPPTRF: Mdle,info = ',Mdle,info
      stop 1
   endif
! 
   call ZTRSM('L',uplo,'C','N',N,NRHS,ZONE,AP,N,STIFF_ALL,N)

   allocate(Zaloc(NRHS,NRHS)) ; Zaloc = ZERO

   call ZHERK('U','C',NRHS,N,1.0d0,STIFF_ALL,N,0.0d0,Zaloc,NRHS)
! 
   deallocate(STIFF_ALL,AP)

   do i=1,NRHS-1
      Zaloc(i+1:NRHS,i) = conjg(Zaloc(i,i+1:NRHS))
   enddo
! 
#else
!
!..factorize the test stiffness matrix
   call DPOTRF(uplo,N,AP,N,info)
   if (info.ne.0) then
      write(*,*) 'elem: AP DPPTRF: Mdle,info = ',Mdle,info
      stop 1
   endif
!
   call DTRSM('L',uplo,'C','N',N,NRHS,ZONE,AP,N,STIFF_ALL,N)
!
   call DSYRK('U','C',NRHS,N,1.0d0,STIFF_ALL,N,0.0d0,Zaloc,NRHS)
! 
   do i=1,NRHS-1
      Zaloc(i+1:NRHS,i) = Zaloc(i,i+1:NRHS)
   enddo
#endif
! 
   BLOC(1)%array(1:j1,1)      = Zaloc(1:j1,j1+j2+1)
   BLOC(2)%array(1:j2,1)      = Zaloc(j1+1:j1+j2,j1+j2+1) 
   ALOC(1,1)%array(1:j1,1:j1) = Zaloc(1:j1,1:j1)
   ALOC(1,2)%array(1:j1,1:j2) = Zaloc(1:j1,j1+1:j1+j2)
   ALOC(2,1)%array(1:j2,1:j1) = Zaloc(j1+1:j1+j2,1:j1)
   ALOC(2,2)%array(1:j2,1:j2) = Zaloc(j1+1:j1+j2,j1+1:j1+j2)
! 
   deallocate(Zaloc)
! 
! 
   end subroutine elem_primal_acoustics


 
