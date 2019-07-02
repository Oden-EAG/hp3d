!----------------------------------------------------------------------
!                                                                     
!     routine name      - elem
!                                                                     
!---------------------------------------------------------------------- 
!                                                                     
!     latest revision:  - Aug 17
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
   case(1)
      Itest(1:NR_PHYSA)=1; Itrial(1:NR_PHYSA)=1
      call elem_acoustics(Mdle)
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
!     routine name      - elem_acoustics
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
   subroutine elem_acoustics(Mdle)
!
   use control, only: INTEGRATION
   use parameters
   use data_structure3D
   use element_data
   use physics   , only : NR_PHYSA
   use assembly  , only : ALOC,BLOC,NR_RHS
   use common_prob_data, only: OMEGA
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
!..3D quadrature data
   dimension xiloc(3,MAX_NINT3),waloc(MAX_NINT3)
!
!..2D quadrature data
   dimension tloc(2,MAX_NINT2),wtloc(MAX_NINT2)
!
!..BC's flags
   dimension ibc(6,NR_PHYSA)
!..workspace for trial and test variables
   dimension dq(3) , u(3), dp(1:3), v(3), vec(3)
!
!..imaginary unit
   complex*16, parameter :: zi = (0.0d0,1.0d0)   
!..needed for the trick for impedance BC   
   complex*16, allocatable :: Zaloc(:,:), Zbloc(:)
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
!..determine edge and face orientations
   call find_orient(Mdle, norient_edge,norient_face)
!  
!..determine nodes coordinates
   call nodcor(Mdle, xnod)
!                                                                   
!..get the element boundary conditions flags
   call find_bc(Mdle, ibc)
!
!..find number of trial dof for each energy space supported by the element
   call celndof(NODES(Mdle)%type,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
!
!..allocate space for matrices
   allocate(Zaloc(nrdofH,nrdofH)) ; Zaloc = ZERO
   allocate(Zbloc(nrdofH))        ; Zbloc = ZERO
!
!----------------------------------------------------------------------
!     E L E M E N T    I N T E G R A L S                              |
!----------------------------------------------------------------------
!
!..use the enriched order to set the quadrature
   call set_3Dint(etype,norder, nint,xiloc,waloc)
!      
!..loop over integration points      
   do l=1,nint
!      
      xi(1:3)=xiloc(1:3,l) ; wa=waloc(l)
!
!  ...H1 shape functions (for geometry)
      call shape3H(etype,xi,norder,norient_edge,norient_face,nrdofH,shapH,gradH)
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
      do k1=1,nrdofH
!
!     ...Piola transformation
         q = shapH(k1)
         dq(1:3) = gradH(1,k1)*dxidx(1,1:3) &
                 + gradH(2,k1)*dxidx(2,1:3) & 
                 + gradH(3,k1)*dxidx(3,1:3) 
!
!     ...accumulate for the load vector        
!
         Zbloc(k1) = Zbloc(k1) + q*zf*weight
!          
!     ...second loop through H1 test functions
         do k2=1,nrdofH
!        ...Piola transformation
            p = shapH(k2)
            dp(1:3) = gradH(1,k2)*dxidx(1,1:3) &
                    + gradH(2,k2)*dxidx(2,1:3) & 
                    + gradH(3,k2)*dxidx(3,1:3) 
! 
!        ...accumulate for the stiffness matrix

            Zaloc(k1,k2) = Zaloc(k1,k2)   &
                         + (dq(1)*dp(1)+dq(2)*dp(2)+dq(3)*dp(3)  &
                         -  OMEGA**2*q*p)*weight
!  
         enddo
!     ...end of loop through H1 test functions         
      enddo
!  ...end of loop through integration points      
   enddo
!
!----------------------------------------------------------------------
!     B O U N D A R Y    I N T E G R A L S                            |
!----------------------------------------------------------------------
!
!..loop through element faces
   do if=1,nrf
!
!  .....quit if not on the impedance boundary
      if (ibc(if,1) .ne. 9) cycle
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
      call set_2Dint(ftype,norderf, nint,tloc,wtloc)
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
!     ...geometry
         call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign,x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
         weight = bjac*wtloc(l)
!
!     ...get boundary data
         call getg(mdle,x,rn,ibc(if,1),zg)   
!
!     ...loop through H1 test functions
         do k1 = 1,nrdofH
!        ...value of the shape function at the point
            q = shapH(k1)
!        ...accumulate for the load vector
            Zbloc(k1) = Zbloc(k1) + q*zg*weight
! 
!        ...loop through H1 trial functions
            do k2 = 1,nrdofH
!           ...value of the shape function at the point
               p = shapH(k2)
!           ...accumulate for the stiffness matrix
               Zaloc(k1,k2) = Zaloc(k1,k2) + ZI*OMEGA*q*p*weight
!           ...end of loop through H1 trial functions
            enddo                                                
!        ...end of loop through H1 test functions                 
         enddo
!     ...end of loop through integration points
      enddo
!  ...end of loop through faces  
   enddo
!
! 
   BLOC(1)%array(1:nrdofH,1) = Zbloc(1:nrdofH)
!
   ALOC(1,1)%array(1:nrdofH,1:nrdofH) = Zaloc(1:nrdofH,1:nrdofH)


   deallocate(Zaloc, Zbloc)
!
! 
!   
   end subroutine elem_acoustics


 
