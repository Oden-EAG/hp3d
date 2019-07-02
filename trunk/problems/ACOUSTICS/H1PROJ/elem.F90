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
      call elem_H1PROJ(Mdle)
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
!     routine name      - elem_H1PROJ
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
   subroutine elem_H1PROJ(Mdle)
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
!..workspace for trial and test variables
   dimension dq(3) , dp(1:3) 
   dimension zdf(3)
!
!..workspace for exact
   dimension ZvalH(MAXEQNH), &
             ZdvalH(MAXEQNH,3),Zd2valH(MAXEQNH,3,3), &
             ZvalE(3,MAXEQNE), &
             ZdvalE(3,MAXEQNE,3),Zd2valE(3,MAXEQNE,3,3), &
             ZvalV(3,MAXEQNV), &
             ZdvalV(3,MAXEQNV,3),Zd2valV(3,MAXEQNV,3,3), &
             ZvalQ(MAXEQNQ), &
             ZdvalQ(MAXEQNQ,3),Zd2valQ(MAXEQNQ,3,3) 
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
!  ...compute RHS
      call exact(x,1,                  &
                 zvalH,zdvalH,zd2valH, &
                 zvalE,zdvalE,zd2valE, &
                 zvalV,zdvalV,zd2valV, &
                 zvalQ,zdvalQ,zd2valQ) 
!  ...exact solution and derivatives      
      zf = zvalH(1) ; zdf(1:3) = zdvalH(1,1:3)
!
!  ...loop through H1 test functions
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
         Zbloc(k1) = Zbloc(k1)        &
                   + (q*zf + dq(1)*zdf(1) + dq(2)*zdf(2) + dq(3)*zdf(3))*weight
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
                         +  q*p)*weight
!  
         enddo
!     ...end of loop through H1 test functions         
      enddo
!  ...end of loop through integration points      
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
   end subroutine elem_H1PROJ


 
