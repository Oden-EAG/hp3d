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
   case(1)
      Itest(1:NR_PHYSA)=1; Itrial(1:NR_PHYSA)=1
      call elem_L2PROJ(Mdle)
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
   subroutine elem_L2PROJ(Mdle)
!
   use uweak_acoustics_module
   use control, only: INTEGRATION
   use parameters
   use parametersDPG
   use data_structure3D
   use element_data
   use physics   , only : NR_PHYSA
   use assembly  , only : ALOC,BLOC,NR_RHS
   use common_prob_data, only: SYMMETRY_TOL, TEST_NORM, OMEGA
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
!  ...workspace for exact
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
   complex*16, allocatable :: ZalocQQ(:,:), ZblocQ(:)
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
!..determine edge and face orientations
   call find_orient(Mdle, norient_edge,norient_face)
!                                                                     
!..determine nodes coordinates 
   call nodcor(Mdle, xnod)
!
!..extended load vector
   call celndof(NODES(Mdle)%type,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
   allocate(ZblocQ(4*nrdofQ))
   allocate(ZalocQQ(4*nrdofQ,4*nrdofQ))
   ZblocQ  = ZERO ;  ZalocQQ = ZERO

!
!----------------------------------------------------------------------
!     E L E M E N T    I N T E G R A L S                              |
!----------------------------------------------------------------------
!
!..use the enriched order to set the quadrature
   INTEGRATION = NORD_ADD
   call set_3Dint(etype,norder, nint,xiloc,waloc)
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
!  ...geometry map
      call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,x,dxdxi,dxidx,rjac,iflag)
!
!  ...integration weight 
      weight = rjac*wa
!
!   ......compute load vector
      call exact(x,1,                  &
                 zvalH,zdvalH,zd2valH, &
                 zvalE,zdvalE,zd2valE, &
                 zvalV,zdvalV,zd2valV, &
                 zvalQ,zdvalQ,zd2valQ) 

        zf(1:4) = zvalQ(1:4)
!
!  .....loop through L2 test functions

        do k1=1,nrdofQ
           n1 =(k1-1)*4+1
           m1 =(k1-1)*4+2
           l1 =(k1-1)*4+3
           ll1=(k1-1)*4+4

           q = shapQ(k1)/rjac; v1 = q; v2 = q ; v3 = q

           ZblocQ(n1) = ZblocQ(n1) + zf(1)*q*weight
           ZblocQ(m1) = ZblocQ(m1) + zf(2)*v1*weight
           ZblocQ(l1) = ZblocQ(l1) + zf(3)*v2*weight
           ZblocQ(ll1) = ZblocQ(ll1) + zf(4)*v3*weight
           do k2=1,nrdofQ
              n2=(k2-1)*4+1
              m2=(k2-1)*4+2
              l2=(k2-1)*4+3
              ll2=(k2-1)*4+4

              p = shapQ(k2)/rjac; u1 = p; u2 = p ; u3 = p

              ZalocQQ(n1,n2) = ZalocQQ(n1,n2) + p*q*weight
              ZalocQQ(m1,m2) = ZalocQQ(m1,m2) + u1*v1*weight 
              ZalocQQ(l1,l2) = ZalocQQ(l1,l2) + u2*v2*weight 
              ZalocQQ(ll1,ll2) = ZalocQQ(ll1,ll2) + u3*v3*weight 

          enddo
        enddo
!
   enddo
!
!
!----------------------------------------------------------------------
!      Alternative construction of normal matrix
!----------------------------------------------------------------------
! 
! 
   j1 = 4*nrdofQ

   BLOC(1)%array(1:j1,1) = ZblocQ(1:j1)
!
   ALOC(1,1)%array(1:j1,1:j1) = ZalocQQ(1:j1,1:j1)

   deallocate(ZalocQQ,ZblocQ)
!


! 
!   
   end subroutine elem_L2PROJ


 
