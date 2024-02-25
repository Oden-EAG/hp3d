!-------------------------------------------------------------------------
!
!     routine name      - elem
!
!-------------------------------------------------------------------------
!
!     latest revision:  - Oct 2021
!
!     purpose:          - driver for the element routine
!
!     arguments:
!        in:
!             Mdle      - an element middle node number, identified
!                         with the element
!        out:
!             Itest     - index for assembly
!             Itrial    - index for assembly
!
!-------------------------------------------------------------------------
subroutine elem(Mdle, Itest,Itrial)
!
   use data_structure3D
   use physics  , only: NR_PHYSA
   use assembly , only: ALOC,BLOC
!
   implicit none
!
   integer,                    intent(in)  :: Mdle
   integer,dimension(NR_PHYSA),intent(out) :: Itest,Itrial
!
   integer :: norder(19)
   integer :: nrdofH,nrdofE,nrdofV,nrdofQ
!
!-------------------------------------------------------------------------
!
!..activate one physics variable (H(curl)) for assembly
   Itest(1) = 1; Itrial(1) = 1
!
!..determine order of approximation
   call find_order(Mdle, norder)
!
!..find number of dof for each energy space supported by the element
   call celndof(NODES(Mdle)%ntype,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
!
!..call element integration routine
   call elem_maxwell(Mdle,nrdofE, ALOC(1,1)%array,BLOC(1)%array)
!
end subroutine elem
!
!-------------------------------------------------------------------------
!
!     routine name      - elem_maxwell
!
!-------------------------------------------------------------------------
!
!     latest revision:  - Oct 2021
!
!     purpose:          - compute element stiffness and load
!
!     arguments:
!        in:
!             Mdle      - an element middle node number, identified
!                         with the element
!             Nrdof     - trial/test degrees of freedom
!        out:
!             Zaloc     - element stiffness matrix
!             Zbloc     - element load vector
!
!-------------------------------------------------------------------------
!
subroutine elem_maxwell(Mdle,Nrdof, Zaloc,Zbloc)
!
   use common_prob_data
   use data_structure3D
   use element_data
   use parameters
   use mpi_param, only: RANK
!
   implicit none
!
   integer   , intent(in)  :: Mdle
   integer   , intent(in)  :: Nrdof
   complex(8), intent(out) :: Zaloc(Nrdof,Nrdof), Zbloc(Nrdof)
!
!-------------------------------------------------------------------------
!
!..aux variables
   complex(8) :: zJ(3), za, zb
!
   real(8) :: rjac, wa, weight
   integer :: iflag, nrdofH, nrdofE, nint, k1, k2, l
!
   integer :: etype
!
!..element order, orientation for edges and faces
   integer :: norder(19), norient_edge(12), norient_face(6)
!
!..geometry dof
   real(8) :: xnod(3,MAXbrickH)
!
!..geometry
! [
!  xi     coordinates in reference element
!  x      coordinates in physical domain
!  dxdxi  d(x) / dxi  Jacobian
!  dxidx  d(xi) / dx  Inverse Jacobian
! ]
   real(8) :: xi(3), x(3), dxdxi(3,3), dxidx(3,3)
!
!..H1 shape functions
   real(8) :: shapH(MAXbrickH), gradH(3,MAXbrickH)
!
!..H(curl) shape functions
   real(8) :: shapE(3,MAXbrickE), curlE(3,MAXbrickE)
!
!..3D quadrature data
!  [xiloc: integration points (local coordinates)]
!  [waloc: integration weights per integration point]
   real(8) :: xiloc(3,MAX_NINT3), waloc(MAX_NINT3)
!
!..workspace for trial and test variables
   real(8) :: E(3), CE(3), F(3), CF(3)
!
!----------------------------------------------------------------------
!
   Zaloc = ZERO; Zbloc = ZERO
!
!..element type
   etype = NODES(Mdle)%ntype
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
!----------------------------------------------------------------------
!     E L E M E N T    I N T E G R A L S                              |
!----------------------------------------------------------------------
!
!..set quadrature points and weights
   call set_3D_int(etype,norder,norient_face, nint,xiloc,waloc)
!
!..loop over integration points
   do l=1,nint
!
!  ...coordinates and weight of this integration point
      xi(1:3)=xiloc(1:3,l)
      wa=waloc(l)
!
!  ...H1 shape functions (for geometry)
      call shape3DH(etype,xi,norder,norient_edge,norient_face, nrdofH,shapH,gradH)
!
!  ...H(curl) shape functions
      call shape3DE(etype,xi,norder,norient_edge,norient_face, nrdofE,shapE,curlE)
!
!  ...geometry map
      call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH, x,dxdxi,dxidx,rjac,iflag)
!
!  ...integration weight
      weight = rjac*wa
!
!  ...get the RHS
      call getf(Mdle,x, zJ)
!
!  ...loop through H(curl) test functions
      do k1=1,nrdofE
!
!     ...Piola transformation
         F(1:3) = shapE(1,k1)*dxidx(1,1:3) &
                + shapE(2,k1)*dxidx(2,1:3) &
                + shapE(3,k1)*dxidx(3,1:3)
         CF(1:3) = dxdxi(1:3,1)*curlE(1,k1) &
                 + dxdxi(1:3,2)*curlE(2,k1) &
                 + dxdxi(1:3,3)*curlE(3,k1)
         CF(1:3) = CF(1:3)/rjac
!
!     ...accumulate for the load vector: (-iωJ,F)
         za = F(1)*zJ(1) + F(2)*zJ(2) + F(3)*zJ(3)
         Zbloc(k1) = Zbloc(k1) - ZI*OMEGA*za*weight
!
!     ...loop through H(curl) trial functions
         do k2=1,nrdofE
!
!        ...Piola transformation
            E(1:3) = shapE(1,k2)*dxidx(1,1:3) &
                   + shapE(2,k2)*dxidx(2,1:3) &
                   + shapE(3,k2)*dxidx(3,1:3)
            CE(1:3) = dxdxi(1:3,1)*curlE(1,k2) &
                    + dxdxi(1:3,2)*curlE(2,k2) &
                    + dxdxi(1:3,3)*curlE(3,k2)
            CE(1:3) = CE(1:3)/rjac
!
!        ...accumulate for the stiffness matrix: ((1/μ)curl E, curl F)-((ω^2ε-iωσ)E, F)
            za = (CE(1)*CF(1) + CE(2)*CF(2) + CE(3)*CF(3)) / MU
            zb = (OMEGA*OMEGA*EPSILON - ZI*OMEGA*SIGMA)*(E(1)*F(1) + E(2)*F(2) + E(3)*F(3))
            Zaloc(k1,k2) = Zaloc(k1,k2) + (za-zb)*weight
!
!     ...end of loop through H(curl) trial functions
         enddo
!  ...end of loop through H(curl) test functions
      enddo
!..end of loop through integration points
   enddo
!
!
end subroutine elem_maxwell


