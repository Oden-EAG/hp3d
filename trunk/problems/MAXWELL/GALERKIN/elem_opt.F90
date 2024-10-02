!-------------------------------------------------------------------------
!
!     routine name      - elem_opt
!
!-------------------------------------------------------------------------
!
!> @date   Apr 2024
!
!> @brief  Compute element stiffness and load using BLAS3
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
subroutine elem_opt(Mdle,Nrdof, Zaloc,Zbloc)
!
   use common_prob_data
   use data_structure3D
   use element_data
   use parameters
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
   real(8)    :: rjac, wa, weight
   integer    :: etype, iflag, nrdofH, nrdofE, nint, noff, nda, k, l
!
!..element order, orientation for edges and faces
   integer :: norder(19), norient_edge(12), norient_face(6)
!
!..geometry dof
   real(8) :: xnod(3,MAXbrickH)
!
!..geometry
   real(8) :: xi(3), x(3), dxdxi(3,3), dxidx(3,3)
!
!..H1 shape functions
   real(8) :: shapH(MAXbrickH), gradH(3,MAXbrickH)
!
!..H(curl) shape functions
   real(8) :: shapE(3,MAXbrickE), curlE(3,MAXbrickE)
!
!..3D quadrature data
   real(8) :: xiloc(3,MAX_NINT3), waloc(MAX_NINT3)
!
!..workspace for trial and test variables
   real(8) :: F(3), CF(3)
!
!..workspace for auxiliary matrices (storing info at integration points)
   complex(8), allocatable :: M_TEST(:,:), A_TEST(:,:)
!
!----------------------------------------------------------------------
!
   Zbloc = ZERO
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
!..allocate auxiliary matrix
   nda = 3*nint
   allocate(A_TEST(Nrdof,nda))
   allocate(M_TEST(Nrdof,nda))
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
!  ...offset in auxiliary matrices
      noff = 3*(l-1)
!
!  ...loop through H(curl) test functions
      do k=1,nrdofE
!
!     ...Piola transformation
         F(1:3) = shapE(1,k)*dxidx(1,1:3) &
                + shapE(2,k)*dxidx(2,1:3) &
                + shapE(3,k)*dxidx(3,1:3)
         CF(1:3) = dxdxi(1:3,1)*curlE(1,k) &
                 + dxdxi(1:3,2)*curlE(2,k) &
                 + dxdxi(1:3,3)*curlE(3,k)
         CF(1:3) = CF(1:3)/rjac
!
!     ...accumulate for the load vector: (-iωJ,F)
         za = F(1)*zJ(1) + F(2)*zJ(2) + F(3)*zJ(3)
         Zbloc(k) = Zbloc(k) - ZI*OMEGA*za*weight
!
!     ...fill auxiliary matrices for stiffness
         zb = OMEGA*OMEGA*EPSILON - ZI*OMEGA*SIGMA
         M_TEST(k,noff+1:noff+3) = F(1:3) * sqrt(zb*weight)
         A_TEST(k,noff+1:noff+3) = CF(1:3) * sqrt(weight / MU)
!
!  ...end of loop through H(curl) test functions
      enddo
!..end of loop through integration points
   enddo
!
!..compute stiffness matrix
   call ZSYRK('U','N',Nrdof,nda, ZONE,A_TEST(1:Nrdof,1:nda) ,Nrdof,  &
                                 ZERO,Zaloc(1:Nrdof,1:Nrdof),Nrdof)
   call ZSYRK('U','N',Nrdof,nda,-ZONE,M_TEST(1:Nrdof,1:nda) ,Nrdof,  &
                                 ZONE,Zaloc(1:Nrdof,1:Nrdof),Nrdof)
!
!..fill lower triangular part of symmetric stiffness matrix
   do k=1,Nrdof
      Zaloc(k+1:Nrdof,k) = Zaloc(k,k+1:Nrdof)
   enddo
!
   deallocate(A_TEST,M_TEST)
!
end subroutine elem_opt
