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
!..ALOC: holds local element stiffness matrices
!..BLOC: holds local element load vectors
   use data_structure3D
   use element_data
   use parameters
!
   implicit none
!
   integer, intent(in)  :: Mdle
   integer, intent(in)  :: Nrdof
   real(8), intent(out) :: Zaloc(Nrdof,Nrdof), Zbloc(Nrdof)
!
!-------------------------------------------------------------------------
!
!..aux variables
   real(8) :: rjac, fval, wa, weight, q, dq(3)
   integer :: etype, iflag, nrdofH, nint, noff, nda, k, l
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
!..3D quadrature data
   real(8) :: xiloc(3,MAX_NINT3), waloc(MAX_NINT3)
!
!..workspace for auxiliary matrix (storing info at integration points)
   real(8), allocatable :: A_TEST(:,:)
!
!-------------------------------------------------------------------------
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
!
!..loop over integration points
   do l=1,nint
!
!  ...coordinates and weight of this integration point
      xi(1:3)=xiloc(1:3,l)
      wa=waloc(l)
!
!  ...H1 shape functions
      call shape3DH(etype,xi,norder,norient_edge,norient_face, nrdofH,shapH,gradH)
!
!  ...geometry map
      call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH, x,dxdxi,dxidx,rjac,iflag)
!
!  ...integration weight
      weight = rjac*wa
!
!  ...get the RHS
      call getf(Mdle,x, fval)
!
!  ...offset in auxiliary matrix
      noff = 3*(l-1)
!
!  ...loop through H1 test functions
      do k=1,nrdofH
!
!     ...Piola transformation
         q = shapH(k)
         dq(1:3) = gradH(1,k)*dxidx(1,1:3) &
                 + gradH(2,k)*dxidx(2,1:3) &
                 + gradH(3,k)*dxidx(3,1:3)
!
!     ...accumulate for the load vector
         Zbloc(k) = Zbloc(k) + q*fval*weight
!
!     ...fill auxiliary matrix for stiffness
         A_TEST(k,noff+1:noff+3) = dq(1:3) * sqrt(weight)
!
!  ...end of loop through H1 test functions
      enddo
!
!..end of loop through integration points
   enddo
!
!..compute stiffness matrix
   call DSYRK('U','N',Nrdof,nda,1.0d0,A_TEST(1:Nrdof,1:nda) ,Nrdof,  &
                                0.0d0,Zaloc(1:Nrdof,1:Nrdof),Nrdof)
!
!..fill lower triangular part of symmetric stiffness matrix
   do k=1,Nrdof
      Zaloc(k+1:Nrdof,k) = Zaloc(k,k+1:Nrdof)
   enddo
!
   deallocate(A_TEST)
!
end subroutine elem_opt

