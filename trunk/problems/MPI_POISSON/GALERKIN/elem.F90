!----------------------------------------------------------------------
!
!     routine name      - elem
!
!----------------------------------------------------------------------
!
!     latest revision:  - Aug 2019
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
!----------------------------------------------------------------------
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
!----------------------------------------------------------------------
!
!..activate one physics variable (H1) for assembly
   Itest(1) = 1; Itrial(1) = 1
!
!..determine order of approximation
   call find_order(Mdle, norder)
!
!..find number of dof for each energy space supported by the element
   call celndof(NODES(Mdle)%type,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
!
!..call element integration routine
   call elem_poisson(Mdle,nrdofH, ALOC(1,1)%array,BLOC(1)%array)
!
end subroutine elem
!
!-------------------------------------------------------------------------
!
!     routine name      - elem_poisson
!
!-------------------------------------------------------------------------
!
!     latest revision:  - July 2019
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
subroutine elem_poisson(Mdle,Nrdof, Zaloc,Zbloc)
!
!..ALOC: holds local element stiffness matrices
!..BLOC: holds local element load vectors
   use assembly, only: ALOC,BLOC
   use control , only: INTEGRATION
   use physics , only: NR_PHYSA
   use data_structure3D
   use element_data
   use parameters
   use mpi_param, only: RANK
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
   real(8) :: rjac, fval, wa, weight, q, p
   integer :: iflag, nrv, nre, nrf
   integer :: nrdofH, nint, k1, k2, l
!
   character(len=4) :: etype
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
!..3D quadrature data
!  [xiloc: integration points (local coordinates)]
!  [waloc: integration weights per integration point]
   real(8) :: xiloc(3,MAX_NINT3), waloc(MAX_NINT3)
!
!..workspace for trial and test variables
   real(8) :: dq(3), u(3), dp(1:3)
!
!-------------------------------------------------------------------------
!
   Zaloc = ZERO; Zbloc = ZERO
!
!..element type
   etype = NODES(Mdle)%type
!  [number of vertices, edges, and faces of this element (type)]
   nrv = nvert(etype)
   nre = nedge(etype)
   nrf = nface(etype)
!
!..determine order of approximation
   !write(*,2050) '[', RANK, '] FIND ORDER'; call pause
   call find_order(Mdle, norder)
!
!..determine edge and face orientations
   !write(*,2050) '[', RANK, '] FIND ORIENT'; call pause
   call find_orient(Mdle, norient_edge,norient_face)
!
!..determine nodes coordinates
   !write(*,2050) '[', RANK, '] FIND NODCOR'; call pause
   call nodcor(Mdle, xnod)
!
!----------------------------------------------------------------------
!     E L E M E N T    I N T E G R A L S                              |
!----------------------------------------------------------------------
!
!..set quadrature points and weights
!  [
!  in:
!   Type      - element type
!   Norder    - order of approximation
!
!  out:
!   Nint      - number of integration points
!   Xiloc     - integration points
!   Waloc     - weights
!  ]
   call set_3D_int(etype,norder,norient_face, nint,xiloc,waloc)
!
!..loop over integration points
   do l=1,nint
!
!  ...coordinates and weight of this integration point
      xi(1:3)=xiloc(1:3,l)
      wa=waloc(l)
!
! H1 shape functions (for geometry)
! [
!     in:
!       Type           - element type
!       Xi             - master element coordinates
!       Norder         - polynomial order for the nodes (H1 sense)
!       Nedge_orient   - edge orientations
!       Nface_orient   - face orientations
!     out:
!       NrdofH         - number of the element shape functions
!       ShapH          - values of shape functions
!       GradH          - values of derivatives of the shape functions
!                        wrt to master element coordinates
!
! ]
      call shape3DH(etype,xi,norder,norient_edge,norient_face, nrdofH,shapH,gradH)
!
! geometry map
! [
!     in:
!             Mdle     - element middle node
!             Xi       - master element coordinates of a point
!             Xnod     - element geometry dof
!             ShapH    - values of shape functions at a point
!             GradH    - derivatives of the shape functions wrt
!                        master coordinates at the point
!             NrdofH   - number of element H1 dof
!     out:
!             X        - physical coordinates of the point
!             Dxdxi    - Jacobian matrix
!             Dxidx    - inverse Jacobian matrix
!             Rjac     - jacobian (determinant of the Jacobian matrix)
!             Iflag    = 0 OK
!                        1 negative jacobian
! ]
      call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH, x,dxdxi,dxidx,rjac,iflag)
!
!  ...integration weight
      weight = rjac*wa
!
!  ...get the RHS
      call getf(Mdle,x, fval)
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
         Zbloc(k1) = Zbloc(k1) + q*fval*weight
!
!     ...loop through H1 trial functions
         do k2=1,nrdofH
!        ...Piola transformation
            p = shapH(k2)
            dp(1:3) = gradH(1,k2)*dxidx(1,1:3) &
                    + gradH(2,k2)*dxidx(2,1:3) &
                    + gradH(3,k2)*dxidx(3,1:3)
!
!        ...accumulate for the stiffness matrix
            Zaloc(k1,k2) = Zaloc(k1,k2) + weight * ( &
                           dq(1)*dp(1) + dq(2)*dp(2) + dq(3)*dp(3))
!
!     ...end of loop through H1 trial functions
         enddo
!  ...end of loop through H1 test functions
      enddo
!..end of loop through integration points
   enddo
!
!
end subroutine elem_poisson


