!----------------------------------------------------------------------
!                                                                     
!     routine name      - elem
!                                                                     
!----------------------------------------------------------------------
!                                                                     
!     latest revision:  - July 2019
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
!
!..select node%case (see data_structure3D.F)
!..[case = 2^NR_PHYSA-1]
   select case (NODES(Mdle)%case)
!  ...adjust case to problem
!  ...to support multiple physics on different elements
      case (1)
         Itest(1:NR_PHYSA) = 1
         Itrial(1:NR_PHYSA) = 1
         call elem_poisson(Mdle)
      case default
         write(*,*) 'elem: Mdle,NODES(Mdle)%case = ', &
                     Mdle,NODES(Mdle)%case
         call logic_error(ERR_INVALID_VALUE, __FILE__, __LINE__)
   end select
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
!
!-------------------------------------------------------------------------
!
subroutine elem_poisson(Mdle)
!
!..ALOC: holds local element stiffness matrices
!..BLOC: holds local element load vectors
   use assembly, only: ALOC,BLOC, NR_RHS
   use control , only: INTEGRATION
   use physics , only: NR_PHYSA
   use data_structure3D
   use element_data
   use parameters
   use mpi_param, only: RANK
!
   implicit none
!
   integer, intent(in) :: Mdle
!
!..aux variables
   real(8) :: rjac, fval, wa, weight, q, p
   integer :: iflag, nrv, nre, nrf
   integer :: nrdofH, nrdofE, nrdofV, nrdofQ, nint, k1, k2, l
!
!----------------------------------------------------------------------
!
   character(len=4) :: etype,ftype
!
!..element order, orientation for edges and faces
   integer :: norder(19), norient_edge(12), norient_face(6)
!
!..face order
   integer :: norderf(5)
!
!..geometry dof
   real(8) :: xnod(3,MAXbrickH)
!
!..geometry
! [
!  xi     coordinates in reference element
!  dxidt  d(xi) / dt
!  x      coordinates in physical domain
!  dxdxi  d(x) / dxi  Jacobian
!  dxidx  d(xi) / dx  Inverse Jacobian
!  dxdt   d(x) / dt
!  rt     .
!  rn     normal vec
!  t      .
! ]
   real(8) :: xi(3), dxidt(3,2), x(3), dxdxi(3,3), dxidx(3,3), dxdt(3,2)
   real(8) :: rt(3,2), rn(3), t(2)
!
!..H1 shape functions
   real(8) :: shapH(MAXbrickH), gradH(3,MAXbrickH)
!
!..3D quadrature data
!  [xiloc: integration points (local coordinates)]
!  [waloc: integration weights per integration point]
   real(8) :: xiloc(3,MAX_NINT3), waloc(MAX_NINT3)
!
!..2D quadrature data
   real(8) :: tloc(2,MAX_NINT2), wtloc(MAX_NINT2)
!
!..BC's flags
!  [for each attribute and face (max #faces is 6 for 3D elements]
   integer :: ibc(6,NR_PHYSA)
!
!..workspace for trial and test variables
   real(8) :: dq(3), u(3), dp(1:3), v(3), vec(3)
!
   real(8), allocatable :: Zaloc(:,:), Zbloc(:)
!
!----------------------------------------------------------------------
!
!..element type
   etype = NODES(Mdle)%type
!  [number of vertices, edges, and faces of this element (type)]
   nrv = nvert(etype)
   nre = nedge(etype)
   nrf = nface(etype)
!
   !write(*,2050) '[', RANK, '] FIND ORDER'; call pause
!..determine order of approximation
   call find_order(Mdle, norder)
!
   !write(*,2050) '[', RANK, '] FIND ORIENT'; call pause
!..determine edge and face orientations
   call find_orient(Mdle, norient_edge, norient_face)
!
!..determine nodes coordinates
   !write(*,2050) '[', RANK, '] FIND NODCOR'; call pause
   call nodcor(Mdle, xnod)
!
!..get the element boundary conditions flags
   !write(*,2050) '[', RANK, '] FIND BC'; call pause
   call find_bc(Mdle, ibc)
!
!..find number of trial dof for each energy space supported by the element
   !write(*,2050) '[', RANK, '] FIND CELNDOF'; call pause
   call celndof(NODES(Mdle)%type,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
!
!..allocate space for matrices
   allocate(Zaloc(nrdofH,nrdofH)); Zaloc = ZERO
   allocate(Zbloc(nrdofH))       ; Zbloc = ZERO
!
   !write(*,2050) '[', RANK, '] START ELEM INTEGRALS'
   !2050 format(A,I2,A)
!
!----------------------------------------------------------------------
!     E L E M E N T    I N T E G R A L S                              |
!----------------------------------------------------------------------
!
!..use the enriched order to set the quadrature
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
   call set_3Dint(etype, norder, nint, xiloc, waloc)
!
!..loop over integration points
   do l=1,nint
!
!  ...coordinates and weight of this integration point
      xi(1:3)=xiloc(1:3,l);
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
      call shape3H(etype,xi,norder,norient_edge,norient_face, nrdofH,shapH,gradH)
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
      call getf(Mdle, x,fval)
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
!     ...second loop through H1 trial functions
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
   BLOC(1)%array(1:nrdofH,1) = Zbloc(1:nrdofH)
   ALOC(1,1)%array(1:nrdofH,1:nrdofH) = Zaloc(1:nrdofH,1:nrdofH)
!
   deallocate(Zaloc, Zbloc)
!
!
end subroutine elem_poisson


