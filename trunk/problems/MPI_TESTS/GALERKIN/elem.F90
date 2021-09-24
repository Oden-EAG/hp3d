!----------------------------------------------------------------------
!
!     routine name      - elem
!
!----------------------------------------------------------------------
!
!     latest revision:  - Jul 21
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
!
      subroutine elem(Mdle, Itest,Itrial)
!
      use data_structure3D
      use physics
      use assembly , only: ALOC,BLOC
!
      implicit none
!
      integer,                    intent(in)  :: Mdle
      integer,dimension(NR_PHYSA),intent(out) :: Itest,Itrial
!
      integer :: nn
!
!----------------------------------------------------------------------
!
!  ...activate just one physics variable for assembly
      Itest = 0; Itrial = 0
!
!  ...determine leading dimension of ALOC(1,1)%array
      nn = ubound(ALOC(1,1)%array,1)
!
!  ...call the actual element routine
!
!  ...H1 projection
      if (PHYSAm(1)) then
        Itest(1) = 1; Itrial(1) = 1
        nn = ubound(ALOC(1,1)%array,1)
        call elem_proj('contin',Mdle,nn, ALOC(1,1)%array,BLOC(1)%array)
      elseif (PHYSAm(2)) then
        Itest(2) = 1; Itrial(2) = 1
        nn = ubound(ALOC(2,2)%array,1)
        call elem_proj('tangen',Mdle,nn, ALOC(2,2)%array,BLOC(2)%array)
      elseif (PHYSAm(3)) then
        Itest(3) = 1; Itrial(3) = 1
        nn = ubound(ALOC(3,3)%array,1)
        call elem_proj('normal',Mdle,nn, ALOC(3,3)%array,BLOC(3)%array)
      elseif (PHYSAm(4)) then
        Itest(4) = 1; Itrial(4) = 1
        nn = ubound(ALOC(4,4)%array,1)
        call elem_proj('discon',Mdle,nn, ALOC(4,4)%array,BLOC(4)%array)
      endif
!
      end subroutine elem
!
!-------------------------------------------------------------------------
!
!     routine name      - elem_proj
!
!-------------------------------------------------------------------------
!
!     latest revision:  - Jul 21
!
!     purpose:          - compute element stiffness and load for 
!                         a projection problem
!
!     arguments:
!        in:
!             Attr      = 'contin'  H1-projection
!                         'tangen'  H(curl)-projection
!                         'normal'  H(div)-projection
!                         'discon'  L2-projection
!             Mdle      - an element middle node number, identified
!                         with the element
!             Nn        - leading dimension in Zaloc and Zbloc
!        out:
!             Zaloc     - element stiffness matrix
!             Zbloc     - element load vector
!
!-------------------------------------------------------------------------
!
      subroutine elem_proj(Attr,Mdle,Nn, Zaloc,Zbloc)
!
      use control , only: INTEGRATION
      use data_structure3D
      use element_data
      use parameters
      use mpi_param, only: RANK
!
      implicit none
!
      character(len=6), intent(in)  :: Attr
      integer, intent(in)  :: Mdle
      integer, intent(in)  :: Nn
      real(8), intent(out) :: Zaloc(Nn,Nn), Zbloc(Nn)
!
!  ...aux variables
      real(8) :: wa, weight
      integer :: nrv,nre,nrf,iflag
      integer :: nrdofH,nrdofE,nrdofV,nrdofQ,nint,k1,k2,l
!
      character(len=4) :: etype
!
!  ...element order, orientation for edges and faces
      integer :: norder(19), norient_edge(12), norient_face(6)
!
!  ...geometry dof
      real(8) :: xnod(3,MAXbrickH)
!
!  ...geometry
      real(8) :: xi(3),x(3),dxdxi(3,3),dxidx(3,3),rjac
!
!  ...shape functions
      real(8) :: shapH(MAXbrickH), gradH(3,MAXbrickH)
      real(8) :: shapE(3,MAXbrickE), curlE(3,MAXbrickE)
      real(8) :: shapV(3,MAXbrickV), divV(MAXbrickV)
      real(8) :: shapQ(MAXbrickQ)
!
!  ...3D quadrature data
      real(8) :: xiloc(3,MAX_NINT3), waloc(MAX_NINT3)
!
!  ...exact solution
      real(8),dimension(  MAXEQNH    ) ::   valH
      real(8),dimension(  MAXEQNH,3  ) ::  dvalH
      real(8),dimension(  MAXEQNH,3,3) :: d2valH
      real(8),dimension(3,MAXEQNE    ) ::   valE
      real(8),dimension(3,MAXEQNE,3  ) ::  dvalE
      real(8),dimension(3,MAXEQNE,3,3) :: d2valE
      real(8),dimension(3,MAXEQNV    ) ::   valV
      real(8),dimension(3,MAXEQNV,3  ) ::  dvalV
      real(8),dimension(3,MAXEQNV,3,3) :: d2valV
      real(8),dimension(  MAXEQNQ    ) ::   valQ
      real(8),dimension(  MAXEQNQ,3  ) ::  dvalQ
      real(8),dimension(  MAXEQNQ,3,3) :: d2valQ
!
!  ...workspace for trial and test variables
      real(8) :: q,dq(3), p,dp(3), &
                 F(3),curlF(3), G(3),curlG(3), &
                 w(3),divw, u(3),divu
!
#if DEBUG_MODE
      integer :: iprint = 0
#endif
!
!----------------------------------------------------------------------
!  
      Zaloc(1:Nn,1:Nn) = .0d0; Zbloc(1:Nn) = .0d0
!
!  ...element type
      etype = NODES(Mdle)%type
!
!  ...[number of vertices, edges, and faces of this element (type)]
      nrv = nvert(etype)
      nre = nedge(etype)
      nrf = nface(etype)
!
!  ...determine order of approximation
      call find_order(Mdle, norder)
!
#if DEBUG_MODE
      if (iprint.eq.1) then
        write(*,7010) norder(1:nre+nrf+1)
 7010   format('elem_proj: norder = ',19i4)
      endif
#endif
!
!  ...determine edge and face orientations
      call find_orient(Mdle, norient_edge,norient_face)
!
!  ...determine nodes coordinates
      call nodcor(Mdle, xnod)
!
!----------------------------------------------------------------------
!     E L E M E N T    I N T E G R A L S                              |
!----------------------------------------------------------------------
!
!  ...use the order to set the quadrature
      call set_3D_int(etype,norder,norient_face, nint,xiloc,waloc)
!
!  ...loop over integration points
      do l=1,nint
!
!  .....coordinates and weight of this integration point
        xi(1:3) = xiloc(1:3,l)
        wa=waloc(l)
!
!  .....H1 shape functions (for geometry and H1-projection)
        call shape3DH(etype,xi,norder,norient_edge,norient_face, nrdofH,shapH,gradH)
        select case(Attr)
!
!  .....H(curl) shape functions 
        case('tangen')
          call shape3DE(etype,xi,norder,norient_edge,norient_face, nrdofE,shapE,curlE)
!
!  .....H(div) shape functions 
        case('normal')
          call shape3DV(etype,xi,norder,norient_face, nrdofV,shapV,divV)
!
!  .....L2 shape functions 
        case('discon')
          call shape3DQ(etype,xi,norder, nrdofQ,shapQ)
        end select
!
!  .....geometry map
        call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH, x,dxdxi,dxidx,rjac,iflag)
!
!  .....integration weight
        weight = rjac*wa
!
!  .....get the exact solution
        call exact(x,Mdle, valH,dvalH,d2valH, valE,dvalE,d2valE, &
                           valV,dvalV,d2valV, valQ,dvalQ,d2valQ)
!
        select case(Attr)
!
!  .....H1 projection
        case('contin')
!
!  .......loop through H1 test functions
          do k1=1,nrdofH
!
!  .........Piola transformation
            q = shapH(k1)
            dq(1:3) = gradH(1,k1)*dxidx(1,1:3) &
                    + gradH(2,k1)*dxidx(2,1:3) &
                    + gradH(3,k1)*dxidx(3,1:3)
!
!  .........accumulate for the load vector
            Zbloc(k1) = Zbloc(k1) &
                      + (dvalH(1,1)*dq(1)+dvalH(1,2)*dq(2)+dvalH(1,3)*dq(3) &
                        + valH(1)*q)*weight
!
!  .........loop through H1 trial functions
            do k2=1,nrdofH
!
!  ...........Piola transformation
              p = shapH(k2)
              dp(1:3) = gradH(1,k2)*dxidx(1,1:3) &
                      + gradH(2,k2)*dxidx(2,1:3) &
                      + gradH(3,k2)*dxidx(3,1:3)
!
!  ...........accumulate for the stiffness matrix
              Zaloc(k1,k2) = Zaloc(k1,k2) &
                           +(dq(1)*dp(1)+dq(2)*dp(2)+dq(3)*dp(3) &
                           + q*p)*weight
!
!  .........end of loop through H1 trial functions
            enddo
!
!  .......end of loop through H1 test functions
          enddo
!
!  .....H(curl) projection
        case('tangen')
!
!  .......loop through H(curl) test functions
          do k1=1,nrdofE
!
!  .........Piola transformation
            F(1:3) = shapE(1,k1)*dxidx(1,1:3) &
                   + shapE(2,k1)*dxidx(2,1:3) &
                   + shapE(3,k1)*dxidx(3,1:3)
            curlF(1:3) = (dxdxi(1:3,1)*curlE(1,k1) &
                        + dxdxi(1:3,2)*curlE(2,k1) &
                        + dxdxi(1:3,3)*curlE(3,k1))/rjac
!
!  .........accumulate for the load vector
            Zbloc(k1) = Zbloc(k1) &
                      + ((dvalE(3,1,2)-dvalE(2,1,3))*curlF(1) &
                       + (dvalE(1,1,3)-dvalE(3,1,1))*curlF(2) &
                       + (dvalE(2,1,1)-dvalE(1,1,2))*curlF(3) &
                       + valE(1,1)*F(1)+valE(2,1)*F(2)+valE(3,1)*F(3))*weight
!
!  .........loop through H(curl) trial functions
            do k2=1,nrdofE
!
!  ...........Piola transformation
              G(1:3) = shapE(1,k2)*dxidx(1,1:3) &
                     + shapE(2,k2)*dxidx(2,1:3) &
                     + shapE(3,k2)*dxidx(3,1:3)
              curlG(1:3) = (dxdxi(1:3,1)*curlE(1,k2) &
                          + dxdxi(1:3,2)*curlE(2,k2) &
                          + dxdxi(1:3,3)*curlE(3,k2))/rjac
!
!  ...........accumulate for the stiffness matrix
              Zaloc(k1,k2) = Zaloc(k1,k2) &
                           +(curlF(1)*curlG(1)+curlF(2)*curlG(2)+curlF(3)*curlG(3) &
                           + F(1)*G(1)+F(2)*G(2)+F(3)*G(3))*weight
!
!  .........end of loop through H(curl) trial functions
            enddo
!
!  .......end of loop through H(curl) test functions
          enddo
!
!  .....H(div) projection
        case('normal')
!
!  .......loop through H(div) test functions
          do k1=1,nrdofV
!
!  .........Piola transformation
            w(1:3) = (dxdxi(1:3,1)*shapV(1,k1) &
                    + dxdxi(1:3,2)*shapV(2,k1) &
                    + dxdxi(1:3,3)*shapV(3,k1))/rjac
            divw = divV(k1)/rjac
!
!  .........accumulate for the load vector
            Zbloc(k1) = Zbloc(k1) &
                      + ((dvalv(1,1,1)+dvalv(2,1,2)+dvalv(3,1,3))*divw &
                       + valV(1,1)*w(1)+valV(2,1)*w(2)+valV(3,1)*w(3))*weight
!
!  .........loop through H(div) trial functions
            do k2=1,nrdofV
!
!  ...........Piola transformation
              u(1:3) = (dxdxi(1:3,1)*shapV(1,k2) &
                      + dxdxi(1:3,2)*shapV(2,k2) &
                      + dxdxi(1:3,3)*shapV(3,k2))/rjac
              divu = divV(k2)/rjac
!
!  ...........accumulate for the stiffness matrix
              Zaloc(k1,k2) = Zaloc(k1,k2) &
                           +(divw*divu &
                           + w(1)*u(1)+w(2)*u(2)+w(3)*u(3))*weight
!
!  .........end of loop through H(div) trial functions
            enddo
!
!  .......end of loop through H(div) test functions
          enddo
!
!  .....L2 projection
        case('discon')
!
!  .......loop through L2 test functions
          do k1=1,nrdofQ
!
!  .........Piola transformation
            q = shapQ(k1)/rjac
!
!  .........accumulate for the load vector
            Zbloc(k1) = Zbloc(k1) + valQ(1)*q*weight
!
!  .........loop through L2 trial functions
            do k2=1,nrdofQ
!
!  ...........Piola transformation
              p = shapQ(k2)/rjac
!
!  ...........accumulate for the stiffness matrix
              Zaloc(k1,k2) = Zaloc(k1,k2) &
                           + q*p*weight
!
!  .........end of loop through L2 trial functions
            enddo
!
!  .......end of loop through L2 test functions
          enddo
!
!
        end select
!
!  ...end of loop through integration points
      enddo
!
#if DEBUG_MODE
      if (iprint.ge.1) then
        write(*,*) 'elem_proj: Zbloc, Zaloc = '
        select case(Attr)
        case('contin')
          do k1=1,nrdofH
            write(*,8090) k1,Zbloc(k1)
 8090       format('k1 = ',i4,3x,e12.5)
            write(*,8100) Zaloc(k1,1:nrdofH)
 8100       format(16e12.5)
          enddo
        case('tangen')
          do k1=1,nrdofE
            write(*,8090) k1,Zbloc(k1)
            write(*,8100) Zaloc(k1,1:nrdofE)
          enddo
        case('normal')
          do k1=1,nrdofV
            write(*,8090) k1,Zbloc(k1)
            write(*,8100) Zaloc(k1,1:nrdofV)
          enddo
        case('discon')
          do k1=1,nrdofQ
            write(*,8090) k1,Zbloc(k1)
            write(*,8100) Zaloc(k1,1:nrdofQ)
          enddo
        end select
        call pause
      endif
#endif
!
      end subroutine elem_proj
