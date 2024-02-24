!----------------------------------------------------------------------
!
!     routine name      - elem
!
!----------------------------------------------------------------------
!
!     latest revision:  - May 2021
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
   use physics  , only: NR_PHYSA
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
!..activate one physics variable (H1) for assembly
   Itest(1) = 1; Itrial(1) = 1
!
!..determine leading dimension of ALOC(1,1)%array
   nn = ubound(ALOC(1,1)%array,1)
!
!..call element integration routine
   call elem_vect_poisson(Mdle,nn, ALOC(1,1)%array,BLOC(1)%array)
!
end subroutine elem
!
!-------------------------------------------------------------------------
!
!     routine name      - elem_vect_poisson
!
!-------------------------------------------------------------------------
!
!     latest revision:  - June 2021
!
!     purpose:          - compute element stiffness and load
!
!     arguments:
!        in:
!             Mdle      - an element middle node number, identified
!                         with the element
!             Nn        - leading dimension in Zaloc and Zbloc
!        out:
!             Zaloc     - element stiffness matrix
!             Zbloc     - element load vector
!
!-------------------------------------------------------------------------
!
subroutine elem_vect_poisson(Mdle,Nn, Zaloc,Zbloc)
!
   use data_structure3D
   use element_data
   use parameters
   use mpi_param, only: RANK
!
   implicit none
!
   integer, intent(in)  :: Mdle
   integer, intent(in)  :: Nn
   real(8), intent(out) :: Zaloc(Nn,Nn), Zbloc(Nn)
!
!..aux variables
   real(8) :: rjac, fval, gval, wa, weight, q, p
   integer :: iflag, nrv, nre, nrf
   integer :: nrdofH,nint,k1,k2,l,ivar,n1,n2
!
!----------------------------------------------------------------------
!
   integer :: etype,ftype
!
!..element order, orientation for edges and faces
   integer :: norder(19), norient_edge(12), norient_face(6)
!
!..face order
   integer :: norderf(5)
!
!..BC flags
   integer :: ibc(6,NRINDEX)
!
!..geometry dof
   real(8) :: xnod(3,MAXbrickH)
!
!..geometry
   real(8) :: xi(3), dxidt(3,2), x(3), dxdxi(3,3), dxidx(3,3), dxdt(3,2)
   real(8) :: rn(3), t(2),bjac
!
!..H1 shape functions
   real(8) :: shapH(MAXbrickH), gradH(3,MAXbrickH)
!
!..3D quadrature data
   real(8) :: xiloc(3,MAX_NINT3), waloc(MAX_NINT3)
!
!..2D quadrature data
   real(8) :: tloc(2,MAX_NINT2), wtloc(MAX_NINT2)
!
!..workspace for trial and test variables
   real(8) :: dq(3), dp(1:3)
!
   integer :: ifc,nsign
!
#if DEBUG_MODE
   integer :: ivar1,ivar2
   integer :: iprint
   iprint=0
#endif
!
!----------------------------------------------------------------------
!
   Zaloc = ZERO; Zbloc = ZERO
!
!..element type
   etype = NODES(Mdle)%ntype
!..[number of vertices, edges, and faces of this element (type)]
   nrv = nvert(etype)
   nre = nedge(etype)
   nrf = nface(etype)
!
!..determine order of approximation
   !write(*,2050) '[', RANK, '] FIND ORDER'; call pause
   call find_order(Mdle, norder)
!
#if DEBUG_MODE
   if (iprint.eq.1) then
     write(*,7010) norder(1:nre+nrf+1)
7010 format('elem_vect_poisson: norder = ',19i4)
   endif
#endif
!
!..determine edge and face orientations
   !write(*,2050) '[', RANK, '] FIND ORIENT'; call pause
   call find_orient(Mdle, norient_edge,norient_face)
!
!..determine BC flags
   call find_bc(Mdle, ibc)
!
#if DEBUG_MODE
   if (iprint.eq.1) then
      do ivar=1,3
         write(*,7020) ivar, ibc(1:6,ivar)
 7020    format('elem_vect_poisson: BC flags for ivar = ',i1,2x,6i2)
      enddo
   endif
#endif
!
!..determine nodes coordinates
   !write(*,2050) '[', RANK, '] FIND NODCOR'; call pause
   call nodcor(Mdle, xnod)
!
!----------------------------------------------------------------------
!     E L E M E N T    I N T E G R A L S                              |
!----------------------------------------------------------------------
!
!..use the order to set the quadrature
   call set_3D_int(etype,norder,norient_face, nint,xiloc,waloc)
!
!..loop over integration points
   do l=1,nint
!
!  ...coordinates and weight of this integration point
      xi(1:3)=xiloc(1:3,l);
      wa=waloc(l)
!
!  ...H1 shape functions (for geometry)
      call shape3DH(etype,xi,norder,norient_edge,norient_face, nrdofH,shapH,gradH)
!
!  ...geometry map
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
         do ivar=1,3
            n1 = (k1-1)*3 + ivar
            Zbloc(n1) = Zbloc(n1) + q*fval*weight
         enddo
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
            do ivar=1,3
               n1 = (k1-1)*3 + ivar
               n2 = (k2-1)*3 + ivar
               Zaloc(n1,n2) = Zaloc(n1,n2) + weight * ( &
                              dq(1)*dp(1) + dq(2)*dp(2) + dq(3)*dp(3))
            enddo
!
!     ...end of loop through H1 trial functions
         enddo
!  ...end of loop through H1 test functions
      enddo
!..end of loop through integration points
   enddo
!
#if DEBUG_MODE
      if (iprint.eq.1) then
  100   write(*,*) 'elem_vect_possion: SET ivar1,ivar2 (1:3)'
        read(*,*) ivar1,ivar2
        if (ivar1.eq.0) goto 200
        write(*,7100) ivar1,ivar2
 7100   format(    '                   ivar1,ivar2 = '2i2)
        do k1=1,nrdofH
          n1 = (k1-1)*3+ivar1
          write(*,7200) k1
 7200     format(' k1 = ',i3)
          write(*,8100) (Zaloc(n1,(k2-1)*3+ivar2),k2=1,nrdofH)
 8100     format(16e12.5)
        enddo
        goto 100
      endif
  200 continue
#endif
!
!---------------------------------------------------------------------
!    B O U N D A R Y    I N T E G R A L S                            |
!---------------------------------------------------------------------
!
!..loop through element faces
   do ifc=1,nrf
!
!  ...sign factor to determine the OUTWARD normal unit vector
      nsign = nsign_param(etype,ifc)
!
!  ...face type
      ftype = face_type(etype,ifc)
!
!  ...face order of approximation
      call face_order(etype,ifc,norder, norderf)
!
!  ...set 2D quadrature
      call set_2D_int(ftype,norderf,norient_face(ifc), nint,tloc,wtloc)
!
!  ...loop through integration points
      do l=1,nint
!
!     ...face coordinates
         t(1:2) = tloc(1:2,l)
!
!     ...face parametrization
         call face_param(etype,ifc,t, xi,dxidt)
!
!     ...determine element H1 shape functions
         call shape3DH(etype,xi,norder,norient_edge,norient_face, &
                       nrdofH,shapH,gradH)
!
!     ...geometry
         call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign, &
                      x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
!
!     ...get the surface load (for a single eqn)
         call getg(Mdle,x,rn, gval)
!
!     ...integration weight
         weight = bjac*wtloc(l)
!
!     ...loop through H1 test functions
         do k1=1,NrdofH
            q = shapH(k1)
            do ivar=1,3
               n1 = (k1-1)*3+ivar
!
!           ...check if BC flag for the Neumann BC
               if (ibc(ifc,ivar).eq.2) then
                  Zbloc(n1) =  Zbloc(n1) + gval*q*weight
               endif
            enddo
         enddo
!  ...end loop through integration points
      enddo
!
!..end loop through element faces
   enddo
!
end subroutine elem_vect_poisson
