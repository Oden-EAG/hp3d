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
   use parametersDPG
   use physics  , only: NR_PHYSA
   use assembly , only: ALOC,BLOC
!
   implicit none
!
   integer,                    intent(in)  :: Mdle
   integer,dimension(NR_PHYSA),intent(out) :: Itest,Itrial
!
!..element order and enriched order
   integer :: norder(19),norderP(19),nordP
!
!..number of test and trial degrees of freedom
   integer :: nrdofH ,nrdofE ,nrdofV ,nrdofQ
   integer :: nrdofHH,nrdofEE,nrdofVV,nrdofQQ
   integer :: nrdofVi_a,nrdofVi_b,nrTest,nrTrial
   integer :: ndofHmdl,ndofEmdl,ndofVmdl,ndofQmdl
   integer :: nrv,nre,nrf
!
!..element type
   character(len=4) :: etype
!
!-------------------------------------------------------------------------
!
!..activate four physics variables (H1 trace + H(div) trace +L2 field + L2 gradient ) for assembly
   Itest(1:NR_PHYSA) = 1; Itrial(1:NR_PHYSA) = 1
!
   norder (1:19) = 0
   norderP(1:19) = 0
!
   etype = NODES(Mdle)%type
   nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
!..determine order of approximation
   call find_order(Mdle, norder)
!
!..set the enriched order of approximation
   select case(etype)
      case('mdlb')
         nordP = NODES(Mdle)%order+NORD_ADD*111
      case('mdlp')
         nordP = NODES(Mdle)%order+NORD_ADD*11
      case('mdln','mdld')
         nordP = NODES(Mdle)%order+NORD_ADD
      case default
         write(*,*) 'elem: invalid etype param. stop.'
         stop
   end select
!
!..note: compute_enriched_order is only provided for hexa
   call compute_enriched_order(etype,nordP, norderP)
!..compute nrdof for trial
   call celndof(etype,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
!..compute nrdof for test
   call celndof(etype,norderP, nrdofHH,nrdofEE,nrdofVV,nrdofQQ)
!..compute number of bubble DOFs
   call ndof_nod(etype,norder(nre+nrf+1), ndofHmdl,ndofEmdl,ndofVmdl,ndofQmdl)
!..compute number of H1 interface dofs
   nrdofVi_a = nrdofH - ndofHmdl
!..calculate number of H(div) interface DOFs
   nrdofVi_b = nrdofV - ndofVmdl
!..calculate total number of trial and test DOFs
   nrTest  = nrdofHH +  nrdofVV
   nrTrial = nrdofQ + 3 * nrdofQ + nrdofVi_a + nrdofVi_b
!    write(*,*) nrdofVi_a,"    ",BLOC(1)%nrow
!    write(*,*) nrdofVi_b,"    ",BLOC(2)%nrow
!    write(*,*) nrdofQ,"    ",BLOC(3)%nrow
!    write(*,*) nrdofQ*3,"    ",BLOC(4)%nrow
! !
!..call element integration routine
   !  write(*,*) Mdle,nrTest,nrdofHH,nrdofVV,nrdofQ,nrdofH,nrdofV,NORD_ADD
   call elem_poisson_UW(Mdle,nrTest,nrTrial,                                                             &
            nrdofHH,nrdofVV,nrdofH,nrdofQ,nrdofQ,3*nrdofQ,nrdofVi_a,nrdofVi_b,                           &
            BLOC(1)%nrow,BLOC(2)%nrow,BLOC(3)%nrow,BLOC(4)%nrow,                                         &
            BLOC(1)%array,ALOC(1,1)%array,ALOC(1,2)%array,ALOC(1,3)%array,ALOC(1,4)%array,               &
            BLOC(2)%array,ALOC(2,1)%array,ALOC(2,2)%array,ALOC(2,3)%array,ALOC(2,4)%array,               &
            BLOC(3)%array,ALOC(3,1)%array,ALOC(3,2)%array,ALOC(3,3)%array,ALOC(3,4)%array,               &
            BLOC(4)%array,ALOC(4,1)%array,ALOC(4,2)%array,ALOC(4,3)%array,ALOC(4,4)%array)
            !     
end subroutine elem
!
!-------------------------------------------------------------------------
!
!     routine name      - elem_poisson
!
!-------------------------------------------------------------------------
!
!     latest revision:  - Oct 2021
!
!     purpose:          - compute element stiffness and load
!
!     arguments:
!        in:
!              Mdle     - element middle node number
!              NrTest   -  number of total test dofs
!              NrTrial  - total number of trial dof
!              NrdofHH  - number of H1 test dof
!              NrdofVV  - number of Hdiv test dof
!              NrdofQ   - number of L2 dofs which will be used for u and sigma_1,2,3 
!              NrdofU   - number of L2 trial dof for u
!              NrdofS   - number of L2 trial dof for sigma  = grad u (Dimension X NrdofQ)
!              NrdofVi_a- number of H1 trial interface dof
!              NrdofVi_b- number of H(div) trial interface dof
!              MdU      - num rows of AlocUU,AlocUS,AlocUVi_a,AlocUVi_b
!              MdS      - num rows of AlocSU,AlocSS,AlocSVi_a,AlocSVi_b
!              MdVi_a   - num rows of AlocVi_aU,AlocVi_aS,AlocVi_aa,Aloc_Viab
!              MdVi_b   - num rows of AlocVi_bU,AlocVi_bS,AlocVi_ba,Aloc_Vibb

!        out:
!              BlocU    - load vectors
!              BlocV
!              BlocVi_a
!              BlocVi_b
!              AlocUU   - stiffness matrices
!              AlocUS
!              AlocUVi_a
!              AlocUVi_b
!              AlocSU
!              AlocSS
!              AlocSVi_a
!              AlocSVi_b
!              AlocVi_aU
!              AlocVi_aS
!              AlocVi_aa
!              AlocVi_ab
!              AlocVi_bU
!              AlocVi_bS
!              AlocVi_ba
!              AlocVi_bb
!
!-------------------------------------------------------------------------
!
subroutine elem_poisson_UW(Mdle,                                          &
                        NrTest,NrTrial,                                   &
                        NrdofHH,NrdofVV,NrdofH,NrdofQ,                    &
                        NrdofU,NrdofS,                                    &
                        NrdofVi_a,NrdofVi_b,                              &
                        MdVi_a,MdVi_b,MdU,MdS,                            &
                        BlocVi_a,AlocVi_aa,AlocVi_ab,AlocVi_aU,AlocVi_aS, &
                        BlocVi_b,AlocVi_ba,AlocVi_bb,AlocVi_bU,AlocVi_bS, & 
                        BlocU,AlocUVi_a,AlocUVi_b,AlocUU,AlocUS,          &
                        BlocS,AlocSVi_a,AlocSVi_b,AlocSU,AlocSS)
!
!..ALOC: holds local element stiffness matrices
!..BLOC: holds local element load vectors
   use control
   use data_structure3D
   use element_data
   use parametersDPG
   use mpi_param, only: RANK
!
   implicit none
!
!..declare input/output variables
   integer,                       intent(in)  :: Mdle
   integer,                       intent(in)  :: NrTest
   integer,                       intent(in)  :: NrTrial
   integer,                       intent(in)  :: NrdofHH
   integer,                       intent(in)  :: NrdofVV
   integer,                       intent(in)  :: NrdofH
   integer,                       intent(in)  :: NrdofQ
   integer,                       intent(in)  :: NrdofU
   integer,                       intent(in)  :: NrdofS
   integer,                       intent(in)  :: NrdofVi_a
   integer,                       intent(in)  :: NrdofVi_b
   integer,                       intent(in)  :: MdU
   integer,                       intent(in)  :: MdS
   integer,                       intent(in)  :: MdVi_a
   integer,                       intent(in)  :: MdVi_b


   real(8), dimension(MdU),       intent(out) :: BlocU
   real(8), dimension(MdU,MdU),   intent(out) :: AlocUU
   real(8), dimension(MdU,MdS),   intent(out) :: AlocUS
   real(8), dimension(MdU,MdVi_a),   intent(out) :: AlocUVi_a
   real(8), dimension(MdU,MdVi_b),   intent(out) :: AlocUVi_b
   
   real(8), dimension(MdS),       intent(out) :: BlocS
   real(8), dimension(MdS,MdU),   intent(out) :: AlocSU
   real(8), dimension(MdS,MdS),   intent(out) :: AlocSS
   real(8), dimension(MdS,MdVi_a),   intent(out) :: AlocSVi_a
   real(8), dimension(MdS,MdVi_b),   intent(out) :: AlocSVi_b

   real(8), dimension(MdVi_a),       intent(out) :: BlocVi_a
   real(8), dimension(MdVi_a,MdU),   intent(out) :: AlocVi_aU
   real(8), dimension(MdVi_a,MdS),   intent(out) :: AlocVi_aS
   real(8), dimension(MdVi_a,MdVi_a),   intent(out) :: AlocVi_aa
   real(8), dimension(MdVi_a,MdVi_b),   intent(out) :: AlocVi_ab

   real(8), dimension(MdVi_b),       intent(out) :: BlocVi_b
   real(8), dimension(MdVi_b,MdU),   intent(out) :: AlocVi_bU
   real(8), dimension(MdVi_b,MdS),   intent(out) :: AlocVi_bS
   real(8), dimension(MdVi_b,MdVi_a),   intent(out) :: AlocVi_ba
   real(8), dimension(MdVi_b,MdVi_b),   intent(out) :: AlocVi_bb
!
!-------------------------------------------------------------------------
!
!..aux variables
   real(8) :: rjac, bjac, fval, wa, weight
   integer :: iflag, nrv, nre, nrf, nint
   integer :: j1, j2, j3, j4, k1, k2, i, i1, i2, k, l,n1,n2,ivar
   integer :: nordP, nrdof, nsign, ifc, info
!
   character(len=4) :: etype,ftype
!
!..element order, orientation for edges and faces
   integer :: norder(19), norient_edge(12), norient_face(6)
!
!..element nodes order (trial) for interfaces
   integer :: norderi(19)
!
!..face order
   integer :: norderf(5)
!
!..geometry dof
   real(8) :: xnod(3,MAXbrickH)
!
!..geometry
   real(8) :: xi(3), x(3), dxdxi(3,3), dxidx(3,3), rn(3)
   real(8) :: dxidt(3,2), dxdt(3,2), t(2)
!
!..H1 shape functions
   real(8) :: shapH(MAXbrickH), gradH(3,MAXbrickH)
   ! real(8) :: shapHB(MAXbrickH), gradHB(3,MAXbrickH)
!
!..H(div) shape functions
   real(8) :: shapV(3,MAXbrickV), divV(MAXbrickV)

!.. L2 shape functions
   real(8) :: shapQ(MAXbrickQ)
!
!..Enriched H1 shape functions
   real(8) :: shapHH(MAXbrickHH), gradHH(3,MAXbrickHH)

!..Enriched H(div) shape functions
   real(8) :: shapVV(3,MAXbrickVV), divVV(MAXbrickVV)

!..load vector for the enriched space
   real(8) :: bload_V(NrTest)
   real(8) :: dummy_check
!
!..Gram matrix in packed format
   ! real(8), allocatable :: gramPV(:)
   ! real(8), allocatable :: gramPT(:)
   real(8), allocatable :: gramP(:)
!
!..stiffness matrices for the enriched test space
   real(8), allocatable :: stiff_UV(:,:),stiff_SV(:,:),stiff_Vi_aV(:,:),stiff_Vi_bV(:,:)
   real(8), allocatable :: stiff_UT(:,:),stiff_ST(:,:),stiff_Vi_aT(:,:),stiff_Vi_bT(:,:)
   real(8), allocatable :: stiff_ALL(:,:),raloc(:,:)
!
!..3D quadrature data
   real(8) :: xiloc(3,MAXNINT3ADD), waloc(MAXNINT3ADD)
!
!..2D quadrature data
   real(8) :: tloc(2,MAXNINT2ADD), wtloc(MAXNINT2ADD)
!
!..workspace for trial and test variables
   real(8) :: u, v, q, divtau_a, divtau_b, tn, sn, lambda, aux
   real(8) :: sig(3), dv(3), dq(3), tau_a(3), tau_b(3), s(3)
!
   integer :: nk
   nk(k1,k2) = (k2-1)*k2/2 + k1
!
!---------------------------------------------------------------------
!--------- INITIALIZE THE ELEMENT ORDER, ORIENTATION ETC. ------------
!---------------------------------------------------------------------
!
!
!..allocate auxiliary matrices
   allocate(gramP((NrTest)*(NrTest+1)/2))
   allocate(stiff_UV(NrdofHH,NrdofU))
   allocate(stiff_SV(NrdofHH,NrdofS))
   allocate(stiff_Vi_aV(NrdofHH,NrdofVi_a))
   allocate(stiff_Vi_bV(NrdofHH,NrdofVi_b))
   allocate(stiff_UT(NrdofVV,NrdofU))
   allocate(stiff_ST(NrdofVV,NrdofS))
   allocate(stiff_Vi_aT(NrdofVV,NrdofVi_a))
   allocate(stiff_Vi_bT(NrdofVV,NrdofVi_b))
  
!
!..element type
   etype = NODES(Mdle)%type
   nrv = nvert(etype)
   nre = nedge(etype)
   nrf = nface(etype)
!
!..determine order of approximation
   call find_order(Mdle, norder)
   norderi(1:nre+nrf) = norder(1:nre+nrf)
!
!..set the enriched order of approximation
   select case(etype)
      case('mdlb')
         nordP = NODES(Mdle)%order+NORD_ADD*111
         norderi(nre+nrf+1) = 111
      case('mdlp')
         nordP = NODES(Mdle)%order+NORD_ADD*11
         norderi(nre+nrf+1) = 11
      case('mdln','mdld')
         nordP = NODES(Mdle)%order+NORD_ADD
         norderi(nre+nrf+1) = 1
      case default
         write(*,*) 'elem_poisson: invalid etype param. stop.'
         stop
   end select
!
!..determine edge and face orientations
   call find_orient(Mdle, norient_edge,norient_face)
!
!..determine nodes coordinates
   call nodcor(Mdle, xnod)
!
!..clear space for stiffness matrix and load vector:
   BlocU = ZERO; AlocUU = ZERO; AlocUS = ZERO; AlocUVi_a = ZERO; AlocUVi_b = ZERO
   BlocS = ZERO; AlocSU = ZERO; AlocSS = ZERO; AlocSVi_a = ZERO; AlocSVi_b = ZERO
   BlocVi_a = ZERO; AlocVi_aU = ZERO; AlocVi_aS = ZERO; AlocVi_aa = ZERO; AlocVi_ab = ZERO
   BlocVi_b = ZERO; AlocVi_bU = ZERO; AlocVi_bS = ZERO; AlocVi_ba = ZERO; AlocVi_bb = ZERO
!
!..clear space for auxiliary matrices
   bload_V   = ZERO
   gramP     = ZERO
   stiff_UV  = ZERO
   stiff_SV  = ZERO
   stiff_Vi_aV = ZERO
   stiff_Vi_bV = ZERO
   stiff_UT = ZERO
   stiff_ST = ZERO
   stiff_Vi_aT = ZERO
   stiff_Vi_bT = ZERO

!
!----------------------------------------------------------------------
!     E L E M E N T    I N T E G R A L S                              |
!----------------------------------------------------------------------
!
!..use the enriched order to set the quadrature
   INTEGRATION = NORD_ADD
   call set_3D_int_DPG(etype,norder,norient_face, nint,xiloc,waloc)
   INTEGRATION = 0
!
!..loop over integration points
   do l=1,nint
!
!  ...coordinates and weight of this integration point
      xi(1:3)=xiloc(1:3,l)
      wa=waloc(l)
!
!  ...H1 shape functions (for geometry)
      call shape3DH(etype,xi,norder,norient_edge,norient_face, nrdof,shapH,gradH)

!  ...L2 shape function calls
      call shape3DQ(etype,xi,norder, nrdof,shapQ)
!
!  ...discontinuous H1 shape functions
      call shape3HH(etype,xi,nordP, nrdof,shapHH,gradHH)
!  ...discontinuous H(div) shape functions
      call shape3VV(etype,xi,nordP, nrdof,shapVV,divVV)

!  ...geometry map
      call geom3D(Mdle,xi,xnod,shapH,gradH,NrdofH, x,dxdxi,dxidx,rjac,iflag)
!
!  ...integration weight
      weight = rjac*wa
!
!  ...get the RHS
      call getf(Mdle,x, fval)
!
!  ...1st loop through enriched H1 test functions
      do k1=1,NrdofHH
!     ...Piola transformation
         v = shapHH(k1)
         dv(1:3) = gradHH(1,k1)*dxidx(1,1:3) &
                 + gradHH(2,k1)*dxidx(2,1:3) &
                 + gradHH(3,k1)*dxidx(3,1:3)
!
!     ...EXTENDED LOAD VECTOR: (f,v)
         bload_V(k1) = bload_V(k1) + fval*v*weight
!
!     ...loop through L2 trial functions
         do k2=1,NrdofQ

            do ivar = 1,3
               n1 = k1
               ! n2 = (ivar - 1)*NrdofQ + k2
               n2 = (k2 - 1)*3 + ivar
               sig = ZERO
!        ...Piola Transform for L2 fields i.e for the ivar th component of sigma
               sig(ivar) = shapQ(k2)/rjac 
               stiff_SV(n1,n2) = stiff_SV(n1,n2) + weight * (sig(1) * dv(1) &
                                                            +sig(2) * dv(2) + sig(3)*dv(3))
            enddo
!
!     ...end of loop through L2 trial functions
         enddo
!
!     ...2nd loop through enriched H1 test functions for Gram matrix
         do k2=k1,NrdofHH
!        ...Piola transformation
            q = shapHH(k2)
            dq(1:3) = gradHH(1,k2)*dxidx(1,1:3) &
                    + gradHH(2,k2)*dxidx(2,1:3) &
                    + gradHH(3,k2)*dxidx(3,1:3)
!
!        ...determine index in triangular packed format            
            k = nk(k1,k2)
            aux = q*v + (dq(1)*dv(1) + dq(2)*dv(2) + dq(3)*dv(3))
            gramP(k) = gramP(k) + aux*weight
!           
!     ...enddo 2nd loop through enriched H1 test functions
         enddo

      ! cross terms for  graph norm
         do k2 = 1,NrdofVV
            divtau_a = divVV(k2)/rjac
            tau_a(1:3) =     dxdxi(1:3,1) * shapVV(1,k2)      &
                           + dxdxi(1:3,2) * shapVV(2,k2)       &
                           + dxdxi(1:3,3) * shapVV(3,k2)
            
            tau_a(1:3) = tau_a(1:3)/rjac
            k = nk(k1,NrdofHH+k2)

            aux = dv(1) * tau_a(1) + dv(2) * tau_a(2) + dv(3) * tau_a(3)
            gramP(k) = gramP(k) + aux * weight
         enddo

!  ...end of 1st loop through enriched H1 test functions
      enddo

!   ... loop over discontinuous H(div) test functions
      do k1=1,NrdofVV
         ! Piola Transform of the divergence
         divtau_a = divVV(k1)/rjac
         tau_a(1:3) =     dxdxi(1:3,1) * shapVV(1,k1)      &
                        + dxdxi(1:3,2) * shapVV(2,k1)       &
                        + dxdxi(1:3,3) * shapVV(3,k1)

         tau_a(1:3) = tau_a(1:3)/rjac   
         
         do k2=1,NrdofQ
            ! Piola transform of L2 variable
            u = shapQ(k2)/rjac
            stiff_UT(k1,k2) = stiff_UT(k1,k2) + weight*(u * divtau_a)

            do ivar = 1,3
               sig = ZERO
               sig(ivar) = shapQ(k2)/rjac  !Piola Transform for the ivar th comp
               n1 = k1
               ! n2  = (ivar-1) * NrdofQ + k2
               n2 = (k2 - 1)*3 + ivar
               stiff_ST(n1,n2) = stiff_ST(n1,n2) + weight*(tau_a(1)*sig(1) + tau_a(2) * sig(2) &
                                                          + tau_a(3)*sig(3))

            enddo
         enddo

       ! ...Gram matrix contribution for H(div) inner product
         do k2=k1,NrdofVV
            !Piola transform for H(div) test function and its divergence. 
            divtau_b = divVV(k2)/rjac
            tau_b(1:3) =     dxdxi(1:3,1) * shapVV(1,k2)      &
                           + dxdxi(1:3,2) * shapVV(2,k2)       &
                           + dxdxi(1:3,3) * shapVV(3,k2)
   
            tau_b(1:3) = tau_b(1:3)/rjac

            k = nk(k1 + NrdofHH,k2 + NrdofHH)
            aux = divtau_a * divtau_b + 2.d0 * (tau_a(1)*tau_b(1) + tau_a(2)*tau_b(2) + tau_a(3)*tau_b(3))
            gramP(k) = gramP(k) +  weight * aux

         enddo   

      enddo

!..end of loop through integration points
   enddo
  


   open(1, file = 'data2.dat', status='replace')  
   do k1=1,NrTest  
      do k2=k1,NrTest
       write(1,*) k1,",",k2,",",gramP(nk(k1,k2))  
      enddo     
   enddo  
   
   ! write(*,*) NrTest,",",NrdofHH,",",NrdofVV
   close(1) 


!
!---------------------------------------------------------------------
!    B O U N D A R Y    I N T E G R A L S                            |
!---------------------------------------------------------------------
! 
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
      INTEGRATION = NORD_ADD
      call set_2D_int_DPG(ftype,norderf,norient_face(ifc), nint,tloc,wtloc)
      INTEGRATION = 0
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
!     ...determine discontinuous H1 shape functions
         call shape3HH(etype,xi,nordP, nrdof,shapHH,gradHH)
!     ...discontinuous H(div) shape functions
         call shape3VV(etype,xi,nordP, nrdof,shapVV,divVV)
!
!     ...determine element H1 shape functions (for geometry)
         call shape3DH(etype,xi,norder,norient_edge,norient_face, &
                       nrdof,shapH,gradH)
!
!     ...determine element H(div) shape functions (for fluxes)
!     ...for interfaces only (no bubbles)
         call shape3DV(etype,xi,norderi,norient_face, &
                       nrdof,shapV,divV)

!     ...geometry
         call bgeom3D(Mdle,xi,xnod,shapH,gradH,NrdofH,dxidt,nsign, &
                  x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
         weight = bjac*wtloc(l)
!
!     ...loop through enriched H1 test functions
         do k1=1,NrdofHH
            v = shapHH(k1)
!
!        ...loop through H(div) trial functions
            do k2=1,NrdofVi_b
!           ...Piola transformation
               s(1:3) = dxdxi(1:3,1)*shapV(1,k2) &
                        + dxdxi(1:3,2)*shapV(2,k2) &
                        + dxdxi(1:3,3)*shapV(3,k2)
               s(1:3) = s(1:3)/rjac
!           ...normal component
               sn = s(1)*rn(1)+s(2)*rn(2)+s(3)*rn(3)
!
!           ...EXTENDED HV STIFFNESS MATRIX: -<dot(σ,n), v>_{Γ_h}
               stiff_Vi_bV(k1,k2) = stiff_Vi_bV(k1,k2) - sn*v*weight
!        ...end loop through H(div) trial functions
            enddo
!     ...end loop through enriched H1 test functions
         enddo

!        ...loop through enriched H(div) test functions
         do k1=1,NrdofVV
            tau_a(1:3) =     dxdxi(1:3,1) * shapVV(1,k1)      &
                           + dxdxi(1:3,2) * shapVV(2,k1)       &
                           + dxdxi(1:3,3) * shapVV(3,k1)

            tau_a(1:3) = tau_a(1:3)/rjac
            tn = tau_a(1)*rn(1) + tau_a(2)*rn(2) + tau_a(3)*rn(3)
            do k2=1,NrdofVi_a
               lambda = shapH(k2)

               stiff_Vi_aT(k1,k2) = stiff_Vi_aT(k1,k2) - weight * tn * lambda

            enddo

         enddo   

!     ...end loop through integration points
      enddo
!..end loop through element faces
   enddo
!
!---------------------------------------------------------------------
!  Construction of DPG system
!---------------------------------------------------------------------
!
   allocate(stiff_ALL(NrTest,NrTrial+1))
   stiff_ALL = ZERO

!  Total test/trial DOFs of the element
   i1 = NrdofHH; i2 = NrdofVV; j1 = NrdofVi_a; j2 = NrdofVi_b; j3 = NrdofU; j4 = NrdofS
   stiff_ALL(1:i1,1:j1) = stiff_Vi_aV
   stiff_ALL(1:i1,j1+1:j1+j2) = stiff_Vi_bV
   stiff_ALL(1:i1,j1+j2+1:j1+j2+j3) = stiff_UV
   stiff_ALL(1:i1,j1+j2+j3+1:j1+j2+j3+j4) = stiff_SV

   stiff_ALL(i1+1:i1+i2,1:j1) = stiff_Vi_aT
   stiff_ALL(i1+1:i1+i2,j1+1:j1+j2) = stiff_Vi_bT
   stiff_ALL(i1+1:i1+i2,j1+j2+1:j1+j2+j3) = stiff_UT
   stiff_ALL(i1+1:i1+i2,j1+j2+j3+1:j1+j2+j3+j4) = stiff_ST

   stiff_ALL(1:i1+i2,j1+j2+j3+j4+1) = bload_V
   


   deallocate(stiff_UV,stiff_SV,stiff_Vi_aV,stiff_Vi_bV)
   deallocate(stiff_UT,stiff_ST,stiff_Vi_aT,stiff_Vi_bT)

!..A. Compute Cholesky factorization of Gram Matrix, G=U^T U (=LL^T)

   call DPPTRF('U',NrTest,gramP,info)
   if (info.ne.0) then
      write(*,*) 'elem_heat: DPPTRF: Mdle,info = ',Mdle,info,'. stop.'
      stop
   endif

!..B. Solve triangular system to obtain B~, (LX=) U^T X = [B|l]
   call DTPTRS('U','T','N',NrTest,NrTrial+1,gramP,stiff_ALL,NrTest,info)
   if (info.ne.0) then
      write(*,*) 'elem_heat: DTPTRS: Mdle,info = ',Mdle,info,'. stop.'
      stop
   endif

   allocate(raloc(NrTrial+1,NrTrial+1)); raloc = ZERO

!..C. Matrix multiply: B^T G^-1 B (=B~^T B~)
   call DSYRK('U','T',NrTrial+1,NrTest,ZONE,stiff_ALL,NrTest,ZERO,raloc,NrTrial+1)

!
   deallocate(stiff_ALL)

!..D. Fill lower triangular part of Hermitian matrix using the upper triangular matrix
   do i=1,NrTrial
      raloc(i+1:NrTrial+1,i) = raloc(i,i+1:NrTrial+1)
   enddo

!  Fill the Aloc and Bloc matrices
   BlocVi_a(1:j1) = raloc(1:j1,j1+j2+j3+j4+1)
   BlocVi_b(1:j2) = raloc(j1+1:j1+j2,j1+j2+j3+j4+1)
   BlocU(1:j3) = raloc(j1+j2+1:j1+j2+j3,j1+j2+j3+j4+1)
   BlocS(1:j4) = raloc(j1+j2+j3+1:j1+j2+j3+j4,j1+j2+j3+j4+1)


   AlocVi_aa(1:j1,1:j1) = raloc(1:j1,1:j1)
   AlocVi_ab(1:j1,1:j2) = raloc(1:j1,j1+1:j1+j2)
   AlocVi_aU(1:j1,1:j3) = raloc(1:j1,j1+j2+1:j1+j2+j3)
   AlocVi_aS(1:j1,1:j4) = raloc(1:j1,j1+j2+j3+1:j1+j2+j3+j4)

   
   AlocVi_ba(1:j2,1:j1) = raloc(j1+1:j1+j2,1:j1)
   AlocVi_bb(1:j2,1:j2) = raloc(j1+1:j1+j2,j1+1:j1+j2)
   AlocVi_bU(1:j2,1:j3) = raloc(j1+1:j1+j2,j1+j2+1:j1+j2+j3)
   AlocVi_bS(1:j2,1:j4) = raloc(j1+1:j1+j2,j1+j2+j3+1:j1+j2+j3+j4)


   AlocUVi_a(1:j3,1:j1) = raloc(j1+j2+1:j1+j2+j3,1:j1)
   AlocUVi_b(1:j3,1:j2) = raloc(j1+j2+1:j1+j2+j3,j1+1:j1+j2)
   AlocUU(1:j3,1:j3) =    raloc(j1+j2+1:j1+j2+j3,j1+j2+1:j1+j2+j3)
   AlocUS(1:j3,1:j4) =    raloc(j1+j2+1:j1+j2+j3,j1+j2+j3+1:j1+j2+j3+j4)


   AlocSVi_a(1:j4,1:j1) = raloc(j1+j2+j3+1:j1+j2+j3+j4,1:j1)
   AlocSVi_b(1:j4,1:j2) = raloc(j1+j2+j3+1:j1+j2+j3+j4,j1+1:j1+j2)
   AlocSU(1:j4,1:j3) =    raloc(j1+j2+j3+1:j1+j2+j3+j4,j1+j2+1:j1+j2+j3)
   AlocSS(1:j4,1:j4) =    raloc(j1+j2+j3+1:j1+j2+j3+j4,j1+j2+j3+1:j1+j2+j3+j4)






   deallocate(raloc)

!
end subroutine elem_poisson_UW

!----------------------------------------------------------------------
!  routine: compute_enriched_order
!----------------------------------------------------------------------
!  purpose: - compute enriched order vector based on input Nord
!             (enriched order based on mdle order and Nord only)
!----------------------------------------------------------------------
subroutine compute_enriched_order(EType,Nord, Norder)
!
   use parameters, only : MODORDER
!
   implicit none
!
   character(len=4), intent(in)  :: Etype
   integer         , intent(in)  :: Nord
   integer         , intent(out) :: Norder(19)
!
   integer :: temp(2)
   integer :: nordF(3),nordB(3)
!
!----------------------------------------------------------------------
!
!..see implementation of BrokenExactSequence in shape functions
   select case(Etype)
      case('mdlb')
         call decod(Nord,MODORDER,2, temp) !xy face, z edge
         nordF(1) = temp(1); nordB(3) = temp(2)
         call decod(nordF(1),MODORDER,2, nordB(1:2)) !x edge, y edge
         call encod((/nordB(1),nordB(3)/),MODORDER,2, nordF(2)) !xz face
         call encod(nordB(2:3),MODORDER,2, nordF(3)) !yz face
!     ...edges
         Norder(1:4)   = (/nordB(1),nordB(2),nordB(1),nordB(2)/) !x,y,x,y edges
         Norder(5:8)   = (/nordB(1),nordB(2),nordB(1),nordB(2)/) !x,y,x,y edges
         Norder(9:12)  = nordB(3) !z edges
!     ...faces
         Norder(13:14) = nordF(1) !xy faces
         Norder(15:18) = (/nordF(2),nordF(3),nordF(2),nordF(3)/) !xz,yz,xz,yz faces
!     ...element interior
         Norder(19)    = Nord
      case default
         write(*,*) 'compute_enriched_order: only implemented for hexa. stop.'
         stop
   end select
!
end subroutine compute_enriched_order
!

