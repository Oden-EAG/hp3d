!-------------------------------------------------------------------------
!
!     routine name      - elem
!
!-------------------------------------------------------------------------
!
!     latest revision:  - May 2023
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
   integer :: etype
!
!-------------------------------------------------------------------------
!
!..activate four physics variables (H1 trace + H(div) trace +L2 (field variable) 
!   + L2 (approximate gradient) ) for assembly
   Itest(1:NR_PHYSA) = 1; Itrial(1:NR_PHYSA) = 1
!
   norder (1:19) = 0
   norderP(1:19) = 0
!
   etype = NODES(Mdle)%ntype
   nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
!..determine order of approximation
   call find_order(Mdle, norder)
!
!..set the enriched order of approximation
   select case(etype)
      case(MDLB)
         nordP = NODES(Mdle)%order+NORD_ADD*111
      case(MDLP)
         nordP = NODES(Mdle)%order+NORD_ADD*11
      case(MDLN,MDLD)
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
!
!..call element integration routine
   call elem_poisson_UW(Mdle,nrTest,nrTrial,                                                 &
            nrdofHH,nrdofVV,nrdofH,nrdofQ,nrdofQ,3*nrdofQ,nrdofVi_a,nrdofVi_b,               &
            BLOC(1)%nrow,BLOC(2)%nrow,BLOC(3)%nrow,BLOC(4)%nrow,                             &
            BLOC(1)%array,ALOC(1,1)%array,ALOC(1,2)%array,ALOC(1,3)%array,ALOC(1,4)%array,   &
            BLOC(2)%array,ALOC(2,1)%array,ALOC(2,2)%array,ALOC(2,3)%array,ALOC(2,4)%array,   &
            BLOC(3)%array,ALOC(3,1)%array,ALOC(3,2)%array,ALOC(3,3)%array,ALOC(3,4)%array,   &
            BLOC(4)%array,ALOC(4,1)%array,ALOC(4,2)%array,ALOC(4,3)%array,ALOC(4,4)%array)
!     
end subroutine elem
!
!-------------------------------------------------------------------------
!
!     routine name      - elem_poisson_UW
!
!-------------------------------------------------------------------------
!
!     latest revision:  - May 2023
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
!              NrdofH   - number of H1 trial dof
!              NrdofQ   - number of L2 dofs which will be used for u and sigma_1,2,3 
!              NrdofU   - number of L2 trial dof for u
!              NrdofS   - number of L2 trial dof for sigma  = grad u (Dimension X NrdofQ)
!              NrdofVi_a- number of H1 trial interface dof
!              NrdofVi_b- number of H(div) trial interface dof
!              MdQ1     - num rows of AlocQ1Q1,AlocQ1Q2,AlocQ1H,AlocQ1V
!              MdQ2     - num rows of AlocQ2Q1,AlocQ2Q2,AlocQ2H,AlocQ2V
!              MdH      - num rows of AlocHQ1,AlocHQ2,AlocHH,AlocHV
!              MdV      - num rows of AlocVQ1,AlocVSQ2,AlocVH,AlocVV
!
!        out:
!              BlocQ1
!              BlocQ2
!              BlocH    - load vectors
!              BlocV
!
!              AlocQ1Q1
!              AlocQ1Q2
!              AlocQ1H
!              AlocQ1V
!              AlocQ2Q1
!              AlocQ2Q2
!              AlocQ2H
!              AlocQ2V
!              AlocHQ1
!              AlocHQ2  - stiffness matrices
!              AlocHH
!              AlocHV
!              AlocVQ1
!              AlocVQ2
!              AlocVH
!              AlocVV
!
!-------------------------------------------------------------------------
!
subroutine elem_poisson_UW(Mdle,                                        &
                           NrTest,NrTrial,                              &
                           NrdofHH,NrdofVV,NrdofH,NrdofQ,               &
                           NrdofU,NrdofS,                               &
                           NrdofVi_a,NrdofVi_b,                         &
                           MdH,MdV,MdQ1,MdQ2,                           &
                           BlocH,AlocHH,AlocHV,AlocHQ1,AlocHQ2,         &
                           BlocV,AlocVH,AlocVV,AlocVQ1,AlocVQ2,         &
                           BlocQ1,AlocQ1H,AlocQ1V,AlocQ1Q1,AlocQ1Q2,    &
                           BlocQ2,AlocQ2H,AlocQ2V,AlocQ2Q1,AlocQ2Q2)
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
   integer, intent(in) :: Mdle
   integer, intent(in) :: NrTest
   integer, intent(in) :: NrTrial
   integer, intent(in) :: NrdofHH
   integer, intent(in) :: NrdofVV
   integer, intent(in) :: NrdofH
   integer, intent(in) :: NrdofQ
   integer, intent(in) :: NrdofU
   integer, intent(in) :: NrdofS
   integer, intent(in) :: NrdofVi_a
   integer, intent(in) :: NrdofVi_b
   integer, intent(in) :: MdQ1
   integer, intent(in) :: MdQ2
   integer, intent(in) :: MdH
   integer, intent(in) :: MdV
!
   real(8), dimension(MdQ1),       intent(out) :: BlocQ1
   real(8), dimension(MdQ1,MdQ1),  intent(out) :: AlocQ1Q1
   real(8), dimension(MdQ1,MdQ2),  intent(out) :: AlocQ1Q2
   real(8), dimension(MdQ1,MdH),   intent(out) :: AlocQ1H
   real(8), dimension(MdQ1,MdV),   intent(out) :: AlocQ1V
!
   real(8), dimension(MdQ2),       intent(out) :: BlocQ2
   real(8), dimension(MdQ2,MdQ1),  intent(out) :: AlocQ2Q1
   real(8), dimension(MdQ2,MdQ2),  intent(out) :: AlocQ2Q2
   real(8), dimension(MdQ2,MdH),   intent(out) :: AlocQ2H
   real(8), dimension(MdQ2,MdV),   intent(out) :: AlocQ2V
!
   real(8), dimension(MdH),        intent(out) :: BlocH
   real(8), dimension(MdH,MdQ1),   intent(out) :: AlocHQ1
   real(8), dimension(MdH,MdQ2),   intent(out) :: AlocHQ2
   real(8), dimension(MdH,MdH),    intent(out) :: AlocHH
   real(8), dimension(MdH,MdV),    intent(out) :: AlocHV
!
   real(8), dimension(MdV),        intent(out) :: BlocV
   real(8), dimension(MdV,MdQ1),   intent(out) :: AlocVQ1
   real(8), dimension(MdV,MdQ2),   intent(out) :: AlocVQ2
   real(8), dimension(MdV,MdH),    intent(out) :: AlocVH
   real(8), dimension(MdV,MdV),    intent(out) :: AlocVV
!
!-------------------------------------------------------------------------
!
!..aux variables
   real(8) :: rjac, bjac, fval, wa, weight
   integer :: iflag, nrv, nre, nrf, nint
   integer :: j1, j2, j3, j4, k1, k2, i, i1, i2, k, l,n1,n2,ivar
   integer :: nordP, nrdof, nsign, ifc, info
!
   integer :: etype,ftype
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
!
!..load vector for the enriched space
   real(8) :: bload_H(NrTest)
!
!..Gram matrix in packed format
   real(8), allocatable :: gramP(:)
!
!..stiffness matrices for the enriched test space
   real(8), allocatable :: stiff_HQ1(:,:),stiff_HQ2(:,:),stiff_HH(:,:),stiff_HV(:,:)
   real(8), allocatable :: stiff_VQ1(:,:),stiff_VQ2(:,:),stiff_VH(:,:),stiff_VV(:,:)
   real(8), allocatable :: stiff_ALL(:,:),raloc(:,:)
!
!..3D quadrature data
   real(8) :: xiloc(3,MAXNINT3ADD), waloc(MAXNINT3ADD)
!
!..2D quadrature data
   real(8) :: tloc(2,MAXNINT2ADD), wtloc(MAXNINT2ADD)
!
!..workspace for trial and test variables
   real(8) :: u, v, q, divtau_a, divtau_b, tn, sn, u_hat, aux
   real(8) :: sig(3), dv(3), dq(3), tau_a(3), tau_b(3), s(3)
!
   integer, external :: ij_upper_to_packed
!
!---------------------------------------------------------------------
!--------- INITIALIZE THE ELEMENT ORDER, ORIENTATION ETC. ------------
!---------------------------------------------------------------------
!
!
!..allocate auxiliary matrices
   allocate(gramP((NrTest)*(NrTest+1)/2))
   allocate(stiff_HQ1(NrdofHH,NrdofU))
   allocate(stiff_HQ2(NrdofHH,NrdofS))
   allocate(stiff_HH(NrdofHH,NrdofVi_a))
   allocate(stiff_HV(NrdofHH,NrdofVi_b))
   allocate(stiff_VQ1(NrdofVV,NrdofU))
   allocate(stiff_VQ2(NrdofVV,NrdofS))
   allocate(stiff_VH(NrdofVV,NrdofVi_a))
   allocate(stiff_VV(NrdofVV,NrdofVi_b))
!
!..element type
   etype = NODES(Mdle)%ntype
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
      case(MDLB)
         nordP = NODES(Mdle)%order+NORD_ADD*111
         norderi(nre+nrf+1) = 111
      case(MDLP)
         nordP = NODES(Mdle)%order+NORD_ADD*11
         norderi(nre+nrf+1) = 11
      case(MDLN,MDLD)
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
   BlocQ1 = ZERO; AlocQ1Q1 = ZERO; AlocQ1Q2 = ZERO; AlocQ1H = ZERO; AlocQ1V = ZERO
   BlocQ2 = ZERO; AlocQ2Q1 = ZERO; AlocQ2Q2 = ZERO; AlocQ2H = ZERO; AlocQ2V = ZERO
   BlocH  = ZERO; AlocHQ1  = ZERO; AlocHQ2  = ZERO; AlocHH  = ZERO; AlocHV  = ZERO
   BlocV  = ZERO; AlocVQ1  = ZERO; AlocVQ2  = ZERO; AlocVH  = ZERO; AlocVV  = ZERO
!
!..clear space for auxiliary matrices
   bload_H   = ZERO
   gramP     = ZERO
   stiff_HQ1 = ZERO
   stiff_HQ2 = ZERO
   stiff_HH  = ZERO
   stiff_HV  = ZERO
   stiff_VQ1 = ZERO
   stiff_VQ2 = ZERO
   stiff_VH  = ZERO
   stiff_VV  = ZERO
!
!----------------------------------------------------------------------
!     E L E M E N T    I N T E G R A L S                              |
!----------------------------------------------------------------------
!
!..use the enriched order to set the quadrature
   INTEGRATION = NORD_ADD
   call set_3D_int_DPG(etype,norder,norient_face, nint,xiloc,waloc)
!
!..loop over integration points
   do l=1,nint
!
!..coordinates and weight of this integration point
      xi(1:3)=xiloc(1:3,l)
      wa=waloc(l)
!
!..H1 shape functions (for geometry)
      call shape3DH(etype,xi,norder,norient_edge,norient_face, nrdof,shapH,gradH)

!..L2 shape functions
      call shape3DQ(etype,xi,norder, nrdof,shapQ)
!
!..discontinuous H1 shape functions
      call shape3HH(etype,xi,nordP, nrdof,shapHH,gradHH)
!..discontinuous H(div) shape functions
      call shape3VV(etype,xi,nordP, nrdof,shapVV,divVV)

!..geometry map
      call geom3D(Mdle,xi,xnod,shapH,gradH,NrdofH, x,dxdxi,dxidx,rjac,iflag)
!
!..integration weight
      weight = rjac*wa
!
!..get the RHS
      call getf(Mdle,x, fval)
!
!..1st loop through enriched H1 test functions
      do k1=1,NrdofHH
!..Piola transformation
         v = shapHH(k1)
         dv(1:3) = gradHH(1,k1)*dxidx(1,1:3) &
                 + gradHH(2,k1)*dxidx(2,1:3) &
                 + gradHH(3,k1)*dxidx(3,1:3)
!
!..accumulate load vector: (f,v)
         bload_H(k1) = bload_H(k1) + fval*v*weight
!
!..loop through L2 trial functions
         do k2=1,NrdofQ
!
            do ivar = 1,3
               n1 = k1
               n2 = (k2 - 1)*3 + ivar
               sig = ZERO
!..Piola Transform for L2 fields i.e for the ivar-th component of sigma
               sig(ivar) = shapQ(k2)/rjac
               stiff_HQ2(n1,n2) = stiff_HQ2(n1,n2) + weight * &
                                  (sig(1)*dv(1) + sig(2)*dv(2) + sig(3)*dv(3))
            enddo
!
!..end of loop through L2 trial functions
         enddo
!
!..2nd loop through enriched H1 test functions for Gram matrix
         do k2=k1,NrdofHH
!..Piola transformation
            q = shapHH(k2)
            dq(1:3) = gradHH(1,k2)*dxidx(1,1:3) &
                    + gradHH(2,k2)*dxidx(2,1:3) &
                    + gradHH(3,k2)*dxidx(3,1:3)
!
!..determine index in triangular packed format            
            k = ij_upper_to_packed(k1,k2)
!..accumlate Gram components of Gram matrix correspoding to v
            aux = q*v + (dq(1)*dv(1) + dq(2)*dv(2) + dq(3)*dv(3))
            gramP(k) = gramP(k) + aux*weight
!           
!..end of 2nd loop through enriched H1 test functions
         enddo
!
!..cross terms for  graph norm
         do k2 = 1,NrdofVV
            tau_a(1:3) = dxdxi(1:3,1) * shapVV(1,k2) &
                       + dxdxi(1:3,2) * shapVV(2,k2) &
                       + dxdxi(1:3,3) * shapVV(3,k2)
            tau_a(1:3) = tau_a(1:3)/rjac
!
            k = ij_upper_to_packed(k1,NrdofHH+k2)
!
            aux = dv(1)*tau_a(1) + dv(2)*tau_a(2) + dv(3)*tau_a(3)
            gramP(k) = gramP(k) + aux * weight
         enddo
!
!..end of 1st loop through enriched H1 test functions
      enddo
!..loop over discontinuous H(div) test functions
      do k1=1,NrdofVV
!..Piola Transform of the H(div) test functions.
         divtau_a = divVV(k1)/rjac
         tau_a(1:3) = dxdxi(1:3,1) * shapVV(1,k1) &
                    + dxdxi(1:3,2) * shapVV(2,k1) &
                    + dxdxi(1:3,3) * shapVV(3,k1)
         tau_a(1:3) = tau_a(1:3)/rjac   
!
!..loop over L2 trial functions
         do k2=1,NrdofQ
!..Piola transform of L2 variable
            u = shapQ(k2)/rjac
            stiff_VQ1(k1,k2) = stiff_VQ1(k1,k2) + weight*(u * divtau_a)
!..loop over L2 components
            do ivar = 1,3
               sig = ZERO
               sig(ivar) = shapQ(k2)/rjac ! Piola Transform for the ivar-th comp
               n1 = k1
               n2 = (k2 - 1)*3 + ivar
               stiff_VQ2(n1,n2) = stiff_VQ2(n1,n2) + weight * &
                                 (tau_a(1)*sig(1) + tau_a(2)*sig(2) + tau_a(3)*sig(3))
!
!..end of loop over L2 components
            enddo
!..end of loop over L2 trial functions
         enddo
!
!..Gram matrix contribution for H(div) inner product
         do k2=k1,NrdofVV
!..Piola transform for H(div) test function and its divergence. 
            divtau_b = divVV(k2)/rjac
            tau_b(1:3) = dxdxi(1:3,1) * shapVV(1,k2) &
                       + dxdxi(1:3,2) * shapVV(2,k2) &
                       + dxdxi(1:3,3) * shapVV(3,k2)
            tau_b(1:3) = tau_b(1:3)/rjac
!
            k = ij_upper_to_packed(k1 + NrdofHH,k2 + NrdofHH)
            aux = divtau_a * divtau_b + 2.d0 * &
                  (tau_a(1)*tau_b(1) + tau_a(2)*tau_b(2) + tau_a(3)*tau_b(3))
            gramP(k) = gramP(k) +  weight * aux
!
         enddo   
!..end of loop over H(div) test functions
      enddo
!..end of loop through integration points
   enddo
!
!
!---------------------------------------------------------------------
!    B O U N D A R Y    I N T E G R A L S                            |
!---------------------------------------------------------------------
! 
!
!..loop through element faces
   do ifc=1,nrf
!
!..sign factor to determine the outward normal unit vector
      nsign = nsign_param(etype,ifc)
!
!..face type
      ftype = face_type(etype,ifc)
!
!..face order of approximation
      call face_order(etype,ifc,norder, norderf)
!
!..set 2D quadrature
      INTEGRATION = NORD_ADD
      call set_2D_int_DPG(ftype,norderf,norient_face(ifc), nint,tloc,wtloc)
      INTEGRATION = 0
!
!..loop through integration points
      do l=1,nint
!
!..face coordinates
         t(1:2) = tloc(1:2,l)
!
!..face parametrization
         call face_param(etype,ifc,t, xi,dxidt)
!
!..determine discontinuous H1 shape functions
         call shape3HH(etype,xi,nordP, nrdof,shapHH,gradHH)
!..discontinuous H(div) shape functions
         call shape3VV(etype,xi,nordP, nrdof,shapVV,divVV)
!
!..determine element H1 shape functions (for geometry)
         call shape3DH(etype,xi,norder,norient_edge,norient_face, &
                       nrdof,shapH,gradH)
!
!..determine H(div) shape functions (for fluxes)
         call shape3DV(etype,xi,norderi,norient_face, &
                       nrdof,shapV,divV)
!
!..geometry map
         call bgeom3D(Mdle,xi,xnod,shapH,gradH,NrdofH,dxidt,nsign, &
                      x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
!..integration weight
         weight = bjac*wtloc(l)
!
!..loop through enriched H1 test functions
         do k1=1,NrdofHH
            v = shapHH(k1)
!
!..loop through H(div) trial functions
            do k2=1,NrdofVi_b
!..Piola transformation
               s(1:3) = dxdxi(1:3,1)*shapV(1,k2) &
                      + dxdxi(1:3,2)*shapV(2,k2) &
                      + dxdxi(1:3,3)*shapV(3,k2)
               s(1:3) = s(1:3)/rjac
!..normal component
               sn = s(1)*rn(1)+s(2)*rn(2)+s(3)*rn(3)
!
!..EXTENDED HV STIFFNESS MATRIX: -<dot(σ,n), v>_{Γ_h}
               stiff_HV(k1,k2) = stiff_HV(k1,k2) - sn*v*weight
!..end loop through H(div) trial functions
            enddo
!..end loop through enriched H1 test functions
         enddo

!..loop through enriched H(div) test functions
         do k1=1,NrdofVV
!..Piola Transform
            tau_a(1:3) = dxdxi(1:3,1) * shapVV(1,k1)  &
                       + dxdxi(1:3,2) * shapVV(2,k1)  &
                       + dxdxi(1:3,3) * shapVV(3,k1)
            tau_a(1:3) = tau_a(1:3)/rjac
            tn = tau_a(1)*rn(1) + tau_a(2)*rn(2) + tau_a(3)*rn(3)
!
            do k2=1,NrdofVi_a
               u_hat = shapH(k2)
               stiff_VH(k1,k2) = stiff_VH(k1,k2) - weight * tn * u_hat
            enddo
!
!..end loop through enriched H(div) test functions
         enddo
!..end loop through integration points
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
!
!  Total test/trial DOFs of the element
   i1 = NrdofHH; i2 = NrdofVV; j1 = NrdofVi_a; j2 = NrdofVi_b; j3 = NrdofU; j4 = NrdofS
   stiff_ALL(1:i1,1:j1) = stiff_HH
   stiff_ALL(1:i1,j1+1:j1+j2) = stiff_HV
   stiff_ALL(1:i1,j1+j2+1:j1+j2+j3) = stiff_HQ1
   stiff_ALL(1:i1,j1+j2+j3+1:j1+j2+j3+j4) = stiff_HQ2
!
   stiff_ALL(i1+1:i1+i2,1:j1) = stiff_VH
   stiff_ALL(i1+1:i1+i2,j1+1:j1+j2) = stiff_VV
   stiff_ALL(i1+1:i1+i2,j1+j2+1:j1+j2+j3) = stiff_VQ1
   stiff_ALL(i1+1:i1+i2,j1+j2+j3+1:j1+j2+j3+j4) = stiff_VQ2
!
   stiff_ALL(1:i1+i2,j1+j2+j3+j4+1) = bload_H
!
   deallocate(stiff_HQ1,stiff_HQ2,stiff_HH,stiff_HV)
   deallocate(stiff_VQ1,stiff_VQ2,stiff_VH,stiff_VV)
!
!..A. Compute Cholesky factorization of Gram Matrix, G=U^T U (=LL^T)
   call DPPTRF('U',NrTest,gramP,info)
   if (info.ne.0) then
      write(*,*) 'elem_heat: DPPTRF: Mdle,info = ',Mdle,info,'. stop.'
      stop
   endif
!
!..B. Solve triangular system to obtain B~, (LX=) U^T X = [B|l]
   call DTPTRS('U','T','N',NrTest,NrTrial+1,gramP,stiff_ALL,NrTest,info)
   if (info.ne.0) then
      write(*,*) 'elem_heat: DTPTRS: Mdle,info = ',Mdle,info,'. stop.'
      stop
   endif
!
   allocate(raloc(NrTrial+1,NrTrial+1)); raloc = ZERO

!..C. Matrix multiply: B^T G^-1 B (=B~^T B~)
   call DSYRK('U','T',NrTrial+1,NrTest,ZONE,stiff_ALL,NrTest,ZERO,raloc,NrTrial+1)
!
   deallocate(stiff_ALL)
!
!..D. Fill lower triangular part of Hermitian matrix using the upper triangular matrix
   do i=1,NrTrial
      raloc(i+1:NrTrial+1,i) = raloc(i,i+1:NrTrial+1)
   enddo
!
!  Fill the Aloc and Bloc matrices
   BlocH(1:j1) = raloc(1:j1,j1+j2+j3+j4+1)
   BlocV(1:j2) = raloc(j1+1:j1+j2,j1+j2+j3+j4+1)
   BlocQ1(1:j3) = raloc(j1+j2+1:j1+j2+j3,j1+j2+j3+j4+1)
   BlocQ2(1:j4) = raloc(j1+j2+j3+1:j1+j2+j3+j4,j1+j2+j3+j4+1)
!
   AlocHH(1:j1,1:j1)   = raloc(1:j1,1:j1)
   AlocHV(1:j1,1:j2)   = raloc(1:j1,j1+1:j1+j2)
   AlocHQ1(1:j1,1:j3)  = raloc(1:j1,j1+j2+1:j1+j2+j3)
   AlocHQ2(1:j1,1:j4)  = raloc(1:j1,j1+j2+j3+1:j1+j2+j3+j4)
!  
   AlocVH(1:j2,1:j1)   = raloc(j1+1:j1+j2,1:j1)
   AlocVV(1:j2,1:j2)   = raloc(j1+1:j1+j2,j1+1:j1+j2)
   AlocVQ1(1:j2,1:j3)  = raloc(j1+1:j1+j2,j1+j2+1:j1+j2+j3)
   AlocVQ2(1:j2,1:j4)  = raloc(j1+1:j1+j2,j1+j2+j3+1:j1+j2+j3+j4)
!
   AlocQ1H(1:j3,1:j1)  = raloc(j1+j2+1:j1+j2+j3,1:j1)
   AlocQ1V(1:j3,1:j2)  = raloc(j1+j2+1:j1+j2+j3,j1+1:j1+j2)
   AlocQ1Q1(1:j3,1:j3) = raloc(j1+j2+1:j1+j2+j3,j1+j2+1:j1+j2+j3)
   AlocQ1Q2(1:j3,1:j4) = raloc(j1+j2+1:j1+j2+j3,j1+j2+j3+1:j1+j2+j3+j4)
!
   AlocQ2H(1:j4,1:j1)  = raloc(j1+j2+j3+1:j1+j2+j3+j4,1:j1)
   AlocQ2V(1:j4,1:j2)  = raloc(j1+j2+j3+1:j1+j2+j3+j4,j1+1:j1+j2)
   AlocQ2Q1(1:j4,1:j3) = raloc(j1+j2+j3+1:j1+j2+j3+j4,j1+j2+1:j1+j2+j3)
   AlocQ2Q2(1:j4,1:j4) = raloc(j1+j2+j3+1:j1+j2+j3+j4,j1+j2+j3+1:j1+j2+j3+j4)
!
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
   use node_types, only : MDLB
!
   implicit none
!
   integer, intent(in)  :: Etype
   integer, intent(in)  :: Nord
   integer, intent(out) :: Norder(19)
!
   integer :: temp(2)
   integer :: nordF(3),nordB(3)
!
!----------------------------------------------------------------------
!
!..see implementation of BrokenExactSequence in shape functions
   select case(Etype)
   case(MDLB)
      call decod(Nord,MODORDER,2, temp) !xy face, z edge
      nordF(1) = temp(1); nordB(3) = temp(2)
      call decod(nordF(1),MODORDER,2, nordB(1:2)) !x edge, y edge
      call encod((/nordB(1),nordB(3)/),MODORDER,2, nordF(2)) !xz face
      call encod(nordB(2:3),MODORDER,2, nordF(3)) !yz face
!..edges
      Norder(1:4)   = (/nordB(1),nordB(2),nordB(1),nordB(2)/) !x,y,x,y edges
      Norder(5:8)   = (/nordB(1),nordB(2),nordB(1),nordB(2)/) !x,y,x,y edges
      Norder(9:12)  = nordB(3) !z edges
!..faces
      Norder(13:14) = nordF(1) !xy faces
      Norder(15:18) = (/nordF(2),nordF(3),nordF(2),nordF(3)/) !xz,yz,xz,yz faces
!..element interior
      Norder(19)    = Nord
   case default
      write(*,*) 'compute_enriched_order: only implemented for hexa. stop.'
      stop
end select
!
end subroutine compute_enriched_order
!

