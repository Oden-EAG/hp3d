!-------------------------------------------------------------------------
!
!     routine name      - elem_opt_UW
!
!-------------------------------------------------------------------------
!
!> @date  Apr 2024
!
!> @brief Compute element stiffness and load using BLAS3
!
!     arguments:
!        in:
!              Mdle     - element middle node number
!              NrTest   - number of total test dofs
!              NrTrial  - total number of trial dof
!              NrdofHH  - number of H1 test dof
!              NrdofVV  - number of Hdiv test dof
!              NrdofH   - number of H1 trial dof
!              NrdofQ   - number of L2 dofs which will be used for u and sigma_1,2,3 
!              NrdofU   - number of L2 trial dof for u
!              NrdofS   - number of L2 trial dof for sigma  = grad u (Dimension X NrdofQ)
!              NrdofHi  - number of H1 trial interface dof
!              NrdofVi  - number of H(div) trial interface dof
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
subroutine elem_opt_UW( Mdle,                                        &
                        NrTest,NrTrial,                              &
                        NrdofHH,NrdofVV,NrdofH,NrdofQ,               &
                        NrdofU,NrdofS,                               &
                        NrdofHi,NrdofVi,                             &
                        MdH,MdV,MdQ1,MdQ2,                           &
                        BlocH,AlocHH,AlocHV,AlocHQ1,AlocHQ2,         &
                        BlocV,AlocVH,AlocVV,AlocVQ1,AlocVQ2,         &
                        BlocQ1,AlocQ1H,AlocQ1V,AlocQ1Q1,AlocQ1Q2,    &
                        BlocQ2,AlocQ2H,AlocQ2V,AlocQ2Q1,AlocQ2Q2)
!
!..ALOC: holds local element stiffness matrices
!..BLOC: holds local element load vectors
   use common_prob_data
   use control
   use data_structure3D
   use element_data
   use parametersDPG
   use mpi_wrapper
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
   integer, intent(in) :: NrdofHi
   integer, intent(in) :: NrdofVi
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
   real(8) :: rjac, bjac, fval, wa, weight, sqrt_weight, tfac, vfac
   integer :: iflag, nrv, nre, nrf, nint
   integer :: j1, j2, j3, j4, k1, k2, i, i1, i2, k, l
   integer :: nordP, nrdof, nrdofH_ifc, nsign, ifc, info
   integer :: etype, ftype
!
!..element order, orientation for edges and faces
   integer :: norder(19), norient_edge(12), norient_face(6)
!
!..element nodes order (trial) for interfaces
   integer :: norder_ifc(19)
!
!..face order
   integer :: norderf(5)
!
!..number of face dofs
   integer :: ndofH_face, ndofE_face, ndofV_face, ndofQ_face
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
   real(8) :: shapH_ifc(MAXbrickH), gradH_ifc(3,MAXbrickH)
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
!..Gram matrix
   real(8), allocatable :: gram(:,:), gram_HH(:,:), gram_HV(:,:), gram_VV(:,:)
!
!..auxiliary matrices
   real(8), allocatable :: test_H(:,:), test_GH(:,:),   &
                           test_V(:,:), test_DV(:,:),   &
                           trial_Q1(:,:), trial_H(:,:), &
                           trial_Q2(:,:), trial_V(:,:)
   integer :: koff, noff, noffH, noffV, nda, nH, nV, nS, nU
!
!..stiffness matrices for the enriched test space
   real(8), allocatable :: stiff_HQ2(:,:),stiff_HV(:,:),stiff_VH(:,:)
   real(8), allocatable :: stiff_VQ1(:,:),stiff_VQ2(:,:)
   real(8), allocatable :: stiff_ALL(:,:),raloc(:,:)
!
!..3D quadrature data
   real(8) :: xiloc(3,MAXNINT3ADD), waloc(MAXNINT3ADD)
!
!..2D quadrature data
   real(8) :: tloc(2,MAXNINT2ADD), wtloc(MAXNINT2ADD)
!
!..workspace for trial and test variables
   real(8) :: u, v, q, tn, sn, u_hat
   real(8) :: dv(3), tau(3), s(3)
!
!..construction of DPG system assumes uplo = 'U'
   character(len=1), parameter :: uplo = 'U'
!
!..TIMER
   real(8) :: start_time, end_time
   logical, parameter :: timer = .false.
!
!---------------------------------------------------------------------
!--------- INITIALIZE THE ELEMENT ORDER, ORIENTATION ETC. ------------
!---------------------------------------------------------------------
!
!..TIMER
   if (timer) start_time = MPI_Wtime()
!
!..element type
   etype = NODES(Mdle)%ntype
   nrv = nvert(etype)
   nre = nedge(etype)
   nrf = nface(etype)
!
!..determine order of approximation
   call find_order(Mdle, norder)
   norder_ifc(1:nre+nrf) = norder(1:nre+nrf)
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
         write(*,*) 'elem_opt: invalid etype param. stop.'
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
!..allocate auxiliary matrices
   nda = 3*nint
   allocate(test_H (NrdofHH,nint))
   allocate(test_GH(NrdofHH,nda))
   allocate(test_V (NrdofVV,nda))
   allocate(test_DV(NrdofVV,nint))
   allocate(trial_Q1(NrdofU,nint))
   allocate(trial_Q2(NrdofS,nda))
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
!
!  ...L2 shape functions
      call shape3DQ(etype,xi,norder, nrdof,shapQ)
!
!  ...discontinuous H1 shape functions
      call shape3HH(etype,xi,nordP, nrdof,shapHH,gradHH)
!
!  ...discontinuous H(div) shape functions
      call shape3VV(etype,xi,nordP, nrdof,shapVV,divVV)
!
!  ...geometry map
      call geom3D(Mdle,xi,xnod,shapH,gradH,NrdofH, x,dxdxi,dxidx,rjac,iflag)
!
!  ...integration weight
      weight = rjac*wa
      sqrt_weight = sqrt(weight)
!
!  ...get the RHS
      call getf(Mdle,x, fval)
!
!  ...offset in auxiliary matrix
      noff = 3*(l-1)
!
!  ...loop through enriched H1 test functions
      do k=1,NrdofHH
!     ...Piola transformation
         v = shapHH(k)
         dv(1:3) = gradHH(1,k)*dxidx(1,1:3) &
                 + gradHH(2,k)*dxidx(2,1:3) &
                 + gradHH(3,k)*dxidx(3,1:3)
!
!     ...accumulate load vector: (f,v)
         bload_H(k) = bload_H(k) + fval*v*weight
!
!     ...fill test H1 auxiliary matrices
         test_H(k,l) = v * sqrt_weight
         test_GH(k,noff+1:noff+3) = dv(1:3) * sqrt_weight
      enddo
!
!  ...loop through enriched H(div) test functions
      do k=1,NrdofVV
!     ...Piola transformation
         tau(1:3) = dxdxi(1:3,1) * shapVV(1,k) &
                  + dxdxi(1:3,2) * shapVV(2,k) &
                  + dxdxi(1:3,3) * shapVV(3,k)
         tau(1:3) = tau(1:3) / rjac
         q = divVV(k) / rjac
!
!     ...fill test H(div) auxiliary matrices
         test_V(k,noff+1:noff+3) = tau(1:3) * sqrt_weight
         test_DV(k,l) = q * sqrt_weight
      enddo
!
!  ...loop through L2 trial functions
      do k=1,NrdofQ
!     ...Piola transformation
         u = shapQ(k) / rjac
!
!     ...fill trial L2 auxiliary matrices
         koff = 3*(k-1)
         trial_Q1(k,l) = u * sqrt_weight
!
         trial_Q2(koff+1:koff+3,noff+1:noff+3) = ZERO
         trial_Q2(koff+1,noff+1) = u * sqrt_weight
         trial_Q2(koff+2,noff+2) = u * sqrt_weight
         trial_Q2(koff+3,noff+3) = u * sqrt_weight
      enddo
!
!..end of loop through integration points
   enddo
!
!..compute Gram matrix
   nH = NrdofHH; nV = NrdofVV
   allocate(gram_HH(nH,nH))
   allocate(gram_VV(nV,nV))
!
   select case(TEST_NORM)
      case(GRAPH_NORM)
         vfac = ALPHA_NORM
         tfac = ALPHA_NORM + 1.0d0
      case( MATH_NORM)
         vfac = 1.0d0
         tfac = 1.0d0
      case default
         write(*,*) 'elem_opt: invalid test norm!'
         stop
   end select
!
!  α*(v,v)
   call DSYRK(uplo,'N',nH,nint,vfac ,test_H(1:nH,1:nint) ,nH, &
                               0.0d0,gram_HH(1:nH,1:nH)  ,nH)
!  (grad v, grad v)
   call DSYRK(uplo,'N',nH,nda ,1.0d0,test_GH(1:nH,1:nda) ,nH, &
                               1.0d0,gram_HH(1:nH,1:nH)  ,nH)
!  (α+1)*(tau, tau)
   call DSYRK(uplo,'N',nV,nda , tfac,test_V(1:nV,1:nda)  ,nV, &
                               0.0d0,gram_VV(1:nV,1:nV)  ,nV)
!  (div tau, div tau)
   call DSYRK(uplo,'N',nV,nint,1.0d0,test_DV(1:nV,1:nint),nV, &
                               1.0d0,gram_VV(1:nV,1:nV)  ,nV)
!
   if (TEST_NORM .eq. GRAPH_NORM) then
   !  (grad v, tau)
      allocate(gram_HV(nH,nV))
      call DGEMM('N','T',nH,nV,nda,1.0d0,test_GH(1:nH,1:nda),nH, &
                                         test_V(1:nV,1:nda) ,nV, &
                                   0.0d0,gram_HV(1:nH,1:nV) ,nH)
!  ...insert blocks into full Gram matrix
      allocate(gram(NrTest,NrTest))
      gram(   1:nH   ,   1:nH   ) = gram_HH(1:nH,1:nH)
      gram(   1:nH   ,nH+1:nH+nV) = gram_HV(1:nH,1:nV)
      gram(nH+1:nH+nV,nH+1:nH+nV) = gram_VV(1:nV,1:nV)
      deallocate(gram_HH,gram_HV,gram_VV)
   endif
!
!..compute stiffness matrix
   nS = NrdofS; nU = NrdofU
   allocate(stiff_HQ2(nH,nS))
   allocate(stiff_VQ1(nV,NrdofU))
   allocate(stiff_VQ2(nV,nS))
!
!  (sigma, grad v)
   call DGEMM('N','T',nH,nS,nda ,1.0d0,test_GH (1:nH,1:nda),nH,   &
                                       trial_Q2(1:nS,1:nda),nS,   &
                                 0.0d0,stiff_HQ2(1:nH,1:nS),nH)
!  (u, div tau)
   call DGEMM('N','T',nV,nU,nint,1.0d0,test_DV (1:nV,1:nint),nV,  &
                                       trial_Q1(1:nU,1:nint),nU,  &
                                 0.0d0,stiff_VQ1(1:nV,1:nU) ,nV)
!  (sigma, tau)
   call DGEMM('N','T',nV,nS,nda ,1.0d0,test_V(1:nV,1:nda)  ,nV,   &
                                       trial_Q2(1:nS,1:nda),nS ,  &
                                 0.0d0,stiff_VQ2(1:nV,1:nS),nV)
!
   deallocate(test_H,test_GH,test_V,test_DV,trial_Q1,trial_Q2)
!
!..TIMER
   if (timer) then
      end_time = MPI_Wtime()
      !$OMP CRITICAL
      write(*,11) 'elem INTEGR Vol: ', end_time-start_time
      !$OMP END CRITICAL
      11 format(A,f12.5,' s')
      start_time = MPI_Wtime()
   endif
!
!---------------------------------------------------------------------
!    B O U N D A R Y    I N T E G R A L S                            |
!---------------------------------------------------------------------
!
   allocate(stiff_HV(NrdofHH,NrdofVi)); stiff_HV = ZERO
   allocate(stiff_VH(NrdofVV,NrdofHi)); stiff_VH = ZERO
!
   allocate( test_H(NrdofHH,MAXNINT2ADD))
   allocate( test_V(NrdofVV,MAXNINT2ADD))
   allocate(trial_H(NrdofHi,MAXNINT2ADD))
   allocate(trial_V(NrdofVi,MAXNINT2ADD))
!
   noffH = 0; noffV = 0
!
!..loop through element faces
   do ifc=1,nrf
!
      trial_H = ZERO; trial_V = ZERO
!
!  ...sign factor to determine the outward normal unit vector
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
!  ...determine #dof for this face node
      call ndof_nod(face_type(etype,ifc),norder(nre+ifc), &
                    ndofH_face,ndofE_face,ndofV_face,ndofQ_face)
!
!  ...set order for all element edge nodes and this face node
!     note: modified norder_ifc returns shape functions with different
!           enumeration of DOFs than the element DOF ordering
      call initiate_order(etype, norder_ifc)
      norder_ifc(1:nre) = norder(1:nre)
      norder_ifc(nre+ifc) = norder(nre+ifc)
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
!
!     ...discontinuous H(div) shape functions
         call shape3VV(etype,xi,nordP, nrdof,shapVV,divVV)
!
!     ...determine element H1 shape functions
!        (vertex, edge, face trial dofs)
         call shape3DH(etype,xi,norder_ifc,norient_edge,norient_face, &
                       nrdofH_ifc,shapH_ifc,gradH_ifc)
!
!     ...determine H(div) shape functions (for fluxes)
         call shape3DV(etype,xi,norder_ifc,norient_face, &
                       nrdof,shapV,divV)
#if DEBUG_MODE
         if (nrdof .ne. ndofV_face+nrf-1) then
            write(*,*) 'elem_opt: nrdof = ',nrdof, &
                               ', ndofV_face = ',ndofV_face
            stop
         endif
#endif
!
!     ...determine element H1 shape functions (for geometry)
         call shape3DH(etype,xi,norder,norient_edge,norient_face, &
                       nrdof,shapH,gradH)
!
!     ...geometry map
         call bgeom3D(Mdle,xi,xnod,shapH,gradH,NrdofH,dxidt,nsign, &
                      x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
!     ...integration weight
         weight = bjac*wtloc(l)
         sqrt_weight = sqrt(weight)
!
!     ...enriched H1 test functions
         test_H(1:NrdofHH,l) = shapHH(1:NrdofHH) * sqrt_weight
!
!     ...number of unused H(div) shape functions from other faces
!        that we need to skip (one per previous face)
         nrdof = ifc-1
!
!     ...loop through H(div) trial functions
         do k=1,ndofV_face
            k1 = nrdof + k
!        ...Piola transformation
            s(1:3) = dxdxi(1:3,1)*shapV(1,k1) &
                   + dxdxi(1:3,2)*shapV(2,k1) &
                   + dxdxi(1:3,3)*shapV(3,k1)
            s(1:3) = s(1:3)/rjac
!        ...normal component
            sn = s(1)*rn(1)+s(2)*rn(2)+s(3)*rn(3)
!
            k2 = noffV + k
            trial_V(k2,l) = sn * sqrt_weight
         enddo
!
!     ...loop through enriched H(div) test functions
         do k=1,NrdofVV
!        ...Piola Transform
            tau(1:3) = dxdxi(1:3,1) * shapVV(1,k)  &
                     + dxdxi(1:3,2) * shapVV(2,k)  &
                     + dxdxi(1:3,3) * shapVV(3,k)
            tau(1:3) = tau(1:3)/rjac
            tn = tau(1)*rn(1) + tau(2)*rn(2) + tau(3)*rn(3)
!
            test_V(k,l) = tn * sqrt_weight
         enddo
!
!     ...loop through H1 trial shape functions
!        - vertex and edge shape functions
         nrdof = nrdofH_ifc - ndofH_face
         trial_H(1:nrdof,l) = shapH_ifc(1:nrdof) * sqrt_weight
!        - face shape functions
         do k=1,ndofH_face
            k1 = nrdof + k
            u_hat = shapH_ifc(k1)
!
            k2 = nrdof + noffH + k
            trial_H(k2,l) = u_hat * sqrt_weight
         enddo
!
!  ...end loop through integration points
      enddo
!
!  ...EXTENDED HV STIFFNESS MATRIX
!     -<dot(σ,n), v>_{Γ_h}
      call DGEMM('N','T',NrdofHH,NrdofVi,nint,-1.0d0, test_H (1:NrdofHH,1:nint   ),NrdofHH,  &
                                                     trial_V (1:NrdofVi,1:nint   ),NrdofVi,  &
                                               1.0d0,stiff_HV(1:NrdofHH,1:NrdofVi),NrdofHH)
!     -<û, dot(τ,n)>_{Γ_h}
      call DGEMM('N','T',NrdofVV,NrdofHi,nint,-1.0d0, test_V (1:NrdofVV,1:nint   ),NrdofVV,  &
                                                     trial_H (1:NrdofHi,1:nint   ),NrdofHi,  &
                                               1.0d0,stiff_VH(1:NrdofVV,1:NrdofHi),NrdofVV)
!
!  ...increment DOF offset for the next face
      noffH = noffH + ndofH_face
      noffV = noffV + ndofV_face
!
!..end loop through element faces
   enddo
!
   deallocate(test_H,test_V,trial_H,trial_V)
!
!..TIMER
   if (timer) then
      end_time = MPI_Wtime()
      !$OMP CRITICAL
      write(*,11) 'elem INTEGR Bdr: ', end_time-start_time
      !$OMP END CRITICAL
      start_time = MPI_Wtime()
   endif
!
!---------------------------------------------------------------------
!  Construction of DPG system
!---------------------------------------------------------------------
!
   allocate(stiff_ALL(NrTest,NrTrial+1))
   stiff_ALL = ZERO
!
!  Total test/trial DOFs of the element
   i1 = NrdofHH; i2 = NrdofVV
   j1 = NrdofHi; j2 = NrdofVi; j3 = NrdofU; j4 = NrdofS
!
   stiff_ALL(   1:i1   ,j1      +1:j1+j2      ) = stiff_HV
   stiff_ALL(   1:i1   ,j1+j2+j3+1:j1+j2+j3+j4) = stiff_HQ2
!
   stiff_ALL(i1+1:i1+i2,         1:j1         ) = stiff_VH
   stiff_ALL(i1+1:i1+i2,j1+j2   +1:j1+j2+j3   ) = stiff_VQ1
   stiff_ALL(i1+1:i1+i2,j1+j2+j3+1:j1+j2+j3+j4) = stiff_VQ2
!
   stiff_ALL(1:i1+i2,j1+j2+j3+j4+1) = bload_H
!
   deallocate(stiff_HQ2,stiff_HV)
   deallocate(stiff_VQ1,stiff_VQ2,stiff_VH)
!
!..A. Compute Cholesky factorization of Gram Matrix, G=U^T U (=LL^T)
   if (TEST_NORM .eq. GRAPH_NORM) then
      call DPOTRF(uplo,NrTest,gram,NrTest,info)
   else
      call DPOTRF(uplo,nH,gram_HH,nH,info); if (info.ne.0) goto 310
      call DPOTRF(uplo,nV,gram_VV,nV,info)
   endif
   310 continue
   if (info.ne.0) then
      write(*,*) 'elem_opt_UW: DPOTRF: Mdle,info = ',Mdle,info,'. stop.'
      stop
   endif
!
!..B. Solve triangular system to obtain B~, (LX=) U^T X = [B|l]
   if (TEST_NORM .eq. GRAPH_NORM) then
      call DTRTRS(uplo,'T','N',NrTest,NrTrial+1,gram,NrTest,stiff_ALL,NrTest,info)
      deallocate(gram)
   else
      call DTRTRS(uplo,'T','N',nH,NrTrial+1,gram_HH,nH,stiff_ALL(   1:nH   ,:),nH,info)
      if (info.ne.0) goto 320
      call DTRTRS(uplo,'T','N',nV,NrTrial+1,gram_VV,nV,stiff_ALL(nH+1:nH+nV,:),nV,info)
      deallocate(gram_HH,gram_VV)
   endif
   320 continue
   if (info.ne.0) then
      write(*,*) 'elem_opt_UW: DTRTRS: Mdle,info = ',Mdle,info,'. stop.'
      stop
   endif
!
   allocate(raloc(NrTrial+1,NrTrial+1)); raloc = ZERO

!..C. Matrix multiply: B^T G^-1 B (=B~^T B~)
   call DSYRK(uplo,'T',NrTrial+1,NrTest,ZONE,stiff_ALL,NrTest,ZERO,raloc,NrTrial+1)
!
   deallocate(stiff_ALL)
!
!..D. Fill lower triangular part of Hermitian matrix using the upper triangular matrix
   do i=1,NrTrial
      raloc(i+1:NrTrial+1,i) = raloc(i,i+1:NrTrial+1)
   enddo
!
!..Fill the ALOC and BLOC matrices
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
!..TIMER
   if (timer) then
      end_time = MPI_Wtime()
      !$OMP CRITICAL
      write(*,11) 'elem DPG LinAlg: ', end_time-start_time
      !$OMP END CRITICAL
   endif
!
end subroutine elem_opt_UW
