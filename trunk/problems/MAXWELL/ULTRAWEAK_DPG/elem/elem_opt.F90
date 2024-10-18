!------------------------------------------------------------------------------
!> @brief      Evaluates unconstrained stiffness matrix and load vector
!!             using BLAS3
!!
!> @param[in]  Mdle     - element middle node number
!> @param[in]  NrTest   - total number of test dof
!> @param[in]  NrTrial  - total number of trial dof
!> @param[in]  NrdofEE  - number of H(curl) test dof
!> @param[in]  NrdofH   - number of H1 trial dof
!> @param[in]  NrdofE   - number of H(curl) trial dof
!> @param[in]  NrdofQ   - number of L2 trial dof
!> @param[in]  NrdofEi  - number of H(curl) trial interface dof
!> @param[in]  MdE      - num rows of ZalocEE,ZalocEQ
!> @param[in]  MdQ      - num rows of ZalocQE,ZalocQQ
!!
!> @param[out] ZblocE   - load vectors
!> @param[out] ZblocQ
!> @param[out] ZalocEE  - stiffness matrices
!> @param[out] ZalocEQ
!> @param[out] ZalocQE
!> @param[out] ZalocQQ
!!
!> @date       Apr 2024
!------------------------------------------------------------------------------
subroutine elem_opt(Mdle,                         &
                    NrTest,NrTrial,               &
                    NrdofEE,                      &
                    NrdofH,NrdofE,NrdofQ,         &
                    NrdofEi,                      &
                    MdE,MdQ,                      &
                    ZblocE,ZalocEE,ZalocEQ,       &
                    ZblocQ,ZalocQE,ZalocQQ)
!
   use control
   use parametersDPG
   use data_structure3D
   use commonParam
   use mpi_wrapper
!
   implicit none
!
!..declare input/output variables
   integer,                        intent(in)  :: Mdle
   integer,                        intent(in)  :: NrTest
   integer,                        intent(in)  :: NrTrial
   integer,                        intent(in)  :: NrdofEE
   integer,                        intent(in)  :: NrdofH
   integer,                        intent(in)  :: NrdofE
   integer,                        intent(in)  :: NrdofQ
   integer,                        intent(in)  :: NrdofEi
   integer,                        intent(in)  :: MdE
   integer,                        intent(in)  :: MdQ
   complex(8), dimension(MdE),     intent(out) :: ZblocE
   complex(8), dimension(MdE,MdE), intent(out) :: ZalocEE
   complex(8), dimension(MdE,MdQ), intent(out) :: ZalocEQ
   complex(8), dimension(MdQ),     intent(out) :: ZblocQ
   complex(8), dimension(MdQ,MdE), intent(out) :: ZalocQE
   complex(8), dimension(MdQ,MdQ), intent(out) :: ZalocQQ
!
!..declare type variables
   integer :: etype, ftype
!
!..declare element order, orientation for edges and faces
   integer :: norder(19), norient_edge(12), norient_face(6)
!
!..element nodes order (trial) for interfaces
   integer :: norderi(19), norder_ifc(19), norderf(5)
!
!..number of face dofs
   integer :: ndofH_face, ndofE_face, ndofV_face, ndofQ_face
!
!..geometry dof (work space for nodcor)
   real(8) :: xnod(3,MAXbrickH)
!
!..geometry
   real(8) :: xi(3), x(3), rn(3)
   real(8) :: dxidt(3,2), dxdt(3,2)
   real(8) :: dxdxi(3,3), dxidx(3,3)
   real(8) :: t(2), rjac, bjac
!
!..H1 shape functions
   real(8) :: shapH(MAXbrickH), gradH(3,MAXbrickH)
!
!..H(curl) shape functions
   real(8) :: shapE (3,MAXbrickE), curlE(3,MAXbrickE)
   real(8) :: shapFi(3,MAXbrickE)
!
   real(8) :: shapF_T(NrdofEE,3)
!
!  enriched
   real(8) :: shapEE(3,MAXbrickEE), curlEE(3,MAXbrickEE)
!
!..L2 shape functions
   real(8) :: shapQ(MAXbrickQ)
!
!..Field and material values
   real(8) :: u
   real(8) :: fldF(3), crlF(3)
   complex(8) :: eps(3,3)
!
!..auxiliary matrices
   real(8), allocatable :: data_rEPS(:,:),data_rMU(:,:)
   real(8), allocatable :: data_iEPS(:,:),data_iMU(:,:)
!
   real(8), allocatable ::  temp_rQ (:,:), temp_iQ (:,:)
   real(8), allocatable ::  temp_rE (:,:), temp_iE (:,:)
   real(8), allocatable :: stiff_rFE(:,:),stiff_iFE(:,:)
   real(8), allocatable :: stiff_rGH(:,:),stiff_iGH(:,:)
!
   real(8), allocatable :: gram_r(:,:), gram_i(:,:)
   real(8), allocatable :: stiff_rFH(:,:), test_rE (:,:),test_rCE(:,:)
   real(8), allocatable :: stiff_rEE(:,:),trial_rnE(:,:),trial_rQ(:,:)
   real(8), allocatable :: stiff_rFE_imp(:,:),stiff_rFH_imp(:,:),trial_rnnE(:,:)
!
   complex(8), allocatable :: temp_E(:,:), gram_FF(:,:), gram_FG(:,:), gram_GG(:,:)
   complex(8), allocatable :: stiff_ALL(:,:), zaloc(:,:), gram(:,:)
!
!..load vector for the enriched space
   complex(8) :: bload_E(NrTest)
!
   integer :: ioff, joff, koff, noff, noffE, nda, nE, nEi, nQ
!
!..quadrature data
   real(8) :: xiloc(3,MAXNINT3ADD), waloc(MAXNINT3ADD)
   real(8) :: tloc(2,MAXNINT2ADD), wtloc(MAXNINT2ADD)
   real(8) :: sqrt_weight, weight, wa
!
!..BC's flags
   integer :: ibc(6,NRINDEX)
!
!..Maxwell load and auxiliary variables
   complex(8) :: zJ(3), zImp(3)
   real(8), dimension(3) :: E2,rntimesE,rn2timesE
!
!..number of vertices, edges, and faces per element type
   integer :: nrv, nre, nrf
!
!..auxiliary variables
   integer :: jE, jQ, k2, i, j, k, l, nint, n, ifc
   integer :: iflag, info, nrdof, nrdofE_ifc, nordP, nsign
   real(8) :: afac
   logical :: imp_bc
!
   complex(8) :: za(3,3), zc(3,3)
!
!..construction of DPG system assumes uplo = 'U'
   character(len=1), parameter :: uplo = 'U'
!
!..define DPG linear system block-wise
   logical, parameter :: blocks = .false.
!
!..TIMER
   real(8) :: start_time, end_time
   logical, parameter :: timer = .false.
!
#if HP3D_DEBUG
!..Set iprint = 0/1 (Non-/VERBOSE)
   integer :: iprint
   iprint = 0
#endif
!
!-------------------------------------------------------------------------------
!
!..TIMER
   if (timer) start_time = MPI_Wtime()
!
#if HP3D_DEBUG
   if (iprint.eq.1) then
      write(*,*) 'elem_opt: Mdle = ', Mdle
   endif
#endif
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
!..get the element boundary conditions flags
   call find_bc(Mdle, ibc)
!
!..clear space for output matrices
   ZblocE(:)    = ZERO; ZblocQ(:)    = ZERO
   ZalocEE(:,:) = ZERO; ZalocEQ(:,:) = ZERO
   ZalocQE(:,:) = ZERO; ZalocQQ(:,:) = ZERO
!
!..clear space for auxiliary matrices
   bload_E(:) = ZERO
!
!-----------------------------------------------------------------------
!              E L E M E N T   I N T E G R A L S
!-----------------------------------------------------------------------
!
!..use the enriched order to set the quadrature
   INTEGRATION = NORD_ADD
   call set_3D_int_DPG(etype,norder,norient_face, nint,xiloc,waloc)
   INTEGRATION = 0
!
!..allocate auxiliary matrices
   nda = 3*nint
   allocate(test_rE (NrdofEE,nda))
   allocate(test_rCE(NrdofEE,nda))
   allocate(trial_rQ(3*NrdofQ,nda))
   allocate(data_rEPS(3,nda))
   allocate(data_iEPS(3,nda))
   allocate(data_rMU (3,nda))
   allocate(data_iMU (3,nda))
!
!..loop over integration points
   do l=1,nint
!
      xi(1:3)=xiloc(1:3,l); wa=waloc(l)
!
!  ...H1 shape functions (for geometry)
      call shape3DH(etype,xi,norder,norient_edge,norient_face, nrdof,shapH,gradH)
!
!  ...L2 shape functions for the trial space
      call shape3DQ(etype,xi,norder, nrdof,shapQ)
!
!  ...broken H(curl) shape functions for the enriched test space
      call shape3EE(etype,xi,nordP, nrdof,shapEE,curlEE)
!
!  ...geometry map
      call geom3D(Mdle,xi,xnod,shapH,gradH,NrdofH, x,dxdxi,dxidx,rjac,iflag)
#if HP3D_DEBUG
      if (iflag .ne. 0) then
         write(*,5999) Mdle,rjac
5999     format('elem_opt: Negative Jacobian. Mdle,rjac=',i8,2x,e12.5)
         stop
      endif
#endif
!
!  ...get permittivity at x
      call get_permittivity(Mdle,x, eps)
!
      za = (ZI*OMEGA*EPSILON) * eps(1:3,1:3)
      zc = (ZI*OMEGA*MU) * IDENTITY(1:3,1:3)
!
!  ...offset in auxiliary matrix
      noff = 3*(l-1)
!
      data_rEPS(1:3,noff+1:noff+3) = real (za(1:3,1:3))
      data_iEPS(1:3,noff+1:noff+3) = aimag(za(1:3,1:3))
      data_rMU (1:3,noff+1:noff+3) = real (zc(1:3,1:3))
      data_iMU (1:3,noff+1:noff+3) = aimag(zc(1:3,1:3))
!
!  ...integration weight
      weight = rjac*wa
      sqrt_weight = sqrt(weight)
!
!  ...get the RHS
      call getf(Mdle,x, zJ)
!
!  ...loop through L2 trial shape functions
      do k=1,NrdofQ
!     ...Piola transformation
         u = shapQ(k)/rjac
!
!     ...fill trial L2 auxiliary matrix
         koff = 3*(k-1)
         trial_rQ(koff+1:koff+3,noff+1:noff+3) = 0.d0
         trial_rQ(koff+1,noff+1) = u * sqrt_weight
         trial_rQ(koff+2,noff+2) = u * sqrt_weight
         trial_rQ(koff+3,noff+3) = u * sqrt_weight
      enddo
!
!  ...loop through enriched H(curl) test functions
      do k=1,NrdofEE
!
!     ...Piola transformation
         fldF(1:3) = shapEE(1,k)*dxidx(1,1:3) &
                   + shapEE(2,k)*dxidx(2,1:3) &
                   + shapEE(3,k)*dxidx(3,1:3)
         crlF(1:3) = dxdxi(1:3,1)*curlEE(1,k) &
                   + dxdxi(1:3,2)*curlEE(2,k) &
                   + dxdxi(1:3,3)*curlEE(3,k)
         crlF(1:3) = crlF(1:3)/rjac
!
!        RHS:
!        (J^imp,F) first  equation RHS (with first H(curl) test function F)
!        (0    ,G) second equation RHS is zero
         if (blocks) then
            n = k
         else
            n = 2*k-1
         endif
         bload_E(n) = bload_E(n)                                   &
                    + (fldF(1)*zJ(1)+fldF(2)*zJ(2)+fldF(3)*zJ(3))  &
                    * weight
!
!     ...fill test H(curl) auxiliary matrices
         test_rE (k,noff+1:noff+3) = fldF(1:3) * sqrt_weight
         test_rCE(k,noff+1:noff+3) = crlF(1:3) * sqrt_weight
!
!  ...end of loop through enriched H(curl) test functions
      enddo
!
!..end of loop through integration points
   enddo
!
!  --- Stiffness matrix ---
   nE =   NrdofEE
   nQ = 3*NrdofQ
   allocate(stiff_rFE(nE,nQ))
   allocate(stiff_iFE(nE,nQ))
   allocate(stiff_rFH(nE,nQ))
   allocate(stiff_rGH(nE,nQ))
   allocate(stiff_iGH(nE,nQ))
   allocate(temp_rQ  (nQ,nda))
   allocate(temp_iQ  (nQ,nda))
!
!     ...iωε * E
   do l=1,nint
      noff = 3*(l-1)
      call DGEMM('N','T',nQ,3,3,1.d0,trial_rQ (1:nQ,noff+1:noff+3),nQ, &
                                     data_rEPS(1:3 ,noff+1:noff+3), 3, &
                                0.d0, temp_rQ (1:nQ,noff+1:noff+3),nQ) ! real part
      call DGEMM('N','T',nQ,3,3,1.d0,trial_rQ (1:nQ,noff+1:noff+3),nQ, &
                                     data_iEPS(1:3 ,noff+1:noff+3), 3, &
                                0.d0, temp_iQ (1:nQ,noff+1:noff+3),nQ) ! imag part
   enddo
!
!     ...-((iωε)E,F)
   call DGEMM('N','T',nE,nQ,nda,-1.d0, test_rE(1:nE,1:nda),nE,  &
                                       temp_rQ(1:nQ,1:nda),nQ,  &
                                 0.d0,stiff_rFE(1:nE,1:nQ),nE) ! real part
   call DGEMM('N','T',nE,nQ,nda,-1.d0, test_rE(1:nE,1:nda),nE,  &
                                       temp_iQ(1:nQ,1:nda),nQ,  &
                                 0.d0,stiff_iFE(1:nE,1:nQ),nE) ! imag part
!
!     ...(H,curl(F))
!     ...(E,curl(G))
   call DGEMM('N','T',nE,nQ,nda,1.d0,test_rCE(1:nE,1:nda),nE,  &
                                     trial_rQ(1:nQ,1:nda),nQ,  &
                                0.d0,stiff_rFH(1:nE,1:nQ),nE) ! real-valued
!
!     ...iωμ * H
   do l=1,nint
      noff = 3*(l-1)
      call DGEMM('N','T',nQ,3,3,1.d0,trial_rQ (1:nQ,noff+1:noff+3),nQ, &
                                      data_rMU(1:3 ,noff+1:noff+3), 3, &
                                0.d0, temp_rQ (1:nQ,noff+1:noff+3),nQ) ! real part
      call DGEMM('N','T',nQ,3,3,1.d0,trial_rQ (1:nQ,noff+1:noff+3),nQ, &
                                      data_iMU(1:3 ,noff+1:noff+3), 3, &
                                0.d0, temp_iQ (1:nQ,noff+1:noff+3),nQ) ! imag part
   enddo
!     ...((iωμ)H,G)
   call DGEMM('N','T',nE,nQ,nda,1.d0, test_rE(1:nE,1:nda),nE, &
                                      temp_rQ(1:nQ,1:nda),nQ, &
                                0.d0,stiff_rGH(1:nE,1:nQ),nE) ! real part
   call DGEMM('N','T',nE,nQ,nda,1.d0, test_rE(1:nE,1:nda),nE, &
                                      temp_iQ(1:nQ,1:nda),nQ, &
                                0.d0,stiff_iGH(1:nE,1:nQ),nE) ! imag part
!
   deallocate(trial_rQ,temp_rQ,temp_iQ)
!
!  --- Gram matrix ---
!
   select case(TEST_NORM)
      case(GRAPH_NORM); afac = ALPHA_NORM
      case(GRAPH_DIAG); afac = ALPHA_NORM
      case( MATH_NORM); afac = 1.d0
      case default
         write(*,*) 'elem_opt: invalid test norm!'
         stop
   end select
!
   allocate(gram_r(nE,nE))
!
!     ...α*(F,F)
!     ...α*(G,G)
   call DSYRK(uplo,'N',nE,nda,afac,test_rE(1:nE,1:nda),nE,  &
                              0.d0,gram_r (1:nE,1:nE ),nE) ! real-valued
!     ...(curl F, curl F)
!     ...(curl G, curl G)
   call DSYRK(uplo,'N',nE,nda,1.d0,test_rCE(1:nE,1:nda),nE,  &
                              1.d0,gram_r  (1:nE,1:nE ),nE) ! real-valued
!
   if (TEST_NORM .eq. MATH_NORM) goto 110
!
   allocate(gram_FF(nE,nE)); gram_FF = cmplx(gram_r,0.d0,8)
   allocate(gram_GG(nE,nE)); gram_GG = cmplx(gram_r,0.d0,8)
!
   allocate(temp_rE(nE,nda))
   allocate(temp_iE(nE,nda))
   allocate(temp_E (nE,nda))
!
!     ...((iωε)^* F)^* = F^T (iωε)
   do l=1,nint
      noff = 3*(l-1)
      call DGEMM('N','N',nE,3,3,1.d0,test_rE  (1:nE,noff+1:noff+3),nE,  &
                                     data_rEPS(1:3 ,noff+1:noff+3), 3,  &
                                0.d0,temp_rE  (1:nE,noff+1:noff+3),nE) ! real part
      call DGEMM('N','N',nE,3,3,1.d0,test_rE  (1:nE,noff+1:noff+3),nE,  &
                                     data_iEPS(1:3 ,noff+1:noff+3), 3,  &
                                0.d0,temp_iE  (1:nE,noff+1:noff+3),nE) ! imag part
   enddo
   temp_E = cmplx(temp_rE,temp_iE,8)
!
!     ...((iωε)^* F,(iωε)^* F) <-- F^T (iωε) (F^T (iωε))^*
   call ZHERK(uplo,'N',nE,nda,ZONE,temp_E (1:nE,1:nda),nE,  &
                              ZONE,gram_FF(1:nE,1:nE ),nE) ! complex-valued
!
!  Cross-term: upper triangular part (F,G) only
!     ...(curl G,-(iωε)^* F) <-- -F^T (iωε) curl G
   if (TEST_NORM .ne. GRAPH_DIAG) then
      allocate(gram_i(nE,nE))
      call DGEMM('N','T',nE,nE,nda,-1.d0,temp_rE (1:nE,1:nda),nE,  &
                                         test_rCE(1:nE,1:nda),nE,  &
                                    0.d0,gram_r  (1:nE,1:nE ),nE) ! real part
      call DGEMM('N','T',nE,nE,nda,-1.d0,temp_iE (1:nE,1:nda),nE,  &
                                         test_rCE(1:nE,1:nda),nE,  &
                                    0.d0,gram_i  (1:nE,1:nE ),nE) ! imag part
   endif
!
!     ...((iωμ)^* G)^* = G^T (iωμ)
   do l=1,nint
      noff = 3*(l-1)
      call DGEMM('N','N',nE,3,3,1.d0,test_rE (1:nE,noff+1:noff+3),nE,  &
                                     data_rMU(1:3 ,noff+1:noff+3), 3,  &
                                0.d0,temp_rE (1:nE,noff+1:noff+3),nE) ! real part
      call DGEMM('N','N',nE,3,3,1.d0,test_rE (1:nE,noff+1:noff+3),nE,  &
                                     data_iMU(1:3 ,noff+1:noff+3), 3,  &
                                0.d0,temp_iE (1:nE,noff+1:noff+3),nE) ! imag part
   enddo
!     ...((iωμ)^* G,(iωμ)^* G) <-- G^T (iωμ) (G^T (iωμ))^*
   temp_E = cmplx(temp_rE,temp_iE,8)
   call ZHERK(uplo,'N',nE,nda,ZONE,temp_E (1:nE,1:nda),nE,  &
                              ZONE,gram_GG(1:nE,1:nE ),nE) ! complex-valued
!
!  Cross-term: upper triangular part (F,G) only
!     ...((iωμ)^* G, curl F)  <-- (curl F)^T (iωμ)^* G = (curl F)^T (G^T (iωμ))^*
   if (TEST_NORM .ne. GRAPH_DIAG) then
      call DGEMM('N','T',nE,nE,nda, 1.d0,test_rCE(1:nE,1:nda),nE,  &
                                         temp_rE (1:nE,1:nda),nE,  &
                                    1.d0,gram_r  (1:nE,1:nE ),nE) ! real part
      call DGEMM('N','T',nE,nE,nda,-1.d0,test_rCE(1:nE,1:nda),nE,  &
                                         temp_iE (1:nE,1:nda),nE,  &
                                    1.d0,gram_i  (1:nE,1:nE ),nE) ! imag part
      allocate(gram_FG(nE,nE)); gram_FG = cmplx(gram_r,gram_i,8)
      deallocate(gram_i)
   endif
!
   deallocate(gram_r,temp_rE,temp_iE,temp_E)
!
   110 continue
!
   deallocate(data_rEPS,data_rMU,data_iEPS,data_iMU)
   deallocate(test_rE,test_rCE)
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
!-----------------------------------------------------------------------
!              B O U N D A R Y   I N T E G R A L S
!-----------------------------------------------------------------------
!
!..allocate auxiliary matrix
   allocate(stiff_rEE(NrdofEE,NrdofEi)); stiff_rEE = 0.d0
!
   nda = 3 * MAXNINT2ADD
   allocate( test_rE (NrdofEE,nda))
   allocate(trial_rnE(NrdofEi,nda))
!
!..check for impedance BC
   imp_bc = .false.
   do ifc=1,nrf
      if (ibc(ifc,2).eq.3) imp_bc = .true.
   enddo
   if (imp_bc) then
      write(*,*) 'elem_opt: Impedance BC needs to be verified.'; stop
      allocate(stiff_rFE_imp(NrdofEE,NrdofEi)); stiff_rFE_imp = 0.d0
      allocate(stiff_rFH_imp(NrdofEE,NrdofEi)); stiff_rFH_imp = 0.d0
      allocate(trial_rnnE(NrdofEi,nda))
   endif
!
   noffE = 0
!
!..loop through element faces
   do ifc=1,nrf
!
      trial_rnE = 0.d0
      if (ibc(ifc,2).eq.3) trial_rnnE = 0.d0
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
!     ...determine discontinuous Hcurl shape functions
         call shape3EE(etype,xi,nordP, nrdof,shapEE,curlEE)
#if HP3D_DEBUG
         if (nrdof .ne. NrdofEE) then
            write(*,*) 'elem_opt: INCONSISTENCY NrdofEE. stop.'
            stop
         endif
#endif
!
!     ...determine element H1 shape functions (for geometry)
         call shape3DH(etype,xi,norder,norient_edge,norient_face, &
                       nrdof,shapH,gradH)
#if HP3D_DEBUG
      if (nrdof .ne. NrdofH) then
         write(*,*) 'elem_opt: INCONSISTENCY NrdofH. stop.'
         stop
      endif
#endif
!
!     ...geometry
         call bgeom3D(Mdle,xi,xnod,shapH,gradH,NrdofH,dxidt,nsign, &
                      x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
         weight = bjac*wtloc(l)
         sqrt_weight = sqrt(weight)
!
!     ...determine element H(curl) shape functions (for fluxes)
!     ...for interfaces only (no bubbles)
         call shape3DE(etype,xi,norder_ifc,norient_edge,norient_face, &
                       nrdofE_ifc,shapE,curlE)
!
!     ...pullback trial and test functions
         call piola_hcurl_trans(dxidx,NrdofEE,shapEE(1:3,1:NrdofEE), &
                                shapF_T(1:NrdofEE,1:3))
         call piola_hcurl(dxidx,nrdofE_ifc,shapE(1:3,1:nrdofE_ifc), &
                          shapFi(1:3,1:nrdofE_ifc))
!
!     ...check for impedance BC
!        ( < n x H , F > = GAMMA*< n x n x E , F > + < zg , F > )
         if (ibc(ifc,2).eq.3) then
!        ...get the boundary source [zImp should return ZERO here]
            call get_bdSource(Mdle,x,rn, zImp)
            if (blocks) then
               bload_E(1:NrdofEE) = bload_E(1:NrdofEE) - weight       &
                                  * ( zImp(1)*shapF_T(1:NrdofEE,1) +  &
                                      zImp(2)*shapF_T(1:NrdofEE,2) +  &
                                      zImp(3)*shapF_T(1:NrdofEE,3) )
            else
               do k=1,NrdofEE
                  k2 = 2*k-1
                  bload_E(k2) = bload_E(k2) - weight      &
                              * ( zImp(1)*shapF_T(k,1) +  &
                                  zImp(2)*shapF_T(k,2) +  &
                                  zImp(3)*shapF_T(k,3) )
               enddo
            endif
         endif
!
         noff = 3*(l-1)
!
!     ...fill enriched test H(curl) auxiliary matrix
         test_rE(1:NrdofEE,noff+1:noff+3) = shapF_T(1:NrdofEE,1:3) * sqrt_weight
!
!     ...loop through H(curl) trial functions
         nrdof = nrdofE_ifc - ndofE_face
         do k=1,nrdofE_ifc
            E2(1:3) = shapFi(:,k)
            call cross_product(rn,E2, rntimesE)
            if (k .le. nrdof) then
               k2 = k ! edge shape function
            else
               k2 = noffE + k ! face shape function
            endif
!        ...fill trial H(curl) auxiliary matrix
            trial_rnE(k2,noff+1:noff+3) = rntimesE(1:3) * sqrt_weight
!
!        ...check for impedance BC (elimination strategy)
!           (impedance constant is GAMMA for TE10 mode in rectangular waveguide)
!           ( < n x H , F > = GAMMA*< n x n x E , F > + < zg , F > )
            if (ibc(ifc,2).eq.3) then
               call cross_product(rn,rntimesE, rn2timesE)
               trial_rnnE(k2,noff+1:noff+3) = rn2timesE(1:3) * GAMMA * sqrt_weight
            endif
!
!     ...end loop through H(curl) trial functions
         enddo
!  ...end loop through integration points
      enddo
!
!  ...increment DOF offset for the next face
      noffE = noffE + ndofE_face
!
!     Submatrices:
!           ZERO     <n x H, F>
!        <n x E, G>     ZERO
      nda = 3*nint
      call DGEMM('N','T',NrdofEE,NrdofEi,nda,1.0d0, test_rE (1:NrdofEE,1:nda    ),NrdofEE,  &
                                                   trial_rnE(1:NrdofEi,1:nda    ),NrdofEi,  &
                                             1.0d0,stiff_rEE(1:NrdofEE,1:NrdofEi),NrdofEE)
!     With impedance BC, we have three different submatrices:
!     <n x E    , G> on all faces
!     <n x n x E, F> on impedance faces
!     <n x H    , F> on non-impedance faces
      if (.not. imp_bc) cycle
      if (ibc(ifc,2).eq.3) then
         call DGEMM('N','T',NrdofEE,NrdofEi,nda,1.0d0, test_rE     (1:NrdofEE,1:nda    ),NrdofEE,  &
                                                      trial_rnnE   (1:NrdofEi,1:nda    ),NrdofEi,  &
                                                1.0d0,stiff_rFE_imp(1:NrdofEE,1:NrdofEi),NrdofEE)
      else
         call DGEMM('N','T',NrdofEE,NrdofEi,nda,1.0d0, test_rE     (1:NrdofEE,1:nda    ),NrdofEE,  &
                                                      trial_rnE    (1:NrdofEi,1:nda    ),NrdofEi,  &
                                                1.0d0,stiff_rFH_imp(1:NrdofEE,1:NrdofEi),NrdofEE)
      endif
!
!..end loop through faces
   enddo
!
   deallocate(test_rE,trial_rnE)
   if (imp_bc) deallocate(trial_rnnE)
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
!-------------------------------------------------------------------------------
!      Construction of the DPG system
!-------------------------------------------------------------------------------
!
   allocate(stiff_ALL(NrTest,NrTrial+1)); stiff_ALL = ZERO
!
!..Total test/trial DOFs of the element
   jE = 2*NrdofEi; jQ = 6*NrdofQ
!
   nEi =  NrdofEi
   nE  =  NrdofEE
   nQ = 3*NrdofQ
!
!..Set up the DPG linear system from the pre-computed blocks
   if (blocks) then
      ! 1st H(curl) comp.
      if (imp_bc) then
         ! GAMMA <F,n x n x E> (impedance faces)
         stiff_ALL(   1:nE   ,1:nEi) = cmplx(stiff_rFE_imp(1:nE,1:nEi),0.d0,8)
         ! <F,n x H> (non-impedance faces)
         stiff_ALL(nE+1:nE+nE,1:nEi) = cmplx(stiff_rFH_imp(1:nE,1:nEi),0.d0,8)
      else
         stiff_ALL(   1:nE   ,1:nEi) = ZERO
         stiff_ALL(nE+1:nE+nE,1:nEi) = cmplx(stiff_rEE(1:nE,1:nEi),0.d0,8)
      endif
      ! 2nd H(curl) comp
      stiff_ALL(   1:nE   ,nEi+1:nEi+nEi) = cmplx(stiff_rEE(1:nE,1:nEi),0.d0,8)
      stiff_ALL(nE+1:nE+nE,nEi+1:nEi+nEi) = ZERO
      ! 1st-3rd L2 comp
      stiff_ALL(   1:nE   ,jE   +1:jE+nQ   ) = cmplx(stiff_rFE(1:nE,1:nQ),stiff_iFE(1:nE,1:nQ),8)
      stiff_ALL(nE+1:nE+nE,jE   +1:jE+nQ   ) = cmplx(stiff_rFH(1:nE,1:nQ),0.d0                ,8)
      ! 4th-6th L2 comp
      stiff_ALL(   1:nE   ,jE+nQ+1:jE+nQ+nQ) = cmplx(stiff_rFH(1:nE,1:nQ),0.d0                ,8)
      stiff_ALL(nE+1:nE+nE,jE+nQ+1:jE+nQ+nQ) = cmplx(stiff_rGH(1:nE,1:nQ),stiff_iGH(1:nE,1:nQ),8)
      !
   else
!  ...We must extract 2x2 submatrices of the form
!         <F,E>, <F,H>
!         <G,E>, <G,H>
      do j = 1,NrdofEi
         joff = 2*(j-1)
         do i = 1,NrdofEE
            ioff = 2*(i-1)
            if (imp_bc) then
               ! GAMMA <F,n x n x E> (impedance faces)
               stiff_ALL(ioff+1,joff+1) = cmplx(stiff_rFE_imp(i,j),0.d0,8)
               ! <F,n x H> (non-impedance faces)
               stiff_ALL(ioff+1,joff+2) = cmplx(stiff_rFH_imp(i,j),0.d0,8)
            else
               stiff_ALL(ioff+1,joff+1) = ZERO
               stiff_ALL(ioff+1,joff+2) = cmplx(stiff_rEE(i,j),0.d0,8) ! <F,n x H>
            endif
            stiff_ALL(ioff+2,joff+1) = cmplx(stiff_rEE(i,j),0.d0,8) ! <G,n x E>
            stiff_ALL(ioff+2,joff+2) = ZERO
         enddo
      enddo
!
!  ...We must extract 2x6 submatrices of the form
!         (F,E_1 E_2 E_3), (F,H_1 H_2 H_3)
!         (G,E_1 E_2 E_3), (G,H_1 H_2 H_3)
!     instead of passing the computed blocks, because the problem
!     is set up with the L2 variables defined component-wise, i.e.,
!     1 variable / 6 components instead of 6 variables / 1 component
      do j = 1,NrdofQ
         joff = jE + 6*(j-1)
         noff = 3*(j-1)
         do i = 1,NrdofEE
            ioff = 2*(i-1)
            !  (F,E), (F,H)
            stiff_ALL(ioff+1,joff+1) = cmplx(stiff_rFE(i,noff+1),stiff_iFE(i,noff+1),8)
            stiff_ALL(ioff+1,joff+2) = cmplx(stiff_rFE(i,noff+2),stiff_iFE(i,noff+2),8)
            stiff_ALL(ioff+1,joff+3) = cmplx(stiff_rFE(i,noff+3),stiff_iFE(i,noff+3),8)
            stiff_ALL(ioff+1,joff+4) = cmplx(stiff_rFH(i,noff+1),0.d0,8)
            stiff_ALL(ioff+1,joff+5) = cmplx(stiff_rFH(i,noff+2),0.d0,8)
            stiff_ALL(ioff+1,joff+6) = cmplx(stiff_rFH(i,noff+3),0.d0,8)
            !  (G,E), (G,H)
            stiff_ALL(ioff+2,joff+1) = cmplx(stiff_rFH(i,noff+1),0.d0,8)
            stiff_ALL(ioff+2,joff+2) = cmplx(stiff_rFH(i,noff+2),0.d0,8)
            stiff_ALL(ioff+2,joff+3) = cmplx(stiff_rFH(i,noff+3),0.d0,8)
            stiff_ALL(ioff+2,joff+4) = cmplx(stiff_rGH(i,noff+1),stiff_iGH(i,noff+1),8)
            stiff_ALL(ioff+2,joff+5) = cmplx(stiff_rGH(i,noff+2),stiff_iGH(i,noff+2),8)
            stiff_ALL(ioff+2,joff+6) = cmplx(stiff_rGH(i,noff+3),stiff_iGH(i,noff+3),8)
         enddo
      enddo
!..end if define linear system block-wise
   endif
!
!..load vector
   stiff_ALL(1:NrTest,jE+jQ+1) = bload_E(1:NrTest)
!
   if (imp_bc) deallocate(stiff_rFE_imp,stiff_rFH_imp)
   deallocate(stiff_rFE,stiff_iFE,stiff_rFH,stiff_rGH,stiff_iGH,stiff_rEE)
!
!..A. Compute Cholesky factorization of Gram Matrix, G=U^*U (=LL^*)
   allocate(gram(NrTest,NrTest)); gram = ZERO
!
   if (blocks) then
!
      select case (TEST_NORM)
         case (MATH_NORM)
            call DPOTRF(uplo,NrdofEE,gram_r,NrdofEE,info); if (info.ne.0) goto 310
            gram(   1:nE   ,   1:nE   ) = cmplx(gram_r(1:nE,1:nE),0.d0,8)
            gram(nE+1:nE+nE,nE+1:nE+nE) = cmplx(gram_r(1:nE,1:nE),0.d0,8)
         case (GRAPH_DIAG)
            call ZPOTRF(uplo,NrdofEE,gram_FF,NrdofEE,info); if (info.ne.0) goto 310
            call ZPOTRF(uplo,NrdofEE,gram_GG,NrdofEE,info); if (info.ne.0) goto 310
            gram(   1:nE   ,   1:nE   ) = gram_FF(1:nE,1:nE)
            gram(nE+1:nE+nE,nE+1:nE+nE) = gram_GG(1:nE,1:nE)
            deallocate(gram_FF,gram_GG)
         case (GRAPH_NORM)
            ! only need the upper triangle
            gram(   1:nE   ,   1:nE   ) = gram_FF(1:nE,1:nE)
            gram(   1:nE   ,nE+1:nE+nE) = gram_FG(1:nE,1:nE)
            gram(nE+1:nE+nE,nE+1:nE+nE) = gram_GG(1:nE,1:nE)
            deallocate(gram_FF,gram_GG,gram_FG)
            call ZPOTRF(uplo,NrTest,gram,NrTest,info); if (info.ne.0) goto 310
      end select
!
   else
!
!  ...A. Compute Cholesky factorization of Gram Matrix, G=U^*U (=LL^*)
      select case (TEST_NORM)
         case (MATH_NORM)
            call DPOTRF(uplo,NrdofEE,gram_r,NrdofEE,info); if (info.ne.0) goto 310
            do j = 1,NrdofEE
               joff = 2*(j-1)
               do i = 1,j
                  ioff = 2*(i-1)
                  ! (F,F), (G,G)
                  gram(ioff+1,joff+1) = cmplx(gram_r(i,j),0.d0,8)
                  gram(ioff+2,joff+2) = cmplx(gram_r(i,j),0.d0,8)
               enddo
            enddo
            deallocate(gram_r)
         case (GRAPH_DIAG)
            call ZPOTRF(uplo,NrdofEE,gram_FF,NrdofEE,info); if (info.ne.0) goto 310
            call ZPOTRF(uplo,NrdofEE,gram_GG,NrdofEE,info); if (info.ne.0) goto 310
            do j = 1,NrdofEE
               joff = 2*(j-1)
               do i = 1,j
                  ioff = 2*(i-1)
                  ! (F,F), (G,G)
                  gram(ioff+1,joff+1) = gram_FF(i,j)
                  gram(ioff+2,joff+2) = gram_GG(i,j)
               enddo
            enddo
            deallocate(gram_FF,gram_GG)
         case (GRAPH_NORM)
            do j = 1,NrdofEE
               joff = 2*(j-1)
               do i = 1,j
                  ioff = 2*(i-1)
                  ! (F,F), (G,G)
                  gram(ioff+1,joff+1) = gram_FF(i,j)
                  gram(ioff+2,joff+2) = gram_GG(i,j)
                  ! (F,G), (G,F)
                  gram(ioff+1,joff+2) = gram_FG(i,j)
                  gram(ioff+2,joff+1) = conjg(gram_FG(j,i))
               enddo
            enddo
            deallocate(gram_FF,gram_GG,gram_FG)
            call ZPOTRF(uplo,NrTest,gram,NrTest,info); if (info.ne.0) goto 310
      end select
!..end if define linear system block-wise
   endif
   310 continue
   if (info.ne.0) then
      write(*,*) 'elem_opt: POTRF: Mdle,info = ',Mdle,info,'. stop.'
      stop
   endif
!
!..B. Solve triangular system to obtain B~, (LX=) U^*X = [B|l]
   call ZTRTRS(uplo,'C','N',NrTest,NrTrial+1,gram,NrTest,stiff_ALL,NrTest,info)
   if (info.ne.0) then
      write(*,*) 'elem_opt: ZTRTRS: Mdle,info = ',Mdle,info,'. stop.'
      stop
   endif
!
   deallocate(gram)
   allocate(zaloc(NrTrial+1,NrTrial+1)); zaloc = ZERO
!
!..C. Matrix multiply: B^* G^-1 B (=B~^* B~)
   call ZHERK(uplo,'C',NrTrial+1,NrTest,ZONE,stiff_ALL,NrTest,ZERO,zaloc,NrTrial+1)
!
   deallocate(stiff_ALL)
!
!..D. Fill lower triangular part of Hermitian matrix
   do i=1,NrTrial
      zaloc(i+1:NrTrial+1,i) = conjg(zaloc(i,i+1:NrTrial+1))
   enddo
!
!..E. Fill ALOC and BLOC matrices
!  - if using the block-wise DPG system,
!    extract components as defined in the physics file
   if (blocks) then
      nE =   NrdofEi; nQ = 3*NrdofQ
      jE = 2*NrdofEi; jQ = 6*NrdofQ
      !
      do i = 1,nE
         ioff = 2*(i-1)
         ZblocE(ioff+1) = zaloc(   i,jE+jQ+1)
         ZblocE(ioff+2) = zaloc(nE+i,jE+jQ+1)
      enddo
      !
      do i = 1,nQ
         ioff = 6*((i-1)/3) + MOD(i-1,3)
         ZblocQ(ioff+1) = zaloc(jE+   i,jE+jQ+1)
         ZblocQ(ioff+4) = zaloc(jE+nQ+i,jE+jQ+1)
      enddo
      !
      do j = 1,nE
         joff = 2*(j-1)
         do i = 1,nE
            ioff = 2*(i-1)
            ZalocEE(ioff+1,joff+1) = zaloc(   i,   j)
            ZalocEE(ioff+2,joff+1) = zaloc(nE+i,   j)
            ZalocEE(ioff+1,joff+2) = zaloc(   i,nE+j)
            ZalocEE(ioff+2,joff+2) = zaloc(nE+i,nE+j)
         enddo
         do i = 1,nQ
            ioff = 6*((i-1)/3) + MOD(i-1,3)
            ZalocQE(ioff+1,joff+1) = zaloc(jE+   i,   j)
            ZalocQE(ioff+4,joff+1) = zaloc(jE+nQ+i,   j)
            ZalocQE(ioff+1,joff+2) = zaloc(jE+   i,nE+j)
            ZalocQE(ioff+4,joff+2) = zaloc(jE+nQ+i,nE+j)
         enddo
      enddo
      !
      do j = 1,nQ
         joff = 6*((j-1)/3) + MOD(j-1,3)
         do i = 1,nE
            ioff = 2*(i-1)
            ZalocEQ(ioff+1,joff+1) = zaloc(   i,jE+   j)
            ZalocEQ(ioff+2,joff+1) = zaloc(nE+i,jE+   j)
            ZalocEQ(ioff+1,joff+4) = zaloc(   i,jE+nQ+j)
            ZalocEQ(ioff+2,joff+4) = zaloc(nE+i,jE+nQ+j)
         enddo
         do i = 1,nQ
            ioff = 6*((i-1)/3) + MOD(i-1,3)
            ZalocQQ(ioff+1,joff+1) = zaloc(jE+   i,jE+   j)
            ZalocQQ(ioff+4,joff+1) = zaloc(jE+nQ+i,jE+   j)
            ZalocQQ(ioff+1,joff+4) = zaloc(jE+   i,jE+nQ+j)
            ZalocQQ(ioff+4,joff+4) = zaloc(jE+nQ+i,jE+nQ+j)
         enddo
      enddo
      !
   else
      jE = 2*NrdofEi; jQ = 6*NrdofQ
      !
      ZblocE(1:jE) = zaloc(   1:jE   ,jE+jQ+1)
      ZblocQ(1:jQ) = zaloc(jE+1:jE+jQ,jE+jQ+1)
      !
      ZalocEE(1:jE,1:jE) = zaloc(1:jE,1:jE)
      ZalocEQ(1:jE,1:jQ) = zaloc(1:jE,jE+1:jE+jQ)
      !
      ZalocQE(1:jQ,1:jE) = zaloc(jE+1:jE+jQ,1:jE)
      ZalocQQ(1:jQ,1:jQ) = zaloc(jE+1:jE+jQ,jE+1:jE+jQ)
   endif
!
   deallocate(zaloc)
!
!..TIMER
   if (timer) then
      end_time = MPI_Wtime()
      !$OMP CRITICAL
      write(*,11) 'elem DPG LinAlg: ', end_time-start_time
      !$OMP END CRITICAL
   endif
!
!-------------------------------------------------------------------------------
!       I M P E D A N C E   B O U N D A R Y
!-------------------------------------------------------------------------------
!
!..Implementation of impedance BC via L2 penalty term
   if (IBCFLAG.eq.2) call imp_penalty(Mdle,NrdofH,NrdofEi,MdE,          &
                                      norder,norderi, ZblocE,ZalocEE)
!
end subroutine elem_opt
