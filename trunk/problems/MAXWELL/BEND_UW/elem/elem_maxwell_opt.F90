!------------------------------------------------------------------------------
!> @brief      Evaluates unconstrained stiffness matrix and load vector for the
!!             envelope Maxwell UW formulation with a bending-related ansatz
!!             USING integration based on BLAS routines
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
!> @date       Oct 2024
!------------------------------------------------------------------------------
subroutine elem_bend_env_maxwell_opt(Mdle,                      &
                                     NrTest,NrTrial,            &
                                     NrdofEE,                   &
                                     NrdofH,NrdofE,NrdofQ,      &
                                     NrdofEi,                   &
                                     MdE,MdQ,                   &
                                     ZblocE,ZalocEE,ZalocEQ,    &
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
!..declare edge/face type variables
   integer :: ntype,ftype
!
!..declare element order, orientation for edges and faces
   integer :: norder(19), norient_edge(12), norient_face(6)
!
!..element nodes order (trial) for interfaces
   integer :: norderi(19), norderf(5)
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
   real(8) :: shapH(MAXbrickH),  gradH(3,MAXbrickH)
!
!..H(curl) shape functions
   real(8) :: shapE (3,MAXbrickE),  curlE(3,MAXbrickE)
   real(8) :: shapFi(3,NrdofEi)
!  ... enriched
   real(8) :: shapEE(3,MAXbrickEE), curlEE(3,MAXbrickEE)
   real(8) :: shapF (3,NrdofEE),  curlF(3,NrdofEE)
!
!..L2 shape functions
   real(8) :: shapQ(MAXbrickQ)
!
!..load vector for the enriched space
   complex(8) :: bload_E(NrTest)
!
!..Gram matrix in full format (although we'll use upper half only)
   real(8),    allocatable :: gram_r(:,:)
   complex(8), allocatable :: gram_FF(:,:),gram_FG(:,:),gram_GG(:,:)
   complex(8), allocatable :: gram(:,:)
!
!..Field and material values
   real(8) :: fldQ, afac, fldF(3), crlF(3)
   complex(8) :: AstarF1(3,NrdofEE),AstarF2(3,NrdofEE)
   complex(8) :: AstarG1(3,NrdofEE),AstarG2(3,NrdofEE)
!
!..dynamic allocation arrays for filling values at quadrature points
   real(8),    allocatable :: all_shapF(:,:),all_curlF(:,:),all_fldE(:,:)
   complex(8), allocatable :: all_zJ(:),all_zL(:),                        &
                              all_AstarF1(:,:),all_AstarF2(:,:),          &
                              all_AstarG1(:,:),all_AstarG2(:,:)
!..stiffness matrices (transposed) for the enriched test space
   complex(8), allocatable :: zload_F(:),zload_G(:)
   complex(8), allocatable :: zstiff_FE(:,:),zstiff_FH(:,:)
   complex(8), allocatable :: zstiff_GE(:,:),zstiff_GH(:,:)
   complex(8), allocatable :: stiff_EEi(:,:),stiff_EQ(:,:)
   complex(8), allocatable :: stiff_ALL(:,:),zBDPG(:,:)
!
!..quadrature data
   real(8) :: xiloc(3,MAXNINT3ADD), waloc(MAXNINT3ADD)
   real(8) :: tloc(2,MAXNINT2ADD), wtloc(MAXNINT2ADD)
   real(8) :: weight,wa,sqrt_weight
!
!..BC's flags
   integer :: ibc(6,NRINDEX)
!
!..for auxiliary computation
   complex(8) :: zaux,zcux
!
!..Maxwell load and auxiliary variables
   complex(8) :: zJ(3), zImp(3), zL(3)
   real(8), dimension(3) :: E1,E2,rntimesE,rn2timesE
!
!..number of edge,faces per element type
   integer :: nre, nrf
!
!..various variables for the problem
   integer :: i1, j1, j2, k1, k2, i, j, k, l, n, m, nint
   integer :: nda, loff, koff, joff, ioff, nE, nQ
   integer :: iflag, ifc, info, nrdof, nordP, nsign
!
!..TIMER
   real(8) :: start_time, end_time
   logical, parameter :: timer = .false.
!
   integer, external :: ij_upper_to_packed
!
#if HP3D_DEBUG
!..Set iprint = 0/1 (Non-/VERBOSE)
   integer :: iprint
   iprint = 0
   ! if (Mdle.eq.118) iprint = 1
#endif
!
!-------------------------------------------------------------------------------
!
!..TIMER
   if (timer) start_time = MPI_Wtime()
!
#if HP3D_DEBUG
   if (iprint.eq.1) then
      write(*,*) 'elem_maxwell: Mdle = ', Mdle
   endif
#endif
!
!..allocate matrices
   nE =   NrdofEE
   nQ = 3*NrdofQ
   allocate(gram(NrTest,NrTest))
   allocate(stiff_EEi(NrTest,2*NrdofEi))
   allocate(stiff_EQ(NrTest,6*NrdofQ))
!
!..element type
   ntype = NODES(Mdle)%ntype
   nre = nedge(ntype); nrf = nface(ntype)
!
!..determine order of approximation
   call find_order(Mdle, norder)
   norderi(1:nre+nrf) = norder(1:nre+nrf)
!
!..set the enriched order of approximation
   select case(ntype)
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
         write(*,*) 'elem_maxwell: invalid ntype param. stop.'
         stop
   end select
!
!..determine edge and face orientations
   norient_edge(:) = 0
   norient_face(:) = 0
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
   bload_E(:)     = ZERO
   gram(:,:)      = ZERO
   stiff_EEi(:,:) = ZERO
   stiff_EQ(:,:)  = ZERO
!
!-----------------------------------------------------------------------
!              E L E M E N T   I N T E G R A L S
!-----------------------------------------------------------------------
!
!..use the enriched order to set the quadrature
   INTEGRATION = NORD_ADD
   call set_3D_int_DPG(ntype,norder,norient_face, nint,xiloc,waloc)
   INTEGRATION = 0
!
!..we allocate the auxiliary arrays with this dimension
   nda = 3*nint
   allocate(all_fldE(nQ,nda))
   allocate(all_zJ(nda),all_zL(nda))
   allocate(all_shapF(nE,nda),all_curlF(nE,nda))
   allocate(all_AstarF1(nE,nda),all_AstarF2(nE,nda))
   allocate(all_AstarG1(nE,nda),all_AstarG2(nE,nda))
   allocate(zstiff_FE(nE,nQ),zstiff_FH(nE,nQ))
   allocate(zstiff_GE(nE,nQ),zstiff_GH(nE,nQ))
!..initialize the allocated arrays
   all_fldE(:,:)    = ZERO
   all_zJ(:)        = ZERO; all_zL(:)        = ZERO
   all_shapF(:,:)   = ZERO; all_curlF(:,:)   = ZERO
   all_AstarF1(:,:) = ZERO; all_AstarF2(:,:) = ZERO
   all_AstarG1(:,:) = ZERO; all_AstarG2(:,:) = ZERO
   zstiff_FE(:,:)   = ZERO; zstiff_FH(:,:)   = ZERO
   zstiff_GE(:,:)   = ZERO; zstiff_GH(:,:)   = ZERO
!
!..loop over integration points
   do l=1,nint
!  ...offset in auxiliary matrices
      loff = 3*(l-1)
!
      xi(1:3)=xiloc(1:3,l); wa=waloc(l)
!
!  ...H1 shape functions (for geometry)
      call shape3DH(ntype,xi,norder,norient_edge,norient_face, nrdof,shapH,gradH)
!
!  ...L2 shape functions for the trial space
      call shape3DQ(ntype,xi,norder, nrdof,shapQ)
!
!  ...broken H(curl) shape functions for the enriched test space
      call shape3EE(ntype,xi,nordP, nrdof,shapEE,curlEE)
!
!  ...geometry map
      call geom3D(Mdle,xi,xnod,shapH,gradH,NrdofH, x,dxdxi,dxidx,rjac,iflag)
!
!
#if HP3D_DEBUG
      if (iflag .ne. 0) then
         write(*,5999) Mdle,rjac
5999     format('elem_maxwell: Negative Jacobian. Mdle,rjac=',i8,2x,e12.5)
         stop
      endif
#endif
!
!  ...integration weight
      weight = rjac*wa
      sqrt_weight = sqrt(weight)
!
!  ...loop through L2 trial shape functions
      do k=1,NrdofQ
!     ...Piola transformation
         fldQ = shapQ(k)/rjac * sqrt_weight
!
!     ...save in the big matrix of size NrdofQ by 3*nint
         koff = 3*(k-1)
         all_fldE(koff+1,loff+1) = fldQ 
         all_fldE(koff+2,loff+2) = fldQ 
         all_fldE(koff+3,loff+3) = fldQ 
      enddo
!
!  ...get the RHS
      call getf(Mdle,x, zJ, zL)
!  ...save in the big vector of length 3*nint
      all_zJ(loff+1:loff+3) = sqrt_weight*zJ(1:3)
      all_zL(loff+1:loff+3) = sqrt_weight*zL(1:3)
!
!    IDEA 1: FIRST MULTIPLY BLOCKS, THEN MERGE IN CORRECT ORDER
!
!  ...apply pullbacks
      call DGEMM('T','N',3,nE,3,1.d0     ,dxidx,3,shapEE,3,0.d0,shapF,3)!
!  ...save in the big matrix of size nE by 3*nint
      all_shapF(1:nE,loff+1:loff+3) = sqrt_weight*transpose(shapF)
!
      call DGEMM('N','N',3,nE,3,1.d0/rjac,dxdxi,3,curlEE,3,0.d0,curlF,3)
!  ...save in the big matrix of size nE by 3*nint
      all_curlF(1:nE,loff+1:loff+3) = sqrt_weight*transpose(curlF)
!
!  ...Evaluate A^* on (F ; 0) and on (0 ; G)
      call get_Astar_multiple(Mdle,x,nE,shapF,curlF,AstarF1,AstarF2,AstarG1,AstarG2)
!  ...save in big array of size nE by 3*nint; apply conjugate now for convenience
      all_AstarF1(1:nE,loff+1:loff+3) = conjg(sqrt_weight*transpose(AstarF1))
      all_AstarF2(1:nE,loff+1:loff+3) = conjg(sqrt_weight*transpose(AstarF2))
      all_AstarG1(1:nE,loff+1:loff+3) = conjg(sqrt_weight*transpose(AstarG1))
      all_AstarG2(1:nE,loff+1:loff+3) = conjg(sqrt_weight*transpose(AstarG2))
!
!..end of loop through integration points
   enddo
!
!..LOAD VECTOR INTEGRATION
!
   ! write(*,*) 'STARTING LOAD VECTOR INTEGRATION'

   allocate(zload_F(nE),zload_G(nE))

   ! write(*,*) 'AFTER AFTER ALLOCATING zload_F, zload_G'
!..Compute (J,F)
   call ZGEMM('N','N',nE,1,nda,ZONE,cmplx(all_shapF,0.d0,8),nE,all_zJ,nda,ZERO,zload_F,nE)
!..Compute (L,G)
   call ZGEMM('N','N',nE,1,nda,ZONE,cmplx(all_shapF,0.d0,8),nE,all_zL,nda,ZERO,zload_G,nE)
!..organize the full load vector
   do j=1,nE
      bload_E(2*j-1) = zload_F(j)
      bload_E( 2*j ) = zload_G(j)
   enddo
!..free memory that will not be used again
   deallocate(all_zJ,all_zL,zload_F,zload_G)

   ! write(*,*) 'FINISHING LOAD VECTOR INTEGRATION'
   ! write(*,*) 'OPT bload_E = ', bload_E
!
!..STIFFNESS MATRIX INTEGRATION
!
!..Compute ( E , AstarF1 )
   call ZGEMM('N','T',nE,nQ,nda,ZONE,all_AstarF1,nE,cmplx(all_fldE,0.d0,8),nQ,ZERO,zstiff_FE,nE)
!..Compute ( E , AstarG1 )
   call ZGEMM('N','T',nE,nQ,nda,ZONE,all_AstarG1,nE,cmplx(all_fldE,0.d0,8),nQ,ZERO,zstiff_GE,nE)
!..Compute ( H , AstarF2 )
   call ZGEMM('N','T',nE,nQ,nda,ZONE,all_AstarF2,nE,cmplx(all_fldE,0.d0,8),nQ,ZERO,zstiff_FH,nE)
!..Compute ( H , AstarG2 )
   call ZGEMM('N','T',nE,nQ,nda,ZONE,all_AstarG2,nE,cmplx(all_fldE,0.d0,8),nQ,ZERO,zstiff_GH,nE)
!..organize the full stiffness matrix
   do k=1,NrdofQ
      koff = (k-1)*3
      do j=1,nE
         joff = (j-1)*2      
         stiff_EQ(joff+1,2*koff+1:2*koff+3) = zstiff_FE(j,koff+1:koff+3)
         stiff_EQ(joff+2,2*koff+1:2*koff+3) = zstiff_GE(j,koff+1:koff+3)
         stiff_EQ(joff+1,2*koff+4:2*koff+6) = zstiff_FH(j,koff+1:koff+3)
         stiff_EQ(joff+2,2*koff+4:2*koff+6) = zstiff_GH(j,koff+1:koff+3)
      enddo
   enddo
!
!  --- Gram matrix: assemble and factorize ---
!
   allocate(gram_r(nE,nE))
   gram_r(:,:) = 0.d0
!
   select case(TEST_NORM)
      case( MATH_NORM)
         afac = 1.d0
!
!     ...α*(F,F)
!     ...α*(G,G)
         call DSYRK('U','N',nE,nda,afac,all_shapF,nE,0.d0,gram_r,nE) ! real-valued
!     ...(curl F, curl F)
!     ...(curl G, curl G)
         call DSYRK('U','N',nE,nda,1.d0,all_curlF,nE,1.d0,gram_r,nE) ! real-valued
!     ...Cholesky factorization of the real symmetric matrix Gram_r 
         call DPOTRF('U',nE,gram_r,nE,info)
         if (info.ne.0) then
            write(*,*) 'elem_maxwell_opt: DPOTRF: Mdle,info = ',Mdle,info,'. stop.'
            stop
         endif
         do j = 1,nE
            joff = 2*(j-1)
            do i = 1,j
               ioff = 2*(i-1)
!           ...Assemble            
!              (F,F) + (CF,CF)         0 
!                    0           (G,G) + (CG,CG)
               gram(ioff+1,joff+1) = cmplx(gram_r(i,j),0.d0,8)
               gram(ioff+2,joff+2) = cmplx(gram_r(i,j),0.d0,8)
            enddo
         enddo
         deallocate(gram_r,all_shapF,all_curlF)
   
      case(GRAPH_NORM)
         afac = ALPHA_NORM
!  
!        ...α*(F,F)
!        ...α*(G,G)
         call DSYRK('U','N',nE,nda,afac,all_shapF,nE,0.d0,gram_r,nE) ! real-valued
         do j = 1,nE
            joff = 2*(j-1)
            do i = 1,j
               ioff = 2*(i-1)
!           ...Accumulate into the complex-valued matrix
!              α*(F,F)     0 
!                 0     α*(G,G)
               gram(ioff+1,joff+1) = cmplx(gram_r(i,j),0.d0,8)
               gram(ioff+2,joff+2) = cmplx(gram_r(i,j),0.d0,8)
            enddo
         enddo
         deallocate(gram_r,all_shapF,all_curlF)
!         
         allocate(gram_FF(nE,nE),gram_GG(nE,nE),gram_FG(nE,nE))
         gram_FF(:,:) = ZERO; gram_GG(:,:) = ZERO; gram_FG(:,:) = ZERO
!  
!     ...accumulate (AstarF1,AstarF1)
         call ZHERK('U','N',nE,nda,ZONE,all_AstarF1,nE,  &
                    ZERO,gram_FF(1:nE,1:nE ),nE)
!     ...accumulate (AstarF2,AstarF2)
         call ZHERK('U','N',nE,nda,ZONE,all_AstarF2,nE,  &
                    ZONE,gram_FF(1:nE,1:nE ),nE)
!  
!     ...accumulate (AstarG1,AstarG1)
         call ZHERK('U','N',nE,nda,ZONE,all_AstarG1,nE,  &
                    ZERO,gram_GG(1:nE,1:nE ),nE)
!     ...accumulate (AstarG2,AstarG2)
         call ZHERK('U','N',nE,nda,ZONE,all_AstarG2,nE,  &
                    ZONE,gram_GG(1:nE,1:nE ),nE)
!  
!     ...accumulate (AstarF1,AstarG1)
         call ZGEMM('N','C',nE,nE,nda,ZONE,all_AstarF1,nE,all_AstarG1,nE,  &
                    ZERO,gram_FG(1:nE,1:nE ),nE)
!     ...accumulate (AstarF2,AstarG2)
         call ZGEMM('N','C',nE,nE,nda,ZONE,all_AstarF2,nE,all_AstarG2,nE,  &
                    ZONE,gram_FG(1:nE,1:nE ),nE)
!  
!     ...organize blocks within upper part of gram
         do j = 1,nE
            joff = 2*(j-1)
            do i = 1,j
               ioff = 2*(i-1)
               ! (F,F), (G,G)
               gram(ioff+1,joff+1) = gram(ioff+1,joff+1) + gram_FF(i,j)
               gram(ioff+2,joff+2) = gram(ioff+2,joff+2) + gram_GG(i,j)
               ! (F,G), (G,F)
               gram(ioff+1,joff+2) = gram_FG(i,j)
               gram(ioff+2,joff+1) = conjg(gram_FG(j,i))
            enddo
         enddo
         deallocate(gram_FF,gram_GG,gram_FG)


         if (Mdle.eq.133) then
            write(*,*) 'elem_maxwell_opt: Gram for Mdle=',Mdle
            do i=1,10
               write(*,6001) (gram(i,l),l=i,10)
            enddo
            6001 format(20e13.5)
         endif
         

         call ZPOTRF('U',NrTest,gram,NrTest,info)
         if (info.ne.0) then
            write(*,*) 'elem_maxwell_opt: ZPOTRF: Mdle,info = ',Mdle,info,'. stop.'
            stop
         endif
      case default
         write(*,*) 'elem_maxwell: invalid test norm!'
         stop
   end select

   deallocate(all_AstarF1,all_AstarF2)
   deallocate(all_AstarG1,all_AstarG2)

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
!..loop through element faces
   do ifc=1,nrf
!
!  ...sign factor to determine the OUTWARD normal unit vector
      nsign = nsign_param(ntype,ifc)
!
!  ...face type
      ftype = face_type(ntype,ifc)
!
!  ...face order of approximation
      call face_order(ntype,ifc,norder, norderf)
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
         call face_param(ntype,ifc,t, xi,dxidt)
!
!     ...determine discontinuous Hcurl shape functions
         call shape3EE(ntype,xi,nordP, nrdof,shapEE,curlEE)
#if HP3D_DEBUG
         if (nrdof .ne. NrdofEE) then
            write(*,*) 'elem_maxwell: INCONSISTENCY NrdofEE. stop.'
            stop
         endif
#endif
!
!     ...determine element H1 shape functions (for geometry)
         call shape3DH(ntype,xi,norder,norient_edge,norient_face, &
                       nrdof,shapH,gradH)
#if HP3D_DEBUG
      if (nrdof .ne. NrdofH) then
         write(*,*) 'elem_maxwell: INCONSISTENCY NrdofH. stop.'
         stop
      endif
#endif
!
!     ...determine element H(curl) shape functions (for fluxes)
!     ...for interfaces only (no bubbles)
         call shape3DE(ntype,xi,norderi,norient_edge,norient_face, &
                       nrdof,shapE,curlE)
#if HP3D_DEBUG
      if (nrdof .ne. NrdofEi) then
         write(*,*) 'elem_maxwell: INCONSISTENCY NrdofEi. stop.'
         stop
      endif
#endif
!
!     ...geometry
         call bgeom3D(Mdle,xi,xnod,shapH,gradH,NrdofH,dxidt,nsign, &
                      x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
         weight = bjac*wtloc(l)
!
!     ...pullback trial and test functions
         call DGEMM('T','N',3,NrdofEE,3,1.d0,dxidx,3,shapEE,3,0.d0,shapF ,3)
         call DGEMM('T','N',3,NrdofEi,3,1.d0,dxidx,3,shapE ,3,0.d0,shapFi,3)
!
!        COMPUTE IMPEDANCE LOAD (elimination strategy)
         if (ibc(ifc,2).eq.3) then
!        ...loop through enriched H(curl) test functions
            do k1=1,NrdofEE
               E1(1:3) = shapF(:,k1)
!           ...check for impedance BC (elimination strategy)
!              (impedance constant is GAMMA for TE10 mode in rectangular waveguide)
!              ( < n x H , F > = GAMMA*< n x n x E , F > + < zg , F > )               
!           ...get the boundary source [zImp should be zero here]
               call get_bdSource(Mdle,x,rn, zImp)
!           ...accumulate for the load vector
               k = 2*k1-1
               bload_E(k) = bload_E(k) &
                          - (zImp(1)*E1(1)+zImp(2)*E1(2)+zImp(3)*E1(3))*weight
            enddo
!     ...end if for impedance BC
         endif
!
!        COMPUTE STIFFNESS CONTRIBUTIONS
!     ...loop through H(curl) trial functions
         do k2=1,NrdofEi
            E2(1:3) = shapFi(:,k2)
            call cross_product(rn,E2, rntimesE)
!        ...loop through enriched H(curl) test functions
            do k1=1,NrdofEE
               E1(1:3) = shapF(:,k1)
!
!           ...check for impedance BC (elimination strategy)
               if (ibc(ifc,2).eq.3) then
!           ...accumulate for the extended stiffness matrix on IBC
                  call cross_product(rn,rntimesE, rn2timesE)
                  stiff_EEi(2*k1-1,2*k2-1) = stiff_EEi(2*k1-1,2*k2-1) &
                                           + (  E1(1)*rn2timesE(1) &
                                              + E1(2)*rn2timesE(2) &
                                              + E1(3)*rn2timesE(3) &
                                             )*GAMMA*weight
               else
!           ...accumulate for the extended stiffness matrix without IBC
                  stiff_EEi(2*k1-1,2*k2) = stiff_EEi(2*k1-1,2*k2) &
                                         + (  E1(1)*rntimesE(1) &
                                            + E1(2)*rntimesE(2) &
                                            + E1(3)*rntimesE(3) &
                                               )*weight
!           ...end if for impedance BC
               endif
               stiff_EEi(2*k1,2*k2-1) = stiff_EEi(2*k1,2*k2-1) &
                                      + (  E1(1)*rntimesE(1) &
                                         + E1(2)*rntimesE(2) &
                                         + E1(3)*rntimesE(3) &
                                            )*weight
!        ...end loop through H(curl) trial functions
            enddo
!     ...end loop through the enriched H(curl) test functions
         enddo
!  ...end loop through integration points
      enddo
!..end loop through faces
   enddo
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
   allocate(stiff_ALL(NrTest,NrTrial+1))
!
!..Total test/trial DOFs of the element
   i1 = NrTest ; j1 = 2*NrdofEi ; j2 = 6*NrdofQ
!
   stiff_ALL(1:i1,1:j1)       = stiff_EEi(1:i1,1:j1)
   stiff_ALL(1:i1,j1+1:j1+j2) = stiff_EQ(1:i1,1:j2)
   stiff_ALL(1:i1,j1+j2+1)    = bload_E(1:i1)
!
   deallocate(stiff_EEi,stiff_EQ)

#if HP3D_DEBUG
   if (iprint.eq.1) then
      do i = 1,NrdofEE
         write(*,*)   'i=',i
         k = 2*(i-1)
        do m = 1,NrdofQ
           l = j1 + 6*(m-1)
           ! if (i.eq.1 .and. m.eq.1) then
              write(*,124) stiff_ALL(k+1,l+1:l+6)
              write(*,124) stiff_ALL(k+2,l+1:l+6)
              124 format(6f14.9,6f14.9)
           ! endif
        enddo
      enddo
      write(*,*) 'bload_E='
      write(*,124) bload_E(1:i1)
      call pause
   endif
#endif
!
!
!..B. Solve triangular system to obtain B~, (LX=) U^*X = [B|l]
   call ZTRTRS('U','C','N',NrTest,NrTrial+1,gram,NrTest,stiff_ALL,NrTest,info)
   if (info.ne.0) then
      write(*,*) 'elem_maxwell_opt: ZTPTRS: Mdle,info = ',Mdle,info,'. stop.'
      stop
   endif
!
   deallocate(gram)
!
   allocate(zBDPG(NrTrial+1,NrTrial+1)); zBDPG = ZERO
!
!..C. Matrix multiply: B^* G^-1 B (=B~^* B~)
   call ZHERK('U','C',NrTrial+1,NrTest,ZONE,stiff_ALL,NrTest,ZERO,zBDPG,NrTrial+1)
!
   deallocate(stiff_ALL)
!
!..D. Fill lower triangular part of Hermitian matrix
   do i=1,NrTrial
      zBDPG(i+1:NrTrial+1,i) = conjg(zBDPG(i,i+1:NrTrial+1))
   enddo
!
!..E. Fill ALOC and BLOC matrices
   ZblocE(1:j1) = zBDPG(1:j1,j1+j2+1)
   ZblocQ(1:j2) = zBDPG(j1+1:j1+j2,j1+j2+1)
!
   ZalocEE(1:j1,1:j1) = zBDPG(1:j1,1:j1)
   ZalocEQ(1:j1,1:j2) = zBDPG(1:j1,j1+1:j1+j2)
!
   ZalocQE(1:j2,1:j1) = zBDPG(j1+1:j1+j2,1:j1)
   ZalocQQ(1:j2,1:j2) = zBDPG(j1+1:j1+j2,j1+1:j1+j2)
!
   deallocate(zBDPG)
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
end subroutine elem_bend_env_maxwell_opt
