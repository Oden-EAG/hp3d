!-------------------------------------------------------------------------
!
!     routine name      - elem_opt
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
!              NrTest   - total number of test dof
!              NrTrial  - total number of trial dof
!              NrdofHH  - number of H1 test dof
!              NrdofH   - number of H1 trial dof
!              NrdofV   - number of H(div) trial dof
!              NrdofVi  - number of H(div) trial interface dof
!              MdH      - num rows of AlocHH,AlocHV
!              MdV      - num rows of AlocVH,AlocVV
!        out:
!              BlocH    - load vectors
!              BlocV
!              AlocHH   - stiffness matrices
!              AlocHV
!              AlocVH
!              AlocVV
!
!-------------------------------------------------------------------------
!
subroutine elem_opt(Mdle,                   &
                    NrTest,NrTrial,         &
                    NrdofHH,                &
                    NrdofH,NrdofV,          &
                    NrdofVi,                &
                    MdH,MdV,                &
                    BlocH,AlocHH,AlocHV,    &
                    BlocV,AlocVH,AlocVV)
!..ALOC: holds local element stiffness matrices
!..BLOC: holds local element load vectors
   use control
   use data_structure3D
   use element_data
   use parametersDPG
   use mpi_wrapper
!
   implicit none
!
!..declare input/output variables
   integer,                       intent(in)  :: Mdle
   integer,                       intent(in)  :: NrTest
   integer,                       intent(in)  :: NrTrial
   integer,                       intent(in)  :: NrdofHH
   integer,                       intent(in)  :: NrdofH
   integer,                       intent(in)  :: NrdofV
   integer,                       intent(in)  :: NrdofVi
   integer,                       intent(in)  :: MdH
   integer,                       intent(in)  :: MdV
   real(8), dimension(MdH),       intent(out) :: BlocH
   real(8), dimension(MdH,MdH),   intent(out) :: AlocHH
   real(8), dimension(MdH,MdV),   intent(out) :: AlocHV
   real(8), dimension(MdV),       intent(out) :: BlocV
   real(8), dimension(MdV,MdH),   intent(out) :: AlocVH
   real(8), dimension(MdV,MdV),   intent(out) :: AlocVV
!
!-------------------------------------------------------------------------
!
!..aux variables
   real(8) :: rjac, bjac, fval, wa, weight, sqrt_weight
   integer :: iflag, nrv, nre, nrf, nint
   integer :: j1, j2, i, k, l
   integer :: nordP, nrdof, nsign, ifc, info
!
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
!
!..H(div) shape functions
   real(8) :: shapV(3,MAXbrickV), divV(MAXbrickV)
!
!..Enriched H1 shape functions
   real(8) :: shapHH(MAXbrickHH), gradHH(3,MAXbrickHH)
!
!..load vector for the enriched space
   real(8) :: bload_H(NrTest)
!
!..Gram matrix
   real(8), allocatable :: gram_RFP(:)
!
!..stiffness matrices for the enriched test space
   real(8), allocatable :: stiff_HV(:,:),stiff_ALL(:,:),raloc(:,:)
!
!..auxiliary matrices (storing info at integration points)
   real(8), allocatable :: test_H(:,:), test_GH(:,:), trial_GH(:,:), trial_V(:,:)
   integer :: noff, nda
!
!..3D quadrature data
   real(8) :: xiloc(3,MAXNINT3ADD), waloc(MAXNINT3ADD)
!
!..2D quadrature data
   real(8) :: tloc(2,MAXNINT2ADD), wtloc(MAXNINT2ADD)
!
!..workspace for trial and test variables
   real(8) :: v, sn
   real(8) :: dp(3), dv(3), s(3)
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
!..clear space for stiffness matrix and load vector
   BlocH = ZERO; AlocHH = ZERO; AlocHV = ZERO
   BlocV = ZERO; AlocVH = ZERO; AlocVV = ZERO
!
!..clear space for auxiliary vector
   bload_H = ZERO
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
   allocate(test_H  (NrdofHH,nint))
   allocate(test_GH (NrdofHH,nda))
   allocate(trial_GH(NrdofH ,nda))
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
!  ...discontinuous H1 shape functions
      call shape3HH(etype,xi,nordP, nrdof,shapHH,gradHH)
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
!
!     ...accumulate load vector: (f,v)
         bload_H(k) = bload_H(k) + fval*v*weight
!
!     ...fill auxiliary matrices
         test_H(k,l) = v * sqrt_weight

         dv(1:3) = gradHH(1,k)*dxidx(1,1:3) &
                 + gradHH(2,k)*dxidx(2,1:3) &
                 + gradHH(3,k)*dxidx(3,1:3)
         test_GH(k,noff+1:noff+3) = dv(1:3) * sqrt_weight
!
!  ...end of loop through enriched H1 test functions
      enddo
!
!  ...loop through H1 trial functions
      do k=1,NrdofH
!     ...Piola transformation
         dp(1:3) = gradH(1,k)*dxidx(1,1:3) &
                 + gradH(2,k)*dxidx(2,1:3) &
                 + gradH(3,k)*dxidx(3,1:3)
!
!     ...fill auxiliary matrix
         trial_GH(k,noff+1:noff+3) = dp(1:3) * sqrt_weight
!
!  ...end of loop through H1 trial functions
      enddo
!
!  ...alternative computation (apply pull-backs to all shape functions via DGEMM)
!      dxidx = dxidx * sqrt_weight
!      call piola_hcurl_trans(dxidx,NrdofHH,gradHH,  test_GH(1:NrdofHH,noff+1:noff+3))
!      call piola_hcurl_trans(dxidx,NrdofH ,gradH , trial_GH(1:NrdofH ,noff+1:noff+3))
!
!..end of loop through integration points
   enddo
!
!..compute Gram matrix using rectangular full packed (RFP) format
   allocate(gram_RFP(NrTest*(NrTest+1)/2))
   call DSFRK('N',uplo,'N',NrdofHH,nda ,1.0d0,test_GH(1:NrdofHH,1:nda),NrdofHH,  &
                                        0.0d0,gram_RFP)
   call DSFRK('N',uplo,'N',NrdofHH,nint,1.0d0,test_H(1:NrdofHH,1:nint),NrdofHH,  &
                                        1.0d0,gram_RFP)
!
!..compute stiffness matrix
   allocate(stiff_ALL(NrTest,NrTrial+1))
   call DGEMM('N','T',NrdofHH,NrdofH,nda,1.0d0,test_GH (1:NrdofHH,1:nda)    ,NrdofHH,  &
                                               trial_GH(1:NrdofH ,1:nda)    ,NrdofH ,  &
                                         0.0d0,stiff_ALL(1:NrdofHH,1:NrdofH),NrdofHH)
!
   deallocate(test_H,test_GH,trial_GH)
!
!..TIMER
   if (timer) then
      end_time = MPI_Wtime()
      !$OMP CRITICAL
      write(*,11) 'elem INTEGR Vol: ', end_time-start_time
      !$OMP END CRITICAL
      11 format(A,f12.7,' s')
      start_time = MPI_Wtime()
   endif
!
!---------------------------------------------------------------------
!    B O U N D A R Y    I N T E G R A L S                            |
!---------------------------------------------------------------------
!
   allocate(stiff_HV(NrTest,NrdofVi)); stiff_HV = ZERO
   allocate( test_H(NrdofHH,MAXNINT2ADD))
   allocate(trial_V(NrdofVi,MAXNINT2ADD))
!
   noff = 0
!
!..loop through element faces
   do ifc=1,nrf
!
      trial_V = ZERO
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
!     ...determine discontinuous H1 shape functions
         call shape3HH(etype,xi,nordP, nrdof,shapHH,gradHH)
!
!     ...determine element H(div) shape functions (for fluxes)
!        for interfaces only (no bubbles)
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
!     ...geometry
         call bgeom3D(Mdle,xi,xnod,shapH,gradH,NrdofH,dxidt,nsign, &
                      x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
         weight = bjac*wtloc(l)
         sqrt_weight = sqrt(weight)
!
!     ...number of unused shape functions from other faces
!        that we need to skip (one per previous face)
         nrdof = ifc-1
!
!     ...loop through enriched H1 test functions
         test_H(1:NrdofHH,l) = shapHH(1:NrdofHH) * sqrt_weight
!
!     ...loop through H(div) trial functions
         do i=1,ndofV_face
            k = nrdof + i
!        ...Piola transformation
            s(1:3) = dxdxi(1:3,1)*shapV(1,k) &
                   + dxdxi(1:3,2)*shapV(2,k) &
                   + dxdxi(1:3,3)*shapV(3,k)
            s(1:3) = s(1:3)/rjac
!        ...normal component
            sn = s(1)*rn(1)+s(2)*rn(2)+s(3)*rn(3)
!
            k = noff + i
            trial_V(k,l) = sn * sqrt_weight
!     ...end loop through H(div) trial functions
         enddo
!  ...end loop through integration points
      enddo
!
!  ...EXTENDED HV STIFFNESS MATRIX: -<dot(σ,n), v>_{Γ_h}
      call DGEMM('N','T',NrdofHH,NrdofVi,nint,-1.0d0, test_H (1:NrdofHH,1:nint   ),NrdofHH,  &
                                                     trial_V (1:NrdofVi,1:nint   ),NrdofVi,  &
                                               1.0d0,stiff_HV(1:NrdofHH,1:NrdofVi),NrdofHH)
!
!  ...increment DOF offset for the next face
      noff = noff + ndofV_face
!
!..end loop through element faces
   enddo
!
   deallocate(test_H,trial_V)
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
!..Total test/trial DOFs of the element
   i = NrTest ; j1 = NrdofH ; j2 = NrdofVi
!
   stiff_ALL(1:i,j1+1:j1+j2) = stiff_HV(1:i,1:j2)
   stiff_ALL(1:i,j1+j2+1)    = bload_H(1:i)
!
   deallocate(stiff_HV)
!
!..A. Compute Cholesky factorization of Gram Matrix, G=U^T U (=LL^T)
   call DPFTRF('N',uplo,NrTest,gram_RFP,info)
   if (info.ne.0) then
      write(*,*) 'elem_opt: DPFTRF: Mdle,info = ',Mdle,info,'. stop.'
      stop
   endif
!
!..B. Solve triangular system to obtain B~, (LX=) U^T X = [B|l]
   call DTFSM('N','L',uplo,'T','N',NrTest,NrTrial+1,1.d0,gram_RFP,stiff_ALL,NrTest)
!
   deallocate(gram_RFP)
   allocate(raloc(NrTrial+1,NrTrial+1))
!
!..C. Matrix multiply: B^T G^-1 B (=B~^T B~)
   call DSYRK(uplo,'T',NrTrial+1,NrTest,ZONE,stiff_ALL,NrTest,ZERO,raloc,NrTrial+1)
!
   deallocate(stiff_ALL)
!
!..D. Fill lower triangular part of symmetric matrix
   do i=1,NrTrial
      raloc(i+1:NrTrial+1,i) = raloc(i,i+1:NrTrial+1)
   enddo
!
!..E. Fill ALOC and BLOC matrices
   BlocH(1:j1) = raloc(1:j1,j1+j2+1)
   BlocV(1:j2) = raloc(j1+1:j1+j2,j1+j2+1)
!
   AlocHH(1:j1,1:j1) = raloc(1:j1,1:j1)
   AlocHV(1:j1,1:j2) = raloc(1:j1,j1+1:j1+j2)
!
   AlocVH(1:j2,1:j1) = raloc(j1+1:j1+j2,1:j1)
   AlocVV(1:j2,1:j2) = raloc(j1+1:j1+j2,j1+1:j1+j2)
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
end subroutine elem_opt
