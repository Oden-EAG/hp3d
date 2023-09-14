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
   integer :: nrdofVi,nrTest,nrTrial
   integer :: ndofHmdl,ndofEmdl,ndofVmdl,ndofQmdl
   integer :: nrv,nre,nrf
!
!..element type
   integer :: etype
!
!-------------------------------------------------------------------------
!
!..activate two physics variables (H1 field + H(div) trace) for assembly
   Itest(1:2) = 1; Itrial(1:2) = 1
!
   norder (1:19) = 0
   norderP(1:19) = 0
!
   etype = NODES(Mdle)%ntype
   nrv = NVERT(etype); nre = NEDGE(etype); nrf = NFACE(etype)
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
!
!..calculate number of H(div) interface DOFs
   nrdofVi = nrdofV-ndofVmdl
!..calculate total number of trial and test DOFs
   nrTest  = nrdofHH
   nrTrial = nrdofH + nrdofVi
!
!..call element integration routine
   call elem_poisson(Mdle,nrTest,nrTrial,                               &
            nrdofHH,nrdofH,nrdofV,nrdofVi,                              &
            BLOC(1)%nrow,BLOC(2)%nrow,                                  &
            BLOC(1)%array,ALOC(1,1)%array,ALOC(1,2)%array,              &
            BLOC(2)%array,ALOC(2,1)%array,ALOC(2,2)%array)
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
subroutine elem_poisson(Mdle,                   &
                        NrTest,NrTrial,         &
                        NrdofHH,                &
                        NrdofH,NrdofV,          &
                        NrdofVi,                &
                        MdH,MdV,                &
                        BlocH,AlocHH,AlocHV,    &
                        BlocV,AlocVH,AlocVV)
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
   real(8) :: rjac, bjac, fval, wa, weight
   integer :: iflag, nrv, nre, nrf, nint
   integer :: j1, j2, k1, k2, i, k, l
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
!
!..Enriched H1 shape functions
   real(8) :: shapHH(MAXbrickHH), gradHH(3,MAXbrickHH)
!
!..load vector for the enriched space
   real(8) :: bload_H(NrTest)
!
!..Gram matrix in packed format
   real(8), allocatable :: gramP(:)
!
!..stiffness matrices for the enriched test space
   real(8), allocatable :: stiff_HH(:,:),stiff_HV(:,:),stiff_ALL(:,:),raloc(:,:)
!
!..3D quadrature data
   real(8) :: xiloc(3,MAXNINT3ADD), waloc(MAXNINT3ADD)
!
!..2D quadrature data
   real(8) :: tloc(2,MAXNINT2ADD), wtloc(MAXNINT2ADD)
!
!..workspace for trial and test variables
   real(8) :: q, v, sn, aux
   real(8) :: dq(3), dp(3), dv(3), s(3)
!
!..for Gram matrix compressed storage format
   integer :: nk
   nk(k1,k2) = (k2-1)*k2/2+k1
!
!---------------------------------------------------------------------
!--------- INITIALIZE THE ELEMENT ORDER, ORIENTATION ETC. ------------
!---------------------------------------------------------------------
!
!
!..allocate auxiliary matrices
   allocate(gramP(NrTest*(NrTest+1)/2))
   allocate(stiff_HH(NrTest,NrdofH))
   allocate(stiff_HV(NrTest,NrdofVi))
!
!..element type
   etype = NODES(Mdle)%ntype
   nrv = NVERT(etype)
   nre = NEDGE(etype)
   nrf = NFACE(etype)
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
   BlocH = ZERO; AlocHH = ZERO; AlocHV = ZERO
   BlocV = ZERO; AlocVH = ZERO; AlocVV = ZERO
!
!..clear space for auxiliary matrices
   bload_H   = ZERO
   gramP     = ZERO
   stiff_HH  = ZERO
   stiff_HV  = ZERO
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
      call shape3DH(etype,xi,norder,norient_edge,norient_face, nrdofH,shapH,gradH)
!
!  ...discontinuous H1 shape functions
      call shape3HH(etype,xi,nordP, nrdof,shapHH,gradHH)
!
!  ...geometry map
      call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH, x,dxdxi,dxidx,rjac,iflag)
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
         bload_H(k1) = bload_H(k1) + fval*v*weight
!
!     ...loop through H1 trial functions
         do k2=1,nrdofH
!        ...Piola transformation
            dp(1:3) = gradH(1,k2)*dxidx(1,1:3) &
                    + gradH(2,k2)*dxidx(2,1:3) &
                    + gradH(3,k2)*dxidx(3,1:3)
!
!        ...EXTENDED HH STIFFNESS MATRIX: (grad u, grad_h v)
            stiff_HH(k1,k2) = stiff_HH(k1,k2) + weight * ( &
                           dv(1)*dp(1) + dv(2)*dp(2) + dv(3)*dp(3))
!
!     ...end of loop through H1 trial functions
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
!
!        ...H1 test inner product: (q,v) + (grad_h q, grad_h v)
            aux = q*v + (dq(1)*dv(1) + dq(2)*dv(2) + dq(3)*dv(3))
            gramP(k) = gramP(k) + aux*weight
!
!     ...enddo 2nd loop through enriched H1 test functions
         enddo
!
!  ...end of 1st loop through enriched H1 test functions
      enddo
!..end of loop through integration points
   enddo
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
!
!     ...determine element H1 shape functions (for geometry)
         call shape3DH(etype,xi,norder,norient_edge,norient_face, &
                       nrdof,shapH,gradH)
!
!     ...determine element H(div) shape functions (for fluxes)
!     ...for interfaces only (no bubbles)
         call shape3DV(etype,xi,norderi,norient_face, &
                       nrdof,shapV,divV)
!
!     ...geometry
         call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign, &
                  x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
         weight = bjac*wtloc(l)
!
!     ...loop through enriched H1 test functions
         do k1=1,NrdofHH
            v = shapHH(k1)
!
!        ...loop through H(div) trial functions
            do k2=1,NrdofVi
!           ...Piola transformation
               s(1:3) = dxdxi(1:3,1)*shapV(1,k2) &
                        + dxdxi(1:3,2)*shapV(2,k2) &
                        + dxdxi(1:3,3)*shapV(3,k2)
               s(1:3) = s(1:3)/rjac
!           ...normal component
               sn = s(1)*rn(1)+s(2)*rn(2)+s(3)*rn(3)
!
!           ...EXTENDED HV STIFFNESS MATRIX: -<dot(σ,n), v>_{Γ_h}
               stiff_HV(k1,k2) = stiff_HV(k1,k2) - sn*v*weight
!        ...end loop through H(div) trial functions
            enddo
!     ...end loop through enriched H1 test functions
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
!
!..Total test/trial DOFs of the element
   i = NrTest ; j1 = NrdofH ; j2 = NrdofVi
!
   stiff_ALL(1:i,1:j1)       = stiff_HH(1:i,1:j1)
   stiff_ALL(1:i,j1+1:j1+j2) = stiff_HV(1:i,1:j2)
   stiff_ALL(1:i,j1+j2+1)    = bload_H(1:i)
!
   deallocate(stiff_HH,stiff_HV)
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
!
!..C. Matrix multiply: B^T G^-1 B (=B~^T B~)
   call DSYRK('U','T',NrTrial+1,NrTest,ZONE,stiff_ALL,NrTest,ZERO,raloc,NrTrial+1)
!
   deallocate(stiff_ALL)
!
!..D. Fill lower triangular part of Hermitian matrix
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
end subroutine elem_poisson

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

