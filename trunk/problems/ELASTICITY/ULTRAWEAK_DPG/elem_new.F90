!--------------------------------------------------------------------------
!> @brief      Return stiffness matrix and load vector for element
!!
!> @param[in]  Mdle   - an element (middle node) number
!> @param[out] Nrdof  - number of dof for a single component
!> @param[out] Itest  - index for assembly
!> @param[out] Itrial - index for assembly
!!
!> @date       July 2023
!--------------------------------------------------------------------------
   subroutine elem(Mdle, Itest,Itrial)
!
      use physics   , only : NR_PHYSA
      use assembly  , only : ALOC,BLOC,NR_RHS
      use data_structure3D
!
      implicit none
!
      integer, intent(in)  :: Mdle
      integer, intent(out) :: Itest(NR_PHYSA), Itrial(NR_PHYSA)
!
!  ...local variables
      integer :: ntype
      integer :: nrv, nre, nrf
      integer :: nord_add_local, nordP
      integer :: norder(19), norderP(19)
      integer :: nrdofH, nrdofE, nrdofV, nrdofQ
      integer :: nrdofHH, nrdofEE, nrdofVV, nrdofQQ
      integer :: ndofHmdl, ndofEmdl, ndofVmdl, ndofQmdl
      integer :: nrdofHi, nrdofVi
      integer :: nrTest, nrTrial
!
!--------------------------------------------------------------------------
!
!  ...get element type and number of vertices, edges, and faces
      ntype = NODES(Mdle)%ntype
      nrv = nvert(ntype); nre = nedge(ntype); nrf = nface(ntype)
!
      Itest (1:NR_PHYSA) = 0
      Itrial(1:NR_PHYSA) = 0
!
      select case(NODES(Mdle)%case)
      case(15)
         Itest(1:NR_PHYSA)=1; Itrial(1:NR_PHYSA)=1
         call elem_DPG_UWEAK_SYMM(Mdle)
!
!     ...get DPG enrichment order
         nord_add_local = NORD_ADD
!
!     ...determine order of approximation
         call find_order(Mdle, norder)
!
 10      continue
!     ...set the enriched order of approximation
         select case(ntype)
            case(MDLB);      nordP = NODES(Mdle)%order + nord_add_local*111
            case(MDLN,MDLD); nordP = NODES(Mdle)%order + nord_add_local*1
            case(MDLP);      nordP = NODES(Mdle)%order + nord_add_local*11
         end select
!
!     ...get enriched order vector
         call compute_enriched_order(ntype,nordP, norderP)
!     ...number of element trial DOFs
         call celndof(ntype,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
!     ...number of element test DOFs
         call celndof(ntype,norderP, nrdofHH,nrdofEE,nrdofVV,nrdofQQ)
!     ...number of element bubble DOFs of each type
         call ndof_nod(ntype,norder(nre+nrf+1), ndofHmdl,ndofEmdl,ndofVmdl,ndofQmdl)
!
!     ...number of test functions
         nrTest = 9*nrdofHH
!
!     ...number of trial functions of each type (H,V,Q)
         nrdofHi = nrdofH - ndofHmdl
         nrdofVi = nrdofV - ndofVmdl
         nrTrial = 3*nrdofHi + 3*nrdofVi + 9*nrdofQ
!
         if (nrTest .le. nrTrial-ndofHmdl-ndofVmdl) then
            nord_add_local = nord_add_local + 1
!$omp critical
            write(*,*) 'elem: WARNING mdle =  ', mdle
            write(*,*) 'elem: NrTest, NrTrial = ', nrTest, nrTrial-ndofHmdl-ndofVmdl
            write(*,*) 'elem: nord_add_local  = ', nord_add_local
!$omp end critical
            go to 10
         endif
!
         call elem_uw_elasticity_DPG(nord_add_local,nrTest,nrTrial,  &
                                     nrdofHH,nrdofH,nrdofV,nrdofQ,   &
                                     ndofHmdl,ndofVmdl)
!
      case default
         write(*,*) 'elem: Mdle,NODES(Mdle)%case = ',  &
                     Mdle,NODES(Mdle)%case
         call logic_error(ERR_INVALID_VALUE, __FILE__,__LINE__)
      end select

!
   end subroutine elem




!--------------------------------------------------------------------------
!> Purpose :   Element stiffness matrix and load vector for UW DPG
!!             elasticity problem
!!
!> @param[in]  Mdle      - middle node number
!> @param[in]  Dp        - DPG enrichment order
!> @param[in]  NrTest    - number of test functions
!> @param[in]  NrTrial   - number of trial functions
!> @param[in]  NrdofHH   - number of enriched H1 shape functions
!> @param[in]  NrdofH    - number of H1 shape functions
!> @param[in]  NrdofV    - number of Hdiv shape functions
!> @param[in]  NrdofQ    - number of L2 shape functions
!> @param[in]  NrdofHmdl - number of H1 bubble shape functions
!> @param[in]  NrdofVmdl - number of Hdiv bubble shape functions
!!
!> @date       July 2023
!--------------------------------------------------------------------------
!
!                |   \tau \in H(div)^3   |       v \in (H1)^3     |
!
!                   - <\hat u,(\tau n)>  +            0
!  +                        0            + - <\hat (\sigma n),v>
!  + \int_\Omega [   u \cdot div(\tau)   +            0           ]
!  + \int_\Omega [    A \sigma : \tau    +    \sigma : grad(v)    ]
!  + \int_\Omega [     \omega : \tau     +            0           ]
!  = \int_\Omega [          0            +        f \cdot v       ]
!
!--------------------------------------------------------------------------
   subroutine elem_uw_elasticity_DPG(Mdle,dp,NrTest,NrTrial,         &
                                     NrdofHH,NrdofH,NrdofV,NrdofQ,   &
                                     NrdofHmdl,NrdofVmdl)
!
      use control,                  only: INTEGRATION
      use parameters
      use parametersDPG
      use data_structure3D
      use element_data
      use isotropic_elast_material
      use physics,                  only: NR_PHYSA
      use assembly,                 only: ALOC,BLOC,NR_RHS
      use common_prob_data,         only: SYMMETRY_TOL, TEST_NORM, ALPHA
!
      implicit none
!
      integer, intent(in) :: Mdle
      integer, intent(in) :: Dp
      integer, intent(in) :: NrTest, NrTrial
      integer, intent(in) :: NrdofHH
      integer, intent(in) :: NrdofH, NrdofV, NrdofQ
      integer, intent(in) :: NrdofHmdl, NrdofVmdl
!
!  ...element and face data
      integer :: ntype,ftype
      integer :: nrv,nre,nrf
      integer :: norder(19), nordf(5), nordP
      integer :: nedge_orient(12),nface_orient(6)
!
!  ...shape functions
      real(8) :: shapH  (  MAXbrickH), gradH(3,MAXbrickH)
      real(8) :: shapV  (3,MAXbrickV), divV (  MAXbrickV)
      real(8) :: shapV_n(  MAXbrickV)
      real(8) :: shapQ  (  MAXbrickQ)
      real(8) :: shapHH (  MAXbrickHH), gradHH(3,MAXbrickHH)
!
      integer :: nrdofHi, nrdofVi
      integer :: ndofH, ndofV, ndofQ, ndofHH
!
!  ...geometry
      real(8) :: xnod(3,MAXbrickH)
      real(8) :: xi(3), x(3), rn(3)
      real(8) :: dxdxi(3,3), dxidx(3,3)
      real(8) :: t(2), rjac, bjac
      real(8) :: dxidt(3,2), dxdt(3,2), rt(3,2)
      integer :: nsign
!
!  ...tensors in physical coordinates
      real(8), dimension(3,3,3,3) :: A,AA,Symm
!
!  ...temporary variables
      real(8) :: tmp
      real(8) :: tmpDivTau1(3), tmpDivTau2(3), tmpTau_n(3)
!
!  ...source term (don't need Neumann term)
      real(8) :: fval(3,NR_RHS)
!
!  ...quadrature data
      real(8) :: xiloc(3,MAXNINT3ADD), wxi(MAXNINT3ADD)
      real(8) :: tloc(2,MAXNINT2ADD),  wt(MAXNINT2ADD)
      real(8) :: weight, wa
      integer :: nint
!
!  ...matrices
      real(8), allocatable :: EnrTraceDispl(:,:)
      real(8), allocatable :: EnrTraceStress(:,:)
      real(8), allocatable :: EnrFieldDispl(:,:)
      real(8), allocatable :: EnrFieldStress(:,:)
      real(8), allocatable :: EnrLoad(:,:)
      real(8), allocatable :: Stiff_ALL(:,:)
      real(8), allocatable :: FullDPG(:,:)
      real(8), allocatable :: Gram(:,:)
!
!  ...miscellaneous
      integer :: i, j, k, l, m, n, mm, nn
      integer :: k1, k2, k3, k4, k5
      integer :: m1, m2, m3, m4, m5
      integer :: n1, n2, n3
      integer :: i1, j1, j2, j3, j4, j1, j12, j123, j1234
      integer :: ic, jc, kc, lc, klc, ijc
      integer :: ipt, ifc, iflag, info
      integer :: kH, kV, kQ, lH, lV, lQ
      integer :: ndofphys(NR_PHYSA)
      real(8) :: diffmax,dmax
!
      integer :: iprint = 0
!
      integer :: nk
      nk(k1,k2) = (k2-1)*k2/2+k1
!
!--------------------------------------------------------------------------
!
!  ...element type
      ntype = NODES(Mdle)%ntype
      nrv = nvert(ntype); nre = nedge(ntype); nrf = nface(ntype)
!
!  ...order of approximation, orientations, geometry dof's (don't need bc flags)
      call find_order (Mdle, norder)
      call find_orient(Mdle, nedge_orient,nface_orient)
      call nodcor     (Mdle, xnod)
!
!  ...set the enriched order of appoximation
      select case(ntype)
      case(MDLB);       nordP = NODES(Mdle)%order + Dp*111
      case(MDLP);       nordP = NODES(Mdle)%order + Dp*11
      case(MDLN,MDLD);  nordP = NODES(Mdle)%order + Dp*1
      end select
!
      nrdofHi = NrdofH - NrdofHmdl
      nrdofVi = NrdofV - NrdofVmdl
!
!  ...initialize the enriched local element stiffness matrices and load vectors
      allocate(EnrTraceDispl (3*nrdofHi,NrTest));  EnrTraceDispl (:,:) = ZERO
      allocate(EnrTraceStress(3*nrdofVi,NrTest));  EnrTraceStress(:,:) = ZERO
      allocate(EnrFieldDispl (3*NrdofQ, NrTest));  EnrFieldDispl (:,:) = ZERO
      allocate(EnrFieldStress(6*NrdofQ, NrTest));  EnrFieldStress(:,:) = ZERO
      allocate(EnrLoad(NrTest,NR_RHS));            EnrLoad(:) = ZERO
!
!  ...initialize the Gram matrix
      allocate(Gram(NrTest,NrTest))
      Gram(:,:) = ZERO
!
!--------------------------------------------------------------------------
!      E L E M E N T    I N T E G R A L
!--------------------------------------------------------------------------
!
!  ...set up the element quadrature
      INTEGRATION = Dp
      call set_3Dint(ntype,norder, nint,xiloc,wxi)
      INTEGRATION = 0
!
!  ...Auxiliary tensors
      call getSymm(Symm)
!
!  ...loop through integration points
      do ipt=1,nint
         xi(1:3) = xiloc(1:3,ipt); wa = wxi(ipt)
!
!     ...geometry shape functions
         call shape3DH(ntype,xi,norder,nedge_orient,nface_orient,  &
                      ndofH,shapH,gradH)
!
!     ...L2 shape functions
         call shape3DQ(ntype,xi,norder, ndofQ,shapQ)
!
!     ...enriched H1 shape functions
         call shape3HH(ntype,xi,nordP, ndofHH,shapHH,gradHH)
!
!     ...geometry map
         call geom3D(Mdle,xi,xnod,shapH,gradH,ndofH,  &
                     x,dxdxi,dxidx,rjac,iflag)
!
         if (iflag.ne.0) then
            write(*,1000) Mdle,rjac
 1000       format(' Negative Jacobian! Mdle,rjac = ',i8,2x,e12.5)
            stop
         endif
!
!     ...Apply pullbacks
!        L2 (trial)
         shapQ(1:NrdofQ) = shapQ(1:NrdofQ)/rjac
!
! TODO: call Lapack routine
!     ...H1 (test)
         do k=1,NrdofHH
            gradHH(1:3,k) = gradHH(1,k) * dxidx(1,1:3)  &
                          + gradHH(2,k) * dxidx(2,1:3)  &
                          + gradHH(3,k) * dxidx(3,1:3)
         enddo
!
!     ...integration weight
         weight = wa*rjac
!
!     ...compute the compliance tensor
         call getA(x, A)
!
!     ...need this for the adjoint graph norm
         call getAA(X, AA)
!
!     ...get the source term
         call getf(Mdle,x, fval)
!
!    P A R T  1 : go through \tau\in H(div)^3 test space
!
!     ...FIRST OUTER loop through enriched stress dofs
         do k1 = 1,NrdofHH
!        ...OUTER loop through components
            do jc=1,3
               do ic=1,jc
                  ijc = nk(ic,jc)
                  m1 = (k1-1)*6+ijc
!
                  tmpDivTau1=0.d0
                  tmpDivTau1(ic) = gradHH(jc,k1)
                  tmpDivTau1(jc) = gradHH(ic,k1)
!
!     --- Gram matrix ---
!
                  select case(TEST_NORM)
!
!   (A:\tau_1,A:\tau_2+\varepsilon(v_2))+(div(\tau_1),div(tau_2))
!
                  case(1)
!
!                    Tau loop
!
                     do m2 = 1,m1
                        k2 = int((m2-1)/6)+1
                        klc = mod(m2-1,6)+1
!
                        select case(klc)
                        case(1); kc=1; lc=1;
                        case(2); kc=1; lc=2;
                        case(3); kc=2; lc=2;
                        case(4); kc=1; lc=3;
                        case(5); kc=2; lc=3;
                        case(6); kc=3; lc=3;
                        end select
!
                        tmpDivTau2=0.d0
                        tmpDivTau2(kc) = gradHH(lc,k2)
                        tmpDivTau2(lc) = gradHH(kc,k2)
!
                        tmp=0.d0
!
!                    ...L2-terms
                        if (ic.ne.jc) then
                           if (kc.ne.lc) then
                              tmp = AA(kc,lc,ic,jc) * shapHH(k2) * shapHH(k1)    &
                                  + AA(kc,lc,jc,ic) * shapHH(k2) * shapHH(k1)    &
                                  + AA(lc,kc,ic,jc) * shapHH(k2) * shapHH(k1)    &
                                  + AA(lc,kc,jc,ic) * shapHH(k2) * shapHH(k1)
                           else
                              tmp = AA(kc,lc,ic,jc) * shapHH(k2) * shapHH(k1)    &
                                  + AA(kc,lc,jc,ic) * shapHH(k2) * shapHH(k1)
                           endif
                        elseif (kc.ne.lc) then
                           tmp = AA(kc,lc,ic,jc) * shapHH(k2) * shapHH(k1)       &
                               + AA(lc,kc,ic,jc) * shapHH(k2) * shapHH(k1)
                        else
                           tmp = AA(kc,lc,ic,jc) * shapHH(k2) * shapHH(k1)
                        endif
!
!                    ...div term
                        do m = 1,3
                           tmp = tmp + tmpDivTau1(m) * tmpDivTau2(m)
                        enddo
!
!                    ...additional L2-terms
                        if (ic.eq.jc) then
                           if ((kc.eq.ic).and.(lc.eq.jc)) then
                              tmp = tmp + ALPHA * shapHH(k1) * shapHH(k2)
                           endif
                        elseif ((kc.eq.ic).and.(lc.eq.jc)) then
                           tmp = tmp + 2 * ALPHA * shapHH(k1) * shapHH(k2)
                        endif
!
                        Gram(m2,m1) = Gram(m2,m1) + tmp*weight
                     enddo
!
!                    v-loop
!
                     do m2 = 6*NrdofHH+1,6*NrdofHH+3*NrdofHH
                        k2 = int((m2-6*NrdofHH-1)/3)+1
                        kc = mod(m2-6*NrdofHH-1,3)+1

                        tmp=0.d0
                        if (ic.eq.jc) then
                           do m=1,3
                              tmp = tmp  &
                                  + A(ic,jc,kc,m) * gradHH(m,k2) * shapHH(k1)
                           enddo
                        else
                           do m=1,3
                              tmp = tmp  &
                                  + A(ic,jc,kc,m) * gradHH(m,k2) * shapHH(k1)  &
                                  + A(jc,ic,kc,m) * gradHH(m,k2) * shapHH(k1)
                           enddo
                        endif
!
                        Gram(m2,m1) = Gram(m2,m1) + tmp * weight
                     enddo
!
!   (\tau_1,\tau_2)+(div(\tau_1),div(tau_2))
!
                  case(2)
                     do m2=1,m1
                        k2 = int((m2-1)/6)+1
                        klc = mod(m2-1,6)+1
!
                        select case(klc)
                        case(1); kc=1; lc=1;
                        case(2); kc=1; lc=2;
                        case(3); kc=2; lc=2;
                        case(4); kc=1; lc=3;
                        case(5); kc=2; lc=3;
                        case(6); kc=3; lc=3;
                        end select
!
                        tmpDivTau2=0.d0
                        tmpDivTau2(kc) = gradHH(lc,k2)
                        tmpDivTau2(lc) = gradHH(kc,k2)
!
                        tmp=0.d0
!
!                    ...L2 terms
                        if (ic.eq.jc) then
                           if ((kc.eq.ic).and.(lc.eq.jc)) then
                              tmp = shapHH(k1) * shapHH(k2)
                           endif
                        elseif ((kc.eq.ic).and.(lc.eq.jc)) then
                           tmp = 2*shapHH(k1) * shapHH(k2)
                        endif
!
!                    ...Div term
                        do m=1,3
                           tmp = tmp + tmpDivTau1(m) * tmpDivTau2(m)
                        enddo
!
                        Gram(m2,m1) = Gram(m2,m1) + tmp*weight
                     enddo
                  end select
!
!     --- cauchy stress stiffness matrix ---
!
!   A \sigma : \tau
!
!          ( sigma1  sigma2  sigma4 )
!  sigma = ( sigma2  sigma3  sigma5 )
!          ( sigma4  sigma5  sigma6 )
!
!              ...INNER loop through trial dofs for Cauchy stress
                  do k3=1,NrdofQ
                     do klc=1,6
                        m3 = (k3-1)*6+klc
!
                        if (ic.ne.jc) then
                           if (klc.eq.1) then
                              tmp = A(1,1,ic,jc) * shapQ(k3) * shapHH(k1)  &
                                  + A(1,1,jc,ic) * shapQ(k3) * shapHH(k1)
                           elseif (klc.eq.2) then
                              tmp = A(1,2,ic,jc) * shapQ(k3) * shapHH(k1)  &
                                  + A(1,2,jc,ic) * shapQ(k3) * shapHH(k1)  &
                                  + A(2,1,ic,jc) * shapQ(k3) * shapHH(k1)  &
                                  + A(2,1,jc,ic) * shapQ(k3) * shapHH(k1)
                           elseif (klc.eq.3) then
                              tmp = A(2,2,ic,jc) * shapQ(k3) * shapHH(k1)  &
                                  + A(2,2,jc,ic) * shapQ(k3) * shapHH(k1)
                           elseif (klc.eq.4) then
                              tmp = A(1,3,ic,jc) * shapQ(k3) * shapHH(k1)  &
                                  + A(1,3,jc,ic) * shapQ(k3) * shapHH(k1)  &
                                  + A(3,1,ic,jc) * shapQ(k3) * shapHH(k1)  &
                                  + A(3,1,jc,ic) * shapQ(k3) * shapHH(k1)
                           elseif (klc.eq.5) then
                              tmp = A(2,3,ic,jc) * shapQ(k3) * shapHH(k1)  &
                                  + A(2,3,jc,ic) * shapQ(k3) * shapHH(k1)  &
                                  + A(3,2,ic,jc) * shapQ(k3) * shapHH(k1)  &
                                  + A(3,2,jc,ic) * shapQ(k3) * shapHH(k1)
                           elseif (klc.eq.6) then
                              tmp = A(3,3,ic,jc) * shapQ(k3) * shapHH(k1)  &
                                  + A(3,3,jc,ic) * shapQ(k3) * shapHH(k1)
                           endif
                        else
                           if (klc.eq.1) then
                              tmp = A(1,1,ic,jc) * shapQ(k3) * shapHH(k1)
                           elseif (klc.eq.2) then
                              tmp = A(1,2,ic,jc) * shapQ(k3) * shapHH(k1)  &
                                  + A(2,1,ic,jc) * shapQ(k3) * shapHH(k1)
                           elseif (klc.eq.3) then
                              tmp = A(2,2,ic,jc) * shapQ(k3) * shapHH(k1)
                           elseif (klc.eq.4) then
                              tmp = A(1,3,ic,jc) * shapQ(k3) * shapHH(k1)  &
                                  + A(3,1,ic,jc) * shapQ(k3) * shapHH(k1)
                           elseif (klc.eq.5) then
                              tmp = A(2,3,ic,jc) * shapQ(k3) * shapHH(k1)  &
                                  + A(3,2,ic,jc) * shapQ(k3) * shapHH(k1)
                           elseif (klc.eq.6) then
                              tmp = A(3,3,ic,jc) * shapQ(k3) * shapHH(k1)
                           endif
                        endif

                        EnrFieldStress(m3,m1) = EnrFieldStress(m3,m1)  &
                                              + tmp * weight
                     enddo
                  enddo
!
!     --- displacement stiffness matrix --
!
!   u \cdot div(\tau)
!
!              ...INNER loop through trial dofs for displacement
                  do k4=1,NrdofQ
                     do kc=1,3
                        m4 = (k4-1)*3 + kc
!
                        EnrFieldDispl(m4,m1) = EnrFieldDispl(m4,m1)  &
                                             + shapQ(k4) * tmpDivTau1(kc) * weight
                     enddo
                  enddo
!           ...END OUTER LOOP through test stresses
               enddo
            enddo
         enddo
!
!    P A R T  2 : go through v\in(H1)^3 test space
!
!     ...SECOND OUTER loop through enriched H1 dofs
         do k1=1,NrdofHH
!        ...OUTER loop through components
            do ic=1,3
!
!           ...counter of row
               m1 = 6*NrdofHH+(k1-1)*3+ic
!
!     --- load vector ---
!
!           ...f \cdot v
               EnrLoad(m1,1:NR_RHS) = EnrLoad(m1,1:NR_RHS) &
                                    + fval(ic,1:NR_RHS) * shapHH(k1) * weight
!
!     --- Gram matrix ---
!
               select case(TEST_NORM)
!
!   (\varepsilon(v_1),A:\tau_2+\varepsilon(v_2))+(v_1,v_2)
!
               case(1)
                  do m2=6*NrdofHH+1,m1
                     k2 = int((m2-6*NrdofHH-1)/3)+1
                     jc = mod(m2-1,3)+1
!
                     tmp=0.d0
!                 ...grad terms
                     do m=1,3
                        do n=1,3
                           tmp = tmp  &
                               + Symm(jc,m,ic,n) * gradHH(m,k2) * gradHH(n,k1)
                        enddo
                     enddo
!                 ...L2-term
                     if (ic.eq.jc) then
                        tmp = tmp + shapHH(k1)*shapHH(k2)
                     endif
!
                     Gram(m2,m1) = Gram(m2,m1) + tmp*weight
                  enddo
!
!   (v_2,v) + (grad(v_2),grad(v))
!
               case(2)
                  do m2=6*NrdofHH+1,m1
                     jc = mod(m2-1,3)+1
                     if (ic.eq.jc) then
                        k2 = int((m2-6*NrdofHH-1)/3)+1
                        Gram(m2,m1) = Gram(m2,m1)                          &
                                + ( shapHH(k1) * shapHH(k2)                &
                                  + gradHH(1,k1) * gradHH(1,k2)            &
                                  + gradHH(2,k1) * gradHH(2,k2)            &
                                  + gradHH(3,k1) * gradHH(3,k2) ) * weight
                     endif
                  enddo
               end select
!
!     --- cauchy stress stiffness matrix ---
!
!   + \sigma : \varepsilon(v)
!
!          ( sigma1  sigma2  sigma4 )
!  sigma = ( sigma2  sigma3  sigma5 )
!          ( sigma4  sigma5  sigma3 )
!
!           ...INNER loop through trial dofs for Cauchy stress
               do k3=1,NrdofQ
                  do klc=1,6
                     m3 = (k3-1)*6 + klc
!
                     select case(klc)
                     case(1); kc=1; lc=1;
                     case(2); kc=1; lc=2;
                     case(3); kc=2; lc=2;
                     case(4); kc=1; lc=3;
                     case(5); kc=2; lc=3;
                     case(6); kc=3; lc=3;
                     end select
!
                     tmp=0.d0
                     if (kc.eq.ic) then
                        tmp = gradHH(lc,k1) * shapQ(k3)
                     elseif (lc.eq.ic) then
                        tmp = gradHH(kc,k1) * shapQ(k3)
                     endif

                     EnrFieldStress(m3,m1) = EnrFieldStress(m3,m1)  &
                                           + tmp * weight
                  enddo
               enddo
!
!        ...end outer loop
            enddo
         enddo
!
!  ...end of loop through integration points
      enddo
!
      j1 = 6*NrdofHH; j2 = 6*NrdofHH+3*NrdofHH
      Gram(1:j1,j1+1:j2) = transpose(Gram(j1+1:j2,1:j1))
!
!--------------------------------------------------------------------------
!     B O U N D A R Y    I N T E G R A L
!--------------------------------------------------------------------------
!
!  ...loop through element faces
      do ifc=1,nrf
!
!     ...sign factor to determine the OUTWARD normal unit vector
         nsign = nsign_param(ntype,ifc)
!
!     ...face type
         ftype = face_type(ntype,ifc)
!
!     ...face order of approximation
         call face_order(ntype,ifc,norder, nordf)
!
!     ...set up the face quadrature
         INTEGRATION = NORD_ADD
         call set_2Dint(ftype,nordf, nint,tloc,wt)
         INTEGRATION = 0
!
!     ...loop through face integration points
         do ipt=1,nint
!
!        ...face coordinates
            t(1:2) = tloc(1:2,ipt); wa = wt(ipt)
!
!        ...master element coordinates using face parameterization
            call face_param(ntype,ifc,t, xi,dxidt)
!
!        ...Compute shape functions needed for test/trial field variables and geometry
!           H1 (geometry and bdry trial)
            call shape3DH(ntype,xi,norder,nedge_orient,nface_orient,  &
                          ndofH,shapH,gradH)
!        ...H(div) (bdry trial)
            call shape3DV(ntype,xi,norder,nface_orient,  &
                          ndofV,shapV,divV)
!
! TODO: norderi and get interface only
!        ...H1 (test)
            call shape3HH(ntype,xi,nordP, ndofHH,shapHH,gradHH)
!
!        ...geometry map
            call bgeom3D(Mdle,xi,xnod,shapH,gradH,ndofH,dxidt,nsign,  &
                         x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
            weight = wa*bjac
!
!        ...apply pullbacks
! TODO: use Lapack
!           H(div) (trial)
            do k=1,nrdofVi
               shapV(1:3,k) = dxdxi(1:3,1)*shapV(1,k)  &
                            + dxdxi(1:3,2)*shapV(2,k)  &
                            + dxdxi(1:3,3)*shapV(3,k)
               shapV_n(k) = shapV(1,k)*rn(1)  &
                          + shapV(2,k)*rn(2)  &
                          + shapV(3,k)*rn(3)
            enddo
            shapV_n(1:nrdofVi) = shapV_n(1:nrdofVi)/rjac
!
!    P A R T  1 : go through \tau test space
!
!        ...OUTER loop through enriched H(div) test functions
            do k1=1,NrdofHH
               do jc=1,3
                  do ic=1,jc
                     ijc = nk(ic,jc)
                     m1 = (k1-1)*6+ijc

                     tmpTau_n(1:3)=0.d0
                     tmpTau_n(ic) = shapHH(k1)*rn(jc)
                     tmpTau_n(jc) = shapHH(k1)*rn(ic)
!
!        --- displacement stiffness matrix ---
!
!   - <\hat u,(\tau n)>
!
!                 ...INNER loop through H1 bdry trial dofs
                     do k2=1,nrdofHi
                        do kc=1,3
                           m2 = (k2-1)*3+kc
!
                           EnrTraceDispl(m2,m1) = EnrTraceDispl(m2,m1)  &
                                                - shapH(k2) * tmpTau_n(kc) * weight
                        enddo
                     enddo
!
!              ...END OUTER LOOP
                  enddo
               enddo
            enddo
!
!
!    P A R T  2 : go through v\in(H1)^3 test space 
!
!        ...SECOND OUTER loop through enriched H1 test function
            do k1=1,NrdofHH
!           ...OUTER loop through components
               do jc=1,3
                  m1 = 6*NrdofHH+(k1-1)*3+jc
!
!        --- cauchy stress stiffness matrix ---
!
!   - <\hat \sigma,v>
!
!              ...INNER loop through H(div) bdry trial dofs
                  do k3=1,nrdofVi
                     do ic=1,3
                        m3 = (k3-1)*3+ic
                        if (ic.eq.jc) then
                           EnrTraceStress(m3,m1) = EnrTraceStress(m3,m1)  &
                                                 - shapV_n(k3) * shapHH(k1) * weight
                        endif
                     enddo
                  enddo
!
!           ...END OUTER LOOP
               enddo
            enddo
!     ...end of loop over integration points
         enddo
!
!  ...end of loop over faces
      enddo
!
#if DEBUG_MODE
      if (iprint.eq.1) then
         write(*,*) 'Gram = '
         do i=1,25
            write(*,6000) i,Gram((i+1)*i/2-(i-1):(i+1)*i/2)
 6000       format('i = ',i3,'  ',25e12.5)
         enddo
!
         call pause
!
         write(*,*) 'EnrFieldDispl = '
         do i=1,6*NrdofHH+1
            write(*,6001) i,EnrFieldDispl(i,1:3*NrdofQ)
 6001       format('i = ',i4,'  ',15(/,10e12.5))
         enddo
!
         call pause
!
         write(*,*) 'EnrFieldStress = '
         do i=1,6*NrdofHH
           write(*,6001) i,EnrFieldStress(1:6*NrdofQ,i)
         enddo
!
         call pause
!
         do i=1+6*NrdofHH,6*NrdofHH+3*NrdofHH
            write(*,6001) i,EnrFieldStress(1:6*NrdofQ,i)
         enddo
!
         call pause
!
         write(*,*) 'EnrTraceDispl = '
         do i=1,6*NrdofHH+1
            write(*,6001) i,EnrTraceDispl(1:3*nrdofHi,i)
         enddo
!
         call pause
!
         write(*,*) 'EnrTraceStress = '
         do i=6*NrdofHH,6*NrdofHH+3*NrdofHH
            write(*,6001) i,EnrTraceStress(1:3*nrdofVi,i)
         enddo
!
         call pause
!
      endif
#endif
!
!--------------------------------------------------------------------------
!     D P G    F I N A L    A S S E M B L Y
!--------------------------------------------------------------------------
!
!  ...Create vector of indices with dof of each physics variable
      ndofphys = (/3*nrdofHi,3*nrdofVi,3*NrdofQ,6*NrdofQ/)
!
!  ...Construct Stiff_ALL by packing all Enriched Matrices
      i1 = NrTest
      j1 = 3*nrdofHi; j2 = 3*nrdofVi; j3 = 3*NrdofQ; j4 = 6*NrdofQ
      j12   = j1 + j2
      j123  = j12 + j3
      j1234 = j123 + j4
!
      allocate(Stiff_ALL(j1,j1234+NR_RHS))
      Stiff_ALL(1:i1,1:j1)                 = transpose(EnrTraceDispl (1:j1,1:i1))
      Stiff_ALL(1:i1,j1+1:j1+j2)           = transpose(EnrTraceStress(1:j2,1:i1))
      Stiff_ALL(1:i1,j12+1:j12+j3)         = transpose(EnrFieldDispl (1:j3,1:i1))
      Stiff_ALL(1:i1,j123+1:j123+j4)       = transpose(EnrFieldStress(1:j4,1:i1))
      Stiff_ALL(1:i1,j1234+1:j1234+NR_RHS) = EnrLoad(1:i1,1:NR_RHS)
!
      deallocate(EnrTraceDispl,EnrTraceStress,EnrFieldDispl,EnrFieldStress)
      deallocate(EnrLoad)
!
      N     = NrTest
      NRHS  = j1234 + NR_RHS
!
!  ...G^-1 * Stiff_ALL
      call DPOTRF('U',N,Gram,N,info)
      if (info.ne.0) then
         write(*,*) 'elem_DPG_UWEAK_SYMM: info = ',info ; stop
      endif
!
!  ...Build full DPG matrix (stiffness + load) in one go
      call DTRSM('L','U','T','N',N,NRHS,ZONE,Gram,N,Stiff_ALL,N)
      deallocate(Gram); allocate(FullDPG(NRHS,NRHS))
      call DHERK('U','T',NRHS,N,ZONE,Stiff_ALL,N,ZERO,FullDPG,NRHS)
      deallocate(Stiff_ALL)
!
!  ...Fill lower triangular part
      do i=1,j1234
         FullDPG(i+1:NRHS,i) = FullDPG(i,i+1:NRHS)
      enddo
!
!  ...Populate the assembly matrices accordingly
      n1 = 0
      do i = 1,NR_PHYSA
         n = ndofphys(i)
         n2 = n1+n
         m1 = 0
         do j = 1,NR_PHYSA
            m = ndofphys(j)
            m2 = m1 + m
!        ...First initialize
            ALOC(i,j)%array(:,:) = ZERO
!        ...Then populate
            ALOC(i,j)%array(1:n,1:m) = FullDPG(n1+1:n2,m1+1:m2)
            m1 = m2
         enddo
         m2 = m1 + NR_RHS
         BLOC(i)%array(1:n,1:NR_RHS) = FullDPG(n1+1:n2,m1+1:m2)
         n1 = n2
      enddo
!
      deallocate(FullDPG)
!
   end subroutine elem_uw_elasticity_DPG





!----------------------------------------------------------------------
!> @brief       Compute enriched order for all element nodes
!!
!> @param[in]   Nord    - enriched element order
!> @param[out]  Norder  - integer array containing order of elem nodes
!!
!> @date        July 2023
!----------------------------------------------------------------------
   subroutine compute_enriched_order(EType,Nord, Norder)
!
      use node_types, only : MDLB, MDLP, MDLN
      use parameters, only : MODORDER
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
!  ...see implementation of BrokenExactSequence in shape functions
      select case(Etype)
         case(MDLB)
            call decod(Nord,MODORDER,2, temp) !xy face, z edge
            nordF(1) = temp(1); nordB(3) = temp(2)
            call decod(nordF(1),MODORDER,2, nordB(1:2)) !x edge, y edge
            call encod((/nordB(1),nordB(3)/),MODORDER,2, nordF(2)) !xz face
            call encod(nordB(2:3),MODORDER,2, nordF(3)) !yz face
!        ...edges
            Norder(1:4)   = (/nordB(1),nordB(2),nordB(1),nordB(2)/) !x,y,x,y edges
            Norder(5:8)   = (/nordB(1),nordB(2),nordB(1),nordB(2)/) !x,y,x,y edges
            Norder(9:12)  = nordB(3) !z edges
!        ...faces
            Norder(13:14) = nordF(1) !xy faces
            Norder(15:18) = (/nordF(2),nordF(3),nordF(2),nordF(3)/) !xz,yz,xz,yz faces
!        ...element interior
            Norder(19)    = Nord
         case(MDLP)
            call decod(Nord,MODORDER,2, nordB(1:2))
!        ...edges (xy)
            norder(1:6)   = nordB(1)
!        ...edges (z)
            norder(7:9)   = nordB(2)
!        ...triangular faces
            norder(10:11) = nordB(1)
!        ...quadrilateral faces (z)
            norder(12:14) = Nord
!        ...element interior
            norder(15)    = Nord
         case(MDLN)
            norder(1:11)  = Nord
         case default
            write(*,*) 'compute_enriched_order: only available for hexa and prism. stop.'
            stop
      end select
!
   end subroutine compute_enriched_order
!

