!
#include "typedefs.h"
!
!------------------------------------------------------------------------------
!> @brief      Evaluates unconstrained stiffness matrix and load vector
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
!> @date       July 2023
!------------------------------------------------------------------------------
   subroutine elem_maxwell(Mdle,                         &
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
      use mpi_param
!
      implicit none
!
!  ...declare input/output variables
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
!  ...declare edge/face type variables
      integer :: ntype,ftype
!
!  ...declare element order, orientation for edges and faces
      integer :: norder(19), norient_edge(12), norient_face(6)
!
!  ...element nodes order (trial) for interfaces
      integer :: norderi(19), norderf(5)
!
!  ...geometry dof (work space for nodcor)
      real(8) :: xnod(3,MAXbrickH)
!
!  ...geometry
      real(8) :: xi(3), x(3), rn(3)
      real(8) :: dxidt(3,2), dxdt(3,2)
      real(8) :: dxdxi(3,3), dxidx(3,3)
      real(8) :: t(2), rjac, bjac
!
!  ...H1 shape functions
      real(8) :: shapH(MAXbrickH),  gradH(3,MAXbrickH)
!
!  ...H(curl) shape functions
      real(8) :: shapE (3,MAXbrickE),  curlE(3,MAXbrickE)
      real(8) :: shapF (3,MAXbrickE),  curlF(3,MAXbrickE)
      real(8) :: shapFi(3,MAXbrickE)
!
!     enriched
      real(8) :: shapEE(3,MAXbrickEE), curlEE(3,MAXbrickEE)
!
!     complex
      complex(8) :: zshapF(3,MAXbrickE), zcurlF(3,MAXbrickE)
      complex(8) :: epsTshapF(3,NrdofEE), epscurlF(3,NrdofEE)
!
!  ...L2 shape functions
      real(8) :: shapQ(MAXbrickQ)
!
!  ...load vector for the enriched space
      complex(8) :: bload_E(NrTest)
!
!  ...Gram matrix in packed format
      complex(8), allocatable :: gramP(:)
!
!  ...Field and material values
      real(8) :: FF
      real(8) :: fldE(3), fldH(3), crlE(3)
      real(8) :: fldF(3), fldG(3), crlF(3), crlG(3)
      complex(8) :: epsTfldE(3), epsTfldF(3), epscrlE(3)
      complex(8) :: eps(3,3)
!
!  ...matrices for transpose filling (swapped loops)
!  ...stiffness matrices (transposed) for the enriched test space
      complex(8), allocatable :: stiff_EE_T(:,:),stiff_EQ_T(:,:)
      complex(8), allocatable :: stiff_ALL(:,:),zaloc(:,:)
!
!  ...quadrature data
      real(8) :: xiloc(3,MAXNINT3ADD), waloc(MAXNINT3ADD)
      real(8) :: tloc(2,MAXNINT2ADD), wtloc(MAXNINT2ADD)
      real(8) :: weight,wa
!
!  ...BC's flags
      integer :: ibc(6,NRINDEX)
!
!  ...for auxiliary computation
      complex(8) :: zaux,zcux
!
!  ...Maxwell load and auxiliary variables
      complex(8) :: zJ(3), zImp(3)
      real(8), dimension(3) :: E1,E2,rntimesE,rn2timesE
!
!  ...number of edge,faces per element type
      integer :: nre, nrf
!
!  ...various variables for the problem
      real(8) :: CC
      integer :: i1, j1, j2, k1, k2, i, k, l, nint, n, m
      integer :: iflag, ifc, info, nrdof, nordP, nsign
!
      complex(8) :: zc1
      complex(8) :: za(3,3), zc(3,3)
!
#if DEBUG_MODE
      integer :: iprint
#endif
!
!  ...for Gram matrix compressed storage format
      integer :: nk
      nk(k1,k2) = (k2-1)*k2/2+k1
!
!-------------------------------------------------------------------------------
!
#if DEBUG_MODE
!  ...Set iprint = 0/1 (Non-/VERBOSE)
      iprint = 0
      if (iprint.eq.1) then
         write(*,*) 'elem_maxwell: Mdle = ', Mdle
      endif
#endif
!
!  ...allocate matrices
      allocate(gramP(NrTest*(NrTest+1)/2))
      allocate(stiff_EE_T(2*NrdofEi,NrTest))
      allocate(stiff_EQ_T(6*NrdofQ ,NrTest))
!
!  ...element type
      ntype = NODES(Mdle)%ntype
      nre = nedge(ntype); nrf = nface(ntype)
!
!  ...determine order of approximation
      call find_order(Mdle, norder)
      norderi(1:nre+nrf) = norder(1:nre+nrf)
!
!  ...set the enriched order of approximation
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
!  ...determine edge and face orientations
      norient_edge(:) = 0
      norient_face(:) = 0
      call find_orient(Mdle, norient_edge,norient_face)
!
!  ...determine nodes coordinates
      call nodcor(Mdle, xnod)
!
!  ...get the element boundary conditions flags
      call find_bc(Mdle, ibc)
!
!  ...clear space for output matrices
      ZblocE(:)    = ZERO; ZblocQ(:)    = ZERO
      ZalocEE(:,:) = ZERO; ZalocEQ(:,:) = ZERO
      ZalocQE(:,:) = ZERO; ZalocQQ(:,:) = ZERO
!
!  ...clear space for auxiliary matrices
      bload_E(:)      = ZERO
      gramP(:)        = ZERO
      stiff_EE_T(:,:) = ZERO
      stiff_EQ_T(:,:) = ZERO
!
!-----------------------------------------------------------------------
!                 E L E M E N T   I N T E G R A L S
!-----------------------------------------------------------------------
!
!  ...use the enriched order to set the quadrature
      INTEGRATION = NORD_ADD
      call set_3D_int_DPG(ntype,norder,norient_face, nint,xiloc,waloc)
      INTEGRATION = 0
!
!
!  ...loop over integration points
      do l=1,nint
!
         xi(1:3)=xiloc(1:3,l); wa=waloc(l)
!
!     ...H1 shape functions (for geometry)
         call shape3DH(ntype,xi,norder,norient_edge,norient_face, nrdofH,shapH,gradH)
!
!     ...L2 shape functions for the trial space
         call shape3DQ(ntype,xi,norder, nrdofQ,shapQ)
!
!     ...broken H(curl) shape functions for the enriched test space
         call shape3EE(ntype,xi,nordP, nrdofEE,shapEE,curlEE)
!
!     ...geometry map
         call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH, x,dxdxi,dxidx,rjac,iflag)
!
!     ...get permittivity at x
         call get_permittivity(mdle,x, eps)
!
#if DEBUG_MODE
         if (iflag .ne. 0) then
            write(*,5999) Mdle,rjac
   5999     format('elem_maxwell: Negative Jacobian. Mdle,rjac=',i8,2x,e12.5)
            stop
         endif
#endif
!
!     ...integration weight
         weight = rjac*wa
!
!     ...get the RHS
         call getf(Mdle,x, zJ)
!
         za = (ZI*OMEGA*EPSILON) * eps(:,:)
         zc = (ZI*OMEGA*MU) * IDENTITY(:,:)
!
!     ...scalar permeability (will occasionally be more convenient)
         zc1 = ZI*OMEGA*MU
!
!     ...apply pullbacks
         call DGEMM('T','N',3,nrdofEE,3,1.d0     ,dxidx,3,shapEE,3,0.d0,shapF,3)
         call DGEMM('N','N',3,nrdofEE,3,1.d0/rjac,dxdxi,3,curlEE,3,0.d0,curlF,3)
!
!     ...apply permittivity
         zshapF = cmplx(shapF,0.d0,8)
         zcurlF = cmplx(curlF,0.d0,8)
         call ZGEMM('C','N',3,nrdofEE,3,ZONE,za,3,zshapF,3,ZERO,epsTshapF,3)
         call ZGEMM('N','N',3,nrdofEE,3,ZONE,za,3,zcurlF,3,ZERO,epscurlF ,3)
!
!     ...loop through enriched H(curl) test functions
         do k1=1,nrdofEE
!
!        ...pickup pulled back test functions
            fldF(:) = shapF(:,k1);  crlF(:) = curlF(:,k1)
            fldG(:) = fldF(:);      crlG(:) = crlF(:)
            epsTfldF(:) = epsTshapF(:,k1)
!
!  --- load ---
!
!           RHS:
!           (J^imp,F) first  equation RHS (with first H(curl) test function F)
!           (0    ,G) second equation RHS is zero
            n = 2*k1-1
            bload_E(n) = bload_E(n)                                   &
                       + (fldF(1)*zJ(1)+fldF(2)*zJ(2)+fldF(3)*zJ(3))  &
                       * weight
!
!
!  --- stiffness matrix ---
!
!        ...loop through L2 trial shape functions
            do k2=1,nrdofQ
!
!           ...first L2 variable
               m = (k2-1)*6
!
!           ...Piola transformation
               fldE(1:3) = shapQ(k2)/rjac; fldH = fldE
!
!           ...-iωε(E,F)
!           ...(H,curl(F))
               n = 2*k1-1
               stiff_EQ_T(m+1:m+3,n) = stiff_EQ_T(m+1:m+3,n) - fldE(:)*conjg(epsTfldF(:))*weight
               stiff_EQ_T(m+4:m+6,n) = stiff_EQ_T(m+4:m+6,n) + fldH(:)*    crlF(:)*weight
!
!           ...(E,curl(G))
!           ...iωμ(H,G)
               n = 2*k1
               stiff_EQ_T(m+1:m+3,n) = stiff_EQ_T(m+1:m+3,n) +     fldE(:)*crlG(:)*weight
               stiff_EQ_T(m+4:m+6,n) = stiff_EQ_T(m+4:m+6,n) + zc1*fldH(:)*fldG(:)*weight
!
!        ...end of loop through L2 trial functions
            enddo
!
!  --- Gram matrix ---
!
!        ...loop through enriched H(curl) test functions
            do k2=k1,nrdofEE
               fldE(:) = shapF(:,k2);          crlE(:) = curlF(:,k2)
               epsTfldE(:) = epsTshapF(:,k2);  epscrlE(:) = epscurlF(:,k2)
!
               call dot_product(fldF,fldE, FF)
               call dot_product(crlF,crlE, CC)
!
!          ...accumulate for the Hermitian Gram matrix
!             (compute upper triangular only)
!          ...testNorm = Scaled Adjoint Graph norm
!                ||v|| = alpha*(v,v) + (A^* v, A^* v)
!             (first eqn multiplied by F, second eqn by G)
!             G_ij=(phi_j,phi_i)_testNorm is 2x2 matrix
!             where (phi_j,phi_i)_l2Norm = Int[phi_i^* phi_j]
!             and phi_i = (F_i,G_i), phi_j = (F_j,G_j).
!             -------------------------
!             | (F_i,F_j)   (F_i,G_j) |
!             | (G_i,F_j)   (G_i,G_j) |
!             -------------------------
!             F_i/G_i are outer loop shape functions (fldF)
!             F_j/G_j are inner loop shape functions (fldE)
!
!             (F_j,F_i) terms = Int[F_^*i F_j] terms (G_11)
               n = 2*k1-1; m = 2*k2-1
!
               k = nk(n,m)
               zaux = conjg(epsTfldF(1))*epsTfldE(1) + &
                      conjg(epsTfldF(2))*epsTfldE(2) + &
                      conjg(epsTfldF(3))*epsTfldE(3)
               gramP(k) = gramP(k) &
                        + (zaux + ALPHA_NORM*FF + CC)*weight
!
!              (G_j,F_i) terms = Int[F_^*i G_j] terms (G_12)
               n = 2*k1-1; m = 2*k2
               k = nk(n,m)
               zaux = - (fldF(1)*epscrlE(1) + &
                         fldF(2)*epscrlE(2) + &
                         fldF(3)*epscrlE(3) )
               zcux = conjg(zc1)*(crlF(1)*fldE(1) + &
                                  crlF(2)*fldE(2) + &
                                  crlF(3)*fldE(3) )
               gramP(k) = gramP(k) + (zaux+zcux)*weight
!
!           ...compute lower triangular part of 2x2 G_ij matrix
!              only if it is not a diagonal element, G_ii
               if (k1 .ne. k2) then
!                 (F_j,G_i) terms = Int[G_^*i F_j] terms (G_21)
                  n = 2*k1; m = 2*k2-1
                  k = nk(n,m)
                  zaux = - (crlF(1)*epsTfldE(1) + &
                            crlF(2)*epsTfldE(2) + &
                            crlF(3)*epsTfldE(3) )
                  zcux = zc1*(fldF(1)*crlE(1) + &
                              fldF(2)*crlE(2) + &
                              fldF(3)*crlE(3) )
                  gramP(k) = gramP(k) + (zaux+zcux)*weight
               endif
!
!              (G_j,G_i) terms = Int[G_^*i G_j] terms (G_22)
               n = 2*k1; m = 2*k2
               k = nk(n,m)
               zcux = abs(zc1)**2*(fldF(1)*fldE(1) + &
                                   fldF(2)*fldE(2) + &
                                   fldF(3)*fldE(3) )
               gramP(k) = gramP(k) &
                        + (zcux + ALPHA_NORM*FF + CC)*weight
!
!        ...end of loop through enriched H(curl) test functions
            enddo
!
!     ...end of loop through enriched H(curl) test functions
         enddo
!
!  ...end of loop through integration points
      enddo
!
!-----------------------------------------------------------------------
!                 B O U N D A R Y   I N T E G R A L S
!-----------------------------------------------------------------------
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
         call face_order(ntype,ifc,norder, norderf)
!
!     ...set 2D quadrature
         INTEGRATION = NORD_ADD
         call set_2D_int_DPG(ftype,norderf,norient_face(ifc), nint,tloc,wtloc)
         INTEGRATION = 0
!
!     ...loop through integration points
         do l=1,nint
!
!        ...face coordinates
            t(1:2) = tloc(1:2,l)
!
!        ...face parametrization
            call face_param(ntype,ifc,t, xi,dxidt)
!
!        ...determine discontinuous Hcurl shape functions
            call shape3EE(ntype,xi,nordP, nrdof,shapEE,curlEE)
#if DEBUG_MODE
            if (nrdof .ne. NrdofEE) then
               write(*,*) 'elem_maxwell: INCONSISTENCY NrdofEE. stop.'
               stop
            endif
#endif
!
!        ...determine element H1 shape functions (for geometry)
            call shape3DH(ntype,xi,norder,norient_edge,norient_face, &
                          nrdof,shapH,gradH)
#if DEBUG_MODE
         if (nrdof .ne. NrdofH) then
            write(*,*) 'elem_maxwell: INCONSISTENCY NrdofH. stop.'
            stop
         endif
#endif
!
!        ...determine element H(curl) shape functions (for fluxes)
!        ...for interfaces only (no bubbles)
            call shape3DE(ntype,xi,norderi,norient_edge,norient_face, &
                          nrdof,shapE,curlE)
#if DEBUG_MODE
         if (nrdof .ne. NrdofEi) then
            write(*,*) 'elem_maxwell: INCONSISTENCY NrdofEi. stop.'
            stop
         endif
#endif
!
!        ...geometry
            call bgeom3D(Mdle,xi,xnod,shapH,gradH,NrdofH,dxidt,nsign, &
                         x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
            weight = bjac*wtloc(l)
!
!        ...pullback trial and test functions
            call DGEMM('T','N',3,nrdofEE,3,1.d0,dxidx,3,shapEE,3,0.d0,shapF ,3)
            call DGEMM('T','N',3,nrdofEi,3,1.d0,dxidx,3,shapE ,3,0.d0,shapFi,3)
!
!        ...loop through enriched H(curl) test functions
            do k1=1,nrdofEE
               E1(1:3) = shapF(:,k1)
!
!           ...check for impedance BC (elimination strategy)
!              (impedance constant is GAMMA for TE10 mode in rectangular waveguide)
!              ( < n x H , G > = GAMMA*< n x n x E , G > + < zg , G > )
               if (ibc(ifc,2).eq.3) then
!              ...get the boundary source [zImp should be zero here]
                  call get_bdSource(Mdle,x,rn, zImp)
!              ...accumulate for the load vector
                  k = 2*k1-1
                  bload_E(k) = bload_E(k) &
                             - (zImp(1)*E1(1)+zImp(2)*E1(2)+zImp(3)*E1(3))*weight
!           ...end if for impedance BC
               endif
!
!           ...loop through H(curl) trial functions
               do k2=1,NrdofEi
                  E2(1:3) = shapFi(:,k2)
                  call cross_product(rn,E2, rntimesE)
!
!              ...check for impedance BC (elimination strategy)
                  if (ibc(ifc,2).eq.3) then
!                 ...accumulate for the extended stiffness matrix on IBC
                     call cross_product(rn,rntimesE, rn2timesE)
                     stiff_EE_T(2*k2-1,2*k1-1) = stiff_EE_T(2*k2-1,2*k1-1) &
                                                  + (E1(1)*rn2timesE(1) &
                                                  +  E1(2)*rn2timesE(2) &
                                                  +  E1(3)*rn2timesE(3) &
                                                    )*GAMMA*weight
                  else
!                 ...accumulate for the extended stiffness matrix without IBC
                     stiff_EE_T(2*k2,2*k1-1) = stiff_EE_T(2*k2,2*k1-1) &
                                                + (E1(1)*rntimesE(1) &
                                                +  E1(2)*rntimesE(2) &
                                                +  E1(3)*rntimesE(3) &
                                                  )*weight
!              ...end if for impedance BC
                  endif
                  stiff_EE_T(2*k2-1,2*k1) = stiff_EE_T(2*k2-1,2*k1) &
                                             + (E1(1)*rntimesE(1) &
                                             +  E1(2)*rntimesE(2) &
                                             +  E1(3)*rntimesE(3) &
                                               )*weight
!           ...end loop through H(curl) trial functions
               enddo
!        ...end loop through the enriched H(curl) test functions
            enddo
!     ...end loop through integration points
         enddo
!  ...end loop through faces
      enddo
!
!-------------------------------------------------------------------------------
!      Construction of the DPG system
!-------------------------------------------------------------------------------
!
      allocate(stiff_ALL(NrTest,NrTrial+1))
!
!  ...Total test/trial DOFs of the element
      i1 = NrTest ; j1 = 2*NrdofEi ; j2 = 6*NrdofQ
!
      stiff_ALL(1:i1,1:j1)       = transpose(stiff_EE_T(1:j1,1:i1))
      stiff_ALL(1:i1,j1+1:j1+j2) = transpose(stiff_EQ_T(1:j2,1:i1))
      stiff_ALL(1:i1,j1+j2+1)    = bload_E(1:i1)
!
      deallocate(stiff_EE_T,stiff_EQ_T)
!
!  ...A. Compute Cholesky factorization of Gram Matrix, G=U^*U (=LL^*)
      call ZPPTRF('U',NrTest,gramP,info)
      if (info.ne.0) then
         write(*,*) 'elem_maxwell: ZPPTRF: Mdle,info = ',Mdle,info,'. stop.'
         stop
      endif
!
!  ...B. Solve triangular system to obtain B~, (LX=) U^*X = [B|l]
      call ZTPTRS('U','C','N',NrTest,NrTrial+1,gramP,stiff_ALL,NrTest,info)
      if (info.ne.0) then
         write(*,*) 'elem_maxwell: ZTPTRS: Mdle,info = ',Mdle,info,'. stop.'
         stop
      endif
!
      deallocate(gramP)
      allocate(zaloc(NrTrial+1,NrTrial+1)); zaloc = ZERO
!
!  ...C. Matrix multiply: B^* G^-1 B (=B~^* B~)
      call ZHERK('U','C',NrTrial+1,NrTest,ZONE,stiff_ALL,NrTest,ZERO,zaloc,NrTrial+1)
!
      deallocate(stiff_ALL)
!
!  ...D. Fill lower triangular part of Hermitian matrix
      do i=1,NrTrial
         zaloc(i+1:NrTrial+1,i) = conjg(zaloc(i,i+1:NrTrial+1))
      enddo
!
!  ...E. Fill ALOC and BLOC matrices
      ZblocE(1:j1) = zaloc(1:j1,j1+j2+1)
      ZblocQ(1:j2) = zaloc(j1+1:j1+j2,j1+j2+1)
!
      ZalocEE(1:j1,1:j1) = zaloc(1:j1,1:j1)
      ZalocEQ(1:j1,1:j2) = zaloc(1:j1,j1+1:j1+j2)
!
      ZalocQE(1:j2,1:j1) = zaloc(j1+1:j1+j2,1:j1)
      ZalocQQ(1:j2,1:j2) = zaloc(j1+1:j1+j2,j1+1:j1+j2)
!
      deallocate(zaloc)
!
!-------------------------------------------------------------------------------
!       I M P E D A N C E   B O U N D A R Y
!-------------------------------------------------------------------------------
!
!  ...Implementation of impedance BC via L2 penalty term
      if (IBCFLAG.eq.2) call imp_penalty(Mdle,NrdofH,NrdofEi,MdE,          &
                                         norder,norderi, ZblocE,ZalocEE)
!
   end subroutine elem_maxwell


!------------------------------------------------------------------------------
!> @brief       Routine adds impedance L2 penalty terms to the stiffness matrix
!               and load vector for the UW DPG Maxwell problem
!!
!> @param[in]   Mdle     - element middle node number
!> @param[in]   NrdofH   - number of H1 trial dof
!> @param[in]   NrdofEi  - number of H(curl) trial interface dof
!> @param[in]   MdE      - num rows of ZalocEE
!> @param[in]   Norder   - element nodes order (trial)
!> @param[in]   Norderi  - element nodes order (trial) for interfaces
!!
!> @param[out]  ZblocE   - load vector
!> @param[out]  ZalocEE  - stiffness matrix
!!
!> @date       July 2023
!------------------------------------------------------------------------------
   subroutine imp_penalty(Mdle,NrdofH,NrdofEi,MdE,Norder,Norderi, &
                        ZblocE,ZalocEE)
!
      use control
      use parametersDPG
      use data_structure3D
      use commonParam
!
      implicit none
!
      integer, intent(in)    :: Mdle
      integer, intent(in)    :: NrdofH
      integer, intent(in)    :: NrdofEi
      integer, intent(in)    :: MdE
      integer, intent(in)    :: Norder(19)
      integer, intent(in)    :: Norderi(19)
      VTYPE,   intent(inout) :: ZblocE(MdE)
      VTYPE,   intent(inout) :: ZalocEE(MdE,MdE)
!
!  ...element and face types
      integer :: ntype,ftype
!
!  ...orientations/order
      integer :: norient_edge(12), norient_face(6)
      integer :: norderf(5)
!
!  ...geometry
      real(8) :: xnod(3,MAXbrickH)
      real(8) :: xi(3), x(3), rn(3)
      real(8) :: dxidt(3,2), dxdt(3,2)
      real(8) :: dxdxi(3,3), dxidx(3,3)
      real(8) :: t(2)
      real(8) :: rjac, bjac
      integer :: nsign
!
!  ...shape functions
      real(8) :: shapH(MAXbrickH),   gradH(3,MAXbrickH)
      real(8) :: shapE(3,MAXbrickE), curlE(3,MAXbrickE)
!
!  ...quadrature
      real(8) :: tloc(2,MAXNINT2ADD), wtloc(MAXNINT2ADD)
      real(8) :: weight
!
!  ...BCs
      integer, dimension(6,NRINDEX)      :: ibc
!
!  ...Maxwell load and auxiliary variables
      VTYPE  , dimension(3) :: zImp
      real(8), dimension(3) :: E2,rntimesE,rn2timesE
      real(8), dimension(3) :: F1,rntimesF,rn2timesF
!
!  ...penalty
      real(8) :: penalty
!
!  ...misc
      integer :: k1,k2,k,l,nint
      integer :: nrdof,ifc,nrf
!
!-------------------------------------------------------------------------------
!
      if (IBCFLAG .ne. 2) then
         write(*,*) 'imp_penalty called for IBCFLAG.ne.2, returning...'
         return
      endif
!
!  ...define the weight of the penalty term
      penalty = 1.d0
!
!  ...no need for overintegration in penalty term
      INTEGRATION = 0
!
!  ...determine element type and number of faces
      ntype = NODES(Mdle)%ntype
      nrf = nface(ntype)
!
!  ...determine edge and face orientations
      call find_orient(Mdle, norient_edge,norient_face)
!
!  ...determine nodes coordinates
      call nodcor(Mdle, xnod)
!
!  ...get the element boundary conditions flags
      call find_bc(Mdle, ibc)
!
!  ...loop through element faces
      do ifc=1,nrf
!
!     ...skip the face if it is not on the impedance boundary
         if (ibc(ifc,2).ne.2) cycle
!
!     ...sign factor to determine the OUTWARD normal unit vector
         nsign = nsign_param(ntype,ifc)
!
!     ...face type
         ftype = face_type(ntype,ifc)
!
!     ...face order of approximation
         call face_order(ntype,ifc,Norder, norderf)
!
!     ...set 2D quadrature
         call set_2D_int(ftype,norderf,norient_face(ifc), nint,tloc,wtloc)
!
!     ...loop through integration points
         do l=1,nint
!
!        ...face coordinates
            t(1:2) = tloc(1:2,l)
!
!        ...face parametrization
            call face_param(ntype,ifc,t, xi,dxidt)
!
!        ...determine element H1 shape functions (for geometry)
            call shape3DH(ntype,xi,Norder,norient_edge,norient_face, &
                          nrdof,shapH,gradH)
!
!        ...determine element H(curl) shape functions (for fluxes)
!        ...for interfaces only (no bubbles)
            call shape3DE(ntype,xi,Norderi,norient_edge,norient_face, &
                          nrdof,shapE,curlE)
!
!        ...geometry
            call bgeom3D(Mdle,xi,xnod,shapH,gradH,NrdofH,dxidt,nsign, &
                         x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
            weight = bjac*wtloc(l)
!
!        ...get the boundary source [zImp should be zero here]
            call get_bdSource(Mdle,x,rn, zImp)
!
!        ...loop through H(curl) test functions
            do k1=1,NrdofEi
               F1(1:3) = shapE(1,k1)*dxidx(1,1:3) &
                       + shapE(2,k1)*dxidx(2,1:3) &
                       + shapE(3,k1)*dxidx(3,1:3)
!
               call cross_product(rn,F1, rntimesF)
               call cross_product(rn,rntimesF, rn2timesF)
!
!           ...compute impedance BC equation for the load
!              (impedance constant is GAMMA for TE10 mode in rectangular waveguide)
!              ( < n x H , G > = GAMMA*< n x n x E , G > + < zg , G > )
!
!           ...accumulate directly into load vector ZblocE(1:2*NrdofEi)
!              (zImp,F): -gamma * (zImp , n x n x F)
               k = 2*k1-1
               ZblocE(k) = ZblocE(k) - (rn2timesF(1)*zImp(1)  &
                                     +  rn2timesF(2)*zImp(2)  &
                                     +  rn2timesF(3)*zImp(3)  &
                                       )*GAMMA*weight/penalty
!              (zImp,G): (zImp , n x G)
               k = 2*k1
               ZblocE(k) = ZblocE(k) + (rntimesF(1)*zImp(1)  &
                                     +  rntimesF(2)*zImp(2)  &
                                     +  rntimesF(3)*zImp(3)  &
                                       )*weight/penalty
!
!           ...loop through H(curl) trial functions
               do k2=1,NrdofEi
                  E2(1:3) = shapE(1,k2)*dxidx(1,1:3) &
                          + shapE(2,k2)*dxidx(2,1:3) &
                          + shapE(3,k2)*dxidx(3,1:3)
!
                  call cross_product(rn,E2, rntimesE)
                  call cross_product(rn,rntimesE, rn2timesE)
!              ...accumulate directly into stiffness ZalocEE(1:2*NrdofEi,1:2*NrdofEi)
!                 4 contributions: trial (E,H), test (F,G)
!                 (E,F):  gamma^2 * (n x n x E , n x n x F)
                  ZalocEE(2*k1-1,2*k2-1) = ZalocEE(2*k1-1,2*k2-1) &
                                           + (rn2timesF(1)*rn2timesE(1)  &
                                           +  rn2timesF(2)*rn2timesE(2)  &
                                           +  rn2timesF(3)*rn2timesE(3)  &
                                             )*GAMMA*GAMMA*weight/penalty
!                 (H,F): -gamma   * (n x H , n x n x F)
                  ZalocEE(2*k1-1,2*k2  ) = ZalocEE(2*k1-1,2*k2  ) &
                                           - (rn2timesF(1)*rntimesE(1)   &
                                           +  rn2timesF(2)*rntimesE(2)   &
                                           +  rn2timesF(3)*rntimesE(3)   &
                                             )*GAMMA*weight/penalty
!                 (E,G): -gamma   * (n x n x E , n x G)
                  ZalocEE(2*k1  ,2*k2-1) = ZalocEE(2*k1  ,2*k2-1) &
                                           - (rntimesF(1)*rn2timesE(1)  &
                                           +  rntimesF(2)*rn2timesE(2)  &
                                           +  rntimesF(3)*rn2timesE(3)  &
                                             )*GAMMA*weight/penalty
!                 (H,G):            (n x H , n x G)
                  ZalocEE(2*k1  ,2*k2  ) = ZalocEE(2*k1  ,2*k2  ) &
                                           + (rntimesF(1)*rntimesE(1)   &
                                           +  rntimesF(2)*rntimesE(2)   &
                                           +  rntimesF(3)*rntimesE(3)   &
                                             )*weight/penalty
!           ...end loop through H(curl) trial functions
               enddo
!        ...end loop through the H(curl) test functions
            enddo
!     ...end loop through integration points
         enddo
!  ...end loop through faces
      enddo
!
   end subroutine imp_penalty
