!
#include "typedefs.h"
!
!------------------------------------------------------------------------------
!> @brief      Evaluates element residual (squared) for UW  Maxwell problem
!!
!> @param[in]  Mdle        - element middle node number
!> @param[in]  NrdofEE     - number of H(curl) test dof
!> @param[in]  NrdofH      - number of H1 trial dof
!> @param[in]  NrdofE      - number of H(curl) trial dof
!> @param[in]  NrdofQ      - number of L2 trial dof
!!
!> @param[out] Resid       - element residual (squared)
!> @param[out] Nref_flag   - suggested h-refinement flag
!!
!> @date       July 2023
!------------------------------------------------------------------------------
   subroutine elem_residual_maxwell(Mdle,NrTest,                    &
                                    NrdofEE,NrdofH,NrdofE,NrdofQ,   &
                                    Resid,Nref_flag)
!
      use control
      use parametersDPG
      use element_data
      use data_structure3D
      use commonParam
      use mpi_param
!
      implicit none
!
!  ...declare input/output variables
      integer, intent(in)  :: Mdle
      integer, intent(in)  :: NrTest
      integer, intent(in)  :: NrdofEE
      integer, intent(in)  :: NrdofH
      integer, intent(in)  :: NrdofE
      integer, intent(in)  :: NrdofQ
      integer, intent(out) :: Nref_flag
      real(8), intent(out) :: Resid
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
!  ...solution dof (work space for solelm)
      complex(8) :: zdofH(MAXeqnH,MAXbrickH)
      complex(8) :: zdofE(MAXeqnE,MAXbrickE)
      complex(8) :: zdofV(MAXeqnV,MAXbrickV)
      complex(8) :: zdofQ(MAXeqnQ,MAXbrickQ)
!
!  ...geometry
      real(8) :: xi(3), x(3), rn(3), daux(3)
      real(8) :: dxidt(3,2), dxdt(3,2), rt(3,2)
      real(8) :: dxdxi(3,3), dxidx(3,3)
      real(8) :: t(2), rjac, bjac
!
!  ...H1 shape functions
      real(8), dimension(MAXbrickH)   :: shapH
      real(8), dimension(3,MAXbrickH) :: gradH
!
!  ...H(curl) shape functions
      real(8), dimension(3,MAXbrickE)  :: shapE, curlE
      real(8), dimension(3,NrdofEE)    :: shapF, curlF
      complex(8), dimension(3,NrdofEE) :: zshapF, zcurlF
      real(8), dimension(3,NrdofE)     :: shapFi
      complex(8), dimension(3,NrdofEE) :: epsTshapE, epsTshapF, epscurlF
!
!  ...L2 shape functions
      real(8), dimension(MAXbrickQ) :: shapQ
!
!  ...enriched Hcurl shape functions
      real(8), dimension(3,MAXbrickEE) :: shapEE
      real(8), dimension(3,MAXbrickEE) :: curlEE
!
!  ...Gram matrix in packed format
      VTYPE, allocatable :: gramP(:)
!
!  ...intermediate values
      real(8) :: FF
      real(8) :: fldE(3), fldH(3), crlE(3), crlH(3)
      real(8) :: fldF(3), fldG(3), crlF(3), crlG(3)
      complex(8) :: epsTfldE(3), epsTfldF(3), epscrlF(3), epscrlE(3)
      complex(8) :: eps(3,3), CF, FC
!
      real(8) :: D_aux(3,3) ,D_za(3,3), D_zc(3,3)
!
!  ...load vector for the enriched space
      VTYPE, dimension(NrTest)   :: bload_E,bload_Ec
      VTYPE, dimension(2*NrdofE) :: bload_Imp
!
!  ...quadrature data
      real(8) :: xiloc(3,MAXNINT3ADD), waloc(MAXNINT3ADD)
      real(8) :: tloc(2,MAXNINT2ADD), wtloc(MAXNINT2ADD)
      real(8) :: weight,wa
!
!  ...BC's flags
      integer :: ibc(6,NRINDEX)
!
!  ...Maxwell load and auxiliary variables
      complex(8) :: zJ(3), zImp(3)
      real(8)    :: E1(3), rntimesE(3), rn2timesE(3)
      real(8)    :: pen
!
!  ...approximate solution
      VTYPE, dimension(3,2) :: zsolExi,zsolE,zflux,zflux2
      VTYPE, dimension(6)   :: zsolQ
!
!  ...auxiliary
      VTYPE :: zresid, zaux, zbux, zcux
!
!  ...number of faces per element type
      integer :: nrf
!
!  ...various variables for the problem
      real(8) :: CC,EE,CE,E,EC,q,h
      integer :: i1,i2,j1,j2,k1,k2,kH,kk,i,j,m,n,nint,kE,k,l,ivar,iflag
      integer :: nordP,nsign,ifc,ndom,info,icomp,nrdof,nrdof_eig,idec
      VTYPE   :: zfval
      VTYPE   :: za(3,3),zc(3,3),zc1
!
!  ...timer
      real(8) :: MPI_Wtime,start_time,end_time
!
!  ...debug variables
#if DEBUG_MODE
      integer :: iprint = 0
#endif
!
!  ...for Gram matrix compressed storage format
      integer :: nk
      nk(k1,k2) = (k2-1)*k2/2+k1
!
!--------------------------------------------------------------------------
!
      allocate(gramP(NrTest*(NrTest+1)/2))
!
!  ...element type
      ntype = NODES(Mdle)%ntype
      nrf = nface(ntype)
!
!  ...determine order of approximation
      call find_order(Mdle, norder)
!  ...set the enriched order of approximation
      select case(ntype)
         case(MDLB);       nordP = NODES(Mdle)%order + NORD_ADD*111
         case(MDLP);       nordP = NODES(Mdle)%order + NORD_ADD*11
         case(MDLD,MDLN);  nordP = NODES(Mdle)%order + NORD_ADD
      end select
!  ...determine edge and face orientations
      call find_orient(Mdle, norient_edge,norient_face)
!  ...determine nodes coordinates
      call nodcor(Mdle, xnod)
!  ...get the element boundary conditions flags
      call find_bc(Mdle, ibc)
!  ...get current solution dofs
      call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
!
!  ...clear space for auxiliary matrices
      bload_E(:) = ZERO;   gramP(:) = ZERO;
      bload_Ec(:) = ZERO;  bload_Imp(:) = ZERO
!
!--------------------------------------------------------------------------
!
!              E L E M E N T   I N T E G R A L S
!
!--------------------------------------------------------------------------
!
!  ...use the enriched order to set the quadrature
      INTEGRATION = NORD_ADD + 1
      call set_3D_int_DPG(ntype,norder,norient_face, nint,xiloc,waloc)
      INTEGRATION = 0
!
!  ...loop over
      do l=1,nint
         xi(1:3) = xiloc(1:3,l)
         wa = waloc(l)
!
!     ...determine element H1 shape functions
         call shape3DH(ntype,xi,norder,norient_edge,norient_face,  &
                       nrdof,shapH,gradH)
#if DEBUG_MODE
         if (nrdof .ne. NrdofH) then
            write(*,*) 'elem_residual_maxwell: INCONSISTENCY NrdofH. stop.'
            stop
         endif
#endif
!   ...determine element H(curl) shape functions
         call shape3DE(ntype,xi,norder,norient_edge,norient_face, &
                       nrdof,shapE,curlE)
#if DEBUG_MODE
         if (nrdof .ne. NrdofE) then
            write(*,*) 'elem_residual_maxwell: INCONSISTENCY NrdofE. stop.'
            stop
         endif
#endif
!     ...determine element L2 shape functions
         call shape3DQ(ntype,xi,norder, nrdof,shapQ)
#if DEBUG_MODE
         if (nrdof .ne. NrdofQ) then
            write(*,*) 'elem_residual_maxwell: INCONSISTENCY NrdofQ. stop.'
            stop
         endif
#endif
!     ...determine discontinuous H(curl) shape functions
         call shape3EE(ntype,xi,nordP, nrdof,shapEE,curlEE)
#if DEBUG_MODE
         if (nrdof .ne. NrdofEE) then
            write(*,*) 'elem_residual_maxwell: INCONSISTENCY NrdofEE. stop.'
            stop
         endif
#endif
!     ...geometry
         call geom3D(Mdle,xi,xnod,shapH,gradH,NrdofH, &
                      x,dxdxi,dxidx,rjac,iflag)
!
!     ...get permittivity tensor
         call get_permittivity(mdle,x, eps)
!
!     ...integration weight
         weight = rjac*wa
!     ...compute the approximate solution
         zsolQ = ZERO
!
         do k=1,NrdofQ
            zsolQ(1:6)  = zsolQ(1:6)  + zdofQ(1:6,k)*shapQ(k)
         enddo
         zsolQ = zsolQ/rjac
!
!     ...get the RHS
!     ...zfval (heat eqn rhs), zJ (maxwell rhs)
         call getf(Mdle,x, zfval,zJ)
!
!     ...set auxiliary constants
         za = (ZI*OMEGA*EPSILON)*eps(:,:)
         zc = (ZI*OMEGA*MU)*IDENTITY
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
         do k1=1,NrdofEE
!
!        ...pick up pulled-back shape functions
            fldF(:) = shapF(:,k1);  crlF(:) = curlF(:,k1)
            fldG(:) = fldF(:);      crlG(:) = crlF(:)
            epsTfldF(:) = epsTshapF(:,k1);
!
!  --- Residual ---
!
!        ...accumulate for the load
!        ...first eqn
!           (J^imp,F)
            k = 2*k1-1
            zaux = zJ(1)*fldF(1) + zJ(2)*fldF(2) + zJ(3)*fldF(3)
            bload_E(k) = bload_E(k) + zaux*weight
!
!        ...first eqn
            k = 2*k1-1
!           - ( (H,curl(F)) - iωε(E,F) )
            zaux = crlF(1)*zsolQ(4) + crlF(2)*zsolQ(5) + crlF(3)*zsolQ(6)
            zcux = conjg(epsTfldF(1))*zsolQ(1) + &
                   conjg(epsTfldF(2))*zsolQ(2) + &
                   conjg(epsTfldF(3))*zsolQ(3)
            bload_E(k) = bload_E(k) - (zaux - zcux)*weight
!
!        ...second eqn
            k = 2*k1
!           - ( (E,curl(G)) + iωμ(H,G) )
            zaux = crlG(1)*zsolQ(1) + crlG(2)*zsolQ(2) + crlG(3)*zsolQ(3)
            zcux = zc1*(fldG(1)*zsolQ(4) + &
                        fldG(2)*zsolQ(5) + &
                        fldG(3)*zsolQ(6) )
            bload_E(k) = bload_E(k) - (zaux + zcux)*weight
!
!  --- Gram matrix ---
!
!        ...loop through enriched H(curl) test functions
            do k2=k1,NrdofEE
               fldE(:) = shapF(:,k2);  crlE(:) = curlF(:,k2)
               epsTfldE(:) = epsTshapF(:,k2);  epscrlE(:) = epscurlF(:,k2);
!
               call dot_product(fldF,fldE, FF)
               call dot_product(crlF,crlE, CC)
!
!           ...accumulate for the Hermitian Gram matrix
!              (compute upper triangular only)
!           ...testNorm = Scaled Adjoint Graph norm
!                 ||v|| = alpha*(v,v) + (A^* v, A^* v)
!              (first eqn multiplied by F, second eqn by G)
!              G_ij=(phi_j,phi_i)_testNorm is 2x2 matrix
!              where (phi_j,phi_i)_l2Norm = Int[phi_i^* phi_j]
!              -------------------------
!              | (F_i,F_j)   (F_i,G_j) |
!              | (G_i,F_j)   (G_i,G_j) |
!              -------------------------
!              (F_j,F_i) terms = Int[F_^*i F_j] terms (G_11)
               n = 2*k1-1; m = 2*k2-1
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
         enddo
      enddo
!
!--------------------------------------------------------------------------
!
!              B O U N D A R Y      I N T E G R A L S
!
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
!        ...determine discontinuous H(curl) shape functions
            call shape3EE(ntype,xi,nordP, NrdofEE,shapEE,curlEE)
!
!        ...determine element H1 shape functions (for geometry)
            call shape3DH(ntype,xi,norder,norient_edge,norient_face, &
                          nrdof,shapH,gradH)
#if DEBUG_MODE
            if (nrdof .ne. NrdofH) then
               write(*,*) 'elem_residual_maxwell: INCONSISTENCY NrdofH. stop.'
               stop
            endif
#endif
!
!        ...determine element H(curl) shape functions (for fluxes)
            call shape3DE(ntype,xi,norder,norient_edge,norient_face, &
                          nrdof,shapE,curlE)
#if DEBUG_MODE
            if (nrdof .ne. NrdofE) then
               write(*,*) 'elem_residual_maxwell: INCONSISTENCY NrdofE. stop.'
               stop
            endif
#endif
!
!        ...geometry
            call bgeom3D(Mdle,xi,xnod,shapH,gradH,NrdofH,dxidt,nsign, &
                         x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
            weight = bjac*wtloc(l)
!
!        ...compute approximate fluxes at the point
            zsolExi = ZERO
!
            do ivar=1,2
               do k=1,NrdofE
                  zsolExi(1:3,ivar) = zsolExi(1:3,ivar) &
                                    + zdofE(ivar,k)*shapE(1:3,k)
               enddo
               zsolE(1:3,ivar) = zsolExi(1,ivar)*dxidx(1,1:3) &
                               + zsolExi(2,ivar)*dxidx(2,1:3) &
                               + zsolExi(3,ivar)*dxidx(3,1:3)
               call zcross_product(rn,zsolE(1:3,ivar), zflux (1:3,ivar))
               call zcross_product(rn,zflux(1:3,ivar), zflux2(1:3,ivar))
            enddo
!
!        ...check for impedance BC (elimination strategy)
            if (ibc(ifc,2).eq.3) then
!           ...impedance surface load [zImp should be zero here]
               call get_bdSource(Mdle,x,rn, zImp)
               zflux2(1:3,1) = GAMMA*zflux2(1:3,1) + zImp
            endif
!
!        ...loop through enriched test functions
            do k1=1,NrdofEE
               E1(1:3) = shapEE(1,k1)*dxidx(1,1:3) &
                       + shapEE(2,k1)*dxidx(2,1:3) &
                       + shapEE(3,k1)*dxidx(3,1:3)
!
               k=2*k1-1
!           ...check for impedance BC (elimination strategy)
               if (ibc(ifc,2).eq.3) then
!              - GAMMA * < n x n x E , G >
                  zaux = E1(1)*zflux2(1,1) + E1(2)*zflux2(2,1) + E1(3)*zflux2(3,1)
                  bload_E(k) = bload_E(k) - zaux * weight
               else
!              - <n x H, F>
                  zaux = E1(1)*zflux(1,2) + E1(2)*zflux(2,2) + E1(3)*zflux(3,2)
                  bload_E(k) = bload_E(k) - zaux * weight
               endif
!              - <n x E, G>
               k = 2*k1
               zaux = E1(1)*zflux(1,1) + E1(2)*zflux(2,1) + E1(3)*zflux(3,1)
               bload_E(k) = bload_E(k) - zaux * weight
            enddo
!
!        ...check for impedance BC (L2 penalty method)
            if (ibc(ifc,2).ne.2) cycle
!
!        ...impedance surface load [zImp should be zero here]
            call get_bdSource(Mdle,x,rn, zImp)
!        ...define the weight of the penalty term
            pen = 1.d0
!        ...compute residual contribution from impedance boundary
            do k1=1,NrdofE
               E1(1:3) = shapE(1,k1)*dxidx(1,1:3) &
                       + shapE(2,k1)*dxidx(2,1:3) &
                       + shapE(3,k1)*dxidx(3,1:3)
!
               call cross_product(rn,E1, rntimesE)
               call cross_product(rn,rntimesE, rn2timesE)
!
!           ...1st test function
               k = 2*k1-1
!              + GAMMA^2 * < n x n x E , n x n x F >
               bload_Imp(k) = bload_Imp(k) + (rn2timesE(1)*zflux2(1,1)  &
                                           +  rn2timesE(2)*zflux2(2,1)  &
                                           +  rn2timesE(3)*zflux2(3,1)  &
                                             )*GAMMA*GAMMA*weight/pen
!              - GAMMA * < n x H, n x n x F >
               bload_Imp(k) = bload_Imp(k) - (rn2timesE(1)*zflux(1,2)   &
                                           +  rn2timesE(2)*zflux(2,2)   &
                                           +  rn2timesE(3)*zflux(3,2)   &
                                             )*GAMMA*weight/pen
!              + GAMMA * < zImp , n x n x F >
               bload_Imp(k) = bload_Imp(k) + (rn2timesE(1)*zImp(1)  &
                                           +  rn2timesE(2)*zImp(2)  &
                                           +  rn2timesE(3)*zImp(3)  &
                                             )*GAMMA*weight/pen
!           ...2nd test function
               k = 2*k1
!              - GAMMA * < n x n x E , n x G >
               bload_Imp(k) = bload_Imp(k) - (rntimesE(1)*zflux2(1,1)  &
                                           +  rntimesE(2)*zflux2(2,1)  &
                                           +  rntimesE(3)*zflux2(3,1)  &
                                             )*GAMMA*weight/pen
!              + < n x H , n x G >
               bload_Imp(k) = bload_Imp(k) + (rntimesE(1)*zflux(1,2)   &
                                           +  rntimesE(2)*zflux(2,2)   &
                                           +  rntimesE(3)*zflux(3,2)   &
                                             )*weight/pen
!              - < zImp , n x G >
               bload_Imp(k) = bload_Imp(k) - (rntimesE(1)*zImp(1)  &
                                           +  rntimesE(2)*zImp(2)  &
                                           +  rntimesE(3)*zImp(3)  &
                                             )*weight/pen
            enddo
!
         enddo
      enddo
!
#if DEBUG_MODE
      if (iprint.gt.0) then
         write(*,7015) bload_E(1:2*NrdofEE)
 7015    format('elem_residual_maxwell: FINAL bload_E = ',10(/,6(2e12.5,2x)))
         call pause
      endif
#endif
!
!--------------------------------------------------------------------------
!
!  ...factorize the test Gram matrix
      call ZPPTRF('U', NrTest, gramP, info)
      if (info.ne.0) then
         write(*,*) 'elem_residual_maxwell: ZPPTRF: Mdle,info = ',Mdle,info,'. stop.'
         stop
      endif
!
!  ...save copies of the RHS to compute later the residual
      bload_Ec = bload_E
!
!  ...compute the product of inverted test Gram matrix with RHS,
!  ...bload_E is overwritten with the solution
      call ZPPTRS('U', NrTest, 1, gramP, bload_E, NrTest, info)
      if (info.ne.0) then
         write(*,*) 'elem_residual_maxwell: ZPPTRS: Mdle,info = ',Mdle,info,'. stop.'
         stop
      endif
!
      deallocate(gramP)
!
!  ...compute the residual
      zresid = ZERO
      do k=1,NrTest
         zresid = zresid + bload_Ec(k)*conjg(bload_E(k))
      enddo
!
!  ...account for impedance BC penalty term (L2 penalty method)
!     test norm for the residual has then two separate contributions:
!     1) Usual DPG residual measured in adjoint test norm ||\psi||_V
!     2) Additional Impedance BC residual measured in L2 norm ||\phi||
      if (IBCFLAG.eq.2) then
         do k=1,2*NrdofE
            zresid = zresid + bload_Imp(k)*conjg(bload_Imp(k))
         enddo
      endif
!
      Resid = real(zresid,8)
!
!  ...set suggested refinement flag
      select case(ntype)
         case(MDLB);      Nref_flag = 111
         case(MDLP);      Nref_flag = 11
         case(MDLN,MDLD); Nref_flag = 1
      end select
!
#if DEBUG_MODE
      if (iprint.eq.1) then
         write(*,7010) Mdle, Resid
 7010    format('elem_residual_maxwell: Mdle, Resid = ',i5,3x,e12.5)
         call pause
      endif
#endif
!
   end subroutine elem_residual_maxwell

