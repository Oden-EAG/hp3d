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
!> @date       March 2025
!------------------------------------------------------------------------------
   subroutine elem_residual_maxwell_opt(Mdle,NrTest,                    &
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
      integer :: norderf(5)
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
      real(8) :: xi(3), x(3), rn(3)
      real(8) :: dxidt(3,2), dxdt(3,2)
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
!
!  ...L2 shape functions
      real(8), dimension(MAXbrickQ) :: shapQ
!
!  ...enriched Hcurl shape functions
      real(8), dimension(3,MAXbrickEE) :: shapEE
      real(8), dimension(3,MAXbrickEE) :: curlEE
!
!   ..Gram matrix in full format (although we'll use upper half only)
      real(8),    allocatable :: gram_r(:,:)
      complex(8), allocatable :: gram_FF(:,:),gram_FG(:,:),gram_GG(:,:)
      complex(8), allocatable :: gram(:,:)
!
!  ...intermediate values
      real(8) :: FF, CC, afac
      real(8) :: fldE(3), fldH(3), crlE(3), crlH(3), rKTF(3)
      real(8) :: fldF(3), fldG(3), crlF(3), crlG(3)
      complex(8) :: AstarF1(3,NrdofEE),AstarF2(3,NrdofEE)
      complex(8) :: AstarG1(3,NrdofEE),AstarG2(3,NrdofEE)
!   ..dynamic allocation arrays for filling values at quadrature points
      real(8),    allocatable :: all_shapF(:,:),all_curlF(:,:)
      complex(8), allocatable :: all_AstarF1(:,:),all_AstarF2(:,:),   &
                                 all_AstarG1(:,:),all_AstarG2(:,:)
!
!  ...load vector for the enriched space
      VTYPE, dimension(NrTest)   :: bload_E
      VTYPE, dimension(2*NrdofE) :: bload_Imp
!
!  ...quadrature data
      real(8) :: xiloc(3,MAXNINT3ADD), waloc(MAXNINT3ADD)
      real(8) :: tloc(2,MAXNINT2ADD), wtloc(MAXNINT2ADD)
      real(8) :: weight,wa,sqrt_weight
!
!  ...BC's flags
      integer :: ibc(6,NRINDEX)
!
!  ...Maxwell load and auxiliary variables
      complex(8) :: zJ(3), zImp(3) , zL(3)
      real(8)    :: E1(3), rntimesE(3), rn2timesE(3)
      real(8)    :: pen
!
!  ...approximate solution
      VTYPE, dimension(3,2) :: zsolExi,zsolE,zflux,zflux2
      VTYPE, dimension(6)   :: zsolQ
!
!  ...auxiliary
      VTYPE :: zresid, zaux, zcux
!
!  ...number of faces per element type
      integer :: nrf
!
!  ...various variables for the problem
      integer :: k1,k2,m,n,nint,k,l,ivar,iflag
      integer :: nordP,nsign,ifc,info,nrdof

      integer :: nda, loff, koff, joff, ioff, nE, nQ, i , j
      VTYPE   :: za(3,3),zc(3,3),zc1
!
      integer, external :: ij_upper_to_packed
!
#if HP3D_DEBUG
      integer :: iprint
      iprint = 0
#endif
!
!--------------------------------------------------------------------------
!
      select case(TEST_NORM)
         case(GRAPH_NORM); afac = ALPHA_NORM
         case(GRAPH_DIAG); afac = ALPHA_NORM
         case( MATH_NORM); afac = 1.d0
         case default
            write(*,*) 'elem_residual_maxwell: invalid test norm!'
            stop
      end select
!
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
      bload_E(:) = ZERO
      bload_Imp(:) = ZERO

!   ..allocate matrices
      nE =   NrdofEE
      nQ = 3*NrdofQ
      allocate(gram(NrTest,NrTest))
!
!--------------------------------------------------------------------------
!
!              E L E M E N T   I N T E G R A L S
!
!--------------------------------------------------------------------------
!
!  ...use the enriched order to set the quadrature
      INTEGRATION = NORD_ADD
      call set_3D_int_DPG(ntype,norder,norient_face, nint,xiloc,waloc)
      INTEGRATION = 0
!      
!   ..we allocate the auxiliary arrays with this dimension
      nda = 3*nint
      allocate(all_shapF(nE,nda),all_curlF(nE,nda))
      allocate(all_AstarF1(nE,nda),all_AstarF2(nE,nda))
      allocate(all_AstarG1(nE,nda),all_AstarG2(nE,nda))
!
!  ...loop over
      do l=1,nint
         xi(1:3) = xiloc(1:3,l)
         wa = waloc(l)
!
!     ...determine element H1 shape functions
         call shape3DH(ntype,xi,norder,norient_edge,norient_face,  &
                       nrdof,shapH,gradH)
#if HP3D_DEBUG
         if (nrdof .ne. NrdofH) then
            write(*,*) 'elem_residual_maxwell: INCONSISTENCY NrdofH. stop.'
            stop
         endif
#endif
!   ...determine element H(curl) shape functions
         call shape3DE(ntype,xi,norder,norient_edge,norient_face, &
                       nrdof,shapE,curlE)
#if HP3D_DEBUG
         if (nrdof .ne. NrdofE) then
            write(*,*) 'elem_residual_maxwell: INCONSISTENCY NrdofE. stop.'
            stop
         endif
#endif
!     ...determine element L2 shape functions
         call shape3DQ(ntype,xi,norder, nrdof,shapQ)
#if HP3D_DEBUG
         if (nrdof .ne. NrdofQ) then
            write(*,*) 'elem_residual_maxwell: INCONSISTENCY NrdofQ. stop.'
            stop
         endif
#endif
!     ...determine discontinuous H(curl) shape functions
         call shape3EE(ntype,xi,nordP, nrdof,shapEE,curlEE)
#if HP3D_DEBUG
         if (nrdof .ne. NrdofEE) then
            write(*,*) 'elem_residual_maxwell: INCONSISTENCY NrdofEE. stop.'
            stop
         endif
#endif
!     ...geometry
         call geom3D(Mdle,xi,xnod,shapH,gradH,NrdofH, &
                      x,dxdxi,dxidx,rjac,iflag)
!
!     ...integration weight
         weight = rjac*wa
         sqrt_weight = sqrt(weight)
!
!     ...compute the approximate solution
         zsolQ = ZERO
!
         do k=1,NrdofQ
            zsolQ(1:6)  = zsolQ(1:6)  + zdofQ(1:6,k)*shapQ(k)
         enddo
         zsolQ = zsolQ/rjac
!
!     ...get the RHS
!     ...zJ (maxwell rhs)
         call getf(Mdle,x, zJ,zL)
!
!         
!     ...STORE ALL THE VALUES OF ADJOINT AT INTEGRATION POINTS
!     ...offset in auxiliary matrices
         loff = 3*(l-1)
!     ...apply pullbacks
         call DGEMM('T','N',3,nE,3,1.d0     ,dxidx,3,shapEE,3,0.d0,shapF,3)!
!     ...save in the big matrix of size nE by 3*nint
         all_shapF(1:nE,loff+1:loff+3) = sqrt_weight*transpose(shapF)
!   
         call DGEMM('N','N',3,nE,3,1.d0/rjac,dxdxi,3,curlEE,3,0.d0,curlF,3)
!     ...save in the big matrix of size nE by 3*nint
         all_curlF(1:nE,loff+1:loff+3) = sqrt_weight*transpose(curlF)
!   
!     ...Evaluate A^* on (F ; 0) and on (0 ; G)
         call get_Astar_multiple(Mdle,x,nE,shapF,curlF,AstarF1,AstarF2,AstarG1,AstarG2)
!     ...save in big array of size nE by 3*nint; apply conjugate now for convenience
         all_AstarF1(1:nE,loff+1:loff+3) = conjg(sqrt_weight*transpose(AstarF1))
         all_AstarF2(1:nE,loff+1:loff+3) = conjg(sqrt_weight*transpose(AstarF2))
         all_AstarG1(1:nE,loff+1:loff+3) = conjg(sqrt_weight*transpose(AstarG1))
         all_AstarG2(1:nE,loff+1:loff+3) = conjg(sqrt_weight*transpose(AstarG2))
!
!     ...loop through enriched H(curl) test functions
         do k1=1,NrdofEE
!
!        ...pick up pulled-back shape functions
            fldF(:) = shapF(:,k1);  crlF(:) = curlF(:,k1)
            fldG(:) = fldF(:);      crlG(:) = crlF(:)
!
!  --- Residual ---
!
!        ...accumulate for the load
!           RHS:
!           (J^imp,F) first  equation RHS (with first H(curl) test function F)
!           (L^fdy,G) second equation RHS (Faraday's eqn artificial load)
            n = 2*k1-1
            bload_E(n) = bload_E(n)                                   &
                       + (fldF(1)*zJ(1)+fldF(2)*zJ(2)+fldF(3)*zJ(3))  &
                       * weight
!
            n = 2*k1
            bload_E(n) = bload_E(n)                                   &
                       + (fldG(1)*zL(1)+fldG(2)*zL(2)+fldG(3)*zL(3))  &
                       * weight
!
!   1ST OPTION TO INTEGRATE: evaluate adjoint operator
!        ...First equation. Test function F
            n = 2*k1-1
!        ...Accumulate    -( E , AstarF1 )   -( H , AstarF2 )
            bload_E(n) = bload_E(n)                                  &
                       - ( SUM(zsolQ(1:3)*conjg(AstarF1(:,k1)))         &
                          +SUM(zsolQ(4:6)*conjg(AstarF2(:,k1))))*weight
!
!        ...Second equation. Test function G
            n = 2*k1
!        ...Accumulate    -( E , AstarG1 )   -( H , AstarG2 )
            bload_E(n) = bload_E(n)                                  &
                       - ( SUM(zsolQ(1:3)*conjg(AstarG1(:,k1)))         &
                          +SUM(zsolQ(4:6)*conjg(AstarG2(:,k1))))*weight
!
!        ...end of loop through enriched H(curl) test functions
         enddo
!   
!   ..end of loop through integration points
      enddo
!
!
!     --- Gram matrix: assemble and factorize ---
!   
      allocate(gram_r(nE,nE))
      gram_r(:,:) = 0.d0
!   
      select case(TEST_NORM)
         case( MATH_NORM)
            afac = 1.d0
!   
!        ...α*(F,F)
!        ...α*(G,G)
            call DSYRK('U','N',nE,nda,afac,all_shapF,nE,0.d0,gram_r,nE) ! real-valued
!        ...(curl F, curl F)
!        ...(curl G, curl G)
            call DSYRK('U','N',nE,nda,1.d0,all_curlF,nE,1.d0,gram_r,nE) ! real-valued
!        ...Cholesky factorization of the real symmetric matrix Gram_r 
            call DPOTRF('U',nE,gram_r,nE,info)
            if (info.ne.0) then
               write(*,*) 'elem_maxwell_opt: DPOTRF: Mdle,info = ',Mdle,info,'. stop.'
               stop
            endif
            do j = 1,nE
               joff = 2*(j-1)
               do i = 1,j
                  ioff = 2*(i-1)
!              ...Assemble            
!                 (F,F) + (CF,CF)         0 
!                       0           (G,G) + (CG,CG)
                  gram(ioff+1,joff+1) = cmplx(gram_r(i,j),0.d0,8)
                  gram(ioff+2,joff+2) = cmplx(gram_r(i,j),0.d0,8)
               enddo
            enddo
            deallocate(gram_r,all_shapF,all_curlF)
      
         case(GRAPH_NORM)
            afac = ALPHA_NORM
!     
!           ...α*(F,F)
!           ...α*(G,G)
            call DSYRK('U','N',nE,nda,afac,all_shapF,nE,0.d0,gram_r,nE) ! real-valued
            do j = 1,nE
               joff = 2*(j-1)
               do i = 1,j
                  ioff = 2*(i-1)
!              ...Accumulate into the complex-valued matrix
!                 α*(F,F)     0 
!                    0     α*(G,G)
                  gram(ioff+1,joff+1) = cmplx(gram_r(i,j),0.d0,8)
                  gram(ioff+2,joff+2) = cmplx(gram_r(i,j),0.d0,8)
               enddo
            enddo
            deallocate(gram_r,all_shapF,all_curlF)
!            
            allocate(gram_FF(nE,nE),gram_GG(nE,nE),gram_FG(nE,nE))
            gram_FF(:,:) = ZERO; gram_GG(:,:) = ZERO; gram_FG(:,:) = ZERO
!     
!        ...accumulate (AstarF1,AstarF1)
            call ZHERK('U','N',nE,nda,ZONE,all_AstarF1,nE,  &
                       ZERO,gram_FF(1:nE,1:nE ),nE)
!        ...accumulate (AstarF2,AstarF2)
            call ZHERK('U','N',nE,nda,ZONE,all_AstarF2,nE,  &
                       ZONE,gram_FF(1:nE,1:nE ),nE)
!     
!        ...accumulate (AstarG1,AstarG1)
            call ZHERK('U','N',nE,nda,ZONE,all_AstarG1,nE,  &
                       ZERO,gram_GG(1:nE,1:nE ),nE)
!        ...accumulate (AstarG2,AstarG2)
            call ZHERK('U','N',nE,nda,ZONE,all_AstarG2,nE,  &
                       ZONE,gram_GG(1:nE,1:nE ),nE)
!     
!        ...accumulate (AstarF1,AstarG1)
            call ZGEMM('N','C',nE,nE,nda,ZONE,all_AstarF1,nE,all_AstarG1,nE,  &
                       ZERO,gram_FG(1:nE,1:nE ),nE)
!        ...accumulate (AstarF2,AstarG2)
            call ZGEMM('N','C',nE,nE,nda,ZONE,all_AstarF2,nE,all_AstarG2,nE,  &
                       ZONE,gram_FG(1:nE,1:nE ),nE)
!     
!        ...organize blocks within upper part of gram
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
            call shape3EE(ntype,xi,nordP, nrdof,shapEE,curlEE)
!
!        ...determine element H1 shape functions (for geometry)
            call shape3DH(ntype,xi,norder,norient_edge,norient_face, &
                          nrdof,shapH,gradH)
#if HP3D_DEBUG
            if (nrdof .ne. NrdofH) then
               write(*,*) 'elem_residual_maxwell: INCONSISTENCY NrdofH. stop.'
               stop
            endif
#endif
!
!        ...determine element H(curl) shape functions (for fluxes)
            call shape3DE(ntype,xi,norder,norient_edge,norient_face, &
                          nrdof,shapE,curlE)
#if HP3D_DEBUG
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
#if HP3D_DEBUG
      if (iprint.gt.0) then
         write(*,7015) bload_E(1:2*NrdofEE)
 7015    format('elem_residual_maxwell: FINAL bload_E = ',10(/,6(2e12.5,2x)))
         call pause
      endif
#endif
!
!--------------------------------------------------------------------------
!   ..Solve triangular system to obtain R~, (LX=) U^*X = [l]
      call ZTRTRS('U','C','N',NrTest,1,gram,NrTest,bload_E,NrTest,info)
      if (info.ne.0) then
         write(*,*) 'elem_residual_maxwell_opt: ZTPTRS: Mdle,info = ',Mdle,info,'. stop.'
         stop
      endif
! !..C. Matrix multiply: B^* G^-1 B (=B~^* B~)
!    call ZHERK('U','C',1,NrTest,ZONE,bload_E,NrTest,ZERO,zBDPG,NrTrial+1)
! !
!    deallocate(stiff_ALL)
!
      deallocate(gram)
!
!  ...compute the residual
      zresid = ZERO
      do k=1,NrTest
         zresid = zresid + bload_E(k)*conjg(bload_E(k))
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
#if HP3D_DEBUG
      if (iprint.eq.1) then
         write(*,7010) Mdle, Resid
 7010    format('elem_residual_maxwell: Mdle, Resid = ',i5,3x,e12.5)
         call pause
      endif
#endif
!
   end subroutine elem_residual_maxwell_opt

