!-------------------------------------------------------------------------------
!
!    routine name      - elem_maxwell
!
!-------------------------------------------------------------------------------
!
!     latest revision:  - June 2021
!
!> @brief         - routine returns unconstrained (ordinary)
!                         stiffness matrix and load vector for the
!                         UW DPG formulation for Maxwell equations
!                       - Uses sum factorization for fast integration
!                         with swapped loops
!
!     arguments:
!        in:
!           Mdle        - element middle node number
!           Fld_flag    - field flag (0: pump, 1: signal)
!           NrTest      - total number of test dof
!           NrTrial     - total number of trial dof
!           NrdofEE     - number of H(curl) test dof
!           NrdofH      - number of H1 trial dof
!           NrdofE      - number of H(curl) trial dof
!           NrdofQ      - number of L2 trial dof
!           NrdofEi     - number of H(curl) trial interface dof
!           MdE         - num rows of ZalocEE,ZalocEQ
!           MdQ         - num rows of ZalocQE,ZalocQQ
!        out:
!           ZblocE      - load vectors
!           ZblocQ
!           ZalocEE     - stiffness matrices
!           ZalocEQ
!           ZalocQE
!           ZalocQQ
!
!-------------------------------------------------------------------------------
!
#include "typedefs.h"
!
subroutine elem_maxwell(Mdle,Fld_flag,                &
                        NrTest,NrTrial,               &
                        NrdofEE,                      &
                        NrdofH,NrdofE,NrdofQ,         &
                        NrdofEi,                      &
                        MdE,MdQ,                      &
                        ZblocE,ZalocEE,ZalocEQ,       &
                        ZblocQ,ZalocQE,ZalocQQ)
!..modules used
   use control
   use parametersDPG
   use data_structure3D
   use laserParam
   use commonParam
   use mpi_wrapper
!..no implicit statements
   implicit none
!..declare input/output variables
   integer,                   intent(in)  :: Mdle
   integer,                   intent(in)  :: Fld_flag
   integer,                   intent(in)  :: NrTest
   integer,                   intent(in)  :: NrTrial
   integer,                   intent(in)  :: NrdofEE
   integer,                   intent(in)  :: NrdofH
   integer,                   intent(in)  :: NrdofE
   integer,                   intent(in)  :: NrdofQ
   integer,                   intent(in)  :: NrdofEi
   integer,                   intent(in)  :: MdE
   integer,                   intent(in)  :: MdQ
   VTYPE, dimension(MdE),     intent(out) :: ZblocE
   VTYPE, dimension(MdE,MdE), intent(out) :: ZalocEE
   VTYPE, dimension(MdE,MdQ), intent(out) :: ZalocEQ
   VTYPE, dimension(MdQ),     intent(out) :: ZblocQ
   VTYPE, dimension(MdQ,MdE), intent(out) :: ZalocQE
   VTYPE, dimension(MdQ,MdQ), intent(out) :: ZalocQQ
!
!..declare edge/face type variables
   integer :: etype,ftype
!
!..declare element order, orientation for edges and faces
   integer, dimension(19)    :: norder
   integer, dimension(12)    :: norient_edge
   integer, dimension(6)     :: norient_face
!
!..element nodes order (trial) for interfaces
   integer, dimension(19)    :: norderi
!
!..face order
   integer, dimension(5)     :: norderf
!
!..geometry dof (work space for nodcor)
   real(8), dimension(3,MAXbrickH) :: xnod
!
!..solution dof (work space for solelm)
   VTYPE, dimension(MAXEQNH,MAXbrickH) :: zdofH
   VTYPE, dimension(MAXEQNE,MAXbrickE) :: zdofE
   VTYPE, dimension(MAXEQNV,MAXbrickV) :: zdofV
   VTYPE, dimension(MAXEQNQ,MAXbrickQ) :: zdofQ
!
!..approximate solution -- using soleval
   integer :: nflag
   VTYPE, dimension(  MAXEQNH  ) ::  zsolH_soleval
   VTYPE, dimension(  MAXEQNH,3) :: zdsolH_soleval
   VTYPE, dimension(3,MAXEQNE  ) ::  zsolE_soleval
   VTYPE, dimension(3,MAXEQNE  ) :: zcurlE_soleval
   VTYPE, dimension(3,MAXEQNV  ) ::  zsolV_soleval
   VTYPE, dimension(  MAXEQNV  ) ::  zdivV_soleval
   VTYPE, dimension(  MAXEQNQ  ) ::  zsolQ_soleval
   real(8) :: rsolH
!
!..variables for geometry
   real(8), dimension(3)    :: xi,x,rn,daux
   real(8), dimension(3,2)  :: dxidt,dxdt
   real(8), dimension(3,3)  :: dxdxi,dxidx
   real(8), dimension(2)    :: t
!
!..H1 shape functions
   real(8), dimension(MAXbrickH)    :: shapH
   real(8), dimension(3,MAXbrickH)  :: gradH
!
!..H(curl) shape functions
   real(8), dimension(3,MAXbrickE)  :: shapE
   real(8), dimension(3,MAXbrickE)  :: curlE
!
!..L2 shape functions
   real(8), dimension(MAXbrickQ)    :: shapQ
!
!..Enriched H1 shape functions
   real(8) , dimension(3,MAXbrickEE) :: shapEE
   real(8) , dimension(3,MAXbrickEE) :: curlEE
!
!..nrdof for interface only (without bubbles)
   !integer :: ik,nrdofEEi
!
!..H(curl) bubble index
   !integer, allocatable :: idxEE(:)
!
!..load vector for the enriched space
   VTYPE :: bload_E(NrTest)
!
!..Gram matrix in packed format
   !VTYPE :: gramP(NrTest*(NrTest+1)/2)
   VTYPE, allocatable :: gramP(:)
   real(8)  :: FF, CF, FC
   real(8)  :: fldE(3), fldH(3), crlE(3), rotE(3)
   real(8)  :: fldF(3), fldG(3), crlF(3), crlG(3), rotF(3), rotG(3)
!
!..matrices for transpose filling (swapped loops)
!..stiffness matrices (transposed) for the enriched test space
   complex(8), allocatable :: stiff_EE_T(:,:),stiff_EQ_T(:,:)
   complex(8), allocatable :: stiff_ALL(:,:),zaloc(:,:)
!
!..3D quadrature data
   real(8), dimension(3,MAXNINT3ADD)  :: xiloc
   real(8), dimension(MAXNINT3ADD)    :: waloc
!
!..2D quadrature data
   real(8), dimension(2,MAXNINT2ADD)  :: tloc
   real(8), dimension(MAXNINT2ADD)    :: wtloc
!
!..BC's flags
   integer, dimension(6,NRINDEX)      :: ibc
!
!..for auxiliary computation
   VTYPE :: zaux,zbux,zcux
!
!..Maxwell load and auxiliary variables
   VTYPE  , dimension(3) :: zJ,zImp
   real(8), dimension(3) :: E1,E2,rntimesE,rn2timesE
!
!..number of edge,faces per element type
   integer :: nre, nrf
!
!..various variables for the problem
   real(8) :: rjac,weight,wa,CC
   real(8) :: bjac,minz,maxz,elem_z
   integer :: i1,j1,j2,k1,k2,i,k,l,nint,n,m
   integer :: iflag,info
   integer :: nrdof,nordP,nsign,ifc,ndom
   VTYPE   :: zfval
   VTYPE   :: za(3,3),zc(3,3)
   VTYPE   :: zaJ(3,3),zcJ(3,3)
!
!..for polarizations function
   VTYPE, dimension(3,3) :: bg_pol,gain_pol,raman_pol
   real(8) :: delta_n
   integer :: dom_flag
!
!..OMEGA_RATIO_SIGNAL or OMEGA_RATIO_PUMP
   real(8) :: OMEGA_RATIO_FLD
!..WAVENUM_SIGNAL or WAVENUM_PUMP
   real(8) :: WAVENUM_FLD
!
!..for PML
   VTYPE :: zbeta,zdbeta,zd2beta,detJstretch
   VTYPE, dimension(3,3) :: Jstretch,invJstretch,JJstretch
!
   integer, external :: ij_upper_to_packed
!
!..timer
!   real(8) :: start_time,end_time
!
#if HP3D_DEBUG
   integer :: iprint
   iprint = 0
#endif
!
!-------------------------------------------------------------------------------
!
#if HP3D_DEBUG
   if (iprint.eq.1) then
      write(*,*) 'elem_maxwell: Mdle = ', Mdle
   endif
#endif
!
!..allocate matrices
   allocate(gramP(NrTest*(NrTest+1)/2))
   allocate(stiff_EE_T(2*NrdofEi,NrTest))
   allocate(stiff_EQ_T(6*NrdofQ ,NrTest))
!
!..element type
   etype = NODES(Mdle)%ntype
   nre = nedge(etype); nrf = nface(etype)
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
         write(*,*) 'elem_maxwell: invalid etype param. stop.'
         stop
   end select
!
!..determine edge and face orientations
   call find_orient(Mdle, norient_edge,norient_face)
!
!..determine nodes coordinates
   call nodcor(Mdle, xnod)
!
!..determine z-coordinate inside the element
   select case(etype)
      case(MDLB)
         maxz = maxval(xnod(3,1:8))
         minz = minval(xnod(3,1:8))
      case(MDLP)
         maxz = maxval(xnod(3,1:6))
         minz = minval(xnod(3,1:6))
      case default
         write(*,*) 'elem_maxwell: unexpected etype=',etype,'. stop.'
         stop
   end select
   elem_z = (minz + maxz) / 2.d0
!
!..get the element boundary conditions flags
   call find_bc(Mdle, ibc)
!
!..clear space for output matrices
   ZblocE  = ZERO; ZblocQ  = ZERO
   ZalocEE = ZERO; ZalocEQ = ZERO
   ZalocQE = ZERO; ZalocQQ = ZERO
!
!..clear space for auxiliary matrices
   bload_E    = ZERO
   gramP      = ZERO
   stiff_EE_T = ZERO
   stiff_EQ_T = ZERO
!
!..set OMEGA_RATIO_FLD and WAVENUM_FLD
   select case(Fld_flag)
      case(0)
         OMEGA_RATIO_FLD = OMEGA_RATIO_PUMP
         WAVENUM_FLD     = WAVENUM_PUMP
      case(1)
         OMEGA_RATIO_FLD = OMEGA_RATIO_SIGNAL ! 1.0d0
         WAVENUM_FLD     = WAVENUM_SIGNAL
      case default
      write(*,*) 'elem_maxwell: invalid Fld_flag param. stop.'
         stop
   end select
!
!..initialize PML matrices
   Jstretch = ZERO
   Jstretch(1,1) = ZONE
   Jstretch(2,2) = ZONE
!
   invJstretch = ZERO
   invJstretch(1,1) = ZONE
   invJstretch(2,2) = ZONE
!
   JJstretch = ZERO
   zaJ = ZERO
   zcJ = ZERO
!
!..initialize the background polarization
   call find_domain(Mdle, ndom)
!..select case of GEOM_NO to set
!..refractive index according to domain
!..Fld_flag = 1 when we are in signal element routine
   select case(GEOM_NO)
      case(1)
         bg_pol = ZERO; gain_pol = ZERO; raman_pol = ZERO
      case(4)
         bg_pol = ZERO; gain_pol = ZERO; raman_pol = ZERO
      case(5)
         gain_pol = ZERO; raman_pol = ZERO
         x(1:3) = 0.d0
         select case(ndom)
            case(1,2)
               dom_flag = 1 ! Fiber core
               call get_bgPol(dom_flag,Fld_flag,0.d0,x, bg_pol)
            case(3,4)
               dom_flag = 0 ! Fiber cladding
               call get_bgPol(dom_flag,Fld_flag,0.d0,x, bg_pol)
            case default
               write(*,*) 'elem_maxwell: unexpected ndom param. stop.'
               stop
         end select
      case default
         write(*,*) 'elem_maxwell: unexpected GEOM_NO: ', GEOM_NO, '. stop.'
         stop
!..end select case of GEOM_NO
   end select
!
   if(NONLINEAR_FLAG.eq.1) then
!  ...get current solution dofs
      call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
   endif
!
!-----------------------------------------------------------------------
!  Element Integration
!
!..use the enriched order to set the quadrature
   INTEGRATION = NORD_ADD
   call set_3D_int_DPG(etype,norder,norient_face, nint,xiloc,waloc)
   INTEGRATION = 0
!
!..start timer
!   start_time = MPI_Wtime()
!
!..loop over integration points
   do l=1,nint
!
      xi(1:3)=xiloc(1:3,l); wa=waloc(l)
!
!  ...H1 shape functions (for geometry)
      call shape3DH(etype,xi,norder,norient_edge,norient_face, nrdofH,shapH,gradH)
!
!  ...L2 shape functions for the trial space
      call shape3DQ(etype,xi,norder, nrdofQ,shapQ)
!
!  ...broken H(curl) shape functions for the enriched test space
      call shape3EE(etype,xi,nordP, nrdofEE,shapEE,curlEE)
!
!  ...geometry map
      call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH, x,dxdxi,dxidx,rjac,iflag)
#if HP3D_DEBUG
      if (iflag .ne. 0) then
         write(*,5999) Mdle,rjac
   5999  format('elem_maxwell: Negative Jacobian. Mdle,rjac=',i8,2x,e12.5)
         stop
      endif
#endif
!
!  ...integration weight
      weight = rjac*wa
!
!  ...get the RHS
      call getf(Mdle,x, zfval,zJ)
!
!  ...set auxiliary constants (updated below if nonlinear)
!     ...update bgpol depending on coordinates if refractive index is varying
!        (other than step-index)
      if (GEOM_NO .eq. 5) then
         call get_bgPol(dom_flag,Fld_flag,0.d0,x, bg_pol)
      endif
      za = (ZI*OMEGA*OMEGA_RATIO_FLD*EPSILON+SIGMA)*IDENTITY+bg_pol
      zc = (ZI*OMEGA*OMEGA_RATIO_FLD*MU)*IDENTITY
!
!  ...FOR THE NONLINEAR_FLAG = 1
      if(NONLINEAR_FLAG.eq.1) then
!     ...compute current solution using soleval
         nflag = 1
         call soleval(Mdle,xi,norient_edge,norient_face,norder,xnod, &
           zdofH,zdofE,zdofV,zdofQ,nflag,x,dxdxi,                    &
           zsolH_soleval,zdsolH_soleval,zsolE_soleval,zcurlE_soleval,&
           zsolV_soleval,zdivV_soleval,zsolQ_soleval)
!     ...initialize material refractive index perturbation
         delta_n = 0.d0
         if(HEAT_FLAG.eq.1) then
!        ...compute thermally induced perturbation to refractive index at x
            rsolH = real(zsolH_soleval(1))
            delta_n = THERMO_OPT_COEFF*rsolH
         endif
!
!     ...update background polarization with thermal perturbation
         if(delta_n .ne. 0.d0) then
            call get_bgPol(dom_flag,Fld_flag,delta_n,x, bg_pol)
         endif
!
!     ...initialize gain polarization, raman polarization
         gain_pol = ZERO; raman_pol = ZERO
!     ...skip nonlinear gain computation if inside PML region
         if ( USE_PML .and. ( (x(3).gt.PML_REGION) .or. &
                              ( (COPUMP.eq.0).and.(x(3).lt.(ZL-PML_REGION)) ) &
                            ) &
            ) goto 190
         if (ACTIVE_GAIN .gt. 0.d0) then
            if (dom_flag .eq. 1) then
               call get_activePol(zsolQ_soleval(1:12),Fld_flag,delta_n,elem_z, gain_pol)
            endif
         endif
         if (RAMAN_GAIN .gt. 0.d0) then
!        ...Fld_flag = 1 when in signal element routine
            if(Fld_flag .eq. 1) then
               call get_ramanPol(zsolQ_soleval(7:9),zsolQ_soleval(10:12), &
                              dom_flag,Fld_flag,delta_n, raman_pol)
            else
               call get_ramanPol(zsolQ_soleval(1:3),zsolQ_soleval(4:6), &
                              dom_flag,Fld_flag,delta_n, raman_pol)
            endif
         endif
 190     continue
!     ...update auxiliary constant za
         za = (ZI*OMEGA*OMEGA_RATIO_FLD*EPSILON+SIGMA)*IDENTITY+bg_pol+gain_pol+raman_pol
!  ...endif NONLINEAR_FLAG
      endif
!
!.....................................................
!...............toggle PML............................
!
      if(.not. USE_PML) then
         JJstretch      = ZERO
         JJstretch(1,1) = ZONE
         JJstretch(2,2) = ZONE
         JJstretch(3,3) = ZONE
         detJstretch    = ZONE
      else
!     ...get PML function
         call get_Beta(x,Fld_flag, zbeta,zdbeta,zd2beta)
         Jstretch(3,3) = zdbeta
!     ...compute det(J) * J^-1 * J^-T
!        (J is a diagonal matrix)
         invJstretch(3,3) = 1.d0/zdbeta
         call ZGEMM('N', 'N', 3, 3, 3, ZONE, invJstretch, 3, &
                       invJstretch, 3, ZERO, JJstretch, 3)
         detJstretch = zdbeta
         JJstretch = detJstretch*JJstretch
      endif
!
!  ...PML stretching
      zaJ(1,1) = JJstretch(1,1)*za(1,1)
      zaJ(2,2) = JJstretch(2,2)*za(2,2)
      zaJ(3,3) = JJstretch(3,3)*za(3,3)
!
      zcJ(1,1) = JJstretch(1,1)*zc(1,1)
      zcJ(2,2) = JJstretch(2,2)*zc(2,2)
      zcJ(3,3) = JJstretch(3,3)*zc(3,3)
!
!...........end toggle PML............................
!.....................................................
!
!  ...loop through enriched H(curl) test functions in the enriched space
      do k1=1,nrdofEE
!
!     ...Piola transformation
         do i = 1,3
            daux = dxdxi(i,1:3)
            call dot_product(shapEE(1:3,k1),dxidx(1:3,i), fldF(i))
            call dot_product(curlEE(1:3,k1),daux        , crlF(i))
         enddo
         crlF(1:3) = crlF(1:3)/rjac
         fldG = fldF; crlG = crlF
!     ...e_z x F, e_z x G
         rotF = 0.d0; rotF(1) = -fldF(2); rotF(2) = fldF(1); rotG = rotF
!
!        RHS:
!        (J^imp,F) first  equation RHS (with first H(curl) test function F)
!        (0    ,G) second equation RHS is zero
         n = 2*k1-1
         bload_E(n) = bload_E(n)                                   &
                    + (fldF(1)*zJ(1)+fldF(2)*zJ(2)+fldF(3)*zJ(3))  &
                    * weight
!
!     ...loop through L2 trial shape functions
         do k2=1,nrdofQ
!
!        ...Piola transformation
            fldE(1:3) = shapQ(k2)/rjac; fldH = fldE
!
!        ...testing with F (first H(curl) test function))
            n = 2*k1-1
!        ...-iωε(E,F)
            m = (k2-1)*6+1
            stiff_EQ_T(m,n) = stiff_EQ_T(m,n) - zaJ(1,1)*fldE(1)*fldF(1)*weight
            m = (k2-1)*6+2
            stiff_EQ_T(m,n) = stiff_EQ_T(m,n) - zaJ(2,2)*fldE(2)*fldF(2)*weight
            m = (k2-1)*6+3
            stiff_EQ_T(m,n) = stiff_EQ_T(m,n) - zaJ(3,3)*fldE(3)*fldF(3)*weight
!        ...(H,curl(F))
            m = (k2-1)*6+4
            stiff_EQ_T(m,n) = stiff_EQ_T(m,n) + fldH(1)*crlF(1)*weight
            m = (k2-1)*6+5
            stiff_EQ_T(m,n) = stiff_EQ_T(m,n) + fldH(2)*crlF(2)*weight
            m = (k2-1)*6+6
            stiff_EQ_T(m,n) = stiff_EQ_T(m,n) + fldH(3)*crlF(3)*weight
!
!        ...additional stiffness contribution if solving vectorial envelope equation
            if (ENVELOPE) then
!           ...-ik(e_z x H,F), where e_z x H = (-H_y,H_x,0), or
!           ...+ik(H,e_z x F), where e_z x F = (-F_y,F_x,0)
               m = (k2-1)*6+4
               stiff_EQ_T(m,n) = stiff_EQ_T(m,n) + ZI*WAVENUM_FLD*detJstretch*fldH(1)*rotF(1)*weight
               m = (k2-1)*6+5
               stiff_EQ_T(m,n) = stiff_EQ_T(m,n) + ZI*WAVENUM_FLD*detJstretch*fldH(2)*rotF(2)*weight
            endif
!
!        ...testing with G (second H(curl) test function))
            n = 2*k1
!        ...(E,curl(G))
            m = (k2-1)*6+1
            stiff_EQ_T(m,n) = stiff_EQ_T(m,n) + fldE(1)*crlG(1)*weight
            m = (k2-1)*6+2
            stiff_EQ_T(m,n) = stiff_EQ_T(m,n) + fldE(2)*crlG(2)*weight
            m = (k2-1)*6+3
            stiff_EQ_T(m,n) = stiff_EQ_T(m,n) + fldE(3)*crlG(3)*weight
!
!        ...iωμ(H,G)
            m = (k2-1)*6+4
            stiff_EQ_T(m,n) = stiff_EQ_T(m,n) + zcJ(1,1)*fldH(1)*fldG(1)*weight
            m = (k2-1)*6+5
            stiff_EQ_T(m,n) = stiff_EQ_T(m,n) + zcJ(2,2)*fldH(2)*fldG(2)*weight
            m = (k2-1)*6+6
            stiff_EQ_T(m,n) = stiff_EQ_T(m,n) + zcJ(3,3)*fldH(3)*fldG(3)*weight
!
!        ...additional stiffness contribution if solving vectorial envelope equation
            if (ENVELOPE) then
!           ...-ik(e_z x E,G), where e_z x E = (-E_y,E_x,0), or
!           ...+ik(E,e_z x G), where e_z x G = (-G_y,G_x,0)
               m = (k2-1)*6+1
               stiff_EQ_T(m,n) = stiff_EQ_T(m,n) + ZI*WAVENUM_FLD*detJstretch*fldE(1)*rotG(1)*weight
               m = (k2-1)*6+2
               stiff_EQ_T(m,n) = stiff_EQ_T(m,n) + ZI*WAVENUM_FLD*detJstretch*fldE(2)*rotG(2)*weight
            endif
!
!     ...end of loop through L2 trial functions
         enddo
!
!     ...loop through enriched H(curl) test functions in the enriched space
         do k2=k1,nrdofEE
!
!        ...Piola transformation
            do i = 1,3
               daux = dxdxi(i,1:3)
               call dot_product(shapEE(1:3,k2),dxidx(1:3,i), fldE(i))
               call dot_product(curlEE(1:3,k2),daux        , crlE(i))
            enddo
            crlE(1:3) = crlE(1:3)/rjac
!        ...e_z x E
            rotE = 0.d0; rotE(1) = -fldE(2); rotE(2) = fldE(1)
!
            call dot_product(fldF,fldE, FF)
            call dot_product(fldF,crlE, FC)
            call dot_product(crlF,fldE, CF)
            call dot_product(crlF,crlE, CC)
!
!        ...accumulate for the Hermitian Gram matrix
!           (compute upper triangular only)
!        ...testNorm = Scaled Adjoint Graph norm
!              ||v|| = alpha*(v,v) + (A^* v, A^* v)
!           (first eqn multiplied by F, second eqn by G)
!           G_ij=(phi_j,phi_i)_testNorm is 2x2 matrix
!           where (phi_j,phi_i)_l2Norm = Int[phi_i^* phi_j]
!           and phi_i = (F_i,G_i), phi_j = (F_j,G_j).
!           -------------------------
!           | (F_i,F_j)   (F_i,G_j) |
!           | (G_i,F_j)   (G_i,G_j) |
!           -------------------------
!           F_i/G_i are outer loop shape functions (fldF)
!           F_j/G_j are inner loop shape functions (fldE)
!
!           (F_j,F_i) terms = Int[F_^*i F_j] terms (G_11)
            n = 2*k1-1; m = 2*k2-1
            k = ij_upper_to_packed(n,m)
            zaux = abs(zaJ(1,1))**2*fldF(1)*fldE(1) + &
                   abs(zaJ(2,2))**2*fldF(2)*fldE(2) + &
                   abs(zaJ(3,3))**2*fldF(3)*fldE(3)
            gramP(k) = gramP(k) &
                     + (zaux + ALPHA_NORM*FF + CC)*weight
!
            if (ENVELOPE) then
!           ...ik(e_z x F_i, curl F_j)
               zaux = ZI*WAVENUM_FLD*detJstretch* &
                     (rotF(1)*crlE(1)+rotF(2)*crlE(2)+rotF(3)*crlE(3))
!           ...k^2(e_z x F_i, e_z x F_j)
               zbux = (WAVENUM_FLD*abs(detJstretch))**2 * &
                     (rotF(1)*rotE(1)+rotF(2)*rotE(2)+rotF(3)*rotE(3))
!           ...-ik(curl F_i, e_z x F_j)
               zcux = -ZI*WAVENUM_FLD*conjg(detJstretch)* &
                     (crlF(1)*rotE(1)+crlF(2)*rotE(2)+crlF(3)*rotE(3))
               gramP(k) = gramP(k) + (zaux+zbux+zcux)*weight
            endif
!
!           (G_j,F_i) terms = Int[F_^*i G_j] terms (G_12)
            n = 2*k1-1; m = 2*k2
            k = ij_upper_to_packed(n,m)
            zaux = - (zaJ(1,1)*fldF(1)*crlE(1) + &
                      zaJ(2,2)*fldF(2)*crlE(2) + &
                      zaJ(3,3)*fldF(3)*crlE(3) )
            zcux = conjg(zcJ(1,1))*crlF(1)*fldE(1) + &
                   conjg(zcJ(2,2))*crlF(2)*fldE(2) + &
                   conjg(zcJ(3,3))*crlF(3)*fldE(3)
            gramP(k) = gramP(k) + (zaux+zcux)*weight
!
            if (ENVELOPE) then
!           ...ik(iωε F_i, e_z x G_j)
               zaux = ZI*WAVENUM_FLD*conjg(detJstretch)* &
                        (zaJ(1,1)*fldF(1)*rotE(1) + &
                         zaJ(2,2)*fldF(2)*rotE(2) + &
                         zaJ(3,3)*fldF(3)*rotE(3) )
!           ...ik(e_z x F_i, (iωμ)^* G_j)
               zcux = ZI*WAVENUM_FLD*detJstretch* &
                        (conjg(zcJ(1,1))*rotF(1)*fldE(1) + &
                         conjg(zcJ(2,2))*rotF(2)*fldE(2) + &
                         conjg(zcJ(3,3))*rotF(3)*fldE(3) )
               gramP(k) = gramP(k) + (zaux+zcux)*weight
            endif
!
!        ...compute lower triangular part of 2x2 G_ij matrix
!           only if it is not a diagonal element, G_ii
            if (k1 .ne. k2) then
!              (F_j,G_i) terms = Int[G_^*i F_j] terms (G_21)
               n = 2*k1; m = 2*k2-1
               k = ij_upper_to_packed(n,m)
               zaux = - (conjg(zaJ(1,1))*crlF(1)*fldE(1) + &
                         conjg(zaJ(2,2))*crlF(2)*fldE(2) + &
                         conjg(zaJ(3,3))*crlF(3)*fldE(3) )
               zcux = zcJ(1,1)*fldF(1)*crlE(1) + &
                      zcJ(2,2)*fldF(2)*crlE(2) + &
                      zcJ(3,3)*fldF(3)*crlE(3)
               gramP(k) = gramP(k) + (zaux+zcux)*weight
!
               if (ENVELOPE) then
!              ...-ik (e_z x G_i, (iωε)^* F_j)
                  zaux = -ZI*WAVENUM_FLD*detJstretch* &
                         (conjg(zaJ(1,1))*rotF(1)*fldE(1) + &
                          conjg(zaJ(2,2))*rotF(2)*fldE(2) + &
                          conjg(zaJ(3,3))*rotF(3)*fldE(3) )
!              ...-ik (iωμ G_i, e_z x F_j)
                  zcux = -ZI*WAVENUM_FLD*conjg(detJstretch)* &
                         (zcJ(1,1)*fldF(1)*rotE(1) + &
                          zcJ(2,2)*fldF(2)*rotE(2) + &
                          zcJ(3,3)*fldF(3)*rotE(3) )
                  gramP(k) = gramP(k) + (zaux+zcux)*weight
               endif
            endif
!
!           (G_j,G_i) terms = Int[G_^*i G_j] terms (G_22)
            n = 2*k1; m = 2*k2
            k = ij_upper_to_packed(n,m)
            zcux = abs(zcJ(1,1))**2*fldF(1)*fldE(1) + &
                   abs(zcJ(2,2))**2*fldF(2)*fldE(2) + &
                   abs(zcJ(3,3))**2*fldF(3)*fldE(3)
            gramP(k) = gramP(k) &
                     + (zcux + ALPHA_NORM*FF + CC)*weight
!
            if (ENVELOPE) then
!           ...ik(e_z x G_i, curl G_j)
               zaux = ZI*WAVENUM_FLD*detJstretch* &
                     (rotF(1)*crlE(1)+rotF(2)*crlE(2)+rotF(3)*crlE(3))
!           ...k^2(e_z x G_i, e_z x G_j)
               zbux = (WAVENUM_FLD*abs(detJstretch))**2 * &
                     (rotF(1)*rotE(1)+rotF(2)*rotE(2)+rotF(3)*rotE(3))
!           ...-ik(curl G_i, e_z x G_j)
               zcux = -ZI*WAVENUM_FLD*conjg(detJstretch)* &
                     (crlF(1)*rotE(1)+crlF(2)*rotE(2)+crlF(3)*rotE(3))
               gramP(k) = gramP(k) + (zaux+zbux+zcux)*weight
            endif
!
!     ...end of loop through enriched H(curl) test functions
         enddo
!
!  ...end of loop through enriched H(curl) test functions
      enddo
!
!..end of loop through integration points
   enddo
!
!..end timer
!   end_time = MPI_Wtime()
!   !$OMP CRITICAL
!      !write(*,10) S_Type(etype), end_time-start_time
!      write(*,11) end_time-start_time
!! 10   format(A,' elem : ',f12.5,'  seconds')
! 11   format(f12.5)
!   !$OMP END CRITICAL
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                 B O U N D A R Y   I N T E G R A L S
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
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
!     ...determine discontinuous Hcurl shape functions
         call shape3EE(etype,xi,nordP, nrdof,shapEE,curlEE)
#if HP3D_DEBUG
         if (nrdof .ne. NrdofEE) then
            write(*,*) 'elem_maxwell: INCONSISTENCY NrdofEE. stop.'
            stop
         endif
#endif
!
!     ...determine element H1 shape functions (for geometry)
         call shape3DH(etype,xi,norder,norient_edge,norient_face, &
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
         call shape3DE(etype,xi,norderi,norient_edge,norient_face, &
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

!     ...loop through Hcurl enriched test functions
         !do ik=1,nrdofEEi
         do k1=1,nrdofEE
            !k1 = idxEE(ik)
            E1(1:3) = shapEE(1,k1)*dxidx(1,1:3) &
                    + shapEE(2,k1)*dxidx(2,1:3) &
                    + shapEE(3,k1)*dxidx(3,1:3)
!        ...check for impedance BC (elimination strategy)
!           (impedance constant is GAMMA for TE10 mode in rectangular waveguide)
!           ( < n x H , G > = GAMMA*< n x n x E , G > + < zg , G > )
            if ((ibc(ifc,3).eq.3 .and. Fld_flag.eq.1) .or. &
                (ibc(ifc,5).eq.3 .and. Fld_flag.eq.0)) then
!           ...get the boundary source [zImp should be zero here]
               call get_bdSource(Mdle,x,rn, zImp)
!           ...accumulate for the load vector
               k = 2*k1-1
               bload_E(k) = bload_E(k) &
                          - (zImp(1)*E1(1)+zImp(2)*E1(2)+zImp(3)*E1(3))*weight
!        ...end if for impedance BC
            endif
!        ...loop through H(curl) trial functions
            do k2=1,NrdofEi
               E2(1:3) = shapE(1,k2)*dxidx(1,1:3) &
                       + shapE(2,k2)*dxidx(2,1:3) &
                       + shapE(3,k2)*dxidx(3,1:3)
!
               call cross_product(rn,E2, rntimesE)
!           ...check for impedance BC (elimination strategy)
               if ((ibc(ifc,3).eq.3 .and. Fld_flag.eq.1) .or. &
                   (ibc(ifc,5).eq.3 .and. Fld_flag.eq.0)) then
!              ...accumulate for the extended stiffness matrix on IBC
                  call cross_product(rn,rntimesE, rn2timesE)
                  stiff_EE_T(2*k2-1,2*k1-1) = stiff_EE_T(2*k2-1,2*k1-1) &
                                               + (E1(1)*rn2timesE(1) &
                                               +  E1(2)*rn2timesE(2) &
                                               +  E1(3)*rn2timesE(3) &
                                                 )*GAMMA*weight
!
               else
!              ...accumulate for the extended stiffness matrix without IBC
                  stiff_EE_T(2*k2,2*k1-1) = stiff_EE_T(2*k2,2*k1-1) &
                                             + (E1(1)*rntimesE(1) &
                                             +  E1(2)*rntimesE(2) &
                                             +  E1(3)*rntimesE(3) &
                                               )*weight
!           ...end if for impedance BC
               endif
               stiff_EE_T(2*k2-1,2*k1) = stiff_EE_T(2*k2-1,2*k1) &
                                          + (E1(1)*rntimesE(1) &
                                          +  E1(2)*rntimesE(2) &
                                          +  E1(3)*rntimesE(3) &
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
!-------------------------------------------------------------------------------
!      Construction of the DPG system
!-------------------------------------------------------------------------------
!
   allocate(stiff_ALL(NrTest,NrTrial+1))
!
!..Total test/trial DOFs of the element
   i1 = NrTest ; j1 = 2*NrdofEi ; j2 = 6*NrdofQ
!
   stiff_ALL(1:i1,1:j1)       = transpose(stiff_EE_T(1:j1,1:i1))
   stiff_ALL(1:i1,j1+1:j1+j2) = transpose(stiff_EQ_T(1:j2,1:i1))
   stiff_ALL(1:i1,j1+j2+1)    = bload_E(1:i1)
!
   deallocate(stiff_EE_T,stiff_EQ_T)
!
!..A. Compute Cholesky factorization of Gram Matrix, G=U^*U (=LL^*)
   call ZPPTRF('U',NrTest,gramP,info)
   if (info.ne.0) then
      write(*,*) 'elem_maxwell: ZPPTRF: Mdle,info = ',Mdle,info,'. stop.'
      stop
   endif
!
!..B. Solve triangular system to obtain B~, (LX=) U^*X = [B|l]
   call ZTPTRS('U','C','N',NrTest,NrTrial+1,gramP,stiff_ALL,NrTest,info)
   if (info.ne.0) then
      write(*,*) 'elem_maxwell: ZTPTRS: Mdle,info = ',Mdle,info,'. stop.'
      stop
   endif
!
   deallocate(gramP)
   allocate(zaloc(NrTrial+1,NrTrial+1)); zaloc = ZERO
!
!..C. Matrix multiply: B^* G^-1 B (=B~^* B~)
   call ZHERK('U','C',NrTrial+1,NrTest,ZONE,stiff_ALL,NrTest,ZERO,zaloc,NrTrial+1)
!
   deallocate(stiff_ALL)
!
!..D. Fill lower triangular part of Hermitian matrix
   do i=1,NrTrial
      zaloc(i+1:NrTrial+1,i) = conjg(zaloc(i,i+1:NrTrial+1))
   enddo
!
!..E. Fill ALOC and BLOC matrices
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
!..Implementation of impedance BC via L2 penalty term
   if (IBCFLAG.eq.2) call imp_penalty(Mdle,Fld_flag,NrdofH,NrdofEi,MdE, &
                                      norder,norderi, ZblocE,ZalocEE)
!
!-------------------------------------------------------------------------------
!       D E B U G   T E S T S
!-------------------------------------------------------------------------------
#if HP3D_DEBUG
   iprint = 0
   if (iprint.ge.1) then
      write(*,7010)
7010  format('elem_maxwell: ZblocE,ZblocQ = ')
      write(*,7011) ZblocE(1:2*NrdofE)
      write(*,7011) ZblocQ(1:6*NrdofQ)
7011  format(12e12.5)
      call pause
      write(*,7012)
7012  format('elem_maxwell: ZalocEE = ')
      do i=1,2*NrdofE
         write(*,7013) i,ZalocEE(i,1:2*NrdofE)
7013     format('i = ',i3,10(/,6(2e12.5,2x)))
      enddo
      call pause
      write(*,7014)
7014  format('elem_maxwell: ZalocEQ = ')
      do i=1,2*NrdofE
         write(*,7013) i,ZalocEQ(i,1:6*NrdofQ)
      enddo
      call pause
      write(*,7015)
7015  format('elem_maxwell: ZalocQQ = ')
      do i=1,6*NrdofQ
         write(*,7013) i,ZalocQQ(i,1:6*NrdofQ)
      enddo
      call pause
   endif
#endif
!
end subroutine elem_maxwell
!
!
!-------------------------------------------------------------------------------
!
!    routine name      - imp_penalty
!
!-------------------------------------------------------------------------------
!
!     latest revision:  - June 2021
!
!> @brief         - routine adds impedance L2 penalty terms to the
!                         stiffness matrix and load vector for the
!                         UW DPG formulation for Maxwell equations
!
!     arguments:
!        in:
!           Mdle        - element middle node number
!           Fld_flag    - field flag (0: pump, 1: signal)
!           NrdofH      - number of H1 trial dof
!           NrdofEi     - number of H(curl) trial interface dof
!           MdE         - num rows of ZalocEE
!           Norder      - element nodes order (trial)
!           Norderi     - element nodes order (trial) for interfaces
!        inout:
!           ZblocE      - load vector
!           ZalocEE     - stiffness matrix
!
!-------------------------------------------------------------------------------
!
#include "typedefs.h"
!
subroutine imp_penalty(Mdle,Fld_flag,NrdofH,NrdofEi,MdE,Norder,Norderi, &
                       ZblocE,ZalocEE)
!
!..modules used
   use control
   use parametersDPG
   use data_structure3D
   use laserParam
   use commonParam
!..no implicit statements
   implicit none
!..declare input/output variables
   integer, intent(in)    :: Mdle
   integer, intent(in)    :: Fld_flag
   integer, intent(in)    :: NrdofH
   integer, intent(in)    :: NrdofEi
   integer, intent(in)    :: MdE
   integer, intent(in)    :: Norder(19)
   integer, intent(in)    :: Norderi(19)
   VTYPE,   intent(inout) :: ZblocE(MdE)
   VTYPE,   intent(inout) :: ZalocEE(MdE,MdE)
!
!..declare edge/face type variables
   integer :: etype,ftype
!
!..declare orientation for edges and faces
   integer, dimension(12)    :: norient_edge
   integer, dimension(6)     :: norient_face
!
!..face order
   integer, dimension(5)     :: norderf
!
!..geometry dof (work space for nodcor)
   real(8), dimension(3,MAXbrickH) :: xnod
!
!..variables for geometry
   real(8), dimension(3)    :: xi,x,rn
   real(8), dimension(3,2)  :: dxidt,dxdt
   real(8), dimension(3,3)  :: dxdxi,dxidx
   real(8), dimension(2)    :: t
!
!..H1 shape functions
   real(8), dimension(MAXbrickH)    :: shapH
   real(8), dimension(3,MAXbrickH)  :: gradH
!
!..H(curl) shape functions
   real(8), dimension(3,MAXbrickE)  :: shapE
   real(8), dimension(3,MAXbrickE)  :: curlE
!
!..2D quadrature data
   real(8), dimension(2,MAXNINT2ADD)  :: tloc
   real(8), dimension(MAXNINT2ADD)    :: wtloc
!
!..BC's flags
   integer, dimension(6,NRINDEX)      :: ibc
!
!..Maxwell load and auxiliary variables
   VTYPE  , dimension(3) :: zImp
   real(8), dimension(3) :: E2,rntimesE,rn2timesE
   real(8), dimension(3) :: F1,rntimesF,rn2timesF
!
!..various variables for the problem
   real(8) :: rjac,weight,eps,bjac
   integer :: k1,k2,k,l,nint
   integer :: nrdof,nsign,ifc,nrf
!
!-------------------------------------------------------------------------------
!
   if (IBCFLAG .ne. 2) then
      write(*,*) 'imp_penalty called for IBCFLAG.ne.2, returning...'
      return
   endif
!
!..define the weight of the penalty term
   eps = 1.d0
!
!..no need for overintegration in penalty term
   INTEGRATION = 0
!
!..determine element type and number of faces
   etype = NODES(Mdle)%ntype
   nrf = nface(etype)
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
!..loop through element faces
   do ifc=1,nrf
!
!  ...skip the face if it is not on the impedance boundary
      if ((ibc(ifc,3).ne.2 .or. Fld_flag.ne.1) .and. &
          (ibc(ifc,5).ne.2 .or. Fld_flag.ne.0)) cycle
!
!  ...sign factor to determine the OUTWARD normal unit vector
      nsign = nsign_param(etype,ifc)
!
!  ...face type
      ftype = face_type(etype,ifc)
!
!  ...face order of approximation
      call face_order(etype,ifc,Norder, norderf)
!
!  ...set 2D quadrature
      call set_2D_int(ftype,norderf,norient_face(ifc), nint,tloc,wtloc)
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
!     ...determine element H1 shape functions (for geometry)
         call shape3DH(etype,xi,Norder,norient_edge,norient_face, &
                       nrdof,shapH,gradH)
!
!     ...determine element H(curl) shape functions (for fluxes)
!     ...for interfaces only (no bubbles)
         call shape3DE(etype,xi,Norderi,norient_edge,norient_face, &
                       nrdof,shapE,curlE)
!
!     ...geometry
         call bgeom3D(Mdle,xi,xnod,shapH,gradH,NrdofH,dxidt,nsign, &
                      x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
         weight = bjac*wtloc(l)
!
!     ...get the boundary source [zImp should be zero here]
         call get_bdSource(Mdle,x,rn, zImp)
!
!     ...loop through H(curl) test functions
         do k1=1,NrdofEi
            F1(1:3) = shapE(1,k1)*dxidx(1,1:3) &
                    + shapE(2,k1)*dxidx(2,1:3) &
                    + shapE(3,k1)*dxidx(3,1:3)
!
            call cross_product(rn,F1, rntimesF)
            call cross_product(rn,rntimesF, rn2timesF)
!
!        ...compute impedance BC equation for the load
!           (impedance constant is GAMMA for TE10 mode in rectangular waveguide)
!           ( < n x H , G > = GAMMA*< n x n x E , G > + < zg , G > )
!
!        ...accumulate directly into load vector ZblocE(1:2*NrdofEi)
!           (zImp,F): -gamma * (zImp , n x n x F)
            k = 2*k1-1
            ZblocE(k) = ZblocE(k) - (rn2timesF(1)*zImp(1)  &
                                  +  rn2timesF(2)*zImp(2)  &
                                  +  rn2timesF(3)*zImp(3)  &
                                    )*GAMMA*weight/eps
!           (zImp,G): (zImp , n x G)
            k = 2*k1
            ZblocE(k) = ZblocE(k) + (rntimesF(1)*zImp(1)  &
                                  +  rntimesF(2)*zImp(2)  &
                                  +  rntimesF(3)*zImp(3)  &
                                    )*weight/eps
!
!        ...loop through H(curl) trial functions
            do k2=1,NrdofEi
               E2(1:3) = shapE(1,k2)*dxidx(1,1:3) &
                       + shapE(2,k2)*dxidx(2,1:3) &
                       + shapE(3,k2)*dxidx(3,1:3)
!
               call cross_product(rn,E2, rntimesE)
               call cross_product(rn,rntimesE, rn2timesE)
!           ...accumulate directly into stiffness ZalocEE(1:2*NrdofEi,1:2*NrdofEi)
!              4 contributions: trial (E,H), test (F,G)
!              (E,F):  gamma^2 * (n x n x E , n x n x F)
               ZalocEE(2*k1-1,2*k2-1) = ZalocEE(2*k1-1,2*k2-1) &
                                        + (rn2timesF(1)*rn2timesE(1)  &
                                        +  rn2timesF(2)*rn2timesE(2)  &
                                        +  rn2timesF(3)*rn2timesE(3)  &
                                          )*GAMMA*GAMMA*weight/eps
!              (H,F): -gamma   * (n x H , n x n x F)
               ZalocEE(2*k1-1,2*k2  ) = ZalocEE(2*k1-1,2*k2  ) &
                                        - (rn2timesF(1)*rntimesE(1)   &
                                        +  rn2timesF(2)*rntimesE(2)   &
                                        +  rn2timesF(3)*rntimesE(3)   &
                                          )*GAMMA*weight/eps
!              (E,G): -gamma   * (n x n x E , n x G)
               ZalocEE(2*k1  ,2*k2-1) = ZalocEE(2*k1  ,2*k2-1) &
                                        - (rntimesF(1)*rn2timesE(1)  &
                                        +  rntimesF(2)*rn2timesE(2)  &
                                        +  rntimesF(3)*rn2timesE(3)  &
                                          )*GAMMA*weight/eps
!              (H,G):            (n x H , n x G)
               ZalocEE(2*k1  ,2*k2  ) = ZalocEE(2*k1  ,2*k2  ) &
                                        + (rntimesF(1)*rntimesE(1)   &
                                        +  rntimesF(2)*rntimesE(2)   &
                                        +  rntimesF(3)*rntimesE(3)   &
                                          )*weight/eps
!        ...end loop through H(curl) trial functions
            enddo
!     ...end loop through the H(curl) test functions
         enddo
!  ...end loop through integration points
      enddo
!..end loop through faces
   enddo
!
!
end subroutine imp_penalty
