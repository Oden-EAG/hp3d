!
!--------------------------------------------------------------------
!
!    routine name      - elem_maxwell_fi_hexa
!
!--------------------------------------------------------------------
!
!     latest revision:  - Oct 2019
!
!     purpose:          - routine returns unconstrained (ordinary)
!                         stiffness matrix and load vector for the
!                         UW DPG formulation for Maxwell equations
!                       - uses sum factorization for fast integration
!                         with swapped loops
!                       - for hexahedral element only
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
!---------------------------------------------------------------------
!
#include "typedefs.h"
!
subroutine elem_maxwell_fi_hexa(Mdle,Fld_flag,                &
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
   character(len=4) :: etype,ftype
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
   real(8) :: xnod(3,MAXbrickH)
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
   real(8), dimension(3)    :: xi,x,rn
   real(8), dimension(3,2)  :: dxidt,dxdt,rt
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
   real(8), dimension(3,MAXbrickEE) :: shapEE
   real(8), dimension(3,MAXbrickEE) :: curlEE
!
!..nrdof for interface only (without bubbles)
   integer :: nrdofEEi
!
!..H(curl) bubble index
   integer, allocatable :: idxEE(:)
!
!..element mdle node dof
   integer :: ndofHHmdl,ndofEEmdl,ndofVVmdl,ndofQQmdl
!
!..load vector for the enriched space
   complex(8) :: bload_E(NrTest)
!
!..Gram matrix in packed format
   !complex(8) :: gramP(NrTest*(NrTest+1)/2)
   complex(8), allocatable :: gramP(:)
   real(8) :: FF, CF, FC
   real(8) :: fldE(3), fldH(3), crlE(3), crlH(3)
   real(8) :: fldF(3), fldG(3), crlF(3), crlG(3)
!
!..matrices for transpose filling (swapped loops)
!..stiffness matrices (transposed) for the enriched test space
   !complex(8) :: stiff_EE_T(2*NrdofEi,NrTest)
   !complex(8) :: stiff_EQ_T(6*NrdofQ ,NrTest)
   !complex(8) :: stiff_ALL(NrTest   ,NrTrial+1)
   !complex(8) :: zaloc    (NrTrial+1,NrTrial+1)
   complex(8), allocatable :: stiff_EE_T(:,:),stiff_EQ_T(:,:)
   complex(8), allocatable :: stiff_ALL(:,:),zaloc(:,:)
!
!..3D quadrature data
   real(8), dimension(3,MAXNINT3ADD)  :: xiloc
   real(8), dimension(  MAXNINT3ADD)  :: waloc
!
!..2D quadrature data
   real(8), dimension(2,MAXNINT2ADD)  :: tloc
   real(8), dimension(  MAXNINT2ADD)  :: wtloc
!
!..BC's flags
   integer, dimension(6,NR_PHYSA)     :: ibc
!
!..for auxiliary computation
   complex(8) :: zaux
!
!..Maxwell load and auxiliary variables
   complex(8), dimension(3) :: zJ,zImp
   real(8)   , dimension(3) :: E1,E2,rntimesE,rn2timesE
!
!..number of edge,faces per element type
   integer :: nre, nrf
!
!..various variables for the problem
   real(8) :: h_elem,rjac,weight,wa,v2n,CC,EE,CE,E,EC,q,h,omeg,alpha_scale
   real(8) :: bjac
   integer :: i1,i2,j1,j2,k1,k2,kH,kk,i,ik,j,k,l,nint,kE,n,m
   integer :: iflag,iprint,itime,iverb
   integer :: nrdof,nordP,nsign,ifc,ndom,info,icomp,idec
   complex(8) :: zfval
   complex(8) :: za(3,3),zc(3,3)
   complex(8) :: zaJ(3,3),zcJ(3,3)
!
!..for polarizations function
   complex(8), dimension(3,3) :: bg_pol,gain_pol,raman_pol
   real(8) :: delta_n
   integer :: dom_flag
!
!..OMEGA_RATIO_SIGNAL or OMEGA_RATIO_PUMP
   real(8) :: OMEGA_RATIO_FLD
!
!..for PML
   VTYPE :: zbeta,zdbeta,zd2beta,detJstretch
   VTYPE, dimension(3,3) :: Jstretch,invJstretch,JJstretch
!
!..added to use fast integration
   VTYPE, allocatable :: AUXEE_A_zb(:,:,:)    , AUXEE_A_zc(:,:,:)
   VTYPE, allocatable :: AUXEE_B_zb(:,:,:,:)  , AUXEE_B_zc(:,:,:,:)
   VTYPE, allocatable :: AUXCC_A(:,:,:,:,:)
   VTYPE, allocatable :: AUXCC_B(:,:,:,:,:,:)
   VTYPE, allocatable :: AUXEC_A_zb(:,:,:,:)  , AUXEC_A_zc(:,:,:,:)
   VTYPE, allocatable :: AUXCE_A_zb(:,:,:,:)  , AUXCE_A_zc(:,:,:,:)
   VTYPE, allocatable :: AUXEC_B_zb(:,:,:,:,:), AUXEC_B_zc(:,:,:,:,:)
   VTYPE, allocatable :: AUXCE_B_zb(:,:,:,:,:), AUXCE_B_zc(:,:,:,:,:)
!
   VTYPE, allocatable :: STIFQE_A(:,:,:)
   VTYPE, allocatable :: STIFQE_B(:,:,:,:)
   VTYPE, allocatable :: STIFQE_ALPHA_A(:,:,:)
   VTYPE, allocatable :: STIFQE_ALPHA_B(:,:,:,:)
   VTYPE, allocatable :: STIFQC_A(:,:,:,:)
   VTYPE, allocatable :: STIFQC_B(:,:,:,:,:)
   VTYPE, allocatable :: LOADE_A(:,:)
   VTYPE, allocatable :: LOADE_B(:,:,:)
!
   integer :: a,b,sa,sb,sc,alph,beta
   integer :: px,py,pz
   integer :: l1,l2,l3,i3,j3,k3,idxbeta,idxalph,idxa,idxa2,idxa3,idxb,idxb2,idxb3,m1,m2
   integer :: nord1,nord2,nord3,nintx,ninty,nintz
   integer :: nrdofH1,nrdofH2,nrdofH3
   integer :: nrdofH1_tr,nrdofH2_tr,nrdofH3_tr
   integer :: nrdofQ1_tr,nrdofQ2_tr,nrdofQ3_tr
   real(8) :: xi1,xi2,xi3,wt1,wt2,wt3,clock1,clock2
   real(8) :: wt123,weighthh,weightvv
   real(8), dimension(MAXPP+1) :: xilocx,xilocy,xilocz
   real(8), dimension(MAXPP+1) :: wlocx,wlocy,wlocz
   real(8), dimension(3,MAXNINT3ADD) :: wloc3
   real(8), dimension(3) :: xip,dHdx,dHHdx
   real(8), dimension(3,3) :: D_za,D_zc,D_aux,C,D
   VTYPE  , dimension(3,3) :: Z_za,Z_zc,Z_aux
   real(8), dimension(MAXPP+1,2) :: shapH1,shapH2,shapH3
   real(8), dimension(MAXPP+1,MAXPP+1) :: sH2p,sH3p,dsH2p,dsH3p
!
!..timer
   real(8) :: MPI_Wtime,start_time,end_time
!
   integer, dimension(3,3) :: deltak
!
!..for Gram matrix compressed storage format
   integer :: nk
   nk(k1,k2) = (k2-1)*k2/2+k1
!
!..Identity/Kronecker delta tensor
   deltak=ZERO
   do a=1,3
     deltak(a,a)=1
   enddo
!
!---------------------------------------------------------------------
!
!..Set iverb = 0/1 (Non-/VERBOSE)
   iverb = 0
!..Set iprint = 0/1 (Non-/VERBOSE)
   iprint = 0
#if DEBUG_MODE
   if (iprint.eq.1) then
      write(*,*) 'elem_maxwell_fi_hexa: Mdle = ', Mdle
   endif
#endif
!
!..allocate matrices
   allocate(gramP(NrTest*(NrTest+1)/2))
   allocate(stiff_EE_T(2*NrdofEi,NrTest))
   allocate(stiff_EQ_T(6*NrdofQ ,NrTest))
!
!..element type
   etype = NODES(Mdle)%type
   nre = nedge(etype); nrf = nface(etype)
!
!..determine order of approximation
   call find_order(Mdle, norder)
   norderi(1:nre+nrf) = norder(1:nre+nrf)
!
!..set the enriched order of approximation
   select case(etype)
      case('mdlb')
         nordP = NODES(Mdle)%order+NORD_ADD*111
         norderi(nre+nrf+1) = 111
      case default
         write(*,*) 'elem_maxwell_fi_hexa: unsupported etype param. stop.'
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
#if DEBUG_MODE
   if (iprint.eq.1) then
      write(*,7001) Mdle
7001  format('elem_maxwell_fi_hexa: BCFLAGS FOR Mdle = ',i5)
      do i=1,NR_PHYSA
         write(*,7002) PHYSA(i), ibc(1:nrf,i)
7002     format('     ATTRIBUTE = ',a6,' FLAGS = ',6i2)
      enddo
   endif
#endif
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
!..set OMEGA_RATIO_FLD
   select case(Fld_flag)
      case(0)
         OMEGA_RATIO_FLD = OMEGA_RATIO_PUMP
      case(1)
         OMEGA_RATIO_FLD = OMEGA_RATIO_SIGNAL ! 1.0d0
      case default
      write(*,*) 'elem_maxwell_fi_hexa: invalid Fld_flag param. stop.'
         stop
   end select
!
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
      case(2,3)
         write(*,*) 'elem_maxwell_fi_hexa: cannot have prism core geometry with fast integration. stop.'
         stop
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
               write(*,*) 'elem_maxwell_fi_hexa: unexpected ndom param. stop.'
               stop
         end select
!..end select case of GEOM_NO
   end select
!
   if(NONLINEAR_FLAG.eq.1) then
!  ...get current solution dofs
      call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
   endif
!
!-----------------------------------------------------------------------
!..here begins the setup for tensorized num quadrature for hexahedra
!..set up the element quadrature
   xiloc=ZERO
   wloc3=ZERO
   xilocx=ZERO
   xilocy=ZERO
   xilocz=ZERO
   wlocx=ZERO
   wlocy=ZERO
   wlocz=ZERO
   sa=ZERO
   sb=ZERO
!
!..use the enriched order to set the quadrature
   INTEGRATION = NORD_ADD
   call set_3Dint_fi(etype,norder, nord1,nord2,nord3,nintx,ninty,nintz,xiloc,wloc3)
   INTEGRATION = 0
!
!..set up # dof for each direction for 1D H1 test functions with order p+dp
   nrdofH1=nord1+1; nrdofH2=nord2+1; nrdofH3=nord3+1
!..set up # dof for each direction for 1D H1 trial functions with order p
   nrdofH1_tr=nrdofH1-NORD_ADD; nrdofH2_tr=nrdofH2-NORD_ADD; nrdofH3_tr=nrdofH3-NORD_ADD
!..set up # dof for each direction for 1D L2 trial functions with order p
   nrdofQ1_tr=nrdofH1_tr-1; nrdofQ2_tr=nrdofH2_tr-1; nrdofQ3_tr=nrdofH3_tr-1
!..Allocate the auxiliary arrays for sum factorization
   allocate(AUXEE_A_zb(3,3,nrdofH3**2))
   allocate(AUXEE_A_zc(3,3,nrdofH3**2))
   allocate(AUXEE_B_zb(3,3,nrdofH2**2,nrdofH3**2))
   allocate(AUXEE_B_zc(3,3,nrdofH2**2,nrdofH3**2))
   allocate(AUXCC_A(2,2,3,3,nrdofH3**2))
   allocate(AUXCC_B(2,2,3,3,nrdofH2**2,nrdofH3**2))
   allocate(AUXEC_A_zb(2,3,3,nrdofH3**2))
   allocate(AUXEC_A_zc(2,3,3,nrdofH3**2))
   allocate(AUXEC_B_zb(2,3,3,nrdofH2**2,nrdofH3**2))
   allocate(AUXEC_B_zc(2,3,3,nrdofH2**2,nrdofH3**2))
   allocate(AUXCE_A_zb(2,3,3,nrdofH3**2))
   allocate(AUXCE_A_zc(2,3,3,nrdofH3**2))
   allocate(AUXCE_B_zb(2,3,3,nrdofH2**2,nrdofH3**2))
   allocate(AUXCE_B_zc(2,3,3,nrdofH2**2,nrdofH3**2))
!
   allocate(STIFQE_A(3,3,nrdofQ3_tr*nrdofH3))
   allocate(STIFQE_B(3,3,nrdofQ2_tr*nrdofH2,nrdofQ3_tr*nrdofH3))
   allocate(STIFQE_ALPHA_A(3,3,nrdofQ3_tr*nrdofH3))
   allocate(STIFQE_ALPHA_B(3,3,nrdofQ2_tr*nrdofH2,nrdofQ3_tr*nrdofH3))
   allocate(STIFQC_A(2,3,3,nrdofQ3_tr*nrdofH3))
   allocate(STIFQC_B(2,3,3,nrdofQ2_tr*nrdofH2,nrdofQ3_tr*nrdofH3))
   allocate(LOADE_A(3,nrdofH3))
   allocate(LOADE_B(3,nrdofH2,nrdofH3))
!
!............consistency check
!-----------------------------------------------------------------------
#if DEBUG_MODE
!
!..total number of test functions (per field)
   nrdof=nord1*nrdofH2*nrdofH3+nrdofH1*nord2*nrdofH3+nrdofH1*nrdofH2*nord3
   if (NrdofEE .ne. nrdof) then
      write(*,*) 'elem_maxwell_fi_hexa: INCONSISTENCY NrdofEE. stop.'
      stop
   endif
!
!..number of trial functions (per component)
   nrdof=nrdofQ1_tr*nrdofQ2_tr*nrdofQ3_tr
   if (NrdofQ .ne. nrdof) then
      write(*,*) 'elem_maxwell_fi_hexa: INCONSISTENCY NrdofQ. stop.'
      stop
   endif
!
#endif
!........end consistency check
!-----------------------------------------------------------------------
!
!..Quadrature points and weights are copied into 1-dim vectors
   l=0
   do l3=1,nintz
      do l2=1,ninty
         do l1=1,nintx
            l=l+1
            xilocx(l1)=xiloc(1,l)
            xilocy(l2)=xiloc(2,l)
            xilocz(l3)=xiloc(3,l)
            wlocx(l1)=wloc3(1,l)
            wlocy(l2)=wloc3(2,l)
            wlocz(l3)=wloc3(3,l)
         enddo
      enddo
   enddo
!
   xip=ZERO
!
!..start timer
   start_time = MPI_Wtime()
!
!..Loop over quadrature points in direction \xi_1
   do px=1,nintx
!  ...read quadrature point location and weight
      xi1=xilocx(px)
      wt1=wlocx(px)
!  ...call 1D shape functions for coordinate 1
!  ...Function values stored in shapH1(:,1), derivative values in shapH1(:,2)
      call shape1HH(xi1,nord1, nrdofH1,shapH1(:,1),shapH1(:,2))
!
!  ...Initialize auxiliary matrices B: Gram matrix
      AUXEE_B_zb = ZERO; AUXEE_B_zc = ZERO
      AUXCC_B    = ZERO
      AUXEC_B_zb = ZERO; AUXEC_B_zc = ZERO
      AUXCE_B_zb = ZERO; AUXCE_B_zc = ZERO
!  ...Initialize auxiliary matrices B: Stiffness and load
      STIFQE_B       = ZERO
      STIFQE_ALPHA_B = ZERO
      STIFQC_B       = ZERO
      LOADE_B        = ZERO
!
!  ...loop over quadrature points in direction \xi_2
      do py=1,ninty
!     ...read quadrature point location and weight
         xi2=xilocy(py)
         wt2=wlocy(py)
!     ...Shape function subroutine is called only once, when
!        px=1 and stored in sH2p(:,py) and dsH2p(:,py)
         if (px.eq.1) then
            sH2p(:,py)=0.d0
            dsH2p(:,py)=0.d0
            call shape1HH(xi2,nord2,nrdofH2,sH2p(:,py),dsH2p(:,py))
         endif
!     ...Copy shape functions in coord. 2 previously evaluated
         shapH2(:,1)=sH2p(:,py)
         shapH2(:,2)=dsH2p(:,py)
!     ...Initialize auxiliary matrices A: Gram matrix
         AUXEE_A_zb = ZERO; AUXEE_A_zc = ZERO
         AUXCC_A    = ZERO
         AUXEC_A_zb = ZERO; AUXEC_A_zc = ZERO
         AUXCE_A_zb = ZERO; AUXCE_A_zc = ZERO
!     ...Initialize auxiliary matrices A: Stiffness matrix and load vector
         STIFQE_A       = ZERO
         STIFQE_ALPHA_A = ZERO
         STIFQC_A       = ZERO
         LOADE_A        = ZERO
!
!     ...loop over quadrature points in direction \xi_3
         do pz=1,nintz
!        ...read quadrature point location and weight
            xi3=xilocz(pz)
            wt3=wlocz(pz)
!        ...store 3D quadrature point
            xip(1)=xi1
            xip(2)=xi2
            xip(3)=xi3
!        ...Shape function subroutine is called only once, when
!           px=py=1 and stored in sH3p(:,pz) and dsH3p(:,pz)
            if (px*py.eq.1) then
               call shape1HH(xi3,nord3,nrdofH3,sH3p(:,pz),dsH3p(:,pz))
            endif
!        ...Copy shape functions in coord. 3 previously evaluated
            shapH3(:,1)=sH3p(:,pz)
            shapH3(:,2)=dsH3p(:,pz)
!
!        ...Compute shape functions needed for geometry - 3D H1 shape functions
            call shape3DH(etype,xip,norder,norient_edge,norient_face, nrdof,shapH,gradH)
#if DEBUG_MODE
            if (nrdof .ne. NrdofH) then
               write(*,*) 'elem_maxwell_fi_hexa: INCONSISTENCY NrdofH. stop.'
               stop
            endif
#endif
!        ...Geometry map
            call geom3D(Mdle,xip,xnod,shapH,gradH,NrdofH, x,dxdxi,dxidx,rjac,iflag)
#if DEBUG_MODE
            if (iflag.ne.0) then
               write(*,5999) Mdle,rjac
 5999          format('elem_maxwell_fi_hexa: Negative Jacobian. Mdle,rjac=',i8,2x,e12.5)
               stop
            endif
#endif
!
!        ...get the RHS
            call getf(Mdle,x, zfval,zJ)
!
!        ...set auxiliary constants (updated below if nonlinear)
!        ...update bgpol depending on coordinates if refractive index is varying
!           (other than step-index)
            if (GEOM_NO .eq. 5) then
               call get_bgPol(dom_flag,Fld_flag,0.d0,x, bg_pol)
            endif
            za = (ZI*OMEGA*OMEGA_RATIO_FLD*EPSILON+SIGMA)*IDENTITY+bg_pol
            zc = (ZI*OMEGA*OMEGA_RATIO_FLD*MU)*IDENTITY
!
!        ...Nonlinear terms
            if(NONLINEAR_FLAG.eq.1) then
!           ...compute current solution using soleval
               nflag = 1
               call soleval(Mdle,xip,norient_edge,norient_face,norder,xnod,&
                 zdofH,zdofE,zdofV,zdofQ,nflag,x,dxdxi, &
                 zsolH_soleval,zdsolH_soleval,zsolE_soleval,zcurlE_soleval,&
                 zsolV_soleval,zdivV_soleval,zsolQ_soleval)
!           ...initialize material refractive index perturbation
               delta_n = 0.d0
               if(HEAT_FLAG.eq.1) then
!              ...compute thermally induced perturbation to refractive index at x
                  rsolH = real(zsolH_soleval(1))
                  delta_n = THERMO_OPT_COEFF*rsolH
               endif
!
!           ...update background polarization with thermal perturbation
               if(delta_n .ne. 0.d0) then
                  call get_bgPol(dom_flag,Fld_flag,delta_n,x, bg_pol)
               endif
!
!           ...initialize gain polarization, raman polarization
               gain_pol = ZERO; raman_pol = ZERO
!          ...skip nonlinear gain computation if inside PML region
               if ( USE_PML .and. ( (x(3).gt.PML_REGION) .or. &
                                    ( (COPUMP.eq.0).and.(x(3).lt.(ZL-PML_REGION)) ) &
                                  ) &
                  ) goto 190
               if (ACTIVE_GAIN .gt. 0.d0) then
                  if (dom_flag .eq. 1) then ! .and. x(3).le.PML_REGION) then
                     call get_activePol(zsolQ_soleval(1:12),Fld_flag,delta_n, gain_pol)
                  endif
               endif
               if (RAMAN_GAIN .gt. 0.d0) then
!              ...Fld_flag = 1 when in signal element routine
                  if(Fld_flag .eq. 1) then
                     call get_ramanPol(zsolQ_soleval(7:9),zsolQ_soleval(10:12), &
                                    dom_flag,Fld_flag,delta_n, raman_pol)
                  else
                     call get_ramanPol(zsolQ_soleval(1:3),zsolQ_soleval(4:6), &
                                    dom_flag,Fld_flag,delta_n, raman_pol)
                  endif
!           ...endif RAMAN_GAIN
               endif
 190           continue
!           ...update auxiliary constant za
               za = (ZI*OMEGA*OMEGA_RATIO_FLD*EPSILON+SIGMA)*IDENTITY+bg_pol+gain_pol+raman_pol
!        ...endif NONLINEAR_FLAG
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
            else
!           ...get PML function
               call get_Beta(x,Fld_flag, zbeta,zdbeta,zd2beta)
               Jstretch(3,3) = zdbeta
!           ...compute det(J) * J^-1 * J^-T
!              (J is a diagonal matrix)
               invJstretch(3,3) = 1.d0/zdbeta
               call ZGEMM('N', 'N', 3, 3, 3, ZONE, invJstretch, 3, &
                             invJstretch, 3, ZERO, JJstretch, 3)
               detJstretch = zdbeta
               JJstretch = detJstretch*JJstretch
            endif
!
!        ...PML stretching
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
!        ...compute total quadrature weight
            wt123=wt1*wt2*wt3
!        ...compute Jacobian determinant * quadrature weight
            weighthh=wt123*rjac
!        ...Determine D = J^-1 * J^-T. Multiply by appropriate weight
            D(1,1) = dxidx(1,1)**2 + dxidx(1,2)**2 + dxidx(1,3)**2
            D(1,2) = dxidx(1,1)*dxidx(2,1) + &
                     dxidx(1,2)*dxidx(2,2) + &
                     dxidx(1,3)*dxidx(2,3)
            D(1,3) = dxidx(1,1)*dxidx(3,1) + &
                     dxidx(1,2)*dxidx(3,2) + &
                     dxidx(1,3)*dxidx(3,3)
            D(2,1) = D(1,2)
            D(2,2) = dxidx(2,1)**2 + dxidx(2,2)**2 + dxidx(2,3)**2
            D(2,3) = dxidx(2,1)*dxidx(3,1) + &
                     dxidx(2,2)*dxidx(3,2) + &
                     dxidx(2,3)*dxidx(3,3)
            D(3,1) = D(1,3)
            D(3,2) = D(2,3)
            D(3,3) = dxidx(3,1)**2 + dxidx(3,2)**2 + dxidx(3,3)**2
            D = D*weighthh
!        ...compute inverse Jacobian determinant * quadrature weight
            weightvv=wt123/rjac
!        ...Determine C = J^T * J.  Multiply by appropriate weight
            C(1,1) = dxdxi(1,1)**2 + dxdxi(2,1)**2 + dxdxi(3,1)**2
            C(1,2) = dxdxi(1,1)*dxdxi(1,2) + &
                     dxdxi(2,1)*dxdxi(2,2) + &
                     dxdxi(3,1)*dxdxi(3,2)
            C(1,3) = dxdxi(1,1)*dxdxi(1,3) + &
                     dxdxi(2,1)*dxdxi(2,3) + &
                     dxdxi(3,1)*dxdxi(3,3)
            C(2,1) = C(1,2)
            C(2,2) = dxdxi(1,2)**2 + dxdxi(2,2)**2 + dxdxi(3,2)**2
            C(2,3) = dxdxi(1,2)*dxdxi(1,3) + &
                     dxdxi(2,2)*dxdxi(2,3) + &
                     dxdxi(3,2)*dxdxi(3,3)
            C(3,1) = C(1,3)
            C(3,2) = C(2,3)
            C(3,3) = dxdxi(1,3)**2 + dxdxi(2,3)**2 + dxdxi(3,3)**2
            C = C*weightvv
!        ...Determine auxiliary matrices to anisotropic refractive index tensor
!           (also needed for PML stretching in adjoint graph norm)
!        ...D_za = J^-1 |zaJ|^2 J^-T
            D_aux(1:3,1) = dxidx(1:3,1) * abs(zaJ(1,1))**2
            D_aux(1:3,2) = dxidx(1:3,2) * abs(zaJ(2,2))**2
            D_aux(1:3,3) = dxidx(1:3,3) * abs(zaJ(3,3))**2
            call DGEMM('N','T',3,3,3,1.0d0,D_aux,3,dxidx,3,0.0d0,D_za,3)
            D_za = D_za*weighthh
!        ...D_zc = J^-1 |zcJ|^2 J^-T
            D_aux(1:3,1) = dxidx(1:3,1) * abs(zcJ(1,1))**2
            D_aux(1:3,2) = dxidx(1:3,2) * abs(zcJ(2,2))**2
            D_aux(1:3,3) = dxidx(1:3,3) * abs(zcJ(3,3))**2
            call DGEMM('N','T',3,3,3,1.0d0,D_aux,3,dxidx,3,0.0d0,D_zc,3)
            D_zc = D_zc*weighthh
!        ...Z_za = J^-1 * zaJ^T * J
            call ZGEMM('T','N',3,3,3,ZONE,zaJ        ,3,ZONE*dxdxi,3,ZERO,Z_aux,3)
            call ZGEMM('N','N',3,3,3,ZONE,ZONE*dxidx ,3,     Z_aux,3,ZERO,Z_za ,3)
!        ...Z_zc = J^-1 * zcJ^T * J
            call ZGEMM('T','N',3,3,3,ZONE,zcJ        ,3,ZONE*dxdxi,3,ZERO,Z_aux,3)
            call ZGEMM('N','N',3,3,3,ZONE,ZONE*dxidx ,3,     Z_aux,3,ZERO,Z_zc ,3)
!        ...Z_zb = conjg(transpose(Z_za)) = J^T * zaJ * J^-T
!        ...Z_zd = conjg(transpose(Z_zc)) = J^T * zcJ * J^-T
!
!        ...put appropriate quadrature weight on Jacobian and its inverse
            dxdxi = dxdxi * weightvv
            dxidx = dxidx * wt123
!
!--------------------------------------------------------------------------
!        ...SUM FACTORIZATION LOOPS START
!--------------------------------------------------------------------------
!        ... For Gram, Stiffness, Load:
!        ... Hcurl TEST  functions will be identified by indices i1,i2,i3,a
!
!        ... For Gram matrix:
!        ... Hcurl TRIAL functions will be identified by indices j1,j2,j3,b
!
!        ... For stiffness matrix:
!        ... L2    TRIAL functions will be identified by indices j1,j2,j3,b
!--------------------------------------------------------------------------
!
!        ...Integration of innermost integrals for Gram matrix
!
!        ...loop over 1D test function, coord 3, Hcurl
            do i3=1,nrdofH3
!           ...loop over 1D DOFs, coord. 3, trial shape func Hcurl
               do j3=1,nrdofH3
!              ... combine indices i3 and j3 into k3
                  k3=(i3-1)*nrdofH3+j3
!              ...loops on vector components a (for test function), b (for trial)
!              ...The way the 3D Hcurl broken functions are organized, allow for
!              ...computing only the upper-triangular blocks of matrices (b>=a)
                  do a=1,3
                     do b=a,3
!                    ...indices sb and sa are 1 or 2, depending on a and b
                        sb=1+deltak(b,3)
                        sa=1+deltak(a,3)
!                    ...indices idxb, idxa, are equal to j3,i3, shifted by 0 or 1,
!                       depending on whether we need the H1 value (j3+0,i3+0)
!                       or the L2 (j3+1,i3+1), which is determined
!                       by the vector component the shape function lies at
                        idxb=j3+deltak(b,3)
                        idxa=i3+deltak(a,3)
!                    ...we need to check that idxb, idxa don't exceed dimension of
!                       1D H1 enriched test space
                        if ((idxb.le.nrdofH3).and.(idxa.le.nrdofH3)) then
!                       ...accumulate innermost 1D integral for EE term in Gram matrix
                           AUXEE_A_zb(b,a,k3) = AUXEE_A_zb(b,a,k3)         &
                                              + (ALPHA_NORM*D(a,b) +       &
                                                 D_za(a,b)           )     &
                                              * shapH3(idxa,sa)            &
                                              * shapH3(idxb,sb)
!
                           AUXEE_A_zc(b,a,k3) = AUXEE_A_zc(b,a,k3)         &
                                              + (ALPHA_NORM*D(a,b) +       &
                                                 D_zc(a,b)           )     &
                                              * shapH3(idxa,sa)            &
                                              * shapH3(idxb,sb)
!
!                       ...loop over components a+alph, b+beta, (modulo 3)
!                          where the curl of the shape functions lie
                           do beta=1,2; do alph=1,2
!                          ...compute a+alph, b+beta, (modulo 3)
                              idxbeta=mod(b+beta-1,3)+1
                              idxalph=mod(a+alph-1,3)+1
!                          ...determine indices sb,sa for the curl components
                              sb=1+1-deltak(idxbeta,3)
                              sa=1+1-deltak(idxalph,3)
!                          ...accumulate innermost 1D integral for CC term in Gram matrix
                              AUXCC_A(alph,beta,b,a,k3) = AUXCC_A(alph,beta,b,a,k3) &
                                                        + shapH3(idxa,sa)           &
                                                        * shapH3(idxb,sb)           &
                                                        * (-1)**(alph+beta)         &
                                                        * C(idxalph,idxbeta)
                           enddo; enddo
!
!                       ...loop over components a+alph, where the curl of
!                          the TEST shape function, for the CE term of Gram matrix
                           do alph=1,2
                              idxalph=mod(a+alph-1,3)+1
                              sb=1+deltak(b,3)
                              sa=1+1-deltak(idxalph,3)
!                          ...accumulate innermost 1D integral for CE term
!                             (Z_za and Z_zc indices are switched here b/c we need the transpose)
                              AUXCE_A_zb(alph,b,a,k3) = AUXCE_A_zb(alph,b,a,k3)     &
                                                      - conjg(Z_za(b,idxalph))      &
                                                      * shapH3(idxa,sa)             &
                                                      * shapH3(idxb,sb)             &
                                                      * (-1)**(alph-1)*wt123
!
                              AUXCE_A_zc(alph,b,a,k3) = AUXCE_A_zc(alph,b,a,k3)     &
                                                      + conjg(Z_zc(b,idxalph))      &
                                                      * shapH3(idxa,sa)             &
                                                      * shapH3(idxb,sb)             &
                                                      * (-1)**(alph-1)*wt123
                           enddo
!                       ...loop over components b+beta, where the curl of
!                          the TRIAL shape function, for the EC term of Gram matrix
                           do beta=1,2
                              idxbeta=mod(b+beta-1,3)+1
                              sb=1+1-deltak(idxbeta,3)
                              sa=1+deltak(a,3)
!                          ...accumulate innermost 1D integral for EC term
                              AUXEC_A_zb(beta,b,a,k3) = AUXEC_A_zb(beta,b,a,k3)        &
                                                     - Z_za(a,idxbeta)                 &
                                                     * shapH3(idxa,sa)                 &
                                                     * shapH3(idxb,sb)                 &
                                                     * (-1)**(beta-1)*wt123
!
                              AUXEC_A_zc(beta,b,a,k3) = AUXEC_A_zc(beta,b,a,k3)        &
                                                     + Z_zc(a,idxbeta)                 &
                                                     * shapH3(idxa,sa)                 &
                                                     * shapH3(idxb,sb)                 &
                                                     * (-1)**(beta-1)*wt123
!
!                          ...loop over beta ends
                           enddo
                        endif
!                    ...loop over b ends
                     enddo
!                 ...loop over a ends
                  enddo
!              ...loop over j3 ends
               enddo
!           ...loop over i3 ends
            enddo
!
!        ...Integration of enriched stiffness matrix and load vector
!
!        ...start new loop over i3 - Hcurl test functions
            do i3=1,nrdofH3
!           ...start new loop over j3 - L2 trial functions
               do j3=1,nrdofQ3_tr
!              ...combine indices i3, j3 into k3
                  k3=(i3-1)*nrdofQ3_tr+j3
!              ...loop over vector components
                  do b=1,3
                     do a=1,3
!                    ...determine indices for test function
                        sa=1+deltak(a,3)
                        idxa=i3+deltak(a,3)
!                    ...check that idxa don't exceed dimension
!                       of 1D H1 enriched test space
                        if (idxa.le.nrdofH3) then
!                       ...accumulate innermost 1D integral for QE term
!                       ...and account for PML stretching
                           STIFQE_A(a,b,k3) = STIFQE_A(a,b,k3)             &
                                            + (dxidx(a,b)*shapH3(j3+1,2))  &
                                            * zcJ(b,b)                     &
                                            * shapH3(idxa,sa)
!                       ...accumulate innermost 1D integral for alpha_scale*QE term
!                       ...and account for PML stretching
                           STIFQE_ALPHA_A(a,b,k3) = STIFQE_ALPHA_A(a,b,k3)       &
                                                  - (dxidx(a,b)*shapH3(j3+1,2))  &
                                                  * zaJ(b,b)                     &
                                                  * shapH3(idxa,sa)
!                       ...loop over components a+alph, required for curl of test f
                           do alph=1,2
                              idxalph=mod(a+alph-1,3)+1
                              sa=1+1-deltak(idxalph,3)
!                          ...accumulate innermost 1D integral for QC term
                              STIFQC_A(alph,a,b,k3) = STIFQC_A(alph,a,b,k3)           &
                                                    + shapH3(j3+1,2)*dxdxi(b,idxalph) &
                                                    * (-1)**(alph-1)*shapH3(idxa,sa)
                           enddo
                        endif
!                    ...loop over a ends
                     enddo
!                 ...loop over b ends
                  enddo
!              ...loop over j3 ends
               enddo
!           ...loop over vector component a of Hcurl test functions
               do a=1,3
                  sa=1+deltak(a,3)
                  idxa=i3+deltak(a,3)
                  if (idxa.le.nrdofH3) then
!                 ...accumulate innermost integral of load vector
                     LOADE_A(a,i3) = LOADE_A(a,i3) &
                                   + ( dxidx(a,1)*zJ(1)+    &
                                       dxidx(a,2)*zJ(2)+    &
                                       dxidx(a,3)*zJ(3) )   &
                                   * shapH3(idxa,sa)*rjac
                  endif
!              ...loop over a ends
               enddo
!           ...loop over i3 ends
            enddo
!        ...loop over pz ends
         enddo
!
!
!     ...Computation of middle 1D integrals
!
!     ...Compute middle integrals in Gram matrix terms
!
!     ...loop over Hcurl i3 test function
         do i3=1,nrdofH3
!        ...loop over Hcurl j3 trial function
            do j3=1,nrdofH3
!           ...combine i3 and j3 into k3
               k3=(i3-1)*nrdofH3+j3
!           ...loop over Hcurl i3 test function
               do i2=1,nrdofH2
!              ...loop over Hcurl j3 trial function
                  do j2=1,nrdofH2
!                 ...combine i2 and j2 into k2
                     k2=(i2-1)*nrdofH2+j2
!                 ...loop over vector components a (test), b(trial)
!                 ...Only upper blocks of matrix ( b>=a )are computed
                     do a=1,3
                        do b=a,3
!                       ...determine indices for shape functions
                           sb=1+deltak(b,2)
                           sa=1+deltak(a,2)
                           idxb=j2+deltak(b,2)
                           idxa=i2+deltak(a,2)
!                       ...check that dimensions are not surpassed by idxa,idxb
                           if((idxb.le.nrdofH2).and.(idxa.le.nrdofH2)) then
!                          ...accumulate middle 1D integral of EE term in Gram
                              AUXEE_B_zb(b,a,k2,k3) = AUXEE_B_zb(b,a,k2,k3) &
                                                    + shapH2(idxa,sa)       &
                                                    * shapH2(idxb,sb)       &
                                                    * AUXEE_A_zb(b,a,k3)
                              AUXEE_B_zc(b,a,k2,k3) = AUXEE_B_zc(b,a,k2,k3) &
                                                    + shapH2(idxa,sa)       &
                                                    * shapH2(idxb,sb)       &
                                                    * AUXEE_A_zc(b,a,k3)
!                          ...loop over b+beta, a+alph, for curl of trial and test
                              do beta=1,2;do alph=1,2
                                 idxbeta=mod(b+beta-1,3)+1
                                 idxalph=mod(a+alph-1,3)+1
                                 sb=1+1-deltak(idxbeta,2)
                                 sa=1+1-deltak(idxalph,2)
!                             ...accumulate middle 1D integral of CC term in Gram
                                 AUXCC_B(alph,beta,b,a,k2,k3) =            &
                                         AUXCC_B(alph,beta,b,a,k2,k3)      &
                                       + shapH2(idxa,sa)*shapH2(idxb,sb)   &
                                       * AUXCC_A(alph,beta,b,a,k3)
                              enddo; enddo
!                          ... accumulate for CE term
                              do alph=1,2
                                 idxalph=mod(a+alph-1,3)+1
                                 sb=1+deltak(b,2)
                                 sa=1+1-deltak(idxalph,2)
                                 AUXCE_B_zb(alph,b,a,k2,k3) =            &
                                        AUXCE_B_zb(alph,b,a,k2,k3)       &
                                      + shapH2(idxa,sa)*shapH2(idxb,sb)  &
                                      * AUXCE_A_zb(alph,b,a,k3)
                                 AUXCE_B_zc(alph,b,a,k2,k3) =            &
                                        AUXCE_B_zc(alph,b,a,k2,k3)       &
                                      + shapH2(idxa,sa)*shapH2(idxb,sb)  &
                                      * AUXCE_A_zc(alph,b,a,k3)
                              enddo
!                          ... accumulate for EC term
                              do beta=1,2
                                 idxbeta=mod(b+beta-1,3)+1
                                 sb=1+1-deltak(idxbeta,2)
                                 sa=1+deltak(a,2)
                                 AUXEC_B_zb(beta,b,a,k2,k3) =            &
                                        AUXEC_B_zb(beta,b,a,k2,k3)       &
                                      + shapH2(idxa,sa)*shapH2(idxb,sb)  &
                                      * AUXEC_A_zb(beta,b,a,k3)
                                 AUXEC_B_zc(beta,b,a,k2,k3) =            &
                                        AUXEC_B_zc(beta,b,a,k2,k3)       &
                                      + shapH2(idxa,sa)*shapH2(idxb,sb)  &
                                      * AUXEC_A_zc(beta,b,a,k3)
                              enddo
                           endif
!                       ...loop over b ends
                        enddo
!                    ...loop over a ends
                     enddo
!                 ...loop over j2 ends
                  enddo
!              ...loop over i2 ends
               enddo
!           ...loop over j3 ends
            enddo
!        ...loop over i3 ends
         enddo
!
!     ...Compute middle integrals for stiffness terms
!
         do i3=1,nrdofH3
            do j3=1,nrdofQ3_tr
               k3=(i3-1)*nrdofQ3_tr+j3
               do i2=1,nrdofH2
                  do j2=1,nrdofQ2_tr
                     k2=(i2-1)*nrdofQ2_tr+j2
                     do b=1,3
                        do a=1,3
                           sa=1+deltak(a,2)
                           idxa=i2+deltak(a,2)
                           if (idxa.le.nrdofH2) then
!                          ...accumulate zc*QE term
                              STIFQE_B(a,b,k2,k3) = STIFQE_B(a,b,k2,k3)    &
                                                  + shapH2(j2+1,2)         &
                                                  * shapH2(idxa,sa)        &
                                                  * STIFQE_A(a,b,k3)
!                          ...accumulate za*QE term
                              STIFQE_ALPHA_B(a,b,k2,k3) = STIFQE_ALPHA_B(a,b,k2,k3) &
                                                        + shapH2(j2+1,2)            &
                                                        * shapH2(idxa,sa)           &
                                                        * STIFQE_ALPHA_A(a,b,k3)
!                          ...accumulate QC term
                              do alph=1,2
                                 idxalph=mod(a+alph-1,3)+1
                                 sa=1+1-deltak(idxalph,2)
                                 STIFQC_B(alph,a,b,k2,k3) = STIFQC_B(alph,a,b,k2,k3)   &
                                                          + shapH2(j2+1,2)             &
                                                          * shapH2(idxa,sa)            &
                                                          * STIFQC_A(alph,a,b,k3)
                              enddo
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
!     ...Compute middle integral of load vector
         do i3=1,nrdofH3
            do i2=1,nrdofH2
               do a=1,3
                  sa=1+deltak(a,2)
                  idxa=i2+deltak(a,2)
                  if (idxa.le.nrdofH2) then
                     LOADE_B(a,i2,i3) = LOADE_B(a,i2,i3)    &
                                      + shapH2(idxa,sa)     &
                                      * LOADE_A(a,i3)
                  endif
               enddo
            enddo
         enddo
!     ...loop over py ends
      enddo
!
! FINAL COMPUTATION OF GRAM MATRIX
!
!  ...loop over Hcurl trial shape function identified by indices j1,j2,j3,b
      do b=1,3
         do j3=1,nrdofH3
            do j2=1,nrdofH2
               do j1=1,nrdofH1
                  sb=1+deltak(b,1)
                  idxb=j1+deltak(b,1)
                  idxb2=j2+deltak(b,2)
                  idxb3=j3+deltak(b,3)
!              ...determine index of 3D Hcurl trial shape function using j1,j2,j3,b
                  if ((idxb.le.nrdofH1).and.(idxb2.le.nrdofH2).and.(idxb3.le.nrdofH3)) then
                     select case(b)
                        case(1)
                        m2=j1+nord1*(j2-1)+nord1*nrdofH2*(j3-1)
                        case(2)
                        m2=nord1*nrdofH2*nrdofH3 &
                           +j1+nrdofH1*(j2-1)+nrdofH1*nord2*(j3-1)
                        case(3)
                        m2=nord1*nrdofH2*nrdofH3 &
                           +nrdofH1*nord2*nrdofH3 &
                           +j1+nrdofH1*(j2-1)+nrdofH1*nrdofH2*(j3-1)
                     end select
!                 ...loop over Hcurl test shape function identified by indices i1,i2,i3,a
                     do a=1,3
                        do i3=1,nrdofH3
                           do i2=1,nrdofH2
                              do i1=1,nrdofH1
                                 sa=1+deltak(a,1)
                                 idxa=i1+deltak(a,1)
                                 idxa2=i2+deltak(a,2)
                                 idxa3=i3+deltak(a,3)
                                 if ((idxa.le.nrdofH1).and.(idxa2.le.nrdofH2).and.(idxa3.le.nrdofH3)) then
!                                ...combine indices i3, j3 into k3
                                    k3=(i3-1)*nrdofH3+j3
!                                ...combine indices i2, j2 into k2
                                    k2=(i2-1)*nrdofH2+j2
!                                ...determine index of 3D Hcurl trial shape function
                                    select case(a)
                                       case(1)
                                       m1=i1+nord1*(i2-1)+nord1*nrdofH2*(i3-1)
                                       case(2)
                                       m1=nord1*nrdofH2*nrdofH3 &
                                         +i1+nrdofH1*(i2-1)+nrdofH1*nord2*(i3-1)
                                       case(3)
                                       m1=nord1*nrdofH2*nrdofH3+nrdofH1*nord2*nrdofH3 &
                                         +i1+nrdofH1*(i2-1)+nrdofH1*nrdofH2*(i3-1)
                                    end select
!                                ...accumulate integrals only if m1<=m2
                                    if (m1.le.m2) then
                                       sa=1+deltak(a,1)
                                       sb=1+deltak(b,1)
                                       kk = nk(2*m1-1,2*m2-1)
!                                   ...sum EE terms
                                       gramP(kk) = gramP(kk)         &
                                                 + shapH1(idxa,sa)   &
                                                 * shapH1(idxb,sb)   &
                                                 * AUXEE_B_zb(b,a,k2,k3)
!                                   ...sum CC terms
                                       do beta=1,2; do alph=1,2
                                          idxbeta=mod(b+beta-1,3)+1
                                          idxalph=mod(a+alph-1,3)+1
                                          sb=1+1-deltak(idxbeta,1)
                                          sa=1+1-deltak(idxalph,1)
                                          gramP(kk) = gramP(kk)          &
                                                    + shapH1(idxa,sa)    &
                                                    * shapH1(idxb,sb)    &
                                                    * AUXCC_B(alph,beta,b,a,k2,k3)

                                       enddo; enddo

                                       kk = nk(2*m1-1,2*m2)
!                                   ...sum CE terms
                                       do alph=1,2
                                          idxalph=mod(a+alph-1,3)+1
                                          sb=1+deltak(b,1)
                                          sa=1+1-deltak(idxalph,1)
                                          gramP(kk) = gramP(kk)                   &
                                                    + AUXCE_B_zc(alph,b,a,k2,k3)  &
                                                    * shapH1(idxa,sa)*shapH1(idxb,sb)
                                       enddo
!                                   ...sum EC terms
                                       do beta=1,2
                                          idxbeta=mod(b+beta-1,3)+1
                                          sb=1+1-deltak(idxbeta,1)
                                          sa=1+deltak(a,1)
                                          gramP(kk) = gramP(kk)                     &
                                                    + AUXEC_B_zb(beta,b,a,k2,k3)    &
                                                    * shapH1(idxa,sa)*shapH1(idxb,sb)
                                       enddo

                                       if (m1.ne.m2) then

                                          kk = nk(2*m1  ,2*m2-1)
!                                      ...sum CE terms
                                          do alph=1,2
                                             idxalph=mod(a+alph-1,3)+1
                                             sb=1+deltak(b,1)
                                             sa=1+1-deltak(idxalph,1)
                                             gramP(kk) = gramP(kk)                     &
                                                       + AUXCE_B_zb(alph,b,a,k2,k3)    &
                                                       * shapH1(idxa,sa)*shapH1(idxb,sb)
                                          enddo
!                                      ...sum EC terms
                                          do beta=1,2
                                             idxbeta=mod(b+beta-1,3)+1
                                             sb=1+1-deltak(idxbeta,1)
                                             sa=1+deltak(a,1)
                                             gramP(kk) = gramP(kk)                   &
                                                       + AUXEC_B_zc(beta,b,a,k2,k3)  &
                                                       * shapH1(idxa,sa)*shapH1(idxb,sb)
                                          enddo

                                       endif

                                       kk = nk(2*m1  ,2*m2  )
!                                   ...sum EE terms
                                       sb=1+deltak(b,1)
                                       sa=1+deltak(a,1)
                                       gramP(kk) = gramP(kk)         &
                                                 + shapH1(idxa,sa)   &
                                                 * shapH1(idxb,sb)   &
                                                 * AUXEE_B_zc(b,a,k2,k3)
!                                   ...sum CC terms
                                       do beta=1,2; do alph=1,2
                                          idxbeta=mod(b+beta-1,3)+1
                                          idxalph=mod(a+alph-1,3)+1
                                          sb=1+1-deltak(idxbeta,1)
                                          sa=1+1-deltak(idxalph,1)
                                          gramP(kk) = gramP(kk)         &
                                                    + shapH1(idxa,sa)   &
                                                    * shapH1(idxb,sb)   &
                                                    * AUXCC_B(alph,beta,b,a,k2,k3)

                                       enddo; enddo
                                     endif
                                 endif
                              enddo
                           enddo
                        enddo
                     enddo
                  endif
               enddo
            enddo
         enddo
      enddo
!
!  ...Compute final integral for Stiffness matrix
!
!  ...loop over trial L2 shape functions
      do j3=1,nrdofQ3_tr
         do j2=1,nrdofQ2_tr
            do j1=1,nrdofQ1_tr
!           ...determine 3D L2 shape function index
               m2=j1+(j2-1)*nrdofQ1_tr+(j3-1)*nrdofQ1_tr*nrdofQ2_tr
!           ...loop over Hcurl test shape function identified by indices i1,i2,i3,a
               do a=1,3
                  do i3=1,nrdofH3-deltak(a,3)
                     do i2=1,nrdofH2-deltak(a,2)
                        do i1=1,nrdofH1-deltak(a,1)
                           sa=1+deltak(a,1)
                           idxa=i1+deltak(a,1)
                           idxa2=i2+deltak(a,2)
                           idxa3=i3+deltak(a,3)
                           if ((idxa.le.nrdofH1).and.(idxa2.le.nrdofH2).and.(idxa3.le.nrdofH3)) then
!                          ...combine indices i3, j3 into k3
                              k3=(i3-1)*nrdofQ3_tr+j3
!                          ...combine indices i2, j2 into k2
                              k2=(i2-1)*nrdofQ2_tr+j2
!                          ...determine index for 3D Hcurl test function
                              select case(a)
                                 case(1)
                                 m1=i1+nord1*(i2-1)+nord1*nrdofH2*(i3-1)
                                 case(2)
                                 m1=nord1*nrdofH2*nrdofH3 &
                                   +i1+nrdofH1*(i2-1)+nrdofH1*nord2*(i3-1)
                                 case(3)
                                 m1=nord1*nrdofH2*nrdofH3+nrdofH1*nord2*nrdofH3 &
                                   +i1+nrdofH1*(i2-1)+nrdofH1*nrdofH2*(i3-1)
                              end select
!                          ...loop over vector components of L2 trial functions
                              do b=1,3
                                 k = (m2-1)*6 + 3+ b
!                             ...accumulate QC term
                                 do alph=1,2
                                    idxalph=mod(a+alph-1,3)+1
                                    sa=1+1-deltak(idxalph,1)
                                    stiff_EQ_T(k,2*m1-1) = stiff_EQ_T(k,2*m1-1) &
                                                      + shapH1(idxa,sa)*shapH1(j1+1,2) &
                                                      * STIFQC_B(alph,a,b,k2,k3)
                                 enddo

                                 sa=1+deltak(a,1)

                                 k = (m2-1)*6 + b
!                             ...accumulate QE term
                                 stiff_EQ_T(k,2*m1-1) = stiff_EQ_T(k,2*m1-1) &
                                                   +shapH1(idxa,sa)*shapH1(j1+1,2) &
                                                   *STIFQE_ALPHA_B(a,b,k2,k3)
                                 k = (m2-1)*6 + 3+ b
!                             ...accumulate QE term
                                 stiff_EQ_T(k,2*m1) = stiff_EQ_T(k,2*m1) &
                                                   +shapH1(idxa,sa)*shapH1(j1+1,2) &
                                                   *STIFQE_B(a,b,k2,k3)
                                 k = (m2-1)*6 + b
!                             ...accumulate QC term
                                 do alph=1,2
                                    idxalph=mod(a+alph-1,3)+1
                                    sa=1+1-deltak(idxalph,1)
                                    stiff_EQ_T(k,2*m1) = stiff_EQ_T(k,2*m1) &
                                                      + shapH1(idxa,sa)*shapH1(j1+1,2) &
                                                      * STIFQC_B(alph,a,b,k2,k3)
                                 enddo
                              enddo
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
!     Final computation of load vector
!
!  ...loop over Hcurl test shape function identified by indices i1,i2,i3,a
      do a=1,3
         do i3=1,nrdofH3
            do i2=1,nrdofH2
               do i1=1,nrdofH1
                  sa=1+deltak(a,1)
                  idxa=i1+deltak(a,1)
                  idxa2=i2+deltak(a,2)
                  idxa3=i3+deltak(a,3)
                  if ((idxa.le.nrdofH1).and.(idxa2.le.nrdofH2).and.(idxa3.le.nrdofH3)) then
!                 ...combine indices i3, j3 into k3
                     k3=(i3-1)*nrdofQ3_tr+j3
!                 ...combine indices i2, j2 into k2
                     k2=(i2-1)*nrdofQ2_tr+j2
!                 ...determine index for 3D Hcurl test function
                     select case(a)
                        case(1)
                        m1=i1+nord1*(i2-1)+nord1*nrdofH2*(i3-1)
                        case(2)
                        m1=nord1*nrdofH2*nrdofH3 &
                          +i1+nrdofH1*(i2-1)+nrdofH1*nord2*(i3-1)
                        case(3)
                        m1=nord1*nrdofH2*nrdofH3+nrdofH1*nord2*nrdofH3 &
                          +i1+nrdofH1*(i2-1)+nrdofH1*nrdofH2*(i3-1)
                     end select
!                 ...accumulate outermost integral
                     bload_E(2*m1-1)=bload_E(2*m1-1)+shapH1(idxa,sa)*LOADE_B(a,i2,i3)
                  endif
               enddo
            enddo
         enddo
      enddo
!  ...loop over px ends
   enddo
!
!..end timer
   end_time = MPI_Wtime()
   !$OMP CRITICAL
      !write(*,10) etype, end_time-start_time
      write(*,11) end_time-start_time
! 10   format(A,' elem : ',f12.5,'  seconds')
 11   format(f12.5)
   !$OMP END CRITICAL
!
   deallocate(AUXEE_A_zb,AUXEE_A_zc)
   deallocate(AUXEE_B_zb,AUXEE_B_zc)
   deallocate(AUXCC_A)
   deallocate(AUXCC_B)
   deallocate(AUXEC_A_zb,AUXEC_A_zc)
   deallocate(AUXEC_B_zb,AUXEC_B_zc)
   deallocate(AUXCE_A_zb,AUXCE_A_zc)
   deallocate(AUXCE_B_zb,AUXCE_B_zc)
!
   deallocate(STIFQE_A)
   deallocate(STIFQE_B)
   deallocate(STIFQE_ALPHA_A)
   deallocate(STIFQE_ALPHA_B)
   deallocate(STIFQC_A)
   deallocate(STIFQC_B)
   deallocate(LOADE_A)
   deallocate(LOADE_B)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                 B O U N D A R Y   I N T E G R A L S
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!..compute index array to skip bubble shape functions in interface variable
!..find ndof associated with the mdle node of the element
!
   call ndof_nod(etype,nordP, ndofHHmdl,ndofEEmdl,ndofVVmdl,ndofQQmdl)
!
   nrdofEEi = NrdofEE-ndofEEmdl
!
   allocate(idxEE(nrdofEEi))
   ik = 0
   do i3 = 1,nrdofH3
      do i2 = 1,nrdofH2
         do i1 = 1,nord1
            if (i2 .lt.3 .or. i3 .lt.3) then
               ik = ik + 1
               m1 = i1+nord1*(i2-1) + nord1*nrdofH2*(i3-1)
               idxEE(ik) = m1
            endif
         enddo
      enddo
   enddo
   do i3 = 1,nrdofH3
      do i2 = 1,nord2
         do i1 = 1,nrdofH1
            if (i1 .lt.3 .or. i3 .lt.3) then
               ik = ik + 1
               m1 = nord1*nrdofH2*nrdofH3 + i1 + nrdofH1*(i2-1) + nrdofH1*nord2*(i3-1)
               idxEE(ik) = m1
            endif
         enddo
      enddo
   enddo
   do i3 = 1,nord3
      do i2 = 1,nrdofH2
         do i1 = 1,nrdofH1
            if (i1 .lt.3 .or. i2 .lt.3) then
               ik = ik + 1
               m1 = nrdofH1*nord2*nrdofH3 + nord1*nrdofH2*nrdofH3   &
                  + i1+nrdofH1*(i2-1)+nrdofH1*nrdofH2*(i3-1)
               idxEE(ik) = m1
            endif
         enddo
      enddo
   enddo
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
#if DEBUG_MODE
         if (nrdof .ne. NrdofEE) then
            write(*,*) 'elem_maxwell_fi_hexa: INCONSISTENCY NrdofEE. stop.'
            stop
         endif
#endif
!
!     ...determine element H1 shape functions (for geometry)
         call shape3DH(etype,xi,norder,norient_edge,norient_face, &
                       nrdof,shapH,gradH)
#if DEBUG_MODE
         if (nrdof .ne. NrdofH) then
            write(*,*) 'elem_maxwell_fi_hexa: INCONSISTENCY NrdofH. stop.'
            stop
         endif
#endif
!
!     ...determine element H(curl) shape functions (for fluxes)
!     ...for interfaces only (no bubbles)
         call shape3DE(etype,xi,norderi,norient_edge,norient_face, &
                       nrdof,shapE,curlE)
#if DEBUG_MODE
         if (nrdof .ne. NrdofEi) then
            write(*,*) 'elem_maxwell_fi_hexa: INCONSISTENCY NrdofEi. stop.'
            stop
         endif
#endif
!
!     ...geometry
         call bgeom3D(Mdle,xi,xnod,shapH,gradH,NrdofH,dxidt,nsign, &
                      x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
         weight = bjac*wtloc(l)

!     ...loop through Hcurl enriched test functions
         do ik=1,nrdofEEi
         !do k1=1,nrdofEE
            k1 = idxEE(ik)
            E1(1:3) = shapEE(1,k1)*dxidx(1,1:3) &
                    + shapEE(2,k1)*dxidx(2,1:3) &
                    + shapEE(3,k1)*dxidx(3,1:3)
!        ...check for impedance BC
!           (impedance constant is GAMMA for TE10 mode in rectangular waveguide)
!           ( < n x H , G > = GAMMA*< n x n x E , G > + < zg , G > )
            if ((ibc(ifc,2).eq.9 .and. Fld_flag .eq. 1) .or. &
                (ibc(ifc,3).eq.9 .and. Fld_flag .eq. 0)) then
!           ...get the boundary source
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
!           ...check for impedance bc
               if ((ibc(ifc,2).eq.9 .and. Fld_flag .eq. 1) .or. &
                   (ibc(ifc,3).eq.9 .and. Fld_flag .eq. 0)) then
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
   deallocate(idxEE)
!
!----------------------------------------------------------------------
!      Construction of the DPG system
!----------------------------------------------------------------------
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
      write(*,*) 'elem_maxwell_fi_hexa: ZPPTRF: Mdle,info = ',Mdle,info,'. stop.'
      stop
   endif
!
!..B. Solve triangular system to obtain B~, (LX=) U^*X = [B|l]
   call ZTPTRS('U','C','N',NrTest,NrTrial+1,gramP,stiff_ALL,NrTest,info)
   if (info.ne.0) then
      write(*,*) 'elem_maxwell_fi_hexa: ZTPTRS: Mdle,info = ',Mdle,info,'. stop.'
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
end subroutine elem_maxwell_fi_hexa

