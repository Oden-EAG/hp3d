!
!--------------------------------------------------------------------
!
!    routine name      - elem_maxwell_fi_pris
!
!--------------------------------------------------------------------
!
!     latest revision:  - Oct 2019
!
!     purpose:          - routine returns unconstrained (ordinary)
!                         stiffness matrix and load vector for the
!                         UW DPG formulation for Maxwell equations
!                       - Uses sum factorization for fast integration
!                         with swapped loops
!                       - for prismatic and hexahedral elements only
!
!     ASSUMES FOR HEXAS ISOTROPIC ORDER IN X,Y
!
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
#include "implicit_none.h"
!
subroutine elem_maxwell_fi_pris(Mdle,Fld_flag,                &
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
   integer,                        intent(in)  :: Mdle
   integer,                        intent(in)  :: Fld_flag
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
   real*8, parameter :: rZero = 0.d0
!
!..declare edge/face type variables
   character(len=4) :: etype, etype1, ftype
!
!..declare element order, orientation for edges and faces
   integer, dimension(19)    :: norder
   integer, dimension(12)    :: norient_edge
   integer, dimension(6)     :: norient_face
   integer, dimension(5)     :: norder_f, norder_fe
!
!..element nodes order (trial) for interfaces
   integer, dimension(19)    :: norderi
!
!..face order
   integer, dimension(5)     :: norderf
!
!..geometry dof (work space for nodcor)
   real*8, dimension(3,MAXbrickH) :: xnod
!
!..solution dof (work space for solelm)
   complex(8), dimension(MAXEQNH,MAXbrickH) :: zdofH
   complex(8), dimension(MAXEQNE,MAXbrickE) :: zdofE
   complex(8), dimension(MAXEQNV,MAXbrickV) :: zdofV
   complex(8), dimension(MAXEQNQ,MAXbrickQ) :: zdofQ
!
!..approximate solution -- using soleval
   integer :: nflag
   complex(8), dimension(  MAXEQNH  ) ::  zsolH_soleval
   complex(8), dimension(  MAXEQNH,3) :: zdsolH_soleval
   complex(8), dimension(3,MAXEQNE  ) ::  zsolE_soleval
   complex(8), dimension(3,MAXEQNE  ) :: zcurlE_soleval
   complex(8), dimension(3,MAXEQNV  ) ::  zsolV_soleval
   complex(8), dimension(  MAXEQNV  ) ::  zdivV_soleval
   complex(8), dimension(  MAXEQNQ  ) ::  zsolQ_soleval
   real*8 :: rsolH
!
!..variables for geometry
   real*8, dimension(3)    :: xi,x,rn
   real*8, dimension(3,2)  :: dxidt,dxdt,rt
   real*8, dimension(3,3)  :: dxdxi,dxidx
   real*8, dimension(2)    :: t
!
!..H1 shape functions
   real*8, dimension(MAXbrickH)    :: shapH
   real*8, dimension(3,MAXbrickH)  :: gradH
!
!..H(curl) shape functions
   real*8, dimension(3,MAXbrickE)  :: shapE
   real*8, dimension(3,MAXbrickE)  :: curlE
!
!..Enriched H1 shape functions
   real*8 , dimension(3,MAXbrickEE) :: shapEE
   real*8 , dimension(3,MAXbrickEE) :: curlEE
!
!..nrdof for interface only (without bubbles)
   !integer :: nrdofEEi
!
!..load vector for the enriched space
   complex(8) :: bload_E(NrTest)
!
!..Gram matrix in packed format
   complex(8), allocatable :: gramP(:)
!
!..matrices for transpose filling (swapped loops)
   complex(8), allocatable :: stiff_EE_T(:,:),stiff_EQ_T(:,:)
   complex(8), allocatable :: stiff_ALL(:,:),zaloc(:,:)
!
!..2D quadrature data
   real*8, dimension(2,MAXNINT2ADD)  :: tloc
   real*8, dimension(  MAXNINT2ADD)  :: wtloc
!
!..BC's flags
   integer, dimension(6,NR_PHYSA)    :: ibc
!
!..for auxiliary computation
   complex(8) :: zaux
!
!..Maxwell load and auxiliary variables
   complex(8), dimension(3) :: zJ,zImp
   real*8    , dimension(3) :: E1,E2,rntimesE,rn2timesE
!
!..number of edge,faces per element type
   integer :: nre, nrf
!
!..various variables for the problem
   real*8  :: h_elem,rjac,weight,wa,v2n,CC,EE,CE,E,EC,q,h,omeg,alpha_scale
   real*8  :: bjac
   integer :: i1,j1,i2,j2
   integer :: i12,j12,k12,k1,k2,fa,fb,i3mod,j3mod,kH,kk,i,ik,j,k,l,nint,kE,n,m,p,pe
   integer :: iflag,iprint,itime,iverb
   integer :: nordP,nsign,ifc,ndom,info,icomp,idec
   complex(8) :: zfval
   complex(8) :: za(3,3),zc(3,3)
   complex(8) :: zaJ(3,3),zcJ(3,3)
!
!..for polarizations function
   complex(8), dimension(3,3) :: bg_pol,gain_pol,raman_pol
   real*8  :: delta_n
   integer :: dom_flag
!
!..OMEGA_RATIO_SIGNAL or OMEGA_RATIO_PUMP
   real*8  :: OMEGA_RATIO_FLD
!
!..for PML
   complex(8) :: zbeta,zdbeta,zd2beta,detJstretch
   complex(8), dimension(3,3) :: Jstretch,invJstretch,JJstretch
!
!..added to use fast integration
!..Gram matrix aux arrays
   complex(8), allocatable :: AUXEE_zb(:,:,:)  , AUXEE_zc(:,:,:)
   complex(8), allocatable :: AUXCC(:,:,:)
   complex(8), allocatable :: AUXEC_zb(:,:,:)  , AUXEC_zc(:,:,:)
   complex(8), allocatable :: AUXCE_zb(:,:,:)  , AUXCE_zc(:,:,:)
!..Stiffness matrix aux arrays
   complex(8), allocatable :: STIFQE(:,:,:)
   complex(8), allocatable :: STIFQE_ALPHA(:,:,:)
   complex(8), allocatable :: STIFQC(:,:,:)
!..Forcing aux array
   complex(8), allocatable :: LOADE(:,:)
!
   integer :: a,b,sa,sb,sc,alph,beta
   integer :: pxy,pz
   integer :: l1,l2,l3,i3,j3,k3,idxbeta,idxalph
   integer :: idxa,idxa2,idxa3,idxb,idxb2,idxb3,m1,m2
   integer :: nord1,nord2,nord3,nintxy,nintz
   integer :: nrdof
   integer :: nrdofH12, nrdofE12, nrdofQ12, nrdofH3, nrdofQ3
   integer :: nrdofH12_tr,nrdofH3_tr
   integer :: nrdofQ12_tr,nrdofQ3_tr
   real*8 :: xi12(2),xi3,wt12,wt3
   real*8 :: wt123,weighthh,weightvv
   real*8 :: xiloc_xy(2,MAXNINT2ADD), xiloc_z(MAXPP+1)
   real*8 :: wloc_xy(MAXNINT2ADD), wloc_z(MAXPP+1)
   real*8, dimension(3,MAXNINT3ADD) :: wloc3
   real*8, dimension(3) :: xip,dHdx,dHHdx
!
   real(8)   , dimension(3,3) :: D_za,D_zc,D_aux,C,D
   complex(8), dimension(3,3) :: Z_za,Z_zc,Z_aux
!
   real(8), allocatable :: shapeH3(:,:),E12(:,:),C12(:,:),Q12(:)
   real(8), allocatable :: sH12(:,:),gH12(:,:,:),sE12(:,:,:),cE12(:,:),sQ12(:,:)
!
   integer, allocatable :: mapEE(:),mapQQ(:)
!
   integer, dimension(3,3) :: deltak
!
!..for Gram matrix compressed storage format
   integer :: nk
   nk(k1,k2) = (k2-1)*k2/2+k1
!
!..Identity/Kronecker delta tensor
   deltak=0
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
      write(*,*) 'elem_maxwell_fi_pris: Mdle = ', Mdle
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
      case('mdlp')
         nord3 = 0
         nord3 = max(nord3,norder(7))
         nord3 = max(nord3,norder(8))
         nord3 = max(nord3,norder(9))
         nord3 = nord3+NORD_ADD
         nordP = NODES(Mdle)%order+NORD_ADD*11
         norderi(nre+nrf+1) = 11
      case('mdlb')
         nord3 = 0
         nord3 = max(nord3,norder(9))
         nord3 = max(nord3,norder(10))
         nord3 = max(nord3,norder(11))
         nord3 = max(nord3,norder(12))
         nord3 = nord3+NORD_ADD
         nordP = NODES(Mdle)%order+NORD_ADD*111
         norderi(nre+nrf+1) = 111
      case default
         write(*,*) 'elem_maxwell_fi_pris: unsupported etype param. stop.'
         stop
   end select
!
   select case(etype)
      case('mdlp') ! prism
         etype1 = 'mdlt'
!     ...calc face order and enriched face order
         norder_f(1:3) = norder(1:3); norder_fe(1:3) = norder(1:3) + NORD_ADD
         norder_f(4) = max(norder(10),norder(11)); norder_fe(4) = norder_f(4) + NORD_ADD
         norder_f(5) = norder_f(4); norder_fe(5) = norder_fe(4)
!
!     ...get element mapping for prism
         p  = NODES(Mdle)%order/10
         pe = nordP/10
         nrdofH12 = (pe+2)*(pe+1)/2    ! test  H1
         nrdofE12 = (pe+2)*(pe)        ! test  H(curl)
         nrdofQ12 = (p+1)*(p)/2        ! trial L2
!
         allocate(mapEE(nrdofE12*(nord3+1) + nrdofH12*nord3)) ! test  dof ordering
         allocate(mapQQ(nrdofQ12*(nord3-NORD_ADD)))           ! trial dof ordering
         call tens_prism_ordEE(pe,nord3         , mapEE)
         call tens_prism_ordQQ(p ,nord3-NORD_ADD, mapQQ)
!
      case('mdlb') ! hexa
         etype1 = 'mdlq'
 !    ...calc face order and enriched face order
         norder_f(1:4) = norder(1:4); norder_fe(1:4) = norder(1:4) + NORD_ADD
!     ...take max order of this face (13) and the opposite face (14)
         norder_f(5) = max(norder(13),norder(14)); norder_fe(5) = norder_f(5) + NORD_ADD*11
!
         p  = NODES(Mdle)%order/110
         pe = nordP/110
         nrdofH12 = (pe+1)*(pe+1)
         nrdofE12 = 2*(pe+1)*(pe)
         nrdofQ12 = (p)*(p)
!
         allocate(mapEE(nrdofE12*(nord3+1) + nrdofH12*nord3)) ! test  dof ordering
         allocate(mapQQ(nrdofQ12*(nord3-NORD_ADD)))           ! trial dof ordering
!        TODO -- implement hexa mapping
         call tens_hexa_ordEE(pe,nord3         , mapEE)
         call tens_hexa_ordQQ(p ,nord3-NORD_ADD, mapQQ)
!
      case default
         write(*,*) 'elem_maxwell_fi_pris: unexpected element type:', etype
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
7001  format('elem_maxwell_fi_pris: BCFLAGS FOR Mdle = ',i5)
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
      write(*,*) 'elem_maxwell_fi_pris: invalid Fld_flag param. stop.'
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
      case(2,3)
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
               write(*,*) 'elem_maxwell_fi_pris: unexpected ndom param. stop.'
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
!..here begins the setup for tensorized num quadrature for prism
!
!  Set triangle integration (TODO etype 1 could be a quad as well, if hexa elem)
   call set_2Dint_DPG(etype1,norder_fe, nintxy,xiloc_xy,wloc_xy)
!
!  Set interval integration
   call set_1Dint_DPG(nord3, nintz,xiloc_z,wloc_z)
! !$OMP CRITICAL
!    write(*,*) 'FAST INT: mdle,nintxy,nintz = ',mdle,nintxy,nintz
! !$OMP END CRITICAL
!
!..Allocate the auxiliary arrays for sum factorization
!
   allocate(shapeH3(2,nord3+1          ))
   allocate(E12    (3,nrdofE12+nrdofH12))
   allocate(C12    (3,nrdofE12+nrdofH12))
   allocate(Q12    (  nrdofQ12         ))
!
   allocate(sH12(  nrdofH12,nintxy))
   allocate(gH12(2,nrdofH12,nintxy))
   allocate(sE12(2,nrdofE12,nintxy))
   allocate(cE12(  nrdofE12,nintxy))
   allocate(sQ12(  nrdofQ12,nintxy))
!
   allocate(AUXEE_zb(3,3,(nrdofE12+nrdofH12)**2))
   allocate(AUXEE_zc(3,3,(nrdofE12+nrdofH12)**2))
   allocate(AUXCC   (3,3,(nrdofE12+nrdofH12)**2))
   allocate(AUXEC_zb(3,3,(nrdofE12+nrdofH12)**2))
   allocate(AUXEC_zc(3,3,(nrdofE12+nrdofH12)**2))
   allocate(AUXCE_zb(3,3,(nrdofE12+nrdofH12)**2))
   allocate(AUXCE_zc(3,3,(nrdofE12+nrdofH12)**2))
!
   allocate(STIFQE      (3,3,nrdofQ12*(nrdofE12+nrdofH12)))
   allocate(STIFQE_ALPHA(3,3,nrdofQ12*(nrdofE12+nrdofH12)))
   allocate(STIFQC      (3,3,nrdofQ12*(nrdofE12+nrdofH12)))
!
   allocate(LOADE(3,(nrdofE12+nrdofH12)))
!
!..Loop over quadrature points in direction \xi_1
   do pz=1,nintz 
!  ...read quadrature point location and weight
      xi3=xiloc_z(pz)
      wt3=wloc_z(pz)
!
!  ...call 1D shape functions for coordinate 3
!  ...Function values stored in shapH3(1,:), derivative values in shapH3(2,:)
      call shape1HH(xi3,nord3, nrdofH3,shapeH3(1,:),shapeH3(2,:))
      nrdofQ3 = nrdofH3 - NORD_ADD - 1
!
!  ...Initialize auxiliary matrices: Gram matrix
      AUXEE_zb = ZERO; AUXEE_zc = ZERO
      AUXCC    = ZERO
      AUXEC_zb = ZERO; AUXEC_zc = ZERO
      AUXCE_zb = ZERO; AUXCE_zc = ZERO
!
!  ...Initialize auxiliary matrices: Stiffness matrix
      STIFQE = ZERO; STIFQE_ALPHA = ZERO
      STIFQC = ZERO;
!
!  ...Initialize auxiliary matrices: Forcing Term
      LOADE = ZERO
!
!  ...loop over xy quadrature points
      do pxy=1,nintxy
!     ...read quadrature point location and weight
         xi12=xiloc_xy(1:2,pxy)
         wt12=wloc_xy(pxy)
!
!     ...Shape function subroutine is called only once, when
!        pz=1 and stored in sH2p(:,py) and dsH2p(:,py)
         if (pz.eq.1) then
            sH12(:,pxy)   = rZero;  gH12(:,:,pxy) = rZero
            sE12(:,:,pxy) = rZero;  cE12(:,pxy)   = rZero
            sQ12(:,pxy)   = rZero;
!
            call shape2HH(etype1,xi12,norder_fe(5), nrdofH12,sH12(:,pxy),gH12(1:2,:,pxy))
            call shape2EE(etype1,xi12,norder_fe(5), nrdofE12,sE12(1:2,:,pxy),cE12(:,pxy))
!
            call shape2QQ(etype1,xi12,norder_f(5), nrdofQ12,sQ12(:,pxy))
         endif
!
!     ...Copy shape functions in coord. 2 previously evaluated
!     ...E12 is for E terms, separated into family 1 and 2 of tri shape fns.
         E12(1:2, 1:NrdofE12 ) = sE12(1:2,1:NrdofE12,pxy)
         E12(3  , 1:NrdofE12 ) = rZero
         E12(1:2, (NrdofE12+1):(NrdofE12+NrdofH12) ) = rZero
         E12(3  , (NrdofE12+1):(NrdofE12+NrdofH12) ) = sH12(1:NrdofH12,pxy)
!
!     ...C12 is for C terms, curl of E terms, same separation as above
         C12(1 , 1:NrdofE12 ) = -sE12(2,1:NrdofE12,pxy)
         C12(2 , 1:NrdofE12 ) =  sE12(1,1:NrdofE12,pxy)
         C12(3 , 1:NrdofE12 ) =  cE12(  1:NrdofE12,pxy)
         C12(1 , (NrdofE12+1):(NrdofE12+NrdofH12) ) =  gH12(2,1:NrdofH12,pxy)
         C12(2 , (NrdofE12+1):(NrdofE12+NrdofH12) ) = -gH12(1,1:NrdofH12,pxy)
         C12(3 , (NrdofE12+1):(NrdofE12+NrdofH12) ) = rZero
!
!     ...Q12 is for L2 trial shape functions
         Q12(1:nrdofQ12) = sQ12(1:nrdofQ12,pxy)
!
!     ...Condense quad pt into single array
         xip(1:2) = xi12(1:2); xip(3) = xi3
!
!     ...Compute shape functions needed for geometry - 3D H1 shape functions
         call shape3H(etype,xip,norder,norient_edge,norient_face, nrdof,shapH,gradH)
!     ...Geometry map
         call geom3D(Mdle,xip,xnod,shapH,gradH,NrdofH, x,dxdxi,dxidx,rjac,iflag)
#if DEBUG_MODE
         if (iflag.ne.0) then
            write(*,5999) Mdle,rjac
5999        format('elem_maxwell_fi_pris: Negative Jacobian. Mdle,rjac=',i8,2x,e12.5)
            stop
         endif
#endif
!
!     ...get the RHS
         call getf(Mdle,x, zfval,zJ)
!
!     ...set auxiliary constants (updated below if nonlinear)
!        ...update bgpol depending on coordinates if refractive index is varying
!           (other than step-index)
         za = (ZI*OMEGA*OMEGA_RATIO_FLD*EPSILON+SIGMA)*IDENTITY+bg_pol
         zc = (ZI*OMEGA*OMEGA_RATIO_FLD*MU)*IDENTITY
         if (GEOM_NO .eq. 5) then
            call get_bgPol(dom_flag,Fld_flag,0.d0,x, bg_pol)
         endif
!
!     ...Nonlinear terms
         if(NONLINEAR_FLAG.eq.1) then
!        ...compute current solution using soleval
            nflag = 1
            call soleval(Mdle,xip,norient_edge,norient_face,norder,xnod,&
              zdofH,zdofE,zdofV,zdofQ,nflag,x,dxdxi, &
              zsolH_soleval,zdsolH_soleval,zsolE_soleval,zcurlE_soleval,&
              zsolV_soleval,zdivV_soleval,zsolQ_soleval)
!        ...initialize material refractive index perturbation
            delta_n = 0.d0
            if(HEAT_FLAG.eq.1) then
!           ...compute thermally induced perturbation to refractive index at x
               rsolH = real(zsolH_soleval(1))
               delta_n = THERMO_OPT_COEFF*rsolH
            endif
!
!        ...update background polarization with thermal perturbation
            if(delta_n .ne. 0.d0) then
               call get_bgPol(dom_flag,Fld_flag,delta_n,x, bg_pol)
            endif
!
!        ...initialize gain polarization, raman polarization
            gain_pol = ZERO; raman_pol = ZERO
            if (ACTIVE_GAIN .gt. 0.d0) then
               if (dom_flag .eq. 1) then ! .and. x(3).le.PML_REGION) then
                  call get_activePol(zsolQ_soleval(1:12),Fld_flag,delta_n, gain_pol)
               endif
            endif
            if (RAMAN_GAIN .gt. 0.d0) then
!           ...Fld_flag = 1 when in signal element routine
               if(Fld_flag .eq. 1) then
                  call get_ramanPol(zsolQ_soleval(7:9),zsolQ_soleval(10:12), &
                                 dom_flag,Fld_flag,delta_n, raman_pol)
               else
                  call get_ramanPol(zsolQ_soleval(1:3),zsolQ_soleval(4:6), &
                                 dom_flag,Fld_flag,delta_n, raman_pol)
               endif
!        ...endif RAMAN_GAIN
            endif
!        ...update auxiliary constants for za,zb: this is for
!        ...Stiffness and Gram matrix that changes with each nonlinear iteration
            za = (ZI*OMEGA*OMEGA_RATIO_FLD*EPSILON+SIGMA)*IDENTITY+bg_pol+gain_pol+raman_pol
!     ...endif NONLINEAR_FLAG
         endif
!
!.....................................................
!...............toggle PML............................
!
         if(USE_PML.eq.0) then
            JJstretch      = ZERO
            JJstretch(1,1) = ZONE
            JJstretch(2,2) = ZONE
            JJstretch(3,3) = ZONE
         else
!        ...get PML function
            call get_Beta(x,Fld_flag, zbeta,zdbeta,zd2beta)
            Jstretch(3,3) = zdbeta
!        ...compute det(J) * J^-1 * J^-T
!           (J is a diagonal matrix)
            invJstretch(3,3) = 1.d0/zdbeta
            call ZGEMM('N', 'N', 3, 3, 3, ZONE, invJstretch, 3, &
                          invJstretch, 3, ZERO, JJstretch, 3)
            detJstretch = zdbeta
            JJstretch = detJstretch*JJstretch
         endif
!
!     ...PML stretching
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
!     ...compute total quadrature weight
         wt123=wt12*wt3
!     ...compute Jacobian determinant * quadrature weight
         weighthh=wt123*rjac
!     ...Determine D = J^-1 * J^-T. Multiply by appropriate weight
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
!     ...compute inverse Jacobian determinant * quadrature weight
         weightvv=wt123/rjac
!     ...Determine C = J^T * J.  Multiply by appropriate weight
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
!
!     ...Determine auxiliary matrices to anisotropic refractive index tensor
!        (also needed for PML stretching in adjoint graph norm)
!     ...D_za = J^-1 |zaJ|^2 J^-T
         D_aux(1:3,1) = dxidx(1:3,1) * (abs(zaJ(1,1))**2)
         D_aux(1:3,2) = dxidx(1:3,2) * (abs(zaJ(2,2))**2)
         D_aux(1:3,3) = dxidx(1:3,3) * (abs(zaJ(3,3))**2)
         call DGEMM('N','T',3,3,3,1.0d0,D_aux,3,dxidx,3,0.0d0,D_za,3)
         D_za = D_za*weighthh
!     ...D_zc = J^-1 |zcJ|^2 J^-T
         D_aux(1:3,1) = dxidx(1:3,1) * (abs(zcJ(1,1))**2)
         D_aux(1:3,2) = dxidx(1:3,2) * (abs(zcJ(2,2))**2)
         D_aux(1:3,3) = dxidx(1:3,3) * (abs(zcJ(3,3))**2)
         call DGEMM('N','T',3,3,3,1.0d0,D_aux,3,dxidx,3,0.0d0,D_zc,3)
         D_zc = D_zc*weighthh
!     ...Z_za = J^-1 * zaJ^T * J
         call ZGEMM('T','N',3,3,3,ZONE,zaJ        ,3,ZONE*dxdxi,3,ZERO,Z_aux,3)
         call ZGEMM('N','N',3,3,3,ZONE,ZONE*dxidx ,3,     Z_aux,3,ZERO,Z_za ,3)
!     ...Z_zc = J^-1 * zcJ^T * J
         call ZGEMM('T','N',3,3,3,ZONE,zcJ        ,3,ZONE*dxdxi,3,ZERO,Z_aux,3)
         call ZGEMM('N','N',3,3,3,ZONE,ZONE*dxidx ,3,     Z_aux,3,ZERO,Z_zc ,3)
!     ...Z_zb = conjg(transpose(Z_za)) = J^T * zaJ * J^-T
!     ...Z_zd = conjg(transpose(Z_zc)) = J^T * zcJ * J^-T
!
!     ...put appropriate quadrature weight on Jacobian and its inverse
         dxdxi = dxdxi * weightvv
         dxidx = dxidx * wt123
!
!--------------------------------------------------------------------------
!        ...SUM FACTORIZATION LOOPS START
!--------------------------------------------------------------------------
!        ... For Gram, Stiffness, Load:
!        ... Hcurl TEST  functions will be identified by indices i12,i3,a
!
!        ... For Gram matrix:
!        ... Hcurl TRIAL functions will be identified by indices j12,j3,b
!
!        ... For stiffness matrix:
!        ... L2    TRIAL functions will be identified by indices j12,j3,b
!--------------------------------------------------------------------------
!
!     ...Integration of innermost integrals for Gram matrix
!     ...loop over 2D tri test function, Hcurl
         do i12=1,(nrdofE12+nrdofH12)
!        ...loop over 2D DOFs, coord. 12, test shape func Hcurl
            do j12=1,(nrdofE12+nrdofH12)
!           ...combine indices i12 and j12 into k12
               k12 = (i12-1)*(nrdofE12+nrdofH12) + j12
!
!           ...loop over vector components
               do b=1,3
                  do a=1,3
!
! Note E terms have a lot of 0s so we only compute if non-zero
!                 ...Accumulate for EE terms
!                    ...For EE11
                        AUXEE_zb(a,b,k12) = AUXEE_zb(a,b,k12)     &
                              + (ALPHA_NORM*D(a,b) + D_za(a,b))   &
                              * E12(a,i12)*E12(b,j12)
!                    ...For EE22
                        AUXEE_zc(a,b,k12) = AUXEE_zc(a,b,k12)     &
                              + (ALPHA_NORM*D(a,b) + D_zc(a,b))   &
                              * E12(a,i12)*E12(b,j12)
!
!                 ...Accumulate for EC terms
                     AUXEC_zb(a,b,k12) = AUXEC_zb(a,b,k12)        &
                           - Z_za(a,b)                            &
                           * E12(a,i12)*C12(b,j12)*wt123
!
                     AUXEC_zc(a,b,k12) = AUXEC_zc(a,b,k12)        &
                           + Z_zc(a,b)                            &
                           * E12(a,i12)*C12(b,j12)*wt123
!
!                 ...Accumulate for CC terms
                     AUXCC(a,b,k12) = AUXCC(a,b,k12)                 &
                           + C12(a,i12)*C12(b,j12)*C(a,b)
!
!                 ...Accumulate for CE terms
!                    (Z_za and Z_zc indices are switched here b/c we need the transpose)
                     AUXCE_zb(a,b,k12) = AUXCE_zb(a,b,k12)        &
                           - conjg(Z_za(b,a))                     &
                           * C12(a,i12)*E12(b,j12)*wt123
!
                     AUXCE_zc(a,b,k12) = AUXCE_zc(a,b,k12)        &
                           + conjg(Z_zc(b,a))                     &
                           * C12(a,i12)*E12(b,j12)*wt123
!
!              ...loop a ends
                  enddo
!           ...loop b ends
               enddo
!        ...loop over 2D test function ends
            enddo
!
! ------------ Assemble Aux Stiffness Matrix --------------
!        ...Loop over 2D trial function
            do j12=1,nrdofQ12
               k12 = (i12-1)*nrdofQ12 + j12
               do b=1,3
                  do a=1,3
                     STIFQE(a,b,k12) = STIFQE(a,b,k12)                  &
                           + (dxidx(a,b)*Q12(j12))                      &
                           * zcJ(b,b)                                   &
                           * E12(a,i12)
                     STIFQE_ALPHA(a,b,k12) = STIFQE_ALPHA(a,b,k12)      &
                           - (dxidx(a,b)*Q12(j12))                      &
                           * zaJ(b,b)                                   &
                           * E12(a,i12)
!                 ...Curl
                     STIFQC(a,b,k12) = STIFQC(a,b,k12)                     &
                           + Q12(j12)*dxdxi(b,a)*C12(a,i12)
                  enddo
               enddo
!
!        ...loop over 2D trial function ends
            enddo
!
! ------------ Assemble Aux Load Vector --------------    
            do a=1,3
               LOADE(a,i12) = LOADE(a,i12)                              &
                     + ( dxidx(a,1)*zJ(1)+                              &
                         dxidx(a,2)*zJ(2)+                              &
                         dxidx(a,3)*zJ(3) )                             &
                     * E12(a,i12)*rjac
            enddo
!
!     ...loop over 2D test function ends
         enddo
!  ...loop over pxy ends
      enddo
!
! --------------------------------------------
! FINAL COMPUTATION OF MATRICES
! --------------------------------------------
!
!  ...Assemble Gram Matrix
!  ...Loop 1 over 1D test functions
      do i3=1,nrdofH3
!     ...Loop 1 over 2D test functions
         do i12=1,(nrdofE12+nrdofH12)

!        ...Switch between shape functions families
            if (i12 .le. nrdofE12) then
!           ...Parameters i3,sa,fa, switch between fn value or deriv
               i3mod = i3; sa = 1; fa = 1
               m1 = mapEE(i3 + (i12-1)*nrdofH3)
            else
               if (i3 .eq. nrdofH3) exit
               i3mod = i3 + 1; sa = 2; fa = 0
               m1 = mapEE(i3 + (i12-nrdofE12-1)*(nrdofH3-1) + nrdofE12*nrdofH3)
            endif
!
            if (i3mod.le.nrdofH3) then
!
!           ...Loop 2 over 1D test functions
               do j3=1,nrdofH3
!              ...Loop 2 over 2D test functions
                  do j12=1,(nrdofE12+nrdofH12)
                     k12 = (i12-1)*(nrdofE12+nrdofH12) + j12
!
!                 ...Switch between shape functions families
                     if (j12.le.nrdofE12) then
!                    ...Parameters j3,sb,fb, switch between fn value or deriv
                        j3mod = j3; sb = 1; fb = 1
                        m2 = mapEE(j3 + (j12-1)*nrdofH3)
                     else
                        if (j3 .eq. nrdofH3) exit
                        j3mod = j3 + 1; sb = 2; fb = 0
                        m2 = mapEE(j3 + (j12-nrdofE12-1)*(nrdofH3-1) + nrdofE12*nrdofH3)
                     endif
!
                     if (m1.le.m2) then 
                        if (j3mod.le.nrdofH3) then
                           do b=1,3
                              do a=1,3
!                             ...Sum EE_11 and CC_11 terms
                                 kk = nk(2*m1-1,2*m2-1)
!
                                 gramP(kk) = gramP(kk)                    &
                                    + AUXEE_zb(a,b,k12)                   &
                                       * shapeH3(sa,i3mod)                & 
                                       * shapeH3(sb,j3mod)                &
                                    + AUXCC(a,b,k12)                      &
                                       * shapeH3(2-deltak(a,3)*fa,i3mod)  &
                                       * shapeH3(2-deltak(b,3)*fb,j3mod) 
!
!                             ...Sum CE and EC terms
                                 kk = nk(2*m1-1,2*m2)
                                 gramP(kk) = gramP(kk)                    &
                                    + AUXCE_zc(a,b,k12)                   &
                                       * shapeH3(2-deltak(a,3)*fa,i3mod)  &
                                       * shapeH3(sb,j3mod)                &
                                    + AUXEC_zb(a,b,k12)                   &
                                       * shapeH3(sa,i3mod)                &
                                       * shapeH3(2-deltak(b,3)*fb,j3mod)  
!
!                             ...Sum other CE and EC terms
                                 if (m1.ne.m2) then
                                    kk = nk(2*m1  ,2*m2-1)
                                    gramP(kk) = gramP(kk)                    &
                                       + AUXCE_zb(a,b,k12)                   &
                                          * shapeH3(2-deltak(a,3)*fa,i3mod)  &
                                          * shapeH3(sb,j3mod)                &
                                       + AUXEC_zc(a,b,k12)                   &
                                          * shapeH3(sa,i3mod)                &
                                          * shapeH3(2-deltak(b,3)*fb,j3mod)   
                                 endif
!
!                             ...Sum EE_22 and CC_22 terms
                                 kk = nk(2*m1  ,2*m2  )
                                 gramP(kk) = gramP(kk)                    &
                                    + AUXEE_zc(a,b,k12)                   &
                                       * shapeH3(sa,i3mod)                & 
                                       * shapeH3(sb,j3mod)                &
                                    + AUXCC(a,b,k12)                      &
                                       * shapeH3(2-deltak(a,3)*fa,i3mod)  &
                                       * shapeH3(2-deltak(b,3)*fb,j3mod) 
!                          ...end loop a
                              enddo
!                       ...end loop b
                           enddo
                        endif
                     endif
!              ...end loop 2 over 2D test functions
                  enddo
!           ...end loop 2 over 1D test functions
               enddo
!
!           ...Assemble Stiffness Matrix
!           ...Loop over 1D trial functions
               do j3=1,nrdofQ3
!              ...Loop over 2D trial functions
                  do j12=1,nrdofQ12
                     k12 = (i12-1)*nrdofQ12 + j12
                     m2 = mapQQ(j3 + (j12-1)*nrdofQ3)
                     j3mod = j3 + 1
                     do b=1,3
                        do a=1,3
!                       ...accumulate QC terms
                           k = (m2-1)*6 + 3 + b
                           stiff_EQ_T(k,2*m1-1) = stiff_EQ_T(k,2*m1-1)     &
                                 + STIFQC(a,b,k12)*shapeH3(2,j3mod)           &
                                 * shapeH3(2-deltak(a,3)*fa,i3mod)
!
!                       ...accumulate QE terms
                           k = (m2-1)*6 + b
                           stiff_EQ_T(k,2*m1-1) = stiff_EQ_T(k,2*m1-1)     &
                                 + STIFQE_ALPHA(a,b,k12)                   &
                                 * shapeH3(sa,i3mod)*shapeH3(2,j3mod)
!
!                       ...accumulate QE terms
                           k = (m2-1)*6 + 3 + b
                           stiff_EQ_T(k,2*m1) = stiff_EQ_T(k,2*m1)         &
                                 + STIFQE(a,b,k12)                         &
                                 * shapeH3(sa,i3mod)*shapeH3(2,j3mod)
!
!                       ...accumulate QC terms
                           k = (m2-1)*6 + b
                           stiff_EQ_T(k,2*m1) = stiff_EQ_T(k,2*m1)         &
                                 + STIFQC(a,b,k12)*shapeH3(2,j3mod)        &
                                 * shapeH3(2-deltak(a,3)*fa,i3mod)
!                    ...end loop a
                        enddo
!                 ...end loop b
                     enddo
!              ...end loop over 2D trial functions
                  enddo
!           ...end loop over 1D trial functions
               enddo
!
!           ...Assemble Forcing Term
               do a=1,3
                  bload_E(2*m1-1) = bload_E(2*m1-1)                     &
                                  + LOADE(a,i12)*shapeH3(sa,i3mod)
               enddo
            endif
!     ...end loop 1 over 2D test functions
         enddo
!  ...end loop 1 over 1D test functions
      enddo
!..end loop over z direction quadrature points
   enddo
!
   deallocate(mapEE,mapQQ)
   deallocate(shapeH3,E12,C12,Q12)
   deallocate(sH12,gH12,sE12,cE12,sQ12)
!
   deallocate(AUXEE_zb,AUXEE_zc)
   deallocate(AUXCC)
   deallocate(AUXEC_zb,AUXEC_zc)
   deallocate(AUXCE_zb,AUXCE_zc)
!
   deallocate(STIFQE)
   deallocate(STIFQE_ALPHA)
   deallocate(STIFQC)
!
   deallocate(LOADE)
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
! !$OMP CRITICAL
!       write(*,*) 'FAST INT: mdle,ftype,nint = ',mdle,ftype,nint
! !$OMP END CRITICAL
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
            write(*,*) 'elem_maxwell: INCONSISTENCY NrdofEE. stop.'
            stop
         endif
#endif
!
!     ...determine element H1 shape functions (for geometry)
         call shape3H(etype,xi,norder,norient_edge,norient_face, &
                      nrdof,shapH,gradH)
#if DEBUG_MODE
         if (nrdof .ne. NrdofH) then
            write(*,*) 'elem_maxwell: INCONSISTENCY NrdofH. stop.'
            stop
         endif
#endif
!
!     ...determine element H(curl) shape functions (for fluxes)
!     ...for interfaces only (no bubbles)
         call shape3E(etype,xi,norderi,norient_edge,norient_face, &
                      nrdof,shapE,curlE)
#if DEBUG_MODE
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
      write(*,*) 'elem_maxwell_fi_pris: ZPPTRF: Mdle,info = ',Mdle,info,'. stop.'
      stop
   endif
!
!..B. Solve triangular system to obtain B~, (LX=) U^*X = [B|l]
   call ZTPTRS('U','C','N',NrTest,NrTrial+1,gramP,stiff_ALL,NrTest,info)
   if (info.ne.0) then
      write(*,*) 'elem_maxwell_fi_pris: ZTPTRS: Mdle,info = ',Mdle,info,'. stop.'
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
end subroutine elem_maxwell_fi_pris
!
!----------------------------------------------------------------------
! routine: tens_prism_ordQQ
!----------------------------------------------------------------------
subroutine tens_prism_ordQQ(p,q, mapQQ)
   use parametersDPG
   integer, intent(in)  :: p
   integer, intent(in)  :: q
   integer, intent(out) :: mapQQ((p+1)*p*q/2)
   map = 0
   k = 0
   ! Bubble functions (Tri Bub. x Unit edges) [1/2*(p+1)*p*q]
   do j = 1,q
      do i=1,(p+1)*(p)/2
         k = k+1
         mapQQ((i-1)*(q) + j) = k
      enddo
   enddo
end subroutine tens_prism_ordQQ
!
!----------------------------------------------------------------------
! routine: tens_prism_ordEE
!----------------------------------------------------------------------
subroutine tens_prism_ordEE(p,q, mapEE)
   use parametersDPG
   integer, intent(in)  :: p
   integer, intent(in)  :: q
   integer, intent(out) :: mapEE((p+2)*p*(q+1)+(p+2)*(p+1)*q/2)
   integer :: i,j,k
!
   mapEE = 0
   k = 0
!
!..Mixed edge functions (Tri Edges x Unit Verts)    [6*p]
!..Have edges from both T(curl) [p] and T(H1) [p-2] to account for
   do i=1,3*p
      k = k+1
      mapEE(i*(q+1) - (q)) = k
   enddo
   do i=1,3*p
      k = k+1
      mapEE(i*(q+1) - (q-1)) = k
   enddo
!
!..Quad edge functions (Tri Verts. x Unit Edges)     [3*q]
   do i=1,3*q
      k = k+1
      mapEE(p*(p+2)*(q+1) + i) = k
   enddo
!
!..Tri face functions (Tri Bub. x Unit Edges)        [2*p*(p-1)]
   do i=1,p*(p-1)
      k = k+1
      mapEE(3*p*(q+1) + i*(q+1) - (q)) = k
   enddo
   do i=1,p*(p-1)
      k = k+1
      mapEE(3*p*(q+1) + i*(q+1) - (q-1)) = k
   enddo
!
!.. Quad face functions (Tri Hcurl Edges x Unit H1 Edges
!                       and Tri H1 Edges x Unit L2 Edges)    [3*(p)*(q-1) + 3*q*(p-1)]
   do j=1,(q-1)
      do i=1,(p)
         k = k+1
         mapEE(i*(q+1) - (q-1) + j) = k
      enddo
   enddo
   do j=1,(q)
      do i=1,(p-1)
         k = k+1
         mapEE(p*(p+2)*(q+1) + 3*q + (i-1)*q + j) = k
      enddo
   enddo
   do j=1,(q-1)
      do i=1,(p)
         k = k+1
         mapEE((i+(p))*(q+1) - (q-1) + j) = k
      enddo
   enddo
   do j=1,(q)
      do i=1,(p-1)
         k = k+1
         mapEE(p*(p+2)*(q+1) + 3*q + (i+(p-1)-1)*q + j) = k
      enddo
   enddo
   do j=1,(q-1)
      do i=1,(p)
         k = k+1
         mapEE((i+2*(p))*(q+1) - (q-1) + j) = k
      enddo
   enddo
   do j=1,(q)
      do i=1,(p-1)
         k = k+1
         mapEE(p*(p+2)*(q+1) + 3*q + (i+2*(p-1)-1)*q + j) = k
      enddo
   enddo
!
!..Bubble functions (Tri Bub. x Unit edges)  [p*(p-1)*(q-1) + 1/2*(p-1)*(p-2)*q]
!..Tri Hcurl Bub. x H1 Edges
   do j=1,(q-1)
      do i=1,p*(p-1)
         k = k+1
         mapEE(3*p*(q+1) + (i-1)*(q+1) + 2 + j) = k
      enddo
   enddo
!
!..Tri H1 Bub. x L2 edges
   do j = 1,q
      do i=1,(p-1)*(p-2)/2
         k = k+1
         mapEE(p*(p+2)*(q+1) + 3*q*p + (i-1)*q + j) = k
      enddo
   enddo
!
end subroutine tens_prism_ordEE
!
!----------------------------------------------------------------------
! routine: tens_hexa_ordEE
!----------------------------------------------------------------------
! TODO
subroutine tens_hexa_ordEE(p,q, mapEE)
   use parametersDPG
   integer, intent(in)  :: p 
   integer, intent(in)  :: q
   integer, intent(out) :: mapEE(2*(p+1)*p*(q+1) + (p+1)*(p+1)*q)
   integer :: i,j,k
   map = 0
   k = 0
   do i = 1,2*(p+1)*p*(q+1) + (p+1)*(p+1)*q
      k = k+1
      mapEE(i) = k
   enddo
end subroutine tens_hexa_ordEE
!
!----------------------------------------------------------------------
! routine: tens_hexa_ordQQ
!----------------------------------------------------------------------
! TODO
subroutine tens_hexa_ordQQ(p,q, mapQQ)
   use parametersDPG
   integer, intent(in)  :: p 
   integer, intent(in)  :: q
   integer, intent(out) :: mapQQ(p*p*q)
   integer :: i,j,k
   map = 0
   k = 0
   do i = 1,p*p*q
      k = k+1
      mapQQ(i) = k
   enddo
end subroutine tens_hexa_ordQQ
