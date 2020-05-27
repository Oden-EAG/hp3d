!--------------------------------------------------------------------
!
!     routine name      - elem_maxwell_gram_hexa
!
!--------------------------------------------------------------------
!
!     latest revision:  - Oct 2019
!
!     purpose:          - routine returns element Gram matrix
!                         for the ultraweak Maxwell formulation
!                         using fast integration for hexahedra
!
!     arguments:
!
!     in:
!              Mdle     - an element middle node number, identified
!                         with the element
!              Fld_flag - field flag (0: pump, 1: signal)
!              NrTest   - total number of test dof
!              NrdofH   - number of H1 trial dof
!     out:
!              GramP    - gram matrix
!
!---------------------------------------------------------------------
!
#include "implicit_none.h"
!
   subroutine elem_maxwell_gram_hexa(Mdle,Fld_flag,NrTest,NrdofH, GramP)
!
   use data_structure3D
   use control, only: INTEGRATION
   use commonParam
   use laserParam
   use parametersDPG
!
   implicit none
!
   integer, intent(in)  :: Mdle
   integer, intent(in)  :: Fld_flag
   integer, intent(in)  :: NrTest
   integer, intent(in)  :: NrdofH
!
   VTYPE  , intent(out) :: GramP(NrTest*(NrTest+1)/2)
!
!..locals
   character(len=4) :: etype
!
!..number of edge,faces per element type
   integer :: nre, nrf
!
!..declare element order, orientation for edges and faces
   integer, dimension(19) :: norder
   integer, dimension(12) :: norient_edge
   integer, dimension(6)  :: norient_face
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
   real(8), dimension(MAXbrickH)   :: shapH
   real(8), dimension(3,MAXbrickH) :: gradH
!
!..3D quadrature data
   real(8), dimension(3,MAXNINT3ADD) :: xiloc
   real(8), dimension(MAXNINT3ADD)   :: waloc
!
!..for auxiliary computation
   VTYPE :: zaux
!
!..various variables for the problem
   real(8) :: h_elem,rjac,weight,wa,CC,EE,CE,E,EC,q,h,omeg
   real(8) :: bjac
   integer :: i1,i2,j1,j2,k1,k2,kH,kk,i,ik,j,k,l,nint,kE,n,m
   integer :: iflag,iprint,itime,iverb
   integer :: nrdof,nordP,nsign,ifc,ndom,info,icomp,idec
   complex(8) :: zfval
   complex(8) :: za(3,3),zc(3,3)
   complex(8) :: zaJ(3,3),zcJ(3,3)
!
!..for PML
   VTYPE :: zbeta,zdbeta,zd2beta,detJstretch
   VTYPE, dimension(3,3) :: Jstretch,invJstretch,JJstretch
!
!..for polarizations function
   VTYPE, dimension(3,3) :: bg_pol,gain_pol,raman_pol,rndotE
   real(8) :: delta_n
   integer :: dom_flag
!
!..OMEGA_RATIO_SIGNAL or OMEGA_RATIO_PUMP
   real(8) :: OMEGA_RATIO_FLD
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
!----------------------------------------------------------------------
!
!..element type
   etype = NODES(Mdle)%type
   nre = nedge(etype); nrf = nface(etype)
!
!..determine order of approximation
   call find_order(Mdle, norder)
!
!..set the enriched order of approximation
   select case(etype)
      case('mdlb')
         nordP = NODES(Mdle)%order+NORD_ADD*111
      case default
         write(*,*) 'elem_maxwell_gram_hexa: unsupported etype param. stop.'
         stop
   end select
!
!..determine edge and face orientations
   call find_orient(Mdle, norient_edge,norient_face)
!
!..determine nodes coordinates
   call nodcor(Mdle, xnod)
!
   GramP = ZERO
!
!..set OMEGA_RATIO_FLD
   select case(Fld_flag)
      case(0)
         OMEGA_RATIO_FLD = OMEGA_RATIO_PUMP
      case(1)
         OMEGA_RATIO_FLD = OMEGA_RATIO_SIGNAL ! 1.0d0
      case default
      write(*,*) 'elem_maxwell_gram_hexa: invalid Fld_flag param. stop.'
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
         write(*,*) 'elem_maxwell_gram_hexa: cannot have prism core geometry with fast integration. stop.'
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
               write(*,*) 'elem_maxwell_gram_hexa: unexpected ndom param. stop.'
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
!
!        ...Geometry map
            call geom3D(Mdle,xip,xnod,shapH,gradH,NrdofH, x,dxdxi,dxidx,rjac,iflag)
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
            D_aux(1:3,1) = dxidx(1:3,1) * (abs(zaJ(1,1))**2)
            D_aux(1:3,2) = dxidx(1:3,2) * (abs(zaJ(2,2))**2)
            D_aux(1:3,3) = dxidx(1:3,3) * (abs(zaJ(3,3))**2)
            call DGEMM('N','T',3,3,3,1.0d0,D_aux,3,dxidx,3,0.0d0,D_za,3)
            D_za = D_za*weighthh
!        ...D_zc = J^-1 |zcJ|^2 J^-T
            D_aux(1:3,1) = dxidx(1:3,1) * (abs(zcJ(1,1))**2)
            D_aux(1:3,2) = dxidx(1:3,2) * (abs(zcJ(2,2))**2)
            D_aux(1:3,3) = dxidx(1:3,3) * (abs(zcJ(3,3))**2)
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

!        ...put appropriate quadrature weight on Jacobian and its inverse
            dxdxi = dxdxi * weightvv
            dxidx = dxidx * wt123
!
!--------------------------------------------------------------------------
!        ...SUM FACTORIZATION LOOPS START
!--------------------------------------------------------------------------
!        ... Hcurl TEST  functions will be identified by indices i1,i2,i3,a
!        ... Hcurl TRIAL functions will be identified by indices j1,j2,j3,b
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
!        ...loop over pz ends
         enddo
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
!  ...loop over py ends
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
                                       GramP(kk) = GramP(kk)         &
                                                 + shapH1(idxa,sa)   &
                                                 * shapH1(idxb,sb)   &
                                                 * AUXEE_B_zb(b,a,k2,k3)
!                                   ...sum CC terms
                                       do beta=1,2; do alph=1,2
                                          idxbeta=mod(b+beta-1,3)+1
                                          idxalph=mod(a+alph-1,3)+1
                                          sb=1+1-deltak(idxbeta,1)
                                          sa=1+1-deltak(idxalph,1)
                                          GramP(kk) = GramP(kk)          &
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
                                          GramP(kk) = GramP(kk)                   &
                                                    + AUXCE_B_zc(alph,b,a,k2,k3)  &
                                                    * shapH1(idxa,sa)*shapH1(idxb,sb)
                                       enddo
!                                   ...sum EC terms
                                       do beta=1,2
                                          idxbeta=mod(b+beta-1,3)+1
                                          sb=1+1-deltak(idxbeta,1)
                                          sa=1+deltak(a,1)
                                          GramP(kk) = GramP(kk)                     &
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
                                             GramP(kk) = GramP(kk)                     &
                                                       + AUXCE_B_zb(alph,b,a,k2,k3)    &
                                                       * shapH1(idxa,sa)*shapH1(idxb,sb)
                                          enddo
!                                      ...sum EC terms
                                          do beta=1,2
                                             idxbeta=mod(b+beta-1,3)+1
                                             sb=1+1-deltak(idxbeta,1)
                                             sa=1+deltak(a,1)
                                             GramP(kk) = GramP(kk)                   &
                                                       + AUXEC_B_zc(beta,b,a,k2,k3)  &
                                                       * shapH1(idxa,sa)*shapH1(idxb,sb)
                                          enddo

                                       endif

                                       kk = nk(2*m1  ,2*m2  )
!                                   ...sum EE terms
                                       sb=1+deltak(b,1)
                                       sa=1+deltak(a,1)
                                       GramP(kk) = GramP(kk)         &
                                                 + shapH1(idxa,sa)   &
                                                 * shapH1(idxb,sb)   &
                                                 * AUXEE_B_zc(b,a,k2,k3)
!                                   ...sum CC terms
                                       do beta=1,2; do alph=1,2
                                          idxbeta=mod(b+beta-1,3)+1
                                          idxalph=mod(a+alph-1,3)+1
                                          sb=1+1-deltak(idxbeta,1)
                                          sa=1+1-deltak(idxalph,1)
                                          GramP(kk) = GramP(kk)         &
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
!..loop over px ends
   enddo
!
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
!
end subroutine elem_maxwell_gram_hexa

