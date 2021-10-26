!-------------------------------------------------------------------------
!
!     routine name         - elem_residual_maxwell
!
!-------------------------------------------------------------------------
!
!     latest revision:     - June 2021
!
!     purpose:             - routine returns element residual (squared)
!                            for UW time-harmonic Maxwell equation
!
!     arguments:
!        in:
!              Mdle        - element middle node number
!              Fld_flag    - field flag (0: pump, 1: signal)
!              NrdofEE     - number of H(curl) test dof
!              NrdofH      - number of H1 trial dof
!              NrdofE      - number of H(curl) trial dof
!              NrdofQ      - number of L2 trial dof
!        out:
!              Resid       - element residual (squared)
!              Nref_flag   - suggested h-refinement flag
!
!--------------------------------------------------------------------------
!
#include "typedefs.h"
!
subroutine elem_residual_maxwell(Mdle,Fld_flag,          &
                                 NrTest,NrdofEE,         &
                                 NrdofH,NrdofE,NrdofQ,   &
                                 Resid,Nref_flag)
!..modules used
   use control
   use parametersDPG
   use element_data
   use data_structure3D
   use commonParam
   use laserParam
!..no implicit statements
   implicit none
!..declare input/output variables
   integer, intent(in)  :: Mdle
   integer, intent(in)  :: Fld_flag
   integer, intent(in)  :: NrTest
   integer, intent(in)  :: NrdofEE
   integer, intent(in)  :: NrdofH
   integer, intent(in)  :: NrdofE
   integer, intent(in)  :: NrdofQ
   integer, intent(out) :: Nref_flag
   real(8), intent(out) :: Resid
!
!..declare edge/face type variables
   character(len=4) :: etype,ftype
!
!..declare element order, orientation for edges and faces
   integer, dimension(19) :: norder
   integer, dimension(12) :: norient_edge
   integer, dimension(6)  :: norient_face
!
!..face order
   integer, dimension(5)  :: norderf
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
   real(8), dimension(3)   :: xi,x,rn,daux
   real(8), dimension(3,2) :: dxidt,dxdt,rt
   real(8), dimension(3,3) :: dxdxi,dxidx
   real(8), dimension(2)   :: t
!
!..H1 shape functions
   real(8), dimension(MAXbrickH)   :: shapH
   real(8), dimension(3,MAXbrickH) :: gradH
!
!..H(curl) shape functions
   real(8), dimension(3,MAXbrickE) :: shapE
   real(8), dimension(3,MAXbrickE) :: curlE
!
!..L2 shape functions
   real(8), dimension(MAXbrickQ) :: shapQ
!
!..enriched Hcurl shape functions
   real(8), dimension(3,MAXbrickEE) :: shapEE
   real(8), dimension(3,MAXbrickEE) :: curlEE
!
!..Gram matrix in packed format
   !VTYPE, dimension(NrTest*(NrTest+1)/2) :: gramTest
   VTYPE, allocatable :: gramP(:), Grfp(:)
   real(8) :: FF, CF, FC
   real(8) :: fldE(3), fldH(3), crlE(3), crlH(3), rotE(3)
   real(8) :: fldF(3), fldG(3), crlF(3), crlG(3), rotF(3)
!
   real(8) :: D_aux(3,3),D_za(3,3),D_zc(3,3)
!
!..load vector for the enriched space
   VTYPE, dimension(NrTest)   :: bload_E,bload_Ec
   VTYPE, dimension(2*NrdofE) :: bload_Imp
!
!..3D quadrature data
   real(8), dimension(3,MAXNINT3ADD) :: xiloc
   real(8), dimension(MAXNINT3ADD)   :: waloc
!
!..2D quadrature data
   real(8), dimension(2,MAXNINT2ADD) :: tloc
   real(8), dimension(MAXNINT2ADD)   :: wtloc
!
!..BC's flags
   integer, dimension(6,NRINDEX)     :: ibc
!
!..Maxwell load and auxiliary variables
   VTYPE  , dimension(3) :: zJ,zImp
   real(8), dimension(3) :: E1,rntimesE,rn2timesE
   real(8)               :: eps
!
!..approximate solution
   VTYPE, dimension(3,2) :: zsolExi,zsolE,zflux,zflux2
   VTYPE, dimension(6)   :: zsolQ
!
!..auxiliary
   VTYPE :: zresid, zaux, zbux, zcux
!
!..number of faces per element type
   integer :: nrf
!
!..various variables for the problem
   real(8) :: rjac,bjac,weight,wa,CC,EE,CE,E,EC,q,h,tol,diff_r,diff_i,max_r,max_i
   integer :: i1,i2,j1,j2,k1,k2,kH,kk,i,j,m,n,nint,kE,k,l,ivar,iflag
   integer :: nordP,nsign,ifc,ndom,info,icomp,nrdof,nrdof_eig,idec
   VTYPE   :: zfval
   VTYPE   :: za(3,3),zc(3,3)
   VTYPE   :: zaJ(3,3),zcJ(3,3)
!
!..for PML
   VTYPE :: zbeta,zdbeta,zd2beta,detJstretch
   VTYPE, dimension(3,3) :: Jstretch,invJstretch,JJstretch
!
!..OMEGA_RATIO_SIGNAL or OMEGA_RATIO_PUMP
   real(8) :: OMEGA_RATIO_FLD
!..WAVENUM_SIGNAL or WAVENUM_PUMP
   real(8) :: WAVENUM_FLD
!
!..for polarizations function
   VTYPE, dimension(3,3) :: bg_pol,gain_pol,raman_pol,rndotE
   real(8) :: delta_n
   integer :: dom_flag
!..timer
   real(8) :: MPI_Wtime,start_time,end_time
!
!..debug variables
#if DEBUG_MODE
   integer :: iprint = 0
#endif
!
!..for Gram matrix compressed storage format
   integer :: nk
   nk(k1,k2) = (k2-1)*k2/2+k1
!
!--------------------------------------------------------------------------
!
   allocate(gramP(NrTest*(NrTest+1)/2))
!
!..element type
   etype = NODES(Mdle)%type
   nrf = nface(etype)
!
!..determine order of approximation
   call find_order(Mdle, norder)
!..set the enriched order of approximation
   select case(etype)
      case('mdlb'); nordP = NODES(Mdle)%order+NORD_ADD*111
      case('mdln','mdld'); nordP = NODES(Mdle)%order+NORD_ADD
      case('mdlp'); nordP = NODES(Mdle)%order+NORD_ADD*11
   end select
!..determine edge and face orientations
   call find_orient(Mdle, norient_edge,norient_face)
!..determine nodes coordinates
   call nodcor(Mdle, xnod)
!..get the element boundary conditions flags
   call find_bc(Mdle, ibc)
!..get current solution dofs
   call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
!
!..clear space for auxiliary matrices
   bload_E = ZERO; gramP = ZERO; bload_Ec = ZERO; bload_Imp = ZERO
!
!..initialize the background polarization
   call find_domain(Mdle, ndom)
!..select case of GEOM_NO to set
!..refractive index according to domain
!..Fld_flag = 1 when we are in signal element routine
   bg_pol = ZERO
   if(GEOM_NO .eq. 5) then
      x(1:3) = 0.d0
      select case(ndom)
         case(1,2)
            dom_flag = 1 ! Fiber core
            call get_bgPol(dom_flag,Fld_flag,0.d0,x, bg_pol)
         case(3,4)
            dom_flag = 0 ! Fiber cladding
            call get_bgPol(dom_flag,Fld_flag,0.d0,x, bg_pol)
         case default
            write(*,*) 'elem_residual_maxwell: unexpected ndom param. stop.'
            stop
      end select
   endif
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
         write(*,*) 'elem_residual_maxwell. invalid Fld_flag param. stop.'
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
!--------------------------------------------------------------------------
!
!..NORMAL element integrals for LOAD (i.e., l-Bu)
!
!--------------------------------------------------------------------------
!
!..use the enriched order to set the quadrature
   INTEGRATION = NORD_ADD
   call set_3D_int_DPG(etype,norder,norient_face, nint,xiloc,waloc)
   INTEGRATION = 0
!
!..loop over
   do l=1,nint
      xi(1:3) = xiloc(1:3,l)
      wa = waloc(l)
!  ...determine element H1 shape functions
      call shape3DH(etype,xi,norder,norient_edge,norient_face,  &
                    nrdof,shapH,gradH)
#if DEBUG_MODE
      if (nrdof .ne. NrdofH) then
         write(*,*) 'elem_residual_maxwell: INCONSISTENCY NrdofH. stop.'
         stop
      endif
#endif
!  ...determine element H(curl) shape functions
      call shape3DE(etype,xi,norder,norient_edge,norient_face, &
                    nrdof,shapE,curlE)
#if DEBUG_MODE
      if (nrdof .ne. NrdofE) then
         write(*,*) 'elem_residual_maxwell: INCONSISTENCY NrdofE. stop.'
         stop
      endif
#endif
!  ...determine element L2 shape functions
      call shape3DQ(etype,xi,norder, nrdof,shapQ)
#if DEBUG_MODE
      if (nrdof .ne. NrdofQ) then
         write(*,*) 'elem_residual_maxwell: INCONSISTENCY NrdofQ. stop.'
         stop
      endif
#endif
!  ...determine discontinuous H(curl) shape functions
      call shape3EE(etype,xi,nordP, nrdof,shapEE,curlEE)
#if DEBUG_MODE
      if (nrdof .ne. NrdofEE) then
         write(*,*) 'elem_residual_maxwell: INCONSISTENCY NrdofEE. stop.'
         stop
      endif
#endif
!  ...geometry
      call geom3D(Mdle,xi,xnod,shapH,gradH,NrdofH, &
                   x,dxdxi,dxidx,rjac,iflag)
!  ...integration weight
      weight = rjac*wa
!  ...compute the approximate solution
      zsolQ = ZERO
!
      do k=1,NrdofQ
         if(NO_PROBLEM.eq.3) then
            zsolQ(1:6)  = zsolQ(1:6)  + zdofQ(1:6,k)*shapQ(k)
         elseif(NO_PROBLEM.eq.4) then
            zsolQ(1:6)  = zsolQ(1:6)  + zdofQ(7:12,k)*shapQ(k)
         else
            write(*,*) 'elem_residual_maxwell: NO_PROBLEM must be 3 or 4. stop.'
         endif
      enddo
      zsolQ = zsolQ/rjac
!
!  ...get the RHS
!  ...zfval (heat eqn rhs), zJ (maxwell rhs)
      call getf(Mdle,x, zfval,zJ)
!
!  ...set auxiliary constants (updated below if nonlinear)
!  ...update bgpol depending on coordinates if refractive index is varying
!     (other than step-index)
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
               call get_activePol(zsolQ_soleval(1:12),Fld_flag,delta_n, gain_pol)
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
!  ...loop through enriched H(curl) test functions
      do k1=1,NrdofEE
!
!     ...Piola transformation
         do i = 1,3
            daux = dxdxi(i,1:3)
            call dot_product(shapEE(1:3,k1),dxidx(1:3,i), fldF(i))
            call dot_product(curlEE(1:3,k1),daux        , crlF(i))
         enddo
         crlF(1:3) = crlF(1:3)/rjac
         fldG = fldF; crlG = crlF
!     ...e_z x F
         rotF = 0.d0; rotF(1) = -fldF(2); rotF(2) = fldF(1)
!
!     ...accumulate for the load
!     ...first eqn
!        (J^imp,F)
         k = 2*k1-1
         zaux = zJ(1)*fldF(1) + zJ(2)*fldF(2) + zJ(3)*fldF(3)
         bload_E(k) = bload_E(k) + zaux*weight
!
!     ...first eqn
         k = 2*k1-1
!        - ( (H,curl(F)) - iωε(E,F) )
         zaux = crlF(1)*zsolQ(4) + crlF(2)*zsolQ(5) + crlF(3)*zsolQ(6)
         zcux = zaJ(1,1)*fldF(1)*zsolQ(1) + &
                zaJ(2,2)*fldF(2)*zsolQ(2) + &
                zaJ(3,3)*fldF(3)*zsolQ(3)
         bload_E(k) = bload_E(k) - (zaux - zcux)*weight
!     ...additional stiffness contribution if solving vectorial envelope equation
         if (ENVELOPE) then
!        ...-( -ik(e_z x H,F) ), where e_z x H = (-H_y,H_x,0)
            zaux = -ZI*WAVENUM_FLD*detJstretch*(-fldF(1)*zsolQ(5) + fldF(2)*zsolQ(4))
            bload_E(k) = bload_E(k) - zaux*weight
         endif
!
!     ...second eqn
         k = 2*k1
!        - ( (E,curl(G)) + iωμ(H,G) )
         zaux = crlG(1)*zsolQ(1) + crlG(2)*zsolQ(2) + crlG(3)*zsolQ(3)
         zcux = zcJ(1,1)*fldG(1)*zsolQ(4) + &
                zcJ(2,2)*fldG(2)*zsolQ(5) + &
                zcJ(3,3)*fldG(3)*zsolQ(6)
         bload_E(k) = bload_E(k) - (zaux + zcux)*weight
!     ...additional stiffness contribution if solving vectorial envelope equation
         if (ENVELOPE) then
!        ...-( -ik(e_z x E,G) ), where e_z x E = (-E_y,E_x,0)
            zaux = -ZI*WAVENUM_FLD*detJstretch*(-fldG(1)*zsolQ(2) + fldG(2)*zsolQ(1))
            bload_E(k) = bload_E(k) - zaux*weight
         endif
!
! ===============================================================================
!     ...Computation of Gram Matrix (w/o fast integration)
!     ...loop through enriched H(curl) test functions
         if (FAST_INT .eq. 1 .and. etype .eq. 'mdlb') cycle
         if (FAST_INT .eq. 1 .and. etype .eq. 'mdlp') cycle
         do k2=k1,NrdofEE
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
!           -------------------------
!           | (F_i,F_j)   (F_i,G_j) |
!           | (G_i,F_j)   (G_i,G_j) |
!           -------------------------
!           (F_i,F_j) terms
            n = 2*k1-1; m = 2*k2-1
            k = nk(n,m)
            zaux = (abs(zaJ(1,1))**2)*fldF(1)*fldE(1) + &
                   (abs(zaJ(2,2))**2)*fldF(2)*fldE(2) + &
                   (abs(zaJ(3,3))**2)*fldF(3)*fldE(3)
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
!           (F_i,G_j) terms
            n = 2*k1-1; m = 2*k2
            k = nk(n,m)
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
!              (G_i,F_j) terms
               n = 2*k1; m = 2*k2-1
               k = nk(n,m)
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
!           (G_i,G_j) terms
            n = 2*k1; m = 2*k2
            k = nk(n,m)
            zcux = (abs(zcJ(1,1))**2)*fldF(1)*fldE(1) + &
                   (abs(zcJ(2,2))**2)*fldF(2)*fldE(2) + &
                   (abs(zcJ(3,3))**2)*fldF(3)*fldE(3)
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
         enddo
      enddo
   enddo
!
!   start_time = MPI_Wtime()
   if (FAST_INT .eq. 1) then
      select case (etype)
         case('mdlb')
            call elem_maxwell_gram_hexa(Mdle,Fld_flag,NrTest,NrdofH, gramP)
         case('mdlp')
            call elem_maxwell_gram_pris(Mdle,Fld_flag,NrTest,NrdofH, gramP)
      end select
!..end timer
!   end_time = MPI_Wtime()
!   !$OMP CRITICAL
!   write(*,10) etype, end_time-start_time
!   10 format(A,' gram : ',f12.5,'  seconds')
!   !$OMP END CRITICAL
!
!!
!      do k=1,NrTest*(NrTest+1)/2
!!      do k1=1,NrTest
!!      do k2=k1,NrTest
!!         k = nk(k1,k2)
!         if (.not. cmp_dc(gramP(k),gramTest(k))) then
!!           if (mod(k1,2) .eq. 1 .and. mod(k2,2) .eq. 1) cycle
!            if (mod(k1,2) .eq. 0 .and. mod(k2,2) .eq. 0) cycle
!            diff_r = abs(real(gramP(k))-real(gramTest(k)))
!            diff_i = abs(imag(gramP(k))-imag(gramTest(k)))
!            !$OMP CRITICAL
!            write(*,*) 'etype    = ', etype
!            write(*,*) 'Mlde     = ', Mdle
!            write(*,*) 'NrTest   = ', NrTest
!            write(*,*) 'gramSize = ', NrTest*(NrTest+1)/2
!!            write(*,1505) 'gramP .ne. gramTest, (k1,k2), k = ', k1,k2,k
!            write(*,1510) 'gramP   : ',real(gramP(k))   ,' , ',imag(gramP(k))
!            write(*,1510) 'gramTest: ',real(gramTest(k)),' , ',imag(gramTest(k))
!            write(*,1510) 'diff r,i: ',diff_r,diff_i
!            write(*,*) ''
!       1505 format(A,'('I4,','I4,'), ',I7)
!       1510 format(A,es15.8,A,es15.8)
!            !$OMP END CRITICAL
!            exit
!         endif
!!      enddo
!      enddo
   endif
!
!--------------------------------------------------------------------------
!
!              B O U N D A R Y      I N T E G R A L S
!
!--------------------------------------------------------------------------
!
!
!..auxiliary constant for flux computation
   if(NO_PROBLEM.eq.3) then
      i = 0
   elseif(NO_PROBLEM.eq.4) then
      i = 2
   else
      write(*,*) 'elem_residual_maxwell: NO_PROBLEM must be 3 or 4. stop.'
      stop
   endif
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
!     ...determine discontinuous H(curl) shape functions
         call shape3EE(etype,xi,nordP, NrdofEE,shapEE,curlEE)
!
!     ...determine element H1 shape functions (for geometry)
         call shape3DH(etype,xi,norder,norient_edge,norient_face, &
                       nrdof,shapH,gradH)
#if DEBUG_MODE
         if (nrdof .ne. NrdofH) then
            write(*,*) 'elem_residual_maxwell: INCONSISTENCY NrdofH. stop.'
            stop
         endif
#endif
!
!     ...determine element H(curl) shape functions (for fluxes)
         call shape3DE(etype,xi,norder,norient_edge,norient_face, &
                       nrdof,shapE,curlE)
#if DEBUG_MODE
         if (nrdof .ne. NrdofE) then
            write(*,*) 'elem_residual_maxwell: INCONSISTENCY NrdofE. stop.'
            stop
         endif
#endif
!
!     ...geometry
         call bgeom3D(Mdle,xi,xnod,shapH,gradH,NrdofH,dxidt,nsign, &
                      x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
         weight = bjac*wtloc(l)
!
!     ...compute approximate fluxes at the point
         zsolExi = ZERO
!
         do ivar=1,2
            do k=1,NrdofE
               zsolExi(1:3,ivar) = zsolExi(1:3,ivar) &
                                 + zdofE(ivar+i,k)*shapE(1:3,k)
            enddo
            zsolE(1:3,ivar) = zsolExi(1,ivar)*dxidx(1,1:3) &
                            + zsolExi(2,ivar)*dxidx(2,1:3) &
                            + zsolExi(3,ivar)*dxidx(3,1:3)
            call zcross_product(rn,zsolE(1:3,ivar), zflux (1:3,ivar))
            call zcross_product(rn,zflux(1:3,ivar), zflux2(1:3,ivar))
         enddo
!
!     ...check for impedance BC (elimination strategy)
         if ((ibc(ifc,3).eq.3 .and. Fld_flag.eq.1) .or. &
             (ibc(ifc,5).eq.3 .and. Fld_flag.eq.0)) then
!        ...impedance surface load [zImp should be zero here]
            call get_bdSource(Mdle,x,rn, zImp)
            zflux2(1:3,1) = GAMMA*zflux2(1:3,1)+zImp
         endif
!
!     ...loop through enriched test functions
         do k1=1,NrdofEE
            E1(1:3) = shapEE(1,k1)*dxidx(1,1:3) &
                    + shapEE(2,k1)*dxidx(2,1:3) &
                    + shapEE(3,k1)*dxidx(3,1:3)
!
            k=2*k1-1
!        ...check for impedance BC (elimination strategy)
            if ((ibc(ifc,3).eq.3 .and. Fld_flag.eq.1) .or. &
                (ibc(ifc,5).eq.3 .and. Fld_flag.eq.0)) then
!           - GAMMA * < n x n x E , G >
               zaux = E1(1)*zflux2(1,1) + E1(2)*zflux2(2,1) + E1(3)*zflux2(3,1)
               bload_E(k) = bload_E(k) - zaux * weight
            else
!           - <n x H, F>
               zaux = E1(1)*zflux(1,2) + E1(2)*zflux(2,2) + E1(3)*zflux(3,2)
               bload_E(k) = bload_E(k) - zaux * weight
            endif
!           - <n x E, G>
            k = 2*k1
            zaux = E1(1)*zflux(1,1) + E1(2)*zflux(2,1) + E1(3)*zflux(3,1)
            bload_E(k) = bload_E(k) - zaux * weight
         enddo
!
!     ...check for impedance BC (L2 penalty method)
         if ((ibc(ifc,3).ne.2 .or. Fld_flag.ne.1) .and. &
             (ibc(ifc,5).ne.2 .or. Fld_flag.ne.0)) cycle
!     ...impedance surface load [zImp should be zero here]
         call get_bdSource(Mdle,x,rn, zImp)
!     ...define the weight of the penalty term
         eps = 1.d0
!     ...compute residual contribution from impedance boundary
         do k1=1,NrdofE
            E1(1:3) = shapE(1,k1)*dxidx(1,1:3) &
                    + shapE(2,k1)*dxidx(2,1:3) &
                    + shapE(3,k1)*dxidx(3,1:3)
!
            call cross_product(rn,E1, rntimesE)
            call cross_product(rn,rntimesE, rn2timesE)
!
!        ...1st test function
            k = 2*k1-1
!           + GAMMA^2 * < n x n x E , n x n x F >
            bload_Imp(k) = bload_Imp(k) + (rn2timesE(1)*zflux2(1,1)  &
                                        +  rn2timesE(2)*zflux2(2,1)  &
                                        +  rn2timesE(3)*zflux2(3,1)  &
                                          )*GAMMA*GAMMA*weight/eps
!           - GAMMA * < n x H, n x n x F >
            bload_Imp(k) = bload_Imp(k) - (rn2timesE(1)*zflux(1,2)   &
                                        +  rn2timesE(2)*zflux(2,2)   &
                                        +  rn2timesE(3)*zflux(3,2)   &
                                          )*GAMMA*weight/eps
!           + GAMMA * < zImp , n x n x F >
            bload_Imp(k) = bload_Imp(k) + (rn2timesE(1)*zImp(1)  &
                                        +  rn2timesE(2)*zImp(2)  &
                                        +  rn2timesE(3)*zImp(3)  &
                                          )*GAMMA*weight/eps
!        ...2nd test function
            k = 2*k
!           - GAMMA * < n x n x E , n x G >
            bload_Imp(k) = bload_Imp(k) - (rntimesE(1)*zflux2(1,1)  &
                                        +  rntimesE(2)*zflux2(2,1)  &
                                        +  rntimesE(3)*zflux2(3,1)  &
                                          )*GAMMA*weight/eps
!           + < n x H , n x G >
            bload_Imp(k) = bload_Imp(k) + (rntimesE(1)*zflux(1,2)   &
                                        +  rntimesE(2)*zflux(2,2)   &
                                        +  rntimesE(3)*zflux(3,2)   &
                                          )*weight/eps
!           - < zImp , n x G >
            bload_Imp(k) = bload_Imp(k) - (rntimesE(1)*zImp(1)  &
                                        +  rntimesE(2)*zImp(2)  &
                                        +  rntimesE(3)*zImp(3)  &
                                          )*weight/eps
         enddo
!
      enddo
   enddo
!
#if DEBUG_MODE
   if (iprint.gt.0) then
      write(*,7015) bload_E(1:2*NrdofEE)
 7015 format('elem_residual_maxwell: FINAL bload_E = ',10(/,6(2e12.5,2x)))
      call pause
   endif
#endif
!
!--------------------------------------------------------------------------
!
!..Convert Gram matrix from packed to RFP format (same storage cost
!  but RFP can use BLAS3 so is faster)
   allocate(Grfp(NrTest*(Nrtest+1)/2))
   call ZTPTTF('N','U',NrTest,gramP,Grfp,info)
   if (info.ne.0) then
      write(*,*) 'elem_residual_maxwell: ZTPTTF: Mdle,info = ',Mdle,info,'. stop.'
      stop
   endif
   deallocate(gramP)
!
!..factorize the test Gram matrix
   call ZPFTRF('N','U',NrTest,Grfp,info)
!   call ZPPTRF('U', NrTest, gramP, info)
   if (info.ne.0) then
      write(*,*) 'elem_residual_maxwell: ZPFTRF: Mdle,info = ',Mdle,info,'. stop.'
      stop
   endif
!
!..save copies of the RHS to compute later the residual
   bload_Ec = bload_E
!
!
!..compute the product of inverted test Gram matrix with RHS,
!..bload_E is overwritten with the solution
   call ZPFTRS('N','U',NrTest,1,Grfp,bload_E,NrTest,info)
!   call ZPPTRS('U', NrTest, 1, gramP, bload_E, NrTest, info)
   if (info.ne.0) then
      write(*,*) 'elem_residual_maxwell: ZPFTRS: Mdle,info = ',Mdle,info,'. stop.'
      stop
   endif
!
   deallocate(Grfp)
!
!..compute the residual
   zresid = ZERO
   do k=1,NrTest
      zresid = zresid + bload_Ec(k)*conjg(bload_E(k))
   enddo
!
!..account for impedance BC penalty term (L2 penalty method)
!  test norm for the residual has then two separate contributions:
!  1) Usual DPG residual measured in adjoint test norm ||\psi||_V
!  2) Additional Impedance BC residual measured in L2 norm ||\phi||
   if (IBCFLAG.eq.2) then
      do k=1,2*NrdofE
         zresid = zresid + bload_Imp(k)*conjg(bload_Imp(k))
      enddo
   endif
!
   Resid = real(zresid,8)
!
!..set suggested refinement flag
   select case(etype)
   !  isotropic
      case('mdlb');        Nref_flag = 111
      case('mdlp');        Nref_flag = 11
      case('mdln','mdld'); Nref_flag = 1
   end select
!
#if DEBUG_MODE
   if (iprint.eq.1) then
      write(*,7010) Mdle, Resid
 7010 format('elem_residual_maxwell: Mdle, Resid = ',i5,3x,e12.5)
      call pause
   endif
#endif
!
end subroutine elem_residual_maxwell

