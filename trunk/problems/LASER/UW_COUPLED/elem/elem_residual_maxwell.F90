!-------------------------------------------------------------------------
!
!     routine name         - elem_residual_maxwell
!
!-------------------------------------------------------------------------
!
!     latest revision:     - Apr 2019
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
#include "implicit_none.h"
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
   real*8 , intent(out) :: Resid
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
   real*8, dimension(3,MAXbrickH) :: xnod
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
   real*8 :: rsolH
!
!..variables for geometry
   real*8, dimension(3)   :: xi,x,rn
   real*8, dimension(3,2) :: dxidt,dxdt,rt
   real*8, dimension(3,3) :: dxdxi,dxidx
   real*8, dimension(2)   :: t
!
!..H1 shape functions
   real*8, dimension(MAXbrickH)   :: shapH
   real*8, dimension(3,MAXbrickH) :: gradH
!
!..H(curl) shape functions
   real*8, dimension(3,MAXbrickE) :: shapE
   real*8, dimension(3,MAXbrickE) :: curlE
!
!..L2 shape functions
   real*8, dimension(MAXbrickQ) :: shapQ
!
!..enriched Hcurl shape functions
   real*8 , dimension(3,MAXbrickEE) :: shapEE
   real*8 , dimension(3,MAXbrickEE) :: curlEE
!
!..Gram matrix in packed format
   VTYPE, dimension(NrTest*(NrTest+1)/2) :: gramP !, gramTest
   real*8  :: FF, CF, FC
   real*8  :: fldE(3), fldH(3), crlE(3), crlH(3)
   real*8  :: fldF(3), fldG(3), crlF(3), crlG(3)
!
!..load vector for the enriched space
   VTYPE, dimension(NrTest) :: bload_E,bload_Ec
!
!..3D quadrature data
   real*8, dimension(3,MAXNINT3ADD) :: xiloc
   real*8, dimension(MAXNINT3ADD)   :: waloc
!
!..2D quadrature data
   real*8, dimension(2,MAXNINT2ADD) :: tloc
   real*8, dimension(MAXNINT2ADD)   :: wtloc
!
!..BC's flags
   integer, dimension(6,NR_PHYSA)   :: ibc
!
!..Maxwell load and auxiliary variables
   VTYPE , dimension(3) :: zJ,zImp
   real*8, dimension(3) :: E1,E2,rntimesE
!
!..approximate solution
   VTYPE, dimension(3,2) :: zsolExi,zsolE,zflux,zflux2
   VTYPE, dimension(6)   :: zsolQ
!
!..
   VTYPE :: zresid, zaux
!
!..number of faces per element type
   integer :: nrf
!
!..various variables for the problem
   real*8  :: rjac,bjac,weight,wa,CC,EE,CE,E,EC,q,h
   integer :: i1,i2,j1,j2,k1,k2,kH,kk,i,j,m,n,nint,kE,k,iprint,l,ivar,iflag
   integer :: nordP,nsign,ifc,ndom,info,icomp,nrdof,nrdof_eig,idec
   VTYPE   :: zfval,zc
   VTYPE   :: za(3,3),zb(3,3)
!
!..OMEGA_RATIO_SIGNAL or OMEGA_RATIO_PUMP
   real*8  :: OMEGA_RATIO_FLD
!
!..for polarizations function
   VTYPE, dimension(3,3) :: bg_pol,gain_pol,raman_pol,rndotE
   real*8  :: delta_n
   integer :: dom_flag
!
!..for Gram matrix compressed storage format
   integer :: nk
   nk(k1,k2) = (k2-1)*k2/2+k1
!
!--------------------------------------------------------------------------
!
   iprint = 0
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
   bload_E = ZERO; gramP = ZERO; bload_Ec = ZERO
!
!..initialize the background polarization
   call find_domain(Mdle, ndom)
!..select case of GEOM_NO to set
!..refractive index according to domain
!..Fld_flag = 1 when we are in signal element routine
   bg_pol = ZERO
   if(GEOM_NO .eq. 5) then
      select case(ndom)
         case(1,2)
            dom_flag = 1 ! Fiber core
            call get_bgPol(dom_flag,Fld_flag,0.d0, bg_pol)
         case(3,4)
            dom_flag = 0 ! Fiber cladding
            call get_bgPol(dom_flag,Fld_flag,0.d0, bg_pol)
         case default
            write(*,*) 'elem_residual_maxwell: unexpected ndom param. stop.'
            stop
      end select
   endif
!
!..set OMEGA_RATIO_FLD
   select case(Fld_flag)
      case(0)
         OMEGA_RATIO_FLD = OMEGA_RATIO_PUMP
      case(1)
         OMEGA_RATIO_FLD = OMEGA_RATIO_SIGNAL ! 1.0d0
      case default
         write(*,*) 'elem_residual_maxwell. invalid Fld_flag param. stop.'
         stop
   end select
!
!..set auxiliary constants (updated if needed in integration routine)
   za = (ZI*OMEGA*OMEGA_RATIO_FLD*EPSILON + SIGMA)*IDENTITY + bg_pol
   zb = -conjg(za)
   zc = ZI*OMEGA*OMEGA_RATIO_FLD*MU
!
!--------------------------------------------------------------------------
!
!..NORMAL element integrals for LOAD (i.e., l-Bu)
!
!--------------------------------------------------------------------------
!
!..use the enriched order to set the quadrature
   INTEGRATION = NORD_ADD
   call set_3Dint_DPG(etype,norder, nint,xiloc,waloc)
   INTEGRATION = 0
!
!..loop over
   do l=1,nint
      xi(1:3) = xiloc(1:3,l)
      wa = waloc(l)
!  ...determine element H1 shape functions
      call shape3H(etype,xi,norder,norient_edge,norient_face,  &
                   nrdof,shapH,gradH)
#if DEBUG_MODE
      if (nrdof .ne. NrdofH) then
         write(*,*) 'elem_residual_maxwell: INCONSISTENCY NrdofH. stop.'
         stop
      endif
#endif
!  ...determine element H(curl) shape functions
      call shape3E(etype,xi,norder,norient_edge,norient_face, &
                    nrdof,shapE,curlE)
#if DEBUG_MODE
      if (nrdof .ne. NrdofE) then
         write(*,*) 'elem_residual_maxwell: INCONSISTENCY NrdofE. stop.'
         stop
      endif
#endif
!  ...determine element L2 shape functions
      call shape3Q(etype,xi,norder, nrdof,shapQ)
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
            call get_bgPol(dom_flag,Fld_flag,delta_n, bg_pol)
         endif
!
!     ...initialize gain polarization, raman polarization
         gain_pol = ZERO; raman_pol = ZERO
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
!     ...endif RAMAN_GAIN
         endif
!     ...update auxiliary constants for za,zb: this is for
!     ...Stiffness and Gram matrix that changes with each nonlinear iteration
         za = (ZI*OMEGA*OMEGA_RATIO_FLD*EPSILON+SIGMA)*IDENTITY+bg_pol+gain_pol+raman_pol
         zb = -conjg(za)
!  ...endif NONLINEAR_FLAG
      endif
!
!  ...loop through enriched H(curl) test functions
      do k1=1,NrdofEE
!
!     ...Piola transformation
         do i = 1,3
            call dot_product(shapEE(1:3,k1),dxidx(1:3,i), fldF(i))
            call dot_product(curlEE(1:3,k1),dxdxi(i,1:3), crlF(i))
         enddo
         crlF(1:3) = crlF(1:3)/rjac
         fldG = fldF; crlG = crlF
         !n = 2*k1-1
         !bload_E(n) = bload_E(n)                                   &
         !           + (fldG(1)*zJ(1)+fldG(2)*zJ(2)+fldG(3)*zJ(3))  &
         !           * weight
!
!     ...accumulate for the load
         k = 2*k1-1
         bload_E(k) = bload_E(k)+(zJ(1)*fldF(1)+zJ(2)*fldF(2)+zJ(3)*fldF(3))*weight
!
!     ...order is modified so that A*=-A
         k = 2*k1
         bload_E(k) = bload_E(k) &
          - ((crlF(1)*zsolQ(1)+crlF(2)*zsolQ(2)+crlF(3)*zsolQ(3)) &
             +zc*(fldF(1)*zsolQ(4)+fldF(2)*zsolQ(5)+fldF(3)*zsolQ(6)) &
          )*weight
!
         k = 2*k1-1
         bload_E(k) = bload_E(k) &
          - ((crlF(1)*zsolQ(4)+crlF(2)*zsolQ(5)+crlF(3)*zsolQ(6)) &
             -(za(1,1)*fldF(1)*zsolQ(1)+za(2,2)*fldF(2)*zsolQ(2)+za(3,3)*fldF(3)*zsolQ(3)) &
          )*weight
!
! ===============================================================================
!     ...Computation of Gram Matrix (w/o fast integration)
!     ...loop through enriched H(curl) trial functions
         if (FAST_INT .eq. 1 .and. etype .eq. 'mdlb') cycle
         do k2=k1,NrdofEE
!        ...Piola transformation
            do i = 1,3
               call dot_product(shapEE(1:3,k2),dxidx(1:3,i), fldE(i))
               call dot_product(curlEE(1:3,k2),dxdxi(i,1:3), crlE(i))
            enddo
            crlE(1:3) = crlE(1:3)/rjac
!
            call dot_product(fldF,fldE, FF)
            call dot_product(fldF,crlE, FC)
            call dot_product(crlF,fldE, CF)
            call dot_product(crlF,crlE, CC)
!
!        ...accumulate for the gram matrix
            n = 2*k1-1; m = 2*k2-1
            k = nk(n,m)
            zaux = abs(zb(1,1))**2*fldF(1)*fldE(1) + &
                   abs(zb(2,2))**2*fldF(2)*fldE(2) + &
                   abs(zb(3,3))**2*fldF(3)*fldE(3)
            gramP(k) = gramP(k)            &
                     + (zaux + ALPHA_NORM*FF + CC)*weight
!
            n = 2*k1-1; m = 2*k2
            k = nk(n,m)
            zaux = zb(1,1)*fldF(1)*crlE(1) + &
                   zb(2,2)*fldF(2)*crlE(2) + &
                   zb(3,3)*fldF(3)*crlE(3)
            gramP(k) = gramP(k)            &
                     + (-zaux + conjg(zc)*CF)*weight
!
            if (k1 .ne. k2) then
               n = 2*k1; m = 2*k2-1
               k = nk(n,m)
               zaux = za(1,1)*crlF(1)*fldE(1) + &
                      za(2,2)*crlF(2)*fldE(2) + &
                      za(3,3)*crlF(3)*fldE(3)
               gramP(k) = gramP(k)         &
                        + (zaux + zc*FC)*weight
            endif
!
            n = 2*k1; m = 2*k2
            k = nk(n,m)
            gramP(k) = gramP(k)            &
                     + ((abs(zc)**2 + ALPHA_NORM)*FF + CC)*weight
         enddo
      enddo
   enddo
!
   if (FAST_INT .eq. 1 .and. etype .eq. 'mdlb') then
       call elem_maxwell_gram_fi(Mdle,Fld_flag,NrTest,NrdofH, gramP)
!
!      do k=1,NrTest*(NrTest+1)/2
!         if (abs(real(gramP(k))-real(gramTest(k))) > 1.d-8 .or. &
!             abs(imag(gramP(k))-imag(gramTest(k))) > 1.d-8)  then
!            write(*,*) 'Mlde     = ', Mdle
!            write(*,*) 'NrTest   = ', NrTest
!            write(*,*) 'gramSize = ', NrTest*(NrTest+1)/2
!            write(*,1505) 'gramP .ne. gramTest, k = ', k
!            write(*,1510) 'gramP   : ',real(gramP(k))   ,' , ',imag(gramP(k))
!            write(*,1510) 'gramTest: ',real(gramTest(k)),' , ',imag(gramTest(k))
!       1505 format(A,I6)
!       1510 format(A,es15.8,A,es15.8)
!            exit
!         endif
!      enddo
!
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
      call set_2Dint_DPG(ftype,norderf, nint,tloc,wtloc)
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
         call shape3H(etype,xi,norder,norient_edge,norient_face, &
                      nrdof,shapH,gradH)
#if DEBUG_MODE
         if (nrdof .ne. NrdofH) then
            write(*,*) 'elem_residual_maxwell: INCONSISTENCY NrdofH. stop.'
            stop
         endif
#endif
!
!     ...determine element H(curl) shape functions (for fluxes)
         call shape3E(etype,xi,norder,norient_edge,norient_face, &
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
         if(ibc(ifc,2).eq.9) then
            call get_bdSource(Mdle,x,rn, zImp)
            zflux2(1:3,1) = GAMMA*zflux2(1:3,1)+zImp
         endif
!
!     ...loop through enriched test functions
         do k1=1,NrdofEE
            E1(1:3) = shapEE(1,k1)*dxidx(1,1:3) &
                    + shapEE(2,k1)*dxidx(2,1:3) &
                    + shapEE(3,k1)*dxidx(3,1:3)
!        ...modifying for A*=-A
            k = 2*k1
            bload_E(k) = bload_E(k) &
             - (E1(1)*zflux(1,1)+E1(2)*zflux(2,1)+E1(3)*zflux(3,1)) &
             *weight
!        ...modifying for A*=-A
            k=2*k1-1
            if(ibc(ifc,2).eq.9) then
               bload_E(k) = bload_E(k) &
                - (E1(1)*zflux2(1,1)+E1(2)*zflux2(2,1)+E1(3)*zflux2(3,1)) &
                *weight
            else
               bload_E(k) = bload_E(k) &
                - (E1(1)*zflux(1,2)+E1(2)*zflux(2,2)+E1(3)*zflux(3,2)) &
                *weight
            endif
         enddo
      enddo
      if (iprint.ne.0) pause
   enddo
!
   if (iprint.gt.0) then
      write(*,7015) bload_E(1:2*NrdofEE)
 7015 format('elem_residual_maxwell: FINAL bload_E = ',10(/,6(2e12.5,2x)))
      call pause
   endif
!
!--------------------------------------------------------------------------
!
!..factorize the test Gram matrix
   call ZPPTRF('U', NrTest, gramP, info)
   if (info.ne.0) then
      !write(*,*) 'elem_residual_maxwell: info = ',info
      !stop
   endif
!
!..save copies of the RHS to compute later the residual
   bload_Ec = bload_E
!
!..compute the product of inverted test Gram matrix with RHS,
!..bload_E is overwritten with the solution
   call ZPPTRS('U', NrTest, 1, gramP, bload_E, NrTest, info)
   if (info.ne.0) then
      !write(*,*) 'elem_residual_maxwell: info = ',info
      !stop
   endif
!
!..compute the residual
   zresid = ZERO
   do k=1,NrTest
      zresid = zresid + bload_Ec(k)*conjg(bload_E(k))
   enddo
   Resid = real(zresid,8)
!
!..set suggested refinement flag (001 = z-ref)
   Nref_flag = 1
!
   if (iprint.eq.1) then
      write(*,7010) Mdle, Resid
 7010 format('elem_residual_maxwell: Mdle, Resid = ',i5,3x,e12.5)
      pause
   endif
!
end subroutine elem_residual_maxwell