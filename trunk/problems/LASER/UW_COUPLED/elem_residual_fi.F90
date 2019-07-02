!--------------------------------------------------------------------
!
!     routine name      - elem_residual
!
!--------------------------------------------------------------------
!
!     latest revision:  - July 18
!
!     purpose:          - routine returns element residual (squared)
!                         for the Primal Poisson and UW Time Harmonic
!                         Maxwell equation
!
!     arguments:
!
!     in:
!             Mdle      - an element middle node number, identified
!                         with the element
!     out:
!             Resid     - element residual (squared)
!             Nref_flag - suggested h-refinement flag
!
!---------------------------------------------------------------------
!
#include "implicit_none.h"
!
subroutine elem_residual(Mdle, Resid,Nref_flag)
!
   use CommonParam
   use control
   use data_structure3D
   use element_data
   use parametersDPG
!
   implicit none
!
!..declare input/output variables
   integer, intent(in)  :: Mdle
   real*8 , intent(out) :: Resid
   integer, intent(out) :: Nref_flag
!
!..number of test and trial degrees of freedom
   integer :: nrdofH ,nrdofE ,nrdofV ,nrdofQ
   integer :: nrdofHH,nrdofEE,nrdofVV,nrdofQQ
!..
   integer :: norder(19),norderP(19),nordP
!..
   character(len=4) :: etype
!
   select case(NO_PROBLEM)
      case(1,2)
         call elem_residual_Heat(Mdle, Resid,Nref_flag)
      case(3,4)
         etype = NODES(Mdle)%type
!     ...determine order of approximation
         call find_order(Mdle, norder)
!     ...set the enriched order of appoximation
         select case(etype)
            case('mdlb')
               nordP = NODES(Mdle)%order+NORD_ADD*111
            case('mdln','mdld')
               nordP = NODES(Mdle)%order+NORD_ADD
            case('mdlp')
               nordP = NODES(Mdle)%order+NORD_ADD*11
            case default
               write(*,*) 'invalid etype param. stop.'
               stop
         end select
!     ...note: compute_enriched_order works only for hexa currently
         call compute_enriched_order(nordP, norderP)
!     ...compute nrdof for trial
         call celndof(etype,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
!     ...compute nrdof for test
         call celndof(etype,norderP, nrdofHH,nrdofEE,nrdofVV,nrdofQQ)
!
         call elem_residual_UWMaxwell(Mdle,nrdofEE,nrdofH,nrdofE,nrdofQ, Resid,Nref_flag)
      case default
         write(*,*) 'error in elem_residual: NO_PROBLEM must be 1,2,3 or 4. stop.'
         stop
   end select
!
end subroutine elem_residual
!
!
!--------------------------------------------------------------------
!
!     routine name      - elem_residual_Heat
!
!--------------------------------------------------------------------
!
subroutine elem_residual_Heat(Mdle, Resid,Nref_flag)
!
   use control
   use parametersDPG
   use element_data
   use data_structure3D
   use CommonParam
   use LaserParam
!
   implicit none
!
!..declare input/output variables
   integer,                     intent(in)  :: Mdle
   real*8,                      intent(out) :: Resid
   integer,                     intent(out) :: Nref_flag
!..declare edge/face type variables
   character(len=4) :: etype,ftype
!
!..declare element order, orientation for edges and faces
   integer, dimension(19)  :: norder
   integer, dimension(12)  :: norient_edge
   integer, dimension(6)   :: norient_face
!..face order
   integer, dimension(5)   :: norderf
!
!..geometry dof (work space for nodcor)
   real*8, dimension(3,MAXbrickH) :: xnod
!
!..solution dof (work space for solelm)
   VTYPE, dimension(MAXEQNH,MAXbrickH) :: zdofH
   VTYPE, dimension(MAXEQNE,MAXbrickE) :: zdofE
   VTYPE, dimension(MAXEQNV,MAXbrickV) :: zdofV
   VTYPE, dimension(MAXEQNQ,MAXbrickQ) :: zdofQ
!..current H1 solution
   VTYPE :: zsolH
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
!..enriched H1 shape functions
   real*8, dimension(MAXbrickHH)   :: shapHH
   real*8, dimension(3,MAXbrickHH) :: gradHH
!
!..H(div) shape functions
   real*8, dimension(3,MAXbrickV)  :: shapV
   real*8, dimension(MAXbrickV)    :: divV
!
!..test functions and gradients
   real*8  :: v2n,v1,v2
   real*8, dimension(3) :: dv1,dv2
!
!..nrdof for various spaces
   integer  :: nrdofH,nrdofE,nrdofV,nrdofQ,nrdofHH
!
!..space for DPG Computations (Gram Matrix, Stiffness etc.)
   integer, parameter  :: MAXtestH = MAXbrickHH
!
!..stiffness matrix for the local Riesz H1 matrix in LAPACK format
   VTYPE, dimension(MAXtestH*(MAXtestH+1)/2) :: AP_Heat
!
!..load vector for the enriched space
   VTYPE, dimension(MAXtestH) :: BLOADH,BLOADHc
!
!..3D quadrature data
   real*8, dimension(3,MAXNINT3ADD)  :: xiloc
   real*8, dimension(MAXNINT3ADD)    :: waloc
!
!..2D quadrature data
   real*8, dimension(2,MAXNINT2ADD)  :: tloc
   real*8, dimension(MAXNINT2ADD)    :: wtloc
!
!..BC's flags
   integer, dimension(6,NR_PHYSA)    :: ibc
!
!..approximate solution
   VTYPE, dimension(3) :: zgradHxi,zgradH,zsolVxi,zsolV
   VTYPE               :: zfval,zsolVn,zvalVn
!
!..for debug printing
   VTYPE, dimension(10)  :: aux
!
!..directional contributions to element residual
   real*8, dimension(4) :: residd
   VTYPE , dimension(4) :: zresidd
   real*8, dimension(3) :: nref,gradpsi
!
!..Maxwell load and auxiliary variables
   VTYPE, dimension(3) :: zJ
!
!..error representation function
   VTYPE, dimension(3,2)    :: zpsi_xi,zcurl_xi_psi,zpsi,zcurl_psi
!
!..metric for L2 norm
   real*8, dimension(3,3):: aa
!
!..number of vertices,edge,faces per element type
   integer :: nrv, nre, nrf
!..various variables for the problem
   real*8  :: h_elem,rjac,weight,wa
   real*8  :: bjac
   integer :: i1,i2,j1,j2,k1,k2,kH,kk,i,j,nint,nint3,iflag,kE,k,iprint,l
   integer :: N,nRHS,nordP,nsign,ifc,ndom,info,icomp,nrdof_eig
!
!..for LAPACK eigensolve
   complex*16, allocatable :: Z(:,:), WORK(:)
   real*8, allocatable     :: W(:),   RWORK(:)
   integer, allocatable    :: IWORK(:)
!
!..for Gram matrix compressed storage format
   integer :: nk
   nk(k1,k2) = (k2-1)*k2/2+k1
!
!-----------------------------------------------------------------------
!
   iprint = 0
!
!..element type
   etype = NODES(Mdle)%type
   nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
!
!..determine order of approximation
   call find_order(Mdle, norder)
!
!  ...set the enriched order of appoximation
   select case(etype)
      case('mdlb'); nordP = NODES(Mdle)%order+NORD_ADD*111
      case('mdln','mdld'); nordP = NODES(Mdle)%order+NORD_ADD
      case('mdlp'); nordP = NODES(Mdle)%order+NORD_ADD*11
   end select
!
!..determine edge and face orientations
   call find_orient( Mdle, norient_edge,norient_face)
!
!..determine nodes coordinates
   call nodcor(Mdle, xnod)
!
!..determine solution dof
   call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
!
   if (iprint.eq.1) then
      write(*,7020) xnod(1,1:8),xnod(2,1:8),xnod(3,1:8)
7020  format('elem_residual: xnod  = ',8(f8.3,2x), &
       2(  /,'                       ',8(f8.3,2x)))
      write(*,7030) zdofH(1,1:8),zdofV(1,1:6)
7030  format('elem_residual: zdofH = ',8(e12.5,2x), &
           /,'               zdofV = ',6(e12.5,2x))
   endif
!
!..clear space for auxiliary matrices
   BLOADH = ZERO; AP_Heat = ZERO
!
!-----------------------------------------------------------------------
!
!..element integrals...
!
!..use the enriched order to set the quadrature
   INTEGRATION = NORD_ADD
   call set_3Dint_DPG(etype,norder, nint3,xiloc,waloc)
   INTEGRATION = 0
!
!..loop over
   do l=1,nint3
      xi(1:3) = xiloc(1:3,l)
      wa = waloc(l)
!
!  .....determine element H1 shape functions
      call shape3H(etype,xi,norder,norient_edge,norient_face, &
                    nrdofH,shapH,gradH)
!
!  .....determine discontinuous H1 shape functions
      call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
!
!  .....geometry
      call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH, &
                  x,dxdxi,dxidx,rjac,iflag)
!
!  .....integration weight
      weight = rjac*wa
!
!  .....compute the approximate solution
      zsolH = 0.d0; zgradHxi(1:3) = 0.d0
      do k=1,nrdofH
         zsolH = zsolH + zdofH(1,k)*shapH(k)
         zgradHxi(1:3) = zgradHxi(1:3) + zdofH(1,k)*gradH(1:3,k)
      enddo
      zgradH(1:3) = zgradHxi(1)*dxidx(1,1:3) &
                  + zgradHxi(2)*dxidx(2,1:3) &
                  + zgradHxi(3)*dxidx(3,1:3)
!
!  ...get the RHS
      call getf(Mdle,x, zfval,zJ)
!
!  ...1st loop through enriched H1 test functions
      do k1=1,nrdofHH
         v1 = shapHH(k1)
         dv1(1:3) = gradHH(1,k1)*dxidx(1,1:3) &
                  + gradHH(2,k1)*dxidx(2,1:3) &
                  + gradHH(3,k1)*dxidx(3,1:3)
!
!
!     ...single step of heat equation
         BLOADH(k1) = BLOADH(k1) &
            + (KAPPA*DELTAT*(zgradH(1)*dv1(1)+zgradH(2)*dv1(2)+ &
                          zgradH(3)*dv1(3)) &
            + zsolH*v1-zfval*v1)*weight
!
!     ...2nd loop through enriched H1 test functions
         do k2=k1,nrdofHH
            v2 = shapHH(k2)
            dv2(1:3) = gradHH(1,k2)*dxidx(1,1:3) &
                     + gradHH(2,k2)*dxidx(2,1:3) &
                     + gradHH(3,k2)*dxidx(3,1:3)
!
!           ...accumulate for the test stiffness matrix
            k = nk(k1,k2)
            select case(INNER_PRODUCT)
               case(1)
                  AP_Heat(k) = AP_Heat(k) &
                             + (dv1(1)*dv2(1)+dv1(2)*dv2(2)+dv1(3)*dv2(3) &
                             +v1*v2)*weight
            end select
!     ...end 1st loop through enriched H1 test functions
         enddo
!  ...end 2nd loop through enriched H1 test functions
      enddo
!..end loop through integration points
   enddo
   if (iprint.eq.2) then
      do i=1,10
         do j=1,i-1
            aux(j) = AP_Heat(nk(j,i))
         enddo
         do j=i,10
            aux(j) = AP_Heat(nk(i,j))
         enddo
         write(*,7011) aux
 7011    format(10e12.5)
      enddo
      pause
   endif
   if (iprint.ge.1) then
      write(*,7014) BLOADH(1:nrdofHH)
 7014 format('elem_residual: BLOADH AFTER VOL INT = ', &
              20(/,10(e12.5,2x)))
      pause
   endif
!
!-----------------------------------------------------------------------
!
!..boundary integrals
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
!     ...determine discontinuous H1 shape functions
         call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
!
!     ...determine element H1 shape functions (for geometry)
         call shape3H(etype,xi,norder,norient_edge,norient_face, &
                      nrdofH,shapH,gradH)
!
!     ...determine element Hdiv shape functions (for fluxes)
         call shape3V(etype,xi,norder,norient_face, &
                      nrdofV,shapV,divV)
!
!     ...geometry
         call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign, &
                      x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
         weight = bjac*wtloc(l)
!
!     ...compute approximate flux at the point
         zsolVxi(1:3) = 0.d0
         do k=1,nrdofV
            zsolVxi(1:3) = zsolVxi(1:3) + zdofV(1,k)*shapV(1:3,k)
         enddo
         zsolV(1:3) = (dxdxi(1:3,1)*zsolVxi(1) &
                      +dxdxi(1:3,2)*zsolVxi(2) &
                      +dxdxi(1:3,3)*zsolVxi(3))/rjac
         zsolVn = ZsolV(1)*rn(1)+ZsolV(2)*rn(2)+ZsolV(3)*rn(3)
!
!     ...loop through enriched test functions
         do k1=1,nrdofHH
            v1 = shapHH(k1)
!        ...accumulate for the load vector
            BLOADH(k1) = BLOADH(k1) - KAPPA*DELTAT*zsolVn*v1*weight
         enddo
      enddo
      if (iprint.eq.1) pause
   enddo
   if (iprint.ge.1) then
      write(*,7015) BLOADH(1:nrdofHH)
7015  format('elem_residual: FINAL BLOADH = ',10(/,10(e12.5,2x)))
      pause
   endif
!
!-----------------------------------------------------------------------
!
!..factorize the test stiffness matrix
   call ZPPTRF('U', nrdofHH, AP_Heat, info)
   if (info.ne.0) then
      write(*,*) 'elem_dpgH1: info = ',info
      stop
   endif
!
!..save copies of the RHS to compute later the residual
   BLOADHc = BLOADH
!
!..compute the product of inverted test Gram matrix with RHS,
!..BLOADH is overwritten with the solution
   call ZPPTRS('U', nrdofHH, 1, AP_Heat, BLOADH, MAXbrickHH, info)
   if (info.ne.0) then
      write(*,*) 'elem_dpgH1: info = ',info
      stop
   endif
!
!..compute the residual
   Resid = 0.d0
   do k=1,NRDOFHH
      Resid = Resid + BLOADHc(k)*BLOADH(k)
   enddo
   Nref_flag = 111
   if (iprint.ge.1) then
      write(*,7010) Mdle, Resid
 7010 format('elem_residual: Mdle, Resid = ',i5,3x,e12.5)
      call pause
   endif
!
end subroutine elem_residual_Heat
!
!
!--------------------------------------------------------------------
!
!     routine name      - elem_residual_UWMaxwell
!
!--------------------------------------------------------------------
!
subroutine elem_residual_UWMaxwell(Mdle,NrdofEE,NrdofH,NrdofE,NrdofQ, &
                              Resid,Nref_flag)
!
   use control
   use parametersDPG
   use element_data
   use data_structure3D
   use CommonParam
   use LaserParam
!
   implicit none
!
!..declare input/output variables
   integer, intent(in)  :: Mdle
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
!..Enriched Hcurl shape functions
   real*8 , dimension(3,MAXbrickEE) :: shapEE
   real*8 , dimension(3,MAXbrickEE) :: curlEE
!
!..space for DPG Computations (Gram Matrix, Stiffness etc.)
   integer, parameter :: MAXtestE = 2*MAXbrickEE
!
!..stiffness matrix for the local Riesz H1 matrix in LAPACK format
   VTYPE, dimension(2*NrdofEE*(2*NrdofEE+1)/2) :: AP_Maxwell
!
!..load vector for the enriched space
   VTYPE, dimension(2*NrdofEE) :: BLOADE,BLOADEc
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
   real*8, dimension(3) :: qq,p,rntimesp,rn2timesp
   real*8, dimension(3) :: E1,curlE1,E2,curlE2,rntimesE
!
!..approximate solution
   VTYPE, dimension(3,2) :: zsolExi,zsolE,zflux,zflux2
   VTYPE, dimension(6)   :: zsolQ
   VTYPE, dimension(3,2) :: zpsi_xi,zcurl_xi_psi,zpsi,zcurl_psi
   VTYPE, dimension(3,3) :: aa
!
!..
   VTYPE :: zresid,Residual
!
!..number of vertices,edge,faces per element type
   integer :: nrv, nre, nrf, iflag
!
!..various variables for the problem
   real*8  :: h_elem,rjac,weight,wa,v2n,CC,EE,CE,E,EC,q,h,omeg,alpha_scale
   real*8  :: bjac,impedanceConstant
   integer :: i1,i2,j1,j2,k1,k2,kH,kk,i,j,nrTest,nint,nint3,kE,k,iprint,l,ivar
   integer :: N,nRHS,nordP,nsign,ifc,ndom,info,icomp,nrdof_eig,idec
   VTYPE   :: zfval,za,zb,zc,zk2
!
!..auxiliary tensors for fast integration
   VTYPE, allocatable :: AUXEE_A(:,:,:)
   VTYPE, allocatable :: AUXEE_B(:,:,:,:)
   VTYPE, allocatable :: AUXCC_A(:,:,:,:,:)
   VTYPE, allocatable :: AUXCC_B(:,:,:,:,:,:)
   VTYPE, allocatable :: AUXEC_A(:,:,:,:)
   VTYPE, allocatable :: AUXCE_A(:,:,:,:)
   VTYPE, allocatable :: AUXEC_B(:,:,:,:,:)
   VTYPE, allocatable :: AUXCE_B(:,:,:,:,:)
!..auxiliary variables for fast integration
   integer :: a,b,sa,sb,sc,alph,beta
   integer :: px,py,pz
   integer :: l1,l2,l3,i3,j3,k3,idxbeta,idxalph,idxa,idxa2,idxa3,idxb,idxb2,idxb3,m1,m2
   integer :: nord1,nord2,nord3,nintx,ninty,nintz
   integer :: nrdofH1,nrdofH2,nrdofH3
   real*8  :: xi1,xi2,xi3,wt1,wt2,wt3
   real*8  :: wt123,weighthh,weightvv
   real*8, dimension(MAXPP+1) :: xilocx,xilocy,xilocz
   real*8, dimension(MAXPP+1) :: wlocx,wlocy,wlocz
   real*8, dimension(3,MAXNINT3ADD) :: wloc3
   real*8, dimension(3) :: xip,dHdx,dHHdx
   real*8, dimension(3,3) :: D,C
   real*8, dimension(MAXPP+1,2) :: shapH1,shapH2,shapH3
   real*8, dimension(MAXPP+1,MAXPP+1) :: sH2p,sH3p,dsH2p,dsH3p
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
   iprint = 0
!
!..element type
   etype = NODES(Mdle)%type
   nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
!
   nrTest = 2*NrdofEE
!
!..determine order of approximation
   call find_order(Mdle, norder)
!..set the enriched order of appoximation
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
!..determine solution dof
   call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
!
   if (iprint.eq.1) then
      write(*,7020) xnod(1,1:nrv)
 7020 format('elem_residual: xnod  = ',8(f8.3,2x))
      write(*,7025) xnod(2,1:nrv)
      write(*,7025) xnod(3,1:nrv)
 7025 format('                       ',8(f8.3,2x))
      write(*,7030) 1,zdofE(1,1:nre)
      write(*,7030) 2,zdofE(2,1:nre)
 7030 format('elem_residual: zdofE(',i1',*) = ',2(/,6(2e12.5,2x)))
      pause
   endif
!
!..clear space for auxiliary matrices
   BLOADE = ZERO; AP_Maxwell = ZERO; BLOADEc = ZERO
!
!..TODO what about the polarization terms in the nonlinear case??
!..see elem_fi for updates of these terms in each nonlinear iteration.
!..auxiliary constants
   select case(NO_PROBLEM)
      case(3)
         za = ZI*OMEGA*OMEGA_RATIO_SIGNAL*EPSILON + SIGMA
         zb = -conjg(za)
         zc = ZI*OMEGA*OMEGA_RATIO_SIGNAL*MU
      case(4)
         za = ZI*OMEGA*OMEGA_RATIO_PUMP*EPSILON + SIGMA
         zb = -conjg(za)
         zc = ZI*OMEGA*OMEGA_RATIO_PUMP*MU
      case default
         write(*,*) 'error in elem_residual_UWMaxwell: NO_PROBLEM must be 3 or 4. stop.'
         stop
   end select
!
!-----------------------------------------------------------------------
!
!..SUM FACTORIZATION element integrals for Gram Matrix
!
!-----------------------------------------------------------------------
!
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
   D=ZERO
   C=ZERO
!
!..use the enriched order to set the quadrature
   INTEGRATION = NORD_ADD
   call set_3Dint_fi(etype,norder,nord1,nord2,nord3,nintx,ninty,nintz,xiloc,wloc3)
   INTEGRATION = 0
!..set up # dof for each direction for 1D H1 test functions with order p+dp
   nrdofH1=nord1+1; nrdofH2=nord2+1; nrdofH3=nord3+1
!..Allocate the auxiliary arrays for sum factorization
   allocate(AUXEE_A(3,3,nrdofH3**2))
   allocate(AUXEE_B(3,3,nrdofH2**2,nrdofH3**2))
   allocate(AUXCC_A(2,2,3,3,nrdofH3**2))
   allocate(AUXCC_B(2,2,3,3,nrdofH2**2,nrdofH3**2))
   allocate(AUXEC_A(2,3,3,nrdofH3**2))
   allocate(AUXEC_B(2,3,3,nrdofH2**2,nrdofH3**2))
   allocate(AUXCE_A(2,3,3,nrdofH3**2))
   allocate(AUXCE_B(2,3,3,nrdofH2**2,nrdofH3**2))
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
      AUXEE_B =ZERO
      AUXCC_B =ZERO
      AUXEC_B =ZERO
      AUXCE_B =ZERO
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
         AUXEE_A =ZERO
         AUXCC_A =ZERO
         AUXEC_A =ZERO
         AUXCE_A =ZERO
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
            call shape3H(etype,xip,norder,norient_edge,norient_face,nrdofH,shapH,gradH)
!        ...Geometry map
            call geom3D(Mdle,xip,xnod,shapH,gradH,nrdofH,x,dxdxi,dxidx,rjac,iflag)
            if (iflag.ne.0) then
                  write(*,5999) Mdle,rjac
 5999             format('Negative JacobiancMdle,rjac=',i8,2x,e12.5)
               stop
            endif
!
!        ...compute total quadrature weight
            wt123=wt1*wt2*wt3
!
!        ...compute Jacobian determinant * quadrature weight
            weighthh=wt123*rjac
!
!        ...Determine D = J^-1 * J^-T. Multiply by appropriate weight
            D(1,1)=weighthh*(dxidx(1,1)**2+dxidx(1,2)**2+dxidx(1,3)**2)
            D(1,2)=weighthh*(dxidx(1,1)*dxidx(2,1)+ &
                           dxidx(1,2)*dxidx(2,2)+dxidx(1,3)*dxidx(2,3))
            D(1,3)=weighthh*(dxidx(1,1)*dxidx(3,1)+ &
                           dxidx(1,2)*dxidx(3,2)+dxidx(1,3)*dxidx(3,3))
            D(2,1)=D(1,2)
            D(2,2)=weighthh*(dxidx(2,1)**2+dxidx(2,2)**2+dxidx(2,3)**2)
            D(2,3)=weighthh*(dxidx(2,1)*dxidx(3,1)+ &
                           dxidx(2,2)*dxidx(3,2)+dxidx(2,3)*dxidx(3,3))
            D(3,1)=D(1,3)
            D(3,2)=D(2,3)
            D(3,3)=weighthh*(dxidx(3,1)**2+dxidx(3,2)**2+dxidx(3,3)**2)
!
!        ...compute inverse Jacobian determinant * quadrature weight
            weightvv=wt123/rjac
!
!        ...Determine C = J^T * J.  Multiply by appropriate weight
            C(1,1)=weightvv*(dxdxi(1,1)**2+dxdxi(2,1)**2+dxdxi(3,1)**2)
            C(1,2)=weightvv*(dxdxi(1,1)*dxdxi(1,2)+ &
                           dxdxi(2,1)*dxdxi(2,2)+dxdxi(3,1)*dxdxi(3,2))
            C(1,3)=weightvv*(dxdxi(1,1)*dxdxi(1,3)+ &
                           dxdxi(2,1)*dxdxi(2,3)+dxdxi(3,1)*dxdxi(3,3))
            C(2,1)=C(1,2)
            C(2,2)=weightvv*(dxdxi(1,2)**2+dxdxi(2,2)**2+dxdxi(3,2)**2)
            C(2,3)=weightvv*(dxdxi(1,2)*dxdxi(1,3)+ &
                           dxdxi(2,2)*dxdxi(2,3)+dxdxi(3,2)*dxdxi(3,3))
            C(3,1)=C(1,3)
            C(3,2)=C(2,3)
            C(3,3)=weightvv*(dxdxi(1,3)**2+dxdxi(2,3)**2+dxdxi(3,3)**2)
!
!        ...put appropriate quadrature weight on Jacobian and its inverse
            dxdxi = dxdxi * weightvv
            dxidx = dxidx * wt123
!
!--------------------------------------------------------------------------
!        ...SUM FACTORIZATION LOOPS START
!--------------------------------------------------------------------------
!        ... For Gram matrix:
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
                           AUXEE_A(b,a,k3) = AUXEE_A(b,a,k3)    &
                                           + (shapH3(idxa,sa)    &
                                           * shapH3(idxb,sb)*D(a,b))
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
                              AUXCC_A(alph,beta,b,a,k3)=                   &
                                          AUXCC_A(alph,beta,b,a,k3)+       &
                                          shapH3(idxa,sa)*shapH3(idxb,sb)* &
                                          (-1)**(alph+beta)*C(idxalph,idxbeta)
                           enddo; enddo
!
!                       ...loop over components a+alph, where the curl of
!                          the TEST shape function, for the CE term of Gram matrix
                           do alph=1,2
                              idxalph=mod(a+alph-1,3)+1
                              sb=1+deltak(b,3)
                              sa=1+1-deltak(idxalph,3)
!                          ...the only nonzero is when  a+alph == b
!                             ...accumulate innermost 1D integral for CE term
                                 AUXCE_A(alph,b,a,k3)=AUXCE_A(alph,b,a,k3) &
                                      +shapH3(idxa,sa)*shapH3(idxb,sb)     &
                                      *(-1)**(alph-1)*wt123
                           enddo
!                       ...loop over components b+beta, where the curl of
!                          the TRIAL shape function, for the CE term of Gram matrix
                           do beta=1,2
                              idxbeta=mod(b+beta-1,3)+1
                              sb=1+1-deltak(idxbeta,3)
                              sa=1+deltak(a,3)
!                          ...the only nonzero is when  b+beta == a
!                             ...accumulate innermost 1D integral for EC term
                                 AUXEC_A(beta,b,a,k3)=AUXEC_A(beta,b,a,k3) &
                                      +shapH3(idxa,sa)*shapH3(idxb,sb)     &
                                      *(-1)**(beta-1)*wt123
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
!        ...loop over pz ends
         enddo
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
                              AUXEE_B(b,a,k2,k3)=AUXEE_B(b,a,k2,k3)+ &
                                                shapH2(idxa,sa)* &
                                                shapH2(idxb,sb)* &
                                                AUXEE_A(b,a,k3)
!                          ...loop over b+beta, a+alph, for curl of trial and test
                              do beta=1,2;do alph=1,2
                                 idxbeta=mod(b+beta-1,3)+1
                                 idxalph=mod(a+alph-1,3)+1
                                 sb=1+1-deltak(idxbeta,2)
                                 sa=1+1-deltak(idxalph,2)
!                             ...accumulate middle 1D integral of CC term in Gram
                                 AUXCC_B(alph,beta,b,a,k2,k3)=          &
                                       AUXCC_B(alph,beta,b,a,k2,k3)     &
                                       +shapH2(idxa,sa)*shapH2(idxb,sb) &
                                       *AUXCC_A(alph,beta,b,a,k3)
                              enddo; enddo
!                          ... accumulate for CE term
                              do alph=1,2
                                 idxalph=mod(a+alph-1,3)+1
                                 sb=1+deltak(b,2)
                                 sa=1+1-deltak(idxalph,2)
                                 if (idxalph.eq.b) then
                                    AUXCE_B(alph,b,a,k2,k3)=              &
                                          AUXCE_B(alph,b,a,k2,k3)         &
                                         +shapH2(idxa,sa)*shapH2(idxb,sb) &
                                         *AUXCE_A(alph,b,a,k3)
                                 endif
                              enddo
!                          ... accumulate for EC term
                              do beta=1,2
                                 idxbeta=mod(b+beta-1,3)+1
                                 sb=1+1-deltak(idxbeta,2)
                                 sa=1+deltak(a,2)
                                 if (idxbeta.eq.a) then
                                    AUXEC_B(beta,b,a,k2,k3)=              &
                                          AUXEC_B(beta,b,a,k2,k3)         &
                                         +shapH2(idxa,sa)*shapH2(idxb,sb) &
                                         *AUXEC_A(beta,b,a,k3)
                                 endif
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
                                       AP_Maxwell(kk)=AP_Maxwell(kk)      &
                                                     +shapH1(idxa,sa)*    &
                                                      shapH1(idxb,sb)*    &
                                                      AUXEE_B(b,a,k2,k3)* &
                                                      (abs(zb)**2 + 1.d0)
!                                   ...sum CC terms
                                       do beta=1,2; do alph=1,2
                                          idxbeta=mod(b+beta-1,3)+1
                                          idxalph=mod(a+alph-1,3)+1
                                          sb=1+1-deltak(idxbeta,1)
                                          sa=1+1-deltak(idxalph,1)
                                          AP_Maxwell(kk)=AP_Maxwell(kk)+    &
                                                         shapH1(idxa,sa)*   &
                                                         shapH1(idxb,sb)*   &
                                                         AUXCC_B(alph,beta,b,a,k2,k3)

                                       enddo; enddo

                                      kk = nk(2*m1-1,2*m2  )
!                                   ...sum CE terms
                                       do alph=1,2
                                          idxalph=mod(a+alph-1,3)+1
                                          sb=1+deltak(b,1)
                                          sa=1+1-deltak(idxalph,1)
                                          if (idxalph.eq.b) then
                                             AP_Maxwell(kk) = AP_Maxwell(kk)       &
                                                +conjg(zc)*AUXCE_B(alph,b,a,k2,k3) &
                                                *shapH1(idxa,sa)*shapH1(idxb,sb)
                                          endif
                                       enddo
!                                   ...sum EC terms
                                       do beta=1,2
                                          idxbeta=mod(b+beta-1,3)+1
                                          sb=1+1-deltak(idxbeta,1)
                                          sa=1+deltak(a,1)
                                          if (idxbeta.eq.a) then
                                             AP_Maxwell(kk) = AP_Maxwell(kk)    &
                                                   +conjg(zb)*AUXEC_B(beta,b,a,k2,k3)  &
                                                   *shapH1(idxa,sa)*shapH1(idxb,sb)
                                          endif
                                       enddo

                                       if (m1.ne.m2) then

                                          kk = nk(2*m1  ,2*m2-1)
!                                      ...sum CE terms
                                          do alph=1,2
                                             idxalph=mod(a+alph-1,3)+1
                                             sb=1+deltak(b,1)
                                             sa=1+1-deltak(idxalph,1)
                                             if (idxalph.eq.b) then
                                                AP_Maxwell(kk) = AP_Maxwell(kk)       &
                                                   +(zb)*AUXCE_B(alph,b,a,k2,k3) &
                                                   *shapH1(idxa,sa)*shapH1(idxb,sb)
                                             endif
                                          enddo
!                                      ...sum EC terms
                                          do beta=1,2
                                             idxbeta=mod(b+beta-1,3)+1
                                             sb=1+1-deltak(idxbeta,1)
                                             sa=1+deltak(a,1)
                                             if (idxbeta.eq.a) then
                                                AP_Maxwell(kk) = AP_Maxwell(kk)    &
                                                      +(zc)*AUXEC_B(beta,b,a,k2,k3)  &
                                                      *shapH1(idxa,sa)*shapH1(idxb,sb)
                                             endif
                                          enddo

                                       endif

                                       kk = nk(2*m1  ,2*m2  )
!                                   ...sum EE terms
                                       sb=1+deltak(b,1)
                                       sa=1+deltak(a,1)
                                       AP_Maxwell(kk)=AP_Maxwell(kk)      &
                                                     +shapH1(idxa,sa)*    &
                                                      shapH1(idxb,sb)*    &
                                                      AUXEE_B(b,a,k2,k3)* &
                                                      (abs(zc)**2 + 1.d0)
!                                   ...sum CC terms
                                       do beta=1,2; do alph=1,2
                                          idxbeta=mod(b+beta-1,3)+1
                                          idxalph=mod(a+alph-1,3)+1
                                          sb=1+1-deltak(idxbeta,1)
                                          sa=1+1-deltak(idxalph,1)
                                          AP_Maxwell(kk)=AP_Maxwell(kk)+    &
                                                         shapH1(idxa,sa)*   &
                                                         shapH1(idxb,sb)*   &
                                                         AUXCC_B(alph,beta,b,a,k2,k3)

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
!  ...loop over px ends
   enddo
!
!-----------------------------------------------------------------------
!
!..NORMAL element integrals for LOAD (i.e., l-Bu)
!
!-----------------------------------------------------------------------
!
!..use the enriched order to set the quadrature
   INTEGRATION = NORD_ADD
   call set_3Dint_DPG(etype,norder, nint3,xiloc,waloc)
   INTEGRATION = 0
!
!..loop over
   do l=1,nint3
      xi(1:3) = xiloc(1:3,l)
      wa = waloc(l)
!  ...determine element H1 shape functions
      call shape3H(etype,xi,norder,norient_edge,norient_face,  &
                   nrdofH,shapH,gradH)
!  ...determine element H(curl) shape functions
      call shape3E(etype,xi,norder,norient_edge,norient_face, &
                    nrdofE,shapE,curlE)
!  ...determine element L2 shape functions
      call shape3Q(etype,xi,norder, nrdofQ,shapQ)
!  ...determine discontinuous H(curl) shape functions
      call shape3EE(etype,xi,nordP, NrdofEE,shapEE,curlEE)
!  ...geometry
      call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH, &
                   x,dxdxi,dxidx,rjac,iflag)
!  ...integration weight
      weight = rjac*wa
!  ...compute the approximate solution
      zsolQ = ZERO
!
      do k=1,nrdofQ
         if(NO_PROBLEM.eq.3) then
            zsolQ(1:6)  = zsolQ(1:6)  + zdofQ(1:6,k)*shapQ(k)
         elseif(NO_PROBLEM.eq.4) then
            zsolQ(1:6)  = zsolQ(1:6)  + zdofQ(7:12,k)*shapQ(k)
         else
            write(*,*) 'error in elem_residual_UWMaxwell: NO_PROBLEM must be 3 or 4. stop.'
         endif
      enddo
      zsolQ = zsolQ/rjac
!
!  ...get the RHS
!  ...zfval (heat eqn rhs), zJ (maxwell rhs)
      call getf(Mdle,x, zfval,zJ)
!
!  ...loop through enriched H(curl) test functions
      do k1=1,NrdofEE
         E1(1:3) = shapEE(1,k1)*dxidx(1,1:3) &
                 + shapEE(2,k1)*dxidx(2,1:3) &
                 + shapEE(3,k1)*dxidx(3,1:3)
         curlE1(1:3) = dxdxi(1:3,1)*curlEE(1,k1) &
                     + dxdxi(1:3,2)*curlEE(2,k1) &
                     + dxdxi(1:3,3)*curlEE(3,k1)
         curlE1(1:3) = curlE1(1:3)/rjac
!
!     ...accumulate for the load
         k = 2*k1-1
         BLOADE(k) = BLOADE(k)+(zJ(1)*E1(1)+zJ(2)*E1(2)+zJ(3)*E1(3))*weight
!
!     ...order is modified so that A*=-A
         k = 2*k1
         BLOADE(k) = BLOADE(k) &
          - ((curlE1(1)*zsolQ(1)+curlE1(2)*zsolQ(2)+curlE1(3)*zsolQ(3)) &
             +zc*(E1(1)*zsolQ(4)+E1(2)*zsolQ(5)+E1(3)*zsolQ(6)) &
          )*weight
!
         k = 2*k1-1
         BLOADE(k) = BLOADE(k) &
          - ((curlE1(1)*zsolQ(4)+curlE1(2)*zsolQ(5)+curlE1(3)*zsolQ(6)) &
             -za*(E1(1)*zsolQ(1)+E1(2)*zsolQ(2)+E1(3)*zsolQ(3)) &
          )*weight
!
!     ...Computation of Gram Matrix
!     ...loop through enriched H(curl) trial functions
!         do k2=k1,NrdofEE
!            E2(1:3) = shapEE(1,k2)*dxidx(1,1:3) &
!                    + shapEE(2,k2)*dxidx(2,1:3) &
!                    + shapEE(3,k2)*dxidx(3,1:3)
!            curlE2(1:3) = dxdxi(1:3,1)*curlEE(1,k2) &
!                        + dxdxi(1:3,2)*curlEE(2,k2) &
!                        + dxdxi(1:3,3)*curlEE(3,k2)
!            curlE2(1:3) = curlE2(1:3)/rjac
!!
!!        ...auxiliary quantities
!            CC = curlE1(1)*curlE2(1) + curlE1(2)*curlE2(2)+ curlE1(3)*curlE2(3)
!            EE = E1(1)*E2(1) + E1(2)*E2(2) + E1(3)*E2(3)
!            CE = curlE1(1)*E2(1) + curlE1(2)*E2(2) + curlE1(3)*E2(3)
!            EC = E1(1)*curlE2(1) + E1(2)*curlE2(2) + E1(3)*curlE2(3)
!!
!!        ...accumulate for the test Gram matrix: adjoint graph
!            k = nk(2*k1-1,2*k2-1)
!            AP_Maxwell(k) = AP_Maxwell(k) &
!                 + (CC+ (abs(zb)**2 + 1.d0)*EE)*weight
!            k = nk(2*k1-1,2*k2  )
!            AP_Maxwell(k) = AP_Maxwell(k) &
!                  + (-za*EC - zc*CE)*weight
!            if (k1.ne.k2) then
!               k = nk(2*k1  ,2*k2-1)
!               AP_Maxwell(k) = AP_Maxwell(k) &
!                  + (zb*CE + zc*EC)*weight
!            endif
!            k = nk(2*k1  ,2*k2  )
!            AP_Maxwell(k) = AP_Maxwell(k) &
!               + (CC+ ((abs(zc))**2 + 1.d0)*EE)*weight
!         enddo
      enddo
   enddo
!
!-----------------------------------------------------------------------
!
!              B O U N D A R Y      I N T E G R A L S
!
!-----------------------------------------------------------------------
!
!
!..auxiliary constant for flux computation
   if(NO_PROBLEM.eq.3) then
      i = 0
   elseif(NO_PROBLEM.eq.4) then
      i = 2
   else
      write(*,*) 'elem_residual_UWMaxwell: NO_PROBLEM must be 3 or 4. stop.'
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
                      nrdofH,shapH,gradH)
!
!     ...determine element H(curl) shape functions (for fluxes)
         call shape3E(etype,xi,norder,norient_edge,norient_face, &
                      nrdofE,shapE,curlE)
!
!     ...geometry
         call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign, &
                      x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
         weight = bjac*wtloc(l)
!
!     ...compute approximate fluxes at the point
         zsolExi = ZERO
!
         do ivar=1,2
            do k=1,nrdofE
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
            BLOADE(k) = BLOADE(k) &
             - (E1(1)*zflux(1,1)+E1(2)*zflux(2,1)+E1(3)*zflux(3,1)) &
             *weight
!        ...modifying for A*=-A
            k=2*k1-1
            if(ibc(ifc,2).eq.9) then
               BLOADE(k) = BLOADE(k) &
                - (E1(1)*zflux2(1,1)+E1(2)*zflux2(2,1)+E1(3)*zflux2(3,1)) &
                *weight
            else
               BLOADE(k) = BLOADE(k) &
                - (E1(1)*zflux(1,2)+E1(2)*zflux(2,2)+E1(3)*zflux(3,2)) &
                *weight
            endif
         enddo
      enddo
      if (iprint.ne.0) pause
   enddo
!
   if (iprint.gt.0) then
      write(*,7015) BLOADE(1:2*NrdofEE)
 7015 format('elem_residual: FINAL BLOADE = ',10(/,6(2e12.5,2x)))
      call pause
   endif
!
!-----------------------------------------------------------------------
!
!..factorize the test stiffness matrix
   call ZPPTRF('U', nrTest, AP_Maxwell, info)
   if (info.ne.0) then
      write(*,*) 'elem_residual: info = ',info
      stop
   endif
!
!..save copies of the RHS to compute later the residual
   BLOADEc = BLOADE
!
!..compute the product of inverted test Gram matrix with RHS,
!..BLOADE is overwritten with the solution
   call ZPPTRS('U', nrTest, 1, AP_Maxwell, BLOADE, nrTest, info)
   if (info.ne.0) then
      write(*,*) 'elem_residual: info = ',info
      stop
   endif
!
!..compute the residual
   zresid = ZERO
   do k=1,nrTest
      zresid = zresid + BLOADEc(k)*conjg(BLOADE(k))
   enddo
   Resid = real(zresid,8)
   Nref_flag = 001
!
   if (iprint.eq.1) then
      write(*,7010) Mdle, Resid
 7010 format('elem_residual: Mdle, Resid = ',i5,3x,e12.5)
      pause
   endif
!
end subroutine elem_residual_UWMaxwell
