!--------------------------------------------------------------------
!
!     routine name      - elem
!
!--------------------------------------------------------------------
!
!     latest revision:  - July 18
!
!     purpose:          - routine returns unconstrained (ordinary)
!                         stiffness matrix and load vector
!                         for either:
!                         1) single step of transient heat problem
!                         2) multiple steps of transient heat problem
!                         3) steady state Maxwell problem
!
!     arguments:
!        in:
!              Mdle     - an element middle node number, identified
!                         with the element
!
!        out:   
!              Itest    - index for assembly
!              Itrial   - index for assembly
!
!
!-----------------------------------------------------------------------
!
#include "implicit_none.h"
!
subroutine elem(Mdle, Itest,Itrial)
!
   use assembly
   use CommonParam
   use control
   use parametersDPG
   use data_structure3D
!
   implicit none
!
!..declare input/output variables
   integer,                    intent(in)  :: Mdle
   integer,dimension(NR_PHYSA),intent(out) :: Itest,Itrial
!..number of test and trial degrees of freedom
   integer :: nrdofH ,nrdofE ,nrdofV ,nrdofQ
   integer :: nrdofHH,nrdofEE,nrdofVV,nrdofQQ
!..
   integer :: norder(19),norderP(19),nordP
!..
   character(len=4) :: etype
!
!..fld_flag refers to either pump (0) or signal (1) field
   integer :: fld_flag
!
   Itest (1:NR_PHYSA) = 0
   Itrial(1:NR_PHYSA) = 0
!
   select case(NODES(Mdle)%case)
!  ...node supports all physical attributes
!  ...6 physical attributes: case = 2^6-1 = 63
      case(63)
!     ...NO_PROBLEM: (1) - single step of transient heat equation
!        ............(2) - multiple steps of transient heat equation
!        ............(3) - time harmonic Maxwell - for signal EH-field
!        ............(4) - time harmonic Maxwell - for pump EH-field
         select case(NO_PROBLEM)
!        ...Heat cases
            case(1,2)
               Itest(1)=1; Itrial(1)=1
               Itest(4)=1; Itrial(4)=1
               call elem_dpgHeat(Mdle, &
                 BLOC(1)%nrow,BLOC(4)%nrow, &
                 BLOC(1)%array,ALOC(1,1)%array,ALOC(1,4)%array, &
                 BLOC(4)%array,ALOC(4,1)%array,ALOC(4,4)%array)
!        ...Maxwell cases - signal EH-field
            case(3)
               Itest(2)=1; Itrial(2)=1
               Itest(5)=1; Itrial(5)=1
               fld_flag = 1;
               etype = NODES(Mdle)%type
!           ...determine order of approximation
               call find_order(Mdle, norder)
!           ...set the enriched order of appoximation
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
!           ...note: compute_enriched_order works only for hexa currently
               call compute_enriched_order(nordP, norderP)
!           ...compute nrdof for trial
               call celndof(etype,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
!           ...compute nrdof for test
               call celndof(etype,norderP, nrdofHH,nrdofEE,nrdofVV,nrdofQQ)
!
               call elem_uwMaxwell_fi_TRANS(Mdle,fld_flag, &
                 nrdofEE,nrdofH,nrdofE,nrdofQ, &
                 BLOC(2)%nrow,BLOC(5)%nrow, &
                 BLOC(2)%array,ALOC(2,2)%array,ALOC(2,5)%array, &
                 BLOC(5)%array,ALOC(5,2)%array,ALOC(5,5)%array)
!        ...Maxwell cases - pump EH-field
            case(4)
               Itest(3)=1; Itrial(3)=1
               Itest(6)=1; Itrial(6)=1
               fld_flag = 0;
               etype = NODES(Mdle)%type
!           ...determine order of approximation
               call find_order(Mdle, norder)
!           ...set the enriched order of appoximation
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
!           ...note: compute_enriched_order works only for hexa currently
               call compute_enriched_order(nordP, norderP)
!           ...compute nrdof for trial
               call celndof(etype,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
!           ...compute nrdof for test
               call celndof(etype,norderP, nrdofHH,nrdofEE,nrdofVV,nrdofQQ)
!
               call elem_uwMaxwell_fi_TRANS(Mdle,fld_flag, &
                 nrdofEE,nrdofH,nrdofE,nrdofQ, &
                 BLOC(3)%nrow,BLOC(6)%nrow, &
                 BLOC(3)%array,ALOC(3,3)%array,ALOC(3,6)%array, &
                 BLOC(6)%array,ALOC(6,3)%array,ALOC(6,6)%array)
!
!     ...end select for select case(NO_PROBLEM)
         end select
!
      case default
         write(*,*) 'elem: Mdle,NODES(Mdle)%case = ', &
           Mdle,NODES(Mdle)%case, '. stop.'
         call logic_error(ERR_INVALID_VALUE, __FILE__,__LINE__)
         stop
!..end select for select case(NODES(Mdle)%case)
   end select
!
!
end subroutine elem

!
!--------------------------------------------------------------------
!
!    routine name      - elem_dpgHeat
!
!--------------------------------------------------------------------
!
!    latest revision:  - April 17
!
!    purpose:          - routine returns unconstrained (ordinary)
!                        stiffness matrix and load vector
!                        for the primal H1 formulation for a single
!                        step of transient heat equation
!
!    arguments:
!
!    in:
!            Mdle      - an element middle node number, identified
!                        with the element
!            MdH       - column length of ZalocHH,ZalocHV
!            MdV       - column length of ZalocVH,ZalocVV
!    out:
!            ZblocH,ZblocV - load vectors
!            ZalocHH,ZalocHV,ZalocVH,ZalocVV - stiffness matrices
!
!---------------------------------------------------------------------
!
subroutine elem_dpgHeat(Mdle,MdH,MdV, &
   ZblocH,ZalocHH,ZalocHV,ZblocV,ZalocVH,ZalocVV)
!
   use control
   use parametersDPG
   use element_data
   use data_structure3D
   use LaserParam
   use CommonParam
!
   implicit none
!..declare input/output variables
   integer,                     intent(in)  :: Mdle
   integer,                     intent(in)  :: MdH
   integer,                     intent(in)  :: MdV
   VTYPE, dimension(MdH),       intent(out) :: ZblocH
   VTYPE, dimension(MdH,MdH),   intent(out) :: ZalocHH
   VTYPE, dimension(MdH,MdV),   intent(out) :: ZalocHV
   VTYPE, dimension(MdV),       intent(out) :: ZblocV
   VTYPE, dimension(MdV,MdH),   intent(out) :: ZalocVH
   VTYPE, dimension(MdV,MdV),   intent(out) :: ZalocVV
!
!..declare edge/face type variables
   character(len=4) :: etype,ftype
!
!..declare element order, orientation for edges and faces
   integer, dimension(19)    :: norder
   integer, dimension(12)    :: norient_edge
   integer, dimension(6)     :: norient_face
!..face order
   integer, dimension(5)     :: norderf
!
!..geometry dof (work space for nodcor)
   real*8, dimension(3,MAXbrickH) :: xnod
!
!..solution dof (work space for solelm)
   VTYPE, dimension(MAXEQNH,MAXbrickH) :: zdofH
   VTYPE, dimension(MAXEQNE,MAXbrickE) :: zdofE
   VTYPE, dimension(MAXEQNV,MAXbrickV) :: zdofV
   VTYPE, dimension(MAXEQNQ,MAXbrickQ) :: zdofQ
!..approximate solution -- using soleval
   integer :: nflag
   VTYPE, dimension(  MAXEQNH  ) ::  zsolH_soleval
   VTYPE, dimension(  MAXEQNH,3) :: zdsolH_soleval
   VTYPE, dimension(3,MAXEQNE  ) ::  zsolE_soleval
   VTYPE, dimension(3,MAXEQNE  ) :: zcurlE_soleval
   VTYPE, dimension(3,MAXEQNV  ) ::  zsolV_soleval
   VTYPE, dimension(  MAXEQNV  ) ::  zdivV_soleval
   VTYPE, dimension(  MAXEQNQ  ) ::  zsolQ_soleval
!
!..variables for geometry
   real*8, dimension(3)      :: xi,x,rn
   real*8, dimension(3,2)    :: dxidt,dxdt,rt
   real*8, dimension(3,3)    :: dxdxi,dxidx
   real*8, dimension(2)      :: t
!
!..H1 shape functions
   real*8, dimension(MAXbrickH)    :: shapH
   real*8, dimension(3,MAXbrickH)  :: gradH
!
!..H(curl) shape functions
   real*8, dimension(3,MAXbrickE)  :: shapE
   real*8, dimension(3,MAXbrickE)  :: curlE
!
!..H(div) shape functions
   real*8, dimension(3,MAXbrickV)  :: shapV
   real*8, dimension(MAXbrickV)    :: divV
!
!..L2 shape functions
   real*8, dimension(MAXbrickQ)    :: shapQ
!
!..Enriched H1 shape functions
   real*8 , dimension(MAXbrickHH)    :: shapHH
   real*8 , dimension(3,MAXbrickHH)  :: gradHH
!
!..nrdof for various spaces
   integer  :: nrdofH,nrdofE,nrdofV,nrdofQ,nrdofHH
!
!..space for DPG Computations (Gram Matrix, Stiffness etc.)
   integer, parameter   :: MAXtestH = MAXbrickHH
!
!..stiffness matrix for the local Riesz H1 matrix in Lapack format
   VTYPE, dimension(MAXtestH*(MAXtestH+1)/2) :: AP_Heat
!
!..load vector for the enriched space
   VTYPE, dimension(MAXtestH) :: BLOADH
!
!..stiffness matrices for the enriched test space
   VTYPE, dimension(MAXtestH,MAXbrickH) :: STIFFHH
   VTYPE, dimension(MAXtestH,MAXbrickV) :: STIFFHV
!
!..STIFF_ALL for alternative computation of stiffness
   VTYPE, dimension(MAXtestH,MAXbrickH+MAXbrickV+1) :: STIFF_ALLH
   VTYPE, allocatable:: SQUARE_STIFF(:,:)
!
!..TODO
   VTYPE, allocatable :: AP_eig(:)
   VTYPE, allocatable :: DIAG_E(:)
   VTYPE, allocatable :: DIAG_H(:)
!
!..dummy for elem_residual
   VTYPE, dimension(MAXtestH*(MAXtestH+1)/2) :: AP
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
!..derivatives wrt physical coordinates, flux
   real*8, dimension(3) :: dv1,dv2,vec
!
!..for debug printing
   VTYPE, dimension(10) :: zaux

!..Maxwell load
   VTYPE, dimension(3) :: zJ
!
   integer,save :: ivis=0
   character    :: uplo, trans,diag
!
!..number of vertices,edge,faces per element type
   integer :: nrv, nre, nrf
!..for Gram matrix
   integer :: nk
!..various variables for the problem
   real*8  :: h_elem,rjac,weight,wa,v2n,v1,v2,bjac
   integer :: i1,i2,j1,j2,k1,k2,kH,kk,i,j,nrTest,nint,iflag,kE,k,iprint
   integer :: N,nRHS,nordP,l,nsign,ifc,ndom,info
   VTYPE   :: zfval,therm_Load
   nk(k1,k2) = (k2-1)*k2/2+k1

!
!------------------INITIALIZE THE ELEMENT ORDER, ORIENTATION ETC.------
!
   ivis=1
   if (ivis.eq.0) then
      write(*,*) 'elem_dpgHeat: NO_PROBLEM    = ',NO_PROBLEM
      write(*,*) 'elem_dpgHeat: INNER_PRODUCT = ',INNER_PRODUCT
      call pause
      ivis=1
   endif
!
!..element type
   etype = NODES(Mdle)%type
   nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
!
!..determine order of approximation
   call find_order(Mdle, norder)
!
!..set the enriched order of appoximation depending on element type
   select case(etype)
      case('mdlb')        ; nordP = NODES(Mdle)%order + NORD_ADD*111
      case('mdln','mdld') ; nordP = NODES(Mdle)%order + NORD_ADD*1
      case('mdlp')        ; nordP = NODES(Mdle)%order + NORD_ADD*11
      case default
         write(*,*) 'invalid etype. stop.'
         stop
   end select
!
!..determine edge and face orientations
   call find_orient(Mdle, norient_edge,norient_face)
!
!..determine nodes coordinates
   call nodcor(Mdle, xnod)
!
!..determine dof for current solution if running multiple
!..steps of heat equation
   select case(NO_PROBLEM)
      case(2)
         call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
      case default
         write(*,*) 'invalid NO_PROBLEM param. stop.'
         stop
   end select
!
!..determine element size and scale correction
   !!! h_elem = min(abs(xnod(1,2)-xnod(1,1)),abs(xnod(2,4)-xnod(2,1)), &
   !!!              abs(xnod(3,5)-xnod(3,1)))
!
!
!-----------
!-----------  CHECK OUT INNER PRODUCT SCALING FOR HEAT STEPPING
!-----------
!     select case(INNER_PRODUCT)
!     case(2)
!        scale = min(1.d0,EPS_DIFF**3/h_elem**2)
!     case(3)
!        scale =  min(1.d0,EPS_DIFF**2/h_elem**2)
!     end select
!     scale = 1.d0
!-----------
!-----------  CHECK OUT INNER PRODUCT SCALING FOR HEAT STEPPING
!-----------
!
!..get the element boundary conditions flags
   call find_bc(Mdle, ibc)
!
   iprint = 0
   if (iprint.ge.1) then
      write(*,7001) Mdle
7001  format('elem_dpgHeat: BCFLAGS FOR Mdle = ',i5)
      do i=1,NR_PHYSA
         write(*,7002) PHYSA(i), ibc(1:nrf,i)
7002     format('          ATTRIBUTE = ',a6,' FLAGS = ',6i2)
      enddo
   endif
!
!..clear space for stiffness matrix and rhsv:
   ZblocH = ZERO; ZblocV = ZERO
   ZalocHH = ZERO; ZalocHV = ZERO; ZalocVH = ZERO; ZalocVV = ZERO
!
!..extended load vector and extended stiffness matrices
   BLOADH=ZERO ; STIFFHH=ZERO ; STIFFHV=ZERO; STIFF_ALLH=ZERO
!
!..Gram matrix
   AP_Heat=ZERO
!
!
!-----------------------------------------------------------------------
!    E L E M E N T    I N T E G R A L S                               |
!-----------------------------------------------------------------------
!
!..use the enriched order to set the quadrature
   INTEGRATION = NORD_ADD
   call set_3Dint_DPG(etype,norder, nint,xiloc,waloc)
   INTEGRATION = 0
!
!..loop over integration points
   do l=1,nint
      xi(1:3)=xiloc(1:3,l) ; wa=waloc(l)
!
!  ...H1 shape functions
      call shape3H(etype,xi,norder,norient_edge,norient_face, &
        nrdofH,shapH,gradH)
!
!  ...discontinuous H1 shape functions
      call shape3HH(etype,xi,nordP, nrdofHH,shapHH,gradHH)
      !!! write(*,*) 'nrdofHH is: ', nrdofHH
      !!! call pause
!
!  ...geometry computations
      call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH, &
        x,dxdxi,dxidx,rjac,iflag)
!  ...compute current solution using soleval
      nflag = 1
      call soleval(Mdle,xi,norient_edge,norient_face,norder,xnod, &
        zdofH,zdofE,zdofV,zdofQ,nflag,x,dxdxi,zsolH_soleval,zdsolH_soleval, &
        zsolE_soleval,zcurlE_soleval,zsolV_soleval,zdivV_soleval,zsolQ_soleval)
!
!..................................................................
!
!  ...FOR THE LASER_MODE = 1
      if(LASER_MODE.eq.1) then
!     ...compute the H1 solution load at the point
         call get_thermLoad(zsolQ_soleval, therm_Load)
      endif
!
!  ...integration weight
      weight = rjac*wa
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
!        -- EXTENDED LOAD VECTOR --
!
         select case(NO_PROBLEM)
            case(1)
               BLOADH(k1) = BLOADH(k1) + zfval*v1*weight
            case(2)
               BLOADH(k1) = BLOADH(k1) + (therm_Load*DELTAT+zsolH_soleval(1))*v1*weight
!        ...do not solve for heat equation in PML
            if(COPUMP.eq.1) then
               if(x(3).ge.PML_REGION) then
                  BLOADH(k1) = ZERO
               endif
            elseif(COPUMP.eq.0) then
               write(*,*) 'heat equation not enabled for counter pump configuration. stop.'
               stop
            endif
         end select
!
!     ...2nd loop through enriched H1 trial functions
!     ...for Gram matrix
         do k2=k1,nrdofHH
!
            v2 = shapHH(k2)
            dv2(1:3) = gradHH(1,k2)*dxidx(1,1:3) &
                     + gradHH(2,k2)*dxidx(2,1:3) &
                     + gradHH(3,k2)*dxidx(3,1:3)
!
! ----------- GRAM MATRIX --------------
!
!        ...determine index in triangular format
            k = nk(k1,k2)
!
!        ...problem-dependent inner product
            select case(INNER_PRODUCT)
               case(1)
                  AP_Heat(k) = AP_Heat(k) &
                             + (dv1(1)*dv2(1) + dv1(2)*dv2(2) + &
                                dv1(3)*dv2(3) +  v1*v2) * weight
               case default
                  write(*,*) 'unexpected INNER_PRODUCT param. stop.'
                  stop
            end select
!
!     ...enddo 2nd loop through enriched H1 trial functions
         enddo
!
!     ...loop through H1 trial functions
!     ...for stiffness matrix
         do k2=1,nrdofH
            v2 = shapH(k2)
            dv2(1:3) = gradH(1,k2)*dxidx(1,1:3) &
                  + gradH(2,k2)*dxidx(2,1:3) &
                  + gradH(3,k2)*dxidx(3,1:3)
!
!          -- EXTENDED HH STIFFNESS MATRIX --
            STIFFHH(k1,k2) = STIFFHH(k1,k2) &
                           + (KAPPA*DELTAT * (&
                              dv1(1)*dv2(1)+dv1(2)*dv2(2)+dv1(3)*dv2(3) &
                             ) + v1*v2)*weight
!     ...end loop through H1 trial functions
         enddo
!  ...end 1st loop through enriched H1 test functions
      enddo
!..end loop through integration points
   enddo
!
!..printing Gram matrix
   iprint = 0
   if (iprint.eq.1) then
      write(*,*) 'elem_dpgHeat: GRAM MATRIX = '
      do i=1,10
         do j=1,i-1
            zaux(j) = AP_Heat(nk(j,i))
         enddo
         do j=i,10
            zaux(j) = AP_Heat(nk(i,j))
         enddo
         write(*,7011) zaux
      enddo
      call pause
   endif
!
!
!-----------------------------------------------------------------------
!    B O U N D A R Y    I N T E G R A L S                             |
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
      call set_2Dint(ftype,norderf, nint,tloc,wtloc)
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
!     ...loop through enriched H1 test functions
         do k1=1,nrdofHH
            v1 = shapHH(k1)
!
!        ...loop through H(div) trial functions
            do k2=1,nrdofV
!
!           ...normal component (Piola transformation at work!)
               vec(1:3) = dxdxi(1:3,1)*shapV(1,k2) &
                        + dxdxi(1:3,2)*shapV(2,k2) &
                        + dxdxi(1:3,3)*shapV(3,k2)
               vec(1:3) = vec(1:3)/rjac
               v2n = vec(1)*rn(1)+vec(2)*rn(2)+vec(3)*rn(3)
!
!            -- EXTENDED HV STIFFNESS MATRIX --
!
               STIFFHV(k1,k2) = STIFFHV(k1,k2) &
                          - KAPPA*DELTAT*v1*v2n*weight
!        ...end loop through H(div) trial functions
            enddo
!     ...end loop through enriched H1 test functions
         enddo
!     ...end loop through integration points
      enddo
!..end loop through element faces
   enddo
!
!-----------------------------------------------------------------------
!
!..factorize the Gram matrix
   uplo = 'U'
   call ZPPTRF(uplo, nrdofHH, AP_Heat, info)
   if (info.ne.0) then
      write(*,*) 'elem_dpgHeat: info = ',info ; stop
   endif

! ....................................................................
!  construction of DPG system
! ....................................................................
!
!..total trial dof for the element
   nrTest = nrdofHH
   call celndof(etype,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
!
   i1 = MAXtestH ; j1 = nrdofH ; j2 = nrdofV
!
   STIFF_ALLH(1:i1,1:j1) = STIFFHH(1:i1,1:j1)
   STIFF_ALLH(1:i1,j1+1:j1+j2) = STIFFHV(1:i1,1:j2)
   STIFF_ALLH(1:i1,j1+j2+1) = BLOADH(1:i1)
!
   uplo  = 'U' ; trans = 'C' ;   diag  = 'N'
   N     = nrTest
   nRHS  = nrdofH + nrdofV + 1

   !!! write(*,*) 'from elem_dpgHeat before ZHERK: nRHS, nrTest, MAXtestH = ',nRHS, nrTest, MAXtestH
   !!! write(*,*) 'from elem_dpgHeat size of  STIFF_ALLH= ', MAXtestH,MAXbrickH+MAXbrickV+1
   !!! call pause
   allocate(SQUARE_STIFF(NRHS,NRHS))
!
   call ZTPTRS(uplo,trans,diag,N,NRHS,AP_heat,STIFF_ALLH, &
              MAXtestH,info)
!
   !!! call ZHERK(uplo,trans,NRHS,nrTest,ZONE,STIFF_ALLH, &
   !!!             MAXtestH,ZERO, &
   !!!             STIFF_ALLH(1:NRHS,1:NRHS),NRHS)

  call ZHERK(uplo,trans,NRHS,nrTest,ZONE,STIFF_ALLH, &
              MAXtestH,ZERO, &
              SQUARE_STIFF(1:NRHS,1:NRHS),NRHS)
!..ZHERK for complex case
   do i=1,NRHS-1
      !!!STIFF_ALLH(i+1:NRHS,i) = conjg(STIFF_ALLH(i,i+1:NRHS))
      SQUARE_STIFF(i+1:NRHS,i) = conjg(SQUARE_STIFF(i,i+1:NRHS))
   enddo
! !
!   ZblocH(1:j1) = STIFF_ALLH(1:j1,j1+j2+1)
!   ZblocV(1:j2) = STIFF_ALLH(j1+1:j1+j2,j1+j2+1)
! !
!   ZalocHH(1:j1,1:j1) = STIFF_ALLH(1:j1,1:j1)
!   ZalocHV(1:j1,1:j2) = STIFF_ALLH(1:j1,j1+1:j1+j2)
! !
!   ZalocVH(1:j2,1:j1) = STIFF_ALLH(j1+1:j1+j2,1:j1)
!   ZalocVV(1:j2,1:j2) = STIFF_ALLH(j1+1:j1+j2,j1+1:j1+j2)

   ZblocH(1:j1) = SQUARE_STIFF(1:j1,j1+j2+1)
   ZblocV(1:j2) = SQUARE_STIFF(j1+1:j1+j2,j1+j2+1)
!
   ZalocHH(1:j1,1:j1) = SQUARE_STIFF(1:j1,1:j1)
   ZalocHV(1:j1,1:j2) = SQUARE_STIFF(1:j1,j1+1:j1+j2)
!
   ZalocVH(1:j2,1:j1) = SQUARE_STIFF(j1+1:j1+j2,1:j1)
   ZalocVV(1:j2,1:j2) = SQUARE_STIFF(j1+1:j1+j2,j1+1:j1+j2)
   deallocate(SQUARE_STIFF)
!
!
#if C_MODE
 7011   format(6(2e10.3,2x))
#else
 7011   format(10e12.5)
#endif
!
end subroutine elem_dpgHeat
!
!-------------------------------------------------------------------
!
!    routine name      - elem_uwMaxwell_fi_TRANS
!
!--------------------------------------------------------------------
!
!    latest revision:  - April 18
!
!    purpose:          - routine returns unconstrained (ordinary)
!                        stiffness matrix and load vector
!                        for the UW DPG formulation for Maxwell
!                        equations
!                      - Uses sum factorization for fast integration
!                        with swapped loops
!
!    arguments:
!
!    in:
!            Mdle      - an element middle node number, identified
!                        with the element
!            NrdofEE
!            NrdofH
!            NrdofE
!            NrdofQ
!            MdE       - column length of ZalocEE,ZalocEQ
!            MdQ       - column length of ZalocQE,ZalocQQ
!    out:
!            ZblocE,ZblocQ - load vectors
!            ZalocEE,ZalocEQ,ZalocQE,ZalocQQ - stiffness matrices
!
!---------------------------------------------------------------------
!
subroutine elem_uwMaxwell_fi_TRANS(Mdle,fld_flag, &
                              NrdofEE,NrdofH,NrdofE,NrdofQ, &
                              MdE,MdQ, &
                              ZblocE,ZalocEE,ZalocEQ, &
                              ZblocQ,ZalocQE,ZalocQQ)
!..modules used
   use control
   use parametersDPG
   use data_structure3D
   use LaserParam
   use CommonParam
!..no implicit statements
   implicit none
!..declare input/output variables
   integer,                   intent(in)  :: Mdle
   integer,                   intent(in)  :: fld_flag
   integer,                   intent(in)  :: NrdofEE
   integer,                   intent(in)  :: NrdofH
   integer,                   intent(in)  :: NrdofE
   integer,                   intent(in)  :: NrdofQ
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
!..set 111 for middle node (last entry)
   integer, dimension(19)    :: norderi
!
!..face order
   integer, dimension(5)     :: norderf
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
!..L2 shape functions
   real*8, dimension(MAXbrickQ)    :: shapQ
!
!..Enriched H1 shape functions
   real*8 , dimension(3,MAXbrickEE) :: shapEE
   real*8 , dimension(3,MAXbrickEE) :: curlEE
!
!..nrdof for interface only (without bubbles)
   integer :: nrdofEi, nrdofEEi
!
!..H(curl) bubble index
   integer, allocatable :: idxEE(:)
!
!..element mdle node dof
   integer :: ndofHHmdl,ndofEEmdl,ndofVVmdl,ndofQQmdl
!
!..element degrees of freedom
   integer :: nrTest
   integer :: nrTrial
!
!..stiffness matrix for the local Riesz H1 matrix in Lapack format
   VTYPE, dimension(2*NrdofEE*(2*NrdofEE+1)/2) :: AP_Maxwell
!
!..load vector for the enriched space
   VTYPE, dimension(2*NrdofEE) :: BLOADE
!
!..STIFF_ALL for alternative computation of stiffness
   VTYPE, dimension(2*NrdofEE,2*NrdofE+6*NrdofQ+1) :: STIFF_ALL
!
!..Matrices for transpose filling (swapped loops)
!..stiffness matrices (transposed) for the enriched test space
   VTYPE, dimension(6*NrdofQ,2*NrdofEE) :: STIFFEQ_T
   VTYPE, dimension(2*NrdofE,2*NrdofEE) :: STIFFEE_T
!
!..aux matrix for stiffness and load computation
   complex*16 :: zaloc(2*NrdofE+6*NrdofQ+1,2*NrdofE+6*NrdofQ+1)
!
!..TODO
   VTYPE, allocatable :: AP_eig(:)
   VTYPE, allocatable :: DIAG_E(:)
   VTYPE, allocatable :: DIAG_H(:)
!
!..dummy for elem_residual
   VTYPE, dimension(2*NrdofEE*(2*NrdofEE+1)/2) :: AP
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
!..for debug printing
   VTYPE, dimension(21) :: zaux
!
!..Maxwell load and auxiliary variables
   VTYPE, dimension(3)  :: zJ,zImp
   real*8, dimension(3) :: E1,curlE1,E2,curlE2,rntimesE,rn2timesE
!
   integer,save :: ivis=0
   character    :: uplo, trans,diag
!
!..number of vertices,edge,faces per element type
   integer :: nrv, nre, nrf
!
!..for Gram matrix
   integer :: nk
!
!..various variables for the problem
   real*8  :: h_elem,rjac,weight,wa,v2n,CC,EE,CE,E,EC,q,h,omeg,alpha_scale
   real*8  :: bjac
   integer :: i1,i2,j1,j2,k1,k2,kH,kk,i,ik,j,k,l,nint,kE
   integer :: iflag,iprint,itime,iverb
   integer :: nordP,nsign,ifc,ndom,info,icomp,nrdof_eig,idec
   VTYPE   :: zfval,za,zb,zc
!
!..for polarizations function
   VTYPE  :: bg_pol,raman_pol,therm_pol,rndotE,alpha
   integer :: dom_flag
!
!..OMEGA_RATIO_SIGNAL or OMEGA_RATIO_PUMP
   real*8  :: OMEGA_RATIO_FLD
!
!..for PML
   VTYPE :: zbeta,zdbeta,zd2beta,detJstretch
   VTYPE, dimension(3,3) :: Jstretch,INVJStretch,JJStretch
!
!..for lapack eigensolve
   complex*16, allocatable :: Z(:,:), WORK(:)
   real*8, allocatable     :: W(:),   RWORK(:)
   integer, allocatable    :: IWORK(:)
!
!..added to use fast integration
   VTYPE, allocatable :: AUXEE_A(:,:,:)
   VTYPE, allocatable :: AUXEE_B(:,:,:,:)
   VTYPE, allocatable :: AUXCC_A(:,:,:,:,:)
   VTYPE, allocatable :: AUXCC_B(:,:,:,:,:,:)
   VTYPE, allocatable :: AUXEC_A(:,:,:,:)
   VTYPE, allocatable :: AUXCE_A(:,:,:,:)
   VTYPE, allocatable :: AUXEC_B(:,:,:,:,:)
   VTYPE, allocatable :: AUXCE_B(:,:,:,:,:)
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
   real*8 :: xi1,xi2,xi3,wt1,wt2,wt3,clock1,clock2
   real*8 :: wt123,weighthh,weightvv
   real*8, dimension(MAXPP+1) :: xilocx,xilocy,xilocz
   real*8, dimension(MAXPP+1) :: wlocx,wlocy,wlocz
   real*8, dimension(3,MAXNINT3ADD) :: wloc3
   real*8, dimension(3) :: xip,dHdx,dHHdx
   real*8, dimension(3,3) :: D,C
   real*8, dimension(MAXPP+1,2) :: shapH1,shapH2,shapH3
   real*8, dimension(MAXPP+1,MAXPP+1) :: sH2p,sH3p,dsH2p,dsH3p
   integer, dimension(3,3) :: deltak
!
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
   if (iprint.eq.1) then
      write(*,*) 'elem: Mdle = ', Mdle
   endif
!..element type
   etype = NODES(Mdle)%type
   nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
!..compute element degrees of freedom
   nrTest = 2*NrdofEE
   nrTrial = 2*NrdofE + 6*NrdofQ
!..determine order of approximation
   call find_order(Mdle, norder)
   norderi(1:nre+nrf) = norder(1:nre+nrf)
!..set the enriched order of appoximation
   select case(etype)
      case('mdlb')
         nordP = NODES(Mdle)%order+NORD_ADD*111
         norderi(nre+nrf+1) = 111
      case('mdln','mdld')
         nordP = NODES(Mdle)%order+NORD_ADD
         norderi(nre+nrf+1) = 1
      case('mdlp')
         nordP = NODES(Mdle)%order+NORD_ADD*11
         norderi(nre+nrf+1) = 11
      case default
         write(*,*) 'invalid etype param. stop.'
         stop
   end select
!..determine edge and face orientations
   call find_orient(Mdle, norient_edge,norient_face)
!..determine nodes coordinates
   call nodcor(Mdle, xnod)
!..get the element boundary conditions flags
   call find_bc(Mdle, ibc)
   iprint = 0
   if (iprint.ge.1) then
      write(*,7001) Mdle
7001  format('elem: BCFLAGS FOR Mdle = ',i5)
      do i=1,NR_PHYSA
         write(*,7002) PHYSA(i), ibc(1:nrf,i)
7002     format('     ATTRIBUTE = ',a6,' FLAGS = ',6i2)
      enddo
   endif
!
! ! !..............FOR LONG WAVEGUIDE: RE-USING STIFFNESS MATRICES
! ...check if one needs to recompute element matrices
  ! idec=0
  ! call copy_elmmat_UW(MdE,MdQ,xnod,ibc, idec, &
  !                            ZalocEE,ZalocEQ, &
  !                            ZalocQE,ZalocQQ,ZblocE,ZblocQ)
  ! if (idec.eq.1) return
!..............FOR LONG WAVEGUIDE: RE-USING STIFFNESS MATRICES
!
!..clear space for stiffness matrix and rhsv:
   ZblocE  = ZERO; ZblocQ  = ZERO
   ZalocEE = ZERO; ZalocEQ = ZERO
   ZalocQE = ZERO; ZalocQQ = ZERO
!..clear space for auxiliary matrices
   BLOADE     = ZERO
   AP_Maxwell = ZERO
   STIFF_ALL  = ZERO
   zaloc      = ZERO
!..clear space for transpose matrices
   STIFFEE_T = ZERO
   STIFFEQ_T = ZERO
!
!..set fld_flag to OMEGA_RATIO_SIGNAL or OMEGA_RATIO_PUMP
   select case(fld_flag)
      case(0)
         OMEGA_RATIO_FLD = OMEGA_RATIO_PUMP
      case(1)
         OMEGA_RATIO_FLD = OMEGA_RATIO_SIGNAL ! 1.0d0
      case default
         write(*,*) 'invalid fld_flag param. stop.'
         stop
   end select
!..auxiliary constant
   za = ZI*OMEGA*OMEGA_RATIO_FLD*EPSILON + SIGMA
   zb = -conjg(za)
   zc = ZI*OMEGA*OMEGA_RATIO_FLD*MU
!
   Jstretch = ZERO
   Jstretch(1,1) = ZONE
   Jstretch(2,2) = ZONE
!
   INVJStretch = ZERO
   INVJStretch(1,1) = ZONE
   INVJStretch(2,2) = ZONE
!
   JJStretch = ZERO
!..initialize the background polarization
   call find_domain(Mdle, ndom)
!..select case of GEOM_NO to set
!..refractive index according to domain
!..fld_flag = 1 when we are in signal element routine
   select case(GEOM_NO)
      case(1)
         bg_pol = ZERO
         raman_pol = ZERO
         therm_pol = ZERO
      case(2,3)
         write(*,*) 'cannot have prism core geometry with fast integration. stop.'
         stop
      case(4)
         bg_pol = ZERO
         raman_pol = ZERO
         therm_pol = ZERO
      case(5)
         raman_pol = ZERO
         therm_pol = ZERO
         select case(ndom)
            case(1,2,3,4,5)
               dom_flag = 1
               call get_bgPol(dom_flag,fld_flag, bg_pol)
            case(6,7,8,9)
               dom_flag = 0
               call get_bgPol(dom_flag,fld_flag, bg_pol)
            case default
               write(*,*) 'unexpected ndom param. stop.'
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
!
!..set up # dof for each direction for 1D H1 test functions with order p+dp
   nrdofH1=nord1+1; nrdofH2=nord2+1; nrdofH3=nord3+1
!..set up # dof for each direction for 1D H1 trial functions with order p
   nrdofH1_tr=nrdofH1-NORD_ADD; nrdofH2_tr=nrdofH2-NORD_ADD; nrdofH3_tr=nrdofH3-NORD_ADD
!..set up # dof for each direction for 1D L2 trial functions with order p
   nrdofQ1_tr=nrdofH1_tr-1; nrdofQ2_tr=nrdofH2_tr-1; nrdofQ3_tr=nrdofH3_tr-1
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
   allocate(STIFQE_A(3,3,nrdofQ3_tr*nrdofH3))
   allocate(STIFQE_B(3,3,nrdofQ2_tr*nrdofH2,nrdofQ3_tr*nrdofH3))
   allocate(STIFQE_ALPHA_A(3,3,nrdofQ3_tr*nrdofH3))
   allocate(STIFQE_ALPHA_B(3,3,nrdofQ2_tr*nrdofH2,nrdofQ3_tr*nrdofH3))
   allocate(STIFQC_A(2,3,3,nrdofQ3_tr*nrdofH3))
   allocate(STIFQC_B(2,3,3,nrdofQ2_tr*nrdofH2,nrdofQ3_tr*nrdofH3))
   allocate(LOADE_A(3,nrdofH3))
   allocate(LOADE_B(3,nrdofH2,nrdofH3))
!
!..TODO stefan [consistency check in debug here instead of reassigning value]
!..total number of test functions (per field)
   !nrdofEE= nord1*nrdofH2*nrdofH3+nrdofH1*nord2*nrdofH3+nrdofH1*nrdofH2*nord3
!..number of trial functions (per component)
   !nrdofQ = nrdofQ1_tr*nrdofQ2_tr*nrdofQ3_tr
!..consistency check
!   if (nrdofEtest .ne. nrdofEE) then
!      write(*,*) 'inconsistent values for:'
!      write(*,*) 'nrdofEtest', nrdofEtest
!      write(*,*) 'nrdofEE   ', nrdofEE
!      write(*,*) 'stop.'
!      stop
!   endif
!   if (nrdofQtrial .ne. nrdofQ) then
!      write(*,*) 'inconsistent values for:'
!      write(*,*) 'nrdofQtrial', nrdofQtrial
!      write(*,*) 'nrdofQ     ', nrdofQ
!      write(*,*) 'stop.'
!      stop
!   endif
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
!  ...Initialize auxiliary matrices B: Stiffness and load
      STIFQE_B      =ZERO
      STIFQE_ALPHA_B=ZERO
      STIFQC_B      =ZERO
      LOADE_B       =ZERO
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
!     ...Initialize auxiliary matrices A: Stiffness matrix and load vector
         STIFQE_A      =ZERO
         STIFQE_ALPHA_A=ZERO
         STIFQC_A      =ZERO
         LOADE_A       =ZERO
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
!   ........compute current solution using soleval
            nflag = 1
            call soleval(Mdle,xip,norient_edge,norient_face,norder,xnod,&
              zdofH,zdofE,zdofV,zdofQ,nflag, x,dxdxi,zsolH_soleval,zdsolH_soleval,&
              zsolE_soleval,zcurlE_soleval,zsolV_soleval,zdivV_soleval,zsolQ_soleval)
!
!        ...compute total quadrature weight
            wt123=wt1*wt2*wt3
!        ...compute Jacobian determinant * quadrature weight
            weighthh=wt123*rjac
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
!        ...compute inverse Jacobian determinant * quadrature weight
            weightvv=wt123/rjac
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
!        ...get the RHS
            call getf(Mdle,x, zfval,zJ)
!        ...get PML function, fld_flag = 1 for signal
            call get_Beta(x,fld_flag, zbeta,zdbeta,zd2beta)
!.....................................................
!...............toggle PML............................
            if(USE_PML.eq.0) then
              zdbeta = ZONE
            endif
            !!!write(*,*) 'from elem_uwMaxwell_fi_TRANS:zdbeta is ', zdbeta
!...............toggle PML............................
!.....................................................
            Jstretch(3,3) = zdbeta
            INVJStretch(3,3) = 1.d0/zdbeta
            call zgemm('N', 'N', 3, 3, 3, ZONE, INVJStretch, 3, INVJStretch, 3, ZERO, JJStretch, 3)
            detJstretch = zdbeta
            JJStretch = detJstretch*JJStretch

!         ..Section to determine coefficients from heat equation solution
!
! ..........FOR THE NONLINEAR_FLAG = 1
            if(NONLINEAR_FLAG.eq.1) then
!           ..initialize raman polarization
              raman_pol = ZERO
!           ..fld_flag = 1 when in signal element routine
              if(fld_flag.eq.1) then
                call get_ramanPol(zsolQ_soleval(7:9),zsolQ_soleval(10:12), &
                                  dom_flag,fld_flag, raman_pol)
              else
                call get_ramanPol(zsolQ_soleval(1:3),zsolQ_soleval(4:6), &
                                  dom_flag,fld_flag, raman_pol)
              endif
              if(LASER_MODE.eq.1) then
!             ..compute the H1 solution at the point
                call get_thermPol(dom_flag,fld_flag,zsolH_soleval(1), therm_pol)
!           ..end if LASER_MODE check
              endif
!         ..end if NONLINEAR_FLAG check
            endif
!         ..update auxiliary constants for za,zb,zc: this is for
!         ..Stiffness and Gram matrix that changes with each nonlinear iteration
            za = ZI*OMEGA*OMEGA_RATIO_FLD*EPSILON+SIGMA+bg_pol+raman_pol+therm_pol
            zb = -conjg(za)
            zc = ZI*OMEGA*OMEGA_RATIO_FLD*MU
            !... writing bg_pol and raman_pol
            !!! write(*,*) 'bg_pol from elem_uwMaxwell_fi_TRANS = ', bg_pol
            !!! write(*,*) 'raman_pol from elem_uwMaxwell_fi_TRANS = ', raman_pol
            !!! write(*,*) 'za from elem_uwMaxwell_fi_TRANS = ', za
            !!! write(*,*) 'fld_flag, NO_PROBLEM = ', fld_flag, NO_PROBLEM
            !!! call pause
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
                           STIFQE_A(a,b,k3) = STIFQE_A(a,b,k3) &
                                            + (dxidx(a,b)*shapH3(j3+1,2)) &
                                            * shapH3(idxa,sa)*zc
!                       ...accumulate innermost 1D integral for alpha_scale*QE term
                           STIFQE_ALPHA_A(a,b,k3) = STIFQE_ALPHA_A(a,b,k3) &
                                             - (dxidx(a,b)*shapH3(j3+1,2)) &
                                             * shapH3(idxa,sa)*(za)
!                       ...loop over components a+alph, required for curl of test f
                           do alph=1,2
                              idxalph=mod(a+alph-1,3)+1
                              sa=1+1-deltak(idxalph,3)
!                          ...accumulate innermost 1D integral for QC term
                              STIFQC_A(alph,a,b,k3) = STIFQC_A(alph,a,b,k3) &
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
                     LOADE_A(a,i3)=LOADE_A(a,i3) &
                              +(dxidx(a,1)*zJ(1)+dxidx(a,2)*zJ(2)+dxidx(a,3)*zJ(3)) &
                              *shapH3(idxa,sa)*rjac
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
                              STIFQE_B(a,b,k2,k3)=STIFQE_B(a,b,k2,k3)             &
                                                 +shapH2(j2+1,2)*shapH2(idxa,sa)  &
                                                 *STIFQE_A(a,b,k3)
!                          ...accumulate za*QE term
                              STIFQE_ALPHA_B(a,b,k2,k3)=STIFQE_ALPHA_B(a,b,k2,k3) &
                                                 +shapH2(j2+1,2)*shapH2(idxa,sa)  &
                                                 *STIFQE_ALPHA_A(a,b,k3)
!                          ...accumulate QC term
                              do alph=1,2
                                 idxalph=mod(a+alph-1,3)+1
                                 sa=1+1-deltak(idxalph,2)
                                 STIFQC_B(alph,a,b,k2,k3)=STIFQC_B(alph,a,b,k2,k3) &
                                                 + shapH2(j2+1,2)*shapH2(idxa,sa)  &
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
                     LOADE_B(a,i2,i3)=LOADE_B(a,i2,i3)+shapH2(idxa,sa)*LOADE_A(a,i3)
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
                                    STIFFEQ_T(k,2*m1-1) = STIFFEQ_T(k,2*m1-1) &
                                                      + shapH1(idxa,sa)*shapH1(j1+1,2) &
                                                      * STIFQC_B(alph,a,b,k2,k3)
                                 enddo

                                 sa=1+deltak(a,1)

                                 k = (m2-1)*6 + b
!                             ...accumulate QE term
!                             ... stretch for PML
                                 STIFFEQ_T(k,2*m1-1) = STIFFEQ_T(k,2*m1-1) &
                                                   +shapH1(idxa,sa)*shapH1(j1+1,2) &
                                                   *JJStretch(a,a)*STIFQE_ALPHA_B(a,b,k2,k3)
                                 k = (m2-1)*6 + 3+ b
!                             ...accumulate QE term
                                 STIFFEQ_T(k,2*m1) = STIFFEQ_T(k,2*m1) &
                                                   +shapH1(idxa,sa)*shapH1(j1+1,2) &
                                                   *JJStretch(a,a)*STIFQE_B(a,b,k2,k3)
                                 k = (m2-1)*6 + b
!                             ...accumulate QC term
                                 do alph=1,2
                                    idxalph=mod(a+alph-1,3)+1
                                    sa=1+1-deltak(idxalph,1)
                                    STIFFEQ_T(k,2*m1) = STIFFEQ_T(k,2*m1) &
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
                     BLOADE(2*m1-1)=BLOADE(2*m1-1)+shapH1(idxa,sa)*LOADE_B(a,i2,i3)
                  endif
               enddo
            enddo
         enddo
      enddo
!  ...loop over px ends
   enddo
!
!..printing Gram matrix
   iprint = 0
   if (iprint.eq.1) then
      write(*,*) 'elem: AP_Maxwell = '
      !!!do i=1,10
      do i=70,90
         do j=70,i-1
            zaux(j-69) = AP_Maxwell(nk(j,i))
         enddo
         do j=i,75
            zaux(j-69) = AP_Maxwell(nk(i,j))
         enddo
         write(*,7011) zaux(1:21)
      enddo
      call pause
   endif
!
   deallocate(AUXEE_A)
   deallocate(AUXEE_B)
   deallocate(AUXCC_A)
   deallocate(AUXCC_B)
   deallocate(AUXEC_A)
   deallocate(AUXEC_B)
   deallocate(AUXCE_A)
   deallocate(AUXCE_B)
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
   call ndof_nod(etype,nordP,ndofHHmdl,ndofEEmdl,ndofVVmdl,ndofQQmdl)
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
      call set_2Dint(ftype,norderf, nint,tloc,wtloc)
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
         call shape3EE(etype,xi,nordP, nrdofEE,shapEE,curlEE)
!
!     ...determine element H1 shape functions (for geometry)
         call shape3H(etype,xi,norder,norient_edge,norient_face, &
                     nrdofH,shapH,gradH)
!
!     ...determine element H(curl) shape functions (for fluxes)
!     ...for interfaces only (no bubbles)
         call shape3E(etype,xi,norderi,norient_edge,norient_face, &
                     nrdofEi,shapE,curlE)
!
!     ...geometry
         call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign, &
                     x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
         weight = bjac*wtloc(l)

!     ...loop through Hcurl enriched test functions
         do ik=1,NrdofEEi
            k1 = idxEE(ik)
            E1(1:3) = shapEE(1,k1)*dxidx(1,1:3) &
                    + shapEE(2,k1)*dxidx(2,1:3) &
                    + shapEE(3,k1)*dxidx(3,1:3)
!        ...check for impedance BC
!        ...impedance boundary
            if ((ibc(ifc,2).eq.9).or.(ibc(ifc,3).eq.9)) then
!           ...get the boundary source
               call get_bdSource(Mdle,x,rn, zImp)
!           ...accumulate for the load vector
               k = 2*k1-1
               BLOADE(k) = BLOADE(k) &
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
               if ((ibc(ifc,2).eq.9).or.(ibc(ifc,3).eq.9)) then
!              ...accumulate for the extended stiffness matrix on IBC
                  call cross_product(rn,rntimesE, rn2timesE)
                  STIFFEE_T(2*k2-1,2*k1-1) = STIFFEE_T(2*k2-1,2*k1-1) &
                                               + (E1(1)*rn2timesE(1) &
                                               +  E1(2)*rn2timesE(2) &
                                               +  E1(3)*rn2timesE(3) &
                                                 )*GAMMA*weight
!
                  STIFFEE_T(2*k2-1,2*k1) = STIFFEE_T(2*k2-1,2*k1) &
                                             + (E1(1)*rntimesE(1) &
                                             +  E1(2)*rntimesE(2) &
                                             +  E1(3)*rntimesE(3) &
                                               )*weight
               else
!              ...accumulate for the extended stiffness matrix without IBC
                  STIFFEE_T(2*k2,2*k1-1) = STIFFEE_T(2*k2,2*k1-1) &
                                             + (E1(1)*rntimesE(1) &
                                             +  E1(2)*rntimesE(2) &
                                             +  E1(3)*rntimesE(3) &
                                               )*weight
!
                  STIFFEE_T(2*k2-1,2*k1) = STIFFEE_T(2*k2-1,2*k1) &
                                             + (E1(1)*rntimesE(1) &
                                             +  E1(2)*rntimesE(2) &
                                             +  E1(3)*rntimesE(3) &
                                               )*weight
!           ...end if for impedance BC
               endif
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
!!..TODO needs to be fixed (use transpose instead)
!!..printing element matrices
!   iprint = 0
!   if (iprint.eq.1) then
!      write(*,*) 'elem: STIFFEE, STIFFEQ, BLOADE = '
!      do j=1,2*nrdofEE
!         write(*,7200) j, BLOADE(j)
!         write(*,*) 'STIFFEE for n x E = '
!         write(*,7210) (STIFFEE(j,2*k-1),k=1,nrdofE)
!         write(*,*) 'STIFFEE for n x H = '
!         write(*,7210) (STIFFEE(j,2*k  ),k=1,nrdofE)
!         write(*,*) 'STIFFEQ = '
!         write(*,7210) STIFFEQ(j,1:6*nrdofQ)
!7200     format('j = ',i5,' BLOAD(j) = ',2e12.5)
!7210     format(6(2e12.5,1x))
!         write(*,*) 'residual = '
!         if (j/2*2.eq.j) then
!            write(*,7210) BLOADE(j) - STIFFEQ(j,1)
!         else
!            write(*,7210) BLOADE(j) - STIFFEQ(j,1) &
!              - (STIFFEE(j,1)+STIFFEE(j,5)+STIFFEE(j,9)+STIFFEE(j,13))
!         endif
!         call pause
!      enddo
!   endif
!
!! ....................................................................
!! ...Construction of DPG system
!! ....................................................................
!
!..Total test/trial DOFs of the element
   i1 = nrTest ; j1 = 2*nrdofE ; j2 = 6*nrdofQ
!
   STIFF_ALL(1:i1,1:j1) = transpose(STIFFEE_T(1:i1,1:j1))
   STIFF_ALL(1:i1,j1+1:j1+j2) = transpose(STIFFEQ_T(1:i1,1:j2))
   STIFF_ALL(1:i1,j1+j2+1) = BLOADE(1:i1)
!
!----------------------------------------------------------------------
!----Preconditioning---------------------------------------------------
!----------------------------------------------------------------------
!   allocate(DIAG_E(nrTest))
!!..preconditioning: getting diagonal entries of Gram matrix
!   do k1=1,nrTest
!      k = nk(k1,k1)
!      DIAG_E(k1) = AP_Maxwell(k)
!   enddo
!!..preconditioning: AP_Maxwell= D^-1/2 * AP_Maxwell * D^-1/2
!   do k1=1,nrTest
!      do k2=k1,nrTest
!         k = nk(k1,k2)
!         AP_Maxwell(k) = AP_Maxwell(k)/sqrt(DIAG_E(k1)*DIAG_E(k2))
!      enddo
!   enddo
!!..preconditioning: STIFF_ALL = D^-1/2 * STIFF_ALL
!   do k2=1,nrTrial+1
!      do k1=1,nrTest
!         STIFF_ALL(k1,k2) = STIFF_ALL(k1,k2)/sqrt(DIAG_E(k1))
!      enddo
!   enddo
!   deallocate(DIAG_E)
!
!----------------------------------------------------------------------
!----Condition number of Gram matrix-----------------------------------
!----------------------------------------------------------------------
!
!!..check condition number
!   nrdof_eig = nrdofEE*2
!   kk = nrdof_eig*(nrdof_eig+1)/2
!   allocate(AP_eig(kk))
!   AP_eig(1:kk) = AP_Maxwell(1:kk)
!   allocate(W(nrdof_eig))
!   allocate(Z(1,nrdof_eig))
!   allocate(WORK(nrdof_eig))
!   allocate(RWORK(nrdof_eig))
!   allocate(IWORK(1))
!   call ZHPEVD('N','U',nrdof_eig, &
!               AP_eig, W,Z,1,WORK,nrdof_eig, &
!               RWORK,nrdof_eig,IWORK,1,info)
!   if (info .ne. 0) then
!      write(*,*) 'eig_solve_sc: info = ', info
!      stop 1
!   endif
!
!   write(*,6999) W(nrdof_eig),W(1)
!6999   format('elem_dpgMaxwell: AP_Maxwell: max_eig, min_eig = ', 2e13.4)
!
!
!   write(*,7000)  W(nrdof_eig)/W(1)
!7000   format('elem_dpgMaxwell: AP_Maxwell condition number = ', 1e13.4)
!   deallocate(IWORK,W,WORK)
!   deallocate(RWORK,Z)
!   deallocate(AP_eig)
!   pause
!
!---------------------------
!  STANDARD APPROACH: B^*B
!---------------------------
!
!!!..Cholesky factorization of Gram matrix: G=U^*U (=LL^*)
!   call ZPPTRF('U',nrTest,AP_Maxwell,info)
!   if (info.ne.0) then
!      write(*,*) 'ZPPTRF: info = ',info,'. stop.'
!      stop
!   endif
!!..Solve triangular system, (LX=) U^*X = [B|l]
!   call ZTPTRS('U','C','N',nrTest,nrTrial+1,AP_Maxwell,STIFF_ALL, &
!               nrTest,info)
!   if (info.ne.0) then
!      write(*,*) 'ZTPTRS: info = ',info,'. stop.'
!      stop
!   endif
!!
!!..ZHERK for complex case, zaloc=X^*X
!   call ZHERK('U','C',nrTrial+1,nrTest,ZONE,STIFF_ALL, &
!              nrTest,ZERO, &
!              zaloc(1:nrTrial+1,1:nrTrial+1),nrTrial+1)
!   if (info.ne.0) then
!      write(*,*) 'ZHERK: info = ',info,'. stop.'
!      stop
!   endif
!!..do conjugate transpose of zaloc
!   do i=1,nrTrial
!      zaloc(i+1:nrTrial+1,i) = conjg(zaloc(i,i+1:nrTrial+1))
!   enddo
!!
!   ZblocE(1:j1) = zaloc(1:j1,j1+j2+1)
!   ZblocQ(1:j2) = zaloc(j1+1:j1+j2,j1+j2+1)
!!
!   ZalocEE(1:j1,1:j1) = zaloc(1:j1,1:j1)
!   ZalocEQ(1:j1,1:j2) = zaloc(1:j1,j1+1:j1+j2)
!!
!   ZalocQE(1:j2,1:j1) = zaloc(j1+1:j1+j2,1:j1)
!   ZalocQQ(1:j2,1:j2) = zaloc(j1+1:j1+j2,j1+1:j1+j2)
!
!
!------------------------
!  BLOCK APPROACH: B^*B
!------------------------
!
!..Cholesky factorization of Gram matrix: G=U^*U (=LL^*)
   call ZPPTRF('U',nrTest,AP_Maxwell,info)
!
   if (info.ne.0) then
      write(*,*) 'elem: gram ZPPTRF: Mdle,info = ',Mdle,info
      stop 1
   endif
!
!..triangular solve left block (1) and right block (3) of STIFF_ALL
   call ZTPTRS('U','C','N',nrTest,2*nrdofEi,AP_Maxwell, &
     STIFF_ALL(1:nrTest,1:2*nrdofEi),nrTest,info)
   call ZTPTRS('U','C','N',nrTest,6*NrdofQ+1,AP_Maxwell, &
     STIFF_ALL(1:nrTest,2*NrdofE+1:nrTrial+1),nrTest,info)
!
!..block (1,1): H(curl) interface x H(curl) interface
   call ZHERK('U','C',2*nrdofEi,nrTest,1.0d0,STIFF_ALL(1:nrTest,1:2*nrdofEi), &
                 nrTest,ZERO,zaloc(1:2*nrdofEi,1:2*nrdofEi),2*nrdofEi)
!
!..block (1,2): H(curl) interface x H(curl) interior
!..THIS BLOCK IS SKIPPED BECAUSE IT'S ZERO ANYWAY
!
!..block (1,3): H(curl) interface x L2 interior
   call ZGEMM('C','N',2*nrdofEi,6*NrdofQ+1,nrTest,ZONE,STIFF_ALL(1:nrTest,1:2*nrdofEi), &
               nrTest, STIFF_ALL(1:nrTest,2*NrdofE+1:nrTrial+1),nrTest,ZERO, &
               zaloc(1:2*nrdofEi,2*NrdofE+1:nrTrial+1),2*nrdofEi)
!
!..block (2,:): H(curl) interior x anything 
!..THESE BLOCKS ARE SKIPPED BECAUSE THEY ARE ZERO ANYWAY
!
!..block (3,1): L2 interior x H(curl) interface
!..already computed = block (1,3)^*
!
!..block (3,2): L2 interior x H(curl) interior
!..THIS BLOCK IS SKIPPED BECAUSE IT'S ZERO ANYWAY
!
!..block (3,3): L2 interior x L2 interior
   call ZHERK('U','C',6*nrdofQ+1,nrTest,1.0d0,STIFF_ALL(1:nrTest,2*NrdofE+1:nrTrial+1), &
               nrTest,ZERO,zaloc(2*NrdofE+1:nrTrial+1,2*NrdofE+1:nrTrial+1),6*NrdofQ+1)
!
   do i=1,nrTrial
      zaloc(i+1:nrTrial+1,i) = conjg(zaloc(i,i+1:nrTrial+1))
   enddo
!
   ZblocE(1:j1) = zaloc(1:j1,j1+j2+1)
   ZblocQ(1:j2) = zaloc(j1+1:j1+j2,j1+j2+1)
!
   ZalocEE(1:j1,1:j1) = zaloc(1:j1,1:j1)
   ZalocEQ(1:j1,1:j2) = zaloc(1:j1,j1+1:j1+j2)
!
   ZalocQE(1:j2,1:j1) = zaloc(j1+1:j1+j2,1:j1)
   ZalocQQ(1:j2,1:j2) = zaloc(j1+1:j1+j2,j1+1:j1+j2)
!
!----------------------
!  END BLOCK APPROACH
!----------------------
!
!.................LONG WAVEGUIDE.............................
  ! if (idec.ne.2) then
  !  write(*,*) 'elem_dpgMaxwell: idec = ',idec
  !  stop 1
  ! endif
  ! call copy_elmmat_UW(MdE,MdQ,xnod,ibc, idec, &
  !                           ZalocEE,ZalocEQ, &
  !                           ZalocQE,ZalocQQ,ZblocE,ZblocQ)
!.................LONG WAVEGUIDE.............................
!
   iprint = 0
   if (iprint.ge.1) then
      write(*,7010)
7010  format('elem_dpgMaxwell: ZblocE,ZblocQ = ')
      write(*,7011) ZblocE(1:2*NrdofE)
      write(*,7011) ZblocQ(1:6*NrdofQ)
7011  format(12e12.5)
      pause
      write(*,7012)
7012  format('elem_dpgMaxwell: ZalocEE = ')
      do i=1,2*NrdofE
         write(*,7013) i,ZalocEE(i,1:2*NrdofE)
7013     format('i = ',i3,10(/,6(2e12.5,2x)))
      enddo
      pause
      write(*,7014)
7014  format('elem_dpgMaxwell: ZalocEQ = ')
      do i=1,2*NrdofE
         write(*,7013) i,ZalocEQ(i,1:6*NrdofQ)
      enddo
      pause
      write(*,7015)
7015  format('elem_dpgMaxwell: ZalocQQ = ')
      do i=1,6*NrdofQ
         write(*,7013) i,ZalocQQ(i,1:6*NrdofQ)
      enddo
      pause
   endif
!
end subroutine elem_uwMaxwell_fi_TRANS


! subroutine works only for hexa currently
subroutine compute_enriched_order(Nord, Norder)
!
   use parameters, only : MODORDER
!
   implicit none
!
   integer, intent(in)  :: Nord
   integer, intent(out) :: Norder(19)
!
   integer :: temp(2)
   integer :: nordF(3),nordB(3),ndofH(3)
!
   call decod(Nord,MODORDER,2, temp)
   nordF(1) = temp(1) ; nordB(3) = temp(2)
   call decod(nordF(1),MODORDER,2, nordB(1:2))
   call encod((/nordB(1),nordB(3)/),MODORDER,2, nordF(2))
   call encod(nordB(2:3),MODORDER,2, nordF(3))
   Norder(1:4)   = (/nordB(1),nordB(2),nordB(1),nordB(2)/)
   Norder(5:8)   = (/nordB(1),nordB(2),nordB(1),nordB(2)/)
   Norder(9:12)  = nordB(3)
   Norder(13:14) = nordF(1)
   Norder(15:18) = (/nordF(2),nordF(3),nordF(2),nordF(3)/)
   Norder(19)    = Nord
!
end subroutine compute_enriched_order

