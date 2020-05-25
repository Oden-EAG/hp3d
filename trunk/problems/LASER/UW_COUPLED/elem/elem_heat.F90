!
!--------------------------------------------------------------------
!
!     routine name      - elem_heat
!
!--------------------------------------------------------------------
!
!     latest revision:  - Mar 2019
!
!     purpose:          - routine returns unconstrained DPG
!                         stiffness matrix and load vector
!                         for the primal DPG formulation for a single
!                         step of the transient heat equation
!
!     arguments:
!        in:
!              Mdle     - element middle node number
!              NrTest   - total number of test dof
!              NrTrial  - total number of trial dof
!              NrdofHH  - number of H1 test dof
!              NrdofH   - number of H1 trial dof
!              NrdofV   - number of H(div) trial dof
!              NrdofVi  - number of H(div) trial interface dof
!              MdH      - num rows of ZalocHH,ZalocHV
!              MdV      - num rows of ZalocVH,ZalocVV
!        out:
!              ZblocH   - load vectors
!              ZblocV
!              ZalocHH  - stiffness matrices
!              ZalocHV
!              ZalocVH
!              ZalocVV
!
!---------------------------------------------------------------------
!
#include "implicit_none.h"
!
subroutine elem_heat(Mdle,                   &
                     NrTest,NrTrial,         &
                     NrdofHH,                &
                     NrdofH,NrdofV,          &
                     NrdofVi,                &
                     MdH,MdV,                &
                     ZblocH,ZalocHH,ZalocHV, &
                     ZblocV,ZalocVH,ZalocVV)
!..modules used
   use control
   use parametersDPG
   use element_data
   use data_structure3D
   use laserParam
   use commonParam
!..no implicit statements
   implicit none
!..declare input/output variables
   integer,                     intent(in)  :: Mdle
   integer,                     intent(in)  :: NrTest
   integer,                     intent(in)  :: NrTrial
   integer,                     intent(in)  :: NrdofHH
   integer,                     intent(in)  :: NrdofH
   integer,                     intent(in)  :: NrdofV
   integer,                     intent(in)  :: NrdofVi
   integer,                     intent(in)  :: MdH
   integer,                     intent(in)  :: MdV
   VTYPE, dimension(MdH),       intent(out) :: ZblocH
   VTYPE, dimension(MdH,MdH),   intent(out) :: ZalocHH
   VTYPE, dimension(MdH,MdV),   intent(out) :: ZalocHV
   VTYPE, dimension(MdV),       intent(out) :: ZblocV
   VTYPE, dimension(MdV,MdH),   intent(out) :: ZalocVH
   VTYPE, dimension(MdV,MdV),   intent(out) :: ZalocVV
!
!..auxiliary parameters
   real(8), parameter :: rZERO = 0.d0
   real(8), parameter :: rONE  = 1.d0
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
   real(8), dimension(3,MAXbrickH) :: xnod
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
!..H(div) shape functions
   real(8), dimension(3,MAXbrickV)  :: shapV
   real(8), dimension(MAXbrickV)    :: divV
!
!..Enriched H1 shape functions
   real(8) , dimension(MAXbrickHH)    :: shapHH
   real(8) , dimension(3,MAXbrickHH)  :: gradHH
!
!..load vector for the enriched space
   real(8) :: bload_H(NrTest)
!
!..gram matrix in packed format
   !real(8) :: gramP(NrTest*(NrTest+1)/2)
   real(8), allocatable :: gramP(:)
   !complex(8), allocatable :: gramEigen(:)
!
!..stiffness matrices for the enriched test space
!   real(8) :: stiff_HH (NrTest   ,NrdofH)
!   real(8) :: stiff_HV (NrTest   ,NrdofVi)
!   real(8) :: stiff_ALL(NrTest   ,NrTrial+1)
!   real(8) :: raloc    (NrTrial+1,NrTrial+1)
   real(8), allocatable :: stiff_HH(:,:),stiff_HV(:,:),stiff_ALL(:,:),raloc(:,:)
!
!..3D quadrature data
   real(8), dimension(3,MAXNINT3ADD)  :: xiloc
   real(8), dimension(  MAXNINT3ADD)  :: waloc
!
!..2D quadrature data
   real(8), dimension(2,MAXNINT2ADD)  :: tloc
   real(8), dimension(  MAXNINT2ADD)  :: wtloc
!
!..derivatives wrt physical coordinates, flux
   real(8), dimension(3) :: dv1,dv2,vec
!
!..auxiliary
   !real(8) :: raux(10)
   complex(8) :: zvoid(3)
!
!..number of vertices,edge,faces per element type
   integer :: nrv,nre,nrf
!
!..various variables for the problem
   real(8)    :: rjac,bjac,weight,wa,v2n,v1,v2
   integer    :: i1,i2,j1,j2,k1,k2,kH,kk,i,j,nint,iflag,kE,k
   integer    :: nordP,nrdof,l,nsign,ifc,info,ndom
   real(8)    :: rfval,therm_Load
   complex(8) :: zfval
!
!..for lapack eigensolve
   ! complex(8), allocatable :: Z(:,:), WORK(:)
   ! real(8), allocatable     :: W(:),   RWORK(:)
   ! integer, allocatable    :: IWORK(:)
!
#if DEBUG_MODE
!..BC's flags
   integer, dimension(6,NR_PHYSA) :: ibc
!..auxiliary
   integer :: iprint
#endif
!
!..for Gram matrix compressed storage format
   integer :: nk
   nk(k1,k2) = (k2-1)*k2/2+k1
!
!---------------------------------------------------------------------
!--------- INITIALIZE THE ELEMENT ORDER, ORIENTATION ETC. ------------
!---------------------------------------------------------------------
!
#if DEBUG_MODE
!..Set iprint = 0/1 (Non-/VERBOSE)
   iprint = 0
   if (iprint.eq.1) then
      write(*,*) 'elem_heat: Mdle = ', Mdle
   endif
#endif
!
!..allocate auxiliary matrices
   allocate(gramP(NrTest*(NrTest+1)/2))
   allocate(stiff_HH (NrTest   ,NrdofH))
   allocate(stiff_HV (NrTest   ,NrdofVi))
!
!..element type
   etype = NODES(Mdle)%type
   nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
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
      case('mdln','mdld')
         nordP = NODES(Mdle)%order+NORD_ADD
         norderi(nre+nrf+1) = 1
      case('mdlp')
         nordP = NODES(Mdle)%order+NORD_ADD*11
         norderi(nre+nrf+1) = 11
      case default
         write(*,*) 'elem_heat: invalid etype param. stop.'
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
      case(1)
!     ...single-step solution assumes zero initial data
      case(2)
         call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
      case default
         write(*,*) 'elem_heat: invalid NO_PROBLEM param. stop.'
         stop
   end select
!
#if DEBUG_MODE
!..get the element boundary conditions flags
   call find_bc(Mdle, ibc)
   if (iprint.eq.1) then
      write(*,7001) Mdle
7001  format('elem_heat: BCFLAGS FOR Mdle = ',i5)
      do i=1,NR_PHYSA
         write(*,7002) PHYSA(i), ibc(1:nrf,i)
7002     format('          ATTRIBUTE = ',a6,' FLAGS = ',6i2)
      enddo
   endif
#endif
!
!..determine element domain
   call find_domain(Mdle, ndom)
!
!..clear space for stiffness matrix and rhsv:
   ZblocH = ZERO; ZalocHH = ZERO; ZalocHV = ZERO
   ZblocV = ZERO; ZalocVH = ZERO; ZalocVV = ZERO
!
!..clear space for auxiliary matrices
   bload_H   = rZERO
   gramP     = rZERO
   stiff_HH  = rZERO
   stiff_HV  = rZERO
!
!---------------------------------------------------------------------
!              E L E M E N T    I N T E G R A L S                    |
!---------------------------------------------------------------------
!
!..use the enriched order to set the quadrature
   INTEGRATION = NORD_ADD
   call set_3D_int_DPG(etype,norder,norient_face, nint,xiloc,waloc)
   INTEGRATION = 0
!
!..loop over integration points
   do l=1,nint
      xi(1:3)=xiloc(1:3,l); wa=waloc(l)
!
!  ...H1 shape functions
      call shape3H(etype,xi,norder,norient_edge,norient_face, &
        nrdof,shapH,gradH)
#if DEBUG_MODE
      if (nrdof .ne. NrdofH) then
         write(*,*) 'elem_heat: INCONSISTENCY NrdofH. stop.'
         stop
      endif
#endif
!
!  ...discontinuous H1 shape functions
      call shape3HH(etype,xi,nordP, nrdof,shapHH,gradHH)
#if DEBUG_MODE
      if (nrdof .ne. NrdofHH) then
         write(*,*) 'elem_heat: INCONSISTENCY NrdofHH. stop.'
         stop
      endif
#endif
!
!  ...geometry computations
      call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH, &
        x,dxdxi,dxidx,rjac,iflag)
!
!  ...compute current solution using soleval
      nflag = 1
      call soleval(Mdle,xi,norient_edge,norient_face,norder,xnod, &
        zdofH,zdofE,zdofV,zdofQ,nflag,x,dxdxi,zsolH_soleval,zdsolH_soleval, &
        zsolE_soleval,zcurlE_soleval,zsolV_soleval,zdivV_soleval,zsolQ_soleval)
      rsolH = real(zsolH_soleval(1))
!
!  ...coupled Maxwell/Heat problem
      therm_Load = rZero
      if(NONLINEAR_FLAG .eq. 1) then
         if (GEOM_NO .ne. 5) then
            write(*,*) 'elem_heat: HEAT coupling GEOM_NO .ne. 5'
            stop
         endif
!     ...heat deposition in fiber core only
         if (ndom .eq. 1 .or. ndom .eq.  2) then
!        ...compute the H1 solution load at the point
            call get_thermLoad(zsolQ_soleval, therm_Load)
         endif
      endif
!
!  ...integration weight
      weight = rjac*wa
!
!  ...get the RHS
!     ...returns zero fval for coupled problem,
!     ...returns non-zero fval for manufactured problem
      call getf(Mdle,x, zfval,zvoid)
      rfval = real(zfval)
!
!  ...1st loop through enriched H1 test functions
      do k1=1,NrdofHH
!
         v1 = shapHH(k1)
         dv1(1:3) = gradHH(1,k1)*dxidx(1,1:3) &
                  + gradHH(2,k1)*dxidx(2,1:3) &
                  + gradHH(3,k1)*dxidx(3,1:3)
!
!     ...EXTENDED LOAD VECTOR
         select case(NO_PROBLEM)
            case(1)
               bload_H(k1) = bload_H(k1) + rfval*v1*weight
            case(2)
               bload_H(k1) = bload_H(k1) +            &
                             ( rsolH + DELTA_T*therm_Load &
                             )*v1*weight
!           ...do not solve for heat equation in PML
               if(COPUMP.eq.1) then
                  if(x(3).ge.PML_REGION) bload_H(k1) = rZERO
               elseif(COPUMP.eq.0) then
                  write(*,*) 'elem_heat: heat equation not enabled for counter pump configuration. stop.'
                  stop
               endif
         end select
!
!     ...2nd loop through enriched H1 trial functions
!     ...for Gram matrix
         do k2=k1,NrdofHH
!
            v2 = shapHH(k2)
            dv2(1:3) = gradHH(1,k2)*dxidx(1,1:3) &
                     + gradHH(2,k2)*dxidx(2,1:3) &
                     + gradHH(3,k2)*dxidx(3,1:3)
!
!------------------ G R A M   M A T R I X --------------------------------
!
!        ...determine index in triangular format
            k = nk(k1,k2)
!
!        ...problem-dependent inner product
!        ...primal DPG heat: H1 inner product
!           (TODO: try out anisotropic inner product, and incl. delta t)
            select case(INNER_PRODUCT)
               case(1)
                  if (ANISO_HEAT .eq. 1) then
                     dv2(3) = ALPHA_Z*ALPHA_Z*dv2(3)
                  endif
                  gramP(k) = gramP(k) &
                           + ( v1*v2 + DELTA_T*ALPHA_0 * &
                                ( dv1(1)*dv2(1) + &
                                  dv1(2)*dv2(2) + &
                                  dv1(3)*dv2(3) ) &
                              ) *weight
               case default
                  write(*,*) 'elem_heat: unexpected INNER_PRODUCT param. stop.'
                  stop
            end select
!
!     ...enddo 2nd loop through enriched H1 trial functions
         enddo
!
!     ...loop through H1 trial functions
!     ...for stiffness matrix
         do k2=1,NrdofH
            v2 = shapH(k2)
            dv2(1:3) = gradH(1,k2)*dxidx(1,1:3) &
                     + gradH(2,k2)*dxidx(2,1:3) &
                     + gradH(3,k2)*dxidx(3,1:3)
!
!        ...account for short fiber scaling (anisotropic diffusion operator)
            if (ANISO_HEAT .eq. 1) then
               dv2(3) = ALPHA_Z*ALPHA_Z*dv2(3)
            endif
!
!        ...EXTENDED HH STIFFNESS MATRIX
            stiff_HH(k1,k2) = stiff_HH(k1,k2) &
                            + ( v1*v2 + DELTA_T*ALPHA_0 * &
                                ( dv1(1)*dv2(1) + &
                                  dv1(2)*dv2(2) + &
                                  dv1(3)*dv2(3) ) &
                              ) *weight
!     ...end loop through H1 trial functions
         enddo
!  ...end 1st loop through enriched H1 test functions
      enddo
!..end loop through integration points
   enddo
!
#if DEBUG_MODE
!..printing Gram matrix
!   iprint = 0
!   if (iprint.eq.1) then
!      write(*,*) 'elem_heat: GRAM MATRIX = '
!      do i=1,10
!         do j=1,i-1
!            raux(j) = gramP(nk(j,i))
!         enddo
!         do j=i,10
!            raux(j) = gramP(nk(i,j))
!         enddo
!         write(*,7011) raux
!      enddo
!      pause
!   endif
#endif
!
!---------------------------------------------------------------------
!    B O U N D A R Y    I N T E G R A L S                            |
!---------------------------------------------------------------------
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
!     ...determine discontinuous H1 shape functions
         call shape3HH(etype,xi,nordP, nrdof,shapHH,gradHH)
#if DEBUG_MODE
         if (nrdof .ne. NrdofHH) then
            write(*,*) 'elem_heat: INCONSISTENCY NrdofHH. stop.'
            stop
         endif
#endif
!
!     ...determine element H1 shape functions (for geometry)
         call shape3H(etype,xi,norder,norient_edge,norient_face, &
                  nrdof,shapH,gradH)
#if DEBUG_MODE
         if (nrdof .ne. NrdofH) then
            write(*,*) 'elem_heat: INCONSISTENCY NrdofH. stop.'
            stop
         endif
#endif
!
!     ...determine element Hdiv shape functions (for fluxes)
!     ...for interfaces only (no bubbles)
         call shape3V(etype,xi,norderi,norient_face, &
                  nrdof,shapV,divV)
#if DEBUG_MODE
         if (nrdof .ne. NrdofVi) then
            write(*,*) 'elem_heat: INCONSISTENCY NrdofV. stop.'
            stop
         endif
#endif
!
!     ...geometry
         call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign, &
                  x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
         weight = bjac*wtloc(l)
!
!     ...loop through enriched H1 test functions
         do k1=1,NrdofHH
            v1 = shapHH(k1)
!
!        ...loop through H(div) trial functions
            do k2=1,NrdofVi
!
!           ...normal component (Piola transformation at work!)
               vec(1:3) = dxdxi(1:3,1)*shapV(1,k2) &
                        + dxdxi(1:3,2)*shapV(2,k2) &
                        + dxdxi(1:3,3)*shapV(3,k2)
               vec(1:3) = vec(1:3)/rjac
               v2n = vec(1)*rn(1)+vec(2)*rn(2)+vec(3)*rn(3)
!
!           ...EXTENDED HV STIFFNESS MATRIX
               stiff_HV(k1,k2) = stiff_HV(k1,k2) &
                               - DELTA_T*ALPHA_0*v1*v2n*weight
!        ...end loop through H(div) trial functions
            enddo
!     ...end loop through enriched H1 test functions
         enddo
!     ...end loop through integration points
      enddo
!..end loop through element faces
   enddo
!
!..Compute condition number of Gram matrix 
      !kk = NrTest*(NrTest+1)/2
      !allocate(gramEigen(kk))
      !gramEigen(1:kk)=gramP(1:kk)
 !      allocate(W(NrTest))
 !      allocate(Z(1,NrTest))
 !      allocate(WORK(NrTest))
 !      allocate(RWORK(NrTest))
 !      allocate(IWORK(1))
 !      call ZHPEVD('N','U',NrTest, &
 !                  gramEigen, W,Z,1,WORK,NrTest, &
 !                  RWORK,NrTest,IWORK,1,info)
 !      if (info .ne. 0) then
 !         write(*,*) 'eig_solve_sc: info = ', info
 !         stop
 !      endif

 !      write(*,6999) W(NrTest),W(1)
 ! 6999 format(' elem_heat: Gram eigenvalues: max_eig, min_eig = ', 2es13.4)

 !      write(*,7000)  W(NrTest)/W(1)
 ! 7000 format(' elem_heat: Gram condition number = ', 1es13.4)
 !      deallocate(IWORK,W,WORK)
 !      deallocate(RWORK,Z)
 !      deallocate(gramEigen)
 !      pause 
!
!---------------------------------------------------------------------
!  Construction of DPG system
!---------------------------------------------------------------------
!
   allocate(stiff_ALL(NrTest,NrTrial+1))
!
!..Total test/trial DOFs of the element
   i1 = NrTest ; j1 = NrdofH ; j2 = NrdofVi
!
   stiff_ALL(1:i1,1:j1)       = stiff_HH(1:i1,1:j1)
   stiff_ALL(1:i1,j1+1:j1+j2) = stiff_HV(1:i1,1:j2)
   stiff_ALL(1:i1,j1+j2+1)    = bload_H(1:i1)
!
   deallocate(stiff_HH,stiff_HV)
!
!..A. Compute Cholesky factorization of Gram Matrix, G=U^T U (=LL^T)
   call DPPTRF('U',NrTest,gramP,info)
   if (info.ne.0) then
      write(*,*) 'elem_heat: DPPTRF: Mdle,info = ',Mdle,info,'. stop.'
      stop
   endif
!
!..B. Solve triangular system to obtain B~, (LX=) U^T X = [B|l]
   call DTPTRS('U','T','N',NrTest,NrTrial+1,gramP,stiff_ALL,NrTest,info)
   if (info.ne.0) then
      write(*,*) 'elem_heat: DTPTRS: Mdle,info = ',Mdle,info,'. stop.'
      stop
   endif
!
   allocate(raloc(NrTrial+1,NrTrial+1)); raloc = rZERO
!
!..C. Matrix multiply: B^T G^-1 B (=B~^T B~)
   call DSYRK('U','T',NrTrial+1,NrTest,rONE,stiff_ALL,NrTest,rZERO,raloc,NrTrial+1)
!
   deallocate(stiff_ALL)
!
!..D. Fill lower triangular part of Hermitian matrix
   do i=1,NrTrial
      raloc(i+1:NrTrial+1,i) = raloc(i,i+1:NrTrial+1)
   enddo
!
!..E. Fill ALOC and BLOC matrices
   ZblocH(1:j1) = raloc(1:j1,j1+j2+1)
   ZblocV(1:j2) = raloc(j1+1:j1+j2,j1+j2+1)
!
   ZalocHH(1:j1,1:j1) = raloc(1:j1,1:j1)
   ZalocHV(1:j1,1:j2) = raloc(1:j1,j1+1:j1+j2)
!
   ZalocVH(1:j2,1:j1) = raloc(j1+1:j1+j2,1:j1)
   ZalocVV(1:j2,1:j2) = raloc(j1+1:j1+j2,j1+1:j1+j2)
!
   deallocate(raloc)
!
end subroutine elem_heat

