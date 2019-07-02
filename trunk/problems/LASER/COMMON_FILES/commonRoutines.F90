!----------------------------------------------------------------------
!
!   routine name       - propagate_flag
!
!----------------------------------------------------------------------
!
!   latest revision    - Sep 17
!
!   purpose            - propagate Nflag from element faces to element
!                        edges and vertices; the flag is passed
!                        provided ALL adjacent faces share the flag
!
!   arguments:
!     in :
!              Icomp   - physics attribute number
!              Nflag   - BC flag
!
!----------------------------------------------------------------------
!
#include "implicit_none.h"
!
subroutine propagate_flag(Icomp,Nflag)
!
      use data_structure3D
      implicit none
      integer :: Icomp,Nflag
!
      character(len=4) :: type
      integer :: iel,mdle,if,nrfn,i,j,nod,iprint
!
!  ...element nodes and orientations, face nodes
      integer :: nodesl(27),norientl(27),nface_nodes(9)
!
!  ...element face BC flags, decoded BC flag for a node
      integer :: ibc(6,NR_PHYSA), nodflag(NR_PHYSA)
!
      iprint=0
!
!  ...loop through active elements
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        type = NODES(mdle)%type
!
!  .....determine element nodes
        call elem_nodes(mdle, nodesl,norientl)
!
!  .....get the element boundary conditions flags
        call find_bc(mdle, ibc)
!
!  .....loop through element faces
        do if=1,nface(type)
!
!  .......determine face node numbers
          call face_nodes(type,if, nface_nodes,nrfn)
          if (iprint.eq.1) then
            write(*,7010) mdle,if, ibc(if,Icomp)
 7010       format('propagate_flag: mdle,if,ibc(if,Icomp) = ',i6,i3,i4)
          endif
!
!  .......loop through the face nodes
          do i=1,nrfn-1
            j = nface_nodes(i)
            nod = nodesl(j)
!
!  .........propagate the flag unless prohibited
            if (ibc(if,Icomp).eq.Nflag) then
              if (NODES(nod)%visit.ne.-Nflag) NODES(nod)%visit = Nflag
!
!  .........prohibit the flag to be passed to the node
            else
              NODES(nod)%visit = -Nflag
            endif
          enddo
        enddo
      enddo
!
!  ...change -Nflag to zero
      do nod=1,NRNODS
        if (iprint.eq.1) then
          write(*,7020) nod,NODES(nod)%visit
 7020     format('propagate_flag: nod,NODES(nod)%visit = ',i8,i3)
        endif
        if (NODES(nod)%visit.eq.0) cycle
        call decod(NODES(nod)%bcond,10,NR_PHYSA, nodflag)
        if (NODES(nod)%visit.eq.-Nflag) then
          nodflag(Icomp)=0
        elseif (NODES(nod)%visit.eq.Nflag) then
          nodflag(Icomp)=Nflag
        endif
        call encod(nodflag,10,NR_PHYSA, NODES(nod)%bcond)
!
!  .....reset the node index
        call set_index(NODES(nod)%case,NODES(nod)%bcond, NODES(nod)%index)
      enddo
      call reset_visit
!
end subroutine propagate_flag
!
!
!----------------------------------------------------------------------
!
!   routine name       - get_power
!
!----------------------------------------------------------------------
!
!   latest revision    - Jan 18
!
!   purpose            - Driver routine for computing power in UW
!                        Maxwell, i.e. the Poynting vector at certain
!                        z-points for both pump and signal.
!                        Z-points are samples in this routine.
!
!   arguments          - none
!
!----------------------------------------------------------------------
!
subroutine get_power
!
   use CommonParam
   use LaserParam
!
   implicit none
!
   integer                           :: fld_flag
   real*8, allocatable, dimension(:) :: zValues
   real*8, allocatable, dimension(:) :: signal_power, pump_power
   real*8, allocatable, dimension(:) :: diff_power, efficiency
   integer                           :: zlength, numPts
!
   real*8  :: a,b
   integer :: i
!
!..Initial value
   a = 0.05d0
   b = 0.005d0
!
   numPts = 200
   if(ZL.gt.1.d0) then
      zlength = numPts*(int(ZL)+1)
   else
      zlength = numPts
   endif
!
   allocate(zValues(zlength), signal_power(zlength), &
            pump_power(zlength), diff_power(zlength), &
            efficiency(zlength))
!
!..get z-values (sample points)
   zValues = (/(((i-1)*b+a),i=1,zlength)/)
!
!..make z-values irrational numbers s.t. they are not on element boundaries
   zValues = zValues*PI*(7.d0/22.d0)
   write(*,*) ' from get_power: zValues'
   do i =1,zlength
      write(*,*) zValues(i)
   enddo
!
!..get signal power (fld_flag=1 means signal, 0 means pump)
   fld_flag = 1
   call compute_power(zValues,zlength,fld_flag, signal_power,diff_power)
   write(*,*) ' from get_power: signal_power'
   do i = 1,zlength
      write(*,*) signal_power(i)
   enddo
!
!..get pump power
   fld_flag = 0
   call compute_power(zValues,zlength,fld_flag, pump_power,diff_power)
   write(*,*) ' from get_power: pump_power'
   do i = 1,zlength
      write(*,*) pump_power(i)
   enddo
!
!..get efficiency
   efficiency(1) = 0.d0
   do i = 2,zlength
      if(COPUMP.eq.1) then
         efficiency(i) = (signal_power(i)-signal_power(1))&
                          /(pump_power(1)-pump_power(i))
      elseif(COPUMP.eq.0) then
         efficiency(i) = (signal_power(i)-signal_power(1))&
                          /(pump_power(zlength)-pump_power(i))
      else
         write(*,*) 'from get_power: COPUMP must be 1 or 0. stop.'
         stop
    endif
   enddo
   write(*,*) ' from get_power: efficiency'
   do i = 1,zlength
      write(*,*) efficiency(i)
   enddo
   deallocate(zValues,signal_power,pump_power,diff_power,efficiency)
!
end subroutine get_power
!
!
!----------------------------------------------------------------------
!
!   routine name       - compute_power
!
!----------------------------------------------------------------------
!
!   latest revision    - Jan 18
!
!   purpose            - Evaluates the electric field power of UW
!                        Maxwell along the cross sections specified by
!                        the vector of zValues in the input
!
!   arguments
!        in:
!                      - zValues     : sample points in z-direction
!                      - length_zpts : number of sample points
!                      - fld_flag    : 1 (signal) or 0 (pump)
!       out:
!                      - power       :
!                      - diff_power  :
!
!----------------------------------------------------------------------
!
subroutine compute_power(zValues,length_zpts,fld_flag, power,diff_power)
!
   use data_structure3D
   use environment, only : QUIET_MODE
!
   implicit none
!
   integer                       , intent(in)  :: length_zpts
   real*8, dimension(length_zpts), intent(in)  :: zValues
   integer                       , intent(in)  :: fld_flag
   real*8, dimension(length_zpts), intent(out) :: power
   real*8, dimension(length_zpts), intent(out) :: diff_power
!
!..auxiliary variables
   real*8 :: facePower, face_diff_power, elemPower
!
!..mdle number
   integer                    :: mdle
   integer, dimension(NRELES) :: mdlea
!
!..element, face order, geometry dof
   real*8 ,dimension(3,MAXbrickH) :: xnod
!
!..number of vertices,edge,faces per element type
   integer :: nrv, nre, nrf
!
!..miscellanea
   integer :: iel, i, iprint
!
!..auxiliary variables for timing
   real*8 :: start, OMP_get_wtime
!
!---------------------------------------------------------------------------------------
!
   iprint=0
!
!..initialize outputs (vector of powers for all z-points)
   power = 0.d0
   diff_power = 0.d0
!
!..initialize running powers computed (elements per z-point)
   facePower = 0.d0
   face_diff_power  = 0.d0
!
!..start timer
   start = OMP_get_wtime()
!
   mdle=0
   do iel=1,NRELES
      call nelcon(mdle, mdle)
      mdlea(iel) = mdle
   enddo
!
!..iterate over elements
! TODO: very inefficient computation currently
!       reduce number of sample points
!       sample once per element (at most)
!
!$OMP PARALLEL DO                                     &
!$OMP PRIVATE(mdle,xnod,i,facePower,face_diff_power)  &
!$OMP REDUCTION(+:power,diff_power)                   &
!$OMP SCHEDULE(DYNAMIC)
   do iel=1,NRELES
      mdle = mdlea(iel)
      call nodcor(mdle, xnod)
      do i=1,length_zpts
         if((zValues(i).lt.xnod(3,8)).and.(zValues(i).gt.xnod(3,1))) then
            call compute_facePower(mdle,2,fld_flag,zValues(i), facePower,face_diff_power)
            power(i) = power(i)+facePower
            diff_power(i) = diff_power(i) + face_diff_power
         endif
      enddo
   enddo
!$OMP END PARALLEL DO
!
!..take absolute value after integration
   do i=1,length_zpts
      power(i) = abs(power(i))
      diff_power(i) = abs(diff_power(i))
   enddo
!
!..end timer
   if (.not. QUIET_MODE) then
      write(*,10) OMP_get_wtime()-start
 10   format(/,' compute_power  : ',f12.5,'  seconds',/)
   endif
!
end subroutine compute_power
!
!
!----------------------------------------------------------------------
!
!   routine name       - compute_face_power
!
!----------------------------------------------------------------------
!
!   latest revision    - Jan 18
!
!   purpose            - Evaluates the electric field power of UW
!                        Maxweel by integrating H(curl) trace solution
!                        on a face of a middle node.
!
!   arguments
!        in:
!                      - mdle       : middle element node
!                      - facenumber :
!                      - fld_flag   : 1 (signal) or 0 (pump)
!                      - zpoint     : value of z at cross section face
!       out:
!                      - facePower       :
!                      - face_diff_power :
!
!----------------------------------------------------------------------
!
subroutine compute_facePower(mdle,facenumber,fld_flag,zpoint, facePower,face_diff_power)
!
   use control
   use data_structure3D
   use environment      , only : L2PROJ
   use physics
   use parametersDPG
   use CommonParam
!
   implicit none
!
   integer, intent(in)  :: mdle
   integer, intent(in)  :: fld_flag
   integer, intent(in)  :: facenumber
   real*8,  intent(in)  :: zpoint
   real*8,  intent(out) :: facePower
   real*8,  intent(out) :: face_diff_power
!
!..element, face order, geometry dof
   integer,dimension(19)          :: norder
   real*8 ,dimension(3,MAXbrickH) :: xnod
   integer,dimension(12)          :: nedge_orient
   integer,dimension(6)           :: nface_orient
!
!..face order
   integer, dimension(5) :: norderf
!
!..number of vertices,edge,faces per element type
   integer :: nrv, nre, nrf
!
!..declare edge/face type varibles
   character(len=4) :: etype,ftype
!
!..variables for geometry
   real*8, dimension(3)      :: xi,x,rn,x_new
   real*8, dimension(3,2)    :: dxidt,dxdt,rt
   real*8, dimension(3,3)    :: dxdxi,dxidx
   real*8, dimension(2)      :: t
   real*8                    :: rjac,bjac
!
!..2D quadrature data
   real*8, dimension(2,MAXNINT2ADD)  :: tloc
   real*8, dimension(MAXNINT2ADD)    :: wtloc
!
!..approximate solution dof's
   VTYPE, dimension(MAXEQNH,MAXbrickH) :: zdofH
   VTYPE, dimension(MAXEQNE,MAXbrickE) :: zdofE
   VTYPE, dimension(MAXEQNV,MAXbrickV) :: zdofV
   VTYPE, dimension(MAXEQNQ,MAXbrickQ) :: zdofQ
!..H1 shape functions
   integer                         :: nrdofH
   real*8, dimension(MAXbrickH)    :: shapH
   real*8, dimension(3,MAXbrickH)  :: gradH
!
!..approximate solution
   VTYPE, dimension(  MAXEQNH  ) ::  zsolH
   VTYPE, dimension(  MAXEQNH,3) :: zdsolH
   VTYPE, dimension(3,MAXEQNE  ) ::  zsolE
   VTYPE, dimension(3,MAXEQNE  ) :: zcurlE
   VTYPE, dimension(3,MAXEQNV  ) ::  zsolV
   VTYPE, dimension(  MAXEQNV  ) ::  zdivV
   VTYPE, dimension(  MAXEQNQ  ) ::  zsolQ
!
!..exact solution
   VTYPE,dimension(  MAXEQNH    )  ::   ValH
   VTYPE,dimension(  MAXEQNH,3  )  ::  DvalH
   VTYPE,dimension(  MAXEQNH,3,3)  :: d2valH
   VTYPE,dimension(3,MAXEQNE    )  ::   ValE
   VTYPE,dimension(3,MAXEQNE,3  )  ::  DvalE
   VTYPE,dimension(3,MAXEQNE,3,3)  :: d2valE
   VTYPE,dimension(3,MAXEQNV    )  ::   ValV
   VTYPE,dimension(3,MAXEQNV,3  )  ::  DvalV
!
!..exact solution (UNUSED)
   VTYPE,dimension(3,MAXEQNV,3,3)  :: d2valV
   VTYPE,dimension(  MAXEQNQ    )  ::   valQ
   VTYPE,dimension(  MAXEQNQ,3  )  ::  dvalQ
   VTYPE,dimension(  MAXEQNQ,3,3)  :: d2valQ
!..for Poynting vector
   VTYPE, dimension(3)          :: EtimesH1,EtimesH2
   VTYPE                        :: FdotN
!
!..miscellanea
   integer :: nint,icase,iattr,l,i,j
   real*8  :: weight,wa
   integer :: iel,nsign
   integer :: iprint,nflag,iload
!
!---------------------------------------------------------------------------------------
!
   iprint=0
!
   facePower = 0.d0
   face_diff_power = 0.0d0
   nflag = 1
!..element type
   etype = NODES(mdle)%type
   nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
   call find_order(mdle, norder)
   call find_orient(mdle, nedge_orient,nface_orient)
   call nodcor(mdle, xnod)
   call solelm(mdle, zdofH,zdofE,zdofV,zdofQ)
!..sign factor to determine the OUTWARD normal unit vector
   nsign = nsign_param(etype,facenumber)
!
!..face type
   ftype = face_type(etype,facenumber)
!
!..face order of approximation
   call face_order(etype,facenumber,norder, norderf)
!
!..set 2D quadrature
   INTEGRATION = NORD_ADD
   call set_2Dint(ftype,norderf, nint,tloc,wtloc)
   INTEGRATION = 0
!
!..loop over integration points
   do l=1,nint
!
!  ...face coordinates
      t(1:2) = tloc(1:2,l)
!
!  ...face parametrization
      call face_param(etype,facenumber,t, xi,dxidt)
!
!  ...determine element H1 shape functions (for geometry)
      call shape3H(etype,xi,norder,nedge_orient,nface_orient, &
                     nrdofH,shapH,gradH)
!
!  ...geometry
      call bgeom3D(mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign, &
                     x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
      weight = bjac*wtloc(l)
!
      call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod,zdofH,zdofE,zdofV,zdofQ,nflag, &
                         x,dxdxi,zsolH,zdsolH,zsolE,zcurlE,zsolV,zdivV,zsolQ)
      if(NEXACT.eq.1) then
         call exact(x,Mdle, ValH,DvalH,d2valH, ValE,DvalE,d2valE, &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      endif
!
!     accumulate Poynting vector power for signal (fld_flag=1) or pump (fld_flag=0),
!     i.e., integrate (Real(n \dot ExH^*)) with:
!                     E/H corresponding to signal if fld_flag = 1
!                     E/H corresponding to pump if fld_flag = 0
!  ...first check for signal, i.e, if fld_flag = 1
      if(fld_flag.eq.1) then
         call zz_cross_product(zsolE(1:3,1),conjg((zsolE(1:3,2))), EtimesH1)
         FdotN = EtimesH1(1)*rn(1)+EtimesH1(2)*rn(2)+EtimesH1(3)*rn(3)
         facePower = facePower + (real(FdotN))*weight
!     ...if we have an exact
         if(NEXACT.eq.1) then
            call zz_cross_product(valE(1:3,1),conjg((valE(1:3,2))), EtimesH2)
            face_diff_power = face_diff_power   &
                           + abs(((EtimesH1(1)*rn(1)+EtimesH1(2)*rn(2)+EtimesH1(3)*rn(3))*weight) - &
                           ((EtimesH2(1)*rn(1)+EtimesH2(2)*rn(2)+EtimesH2(3)*rn(3))*weight))
         endif
!  ...next check for pump, i.e, if fld_flag = 0
      else if(fld_flag.eq.0) then
         call zz_cross_product(zsolE(1:3,3),conjg((zsolE(1:3,4))), EtimesH1)
         FdotN = EtimesH1(1)*rn(1)+EtimesH1(2)*rn(2)+EtimesH1(3)*rn(3)
         facePower = facePower + (real(FdotN))*weight
!     ...if we have an exact
         if(NEXACT.eq.1) then
            call zz_cross_product(valE(1:3,3),conjg((valE(1:3,4))), EtimesH2)
            face_diff_power = face_diff_power   &
                            + abs(((EtimesH1(1)*rn(1)+EtimesH1(2)*rn(2)+EtimesH1(3)*rn(3))*weight) - &
                                 ((EtimesH2(1)*rn(1)+EtimesH2(2)*rn(2)+EtimesH2(3)*rn(3))*weight))
         endif
      else
         write(*,*) 'compute_facePower: fld_flag must be 0 or 1. stop.'
         stop
      endif
!..end loop over integration points
   enddo
!
end subroutine compute_facePower
!
!
!-------------------------------------------------------------------
!Routine to evaluate the L2 norm of difference of component No1 and
! No2 (specified in NRCOMS) of L2 variable
!Last modified : December 17
!input: Flag,No1,No2
!output: L2NormDiff
!-------------------------------------------------------------------
!
subroutine get_L2NormCOMS(Flag,No1,No2, L2NormDiff)
!
      use control          , only : NEXACT
      use data_structure3D
      use environment      , only : QUIET_MODE,L2PROJ,FILE_ERR
      use physics
!
      implicit none
      integer, dimension(NR_PHYSA),intent(in) :: Flag
      integer, intent(in) :: No1,No2
!
      real*8,  intent(out) :: L2NormDiff
      real*8               :: CurrL2NormDiff
!
      integer, parameter :: nin = 13
      integer, parameter :: maxvis =2000
!
!
!     miscellanea
      integer :: mdle,i,nint,iattr,nrdof_tot,ic
!
!     printing flag
      integer :: iprint
!
!---------------------------------------------------------------------------------------
!
      iprint=0
!
!     initialize global quantities
      L2NormDiff=0.d0
      CurrL2NormDiff = 0.d0
!
!     loop over active elements
      mdle=0
      do i=1,NRELES
        call nelcon(mdle,mdle)
        call get_elem_L2NormCOMS(mdle,Flag,No1,No2, CurrL2NormDiff)
!   ... accumulate L2NormDiff
        L2NormDiff = L2NormDiff + CurrL2NormDiff
      enddo
!     compute sqrt of total difference
!
      L2NormDiff   =sqrt(L2NormDiff)
end subroutine get_L2NormCOMS



!-------------------------------------------------------------------
!Routine to evaluate the L2 norm of difference between components
!No1 and No2 with components (set using NRCOMS) of L2 variables
!Last modified : December 17
!input: mdle,Flag,No1,No2
!output: L2 norm of difference (FieldNormQ)
!-------------------------------------------------------------------
subroutine get_elem_L2NormCOMS(mdle,Flag,No1,No2, FieldNormQ)
      use LaserParam
      use CommonParam
    !
      use control          , only : INTEGRATION
      use data_structure3D
      use environment      , only : L2PROJ
      use physics
      implicit none

      integer, dimension(NR_PHYSA),intent(in ) :: Flag
      integer,                     intent(in ) :: mdle,No1,No2
      real*8,                      intent(out) :: FieldNormQ
!
!     node case (decimal form)
      integer,dimension(NR_PHYSA) :: icased
!
!     element, face order, geometry dof
      integer,dimension(19)          :: norder
      real*8 ,dimension(3,MAXbrickH) :: xnod
      integer,dimension(12)          :: nedge_orient
      integer,dimension(6)           :: nface_orient
!
!     geometry
      real*8,dimension(3)   :: xi,x
      real*8,dimension(3,3) :: dxidx,dxdxi
      real*8                :: rjac
!
!     3D quadrature data
      real*8,dimension(3,MAX_NINT3) :: xiloc
      real*8,dimension(  MAX_NINT3) :: wxi
!
!     approximate solution dof's
      VTYPE, dimension(MAXEQNH,MAXbrickH) :: zdofH
      VTYPE, dimension(MAXEQNE,MAXbrickE) :: zdofE
      VTYPE, dimension(MAXEQNV,MAXbrickV) :: zdofV
      VTYPE, dimension(MAXEQNQ,MAXbrickQ) :: zdofQ
!
!     approximate solution
      VTYPE, dimension(  MAXEQNH  ) ::  zsolH
      VTYPE, dimension(  MAXEQNH,3) :: zdsolH
      VTYPE, dimension(3,MAXEQNE  ) ::  zsolE
      VTYPE, dimension(3,MAXEQNE  ) :: zcurlE
      VTYPE, dimension(3,MAXEQNV  ) ::  zsolV
      VTYPE, dimension(  MAXEQNV  ) ::  zdivV
      VTYPE, dimension(  MAXEQNQ  ) ::  zsolQ
!
!     miscellanea
      integer :: nint,icase,iattr,l,i,j,ibeg,iflag,iload,icomp,ndom,ivar,nflag
      real*8  :: weight,wa
!
!     printing flag
      integer :: iprint
!
!---------------------------------------------------------------------------------------
!
      iprint=0
!
!     initialize global quantities

      FieldNormQ=0.d0
!
!     order of approx, orientations, geometry dof's, solution dof's
      call find_order( mdle, norder)
      call find_orient(mdle, nedge_orient,nface_orient)
      call nodcor(     mdle, xnod)
      call solelm(     mdle, zdofH,zdofE,zdofV,zdofQ)
!
!     set up the element quadrature
      INTEGRATION=0
      call set_3Dint_DPG(NODES(mdle)%type,norder, nint,xiloc,wxi)
      INTEGRATION=0
!
!     supported physical attributes
      icase=NODES(mdle)%case
      call decod(icase,2,NR_PHYSA, icased)
!
!     loop over physical attributes
      do iattr=1,NR_PHYSA
!
!       if the error not needed, skip
        if (Flag(iattr) == 0) cycle
!
!       if attribute is absent, skip
        if (icased(iattr) == 0) cycle
!
!       address of the 1st component for the attribute
        ibeg=ADRES(iattr)
!
        select case(DTYPE(iattr))
!===================================================================================
!  L2 ATTRIBUTE                                                                    |
!===================================================================================
          case('discon')
!
!         loop through integration points
          do l=1,nint
!
!           Gauss point and weight
            xi(1:3)=xiloc(1:3,l) ; wa=wxi(l)
!
!           -- APPROXIMATE SOLUTION --
            nflag=1
            call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod,zdofH,zdofE,zdofV,zdofQ,nflag, &
                         x,dxdxi,zsolH,zdsolH,zsolE,zcurlE,zsolV,zdivV,zsolQ                          )
!
!           Jacobian
            call geom(dxdxi, dxidx,rjac,iflag)
            if (iflag /= 0) then
              call find_domain(mdle, ndom)
              write(*,*)'from get_elem_L2NormCOMS: mdle,ndom,rjac = ', mdle,ndom,rjac
            endif

!           total weight
            weight=wa*rjac
!
!             loop over components of the physical attribute
            do icomp=1,NR_COMP(iattr)
!
              i=(No1-1)*NRQVAR+ibeg+icomp
              j=(No2-1)*NRQVAR+ibeg+icomp
!
!               accumulate L2 norm
              FieldNormQ = FieldNormQ+(abs(zsolQ(i)-zsolQ(j))**2)*weight
!
            enddo

!         loop over integration points
          enddo
!
          case default
          write(*,*)'get_elem_L2NormCOMS: PHYSICAL ATTRIBUTE MUST BE L2!'
          stop
!
          endselect
!
!       loop over physical attributes
        enddo
end subroutine get_elem_L2NormCOMS

!-------------------------------------------------------------------
!Routine to evaluate the Norm of specified physical attributes
!and specified component
!by integrating the current solution over the entire domain
!Last modified : September 17
!input: Field, No
!output: FieldNorm (norm of field variable) = \int_{\Omega} |X|^2
!        FieldNormH,FieldNormE,FieldNormV,FieldNormQ
!-------------------------------------------------------------------
!
subroutine get_Norm(Flag,No, FieldNormH,FieldNormE,FieldNormV,FieldNormQ)
!
      use control          , only : NEXACT
      use data_structure3D
      use environment      , only : QUIET_MODE,L2PROJ,FILE_ERR
      use physics
!
      implicit none
      integer, dimension(NR_PHYSA),intent(in) :: Flag
      integer, intent(in)					  :: No
!
      real*8,  intent(out) :: FieldNormH,FieldNormE,FieldNormV,FieldNormQ
      real*8               :: CurrFieldNormH,CurrFieldNormE,CurrFieldNormV,CurrFieldNormQ
!
      integer, parameter :: nin = 13
      integer, parameter :: maxvis =2000
!
!
!     miscellanea
      integer :: mdle,i,nint,iattr,nrdof_tot,ic
!
!     printing flag
      integer :: iprint
!
!---------------------------------------------------------------------------------------
!
      iprint=0
!
!
!     initialize global quantities
      FieldNormH=0.d0
      FieldNormE=0.d0
      FieldNormV=0.d0
      FieldNormQ=0.d0

      CurrFieldNormH = 0.d0
      CurrFieldNormE = 0.d0
      CurrFieldNormV = 0.d0
      CurrFieldNormQ = 0.d0


!
!     loop over active elements
      mdle=0
      do i=1,NRELES
        call nelcon(mdle,mdle)
        call get_elem_Norm(mdle,Flag,No, CurrFieldNormH,CurrFieldNormE,CurrFieldNormV,CurrFieldNormQ)
!
!       accumulate
        FieldNormH = FieldNormH + CurrFieldNormH
        FieldNormE = FieldNormE + CurrFieldNormE
        FieldNormV = FieldNormV + CurrFieldNormV
        FieldNormQ = FieldNormQ + CurrFieldNormQ
      enddo
!
!     compute total error
      FieldNormH   =sqrt(FieldNormH)
      FieldNormE   =sqrt(FieldNormE)
      FieldNormV   =sqrt(FieldNormV)
      FieldNormQ   =sqrt(FieldNormQ)
!
end subroutine get_Norm
!
!-------------------------------------------------------------------
!Routine to evaluate the magnitude of a physical attribute
!and component by integrating the current solution over a given middle node
!Last modified : September 17
!input: Mdle,Flag,No
!output: FieldNorm (norm of field variable) consisting of:
!        FieldNormH,FieldNormE,FieldNormV,FieldNormQ
!-------------------------------------------------------------------
subroutine get_elem_Norm(Mdle,Flag,No, FieldNormH,FieldNormE,FieldNormV,FieldNormQ)
      use LaserParam
      use CommonParam
    !
      use control          , only : INTEGRATION
      use data_structure3D
      use environment      , only : L2PROJ
      use physics
      implicit none

      integer, dimension(NR_PHYSA),intent(in ) :: Flag
      integer,                     intent(in ) :: Mdle,No
      real*8,                      intent(out) :: FieldNormH,FieldNormE,FieldNormV,FieldNormQ
!
!     node case (decimal form)
      integer,dimension(NR_PHYSA) :: icased
!
!     element, face order, geometry dof
      integer,dimension(19)          :: norder
      real*8 ,dimension(3,MAXbrickH) :: xnod
      integer,dimension(12)          :: nedge_orient
      integer,dimension(6)           :: nface_orient
!
!     geometry
      real*8,dimension(3)   :: xi,x
      real*8,dimension(3,3) :: dxidx,dxdxi
      real*8                :: rjac
!
!     3D quadrature data
      real*8,dimension(3,MAX_NINT3) :: xiloc
      real*8,dimension(  MAX_NINT3) :: wxi
!
!     approximate solution dof's
      VTYPE, dimension(MAXEQNH,MAXbrickH) :: zdofH
      VTYPE, dimension(MAXEQNE,MAXbrickE) :: zdofE
      VTYPE, dimension(MAXEQNV,MAXbrickV) :: zdofV
      VTYPE, dimension(MAXEQNQ,MAXbrickQ) :: zdofQ
!
!     approximate solution
      VTYPE, dimension(  MAXEQNH  ) ::  zsolH
      VTYPE, dimension(  MAXEQNH,3) :: zdsolH
      VTYPE, dimension(3,MAXEQNE  ) ::  zsolE
      VTYPE, dimension(3,MAXEQNE  ) :: zcurlE
      VTYPE, dimension(3,MAXEQNV  ) ::  zsolV
      VTYPE, dimension(  MAXEQNV  ) ::  zdivV
      VTYPE, dimension(  MAXEQNQ  ) ::  zsolQ
!
!     miscellanea
      integer :: nint,icase,iattr,l,i,j,ibeg,iflag,iload,icomp,ndom,ivar,nflag
      real*8  :: weight,wa
!
!     printing flag
      integer :: iprint
!
!---------------------------------------------------------------------------------------
!
      iprint=0
!
!     initialize global quantities
      FieldNormH=0.d0
      FieldNormE=0.d0
      FieldNormV=0.d0
      FieldNormQ=0.d0
!
!     order of approx, orientations, geometry dof's, solution dof's
      call find_order( mdle, norder)
      call find_orient(mdle, nedge_orient,nface_orient)
      call nodcor(     mdle, xnod)
      call solelm(     mdle, zdofH,zdofE,zdofV,zdofQ)
!
!     set up the element quadrature
      INTEGRATION=0
      call set_3Dint_DPG(NODES(mdle)%type,norder, nint,xiloc,wxi)
      INTEGRATION=0
!
!     supported physical attributes
      icase=NODES(mdle)%case
      call decod(icase,2,NR_PHYSA, icased)
!
!     loop over physical attributes
      do iattr=1,NR_PHYSA
!
!       if the error not needed, skip
        if (Flag(iattr) == 0) cycle
!
!       if attribute is absent, skip
        if (icased(iattr) == 0) cycle
!
!       address of the 1st component for the attribute
        ibeg=ADRES(iattr)
!
        select case(DTYPE(iattr))
!
!===================================================================================
!  H1 ATTRIBUTE                                                                    |
!===================================================================================
          case('contin')
!
!         loop through integration points
          do l=1,nint
!
!           Gauss point and weight
            xi(1:3)=xiloc(1:3,l) ; wa=wxi(l)
!
!           -- APPROXIMATE SOLUTION --
            nflag=1
            call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod,zdofH,zdofE,zdofV,zdofQ,nflag, &
                         x,dxdxi,zsolH,zdsolH,zsolE,zcurlE,zsolV,zdivV,zsolQ                          )
!           Jacobian
            call geom(dxdxi, dxidx,rjac,iflag)
            if (iflag /= 0) then
              call find_domain(mdle, ndom)
              write(*,9997) mdle,ndom,rjac
9997          format(' element_error: mdle,ndom,rjac = ',i8,2x,i2,2x,e12.5)
            endif
!
!           total weight
            weight=wa*rjac
!
!
!             loop over components of physical attribute
              do icomp=1,NR_COMP(iattr)
!
                i=(No-1)*NRHVAR+ibeg+icomp
!
!               accumulate H1 seminorm
                IF (.NOT. L2PROJ) THEN
                do j=1,3
                  FieldNormH = FieldNormH + abs(zdsolH(i,j)              )**2 * weight
                enddo
                ENDIF
!
!               accumulate L2 norm
                FieldNormH = FieldNormH + abs(zsolH(i)           )**2 * weight
!
!             loop over components
              enddo
!
!         loop over integration points
          enddo
!
!===================================================================================
!  H(curl) ATTRIBUTE                                                               |
!===================================================================================
          case('tangen')
!
!         loop through integration points
          do l=1,nint
!
!           Gauss point and weight
            xi(1:3)=xiloc(1:3,l) ; wa=wxi(l)
!
!           -- APPROXIMATE SOLUTION --
            nflag=1
            call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod,zdofH,zdofE,zdofV,zdofQ,nflag, &
                         x,dxdxi,zsolH,zdsolH,zsolE,zcurlE,zsolV,zdivV,zsolQ                          )
!
!           Jacobian
            call geom(dxdxi, dxidx,rjac,iflag)
            if (iflag /= 0) then
              call find_domain(mdle, ndom)
              write(*,9997) mdle,ndom,rjac
            endif
!
!           total weight
            weight=wa*rjac
!
!
!             loop over the components of the physical attribute
              do icomp=1,NR_COMP(iattr)
!
                i=(No-1)*NREVAR+ibeg+icomp
!
!               accumulate for the error and norm
                do ivar=1,3
                  FieldNormE = FieldNormE + abs(zsolE(ivar,i))**2 * weight
                  IF (.NOT. L2PROJ) THEN
                  FieldNormE = FieldNormE + abs(zcurlE(ivar,i))**2 * weight
                  ENDIF
                enddo
!
!             loop over components
              enddo
!
!         loop over integration points
          enddo
!
!===================================================================================
!  H(div) ATTRIBUTE                                                                |
!===================================================================================
          case('normal')
!
!         loop over integration points
          do l=1,nint
!
!           Gauss point and weight
            xi(1:3)=xiloc(1:3,l) ; wa=wxi(l)
!
!           -- APPROXIMATE SOLUTION --
            nflag=1
            call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod,zdofH,zdofE,zdofV,zdofQ,nflag, &
                         x,dxdxi,zsolH,zdsolH,zsolE,zcurlE,zsolV,zdivV,zsolQ                          )

            call geom(dxdxi, dxidx,rjac,iflag)
            if (iflag /= 0) then
              call find_domain(mdle, ndom)
              write(*,9997) mdle,ndom,rjac
            endif
!
!           total weight
            weight=wa*rjac
!
!
!             loop over components of the physical attribute
              do icomp=1,NR_COMP(iattr)
!
                i=(No-1)*NRVVAR+ibeg+icomp
!
!               accumulate H(div) seminorm
                IF (.NOT. L2PROJ) THEN
                FieldNormV = FieldNormV + abs(zdivV(i))**2 * weight
                ENDIF
!
!               accumulate L2 norm
                do ivar=1,3
                  FieldNormV = FieldNormV + abs(zsolV(ivar,i))**2 * weight
                enddo
!
!             end loop over components
              enddo
!
!         loop over integration points
          enddo
!
!===================================================================================
!  L2 ATTRIBUTE                                                                    |
!===================================================================================
          case('discon')
!
!         loop through integration points
          do l=1,nint
!
!           Gauss point and weight
            xi(1:3)=xiloc(1:3,l) ; wa=wxi(l)
!
!           -- APPROXIMATE SOLUTION --
            nflag=1
            call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod,zdofH,zdofE,zdofV,zdofQ,nflag, &
                         x,dxdxi,zsolH,zdsolH,zsolE,zcurlE,zsolV,zdivV,zsolQ                          )
!
!           Jacobian
            call geom(dxdxi, dxidx,rjac,iflag)
            if (iflag /= 0) then
              call find_domain(mdle, ndom)
              write(*,9997) mdle,ndom,rjac
            endif

!           total weight
            weight=wa*rjac
!
!
!             loop over components of the physical attribute
              do icomp=1,NR_COMP(iattr)
!
                i=(No-1)*NRQVAR+ibeg+icomp
!
!               accumulate L2 norm
                FieldNormQ = FieldNormQ + abs(zsolQ(i))**2 * weight
!
!             end loop over components
              enddo
!
!         loop over integration points
          enddo
!
          case default
          write(*,*)'element_error: UNKNOWN PHYSICAL ATTRIBUTE TYPE!'
          stop
!
          endselect
!
!       loop over physical attributes
        enddo
end subroutine get_elem_Norm
!
!----------------------------------------------------------------------
!
!   routine name       - setAnisoRef
!
!----------------------------------------------------------------------
!
!   latest revision    - Sep 17
!
!   purpose            -  This routine anisotropically refines prism and hexa elements
!   inputs             -
!                        ref_xyz: refinement flag to set kref
!                         1 for xy direction
!                         2 for z direction
!
!
!----------------------------------------------------------------------
subroutine setAnisoRef(ref_xyz)
!
    use control
    use data_structure3D
    use parameters, only : ZERO,ZONE,ZIMG
      implicit none
      integer,                       intent(in)  :: ref_xyz
      integer, allocatable, dimension(:) :: list_elem
      integer :: i,iprint,ic,mdle,iel,kref,nr_elem_to_refine
      allocate (list_elem(NRELES))
      write(*,*) 'NRELES is: ', NRELES
        ic=0
        mdle=0
        do iel=1,NRELES
            call nelcon(mdle, mdle)
            ic=ic+1
            list_elem(ic) = mdle
        enddo
        nr_elem_to_refine = NRELES
        if (nr_elem_to_refine.gt.0) then
            !      ...refine the elements from the list
            do iel=1,nr_elem_to_refine
                mdle = list_elem(iel)
                select case(NODES(mdle)%type)
!               PRISM
                case('mdlp')
                  if(ref_xyz.eq.1) then
                    kref=10
                  elseif(ref_xyz.eq.2) then
                    kref=1
                  else
                    write(*,*) 'error from IANISOREF: invalid ref_xyz in prism'
                    stop
                  endif
!
!               BRICK
                case('mdlb')
                  if(ref_xyz.eq.1) then
                    kref=110
                  elseif(ref_xyz.eq.2) then
                    kref=1
                  else
                    write(*,*) 'error from IANISOREF: invalid ref_xyz in bric'
                    stop
                  endif
                end select
                call refine(mdle,kref)
            enddo
            !      ...close the mesh
            call close
        endif
   end subroutine setAnisoRef

!  just to display the current size of the data structure NODES (in bytes)
   subroutine my_sizetest
        use data_structure3D
        implicit none
        integer :: nod, nH,nE,nV,nQ
        integer*8 :: size_coord, size_h1, size_hcurl, size_hdiv, size_l2, size_tot
        size_coord = 0
        size_h1 = 0
        size_hcurl = 0
        size_hdiv = 0
        size_l2 = 0

         do nod=1,NRNODS
             if (associated(NODES(nod)%coord).eq. .true.) then
               size_coord = size_coord + STORAGE_SIZE(NODES(nod)%coord)
             endif
             if (associated(NODES(nod)%zdofH).eq. .true.) then
               size_h1 = size_h1 + STORAGE_SIZE(NODES(nod)%zdofH)
             endif
             if (associated(NODES(nod)%zdofE).eq. .true.) then
               size_hcurl = size_hcurl + STORAGE_SIZE(NODES(nod)%zdofE)
             endif
             if (associated(NODES(nod)%zdofV).eq. .true.) then
               size_hdiv = size_hdiv + STORAGE_SIZE(NODES(nod)%zdofV)
             endif
             if (associated(NODES(nod)%zdofQ).eq. .true.) then
               size_l2 = size_l2 + STORAGE_SIZE(NODES(nod)%zdofQ)
             endif
         size_tot = size_coord + size_h1 + size_hcurl + size_hdiv + size_l2
        enddo

         write(*,*) '1: total DOF size is: ', size_tot
         ! call pause

         nH=0; nE=0; nV=0; nQ=0

         write(*,*) 'my_sizetest: NRHVAR, NREVAR, NRVVAR, NRQVAR = ', &
                                  NRHVAR, NREVAR, NRVVAR, NRQVAR
        do nod = 1, NRNODS
          if (NODES(nod)%act .eq. 0) cycle
          select case(NODES(nod)%type)
          case('vert')
            nH = nH + 1
          case('medg')
            nH = nH + MAXP-1
            nE = nE + NREVAR*MAXP
          case('mdlq')
            nH = nH + NRHVAR*MAXmdlqH
            nE = nE + NREVAR*MAXmdlqE
            nV = nV + NRVVAR*MAXmdlqV
          case('mdlb')
            nH = nH + NRHVAR*MAXmdlbH
            nE = nE + NREVAR*MAXmdlbE
            nV = nV + NRVVAR*MAXmdlbV
            nQ = nQ + NRQVAR*MAXmdlbQ
          end select
        enddo

        size_tot = (nH + nE + nV + nQ)*16
        write(*,*) '2: total DOF size is: ', size_tot
        call pause

end subroutine my_sizetest
!
!-------------------------
!Routine set_PML
!sets the PML data
!last Modified - Jan 2018
!-------------------------
subroutine set_PML
   use CommonParam
   implicit none
   PML_REGION=ZL-PML_FRAC*ZL
end subroutine set_PML
!
!-------------------------
!Routine get_Beta
!returns the PML function at a point Xp
! input : Xp (coordinates x,y,z at a physical point)
!         fld_flag: 1 for signal 0 for pump
! output:
!       zbeta - PML stretch function that stretches only in z-direction
!       zdbeta - z-derivative of PML stretch function
!       zd2beta - second z-derivative of PML stretch function
!last Modified - Jan 2018
!-------------------------
subroutine get_Beta(Xp,fld_flag, zbeta,zdbeta,zd2beta)
!
   use CommonParam
   use LaserParam
      implicit none
      real*8, dimension(3), intent(in) :: Xp
      integer, intent(in)              :: fld_flag
      real*8 :: z,a,b,c,L,n,rho,drho,d2rho
      VTYPE :: zbeta,zdbeta,zd2beta
      z=Xp(3)
      a = EXP_COEFF
      b = PML_REGION
      L = ZL
      c = 5.0d0
      n=9.d0
      zbeta = z
      zdbeta = 1.d0
! ... check if the wave is exponential growth
    if(SIGMA.ne.ZERO) then
        if(z.gt.b) then
          rho = a*c*(z-b)**n*exp(z/L)
          drho = a*c*(z-b)**(n-1.d0)*dexp(z/L)*(n+(z-b)/L)
          d2rho = a*c*dexp(z/L)*(z-b)**(n-2.d0)*(((z-b)/L)**2+2.d0*n*((z-b)/L)+n*(n-1))
          zbeta = z - ZI*rho/OMEGA
          zdbeta = 1.d0 - ZI*drho/OMEGA
          zd2beta = -ZI*d2rho/OMEGA
        endif
    else
      c = 25.0d0
      n=3.d0
! ...  check for co-pumping: PML on same side for both signal and pump
      if(COPUMP.eq.1) then
        select case(fld_flag)
        case(1)
          if(z.gt.PML_REGION) then
            rho = c*((z-b)/(ZL-PML_REGION))**n
            drho = c*n*((z-b)/(ZL-PML_REGION))**(n-1.d0)*(1.d0/(ZL-PML_REGION))
            d2rho = c*n*(n-1.d0)*((z-b)/(ZL-PML_REGION))**(n-2.d0)*(1.d0/(ZL-PML_REGION)**2)
            if((rho.le.0.d0).or.(drho.le.0.d0).or.(d2rho.le.0)) then
              write(*,*) ' from get_Beta: rho, drho,d2rho are negative'
              stop 1
            endif
            zbeta = z - ZI*rho/(OMEGA*OMEGA_RATIO_SIGNAL)
            zdbeta = 1.d0 - ZI*drho/(OMEGA*OMEGA_RATIO_SIGNAL)
            zd2beta = -ZI*d2rho/(OMEGA*OMEGA_RATIO_SIGNAL)
          endif
        case(0)
          if(z.gt.PML_REGION) then
            rho = c*((z-b)/(ZL-PML_REGION))**n
            drho = c*n*((z-b)/(ZL-PML_REGION))**(n-1.d0)*(1.d0/(ZL-PML_REGION))
            d2rho = c*n*(n-1.d0)*((z-b)/(ZL-PML_REGION))**(n-2.d0)*(1.d0/(ZL-PML_REGION)**2)
            if((rho.le.0.d0).or.(drho.le.0.d0).or.(d2rho.le.0)) then
              write(*,*) ' from get_Beta: rho, drho,d2rho are negative'
              stop 1
            endif
            zbeta = z - ZI*rho/(OMEGA*OMEGA_RATIO_PUMP)
            zdbeta = 1.d0 - ZI*drho/(OMEGA*OMEGA_RATIO_PUMP)
            zd2beta = -ZI*d2rho/(OMEGA*OMEGA_RATIO_PUMP)
          endif
        endselect
! ...  next check for counter-pumping: PML on opposite side for signal and pump
! ...  PML @ z.gt.PML_REGION for signal and @ z.lt.PML_FRAC*ZL for pump
      elseif(COPUMP.eq.0) then
        select case(fld_flag)
        case(1)
          if(z.gt.PML_REGION) then
            rho = c*((z-b)/(ZL-PML_REGION))**n
            drho = c*n*((z-b)/(ZL-PML_REGION))**(n-1.d0)*(1.d0/(ZL-PML_REGION))
            d2rho = c*n*(n-1.d0)*((z-b)/(ZL-PML_REGION))**(n-2.d0)*(1.d0/(ZL-PML_REGION)**2)
            if((rho.le.0.d0).or.(drho.le.0.d0).or.(d2rho.le.0)) then
              write(*,*) ' from get_Beta: rho, drho,d2rho are negative'
              stop 1
            endif
            zbeta = z - ZI*rho/(OMEGA*OMEGA_RATIO_SIGNAL)
            zdbeta = 1.d0 - ZI*drho/(OMEGA*OMEGA_RATIO_SIGNAL)
            zd2beta = -ZI*d2rho/(OMEGA*OMEGA_RATIO_SIGNAL)
          endif
        case(0)
          if(z.lt.PML_FRAC*ZL) then
            rho = c*((z)/(ZL-PML_REGION))**n
            drho = c*n*((z)/(ZL-PML_REGION))**(n-1.d0)*(1.d0/(ZL-PML_REGION))
            d2rho = c*n*(n-1.d0)*((z)/(ZL-PML_REGION))**(n-2.d0)*(1.d0/(ZL-PML_REGION)**2)
            if((rho.le.0.d0).or.(drho.le.0.d0).or.(d2rho.le.0)) then
              write(*,*) ' from get_Beta: rho, drho,d2rho are negative'
              stop 1
            endif
            zbeta = z - ZI*rho/(OMEGA*OMEGA_RATIO_PUMP)
            zdbeta = 1.d0 - ZI*drho/(OMEGA*OMEGA_RATIO_PUMP)
            zd2beta = -ZI*d2rho/(OMEGA*OMEGA_RATIO_PUMP)
          endif
        endselect
!  ... endif for COPUMP check
      endif
! ... endif for the wave is exponential growth
    endif

   end subroutine get_Beta


!-------------------------
!Routine get_bgPol
!returns the (linear) background polarization
! input : domain flag: core - 1, cladding - 0
!         fld_flag: signal - 1, pump - 0
! output:
!       bg_pol - value of background polarization
!last Modified - April 2018
!-------------------------
subroutine get_bgPol(dom_flag,fld_flag, bg_pol)
!
      use CommonParam
      use LaserParam
!
      implicit none
      integer, intent(in) :: dom_flag,fld_flag
      VTYPE, intent(out) :: bg_pol
      if (dom_flag.eq.1) then
        bg_pol = (REF_INDEX_CORE**2-1.d0)
      elseif (dom_flag.eq.0) then
        bg_pol = (REF_INDEX_CLAD**2-1.d0)
      else
        write(*,*) ' error in get_bgPol: dom_flag must be 0 or 1'
        stop
      endif
!
      select case(fld_flag)
      case(1)
        bg_pol = ZI*OMEGA*OMEGA_RATIO_SIGNAL*bg_pol
      case(0)
        bg_pol = ZI*OMEGA*OMEGA_RATIO_PUMP*bg_pol
      case default
        write(*,*) ' error in get_bgPol: fld_flag must be 0 or 1'
        stop
      endselect

      if(real(bg_pol).ne.0.d0) then
        write(*,*) 'error from get_bgPol: bg_pol must be purely imaginary'
        stop 1
      endif

end subroutine get_bgPol
!
!-------------------------
!Routine get_ramanPol
!returns the Raman polarization
! input : domain flag: core - 1, cladding - 0
!         field flag: signal - 1, pump - 0
!         E,H  - electric and magnetic fields (L2 variables)
! output:
!       raman_pol - value of background polarization
!last Modified - Jan 2018
!-------------------------
subroutine get_ramanPol(E,H,dom_flag,fld_flag, raman_pol)
!
      use CommonParam
      use LaserParam
!
      implicit none
      VTYPE, dimension(3), intent(in) :: E,H
      integer, intent(in) :: dom_flag
      integer, intent(in) :: fld_flag
      VTYPE, intent(out) :: raman_pol
      VTYPE,  dimension(3) :: EtimesH
      VTYPE :: g_lR
      call zz_cross_product(E,conjg(H), EtimesH)
      g_lR = sqrt((real(EtimesH(1))**2+real(EtimesH(2))**2+real(EtimesH(3))**2))*RAMAN_GAIN
      select case(fld_flag)
      case(1)
        if (dom_flag.eq.1) then
          raman_pol = ZI*REF_INDEX_CORE*UPSILON_RAMAN_SIGNAL*g_lR/(OMEGA*OMEGA_RATIO_SIGNAL)
          raman_pol = ZI*OMEGA*OMEGA_RATIO_SIGNAL*raman_pol
        elseif (dom_flag.eq.0) then
          raman_pol = ZI*REF_INDEX_CLAD*UPSILON_RAMAN_SIGNAL*g_lR/(OMEGA*OMEGA_RATIO_SIGNAL)
          raman_pol = ZI*OMEGA*OMEGA_RATIO_SIGNAL*raman_pol
        else
          write(*,*) ' error in get_ramanPol: dom_flag must be 0 or 1'
          stop 1
        endif
        if((real(raman_pol).gt.0.d0).or.(aimag(raman_pol).ne.0.d0)) then
          write(*,*) 'error in get_ramanPol: for signal, Raman polarization must be purely real with negative real part'
          stop 1
        endif
      case(0)
        if (dom_flag.eq.1) then
          raman_pol = ZI*REF_INDEX_CORE*UPSILON_RAMAN_PUMP*g_lR/(OMEGA*OMEGA_RATIO_PUMP)
          raman_pol = ZI*OMEGA*OMEGA_RATIO_PUMP*raman_pol
        elseif (dom_flag.eq.0) then
          raman_pol = ZI*REF_INDEX_CLAD*UPSILON_RAMAN_PUMP*g_lR/(OMEGA*OMEGA_RATIO_PUMP)
          raman_pol = ZI*OMEGA*OMEGA_RATIO_PUMP*raman_pol
        else
          write(*,*) ' error in get_ramanPol: dom_flag must be 0 or 1'
          stop 1
        endif
        if((real(raman_pol).lt.0.d0).or.(aimag(raman_pol).ne.0.d0)) then
          write(*,*) 'error in get_ramanPol: for pump, Raman polarization must be purely real with positive real part'
          stop 1
        endif
      case default
        write(*,*) 'error in get_ramanPol: fld_flag must be 0 or 1'
        stop 1
      endselect
end subroutine get_ramanPol
!
!-------------------------
!Routine get_NewtonLoad
!returns the rhs for Newton iteration
! input : domain flag: core - 1, cladding - 0
!         field flag: signal - 1, pump - 0
! output:
!       zNewtonLoad - value of Newton load at
!                     iteration n
!last Modified - Feb 2018
!
! The Newton load at step n corresponds to the semi-implicit
! Newton method used for the nonlinear iterations.
!
!  zNewtonLoad = - (n*gR*Upsilon*Real((E_k(n) \times conjg(H_k(n)), E_k(n) \times conjg(\delta H_k(n-1))))E(n) &
!               + n*gR*Upsilon*Real((E_k(n) \times conjg(H_k(n)), \delta E_k(n-1) \times conjg(H_k(n))))E(n)
!                       )/|E_k(n) \times conjg(H_k(n)|
!
!-------------------------
subroutine get_NewtonLoad(dom_flag,fld_flag,En,Hn,deltaE,deltaH, zNewtonLoad)
!
      use CommonParam
      use LaserParam
!
      implicit none
      integer, intent(in) :: dom_flag
      integer, intent(in) :: fld_flag
      VTYPE, dimension(3), intent(in)  :: En,Hn,deltaE,deltaH
      VTYPE, intent(out) :: zNewtonLoad
      VTYPE,  dimension(3) :: EtimesH, deltaEtimesH, EtimesdeltaH
      VTYPE                :: normEtimesH,dotprod1,dotprod2

      if(RAMAN_GAIN.eq.0.d0) then
        zNewtonLoad = 0.d0
        go to 777
      endif
      call zz_cross_product(En,conjg(Hn), EtimesH)
      call zz_cross_product(deltaE,conjg(Hn), deltaEtimesH)
      call zz_cross_product(En,conjg(deltaH), EtimesdeltaH)
      normEtimesH = sqrt(abs(EtimesH(1))**2+abs(EtimesH(2))**2+abs(EtimesH(3))**2)
      dotprod1 = dot_product(EtimesH,deltaEtimesH)
      dotprod2 = dot_product(EtimesH,EtimesdeltaH)
      if(normEtimesH.eq.0.d0) then
        zNewtonLoad = 0.d0
        go to 777
      endif
      select case(fld_flag)
      case(1)
        if (dom_flag.eq.1) then
          zNewtonLoad = -REF_INDEX_CORE*RAMAN_GAIN*UPSILON_RAMAN_SIGNAL*(real(dotprod1)+real(dotprod2))/(normEtimesH)
        elseif (dom_flag.eq.0) then
          zNewtonLoad = -REF_INDEX_CLAD*RAMAN_GAIN*UPSILON_RAMAN_SIGNAL*(real(dotprod1)+real(dotprod2))/(normEtimesH)
        else
          write(*,*) ' error in get_NewtonLoad: dom_flag must be 0 or 1'
          stop
        endif
      case(0)
        if (dom_flag.eq.1) then
          zNewtonLoad = -REF_INDEX_CORE*RAMAN_GAIN*UPSILON_RAMAN_PUMP*(real(dotprod1)+real(dotprod2))/(normEtimesH)
        elseif (dom_flag.eq.0) then
          zNewtonLoad = -REF_INDEX_CLAD*RAMAN_GAIN*UPSILON_RAMAN_PUMP*(real(dotprod1)+real(dotprod2))/(normEtimesH)
        else
          write(*,*) ' error in get_NewtonLoad: dom_flag must be 0 or 1'
          stop
        endif
      case default
        write(*,*) 'error in get_NewtonLoad: fld_flag must be 0 or 1'
        stop
      endselect
777   continue
end subroutine get_NewtonLoad


!-------------------------
!Routine get_thermPol
!returns the thermal polarization
! input : dom_flag: core - 1, cladding - 0
!         fld_flag: signal - 1, pump - 0
!         theta: value of solution of heat equation (H1 variable)
! output:
!       therm_pol - value of background polarization
!last Modified - May 2018
!-------------------------
subroutine get_thermPol(dom_flag,fld_flag,theta, therm_pol)
!
      use CommonParam
      use LaserParam
!
      implicit none
      integer, intent(in) :: dom_flag,fld_flag
      VTYPE, intent(in) :: theta
      VTYPE, intent(out) :: therm_pol
      VTYPE :: dn
      dn = ALPHA_HEAT*theta
      select case(fld_flag)
      case(1)
        if (dom_flag.eq.1) then
          therm_pol = ZI*OMEGA*OMEGA_RATIO_SIGNAL*2.d0*dn*REF_INDEX_CORE
        elseif (dom_flag.eq.0) then
          therm_pol = ZI*OMEGA*OMEGA_RATIO_SIGNAL*2.d0*dn*REF_INDEX_CLAD
        else
          write(*,*) ' error in get_thermPol: dom_flag must be 0 or 1'
          stop
        endif
      case(0)
        if (dom_flag.eq.1) then
          therm_pol = ZI*OMEGA*OMEGA_RATIO_PUMP*2.d0*dn*REF_INDEX_CORE
        elseif (dom_flag.eq.0) then
          therm_pol = ZI*OMEGA*OMEGA_RATIO_PUMP*2.d0*dn*REF_INDEX_CLAD
        else
          write(*,*) ' error in get_thermPol: dom_flag must be 0 or 1'
          stop
        endif
      case default
        write(*,*) 'error in get_thermPol: fld_flag must be 0 or 1'
        stop
      endselect
      if((real(therm_pol).ne.0.d0)) then
        write(*,*) 'error in get_thermPol: therm_pol must be purely imaginary'
        stop
      endif
end subroutine get_thermPol
!
!-------------------------
!Routine get_thermLoad
!returns the thermal polarization
! input : zsolQ - all 12 EH fields (6 EH- for signal, 6 EH - for pump)
! output:
!       therm_Load - value of thermal load
!last Modified - May 2018
!-------------------------
subroutine get_thermLoad(zsolQ, therm_Load)
!
      use CommonParam
      use LaserParam
!
      implicit none
      VTYPE, dimension(12), intent(in) :: zsolQ
      VTYPE, intent(out) :: therm_Load
      VTYPE, dimension(3) :: Es,Hs,Ep,Hp,EsTimesHs,EpTimesHp
      VTYPE :: Ip,Is
      Es = zsolQ(1:3)
      Hs = zsolQ(4:6)
      Ep = zsolQ(7:9)
      Hp = zsolQ(10:12)
      call zz_cross_product(Es,conjg(Hs), EsTimesHs)
      call zz_cross_product(Ep,conjg(Hp), EpTimesHp)
      Is = sqrt((real(EstimesHs(1))**2+real(EstimesHs(2))**2+real(EstimesHs(3))**2))
      Ip = sqrt((real(EptimesHp(1))**2+real(EptimesHp(2))**2+real(EptimesHp(3))**2))
      therm_Load =-KAPPA*(RAMAN_GAIN**2)*Is*Ip*(UPSILON_RAMAN_SIGNAL+UPSILON_RAMAN_PUMP)
end subroutine get_thermLoad
!
!-------------------------
!Routine get_activePol
!returns the active gain polarization
! input : zsolQ - all 12 EH fields (6 EH- for signal, 6 EH - for pump)
!         dom_flag - 1 for core, 0 for cladding
!         fld_flag - 1 for signal, 0 for pump
! output:
!       active_pol - value of active polarization
!last Modified - Jan 2018
!-------------------------
subroutine get_activePol(zsolQ,dom_flag,fld_flag, active_pol)
!
   use CommonParam
   use LaserParam
!
   implicit none
    VTYPE, dimension(12), intent(in) :: zsolQ
    integer, intent(in) :: dom_flag,fld_flag
    VTYPE, intent(out) :: active_pol
    VTYPE, dimension(3) :: Es,Hs,Ep,Hp,EsTimesHs,EpTimesHp
    real*8 :: Nex,Ngd
    real*8 :: sigma_s_abs_scaled,sigma_s_em_scaled
    real*8 :: sigma_p_abs_scaled, sigma_p_em_scaled
    real*8 :: sum1,sum2,Is,Ip

    sigma_s_abs_scaled = SIGMA_S_ABS*N0*CORE_RAD
    sigma_s_em_scaled = SIGMA_S_EM*N0*CORE_RAD
    sigma_p_abs_scaled = SIGMA_P_ABS*N0*CORE_RAD
    sigma_p_em_scaled = SIGMA_P_EM*N0*CORE_RAD
    Es = zsolQ(1:3)
    Hs = zsolQ(4:6)
    Ep = zsolQ(7:9)
    Hp = zsolQ(10:12)
    call zz_cross_product(Es,conjg(Hs), EsTimesHs)
    call zz_cross_product(Ep,conjg(Hp), EpTimesHp)
    Is = sqrt(abs(EsTimesHs(1))**2+abs(EsTimesHs(2))**2+abs(EsTimesHs(3))**2)
    Ip = sqrt(abs(EpTimesHp(1))**2+abs(EpTimesHp(2))**2+abs(EpTimesHp(3))**2)

    sum1 = (sigma_s_abs_scaled*Is/(OMEGA_RATIO_SIGNAL*OMEGA))+ &
           (sigma_p_abs_scaled*Ip/(OMEGA_RATIO_PUMP*OMEGA))
    sum2 = (Is/(OMEGA_RATIO_SIGNAL*OMEGA)*(sigma_s_abs_scaled+sigma_s_em_scaled)) + &
          (Ip/(OMEGA_RATIO_PUMP*OMEGA)*(sigma_p_abs_scaled+sigma_p_em_scaled))
    Nex = sum1/sum2
    if(Nex.gt.1.d0) then
      write(*,*) 'error from get_activePol: Nex must be <= 1'
    endif
    Ngd = 1.d0-Nex
    if(fld_flag.eq.1) then
      if(dom_flag.eq.1) then
        active_pol = (-ZI*REF_INDEX_CORE/(OMEGA_RATIO_SIGNAL*OMEGA))&
        *(sigma_s_em_scaled*Nex-sigma_s_abs_scaled*Ngd)
        active_pol = ZI*OMEGA_RATIO_SIGNAL*OMEGA*active_pol
      elseif(dom_flag.eq.0) then
        active_pol = (-ZI*REF_INDEX_CLAD/(OMEGA_RATIO_SIGNAL*OMEGA))&
        *(sigma_s_em_scaled*Nex-sigma_s_abs_scaled*Ngd)
        active_pol = ZI*OMEGA_RATIO_SIGNAL*OMEGA*active_pol
      else
        write(*,*) 'error from get_activePol: dom_flag must be 1 or 0'
      endif
    elseif(fld_flag.eq.0) then
      if(dom_flag.eq.1) then
        active_pol = (-ZI*REF_INDEX_CORE/(OMEGA_RATIO_PUMP*OMEGA))&
        *(sigma_p_em_scaled*Nex-sigma_p_abs_scaled*Ngd)
        active_pol = ZI*OMEGA_RATIO_PUMP*OMEGA*active_pol
      elseif(dom_flag.eq.0) then
        active_pol = (-ZI*REF_INDEX_CLAD/(OMEGA_RATIO_PUMP*OMEGA))&
        *(sigma_p_em_scaled*Nex-sigma_p_abs_scaled*Ngd)
        active_pol = ZI*OMEGA_RATIO_PUMP*OMEGA*active_pol
      else
        write(*,*) 'error from get_activePol: dom_flag must be 1 or 0'
      endif
    else
      write(*,*) 'error from get_activePol: fld_flag must be 1 or 0'
    endif

   end subroutine get_activePol
!
!
! from Jake's testing
! not used currently (never actually)
subroutine test_Bessel
   use CommonParam
   use LaserParam
   implicit none
!
!.....for bessel function modes
    real*8 :: zeta,rho,xi,theta,rcore,rcladding,alpha
    real*8 :: order, bessJ, bessY, dbessJ, dbessY,bessJc, bessYc, dbessJc, dbessYc
    real*8 :: bessI, bessK, dbessI, dbessK, bessIc, bessKc, dbessIc, dbessKc
    real*8 :: d2bessJ,d2bessY,d2bessI,d2bessK
    integer :: n

    order = ORDER_BESSEL

    do n=1,100
      xi = 10.d0/n
    call dbessJY(xi, order, bessJ, bessY, dbessJ, dbessY)
    call d2bessJY(xi, order, d2bessJ, d2bessY)
    call dbessIK(xi, order, bessI, bessK, dbessI, dbessK)
    call d2bessIK(xi, order, d2bessI, d2bessK)
    write(*,*) 'from test_Bessel : xi = ', xi
    write(*,*) 'from test_Bessel : bessJ, bessY, dbessJ, dbessY, d2bessJ, d2bessY= ', bessJ, bessY, dbessJ, dbessY, d2bessJ, d2bessY
    write(*,*) 'from test_Bessel : bessI, bessK, dbessI, dbessK,d2bessI, d2bessK  = ', bessI, bessK, dbessI, dbessK,d2bessI, d2bessK
   enddo
end subroutine test_Bessel
!
!-----------------------------------------------------------------------------------------
!> Purpose : routine copies all solution dofs from one component set to another
!!
!! @param[in]  No1  - component to copy from
!! @param[in]  No2  - component to copy to
!!
!! @revision Dec 17
!-----------------------------------------------------------------------------------------
!
subroutine copy_coms(No1,No2)
!
      use parameters
      use data_structure3D
      !
      IMPLICIT NONE
!
!     INPUT ARGUMENTS
      integer, intent(in)  :: No1,No2
!
!     LOCAL VARIABLES
      integer :: nod, nf, nt, nn2, i
!
!-----------------------------------------------------------------------------------------
!
!  ...check consistency
      if ((No1.lt.0).or.(No2.lt.0).or.(No1.gt.NRCOMS).or.(No2.gt.NRCOMS)) then
        write(*,*) 'copy_coms: No1,No2,NRCOMS = ', No1,No2,NRCOMS
        stop 1
      endif
      if (No1.eq.No2) return
!
!  ...loop through active nodes
      do nod=1,NRNODS
        if (NODES(nod)%act.eq.0) cycle
!
!  .....H1 dof
        nf = (No1-1)*NRHVAR
        nt = (No2-1)*NRHVAR
        nn2 = ubound(NODES(nod)%zdofH,2)
        if(nn2.gt.0) then
          do i=1,NRHVAR
            NODES(nod)%zdofH(nt+i,1:nn2) = NODES(nod)%zdofH(nf+i,1:nn2)
          enddo
        endif
!
!  .....H(curl) dof
        nf = (No1-1)*NREVAR
        nt = (No2-1)*NREVAR
        nn2 = ubound(NODES(nod)%zdofE,2)
        if(nn2.gt.0) then
          do i=1,NREVAR
            NODES(nod)%zdofE(nt+i,1:nn2) = NODES(nod)%zdofE(nf+i,1:nn2)
          enddo
        endif
!
!  .....H(div) dof
        nf = (No1-1)*NRVVAR
        nt = (No2-1)*NRVVAR
        nn2 = ubound(NODES(nod)%zdofV,2)
        if(nn2.gt.0) then
          do i=1,NRVVAR
            NODES(nod)%zdofV(nt+i,1:nn2) = NODES(nod)%zdofV(nf+i,1:nn2)
          enddo
        endif
!
!  .....L2 dof
        nf = (No1-1)*NRQVAR
        nt = (No2-1)*NRQVAR
        nn2 = ubound(NODES(nod)%zdofQ,2)
        if(nn2.gt.0) then
          do i=1,NRQVAR
            NODES(nod)%zdofQ(nt+i,1:nn2) = NODES(nod)%zdofQ(nf+i,1:nn2)
          enddo
        endif
!
!  ...end of loop through nodes
      enddo
!
!
end subroutine copy_coms
!
!-------------------------------------------------------------------
!Routine to evaluate the electric field norm of UW-Maxwell
!by integrating H(curl) trace solution on a face of middle node.
!the input variable fld_flag specifies signal or pump
!fld_flag = 1 - signal
!fld_flag = 0 - pump
!zpoint - value of z at cross section face
!Last modified : Feb 18
!input: Mdle,facenumber,fld_flag
!output: faceNorm,face_diff_norm
!-------------------------------------------------------------------
!face norms were used before 'get power' integration (probably not used anymore)
!
subroutine compute_element_faceNorm(mdle,facenumber,fld_flag,zpoint, faceNorm,face_diff_norm)
!
      use control
      use data_structure3D
      use environment      , only : L2PROJ
      use physics
      use parametersDPG
      use CommonParam
      implicit none
      integer, intent(in)                         :: mdle
      integer, intent(in)                         :: fld_flag
      integer, intent(in)                         :: facenumber
      real*8,  intent(in)                         :: zpoint
      real*8,  intent(out)                        :: faceNorm,face_diff_norm
!
!     element, face order, geometry dof
      integer,dimension(19)          :: norder
      real*8 ,dimension(3,MAXbrickH) :: xnod
      integer,dimension(12)          :: nedge_orient
      integer,dimension(6)           :: nface_orient
! ...face order
      integer, dimension(5)     :: norderf
! .... number of vertices,edge,faces per element type
      integer :: nrv, nre, nrf
!.......declare edge/face type varibles
      character(len=4) :: etype,ftype
!
!     geometry
! ... variables for geometry
      real*8, dimension(3)      :: xi,x,rn,x_new
      real*8, dimension(3,2)    :: dxidt,dxdt,rt
      real*8, dimension(3,3)    :: dxdxi,dxidx
      real*8, dimension(2)      :: t
      real*8                    :: rjac,bjac
!
!     2D quadrature data
      real*8, dimension(2,MAXNINT2ADD)  :: tloc
      real*8, dimension(MAXNINT2ADD)    :: wtloc
!
!     approximate solution dof's
      VTYPE, dimension(MAXEQNH,MAXbrickH) :: zdofH
      VTYPE, dimension(MAXEQNE,MAXbrickE) :: zdofE
      VTYPE, dimension(MAXEQNV,MAXbrickV) :: zdofV
      VTYPE, dimension(MAXEQNQ,MAXbrickQ) :: zdofQ
! ...H1 shape functions
      integer                         :: nrdofH
      real*8, dimension(MAXbrickH)    :: shapH
      real*8, dimension(3,MAXbrickH)  :: gradH
!
!     approximate solution
      VTYPE, dimension(  MAXEQNH  )   ::  zsolH
      VTYPE, dimension(  MAXEQNH,3)   :: zdsolH
      VTYPE, dimension(3,MAXEQNE  )   ::  zsolE
      VTYPE, dimension(3,MAXEQNE  )   :: zcurlE
      VTYPE, dimension(3,MAXEQNV  )   ::  zsolV
      VTYPE, dimension(  MAXEQNV  )   ::  zdivV
      VTYPE, dimension(  MAXEQNQ  )   ::  zsolQ

!     exact solution
      VTYPE,dimension(  MAXEQNH    )  ::   ValH
      VTYPE,dimension(  MAXEQNH,3  )  ::  DvalH
      VTYPE,dimension(  MAXEQNH,3,3)  :: d2valH
      VTYPE,dimension(3,MAXEQNE    )  ::   ValE
      VTYPE,dimension(3,MAXEQNE,3  )  ::  DvalE
      VTYPE,dimension(3,MAXEQNE,3,3)  :: d2valE
      VTYPE,dimension(3,MAXEQNV    )  ::   ValV
      VTYPE,dimension(3,MAXEQNV,3  )  ::  DvalV
!
!     exact solution (UNUSED)
      VTYPE,dimension(3,MAXEQNV,3,3)  :: d2valV
      VTYPE,dimension(  MAXEQNQ    )  ::   valQ
      VTYPE,dimension(  MAXEQNQ,3  )  ::  dvalQ
      VTYPE,dimension(  MAXEQNQ,3,3)  :: d2valQ
!
!     miscellanea
      integer :: nint,icase,iattr,l,i,j
      real*8  :: weight,wa
      integer :: iel,nsign
!
!     printing flag
      integer :: iprint,nflag,iload
!
!---------------------------------------------------------------------------------------
!
      iprint=0
!
      faceNorm = 0.d0
      face_diff_norm = 0.0d0
      nflag = 1
! ...element type
      etype = NODES(mdle)%type
      nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
      call find_order(mdle, norder)
      call find_orient(mdle, nedge_orient,nface_orient)
      call nodcor(mdle, xnod)
      call solelm(mdle, zdofH,zdofE,zdofV,zdofQ)
! .....sign factor to determine the OUTWARD normal unit vector
      nsign = nsign_param(etype,facenumber)
!
! .....face type
      ftype = face_type(etype,facenumber)
!
! .....face order of approximation
      call face_order(etype,facenumber,norder, norderf)
!
! .....set 2D quadrature
      INTEGRATION = NORD_ADD
      call set_2Dint(ftype,norderf, nint,tloc,wtloc)
      INTEGRATION = 0

      do l=1,nint
!
! .......face coordinates
        t(1:2) = tloc(1:2,l)
!
! .......face parametrization
        call face_param(etype,facenumber,t, xi,dxidt)
!
! .......determine element H1 shape functions (for geometry)
        call shape3H(etype,xi,norder,nedge_orient,nface_orient, &
                     nrdofH,shapH,gradH)
!
! .......geometry
        call bgeom3D(mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign, &
                     x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
        weight = bjac*wtloc(l)

        call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod,zdofH,zdofE,zdofV,zdofQ,nflag, &
                         x,dxdxi,zsolH,zdsolH,zsolE,zcurlE,zsolV,zdivV,zsolQ)
        if(NEXACT.eq.1) then
          call exact(x,Mdle, ValH,DvalH,d2valH, ValE,DvalE,d2valE,  &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
        endif
!
!       accumulate Norm of E - trace variable for signal (fld_flag=1) or pump (fld_flag=0)
!       i.e., integrate (abs(E(1))**2+abs(E(2))**2+abs(E(3))**2) with:
!                     E corresponding to signal if fld_flag = 1
!                     E corresponding to pump if fld_flag = 0
! .... first check for signal, i.e, if fld_flag = 1
          if(fld_flag.eq.1) then
            faceNorm = faceNorm+ &
            sqrt(abs(zsolE(1,1))**2+abs(zsolE(2,1))**2+abs(zsolE(3,1))**2)*weight
!   ..... if we have an exact
            if(NEXACT.eq.1) then
              face_diff_norm = face_diff_norm   &
                           + sqrt(abs(valE(1,1)-zsolE(1,1))**2+ &
                           + abs(valE(2,1)-zsolE(2,1))**2+abs(valE(3,1)-zsolE(3,1))**2)*weight
            endif
! .... next check for pump, i.e, if fld_flag = 0
          else if(fld_flag.eq.0) then
            faceNorm = faceNorm+ &
            sqrt(abs(zsolE(1,3))**2+abs(zsolE(2,3))**2+abs(zsolE(3,3))**2)*weight
!   ..... if we have an exact
            if(NEXACT.eq.1) then
              face_diff_norm = face_diff_norm   &
                           + sqrt(abs(valE(1,3)-zsolE(1,3))**2+ &
                           + abs(valE(2,3)-zsolE(2,3))**2+abs(valE(3,3)-zsolE(3,3))**2)*weight
            endif
          else
            write(*,*) 'from compute_element_faceNorm: fld_flag must be 0 or 1'
            stop
          endif

!       end loop over integration points
      enddo

end subroutine compute_element_faceNorm
!
!-------------------------------------------------------------------
!Routine to evaluate the electric field norm of UW-Maxwell
!by along the cross sections specified by the vector of z-values
!zValues of length length_zpts
!the input variables:
!input: zValues,length_zpts,fld_flag
!       zValues: vector of z-values along z-axis to evaluate
!       length_zpts: length of zValues vector
!       fld_flag:
!          fld_flag = 1 - signal
!          fld_flag = 0 - pump
!Last modified : Jan 18
!output: Norm,diff_norm
!-------------------------------------------------------------------
!
subroutine compute_FaceNorm(zValues,length_zpts,fld_flag, norm,diff_norm)
!
      use data_structure3D
      implicit none

      integer, intent(in)                         :: length_zpts
      real*8, dimension(length_zpts), intent(in)  :: zValues
      integer, intent(in)                         :: fld_flag
      real*8, dimension(length_zpts), intent(out) :: norm,diff_norm
      real*8 :: faceNorm, face_diff_norm, elemNorm
! .... Mdle number
      integer                        :: mdle
!     element, face order, geometry dof
      real*8 ,dimension(3,MAXbrickH) :: xnod
! .... number of vertices,edge,faces per element type
      integer :: nrv, nre, nrf
!
!     miscellanea
      integer :: iel,i
!
!     printing flag
      integer :: iprint
!
!---------------------------------------------------------------------------------------
!
      iprint=0
!
! ... initialize outputs (vector of powers for all z-points)
      norm = 0.d0
      diff_norm = 0.d0

! ... initialize running powers computed (elements per z-point)
      faceNorm = 0.d0
      face_diff_norm  = 0.d0

      mdle=0
!
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        call nodcor(mdle, xnod)
        do i=1,length_zpts
          if((zValues(i).lt.xnod(3,8)).and.(zValues(i).gt.xnod(3,1))) then
            call compute_element_faceNorm(mdle,2,fld_flag,zValues(i), faceNorm,face_diff_norm)
            norm(i) = norm(i)+faceNorm
            diff_norm(i)  = diff_norm(i) +face_diff_norm
          endif
        enddo
      enddo
end subroutine compute_FaceNorm
!
!-------------------------------------------------------------------
!Driver routine to specify values of z-points in order to
! evaluate the face norm of E-trace variables at those points corresponding
!to linear problem
!zValues of length length_zpts
!the input variable:  NONE
!            output: NONE
!Last modified : Feb 18
!
!-------------------------------------------------------------------
subroutine get_faceNorm
!
  use CommonParam
  use LaserParam
  implicit none
  integer                           :: fld_flag
  real*8, allocatable, dimension(:) :: zValues
  real*8, allocatable, dimension(:) :: face_norm
  real*8, allocatable, dimension(:) :: diff_face_norm
  integer              :: zlength,numPts
  real*8 a,b
  integer i
  a = 0.05d0 ! Initial value
  b = 0.005d0
  numPts = 200
  if(ZL.gt.1.d0) then
    zlength = numPts*int(ZL)
  else
    zlength = numPts
  endif

   allocate(zValues(zlength), face_norm(zlength), &
            diff_face_norm(zlength))
! ... get z-values
   zValues = (/(((i-1)*b+a),i=1,zlength)/)
! ... make z-values irrational numbers so that they are not on element
! ... boundaries
   zValues = zValues*PI*(7.d0/22.d0)
   write(*,*) ' from get_faceNorm: zValues '
   do i =1,zlength
    write(*,*) zValues(i)
   enddo
! ... get face norm and difference of face norm for linear problem
   fld_flag = 1
   call compute_FaceNorm(zValues,zlength,fld_flag, face_norm,diff_face_norm)
   write(*,*) ' from get_faceNorm: face_norm'
   do i =1,zlength
    write(*,*) face_norm(i)
   enddo
   write(*,*) ' from get_faceNorm: diff_face_norm'
   do i =1,zlength
    write(*,*) diff_face_norm(i)
   enddo
   deallocate(zValues,face_norm,diff_face_norm)
end subroutine get_faceNorm
!
