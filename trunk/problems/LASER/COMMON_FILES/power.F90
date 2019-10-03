!
#include "implicit_none.h"
!
!----------------------------------------------------------------------
!
!   routine name       - get_power
!
!----------------------------------------------------------------------
!
!   latest revision - Aug 2019
!
!   purpose         - Driver routine for computing power in UW
!                     Maxwell, i.e. the Poynting vector at certain
!                     z-points for pump or signal.
!        ....Z-points are samples in this routine....
!
!   arguments       - Fld: 0 - pump
!                          1 - signal
!                          2 - both
!                   - NumPts
!                   - FileIter: -1: print to stdout
!                              >=0: print to file with suffix=FileIter
!
!----------------------------------------------------------------------
!
subroutine get_power(Fld,NumPts,FileIter)
!
   use commonParam
   use laserParam
   use mpi_param, only: RANK,ROOT
   use MPI      , only: MPI_COMM_WORLD,MPI_IN_PLACE,MPI_REAL8,MPI_SUM
!
   implicit none
!
   integer, intent(in)    :: Fld
   integer, intent(in)    :: FileIter
   integer, intent(inout) :: NumPts
!
   real*8, allocatable :: zValues(:)
   real*8, allocatable :: sign_power(:),pump_power(:)
   real*8, allocatable :: diff_power(:),efficiency(:)
   real*8, allocatable :: core_power(:),clad_power(:)
!
   real*8  :: a,b
   integer :: i
!
   character*8  :: fmt,suffix
   character*64 :: filename
!
   integer :: count,ierr
!
!----------------------------------------------------------------------
!
   if (Fld .ne. 0 .and. Fld .ne. 1 .and. Fld .ne. 2) then
      if (RANK.eq.ROOT) write(*,*) ' get_power: invalid Fld param. returning.'
      return
   endif
!
   if (NumPts.le.0) NumPts = 4
   if (RANK .eq. ROOT) then
      write(*,2001) '  get_power: Number of sample points: ', NumPts
 2001 format(A,i5)
   endif
!
   allocate(zValues(NumPts+1)   , sign_power(NumPts+1), &
            pump_power(NumPts+1), diff_power(NumPts+1), &
            core_power(NumPts+1), clad_power(NumPts+1)  )
!
!..distributing sample points uniformly
   if (RANK .eq. ROOT) then
      write(*,*) ' get_power: Distributing sample points uniformly along waveguide.'
      write(*,2002) ' ZL = ', ZL
 2002 format(A,F7.2,/)
   endif
   b = ZL/NumPts
   a = b/2.d0
   do i=1,NumPts
      zValues(i) = (i-1)*b+a
   enddo
!..irrationalize z values to avoid points on element interfaces
   zValues = zValues*PI*(7.d0/22.d0)
!
!..init arrays
   sign_power(1:NumPts+1) = 0.d0
   pump_power(1:NumPts+1) = 0.d0
   diff_power(1:NumPts+1) = 0.d0
   core_power(1:NumPts+1) = 0.d0
   clad_power(1:NumPts+1) = 0.d0
!
!..get power
   select case (Fld)
      case(0)
         if (RANK.eq.ROOT) write(*,*) ' get_power: computing pump_power..'
         call compute_power(zValues,NumPts,Fld, pump_power,diff_power,core_power,clad_power)
      case(1)
         if (RANK.eq.ROOT) write(*,*) ' get_power: computing sign_power..'
         call compute_power(zValues,NumPts,Fld, sign_power,diff_power,core_power,clad_power)
      case(2)
         if (RANK.eq.ROOT) write(*,*) ' get_power: computing sign_power and pump_power..'
         call compute_power(zValues,NumPts,0, pump_power,diff_power,core_power,clad_power)
         call compute_power(zValues,NumPts,1, sign_power,diff_power,core_power,clad_power)
      case default
         if (RANK.eq.ROOT) write(*,*) ' get_power: invalid Fld param. stop.'
         stop 1
   end select
   !write(*,*) 'finished..'
!
!..gather RHS vector information on host
   count = NumPts+1
   if (RANK .eq. ROOT) then
      call MPI_REDUCE(MPI_IN_PLACE,sign_power,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE,pump_power,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE,diff_power,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE,core_power,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE,clad_power,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
   else
      call MPI_REDUCE(sign_power,sign_power,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(pump_power,pump_power,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(diff_power,diff_power,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(core_power,core_power,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(clad_power,clad_power,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      goto 90
   endif
!
!..Print signal power output values
   if (Fld .eq. 1 .or. Fld .eq. 2) then
      if (FileIter .eq. -1) then
         write(*,*) ' get_power: printing power values (signal):'
         do i = 1,NumPts+1
            write(*,2020) sign_power(i)
       2020 format('    ',es12.5)
         enddo
         if (NONLINEAR_FLAG .eq. 0) then
            i = NumPts
            if (USE_PML) then
               i = (1.0d0 - PML_FRAC) * NumPts
            endif
            write(*,2021) (sign_power(1)-sign_power(i))/sign_power(1) * 100.d0
       2021 format(' Power loss: ',f6.2,' %',/)
         endif
      endif
      if (FileIter .ge. 0) then
         !WRITE TO FILE
         write(*,*) ' get_power: printing power values (signal) to file..'
         fmt = '(I5.5)'
         write (suffix,fmt) FileIter
         filename=trim(OUTPUT_DIR)//'power/signal_'//trim(suffix)//'.dat'
         open(UNIT=9,FILE=filename,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
         do i = 1,NumPts+1
            write(UNIT=9, FMT="(es12.5)") sign_power(i)
         enddo
         close(UNIT=9)
      endif
   endif
   !
   !..Print pump power output values
   if (Fld .eq. 0 .or. Fld .eq. 2) then
      if (FileIter .eq. -1) then
         write(*,*) ' get_power: printing power values (pump):'
         do i = 1,NumPts+1
            write(*,2020) pump_power(i)
         enddo
      endif
      if (FileIter .ge. 0) then
         !WRITE TO FILE
         write(*,*) ' get_power: printing power values (pump) to file..'
         fmt = '(I5.5)'
         write (suffix,fmt) FileIter
         filename=trim(OUTPUT_DIR)//'power/pump_'//trim(suffix)//'.dat'
         open(UNIT=9,FILE=filename,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
         do i = 1,NumPts+1
            write(UNIT=9, FMT="(es12.5)") pump_power(i)
         enddo
         close(UNIT=9)
      endif
   endif
   !
   !..Print fiber core power ratio
   if (GEOM_NO .eq. 5 .and. (Fld .eq. 1 .or. Fld .eq. 2)) then
      if (FileIter .eq. -1) then
         write(*,*) ' get_power: printing fiber core power ratio (signal):'
         do i = 1,NumPts+1
            write(*,2030) core_power(i)/sign_power(i)
       2030 format('    ',f8.4)
         enddo
      endif
      if (FileIter .ge. 0) then
         !WRITE TO FILE
         write(*,*) ' get_power: printing fiber core power ratio (signal) to file..'
         fmt = '(I5.5)'
         write (suffix,fmt) FileIter
         filename=trim(OUTPUT_DIR)//'power/ratio_'//trim(suffix)//'.dat'
         open(UNIT=9,FILE=filename,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
         do i = 1,NumPts+1
            write(UNIT=9, FMT="(f8.4)") core_power(i)/sign_power(i)
         enddo
         close(UNIT=9)
      endif
   endif
!
!..get efficiency
   if (Fld .eq. 2) then
      allocate(efficiency(NumPts))
      write(*,*) ' get_power: computing efficiency..'
      efficiency(1) = 0.d0
      do i = 2,NumPts+1
         if(COPUMP.eq.1) then
            efficiency(i) = (sign_power(i)-sign_power(1))&
                             /(pump_power(1)-pump_power(i))
         elseif(COPUMP.eq.0) then
            efficiency(i) = (sign_power(i)-sign_power(1))&
                             /(pump_power(NumPts)-pump_power(i))
         else
            write(*,*) ' get_power: COPUMP must be 1 or 0. stop.'
            stop
         endif
      enddo
      if (FileIter .eq. -1) then
         write(*,*) ' get_power: printing efficiency:'
         do i = 1,NumPts+1
            write(*,2030) efficiency(i)
         enddo
      endif
      if (FileIter .ge. 0) then
         !WRITE TO FILE
         write(*,*) ' get_power: printing efficiency to file..'
         fmt = '(I5.5)'
         write (suffix,fmt) FileIter
         filename=trim(OUTPUT_DIR)//'power/efficiency_'//trim(suffix)//'.dat'
         open(UNIT=9,FILE=filename,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
         do i = 1,NumPts+1
            write(UNIT=9, FMT="(f8.4)") efficiency(i)
         enddo
         close(UNIT=9)
      endif
      deallocate(efficiency)
   endif
!
   90 continue
   deallocate(zValues,sign_power,pump_power,diff_power,core_power,clad_power)
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
!   latest revision    - Aug 2019
!
!   purpose            - Evaluates the electric field power of UW
!                        Maxwell along the cross sections specified by
!                        the vector of zValues in the input
!
!   arguments
!        in:
!                      - ZValues     : sample points in z-direction
!                      - Num_zpts    : number of sample points
!                      - Fld         : 1 (signal) or 0 (pump)
!       out:
!                      - Power       : Absolute value of power
!                      - DiffPower   : Diff exact to computed power
!                                      (available if NEXAXT=1)
!                      - CorePower   : (available if GEOM_NO=5)
!                      - CladPower   : (available if GEOM_NO=5)
!
!----------------------------------------------------------------------
!
subroutine compute_power(ZValues,Num_zpts,Fld, Power,DiffPower,CorePower,CladPower)
!
   use commonParam
   use data_structure3D
   use control    , only : GEOM_TOL
   use environment, only : QUIET_MODE
   use mpi_param  , only : RANK,ROOT
   use MPI        , only : MPI_COMM_WORLD
   use par_mesh   , only : DISTRIBUTED
!
   implicit none
!
   integer, intent(in)  :: Num_zpts
   real*8,  intent(in)  :: ZValues(Num_zpts)
   integer, intent(in)  :: Fld
   real*8,  intent(out) :: Power(Num_zpts+1)
   real*8,  intent(out) :: DiffPower(Num_zpts+1)
   real*8,  intent(out) :: CorePower(Num_zpts+1)
   real*8,  intent(out) :: CladPower(Num_zpts+1)
!
!..auxiliary variables
   real*8 :: facePower, faceDiffPower, elemPower
!
!..mdle number
   integer :: mdle
!
!..element, face order, geometry dof
   real*8 :: xnod (3,MAXbrickH)
   real*8 :: maxz,minz
!
!..miscellanea
   integer :: iel, i, ndom
!
!..element type
   character(len=4) :: etype
!
!..face number over which power is computed
!  (in brick and prism, face 2 is face normal to xi3, at xi3=1)
   integer, parameter :: faceNum = 2
!  (in brick and prism, face 1 is face normal to xi3, at xi3=0)
   integer, parameter :: faceIn  = 1
!
!..timer
   real(8) :: MPI_Wtime,start_time,end_time
   integer :: ierr
!
!---------------------------------------------------------------------------------------
!
!..initialize outputs (vector of powers for all z-points)
   Power = 0.d0
   DiffPower = 0.d0
!
!..initialize core/clad power for fiber geometry
   CorePower = 0.d0
   CladPower = 0.d0
!
!..initialize running powers computed (elements per z-point)
   facePower = 0.d0
   faceDiffPower  = 0.d0
!
!..start timer
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime
!
   if (.not. DISTRIBUTED) then
      ELEM_SUBD(1:NRELES) = ELEM_ORDER(1:NRELES)
      NRELES_SUBD = NRELES
   endif
!
!..iterate over elements
!
!$OMP PARALLEL DO                                        &
!$OMP PRIVATE(mdle,etype,xnod,maxz,minz,i,ndom,          &
!$OMP         facePower,faceDiffPower)                   &
!$OMP REDUCTION(+:Power,DiffPower,corePower,cladPower)   &
!$OMP SCHEDULE(DYNAMIC)
   do iel=1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
      if (GEOM_NO .eq. 5) call find_domain(mdle, ndom)
      call nodcor(mdle, xnod)
      etype = NODES(Mdle)%type
      select case(etype)
         case('mdlb')
            maxz = maxval(xnod(3,1:8))
            minz = minval(xnod(3,1:8))
         case('mdlp')
            maxz = maxval(xnod(3,1:6))
            minz = minval(xnod(3,1:6))
         case default
            write(*,*) 'compute_power: invalid etype param. stop.'
            stop
      end select
      if (minz .lt. GEOM_TOL) then
         call compute_facePower(mdle,faceIn,Fld, facePower,faceDiffPower)
         Power(1) = Power(1) + abs(facePower)
         DiffPower(1) = DiffPower(1) + abs(faceDiffPower)
         if (GEOM_NO .eq. 5) then
            select case(ndom)
               case(1,2)
                  CorePower(1) = CorePower(1) + abs(facePower)
               case(3,4)
                  CladPower(1) = CladPower(1) + abs(facePower)
            end select
         endif
      endif
      do i=1,Num_zpts
         if((ZValues(i).le.maxz).and.(ZValues(i).gt.minz)) then
            call compute_facePower(mdle,faceNum,Fld, facePower,faceDiffPower)
            Power(i+1) = Power(i+1) + abs(facePower)
            DiffPower(i+1) = DiffPower(i+1) + abs(faceDiffPower)
            if (GEOM_NO .eq. 5) then
               select case(ndom)
                  case(1,2)
                     CorePower(i+1) = CorePower(i+1) + abs(facePower)
                  case(3,4)
                     CladPower(i+1) = CladPower(i+1) + abs(facePower)
               end select
            endif
         endif
      enddo
   enddo
!$OMP END PARALLEL DO
!
!..TODO check (why not abs(facePower) above??)
!..maybe relevant for counter pump configuration
!..take absolute value after integration
!   do i=1,Num_zpts
!      Power(i) = abs(Power(i))
!      DiffPower(i) = abs(DiffPower(i))
!      if (GEOM_NO .eq. 5) then
!         CorePower(i) = abs(CorePower(i))
!         CladPower(i) = abs(CladPower(i))
!      endif
!   enddo
!
!..end timer
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime
   if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) then
      write(*,3010) end_time-start_time
 3010 format('  compute_power : ',f12.5,'  seconds',/)
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
!   latest revision    - Oct 2018
!
!   purpose            - Evaluates the electric field power of UW
!                        Maxwell by integrating H(curl) trace solution
!                        on a face of a middle node.
!
!   arguments
!        in:
!                      - Mdle       : middle element node
!                      - Facenumber :
!                      - Fld        : 1 (signal) or 0 (pump)
!       out:
!                      - FacePower     :
!                      - FaceDiffPower :
!
!----------------------------------------------------------------------
!
subroutine compute_facePower(Mdle,Facenumber,Fld, FacePower,FaceDiffPower)
!
   use control
   use data_structure3D
   use environment, only : L2PROJ
   use physics
   use parametersDPG
   use commonParam
!
   implicit none
!
   integer, intent(in)  :: Mdle
   integer, intent(in)  :: Fld
   integer, intent(in)  :: Facenumber
   real*8,  intent(out) :: FacePower
   real*8,  intent(out) :: FaceDiffPower
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
!
!..for Poynting vector
   VTYPE, dimension(3)  :: EtimesH1,EtimesH2
   VTYPE                :: FdotN
!
!..miscellanea
   integer :: nint,icase,iattr,l,i,j
   real*8  :: weight,wa
   integer :: iel,nsign
   integer :: nflag,iload
!
!---------------------------------------------------------------------------------------
!
   facePower = 0.d0
   faceDiffPower = 0.0d0
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
      call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod, &
                   zdofH,zdofE,zdofV,zdofQ,nflag,x,dxdxi, &
                   zsolH,zdsolH,zsolE,zcurlE,zsolV,zdivV,zsolQ)
      if(NEXACT.eq.1) then
         call exact(x,Mdle, ValH,DvalH,d2valH, ValE,DvalE,d2valE, &
                            ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      endif
!
!     accumulate Poynting vector power for signal (Fld=1) or pump (Fld=0),
!     i.e., integrate (Real(n \dot ExH^*)) with:
!                     E/H corresponding to signal if Fld = 1
!                     E/H corresponding to pump   if Fld = 0
!  ...first check for signal, i.e, if Fld = 1
      if(Fld.eq.1) then
         call zz_cross_product(zsolE(1:3,1),conjg((zsolE(1:3,2))), EtimesH1)
         FdotN = EtimesH1(1)*rn(1)+EtimesH1(2)*rn(2)+EtimesH1(3)*rn(3)
         facePower = facePower + (real(FdotN))*weight
!     ...if we have an exact
         if(NEXACT.eq.1) then
            call zz_cross_product(valE(1:3,1),conjg((valE(1:3,2))), EtimesH2)
            faceDiffPower = faceDiffPower   &
                           + abs(((EtimesH1(1)*rn(1)+EtimesH1(2)*rn(2)+EtimesH1(3)*rn(3))*weight) - &
                           ((EtimesH2(1)*rn(1)+EtimesH2(2)*rn(2)+EtimesH2(3)*rn(3))*weight))
         endif
!  ...next check for pump, i.e, if Fld = 0
      else if(Fld.eq.0) then
         call zz_cross_product(zsolE(1:3,3),conjg((zsolE(1:3,4))), EtimesH1)
         FdotN = EtimesH1(1)*rn(1)+EtimesH1(2)*rn(2)+EtimesH1(3)*rn(3)
         facePower = facePower + (real(FdotN))*weight
!     ...if we have an exact
         if(NEXACT.eq.1) then
            call zz_cross_product(valE(1:3,3),conjg((valE(1:3,4))), EtimesH2)
            faceDiffPower = faceDiffPower   &
                            + abs(((EtimesH1(1)*rn(1)+EtimesH1(2)*rn(2)+EtimesH1(3)*rn(3))*weight) - &
                                 ((EtimesH2(1)*rn(1)+EtimesH2(2)*rn(2)+EtimesH2(3)*rn(3))*weight))
         endif
      else
         write(*,*) 'compute_facePower: Fld must be 0 or 1. stop.'
         stop
      endif
!..end loop over integration points
   enddo
!
end subroutine compute_facePower
