!
#include "implicit_none.h"
!
!--------------------------------------------------------------------
!
!     routine name      - get_avgTemp
!
!--------------------------------------------------------------------
!
!     latest revision:  - Mar 2019
!
!     purpose:          - routine computes average temperature in
!                         the fiber at certain z locations
!
!     arguments:
!        inout:
!             NumPts    - Num points at which temperature is computed
!             FileIter  -  -1: print to stdout
!                         >=0: print to file with suffix=FileIter
!
!---------------------------------------------------------------------
subroutine get_avgTemp(NumPts,FileIter)
!
   use commonParam
   use laserParam
!
   implicit none
!
   integer, intent(in)    :: FileIter
   integer, intent(inout) :: NumPts
!
   real*8, allocatable, dimension(:) :: zValues
   real*8, allocatable, dimension(:) :: coreTemp
!
   real*8  :: a,b
   integer :: i
!
   character*8  :: fmt,suffix
   character*64 :: filename
!
!----------------------------------------------------------------------
!
   if (NumPts.le.0) NumPts = 4
   write(*,2001) '  get_avgTemp: Number of sample points: ', NumPts
 2001 format(A,i3)
!
   allocate(zValues(NumPts),coreTemp(NumPts))
!
!..distributing sample points uniformly
   write(*,*) ' get_avgTemp: Distributing sample points uniformly along waveguide.'
   write(*,2002) '  ZL = ', ZL
 2002 format(A,f5.2)
   b = ZL/NumPts
   a = b/2.d0
   do i=1,NumPts
      zValues(i) = (i-1)*b+a
   enddo
!..irrationalize z values to avoid points on element interfaces
   zValues = zValues*PI*(7.d0/22.d0)
   write(*,*)
!
   write(*,*) ' get_avgTemp: computing core temperature values..'
   call comp_avgTemp(zValues,NumPts, coreTemp)
!
   if (FileIter .eq. -1) then
      write(*,*) ' get_avgTemp: printing core temperature values..'
      do i = 1,NumPts
         write(*,2020) coreTemp(i)
 2020    format('    ',es12.5)
      enddo
   endif
!
   if (FileIter .ge. 0) then
      write(*,*) ' get_avgTemp: printing core temperature values to file..'
      fmt = '(I5.5)'
      write (suffix,fmt) FileIter
      filename=trim(OUTPUT_DIR)//'temp/temp_'//trim(suffix)//'.dat'
      open(UNIT=9,FILE=filename,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
      do i = 1,NumPts
         write(UNIT=9, FMT="(es12.5)") coreTemp(i)
      enddo
      close(UNIT=9)
   endif
!
   deallocate(zValues,coreTemp)
!
end subroutine get_avgTemp
!
!
!--------------------------------------------------------------------
!
!     routine name        - compute_avgTemp
!
!--------------------------------------------------------------------
!
!     latest revision:    - Mar 2019
!
!     purpose:            - routine ...
!
!     arguments:
!        in:
!              ZValues    - sample points in z-direction
!              NumPts     - number of sample points
!        out:
!              CoreTemp   - average temperature in fiber core
!                           computed for elements at ZValues
!
!---------------------------------------------------------------------
subroutine comp_avgTemp(ZValues,NumPts, CoreTemp)
!
   use data_structure3D
   use environment, only : QUIET_MODE
   use commonParam
!
   implicit none
!
   integer, intent(in)  :: NumPts
   real*8 , intent(in)  :: ZValues(NumPts)
   real*8 , intent(out) :: CoreTemp(NumPts)
!
   real*8, parameter :: rZero = 0.d0
!
!..mdle number
   integer  :: mdle
   integer  :: mdlea(NRELES)
!
!..element, face order, geometry dof
   real*8  :: xnod (3,MAXbrickH)
   real*8 :: maxz,minz
!
!..miscellanea
   integer :: iel, i, ndom
   real*8  :: elemTemp,elemVol
   real*8  :: coreVol(NumPts)
!
!..element type
   character(len=4) :: etype
!
!..auxiliary variables for timing
   real*8 :: start, OMP_get_wtime
!
!----------------------------------------------------------------------
!
   CoreTemp = rZero
   coreVol  = rZero
!
!..start timer
   start = OMP_get_wtime()
!
!..mdle node list
   mdle=0
   do iel=1,NRELES
      call nelcon(mdle, mdle)
      mdlea(iel) = mdle
   enddo
!
!..iterate over elements
!$OMP PARALLEL DO                                                 &
!$OMP PRIVATE(mdle,etype,xnod,maxz,minz,i,ndom,elemTemp,elemVol)  &
!$OMP REDUCTION(+:coreTemp,coreVol)                               &
!$OMP SCHEDULE(DYNAMIC)
   do iel=1,NRELES
      mdle = mdlea(iel)
      if (GEOM_NO .eq. 5) then
         call find_domain(mdle, ndom)
         select case(ndom)
            case(1,2) ! core
            case default
               cycle ! skip cladding
         end select
      endif
      call nodcor(mdle, xnod)
      etype = NODES(mdle)%type
      select case(etype)
         case('mdlb')
            maxz = maxval(xnod(3,1:8))
            minz = minval(xnod(3,1:8))
         case('mdlp')
            maxz = maxval(xnod(3,1:6))
            minz = minval(xnod(3,1:6))
         case default
            write(*,*) 'comp_avgTemp: invalid etype param. stop.'
            stop
      end select
      do i=1,NumPts
         if((ZValues(i).le.maxz).and.(ZValues(i).gt.minz)) then
            call comp_elem_avgTemp(mdle, elemTemp,elemVol)
            CoreTemp(i) = CoreTemp(i) + elemTemp*elemVol
            coreVol(i)  = coreVol(i)  + elemVol
         endif
      enddo
   enddo
!$OMP END PARALLEL DO
!
   do i=1,NumPts
      CoreTemp(i) = CoreTemp(i)/coreVol(i)
!      write(*,3005) 'i = ',i, ', CoreTemp = ',CoreTemp(i),', CoreVol = ',coreVol(i)
! 3005 format(A,I3,A,F6.2,A,F6.2)
   enddo
!
!..end timer
   if (.not. QUIET_MODE) then
      write(*,3010) OMP_get_wtime()-start
 3010 format('  compute_avgTemp : ',f12.5,'  seconds',/)
   endif
!
end subroutine comp_avgTemp
!
!
!--------------------------------------------------------------------
!
!     routine name        - comp_elem_avgTemp
!
!--------------------------------------------------------------------
!
!     latest revision:    - Mar 2019
!
!     purpose:            - routine computes average temperature in
!                           an element
!
!     arguments:
!        in:
!              Mdle       - Element mdle node
!        out:
!              ElemTemp   - Average temperature in element
!              ElemVol    - Volume of element in physical space
!
!---------------------------------------------------------------------
subroutine comp_elem_avgTemp(Mdle, ElemTemp,ElemVol)
!
   use data_structure3D
   use control, only: INTEGRATION
!
   implicit none
!
   integer, intent(in)  :: Mdle
   real*8 , intent(out) :: ElemTemp,ElemVol
!
!..element, face order, geometry dof
   integer :: norder(19)
   real*8  :: xnod(3,MAXbrickH)
   integer :: nedge_orient(12)
   integer :: nface_orient(6)
!
!..geometry
   real*8,dimension(3)   :: xi,x
   real*8,dimension(3,3) :: dxidx,dxdxi
   real*8                :: rjac,wa,weight
!
!..3D quadrature data
   real*8 :: xiloc(3,MAX_NINT3)
   real*8 :: wxi(MAX_NINT3)
!
!..approximate solution dof's
   VTYPE :: zdofH(MAXEQNH,MAXbrickH)
   VTYPE :: zdofE(MAXEQNE,MAXbrickE)
   VTYPE :: zdofV(MAXEQNV,MAXbrickV)
   VTYPE :: zdofQ(MAXEQNQ,MAXbrickQ)
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
   integer :: l,nint,iflag,nflag
   real*8  :: rsolH
!
!---------------------------------------------------------------------
!
   nflag=1
!
!..order of approx, orientations, geometry dof's, solution dof's
   call find_order( mdle, norder)
   call find_orient(mdle, nedge_orient,nface_orient)
   call nodcor(     mdle, xnod)
   call solelm(     mdle, zdofH,zdofE,zdofV,zdofQ)
!
!..set up the element quadrature
   INTEGRATION=0
   call set_3Dint_DPG(NODES(mdle)%type,norder, nint,xiloc,wxi)
   INTEGRATION=0
!
   ElemVol  = 0.d0
   ElemTemp = 0.d0
!
!..loop through integration points
   do l=1,nint
!
!  ...Gauss point and weight
      xi(1:3)=xiloc(1:3,l); wa=wxi(l)
!
      call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod, &
                   zdofH,zdofE,zdofV,zdofQ,nflag, &
                   x,dxdxi,zsolH,zdsolH,zsolE,zcurlE,zsolV,zdivV,zsolQ)
      rsolH = real(zsolH(1))
!
!  ...Jacobian
      call geom(dxdxi, dxidx,rjac,iflag)
      weight = wa*rjac
!
!  ...accumulate element volume
      ElemVol = ElemVol + weight
!
!  ...accumulate element temperature
      ElemTemp = ElemTemp + rsolH * weight
!
!..end loop over integration points
   enddo
!
!..compute average temperature in element
   ElemTemp = ElemTemp / ElemVol
!
end subroutine comp_elem_avgTemp
!
!
!---------------------------------------------------------------------
!
!  routine: get_thermLoad
!
!  last modified: Mar 2019
!
!  purpose: returns the thermal load (heat deposition)
!
!  input:   - ZsolQ: EH fields (6 EH - for signal, 6 EH - for pump)
!
!  output:  - Therm_load
!
!---------------------------------------------------------------------
subroutine get_thermLoad(ZsolQ, Therm_load)
!
   use commonParam
   use laserParam
!
   implicit none
!
   VTYPE,  intent(in)   :: ZsolQ(12)
   real*8, intent(out)  :: Therm_load
!
   VTYPE, dimension(3) :: Es,Hs,Ep,Hp,ETimesHs,ETimesHp
!
!..modified irradiance experiment (birefringent fiber) 
   VTYPE, dimension(3) :: Es_mod,Hs_mod
   VTYPE, dimension(3) :: Ep_mod,Hp_mod
   integer :: modified
!
   real*8 :: Ip,Is,gp,gs
   real*8 :: eta,Nex,Ngd,sum1,sum2
!
!---------------------------------------------------------------------
!
!..get fields
   Es = ZsolQ(1:3)
   Hs = ZsolQ(4:6)
   Ep = ZsolQ(7:9)
   Hp = ZsolQ(10:12)
!
   Es_mod(1) = Es(1)+Es(2)
   Es_mod(2) = ZERO
   Es_mod(3) = Es(3)
!
   Hs_mod(1) = ZERO
   Hs_mod(2) = Hs(1)+Hs(2)
   Hs_mod(3) = Hs(3)
!
   Ep_mod(1) = Ep(1)+Ep(2)
   Ep_mod(2) = ZERO
   Ep_mod(3) = Ep(3)
!
   Hp_mod(1) = ZERO
   Hp_mod(2) = Hp(1)+Hp(2)
   Hp_mod(3) = Hp(3)
!
!..compute irradiance
   modified = 0
   if (modified .eq. 1) then
      call zz_cross_product(Es_mod,conjg(Hs_mod), ETimesHs)
      call zz_cross_product(Ep_mod,conjg(Hp_mod), ETimesHp)
   else
      call zz_cross_product(Es,conjg(Hs), ETimesHs)
      call zz_cross_product(Ep,conjg(Hp), ETimesHp)
   endif
   Is = sqrt((real(EtimesHs(1))**2+real(EtimesHs(2))**2+real(EtimesHs(3))**2))
   Ip = sqrt((real(EtimesHp(1))**2+real(EtimesHp(2))**2+real(EtimesHp(3))**2))
!
   sum1 = (SIGMA_S_ABS/OMEGA_SIGNAL)*Is+(SIGMA_P_ABS/OMEGA_PUMP)*Ip
   sum2 = ((SIGMA_S_ABS+SIGMA_S_EMS)/OMEGA_SIGNAL)*Is + &
          ((SIGMA_P_ABS+SIGMA_P_EMS)/OMEGA_PUMP)*Ip
!
   eta = sum1/(TAU_0+sum2)
!
!..signal gain
   gs = -SIGMA_S_ABS + (SIGMA_S_ABS+SIGMA_S_EMS)*eta
   gs = gs * N_TOTAL
!
!..pump gain
   gp = -SIGMA_P_ABS + (SIGMA_P_ABS+SIGMA_P_EMS)*eta
   gp = gp * N_TOTAL
!
   Therm_load = -Q_0*(gs*Is+gp*Ip)
!
end subroutine get_thermLoad
!
!---------------------------------------------------------------------
