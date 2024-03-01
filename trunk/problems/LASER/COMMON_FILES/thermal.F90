!
#include "typedefs.h"
!
!-------------------------------------------------------------------------------
!
!     routine name      - get_avgTemp
!
!-------------------------------------------------------------------------------
!
!     latest revision:  - Jan 2024
!
!> @brief         - routine sets up computation of average and peak
!                         temperature in the fiber at certain z-locations
!                         (peak temperature assumed at R=0 fiber axis)
!
!     arguments:
!        inout:
!             NumPts    - Num points at which temperature is computed
!             FileIter  -  -1: print to stdout
!                         >=0: print to file with suffix=FileIter
!
!-------------------------------------------------------------------------------
subroutine get_avgTemp(NumPts,FileIter)
!
   use commonParam
   use laserParam
   use mpi_wrapper
   use par_mesh , only: DISTRIBUTED,HOST_MESH
!
   implicit none
!
   integer, intent(in)    :: FileIter
   integer, intent(inout) :: NumPts
!
   real(8), allocatable :: zValues(:),coreTemp(:),peakTemp(:)
!
   real(8) :: a,b
   integer :: i
!
   character(8)  :: fmt,suffix
   character(64) :: filename
!
!-------------------------------------------------------------------------------
!
   if (NumPts.le.0) NumPts = 4
   if (RANK .eq. ROOT) then
      write(*,2001) '  get_avgTemp: Number of sample points: ', NumPts
 2001 format(A,i5)
   endif
!
   if ((.not. DISTRIBUTED .or. HOST_MESH) .and. RANK .ne. ROOT) goto 99
!
   allocate(zValues(NumPts),coreTemp(NumPts),peakTemp(NumPts))
!
!..distributing sample points uniformly
   if (RANK .eq. ROOT) then
      write(*,*) ' get_avgTemp: Distributing sample points uniformly along waveguide.'
      write(*,2002) '  ZL = ', ZL
 2002 format(A,f10.2,/)
   endif
   b = ZL/NumPts
   a = b/2.d0
   do i=1,NumPts
      zValues(i) = (i-1)*b+a
   enddo
!..irrationalize z values to avoid points on element interfaces
   zValues = zValues*PI*(113.d0/355.d0)
!
   if (RANK .eq. ROOT) write(*,*) ' get_avgTemp: computing core temperature values..'
   call comp_avgTemp(zValues,NumPts, coreTemp,peakTemp)
!
!..only ROOT proc has valid coreTemp and peakTemp values
   if (RANK .ne. ROOT) goto 90
!
   if (FileIter .eq. -1) then
      write(*,*) ' get_avgTemp: printing average core temperature values..'
      do i = 1,NumPts
         write(*,2020) coreTemp(i)
 2020    format('    ',es12.5)
      enddo
!
      write(*,*) ' get_avgTemp: printing peak temperature values..'
      do i = 1,NumPts
         write(*,2020) peakTemp(i)
      enddo
   endif
!
   if (FileIter .ge. 0) then
      write(*,*) ' get_avgTemp: printing average core temperature values to file..'
      fmt = '(I5.5)'
      write (suffix,fmt) FileIter
      filename=trim(OUTPUT_DIR)//'temp/avg_temp_'//trim(suffix)//'.dat'
      open(UNIT=9,FILE=filename,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
      do i = 1,NumPts
         write(UNIT=9, FMT="(es12.5)") coreTemp(i)
      enddo
      close(UNIT=9)
!
      write(*,*) ' get_avgTemp: printing peak temperature values to file..'
      write (suffix,fmt) FileIter
      filename=trim(OUTPUT_DIR)//'temp/peak_temp_'//trim(suffix)//'.dat'
      open(UNIT=9,FILE=filename,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
      do i = 1,NumPts
         write(UNIT=9, FMT="(es12.5)") peakTemp(i)
      enddo
      close(UNIT=9)
   endif
!
   90 continue
   deallocate(zValues,coreTemp,peakTemp)
!
   99 continue
!
end subroutine get_avgTemp
!
!
!-------------------------------------------------------------------------------
!
!     routine name        - compute_avgTemp
!
!-------------------------------------------------------------------------------
!
!     latest revision:    - Jan 2024
!
!> @brief           - routine computes average and peak temperature in the
!                           fiber at certain z-locations (sample points)
!
!     arguments:
!        in:
!              ZValues    - sample points in z-direction
!              NumPts     - number of sample points
!        out:
!              CoreTemp   - average temperature in fiber core
!                           computed for elements at ZValues
!              PeakTemp   - peak temperature in fiber core
!                           computed for elements at ZValues and R=0
!
!-------------------------------------------------------------------------------
subroutine comp_avgTemp(ZValues,NumPts, CoreTemp,PeakTemp)
!
   use commonParam
   use data_structure3D
   use control    , only : GEOM_TOL
   use environment, only : QUIET_MODE
   use mpi_wrapper
   use par_mesh   , only : DISTRIBUTED,HOST_MESH
!
   implicit none
!
   integer , intent(in)  :: NumPts
   real(8) , intent(in)  :: ZValues(NumPts)
   real(8) , intent(out) :: CoreTemp(NumPts)
   real(8) , intent(out) :: PeakTemp(NumPts)
!
!..mdle number
   integer :: mdle
!
!..element, face order, geometry dof
   real(8) :: xnod(3,8)
   real(8) :: maxz,minz,minr,rval
!
!..miscellanea
   integer :: iel,i,j,minj,ndom,nod,nv
   integer :: nodesl(27),norientl(27)
   real(8) :: elemTemp,elemVol
   real(8) :: coreVol(NumPts)
!
!..timer
   real(8) :: start_time,end_time
   integer :: count,ierr
!
!-------------------------------------------------------------------------------
!
   CoreTemp = rZERO
   PeakTemp = rZERO
   coreVol  = rZERO
!
!..start timer
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!
   if (.not. DISTRIBUTED) then
      ELEM_SUBD(1:NRELES) = ELEM_ORDER(1:NRELES)
      NRELES_SUBD = NRELES
   endif
!
!..iterate over elements
!$OMP PARALLEL DO                                              &
!$OMP PRIVATE(mdle,nv,xnod,maxz,minz,minr,rval,i,j,minj,ndom,  &
!$OMP         nod,nodesl,norientl,elemTemp,elemVol)            &
!$OMP REDUCTION(+:CoreTemp,coreVol)                            &
!$OMP REDUCTION(max:PeakTemp)                                  &
!$OMP SCHEDULE(DYNAMIC)
   do iel=1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
      if (GEOM_NO .eq. 5) then
         call find_domain(mdle, ndom)
         select case(ndom)
            case(1,2) ! core
            case default
               cycle ! skip cladding
         end select
      endif
      call nodcor_vert(mdle, xnod)
      nv = nvert(NODES(mdle)%ntype)
      maxz = maxval(xnod(3,1:nv))
      minz = minval(xnod(3,1:nv))
      do i=1,NumPts
         if((ZValues(i).le.maxz).and.(ZValues(i).gt.minz)) then
!        ...compute avg temp
            call comp_elem_avgTemp(mdle, elemTemp,elemVol)
            CoreTemp(i) = CoreTemp(i) + elemTemp*elemVol
            coreVol(i)  = coreVol(i)  + elemVol
!        ...compute peak temp
            minr = 1.d0; minj = 0
            do j=1,nv
               rval = sqrt(xnod(1,j)**2+xnod(2,j)**2)
               if (rval .lt. minr) then
                  minr = rval; minj = j
               endif
            enddo
            if (minr .lt. GEOM_TOL) then
               if (minj .eq. 0) then
                  !$OMP CRITICAL
                  write(*,*) 'compute_avgTemp: minj = 0.'
                  !$OMP END CRITICAL
                  cycle
               endif
               call elem_nodes(mdle, nodesl,norientl)
               nod = nodesl(minj)
               if (associated(NODES(nod)%dof)) then
                  if (associated(NODES(nod)%dof%zdofH)) then
                     PeakTemp(i) = real(NODES(nod)%dof%zdofH(1,1,N_COMS))
                  else
                     !$OMP CRITICAL
                     write(*,*) 'compute_avgTemp: NOT associated zdofH nod = ', nod
                     !$OMP END CRITICAL
                  endif
               else
                  !$OMP CRITICAL
                  write(*,*) 'compute_avgTemp: NOT associated dof nod = ', nod
                  !$OMP END CRITICAL
               endif
            endif
         endif
      enddo
   enddo
!$OMP END PARALLEL DO
!
!..gather information on host
   if (.not. DISTRIBUTED .or. HOST_MESH) goto 50
   count = NumPts
   if (RANK .eq. ROOT) then
      call MPI_REDUCE(MPI_IN_PLACE,CoreTemp,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE,coreVol ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE,PeakTemp,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
   else
      call MPI_REDUCE(CoreTemp,CoreTemp,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(coreVol ,coreVol ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(PeakTemp,PeakTemp,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
   endif
   50 continue
!
   if (RANK .ne. ROOT) goto 90
!
!..compute average temperature in fiber core at the sample points
   do i=1,NumPts
      CoreTemp(i) = CoreTemp(i)/coreVol(i)
      !write(*,3005) 'i = ',i, ', CoreTemp = ',CoreTemp(i),', CoreVol = ',coreVol(i)
 3005 format(A,I4,A,F6.2,A,F6.2)
   enddo
!
   90 continue
!
!..end timer
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
   if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) then
      write(*,3010) end_time-start_time
 3010 format('  compute_avgTemp : ',f12.5,'  seconds',/)
   endif
!
end subroutine comp_avgTemp
!
!
!-------------------------------------------------------------------------------
!
!     routine name        - comp_elem_avgTemp
!
!-------------------------------------------------------------------------------
!
!     latest revision:    - Mar 2019
!
!> @brief           - routine computes average temperature in
!                           an element
!
!     arguments:
!        in:
!              Mdle       - Element mdle node
!        out:
!              ElemTemp   - Average temperature in element
!              ElemVol    - Volume of element in physical space
!
!-------------------------------------------------------------------------------
subroutine comp_elem_avgTemp(Mdle, ElemTemp,ElemVol)
!
   use data_structure3D
   use control, only: INTEGRATION
!
   implicit none
!
   integer, intent(in)  :: Mdle
   real(8), intent(out) :: ElemTemp,ElemVol
!
!..element, face order, geometry dof
   integer :: norder(19)
   real(8) :: xnod(3,MAXbrickH)
   integer :: nedge_orient(12)
   integer :: nface_orient(6)
!
!..geometry
   real(8),dimension(3)   :: xi,x
   real(8),dimension(3,3) :: dxidx,dxdxi
   real(8)                :: rjac,wa,weight
!
!..3D quadrature data
   real(8) :: xiloc(3,MAX_NINT3)
   real(8) :: wxi(MAX_NINT3)
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
   real(8) :: rsolH
!
!-------------------------------------------------------------------------------
!
   nflag=1
!
!..order of approx, orientations, geometry dof's, solution dof's
   call find_order( Mdle, norder)
   call find_orient(Mdle, nedge_orient,nface_orient)
   call nodcor(     Mdle, xnod)
   call solelm(     Mdle, zdofH,zdofE,zdofV,zdofQ)
!
!..set up the element quadrature
   INTEGRATION=0
   call set_3D_int_DPG(NODES(Mdle)%ntype,norder,nface_orient, nint,xiloc,wxi)
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
      call soleval(Mdle,xi,nedge_orient,nface_orient,norder,xnod, &
                   zdofH,zdofE,zdofV,zdofQ,nflag, &
                   x,dxdxi,zsolH,zdsolH,zsolE,zcurlE,zsolV,zdivV,zsolQ)
      rsolH = real(zsolH(1))
!
!  ...Jacobian
      call geom(dxdxi, dxidx,rjac,iflag)
      weight = wa*rjac
#if HP3D_DEBUG
      if (iflag .ne. 0) then
        write(*,*) 'comp_elem_avgTemp: geom iflag != 0. Mdle,rjac = ',Mdle,rjac
      endif
#endif
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
!-------------------------------------------------------------------------------
!
!> @name get_thermLoad
!
!> @date Mar 2019
!
!> @brief returns the thermal load (heat deposition)
!
!  input:   - ZsolQ  : EH fields (6 EH - for signal, 6 EH - for pump)
!           - Coord_z: Element z-coordinate in physical coordinates
!
!  output:  - Therm_load
!
!-------------------------------------------------------------------------------
subroutine get_thermLoad(ZsolQ,Coord_z, Therm_load)
!
   use commonParam
   use laserParam
!
   implicit none
!
   VTYPE  , intent(in)  :: ZsolQ(12)
   real(8), intent(in)  :: Coord_z
   real(8), intent(out) :: Therm_load
!
   VTYPE, dimension(3) :: Es,Hs,Ep,Hp,ETimesHs,ETimesHp
!
   integer :: jz,numPts
   real(8) :: gs,gp,Is,Ip
   real(8) :: eta,sum1,sum2,dz
!
!-------------------------------------------------------------------------------
!
!..get fields
   Es = ZsolQ(1:3)
   Hs = ZsolQ(4:6)
   Ep = ZsolQ(7:9)
   Hp = ZsolQ(10:12)
!
!..compute signal irradiance
   call zz_cross_product(Es,conjg(Hs), ETimesHs)
   Is = sqrt((real(EtimesHs(1))**2+real(EtimesHs(2))**2+real(EtimesHs(3))**2))
!
!..compute pump irradiance
   Ip = 0.d0
   if (PLANE_PUMP .eq. 1) then
!     set plane pump power the same here and in active gain computation
!  ...assume pump is a plane wave in fiber cladding (cladding-pumped)
      Ip = PLANE_PUMP_POWER / (PI*R_CLAD*R_CLAD) ! calculate non-dimensional irradiance
   elseif (PLANE_PUMP .eq. 2) then
      numPts = size(PUMP_VAL)
      dz = ZL / numPts
      jz = min(INT(Coord_z/dz), numPts-1) + 1
      Ip = PUMP_VAL(jz) / (PI*R_CLAD*R_CLAD) ! calculate non-dimensional irradiance
   else
      call zz_cross_product(Ep,conjg(Hp), ETimesHp)
      Ip = sqrt((real(EtimesHp(1))**2+real(EtimesHp(2))**2+real(EtimesHp(3))**2))
   endif
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

