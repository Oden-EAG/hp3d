!
#include "typedefs.h"
!
!----------------------------------------------------------------------
!
!   routine name       - get_power
!
!----------------------------------------------------------------------
!
!   latest revision - Dec 2025
!
!   purpose         - Driver routine for computing power in Bend_UW
!                     Maxwell, i.e. the Poynting vector at certain
!                     theta-points.
!        ....theta-points are samples in this routine....
!
!   arguments       - NumPts
!                   - ModeProj
!                   - FileIter: -1: print to stdout
!                              >=0: print to file with suffix=FileIter
!
!----------------------------------------------------------------------
!
subroutine get_power(NumPts,ModeProj,FileIter,InputPower,FinalPower,FinalLossExp)
!
   use commonParam
   use mpi_wrapper
   use par_mesh , only: DISTRIBUTED,HOST_MESH
   use control, only: GEOM_TOL
   use bessel_evaluation, only: initialize_bessel_mode_parameters
!
   implicit none
!
   integer, intent(in)    :: FileIter
   logical, intent(in)    :: ModeProj
   integer, intent(inout) :: NumPts
   real(8), intent(out)   :: InputPower,FinalPower,FinalLossExp
!
   real(8), allocatable :: thValues(:)
   real(8), allocatable :: sign_power(:)
   real(8), allocatable :: diff_power(:),power_loss_exp(:)
   real(8), allocatable :: core_power(:),clad_power(:)
!
   real(8), allocatable :: power_M00(:),power_M01(:),power_M02(:)
   real(8), allocatable :: norm_M00(:),norm_M01(:),norm_M02(:)
   real(8), allocatable :: coef_M00_r(:),coef_M01_r(:),coef_M02_r(:)
   real(8), allocatable :: coef_M00_c(:),coef_M01_c(:),coef_M02_c(:)
!
   real(8) :: a,b,gain,loss
   integer :: i,j
!
   character(8)  :: fmt,suffix
   character(64) :: filename
!
   integer :: count,ierr
!
!----------------------------------------------------------------------
!
   if (RANK .eq. ROOT) then
      write(*,*) ' get_power: Starting power computation...'
   endif
!
   if (NumPts.le.3) NumPts = 4
   if (RANK .eq. ROOT) then
      write(*,2001) '  get_power: Number of sample points: ', NumPts
 2001 format(A,i5)
   endif
!
   if ((.not. DISTRIBUTED .or. HOST_MESH) .and. RANK .ne. ROOT) goto 99
!
   allocate(thValues(NumPts)   , sign_power(NumPts), &
            diff_power(NumPts), &
            core_power(NumPts), clad_power(NumPts)  )
!
   if (ModeProj) then
      allocate(power_M00(NumPts),norm_M00(NumPts),coef_M00_r(NumPts),coef_M00_c(NumPts))
      allocate(power_M01(NumPts),norm_M01(NumPts),coef_M01_r(NumPts),coef_M01_c(NumPts))
      allocate(power_M02(NumPts),norm_M02(NumPts),coef_M02_r(NumPts),coef_M02_c(NumPts))
   endif
!
!..distributing sample points uniformly
   if (RANK .eq. ROOT) then
      write(*,*) ' get_power: Distributing sample points uniformly along waveguide.'
      write(*,2002) ' Upper bound of theta, THUP = ', THUP
 2002 format(A,F8.6,/)
   endif
   b = THUP*(1.d0-PMLTHUP)/(NumPts-1)
   a = THUP*(1.d0-PMLTHUP)* 3.d0*GEOM_TOL ! correction to be slightly above the lower face
   do i=1,NumPts
      thValues(i) = (i-1)*b + a
   enddo
   if (RANK .eq. ROOT) write(*,*) ' get_power: Sample points (theta values) :',thValues(1:NumPts)
!
!..get power
   if (RANK.eq.ROOT) write(*,*) ' get_power: computing sign_power..'
   call compute_power_bent_slab(thValues,NumPts, 1 , sign_power,diff_power,core_power,clad_power)
!
   if (ModeProj) then
      if (RANK.eq.ROOT) write(*,*) ' get_power: computing signal mode_power..'
      !
      !store current ISOL
      i = ISOL
      !
      ISOL = 300 ! M00 projection
      call initialize_bessel_mode_parameters(ISOL,RBEND)
      call compute_power_bent_slab(thValues,NumPts,ISOL, power_M00,norm_M00,coef_M00_r,coef_M00_c)
      ISOL = 301 ! M01 projection
      call initialize_bessel_mode_parameters(ISOL,RBEND)
      call compute_power_bent_slab(thValues,NumPts,ISOL, power_M01,norm_M01,coef_M01_r,coef_M01_c)
      ISOL = 302 ! M02 projection
      call initialize_bessel_mode_parameters(ISOL,RBEND)
      call compute_power_bent_slab(thValues,NumPts,ISOL, power_M02,norm_M02,coef_M02_r,coef_M02_c)
      !
      !restore  ISOL
      ISOL = i
      call initialize_bessel_mode_parameters(ISOL,RBEND)
      !
   endif
!
!..gather all values on host
   if (.not. DISTRIBUTED .or. HOST_MESH) goto 50
   count = NumPts
   if (RANK .eq. ROOT) then
      call MPI_REDUCE(MPI_IN_PLACE,sign_power,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE,diff_power,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE,core_power,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE,clad_power,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      !
      if (ModeProj) then
         call MPI_REDUCE(MPI_IN_PLACE,power_M00  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(MPI_IN_PLACE, norm_M00  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(MPI_IN_PLACE, coef_M00_r,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(MPI_IN_PLACE, coef_M00_c,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         !
         call MPI_REDUCE(MPI_IN_PLACE,power_M01  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(MPI_IN_PLACE, norm_M01  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(MPI_IN_PLACE, coef_M01_r,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(MPI_IN_PLACE, coef_M01_c,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         !
         call MPI_REDUCE(MPI_IN_PLACE,power_M02  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(MPI_IN_PLACE, norm_M02  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(MPI_IN_PLACE, coef_M02_r,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE(MPI_IN_PLACE, coef_M02_c,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         !
      endif
   else
      call MPI_REDUCE(sign_power,sign_power,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(diff_power,diff_power,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(core_power,core_power,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(clad_power,clad_power,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      !
      if (ModeProj) then
         call MPI_REDUCE(power_M00  ,power_M00  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE( norm_M00  , norm_M00  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE( coef_M00_r, coef_M00_r,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE( coef_M00_c, coef_M00_c,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         !
         call MPI_REDUCE(power_M01  ,power_M01  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE( norm_M01  , norm_M01  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE( coef_M01_r, coef_M01_r,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE( coef_M01_c, coef_M01_c,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         !
         call MPI_REDUCE(power_M02  ,power_M02  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE( norm_M02  , norm_M02  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE( coef_M02_r, coef_M02_r,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         call MPI_REDUCE( coef_M02_c, coef_M02_c,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         !
      endif
      goto 90
   endif
!
   50 continue
!
   if (ModeProj) then
!$OMP PARALLEL DO
      do i = 1,NumPts
         norm_M00(i)   = sqrt(norm_M00(i))
         coef_M00_r(i) = sqrt(coef_M00_r(i)**2.d0+coef_M00_c(i)**2.d0) ! / norm_M00(i)
         power_M00(i)  = power_M00(i) * ((coef_M00_r(i) / norm_M00(i))**2.d0)
         !
         norm_M01(i)   = sqrt(norm_M01(i))
         coef_M01_r(i) = sqrt(coef_M01_r(i)**2.d0+coef_M01_c(i)**2.d0) ! / norm_M01(i)
         power_M01(i)  = power_M01(i) * ((coef_M01_r(i) / norm_M01(i))**2.d0)
         !
         norm_M02(i)   = sqrt(norm_M02(i))
         coef_M02_r(i) = sqrt(coef_M02_r(i)**2.d0+coef_M02_c(i)**2.d0) ! / norm_M02(i)
         power_M02(i)  = power_M02(i) * ((coef_M02_r(i) / norm_M02(i))**2.d0)
         !
      enddo
!$OMP END PARALLEL DO
   endif
!
!..Print signal power output values
   if (FileIter .eq. -1) then
      write(*,*) ' get_power: printing power values (signal):'
      do i = 1,NumPts
         write(*,2020) thValues(i),sign_power(i)
      2020 format(f8.5,' ',es12.5)
      enddo
   elseif (FileIter .ge. 0) then
      !WRITE TO FILE
      write(*,*) ' get_power: printing power values (signal) to file..'
      fmt = '(I5.5)'
      write (suffix,fmt) FileIter
      filename=trim(OUTPUT_DIR)//'power/signal_'//trim(suffix)//'.dat'
      open(UNIT=9,FILE=filename,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
      do i = 1,NumPts
         write(UNIT=9, FMT=2020)  thValues(i), sign_power(i)
      enddo
      close(UNIT=9)
   endif
!
!..Print fiber core power ratio
   if (FileIter .eq. -1) then
      write(*,*) ' get_power: printing fiber core power ratio (signal):'
      do i = 1,NumPts
         write(*,2030) thValues(i), core_power(i)/sign_power(i)
      2030 format(f8.5,' ',f8.4)
      enddo
   elseif (FileIter .ge. 0) then
      !WRITE TO FILE
      write(*,*) ' get_power: printing fiber core power ratio (signal) to file..'
      fmt = '(I5.5)'
      write (suffix,fmt) FileIter
      filename=trim(OUTPUT_DIR)//'power/ratio_'//trim(suffix)//'.dat'
      open(UNIT=9,FILE=filename,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
      do i = 1,NumPts
         write(UNIT=9, FMT=2030)  thValues(i), core_power(i)/sign_power(i)
      enddo
      close(UNIT=9)
   endif
!
!..get power_loss_exp
   allocate(power_loss_exp(NumPts))
   write(*,*) ' get_power: computing power_loss_exp..'
   power_loss_exp(1) = 0.d0
   do i = 2,NumPts
      power_loss_exp(i) = log(sign_power(i)/sign_power(1)) / (2.d0*thValues(i))
   enddo
!..assign values to output variables
   InputPower = sign_power(1)
   FinalPower = sign_power(NumPts)
   FinalLossExp = power_loss_exp(NumPts)
   
   if (FileIter .eq. -1) then
      write(*,*) ' get_power: printing power_loss_exp:'
      do i = 1,NumPts
!  2040   format('    ',es12.5)
         write(*,2020) thValues(i), power_loss_exp(i)
      enddo
   elseif (FileIter .ge. 0) then
      !WRITE TO FILE
      write(*,*) ' get_power: printing power_loss_exp to file..'
      fmt = '(I5.5)'
      write (suffix,fmt) FileIter
      filename=trim(OUTPUT_DIR)//'power/power_loss_exp_'//trim(suffix)//'.dat'
      open(UNIT=9,FILE=filename,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
      do i = 1,NumPts
         write(UNIT=9, FMT=2020)  thValues(i), power_loss_exp(i)
      enddo
      close(UNIT=9)
   endif
   deallocate(power_loss_exp)
!
!..Print mode power output values
   if (ModeProj) then
      if (FileIter .eq. -1) then
         write(*,*)
         write(*,*) ' get_power: printing M00 (x) mode power (signal):'
         do i = 1,NumPts
!            write(*,*) 'norm_M00: ', norm_M00(i), ', coef_M00: ', coef_M00_r(i)
!            write(*,2020) power_M00(i)
            write(*,2020) thValues(i), power_M00(i)/sign_power(i)
         enddo
         write(*,*)
         write(*,*) ' get_power: printing M01 (x) mode power (signal):'
         do i = 1,NumPts
!            write(*,*) 'norm_M01: ', norm_M01(i), ', coef_M01: ', coef_M01_r(i)
!            write(*,2020) power_M01(i)
            write(*,2020) thValues(i), power_M01(i)/sign_power(i)
         enddo
         write(*,*)
         write(*,*) ' get_power: printing M02 (x) mode power (signal):'
         do i = 1,NumPts
!            write(*,*) 'norm_M02: ', norm_M02(i), ', coef_M02: ', coef_M02_r(i)
!            write(*,2020) power_M02(i)
            write(*,2020) thValues(i), power_M02(i)/sign_power(i)
         enddo
!
      elseif (FileIter .ge. 0) then
         !WRITE TO FILE
         write(*,*) ' get_power: printing M00 (x) mode power (signal) to file..'
         fmt = '(I5.5)'
         write (suffix,fmt) FileIter
         filename=trim(OUTPUT_DIR)//'power/powerM00_'//trim(suffix)//'.dat'
         open(UNIT=9,FILE=filename,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
         do i = 1,NumPts
            write(UNIT=9, FMT=2020)  thValues(i), power_M00(i)
         enddo
         close(UNIT=9)
         !WRITE TO FILE
         write(*,*) ' get_power: printing M01 (x) mode power (signal) to file..'
         fmt = '(I5.5)'
         write (suffix,fmt) FileIter
         filename=trim(OUTPUT_DIR)//'power/powerM01_'//trim(suffix)//'.dat'
         open(UNIT=9,FILE=filename,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
         do i = 1,NumPts
            write(UNIT=9, FMT=2020)  thValues(i), power_M01(i)
         enddo
         close(UNIT=9)
         !WRITE TO FILE
         write(*,*) ' get_power: printing M02 (x) mode power (signal) to file..'
         fmt = '(I5.5)'
         write (suffix,fmt) FileIter
         filename=trim(OUTPUT_DIR)//'power/powerM02_'//trim(suffix)//'.dat'
         open(UNIT=9,FILE=filename,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
         do i = 1,NumPts
            write(UNIT=9, FMT=2020)  thValues(i), power_M02(i)
         enddo
         close(UNIT=9)
      endif
   endif
!
   90 continue
   deallocate(thValues,sign_power,diff_power,core_power,clad_power)
!
   if (ModeProj) then
      deallocate(power_M00,norm_M00,coef_M00_r,coef_M00_c)
      deallocate(power_M01,norm_M01,coef_M01_r,coef_M01_c)
      deallocate(power_M02,norm_M02,coef_M02_r,coef_M02_c)
!
   endif
!
   99 continue
!
end subroutine get_power
!
!
!----------------------------------------------------------------------
!
!   routine name       - compute_power_bent_slab
!
!----------------------------------------------------------------------
!
!   latest revision    - Dec 2025
!
!   purpose            - Evaluates the electric field power of UW
!                        Maxwell along the cross sections specified by
!                        the vector of thValues in the input
!
!   arguments
!        in:
!                      - ThValues    : sample points in theta-direction
!                      - Num_thpts   : number of sample points
!                      - Fld         : 1   signal power
!                                      300 M00 Mode projection
!                                      301 M01 Mode projection
!                                      302 M02 Mode projection
!       out:
!                      - Power       : Absolute value of power at each theta-point
!                      - DiffPower   : Diff exact to computed power (alt: Norm)
!                                      (available if NEXACT>0)
!                      - CorePower   : Power in core subdomains (alt: Coef_r)
!                      - CladPower   : Power in cladding subdomains (alt: Coef_c)
!
!----------------------------------------------------------------------
!
subroutine compute_power_bent_slab(ThValues,Num_thpts,Fld, Power,DiffPower,CorePower,CladPower)
!
   use commonParam
   use data_structure3D
   use environment, only : QUIET_MODE
   use mpi_wrapper
   use par_mesh   , only : DISTRIBUTED
!
   implicit none
!
   integer, intent(in)  :: Num_thpts
   real(8), intent(in)  :: ThValues(Num_thpts)
   integer, intent(in)  :: Fld
   real(8), intent(out) :: Power(Num_thpts)
   real(8), intent(out) :: DiffPower(Num_thpts)
   real(8), intent(out) :: CorePower(Num_thpts)
   real(8), intent(out) :: CladPower(Num_thpts)
!
!..auxiliary variables
   real(8)    :: facePower, faceDiffPower
   real(8)    :: modeNorm
   complex(8) :: modeCoef
!
!..mdle number
   integer :: mdle
!
!..element, face order, geometry dof
   real(8) :: xnod(3,8),xnod_th(8)
   real(8) :: maxth,minth
!
!..miscellanea
   integer :: iel, i, ndom, nv, iv
!
!..face number over which power is computed
!  (in brick and prism, face 1 is face normal to xi3, at xi3=0)
   integer, parameter :: faceNum = 1
!
!..timer
   real(8) :: start_time,end_time
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
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!
   if (.not. DISTRIBUTED) then
      ELEM_SUBD(1:NRELES) = ELEM_ORDER(1:NRELES)
      NRELES_SUBD = NRELES
   endif
!
!..iterate over elements
!
!$OMP PARALLEL DO                                        &
!$OMP PRIVATE(mdle,nv,xnod,maxth,minth,i,ndom,           &
!$OMP         facePower,faceDiffPower,modeNorm,modeCoef) &
!$OMP REDUCTION(+:Power,DiffPower,corePower,cladPower)   &
!$OMP SCHEDULE(DYNAMIC)
   do iel=1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
      call find_domain(mdle, ndom)
      call nodcor_vert(mdle, xnod)
      nv = nvert(NODES(mdle)%ntype)
      do iv=1,nv
         xnod_th(iv) = atan2(xnod(3,iv),xnod(2,iv))
      enddo
      maxth = maxval(xnod_th(1:nv))
      minth = minval(xnod_th(1:nv))
      do i=1,Num_thpts
         if((ThValues(i).le.maxth).and.(ThValues(i).gt.minth)) then
            ! write(*,*) ' compute_power_bent_slab: element ',mdle,' contains thValue ',ThValues(i), &
            !            ' (minth=',minth,', maxth=',maxth,')'
            if (Fld .eq. 1) then
               call compute_facePower(mdle,faceNum, facePower,faceDiffPower)
               DiffPower(i) = DiffPower(i) + abs(faceDiffPower)
               select case(ndom)
                  case(1);   CorePower(i) = CorePower(i) + abs(facePower)
                  case(2,3); CladPower(i) = CladPower(i) + abs(facePower)
               end select
            elseif (Fld.eq.300 .or. Fld.eq.301 .or. Fld.eq.302) then
               call compute_mode_power(mdle,faceNum, facePower,modeNorm,modeCoef)
               DiffPower(i) = DiffPower(i) + modeNorm       ! false name (calc norm)
               CorePower(i) = CorePower(i) + real(modeCoef) ! false name (calc coef_r)
               CladPower(i) = CladPower(i) + imag(modeCoef) ! false name (calc coef_c)
            endif
            Power(i) = Power(i) + abs(facePower)
         endif
      enddo
   enddo
!$OMP END PARALLEL DO
!
   90 continue
!
!..end timer
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
   if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) then
      write(*,3010) end_time-start_time
 3010 format('  compute_power : ',f12.5,'  seconds')
   endif
!
end subroutine compute_power_bent_slab
!
!
!----------------------------------------------------------------------
!
!   routine name       - compute_face_power
!
!----------------------------------------------------------------------
!
!   latest revision    - Oct 2019
!
!   purpose            - Evaluates the electric field power of UW
!                        Maxwell by integrating H(curl) trace solution
!                        on a face of a middle node.
!
!   arguments
!        in:
!                      - Mdle       : middle element node
!                      - Facenumber : element face used for integration
!       out:
!                      - FacePower     :
!                      - FaceDiffPower :
!
!----------------------------------------------------------------------
!
subroutine compute_facePower(Mdle,Facenumber,FacePower,FaceDiffPower)
!
   use control
   use data_structure3D
   use physics
   use parametersDPG
   use commonParam
!
   implicit none
!
   integer, intent(in)  :: Mdle
   integer, intent(in)  :: Facenumber
   real(8), intent(out) :: FacePower
   real(8), intent(out) :: FaceDiffPower
!
!..element, face order, geometry dof
   integer, dimension(19)          :: norder
   real(8), dimension(3,MAXbrickH) :: xnod
   integer, dimension(12)          :: nedge_orient
   integer, dimension(6)           :: nface_orient
!
!..face order
   integer :: norderf(5)
!
!..number of vertices,edge,faces per element type
   integer :: nrv, nre, nrf
!
!..declare edge/face type varibles
   integer :: etype,ftype
!
!..variables for geometry
   real(8), dimension(3)   :: xi,x,rn
   real(8), dimension(3,2) :: dxidt,dxdt
   real(8), dimension(3,3) :: dxdxi,dxidx
   real(8), dimension(2)   :: t
   real(8)                 :: rjac,bjac
!
!..2D quadrature data
   real(8), dimension(2,MAXNINT2ADD) :: tloc
   real(8), dimension(MAXNINT2ADD)   :: wtloc
!
!..approximate solution dof's
   VTYPE, dimension(MAXEQNH,MAXbrickH) :: zdofH
   VTYPE, dimension(MAXEQNE,MAXbrickE) :: zdofE
   VTYPE, dimension(MAXEQNV,MAXbrickV) :: zdofV
   VTYPE, dimension(MAXEQNQ,MAXbrickQ) :: zdofQ
!..H1 shape functions
   integer                         :: nrdofH
   real(8), dimension(MAXbrickH)   :: shapH
   real(8), dimension(3,MAXbrickH) :: gradH
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
   VTYPE,dimension(  MAXEQNH    ) ::   ValH
   VTYPE,dimension(  MAXEQNH,3  ) ::  DvalH
   VTYPE,dimension(  MAXEQNH,3,3) :: d2valH
   VTYPE,dimension(3,MAXEQNE    ) ::   ValE
   VTYPE,dimension(3,MAXEQNE,3  ) ::  DvalE
   VTYPE,dimension(3,MAXEQNE,3,3) :: d2valE
   VTYPE,dimension(3,MAXEQNV    ) ::   ValV
   VTYPE,dimension(3,MAXEQNV,3  ) ::  DvalV
!
!..exact solution (UNUSED)
   VTYPE,dimension(3,MAXEQNV,3,3) :: d2valV
   VTYPE,dimension(  MAXEQNQ    ) ::   valQ
   VTYPE,dimension(  MAXEQNQ,3  ) ::  dvalQ
   VTYPE,dimension(  MAXEQNQ,3,3) :: d2valQ
!
!..for Poynting vector
   VTYPE, dimension(3) :: EtimesH1,EtimesH2
   VTYPE               :: FdotN
!
!..miscellanea
   integer :: nint,l
   real(8) :: weight
   integer :: nsign
   integer :: nflag
!
!---------------------------------------------------------------------------------------
!
   FacePower = 0.d0
   FaceDiffPower = 0.0d0
   nflag = 1
!..element type
   etype = NODES(Mdle)%ntype
   nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
   call find_order(Mdle, norder)
   call find_orient(Mdle, nedge_orient,nface_orient)
   call nodcor(Mdle, xnod)
   call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
!..sign factor to determine the OUTWARD normal unit vector
   nsign = nsign_param(etype,Facenumber)
!
!..face type
   ftype = face_type(etype,Facenumber)
!
!..face order of approximation
   call face_order(etype,Facenumber,norder, norderf)
!
!..set 2D quadrature
   INTEGRATION = NORD_ADD ! why ?
   call set_2D_int(ftype,norderf,nface_orient(Facenumber), nint,tloc,wtloc)
   INTEGRATION = 0
!
!..loop over integration points
   do l=1,nint
!
!  ...face coordinates
      t(1:2) = tloc(1:2,l)
!
!  ...face parametrization
      call face_param(etype,Facenumber,t, xi,dxidt)
!
!  ...determine element H1 shape functions (for geometry)
      call shape3DH(etype,xi,norder,nedge_orient,nface_orient, &
                    nrdofH,shapH,gradH)
!
!  ...geometry
      call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign, &
                   x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
      weight = bjac*wtloc(l)
!
      call soleval(Mdle,xi,nedge_orient,nface_orient,norder,xnod, &
                   zdofH,zdofE,zdofV,zdofQ,nflag,x,dxdxi, &
                   zsolH,zdsolH,zsolE,zcurlE,zsolV,zdivV,zsolQ)
      if(NEXACT.ne.0) then
         call exact(x,Mdle, ValH,DvalH,d2valH, ValE,DvalE,d2valE, &
                            ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      endif

      ! write(*,*) ' compute_facePower: Mdle=',Mdle, &
      !              ', x=(',x(1),',',x(2),',',x(3),')', &
      !              'theta=',atan2(x(3),x(2))
                   !
!     accumulate Poynting vector power for signal,
!     i.e., integrate (Real(n \dot ExH^*))
      ! from traces
      call zz_cross_product(zsolE(1:3,1),conjg(zsolE(1:3,2)), EtimesH1)
      ! ! from fields
      ! call zz_cross_product(zsolQ(1:3),conjg(zsolQ(4:6)), EtimesH1)

      FdotN = EtimesH1(1)*rn(1)+EtimesH1(2)*rn(2)+EtimesH1(3)*rn(3)
      FacePower = FacePower + (real(FdotN))*weight

      ! FacePower = FacePower + sqrt(real(EtimesH1(1)*conjg(EtimesH1(1))+          &
      !                                   EtimesH1(2)*conjg(EtimesH1(2))+          &
      !                                   EtimesH1(3)*conjg(EtimesH1(3)) ))*weight
!     ...if we have an exact
      if(NEXACT.ne.0) then
         call zz_cross_product(valE(1:3,1),conjg(valE(1:3,2)), EtimesH2)
         FaceDiffPower = FaceDiffPower   &
                        + abs(((EtimesH1(1)*rn(1)+EtimesH1(2)*rn(2)+EtimesH1(3)*rn(3))*weight) - &
                              ((EtimesH2(1)*rn(1)+EtimesH2(2)*rn(2)+EtimesH2(3)*rn(3))*weight))
      endif
!..end loop over integration points
   enddo
!
end subroutine compute_facePower
!
!
!..purpose:
!

!
!
!----------------------------------------------------------------------
!
!   routine name       - compute_mode_power
!
!----------------------------------------------------------------------
!
!   latest revision    - Oct 2019
!
!   purpose            - Compute projection of field onto bent 3-layer waveguide modes
!
!   arguments
!        in:
!                      - Mdle       : middle element node
!                      - Facenumber : element face used for integration
!       out:
!                      - ModeNorm   : norm of the mode (for normalization)
!                      - ModeCoef   : coefficient in the projection on mode
!
!----------------------------------------------------------------------
subroutine compute_mode_power(Mdle,Facenumber,ModePower,ModeNorm,ModeCoef)
!
   use control
   use data_structure3D
   use physics
   use parametersDPG
   use commonParam
!
   implicit none
!
   integer   , intent(in)  :: Mdle
   integer   , intent(in)  :: Facenumber
   real(8)   , intent(out) :: ModePower
   real(8)   , intent(out) :: ModeNorm
   complex(8), intent(out) :: ModeCoef
!
!..element, face order, geometry dof
   integer,dimension(19)          :: norder
   real(8),dimension(3,MAXbrickH) :: xnod
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
   integer :: etype,ftype
!
!..variables for geometry
   real(8), dimension(3)   :: xi,x,rn
   real(8), dimension(3,2) :: dxidt,dxdt
   real(8), dimension(3,3) :: dxdxi,dxidx
   real(8), dimension(2)   :: t
   real(8)                 :: rjac,bjac
!
!..2D quadrature data
   real(8), dimension(2,MAXNINT2ADD) :: tloc
   real(8), dimension(MAXNINT2ADD)   :: wtloc
!
!..approximate solution dof's
   VTYPE, dimension(MAXEQNH,MAXbrickH) :: zdofH
   VTYPE, dimension(MAXEQNE,MAXbrickE) :: zdofE
   VTYPE, dimension(MAXEQNV,MAXbrickV) :: zdofV
   VTYPE, dimension(MAXEQNQ,MAXbrickQ) :: zdofQ
!..H1 shape functions
   integer                         :: nrdofH
   real(8), dimension(MAXbrickH)   :: shapH
   real(8), dimension(3,MAXbrickH) :: gradH
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
   VTYPE,dimension(  MAXEQNH    ) ::   ValH
   VTYPE,dimension(  MAXEQNH,3  ) ::  DvalH
   VTYPE,dimension(  MAXEQNH,3,3) :: d2valH
   VTYPE,dimension(3,MAXEQNE    ) ::   ValE
   VTYPE,dimension(3,MAXEQNE,3  ) ::  DvalE
   VTYPE,dimension(3,MAXEQNE,3,3) :: d2valE
   VTYPE,dimension(3,MAXEQNV    ) ::   ValV
   VTYPE,dimension(3,MAXEQNV,3  ) ::  DvalV
!
!..exact solution (UNUSED)
   VTYPE,dimension(3,MAXEQNV,3,3) :: d2valV
   VTYPE,dimension(  MAXEQNQ    ) ::   valQ
   VTYPE,dimension(  MAXEQNQ,3  ) ::  dvalQ
   VTYPE,dimension(  MAXEQNQ,3,3) :: d2valQ
!
!..for Poynting vector
   VTYPE :: EtimesH(3)
   VTYPE :: FdotN
!
!..miscellanea
   integer :: nint,l
   real(8) :: weight
   integer :: nsign
   integer :: nflag
!
!---------------------------------------------------------------------------------------
!
   ModePower = 0.d0
   nflag = 1
!..element type
   etype = NODES(Mdle)%ntype
   nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
   call find_order(Mdle, norder)
   call find_orient(Mdle, nedge_orient,nface_orient)
   call nodcor(Mdle, xnod)
   call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
!..sign factor to determine the OUTWARD normal unit vector
   nsign = nsign_param(etype,Facenumber)
!
!..face type
   ftype = face_type(etype,Facenumber)
!
!..face order of approximation
   call face_order(etype,Facenumber,norder, norderf)
!
!..set 2D quadrature
   INTEGRATION = NORD_ADD ! why ?
   call set_2D_int(ftype,norderf,nface_orient(Facenumber), nint,tloc,wtloc)
   INTEGRATION = 0
!
!..first loop over integration points to find projection coefficients
   ModeNorm = 0.d0; ModeCoef = 0.d0
   do l=1,nint
!
!  ...face coordinates
      t(1:2) = tloc(1:2,l)
!
!  ...face parametrization
      call face_param(etype,Facenumber,t, xi,dxidt)
!
!  ...determine element H1 shape functions (for geometry)
      call shape3DH(etype,xi,norder,nedge_orient,nface_orient, &
                    nrdofH,shapH,gradH)
!
!  ...geometry
      call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign, &
                     x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
      weight = bjac*wtloc(l)
!
      call soleval(Mdle,xi,nedge_orient,nface_orient,norder,xnod, &
                   zdofH,zdofE,zdofV,zdofQ,nflag,x,dxdxi, &
                   zsolH,zdsolH,zsolE,zcurlE,zsolV,zdivV,zsolQ)
!
!  ...compute field of the mode
      call exact(x,Mdle, ValH,DvalH,d2valH, ValE,DvalE,d2valE, &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
!
!     accumulate L2 inner product (signal),
!     i.e., integrate (E     \dot phi_m^*) for m-th mode,
!       and integrate (phi_m \dot phi_m^*) for m-th mode
         ModeCoef = ModeCoef +     (zsolE(1,1) * conjg(valE(1,1)) +    &
                                    zsolE(2,1) * conjg(valE(2,1)) +    &
                                    zsolE(3,1) * conjg(valE(3,1))) * weight
         ModeNorm = ModeNorm + real (valE(1,1) * conjg(valE(1,1)) +   &
                                     valE(2,1) * conjg(valE(2,1)) +   &
                                     valE(3,1) * conjg(valE(3,1))) * weight
!..end loop over integration points
   enddo
!
!..second loop over integration points to calculate power of mode projection
   do l=1,nint
!
!  ...face coordinates
      t(1:2) = tloc(1:2,l)
!
!  ...face parametrization
      call face_param(etype,Facenumber,t, xi,dxidt)
!
!  ...determine element H1 shape functions (for geometry)
      call shape3DH(etype,xi,norder,nedge_orient,nface_orient, &
                    nrdofH,shapH,gradH)
!
!  ...geometry
      call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign, &
                     x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
      weight = bjac*wtloc(l)
!
!  ...compute field of the mode
      call exact(x,Mdle, ValH,DvalH,d2valH, ValE,DvalE,d2valE, &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
!
!     accumulate Poynting vector power for mode in signal
!     i.e., integrate [Real{n \dot (phi_m x conjg(beta*phi_m))}]
      call zz_cross_product(ValE(1:3,1),conjg(ValE(1:3,2)), EtimesH)
      FdotN = EtimesH(1)*rn(1)+EtimesH(2)*rn(2)+EtimesH(3)*rn(3)
      ModePower = ModePower + weight*real(FdotN)
!
!..end loop over integration points
   enddo
!
end subroutine compute_mode_power



