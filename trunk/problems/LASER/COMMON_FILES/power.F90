!
#include "typedefs.h"
!
!----------------------------------------------------------------------
!
!   routine name       - get_power
!
!----------------------------------------------------------------------
!
!   latest revision - Mar 2022
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
   use par_mesh , only: DISTRIBUTED,HOST_MESH
!
   implicit none
!
   integer, intent(in)    :: Fld
   integer, intent(in)    :: FileIter
   integer, intent(inout) :: NumPts
!
   real(8), allocatable :: zValues(:)
   real(8), allocatable :: sign_power(:),pump_power(:)
   real(8), allocatable :: diff_power(:),efficiency(:)
   real(8), allocatable :: core_power(:),clad_power(:)
!
   real(8), allocatable :: power_LP01_x(:),power_LP11a_x(:),power_LP11b_x(:)
   real(8), allocatable :: power_LP21a_x(:),power_LP21b_x(:),power_LP02_x(:)
   real(8), allocatable :: norm_LP01_x(:),norm_LP11a_x(:),norm_LP11b_x(:)
   real(8), allocatable :: norm_LP21a_x(:),norm_LP21b_x(:),norm_LP02_x(:)
   real(8), allocatable :: coef_LP01_r_x(:),coef_LP11a_r_x(:),coef_LP11b_r_x(:)
   real(8), allocatable :: coef_LP21a_r_x(:),coef_LP21b_r_x(:),coef_LP02_r_x(:)
   real(8), allocatable :: coef_LP01_c_x(:),coef_LP11a_c_x(:),coef_LP11b_c_x(:)
   real(8), allocatable :: coef_LP21a_c_x(:),coef_LP21b_c_x(:),coef_LP02_c_x(:)
!
   real(8), allocatable :: power_LP01_y(:),power_LP11a_y(:),power_LP11b_y(:)
   real(8), allocatable :: power_LP21a_y(:),power_LP21b_y(:),power_LP02_y(:)
   real(8), allocatable :: norm_LP01_y(:),norm_LP11a_y(:),norm_LP11b_y(:)
   real(8), allocatable :: norm_LP21a_y(:),norm_LP21b_y(:),norm_LP02_y(:)
   real(8), allocatable :: coef_LP01_r_y(:),coef_LP11a_r_y(:),coef_LP11b_r_y(:)
   real(8), allocatable :: coef_LP21a_r_y(:),coef_LP21b_r_y(:),coef_LP02_r_y(:)
   real(8), allocatable :: coef_LP01_c_y(:),coef_LP11a_c_y(:),coef_LP11b_c_y(:)
   real(8), allocatable :: coef_LP21a_c_y(:),coef_LP21b_c_y(:),coef_LP02_c_y(:)
!
   real(8) :: a,b,gain,loss
   integer :: i,j
!
   character(8)  :: fmt,suffix
   character(64) :: filename
!
   integer :: count,ierr
!
!..activate to calculate power in each signal LP mode (by projection)
   logical, parameter :: modeProj = .false.
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
   if ((.not. DISTRIBUTED .or. HOST_MESH) .and. RANK .ne. ROOT) goto 99
!
   allocate(zValues(NumPts)   , sign_power(NumPts), &
            pump_power(NumPts), diff_power(NumPts), &
            core_power(NumPts), clad_power(NumPts)  )
!
   if (modeProj) then
   allocate(power_LP01_x(NumPts),norm_LP01_x(NumPts),coef_LP01_r_x(NumPts),coef_LP01_c_x(NumPts))
   allocate(power_LP11a_x(NumPts),norm_LP11a_x(NumPts),coef_LP11a_r_x(NumPts),coef_LP11a_c_x(NumPts))
   allocate(power_LP11b_x(NumPts),norm_LP11b_x(NumPts),coef_LP11b_r_x(NumPts),coef_LP11b_c_x(NumPts))
   allocate(power_LP21a_x(NumPts),norm_LP21a_x(NumPts),coef_LP21a_r_x(NumPts),coef_LP21a_c_x(NumPts))
   allocate(power_LP21b_x(NumPts),norm_LP21b_x(NumPts),coef_LP21b_r_x(NumPts),coef_LP21b_c_x(NumPts))
   allocate(power_LP02_x(NumPts),norm_LP02_x(NumPts),coef_LP02_r_x(NumPts),coef_LP02_c_x(NumPts))
!
   allocate(power_LP01_y(NumPts),norm_LP01_y(NumPts),coef_LP01_r_y(NumPts),coef_LP01_c_y(NumPts))
   allocate(power_LP11a_y(NumPts),norm_LP11a_y(NumPts),coef_LP11a_r_y(NumPts),coef_LP11a_c_y(NumPts))
   allocate(power_LP11b_y(NumPts),norm_LP11b_y(NumPts),coef_LP11b_r_y(NumPts),coef_LP11b_c_y(NumPts))
   allocate(power_LP21a_y(NumPts),norm_LP21a_y(NumPts),coef_LP21a_r_y(NumPts),coef_LP21a_c_y(NumPts))
   allocate(power_LP21b_y(NumPts),norm_LP21b_y(NumPts),coef_LP21b_r_y(NumPts),coef_LP21b_c_y(NumPts))
   allocate(power_LP02_y(NumPts),norm_LP02_y(NumPts),coef_LP02_r_y(NumPts),coef_LP02_c_y(NumPts))
   endif
!
!..distributing sample points uniformly
   if (RANK .eq. ROOT) then
      write(*,*) ' get_power: Distributing sample points uniformly along waveguide.'
      write(*,2002) ' ZL = ', ZL
 2002 format(A,F10.2,/)
   endif
   b = ZL/NumPts
   a = b/2.d0
   do i=1,NumPts
      zValues(i) = (i-1)*b+a
   enddo
!..irrationalize z values to avoid points on element interfaces
   zValues = zValues*PI*(113.d0/355.d0)
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
         stop
   end select
!
   if (modeProj) then
   if (RANK.eq.ROOT) write(*,*) ' get_power: computing signal mode_power..'
   i = ISOL; j = ICOMP_EXACT
   ISOL = 13; ICOMP_EXACT = 1 ! LP01 projection (x-polarized)
   call compute_power(zValues,NumPts,ISOL, power_LP01_x,norm_LP01_x,coef_LP01_r_x,coef_LP01_c_x)
   ISOL = 13; ICOMP_EXACT = 2 ! LP01 projection (y-polarized)
   call compute_power(zValues,NumPts,ISOL, power_LP01_y,norm_LP01_y,coef_LP01_r_y,coef_LP01_c_y)
!
   ISOL = 14 ; ICOMP_EXACT = 1 ! LP11a projection (x-polarized)
   call compute_power(zValues,NumPts,ISOL, power_LP11a_x,norm_LP11a_x,coef_LP11a_r_x,coef_LP11a_c_x)
   ISOL = 140; ICOMP_EXACT = 1 ! LP11b projection (x-polarized)
   call compute_power(zValues,NumPts,ISOL, power_LP11b_x,norm_LP11b_x,coef_LP11b_r_x,coef_LP11b_c_x)
   ISOL = 14 ; ICOMP_EXACT = 2 ! LP11a projection (y-polarized)
   call compute_power(zValues,NumPts,ISOL, power_LP11a_y,norm_LP11a_y,coef_LP11a_r_y,coef_LP11a_c_y)
   ISOL = 140; ICOMP_EXACT = 2 ! LP11b projection (y-polarized)
   call compute_power(zValues,NumPts,ISOL, power_LP11b_y,norm_LP11b_y,coef_LP11b_r_y,coef_LP11b_c_y)
!
   ISOL = 15 ; ICOMP_EXACT = 1 ! LP21a projection (x-polarized)
   call compute_power(zValues,NumPts,ISOL, power_LP21a_x,norm_LP21a_x,coef_LP21a_r_x,coef_LP21a_c_x)
   ISOL = 150; ICOMP_EXACT = 1 ! LP21b projection (x-polarized)
   call compute_power(zValues,NumPts,ISOL, power_LP21b_x,norm_LP21b_x,coef_LP21b_r_x,coef_LP21b_c_x)
   ISOL = 15 ; ICOMP_EXACT = 2 ! LP21a projection (y-polarized)
   call compute_power(zValues,NumPts,ISOL, power_LP21a_y,norm_LP21a_y,coef_LP21a_r_y,coef_LP21a_c_y)
   ISOL = 150; ICOMP_EXACT = 2 ! LP21b projection (y-polarized)
   call compute_power(zValues,NumPts,ISOL, power_LP21b_y,norm_LP21b_y,coef_LP21b_r_y,coef_LP21b_c_y)
!
   ISOL = 16; ICOMP_EXACT = 1 ! LP02 projection (x-polarized)
   call compute_power(zValues,NumPts,ISOL, power_LP02_x,norm_LP02_x,coef_LP02_r_x,coef_LP02_c_x)
   ISOL = 16; ICOMP_EXACT = 2 ! LP02 projection (y-polarized)
   call compute_power(zValues,NumPts,ISOL, power_LP02_y,norm_LP02_y,coef_LP02_r_y,coef_LP02_c_y)
   ISOL = i; ICOMP_EXACT = j
   endif
!
!..gather all values on host
   if (.not. DISTRIBUTED .or. HOST_MESH) goto 50
   if (PLANE_PUMP.ne.0 .and. RANK.ne.ROOT) then
      pump_power = 0.d0
   endif
   count = NumPts
   if (RANK .eq. ROOT) then
      call MPI_REDUCE(MPI_IN_PLACE,sign_power,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE,pump_power,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE,diff_power,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE,core_power,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE,clad_power,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      !
      if (modeProj) then
      call MPI_REDUCE(MPI_IN_PLACE,power_LP01_x  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, norm_LP01_x  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, coef_LP01_r_x,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, coef_LP01_c_x,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      !
      call MPI_REDUCE(MPI_IN_PLACE,power_LP01_y  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, norm_LP01_y  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, coef_LP01_r_y,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, coef_LP01_c_y,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      !
      call MPI_REDUCE(MPI_IN_PLACE,power_LP11a_x  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, norm_LP11a_x  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, coef_LP11a_r_x,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, coef_LP11a_c_x,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      !
      call MPI_REDUCE(MPI_IN_PLACE,power_LP11b_x  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, norm_LP11b_x  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, coef_LP11b_r_x,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, coef_LP11b_c_x,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      !
      call MPI_REDUCE(MPI_IN_PLACE,power_LP11a_y  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, norm_LP11a_y  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, coef_LP11a_r_y,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, coef_LP11a_c_y,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      !
      call MPI_REDUCE(MPI_IN_PLACE,power_LP11b_y  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, norm_LP11b_y  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, coef_LP11b_r_y,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, coef_LP11b_c_y,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      !
      call MPI_REDUCE(MPI_IN_PLACE,power_LP21a_x  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, norm_LP21a_x  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, coef_LP21a_r_x,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, coef_LP21a_c_x,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      !
      call MPI_REDUCE(MPI_IN_PLACE,power_LP21b_x  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, norm_LP21b_x  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, coef_LP21b_r_x,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, coef_LP21b_c_x,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      !
      call MPI_REDUCE(MPI_IN_PLACE,power_LP21a_y  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, norm_LP21a_y  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, coef_LP21a_r_y,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, coef_LP21a_c_y,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      !
      call MPI_REDUCE(MPI_IN_PLACE,power_LP21b_y  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, norm_LP21b_y  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, coef_LP21b_r_y,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, coef_LP21b_c_y,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      !
      call MPI_REDUCE(MPI_IN_PLACE,power_LP02_x  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, norm_LP02_x  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, coef_LP02_r_x,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, coef_LP02_c_x,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      !
      call MPI_REDUCE(MPI_IN_PLACE,power_LP02_y  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, norm_LP02_y  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, coef_LP02_r_y,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(MPI_IN_PLACE, coef_LP02_c_y,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      endif
   else
      call MPI_REDUCE(sign_power,sign_power,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(pump_power,pump_power,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(diff_power,diff_power,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(core_power,core_power,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(clad_power,clad_power,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      !
      if (modeProj) then
      call MPI_REDUCE(power_LP01_x  ,power_LP01_x  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( norm_LP01_x  , norm_LP01_x  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( coef_LP01_r_x, coef_LP01_r_x,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( coef_LP01_c_x, coef_LP01_c_x,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      !
      call MPI_REDUCE(power_LP01_y  ,power_LP01_y  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( norm_LP01_y  , norm_LP01_y  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( coef_LP01_r_y, coef_LP01_r_y,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( coef_LP01_c_y, coef_LP01_c_y,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      !
      call MPI_REDUCE(power_LP11a_x  ,power_LP11a_x  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( norm_LP11a_x  , norm_LP11a_x  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( coef_LP11a_r_x, coef_LP11a_r_x,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( coef_LP11a_c_x, coef_LP11a_c_x,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      !
      call MPI_REDUCE(power_LP11b_x  ,power_LP11b_x  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( norm_LP11b_x  , norm_LP11b_x  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( coef_LP11b_r_x, coef_LP11b_r_x,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( coef_LP11b_c_x, coef_LP11b_c_x,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      !
      call MPI_REDUCE(power_LP11a_y  ,power_LP11a_y  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( norm_LP11a_y  , norm_LP11a_y  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( coef_LP11a_r_y, coef_LP11a_r_y,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( coef_LP11a_c_y, coef_LP11a_c_y,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      !
      call MPI_REDUCE(power_LP11b_y  ,power_LP11b_y  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( norm_LP11b_y  , norm_LP11b_y  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( coef_LP11b_r_y, coef_LP11b_r_y,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( coef_LP11b_c_y, coef_LP11b_c_y,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      !
      call MPI_REDUCE(power_LP21a_x  ,power_LP21a_x  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( norm_LP21a_x  , norm_LP21a_x  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( coef_LP21a_r_x, coef_LP21a_r_x,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( coef_LP21a_c_x, coef_LP21a_c_x,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      !
      call MPI_REDUCE(power_LP21b_x  ,power_LP21b_x  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( norm_LP21b_x  , norm_LP21b_x  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( coef_LP21b_r_x, coef_LP21b_r_x,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( coef_LP21b_c_x, coef_LP21b_c_x,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      !
      call MPI_REDUCE(power_LP21a_y  ,power_LP21a_y  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( norm_LP21a_y  , norm_LP21a_y  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( coef_LP21a_r_y, coef_LP21a_r_y,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( coef_LP21a_c_y, coef_LP21a_c_y,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      !
      call MPI_REDUCE(power_LP21b_y  ,power_LP21b_y  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( norm_LP21b_y  , norm_LP21b_y  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( coef_LP21b_r_y, coef_LP21b_r_y,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( coef_LP21b_c_y, coef_LP21b_c_y,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      !
      call MPI_REDUCE(power_LP02_x  ,power_LP02_x  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( norm_LP02_x  , norm_LP02_x  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( coef_LP02_r_x, coef_LP02_r_x,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( coef_LP02_c_x, coef_LP02_c_x,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      !
      call MPI_REDUCE(power_LP02_y  ,power_LP02_y  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( norm_LP02_y  , norm_LP02_y  ,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( coef_LP02_r_y, coef_LP02_r_y,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE( coef_LP02_c_y, coef_LP02_c_y,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      endif
      goto 90
   endif
!
   50 continue
!
   if (modeProj) then
!$OMP PARALLEL DO
   do i = 1,NumPts
      norm_LP01_x(i)   = sqrt(norm_LP01_x(i))
      coef_LP01_r_x(i) = sqrt(coef_LP01_r_x(i)**2.d0+coef_LP01_c_x(i)**2.d0) / norm_LP01_x(i)
      power_LP01_x(i)  = power_LP01_x(i) * ((coef_LP01_r_x(i) / norm_LP01_x(i))**2.d0)
      !
      norm_LP11a_x(i)   = sqrt(norm_LP11a_x(i))
      coef_LP11a_r_x(i) = sqrt(coef_LP11a_r_x(i)**2.d0+coef_LP11a_c_x(i)**2.d0) / norm_LP11a_x(i)
      power_LP11a_x(i)  = power_LP11a_x(i) * ((coef_LP11a_r_x(i) / norm_LP11a_x(i))**2.d0)
      !
      norm_LP11b_x(i)   = sqrt(norm_LP11b_x(i))
      coef_LP11b_r_x(i) = sqrt(coef_LP11b_r_x(i)**2.d0+coef_LP11b_c_x(i)**2.d0) / norm_LP11b_x(i)
      power_LP11b_x(i)  = power_LP11b_x(i) * ((coef_LP11b_r_x(i) / norm_LP11b_x(i))**2.d0)
      !
      norm_LP21a_x(i)   = sqrt(norm_LP21a_x(i))
      coef_LP21a_r_x(i) = sqrt(coef_LP21a_r_x(i)**2.d0+coef_LP21a_c_x(i)**2.d0) / norm_LP21a_x(i)
      power_LP21a_x(i)  = power_LP21a_x(i) * ((coef_LP21a_r_x(i) / norm_LP21a_x(i))**2.d0)
      !
      norm_LP21b_x(i)   = sqrt(norm_LP21b_x(i))
      coef_LP21b_r_x(i) = sqrt(coef_LP21b_r_x(i)**2.d0+coef_LP21b_c_x(i)**2.d0) / norm_LP21b_x(i)
      power_LP21b_x(i)  = power_LP21b_x(i) * ((coef_LP21b_r_x(i) / norm_LP21b_x(i))**2.d0)
      !
      norm_LP02_x(i)   = sqrt(norm_LP02_x(i))
      coef_LP02_r_x(i) = sqrt(coef_LP02_r_x(i)**2.d0+coef_LP02_c_x(i)**2.d0) / norm_LP02_x(i)
      power_LP02_x(i)  = power_LP02_x(i) * ((coef_LP02_r_x(i) / norm_LP02_x(i))**2.d0)
      !
      norm_LP01_y(i)   = sqrt(norm_LP01_y(i))
      coef_LP01_r_y(i) = sqrt(coef_LP01_r_y(i)**2.d0+coef_LP01_c_y(i)**2.d0) / norm_LP01_y(i)
      power_LP01_y(i)  = power_LP01_y(i) * ((coef_LP01_r_y(i) / norm_LP01_y(i))**2.d0)
      !
      norm_LP11a_y(i)   = sqrt(norm_LP11a_y(i))
      coef_LP11a_r_y(i) = sqrt(coef_LP11a_r_y(i)**2.d0+coef_LP11a_c_y(i)**2.d0) / norm_LP11a_y(i)
      power_LP11a_y(i)  = power_LP11a_y(i) * ((coef_LP11a_r_y(i) / norm_LP11a_y(i))**2.d0)
      !
      norm_LP11b_y(i)   = sqrt(norm_LP11b_y(i))
      coef_LP11b_r_y(i) = sqrt(coef_LP11b_r_y(i)**2.d0+coef_LP11b_c_y(i)**2.d0) / norm_LP11b_y(i)
      power_LP11b_y(i)  = power_LP11b_y(i) * ((coef_LP11b_r_y(i) / norm_LP11b_y(i))**2.d0)
      !
      norm_LP21a_y(i)   = sqrt(norm_LP21a_y(i))
      coef_LP21a_r_y(i) = sqrt(coef_LP21a_r_y(i)**2.d0+coef_LP21a_c_y(i)**2.d0) / norm_LP21a_y(i)
      power_LP21a_y(i)  = power_LP21a_y(i) * ((coef_LP21a_r_y(i) / norm_LP21a_y(i))**2.d0)
      !
      norm_LP21b_y(i)   = sqrt(norm_LP21b_y(i))
      coef_LP21b_r_y(i) = sqrt(coef_LP21b_r_y(i)**2.d0+coef_LP21b_c_y(i)**2.d0) / norm_LP21b_y(i)
      power_LP21b_y(i)  = power_LP21b_y(i) * ((coef_LP21b_r_y(i) / norm_LP21b_y(i))**2.d0)
      !
      norm_LP02_y(i)   = sqrt(norm_LP02_y(i))
      coef_LP02_r_y(i) = sqrt(coef_LP02_r_y(i)**2.d0+coef_LP02_c_y(i)**2.d0) / norm_LP02_y(i)
      power_LP02_y(i)  = power_LP02_y(i) * ((coef_LP02_r_y(i) / norm_LP02_y(i))**2.d0)
   enddo
!$OMP END PARALLEL DO
   endif
!
!..Print signal power output values
   if (Fld .eq. 1 .or. Fld .eq. 2) then
      if (FileIter .eq. -1) then
         write(*,*) ' get_power: printing power values (signal):'
         do i = 1,NumPts
            write(*,2020) sign_power(i)
       2020 format('    ',es12.5)
         enddo
      endif
      if (FileIter .ge. 0) then
         !WRITE TO FILE
         write(*,*) ' get_power: printing power values (signal) to file..'
         fmt = '(I5.5)'
         write (suffix,fmt) FileIter
         filename=trim(OUTPUT_DIR)//'power/signal_'//trim(suffix)//'.dat'
         open(UNIT=9,FILE=filename,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
         do i = 1,NumPts
            write(UNIT=9, FMT="(es12.5)") sign_power(i)
         enddo
         close(UNIT=9)
      endif
      if (NONLINEAR_FLAG .eq. 0) then
        i = NumPts
        if (USE_PML) then
           i = (1.0d0 - PML_FRAC) * NumPts
        endif
        write(*,2021) (sign_power(1)-sign_power(i))/sign_power(1) * 100.d0
   2021 format(' Power loss: ',f6.2,' %',/)
     endif
   endif
!
!..Print pump power output values
   if (Fld .eq. 0 .or. Fld .eq. 2) then
      if (FileIter .eq. -1) then
         write(*,*) ' get_power: printing power values (pump):'
         do i = 1,NumPts
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
         do i = 1,NumPts
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
         do i = 1,NumPts
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
         do i = 1,NumPts
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
      if (COPUMP .eq. 1) then
         efficiency(1) = 0.d0
         do i = 2,NumPts
            if (zValues(i) .gt. PML_REGION) then
               efficiency(i) = 0.d0
            else
               efficiency(i) = (sign_power(i)-sign_power(1)) &
                              /(pump_power(1)-pump_power(i))
            endif
         enddo
      elseif (COPUMP.eq.0) then
         ! define efficiency not point-wise but over the whole fiber
         ! compute signal gain and pump loss over amplifier (exclude pml)
         gain = maxval(sign_power) - sign_power(1) ! signal gain
         do i = 1,NumPts
            if (zValues(i).gt.(ZL-PML_REGION)) then
               loss = pump_power(NumPts) - pump_power(i) ! pump loss
               exit
            endif
         enddo
         efficiency(1:NumPts) = gain / loss
      else
         write(*,*) ' get_power: COPUMP must be 1 or 0. stop.'
         stop
      endif
      if (FileIter .eq. -1) then
         write(*,*) ' get_power: printing efficiency:'
         do i = 1,NumPts
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
         do i = 1,NumPts
            write(UNIT=9, FMT="(f8.4)") efficiency(i)
         enddo
         close(UNIT=9)
      endif
      deallocate(efficiency)
   endif
!
!..Print mode power output values
   if (modeProj) then
      if (FileIter .eq. -1) then
         write(*,*)
         write(*,*) ' get_power: printing LP01 (x) mode power (signal):'
         do i = 1,NumPts
!            write(*,*) 'norm_LP01_x: ', norm_LP01_x(i), ', coef_LP01_x: ', coef_LP01_r_x(i)
!            write(*,2020) power_LP01_x(i)
            write(*,2030) power_LP01_x(i)/sign_power(i)
         enddo
         write(*,*)
         write(*,*) ' get_power: printing LP11a (x) mode power (signal):'
         do i = 1,NumPts
!            write(*,*) 'norm_LP11a_x: ', norm_LP11a_x(i), ', coef_LP11a_x: ', coef_LP11a_r_x(i)
!            write(*,2020) power_LP11a_x(i)
            write(*,2030) power_LP11a_x(i)/sign_power(i)
         enddo
         write(*,*)
         write(*,*) ' get_power: printing LP11b (x) mode power (signal):'
         do i = 1,NumPts
!            write(*,*) 'norm_LP11b_x: ', norm_LP11b_x(i), ', coef_LP11b_x: ', coef_LP11b_r_x(i)
!            write(*,2020) power_LP11b_x(i)
            write(*,2030) power_LP11b_x(i)/sign_power(i)
         enddo
         write(*,*)
         write(*,*) ' get_power: printing LP21a (x) mode power (signal):'
         do i = 1,NumPts
!            write(*,*) 'norm_LP21a_x: ', norm_LP21a_x(i), ', coef_LP21a_x: ', coef_LP21a_r_x(i)
!            write(*,2020) power_LP21a_x(i)
            write(*,2030) power_LP21a_x(i)/sign_power(i)
         enddo
         write(*,*)
         write(*,*) ' get_power: printing LP21b (x) mode power (signal):'
         do i = 1,NumPts
!            write(*,*) 'norm_LP21b_x: ', norm_LP21b_x(i), ', coef_LP21b_x: ', coef_LP21b_r_x(i)
!            write(*,2020) power_LP21b_x(i)
            write(*,2030) power_LP21b_x(i)/sign_power(i)
         enddo
         write(*,*)
         write(*,*) ' get_power: printing LP02 (x) mode power (signal):'
         do i = 1,NumPts
!            write(*,*) 'norm_LP02_x: ', norm_LP02_x(i), ', coef_LP02_x: ', coef_LP02_r_x(i)
!            write(*,2020) power_LP02_x(i)
            write(*,2030) power_LP02_x(i)/sign_power(i)
         enddo
!
         write(*,*)
         write(*,*) ' get_power: printing LP01 (y) mode power (signal):'
         do i = 1,NumPts
!            write(*,*) 'norm_LP01_y: ', norm_LP01_y(i), ', coef_LP01_y: ', coef_LP01_r_y(i)
!            write(*,2020) power_LP01_y(i)
            write(*,2030) power_LP01_y(i)/sign_power(i)
         enddo
         write(*,*)
         write(*,*) ' get_power: printing LP11a (y) mode power (signal):'
         do i = 1,NumPts
!            write(*,*) 'norm_LP11a_y: ', norm_LP11a_y(i), ', coef_LP11a_y: ', coef_LP11a_r_y(i)
!            write(*,2020) power_LP11a_y(i)
            write(*,2030) power_LP11a_y(i)/sign_power(i)
         enddo
         write(*,*)
         write(*,*) ' get_power: printing LP11b (y) mode power (signal):'
         do i = 1,NumPts
!            write(*,*) 'norm_LP11b_y: ', norm_LP11b_y(i), ', coef_LP11b_y: ', coef_LP11b_r_y(i)
!            write(*,2020) power_LP11b_y(i)
            write(*,2030) power_LP11b_y(i)/sign_power(i)
         enddo
         write(*,*)
         write(*,*) ' get_power: printing LP21a (y) mode power (signal):'
         do i = 1,NumPts
!            write(*,*) 'norm_LP21a_y: ', norm_LP21a_y(i), ', coef_LP21a_y: ', coef_LP21a_r_y(i)
!            write(*,2020) power_LP21a(i)
            write(*,2030) power_LP21a_y(i)/sign_power(i)
         enddo
         write(*,*)
         write(*,*) ' get_power: printing LP21b (y) mode power (signal):'
         do i = 1,NumPts
!            write(*,*) 'norm_LP21b_y: ', norm_LP21b_y(i), ', coef_LP21b_y: ', coef_LP21b_r_y(i)
!            write(*,2020) power_LP21b(i)
            write(*,2030) power_LP21b_y(i)/sign_power(i)
         enddo
         write(*,*)
         write(*,*) ' get_power: printing LP02 (y) mode power (signal):'
         do i = 1,NumPts
!            write(*,*) 'norm_LP02_y: ', norm_LP02_y(i), ', coef_LP02_y: ', coef_LP02_r_y(i)
!            write(*,2020) power_LP02_y(i)
            write(*,2030) power_LP02_y(i)/sign_power(i)
         enddo
      endif
      if (FileIter .ge. 0) then
         !WRITE TO FILE
         write(*,*) ' get_power: printing LP01 (x) mode power (signal) to file..'
         fmt = '(I5.5)'
         write (suffix,fmt) FileIter
         filename=trim(OUTPUT_DIR)//'power/powerLP01_x_'//trim(suffix)//'.dat'
         open(UNIT=9,FILE=filename,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
         do i = 1,NumPts
            write(UNIT=9, FMT="(es12.5)") power_LP01_x(i)
         enddo
         close(UNIT=9)
         !WRITE TO FILE
         write(*,*) ' get_power: printing LP11a (x) mode power (signal) to file..'
         fmt = '(I5.5)'
         write (suffix,fmt) FileIter
         filename=trim(OUTPUT_DIR)//'power/powerLP11a_x_'//trim(suffix)//'.dat'
         open(UNIT=9,FILE=filename,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
         do i = 1,NumPts
            write(UNIT=9, FMT="(es12.5)") power_LP11a_x(i)
         enddo
         close(UNIT=9)
         !WRITE TO FILE
         write(*,*) ' get_power: printing LP11b (x) mode power (signal) to file..'
         fmt = '(I5.5)'
         write (suffix,fmt) FileIter
         filename=trim(OUTPUT_DIR)//'power/powerLP11b_x_'//trim(suffix)//'.dat'
         open(UNIT=9,FILE=filename,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
         do i = 1,NumPts
            write(UNIT=9, FMT="(es12.5)") power_LP11b_x(i)
         enddo
         close(UNIT=9)
         !WRITE TO FILE
         write(*,*) ' get_power: printing LP21a (x) mode power (signal) to file..'
         fmt = '(I5.5)'
         write (suffix,fmt) FileIter
         filename=trim(OUTPUT_DIR)//'power/powerLP21a_x_'//trim(suffix)//'.dat'
         open(UNIT=9,FILE=filename,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
         do i = 1,NumPts
            write(UNIT=9, FMT="(es12.5)") power_LP21a_x(i)
         enddo
         close(UNIT=9)
         !WRITE TO FILE
         write(*,*) ' get_power: printing LP21b (x) mode power (signal) to file..'
         fmt = '(I5.5)'
         write (suffix,fmt) FileIter
         filename=trim(OUTPUT_DIR)//'power/powerLP21b_x_'//trim(suffix)//'.dat'
         open(UNIT=9,FILE=filename,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
         do i = 1,NumPts
            write(UNIT=9, FMT="(es12.5)") power_LP21b_x(i)
         enddo
         close(UNIT=9)
         !WRITE TO FILE
         write(*,*) ' get_power: printing LP02 (x) mode power (signal) to file..'
         fmt = '(I5.5)'
         write (suffix,fmt) FileIter
         filename=trim(OUTPUT_DIR)//'power/powerLP02_x_'//trim(suffix)//'.dat'
         open(UNIT=9,FILE=filename,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
         do i = 1,NumPts
            write(UNIT=9, FMT="(es12.5)") power_LP02_x(i)
         enddo
         close(UNIT=9)
         !WRITE TO FILE
         write(*,*) ' get_power: printing LP01 (y) mode power (signal) to file..'
         fmt = '(I5.5)'
         write (suffix,fmt) FileIter
         filename=trim(OUTPUT_DIR)//'power/powerLP01_y_'//trim(suffix)//'.dat'
         open(UNIT=9,FILE=filename,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
         do i = 1,NumPts
            write(UNIT=9, FMT="(es12.5)") power_LP01_y(i)
         enddo
         close(UNIT=9)
         !WRITE TO FILE
         write(*,*) ' get_power: printing LP11a (y) mode power (signal) to file..'
         fmt = '(I5.5)'
         write (suffix,fmt) FileIter
         filename=trim(OUTPUT_DIR)//'power/powerLP11a_y_'//trim(suffix)//'.dat'
         open(UNIT=9,FILE=filename,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
         do i = 1,NumPts
            write(UNIT=9, FMT="(es12.5)") power_LP11a_y(i)
         enddo
         close(UNIT=9)
         !WRITE TO FILE
         write(*,*) ' get_power: printing LP11b (y) mode power (signal) to file..'
         fmt = '(I5.5)'
         write (suffix,fmt) FileIter
         filename=trim(OUTPUT_DIR)//'power/powerLP11b_y_'//trim(suffix)//'.dat'
         open(UNIT=9,FILE=filename,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
         do i = 1,NumPts
            write(UNIT=9, FMT="(es12.5)") power_LP11b_y(i)
         enddo
         close(UNIT=9)
         !WRITE TO FILE
         write(*,*) ' get_power: printing LP21a (y) mode power (signal) to file..'
         fmt = '(I5.5)'
         write (suffix,fmt) FileIter
         filename=trim(OUTPUT_DIR)//'power/powerLP21a_y_'//trim(suffix)//'.dat'
         open(UNIT=9,FILE=filename,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
         do i = 1,NumPts
            write(UNIT=9, FMT="(es12.5)") power_LP21a_y(i)
         enddo
         close(UNIT=9)
         !WRITE TO FILE
         write(*,*) ' get_power: printing LP21b (y) mode power (signal) to file..'
         fmt = '(I5.5)'
         write (suffix,fmt) FileIter
         filename=trim(OUTPUT_DIR)//'power/powerLP21b_y_'//trim(suffix)//'.dat'
         open(UNIT=9,FILE=filename,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
         do i = 1,NumPts
            write(UNIT=9, FMT="(es12.5)") power_LP21b_y(i)
         enddo
         close(UNIT=9)
         !WRITE TO FILE
         write(*,*) ' get_power: printing LP02 (y) mode power (signal) to file..'
         fmt = '(I5.5)'
         write (suffix,fmt) FileIter
         filename=trim(OUTPUT_DIR)//'power/powerLP02_y_'//trim(suffix)//'.dat'
         open(UNIT=9,FILE=filename,FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
         do i = 1,NumPts
            write(UNIT=9, FMT="(es12.5)") power_LP02_y(i)
         enddo
         close(UNIT=9)
      endif
   endif
!
   90 continue
   deallocate(zValues,sign_power,pump_power,diff_power,core_power,clad_power)
!
   if (modeProj) then
   deallocate(power_LP01_x,norm_LP01_x,coef_LP01_r_x,coef_LP01_c_x)
   deallocate(power_LP11a_x,norm_LP11a_x,coef_LP11a_r_x,coef_LP11a_c_x)
   deallocate(power_LP11b_x,norm_LP11b_x,coef_LP11b_r_x,coef_LP11b_c_x)
   deallocate(power_LP21a_x,norm_LP21a_x,coef_LP21a_r_x,coef_LP21a_c_x)
   deallocate(power_LP21b_x,norm_LP21b_x,coef_LP21b_r_x,coef_LP21b_c_x)
   deallocate(power_LP02_x,norm_LP02_x,coef_LP02_r_x,coef_LP02_c_x)
!
   deallocate(power_LP01_y,norm_LP01_y,coef_LP01_r_y,coef_LP01_c_y)
   deallocate(power_LP11a_y,norm_LP11a_y,coef_LP11a_r_y,coef_LP11a_c_y)
   deallocate(power_LP11b_y,norm_LP11b_y,coef_LP11b_r_y,coef_LP11b_c_y)
   deallocate(power_LP21a_y,norm_LP21a_y,coef_LP21a_r_y,coef_LP21a_c_y)
   deallocate(power_LP21b_y,norm_LP21b_y,coef_LP21b_r_y,coef_LP21b_c_y)
   deallocate(power_LP02_y,norm_LP02_y,coef_LP02_r_y,coef_LP02_c_y)
   endif
!
   99 continue
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
!   latest revision    - Mar 2022
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
!                      - DiffPower   : Diff exact to computed power ! alt: Norm
!                                      (available if NEXAXT=1)
!                      - CorePower   : (available if GEOM_NO=5)     ! alt: Coef_r
!                      - CladPower   : (available if GEOM_NO=5)     ! alt: Coef_c
!
!----------------------------------------------------------------------
!
subroutine compute_power(ZValues,Num_zpts,Fld, Power,DiffPower,CorePower,CladPower)
!
   use commonParam
   use laserParam
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
   real(8), intent(in)  :: ZValues(Num_zpts)
   integer, intent(in)  :: Fld
   real(8), intent(out) :: Power(Num_zpts)
   real(8), intent(out) :: DiffPower(Num_zpts)
   real(8), intent(out) :: CorePower(Num_zpts)
   real(8), intent(out) :: CladPower(Num_zpts)
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
   real*8 :: xnod (3,8)
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
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!
!..return pre-computed pump values if using constant pump or pump ODE model
   if (Fld.eq.0 .and. PLANE_PUMP.eq.1) then
      Power(1:Num_zpts) = PLANE_PUMP_POWER
      goto 90
   elseif (Fld.eq.0 .and. PLANE_PUMP.eq.2) then
      if (size(PUMP_VAL) .ne. Num_zpts) then
         write(*,*) 'compute_power: size(PUMP_VAL) .ne. Num_zpts. skipping.'
         goto 90
      endif
      Power(1:Num_zpts) = PUMP_VAL(1:Num_zpts)
      goto 90
   endif
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
!$OMP         facePower,faceDiffPower,modeNorm,modeCoef) &
!$OMP REDUCTION(+:Power,DiffPower,corePower,cladPower)   &
!$OMP SCHEDULE(DYNAMIC)
   do iel=1,NRELES_SUBD
      mdle = ELEM_SUBD(iel)
      if (GEOM_NO .eq. 5) call find_domain(mdle, ndom)
      call nodcor_vert(mdle, xnod)
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
      do i=1,Num_zpts
         if((ZValues(i).le.maxz).and.(ZValues(i).gt.minz)) then
            if (Fld .le. 9) then
               call compute_facePower(mdle,faceNum,Fld, facePower,faceDiffPower)
               DiffPower(i) = DiffPower(i) + abs(faceDiffPower)
               if (GEOM_NO .eq. 5) then
                  select case(ndom)
                  case(1,2); CorePower(i) = CorePower(i) + abs(facePower)
                  case(3,4); CladPower(i) = CladPower(i) + abs(facePower)
                  end select
               endif
            elseif ((Fld .ge. 13 .and. Fld .le. 16) .or. Fld .eq. 140 .or. Fld .eq. 150) then
               call compute_mode_power(mdle,faceNum,Fld, facePower,modeNorm,modeCoef)
               DiffPower(i) = DiffPower(i) + modeNorm       ! false name (calc norm)
               CorePower(i) = CorePower(i) + real(modeCoef) ! false name (calc coef_r)
               CladPower(i) = CladPower(i) + imag(modeCoef) ! false name (calc coef_r)
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
end subroutine compute_power
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
   character(len=4) :: etype,ftype
!
!..variables for geometry
   real(8), dimension(3)   :: xi,x,rn,x_new
   real(8), dimension(3,2) :: dxidt,dxdt,rt
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
   integer :: nint,icase,iattr,l,i,j
   real(8) :: weight,wa
   integer :: iel,nsign
   integer :: nflag,iload
!
!---------------------------------------------------------------------------------------
!
   FacePower = 0.d0
   FaceDiffPower = 0.0d0
   nflag = 1
!..element type
   etype = NODES(Mdle)%type
   nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
   call find_order(Mdle, norder)
   call find_orient(Mdle, nedge_orient,nface_orient)
   call nodcor(mdle, xnod)
   call solelm(mdle, zdofH,zdofE,zdofV,zdofQ)
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
         call zz_cross_product(zsolE(1:3,1),conjg(zsolE(1:3,2)), EtimesH1)
         FdotN = EtimesH1(1)*rn(1)+EtimesH1(2)*rn(2)+EtimesH1(3)*rn(3)
         FacePower = FacePower + (real(FdotN))*weight
!     ...if we have an exact
         if(NEXACT.eq.1) then
            call zz_cross_product(valE(1:3,1),conjg(valE(1:3,2)), EtimesH2)
            FaceDiffPower = FaceDiffPower   &
                           + abs(((EtimesH1(1)*rn(1)+EtimesH1(2)*rn(2)+EtimesH1(3)*rn(3))*weight) - &
                                 ((EtimesH2(1)*rn(1)+EtimesH2(2)*rn(2)+EtimesH2(3)*rn(3))*weight))
         endif
!  ...next check for pump, i.e, if Fld = 0
      else if(Fld.eq.0) then
         call zz_cross_product(zsolE(1:3,3),conjg(zsolE(1:3,4)), EtimesH1)
         FdotN = EtimesH1(1)*rn(1)+EtimesH1(2)*rn(2)+EtimesH1(3)*rn(3)
         FacePower = FacePower + (real(FdotN))*weight
!     ...if we have an exact
         if(NEXACT.eq.1) then
            call zz_cross_product(valE(1:3,3),conjg(valE(1:3,4)), EtimesH2)
            FaceDiffPower = FaceDiffPower   &
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
!   purpose            - Compute projection of field onto LP modes
!
!   arguments
!        in:
!                      - Mdle       : middle element node
!                      - Facenumber : element face used for integration
!                      - Fld        : 13  LP01  Mode
!                                     14  LP11a Mode
!                                     140 LP11b Mode
!                                     15  LP21a Mode
!                                     150 LP21b Mode
!                                     16  LP02  Mode
!       out:
!                      - ModeNorm   : norm of the mode (for normalization)
!                      - ModeCoef   : coefficient in the projection on mode
!
!----------------------------------------------------------------------
subroutine compute_mode_power(Mdle,Facenumber,Fld, ModePower,ModeNorm,ModeCoef)
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
   integer   , intent(in)  :: Mdle
   integer   , intent(in)  :: Fld
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
   character(len=4) :: etype,ftype
!
!..variables for geometry
   real(8), dimension(3)   :: xi,x,rn,x_new
   real(8), dimension(3,2) :: dxidt,dxdt,rt
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
   integer :: nint,icase,iattr,l,i,j
   real(8) :: weight,wa
   integer :: iel,nsign
   integer :: nflag,iload
!
!---------------------------------------------------------------------------------------
!
   ModePower = 0.d0
   nflag = 1
!..element type
   etype = NODES(Mdle)%type
   nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
   call find_order(Mdle, norder)
   call find_orient(Mdle, nedge_orient,nface_orient)
   call nodcor(mdle, xnod)
   call solelm(mdle, zdofH,zdofE,zdofV,zdofQ)
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
!  ...compute field of the LP mode
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
!  ...compute field of the LP mode
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



