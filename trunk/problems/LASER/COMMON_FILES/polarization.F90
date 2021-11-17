!
#include "typedefs.h"
!
!-------------------------------------------------------------------------------
! Routine: get_bgPol
!
! last modified: June 2021
!
! purpose: returns the background polarization
!
! input:
!           Dom_flag - 1 for core, 0 for cladding
!           Fld_flag - 1 for signal, 0 for pump
!           Delta_n  - thermally induced refractive index perturbation
!           X        - coordinate
! output:
!           Bg_pol - value of background polarization
!
!-------------------------------------------------------------------------------
subroutine get_bgPol(Dom_flag,Fld_flag,Delta_n,X, Bg_pol)
!
   use commonParam
   use laserParam
!
   implicit none
!
   integer, intent(in)  :: Dom_flag,Fld_flag
   real(8), intent(in)  :: Delta_n
   real(8), intent(in)  :: X(3)
   VTYPE  , intent(out) :: Bg_pol(3,3)
!
   real(8) :: aux(3,3)
!
   real(8), parameter :: spatialFreq = 1.0d0
   real(8), parameter :: temporalFreq = 1.0d0
!..perturbation amplitude
   !real(8), parameter :: gratingCore = 1.25d-4
   real(8), parameter :: gratingCore = 2.50d-4
   !real(8), parameter :: gratingCore = 5.00d-4
   !real(8), parameter :: gratingCore = 1.00d-3
!
   real(8), parameter :: gratingClad = 0.0d0
!
!..Grating frequency depends on beatLength
   real(8), parameter :: beatLength =  0.0511d0 ! LP01 (x), LP02 (x); nx = 1.45
!   real(8), parameter :: beatLength =  0.0203d0 ! LP01 (x), LP11 (x); nx = 1.45
!   real(8), parameter :: beatLength =  0.0453d0 ! LP01 (x), LP21 (x); nx = 1.45
!
   !real(8), parameter :: beatLength =  0.0719d0 ! LP01 (x), LP02 (x); nx = 1.15
   !real(8), parameter :: beatLength = 11.8005d0 ! LP01 (x), LP01 (y); nx = 1.45, ny = 1.65
!
   real(8), parameter :: gratingFreq =  beatLength * spatialFreq
!
!..define asymmetric perturbation region
   real(8) :: x_perturb, y_perturb ! coordinates of the center of perturbation region
   real(8) :: r_perturb ! radius of the perturbation region (relative to core size)
   real(8) :: r ! radial distance of X from perturbation center
!
!..define a phase shift for the sin(z) perturbation function
   !real(8) :: phaseShift = 0.d0
   !real(8) :: phaseShift = -PI / 4.d0
   !real(8) :: phaseShift = -PI / 2.d0
   !real(8) :: phaseShift = -PI * 0.75d0
   !real(8) :: phaseShift = -PI
   real(8) :: phaseShift = -PI * 0.67d0
!
!-------------------------------------------------------------------------------
!
   aux = 0.d0
   if (Dom_flag.eq.1) then
      aux = CORE_N+Delta_n*IDENTITY
      if (ART_GRATING .eq. 1) then
      !  symmetric grating
      !..EXP 0 (LP01 to LP02)
         x_perturb = 0.0d0; y_perturb = 0.0d0
         r_perturb = 0.5d0*R_CORE
      !  asymmetric grating
      !..EXP 1 (LP01 to LP02)
!         x_perturb = 0.3d0*R_CORE; y_perturb = 0.0d0
!         r_perturb = 0.7d0*R_CORE
      !..EXP 2 (LP01 to LP11a)
!         x_perturb = -0.4d0*R_CORE; y_perturb = 0.0d0
!         r_perturb =  0.6d0*R_CORE
      !..EXP 3 (LP01 to LP21a (and LP11b))
!         x_perturb = 0.0d0; y_perturb = 0.5d0*R_CORE
!         r_perturb = 0.5d0*R_CORE
      !..EXP 4 (LP01 to LP21b)
!         x_perturb = 0.35d0*R_CORE; y_perturb = -0.35d0*R_CORE
!         r_perturb = 0.5d0*R_CORE
!
         r = sqrt((X(1)-x_perturb)**2.d0+(X(2)-y_perturb)**2.d0)
!     ...circular perturbation region
         if (r .le. r_perturb) then
            aux = aux + gratingCore*sin(gratingFreq*X(3) + phaseShift)*IDENTITY
!     ...annulus perturbation region
         else
            !aux = aux + gratingCore*sin(gratingFreq*X(3) + phaseShift)*IDENTITY
         endif
      endif
      !if (ART_GRATING .eq. 1) aux(1,2) = aux(1,2) + gratingCore*sin(gratingFreq*X(3))
      !if (ART_GRATING .eq. 1) aux(2,1) = aux(2,1) + gratingCore*sin(gratingFreq*X(3))
      !if (ART_GRATING .eq. 1) aux(1,1) = aux(1,1) + gratingCore*sin(gratingFreq*X(3)) ! grating in x-pol. only
!     TODO: double check if eps*eps is correct here, or if it is eps^T * eps, etc.
      !call DGEMM('N', 'N', 3, 3, 3, 1.d0, aux, 3, aux, 3, 0.d0, aux, 3)
      aux(1,1) = aux(1,1)*aux(1,1)
      aux(2,2) = aux(2,2)*aux(2,2)
      aux(3,3) = aux(3,3)*aux(3,3)
      Bg_pol = aux-IDENTITY
   elseif (Dom_flag.eq.0) then
      aux = CLAD_N+Delta_n*IDENTITY
      !if (ART_GRATING .eq. 1) aux = aux + gratingClad*sin(gratingFreq*X(3))*IDENTITY
      !if (ART_GRATING .eq. 1) aux(1,2) = aux(1,2) + gratingClad*sin(gratingFreq*X(3))
      !if (ART_GRATING .eq. 1) aux(2,1) = aux(2,1) + gratingClad*sin(gratingFreq*X(3))
      !if (ART_GRATING .eq. 1) aux(1,1) = aux(1,1) + gratingClad*sin(gratingFreq*X(3)) ! grating in x-pol. only
      !call DGEMM('N', 'N', 3, 3, 3, 1.d0, aux, 3, aux, 3, 0.d0, aux, 3)
      aux(1,1) = aux(1,1)*aux(1,1)
      aux(2,2) = aux(2,2)*aux(2,2)
      aux(3,3) = aux(3,3)*aux(3,3)
      Bg_pol = aux-IDENTITY
   else
      write(*,*) ' get_bgPol: Dom_flag must be 0 or 1. stop.'
      stop
   endif
!
   select case(Fld_flag)
      case(1)
         Bg_pol = ZI*OMEGA*OMEGA_RATIO_SIGNAL*Bg_pol
      case(0)
         Bg_pol = ZI*OMEGA*OMEGA_RATIO_PUMP*Bg_pol
      case default
        write(*,*) ' get_bgPol: Fld_flag must be 0 or 1. stop.'
        stop
   endselect

   if(real(Bg_pol(1,1)).ne.0.d0 .or. &
      real(Bg_pol(2,2)).ne.0.d0 .or. &
      real(Bg_pol(3,3)).ne.0.d0) then
      write(*,*) ' get_bgPol: Bg_pol must be purely imaginary. stop.'
      stop
   endif

end subroutine get_bgPol
!
!
!-------------------------------------------------------------------------------
!
! Routine: get_activePol
!
! last modified: June 2021
!
! purpose:  returns the active gain polarization
!
! input:
!           ZsolQ    - all 12 EH fields (6 EH for signal, 6 EH for pump)
!           Dom_flag - 1 for core, 0 for cladding
!           Fld_flag - 1 for signal, 0 for pump
!           Delta_n  - thermally induced refractive index perturbation
!
! output:
!           active_pol - value of active polarization
!
!-------------------------------------------------------------------------------
subroutine get_activePol(ZsolQ,Fld_flag,Delta_n, Active_pol)
!
   use commonParam
   use laserParam
!
   implicit none
!
   VTYPE  , intent(in)  :: ZsolQ(12)
   integer, intent(in)  :: Fld_flag
   real(8), intent(in)  :: Delta_n
   VTYPE  , intent(out) :: Active_pol(3,3)
!
   VTYPE, dimension(3) :: Es,Hs,Ep,Hp,ETimesHs,ETimesHp
!
   real(8) :: eta,Nex,Ngd,sum1,sum2,Is,Ip,g0,gain_ampl
   VTYPE   :: gain
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
   Is = sqrt(real(EtimesHs(1))**2+real(EtimesHs(2))**2+real(EtimesHs(3))**2)
!
!..compute pump irradiance
   Ip = 0.d0
   if (PLANE_PUMP .eq. 1) then
!     set fake pump power the same here and in thermal polarization computation
!  ...assume pump is a plane wave in fiber cladding (cladding-pumped)
      Ip = PLANE_PUMP_POWER / (PI*R_CLAD*R_CLAD) ! calculate non-dimensional irradiance
   else
      call zz_cross_product(Ep,conjg(Hp), ETimesHp)
      Ip = sqrt(real(EtimesHp(1))**2+real(EtimesHp(2))**2+real(EtimesHp(3))**2)
   endif
!
   if (Is .eq. 0.d0 .or. Ip .eq. 0.d0) then
      active_pol = ZERO
      goto 60
   endif
!
   sum1 = (SIGMA_S_ABS/OMEGA_SIGNAL)*Is+(SIGMA_P_ABS/OMEGA_PUMP)*Ip
   sum2 = ((SIGMA_S_ABS+SIGMA_S_EMS)/OMEGA_SIGNAL)*Is + &
          ((SIGMA_P_ABS+SIGMA_P_EMS)/OMEGA_PUMP)*Ip
!
   eta = sum1/(TAU_0+sum2)
!
   if(eta.gt.1.d0) then
      write(*,*) 'error from get_activePol: eta must be <= 1. stop.'
      stop
   endif
!
!..signal gain
   if(Fld_flag.eq.1) then
      gain = -SIGMA_S_ABS + (SIGMA_S_ABS+SIGMA_S_EMS)*eta
!
!..pump gain
   elseif(Fld_flag.eq.0) then
      gain = -SIGMA_P_ABS + (SIGMA_P_ABS+SIGMA_P_EMS)*eta
   else
      write(*,*) 'error from get_activePol: fld_flag must be 1 or 0. stop.'
      stop
   endif
!
!..non-dimensional scaling factor for gain function
   g0 = ACTIVE_GAIN*L_0*SIGMA_0*NU_0
   gain = g0 * gain * N_TOTAL
!
 1500 format(I2,A,F10.6)
!
!..compute active polarization term from gain function
   active_pol = -(CORE_N+Delta_n*IDENTITY)*gain
!
  60 continue
!
end subroutine get_activePol
!
!
!-------------------------------------------------------------------------------
!  subroutine: get_ramanPol
!
!  last modified: Mar 2019
!
!  purpose:    returns the Raman polarization
!
!  input:      E,H        : electric and magnetic fields (L2 variables)
!              domain flag: core - 1, cladding - 0
!              field  flag: signal - 1, pump - 0
!              Delta_n    : thermally-induced refractive index perturbation
!
!  output:
!              raman_pol - value of Raman polarization
!-------------------------------------------------------------------------------
subroutine get_ramanPol(E,H,Dom_flag,Fld_flag,Delta_n, Raman_pol)
!
   use commonParam
   use laserParam
!
   implicit none
!
   VTYPE  , intent(in)  :: E(3), H(3)
   integer, intent(in)  :: Dom_flag, Fld_flag
   real(8), intent(in)  :: Delta_n
   VTYPE  , intent(out) :: Raman_pol(3,3)
!..irradiance: |Re{ExH^*}|
   VTYPE :: EtimesH(3)
   VTYPE :: I
   integer :: j,k
!
!-------------------------------------------------------------------------------
!
!..compute irradiance
   call zz_cross_product(E,conjg(H), EtimesH)
   I = sqrt(real(EtimesH(1))**2+real(EtimesH(2))**2+real(EtimesH(3))**2)
!
   select case(Fld_flag)
      case(1)
         if (Dom_flag.eq.1) then
            Raman_pol = -(CORE_N+Delta_n*IDENTITY)*UPSILON_RAMAN_SIGNAL*I*RAMAN_GAIN
         elseif (dom_flag.eq.0) then
            Raman_pol = -(CLAD_N+Delta_n*IDENTITY)*UPSILON_RAMAN_SIGNAL*I*RAMAN_GAIN
         else
            write(*,*) ' get_ramanPol: dom_flag must be 0 or 1. stop.'
            stop
         endif
         do j = 1,3; do k = 1,3
            if(real(Raman_pol(k,j)).gt.0.d0) then
               write(*,*) ' get_ramanPol: for signal, Raman polarization', &
                          ' must be purely real with negative real part. stop.'
               stop
            endif
         enddo; enddo
      case(0)
         if (Dom_flag.eq.1) then
            Raman_pol = -(CORE_N+Delta_n*IDENTITY)*UPSILON_RAMAN_PUMP*I*RAMAN_GAIN
         elseif (Dom_flag.eq.0) then
            Raman_pol = -(CLAD_N+Delta_n*IDENTITY)*UPSILON_RAMAN_PUMP*I*RAMAN_GAIN
         else
            write(*,*) ' get_ramanPol: dom_flag must be 0 or 1. stop.'
            stop
         endif
         do j = 1,3; do k = 1,3
            if(real(Raman_pol(k,j)).lt.0.d0) then
               write(*,*) ' get_ramanPol: for pump, Raman polarization', &
                          ' must be purely real with positive real part. stop.'
               stop
            endif
         enddo; enddo
      case default
         write(*,*) ' get_ramanPol: fld_flag must be 0 or 1. stop.'
         stop
   end select
   do j = 1,3
      do k = 1,3
         if (aimag(Raman_pol(k,j)).ne.0.d0) then
            write(*,*) ' get_ramanPol: aimag(Raman_pol) not equal 0. stop.'
            stop
         endif
      enddo
   enddo
!
end subroutine get_ramanPol
