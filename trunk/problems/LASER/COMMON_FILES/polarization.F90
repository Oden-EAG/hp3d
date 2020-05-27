!
#include "implicit_none.h"
!
!-------------------------
! Routine: get_bgPol
!
! last modified: Apr 2019
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
!-------------------------
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
   real(8), parameter :: gratingCore = 1.0d-3
   real(8), parameter :: gratingClad = 1.0d-6
   real(8), parameter :: gratingFreq = 0.0511d0
!
   aux = 0.d0
   if (Dom_flag.eq.1) then
      aux = CORE_N+Delta_n*IDENTITY
      if (ART_GRATING .eq. 1) aux = aux + gratingCore*sin(gratingFreq*X(3))*IDENTITY
      aux(1,1) = aux(1,1)*aux(1,1)
      aux(2,2) = aux(2,2)*aux(2,2)
      aux(3,3) = aux(3,3)*aux(3,3)
      Bg_pol = aux-IDENTITY
   elseif (Dom_flag.eq.0) then
      aux = CLAD_N+Delta_n*IDENTITY
      if (ART_GRATING .eq. 1) aux = aux + gratingClad*sin(gratingFreq*X(3))*IDENTITY
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
!---------------------------------------------------------------------------
!
! Routine: get_activePol
!
! last modified: Apr 2019
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
!---------------------------------------------------------------------------
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
!..modified irradiance experiment (birefringent fiber)
   VTYPE, dimension(3) :: Es_mod,Hs_mod
   integer :: modified
!
   real(8) :: eta,Nex,Ngd,sum1,sum2,Is,Ip,g0,gain_ampl
   VTYPE   :: gain
!
!---------------------------------------------------------------------------
!
!..get fields
   Es = ZsolQ(1:3)
   Hs = ZsolQ(4:6)
   Ep = ZsolQ(7:9)
   Hp = ZsolQ(10:12)
!
!..modified irradiance experiment
!   Es_mod(1) = Es(1)+Es(2)
!   Es_mod(2) = ZERO
!   Es_mod(3) = Es(3)
!!
!   Hs_mod(1) = ZERO
!   Hs_mod(2) = Hs(1)+Hs(2)
!   Hs_mod(3) = Hs(3)
!
!..compute irradiance
! if modified=1, then gain model is not working properly (yields wrong efficiency)
   modified = 0
   if (modified .eq. 1) then
      !call zz_cross_product(Es_mod,conjg(Hs_mod), ETimesHs)
      write(*,*) 'do not use modified=1 in polarization! stop.'
      stop
   else
      call zz_cross_product(Es,conjg(Hs), ETimesHs)
   endif
   call zz_cross_product(Ep,conjg(Hp), ETimesHp)
   Is = sqrt((real(EtimesHs(1))**2+real(EtimesHs(2))**2+real(EtimesHs(3))**2))
   Ip = sqrt((real(EtimesHp(1))**2+real(EtimesHp(2))**2+real(EtimesHp(3))**2))
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
!..Non-dimensional scaling factor for gain function
   g0 = ACTIVE_GAIN*L_0*SIGMA_0*NU_0
   gain = g0 * gain * N_TOTAL
!
 1500 format(I2,A,F10.6)
!
!..Compute active polarization term from gain function
   active_pol = -(CORE_N+Delta_n*IDENTITY)*gain
!
  60 continue
!
end subroutine get_activePol
!
!
!--------------------------------------------------
!  subroutine: get_ramanPol
!
!  last modified: Mar 2019
!
!  purpose:    returns the Raman polarization
!
!  input:      E,H        : electric and magnetic fields (L2 variables)
!              domain flag: core - 1, cladding - 0
!              field  flag: signal - 1, pump - 0
!              Delta_n    : thermally induced refractive index perturbation
!
!  output:
!              raman_pol - value of Raman polarization
!--------------------------------------------------
subroutine get_ramanPol(E,H,Dom_flag,Fld_flag,Delta_n, Raman_pol)
!
   use commonParam
   use laserParam
!
   implicit none
!--------------------------------------------------
   VTYPE  , intent(in)  :: E(3), H(3)
   integer, intent(in)  :: Dom_flag, Fld_flag
   real(8), intent(in)  :: Delta_n
   VTYPE  , intent(out) :: Raman_pol(3,3)
!..irradiance: |Re{ExH^*}|
   VTYPE :: EtimesH(3)
   VTYPE :: I
   integer :: j,k
!--------------------------------------------------
!
!..Compute irradiance
   call zz_cross_product(E,conjg(H), EtimesH)
   I = sqrt((real(EtimesH(1))**2+real(EtimesH(2))**2+real(EtimesH(3))**2))
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
               write(*,*) ' get_ramanPol: for signal, Raman polarization must be purely real with negative real part. stop.'
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
               write(*,*) ' get_ramanPol: for pump, Raman polarization must be purely real with positive real part. stop.'
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
