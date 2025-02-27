!
#include "typedefs.h"
!
      module bessel_evaluation
      implicit none
!     Define support points of AAA rational approximant as parameters
! !     Number of support points
!       integer, parameter :: NSP= 11 
! !
! !     Coordinates of support points       ~~~~~~~~~~ADD PREFIX OF DOUBLE PRECISION NUMBBERS!!!!
!       real(8), parameter, dimension(NSP) :: ZSP = (/ -12.69999999999709, 2.7585880076148896, 12.69999999999709, -5.034233036589285, 8.242165636875143, -10.630053141212557, -2.1866669868904864, 10.79051887276728, -6.991238972761494, 11.966732838453026, 6.351785532096983 /)
! !
! !     Weights of support points
!       real(8), parameter, dimension(NSP) :: WSP = (/ -0.0056571080833777996, 0.2043089474418977, -0.1352973669184051, 0.24975116015380383, 0.5071715035472382, 0.025722484848089443, -0.19892940036765372, -0.4580279450360627, -0.14416531064415672, 0.39499840940114656, -0.4398753987276523 /)
! !
! !     Images of support points
!       real(8), parameter, dimension(NSP) :: FSP = (/ 0.0, 1.6378108732703571, 2.4868093695396e-12, 0.9869115657366958, 1.113411106543641, 0.2577944628424862, 1.3309707657234617, 0.5239226938660501, 0.7296267686829699, 0.20507305782356977, 1.417443896075109 /)
! !
! !     Product of weight times image for all support points
!       real(8), parameter, dimension(NSP) :: WFSP =(/ -0.0, 0.33461941562676195, -3.364587597267269e-13, 0.2464823085119468, 0.5646903849719326, 0.006631114164387208, -0.26476921633224515, -0.23997123482922508, -0.10518686976147264, 0.08100353165133942, -0.6234986989601156 /)
! !
!     Number of support points
      integer, parameter :: NSP= 10 
!
!     Coordinates of support points       ~~~~~~~~~~ADD PREFIX OF DOUBLE PRECISION NUMBBERS!!!!
      real(8), parameter, dimension(NSP) :: ZSP = (/ 0.5, 0.016625772516048798, -0.5, 0.1931110108031051, 0.3692679500356233, -0.3129842406718888, -0.3927325703743678, 0.4361634823741394, -0.11615909106694211, -0.46708058752803794 /)
!
!     Weights of support points
      real(8), parameter, dimension(NSP) :: WSP = (/ -0.02404793391416352, -0.23160888858301748, 0.2000138176094943, 0.13178127690978847, -0.12807858973391095, -0.44790648849497333, 0.6030318842024355, 0.10139573334133319, 0.2695480283255378, -0.47412883498288194 /)
!
!     Images of support points
      real(8), parameter, dimension(NSP) :: FSP = (/ -1.2836953722228372e-15, 1.0894006495434743, 0.0, 0.935558981529397, 0.46621754465639015, 0.5583036212407708, 0.33004379380205506, 0.23309930928360684, 0.9829825111800042, 0.10296474487346759 /)
!
!     Product of weight times image for all support points
      real(8), parameter, dimension(NSP) :: WFSP =(/ 3.087022147713233e-17, -0.2523148736623814, 0.0, 0.12328915721036514, -0.0597124856287971, -0.25006781450398125, 0.19902693084577336, 0.02363527540616955, 0.26496099776705606, -0.048818554531166854 /)
!
!
!
!     Mode number, Bessel order, Coefficients for the linear combination of the two l.i. solutions
!     mode 1 for k0 = 149.993333460866, r0 = 1300.d0, a=0.5d0
      integer, parameter    :: Imode         =            1
       complex(8), parameter :: ZLAMBDA_MODE =                  (38007107126.7520380677389983893284821d0,0.00000000000000000000000000000000000)
       complex(8), parameter :: ZCOEF10      =                 (1.57184336725612745203169806371674874d0,-0.00000000000000000000000000000000000)
       complex(8), parameter :: ZCOEF01      =                  (2588.85382397807097836830171213751471d0,0.00000000000000000000000000000000000)
! !     mode 2 for k0 = 149.993333460866, r0 = 1300.d0, a=0.5d0       ~~~~~~~~~~ADD PREFIX OF DOUBLE PRECISION NUMBBERS!!!!
      ! integer, parameter    :: Imode        =            2
      ! complex(8), parameter :: ZLAMBDA_MODE =                  (37954285745.8141508857519816549451173d0,0.00000000000000000000000000000000000)
      ! complex(8), parameter :: ZCOEF10      =                (0.380494806890728985773127546667222389d0,-0.00000000000000000000000000000000000)
      ! complex(8), parameter :: ZCOEF01      =                 (-8933.70820961992077300315636080164352d0,0.00000000000000000000000000000000000)
! !    mode 3 for k0 = 149.993333460866, r0 = 1300.d0, a=0.5d0       ~~~~~~~~~~ADD PREFIX OF DOUBLE PRECISION NUMBBERS!!!!
      ! integer, parameter    :: Imode        =            3
      ! complex(8), parameter :: ZLAMBDA_MODE =                  (37871108378.4249748433402696483647618d0,0.00000000000000000000000000000000000)
      ! complex(8), parameter :: ZCOEF10      =                (-1.02590650828264444059353639246478155d0,-0.00000000000000000000000000000000000)
      ! complex(8), parameter :: ZCOEF01      =                 (-2799.17121012013890047201708963238776d0,0.00000000000000000000000000000000000)
! !      mode 4 for k0 = 149.993333460866, r0 = 1300.d0, a=0.5d0
!       integer, parameter    :: Imode        =            4
!       complex(8), parameter :: ZLAMBDA_MODE =                  (37754502379.3769398924868581905513854d0,0.00000000000000000000000000000000000)
!       complex(8), parameter :: ZCOEF10      =               (-0.174962920969328175105637042258464969d0,-0.00000000000000000000000000000000000)
!       complex(8), parameter :: ZCOEF01      =                  (16577.5092176914679245858580875585293d0,0.00000000000000000000000000000000000)
!

! ! IMPEDANCE MODES WITH k0 = 149.993333460866, r0 = 1300.d0, a=0.5d0, d = -1/c = -8.47252801803306E-01
       ! integer, parameter    :: Imode        =            1
       ! complex(8), parameter :: ZLAMBDA_MODE =                 (38007109531.0212535832934389084965127d0,-558800.565253144981038272714515288780d0)
       ! complex(8), parameter :: ZCOEF10      =             (1.57192013800750286812492090876574160d0,-0.009295845543320480213508651202459053779d0)
       ! complex(8), parameter :: ZCOEF01      =                 (2589.47319818030213697897003400320251d0,-136.384737976323732024208621739425211d0)

       ! integer, parameter    :: Imode        =            2
       ! complex(8), parameter :: ZLAMBDA_MODE =                 (37954298296.9526978385617700680371866d0,-1263945.84638945095388623562760170783d0)
       ! complex(8), parameter :: ZCOEF10      =             (0.380678894620398692732367563928460725d0,-0.03696625027400351848227028422406469485d0)
       ! complex(8), parameter :: ZCOEF01      =                (-8936.66859328994639309431793099814257d0,-132.634800034923831092032632738758715d0)

       ! integer, parameter    :: Imode        =            3
       ! complex(8), parameter :: ZLAMBDA_MODE =                 (37871136582.3208965876604975171195350d0,-2584813.69792156068857807909590339654d0)
       ! complex(8), parameter :: ZCOEF10      =            (-1.02677231006622288890435081016054939d0,-0.008426457738186929433364246017511300460d0)
       ! complex(8), parameter :: ZCOEF01      =                 (-2808.53577462661293111495151099351989d0,515.598170754458927381134060433732170d0)

       ! integer, parameter    :: Imode        =            4
       ! complex(8), parameter :: ZLAMBDA_MODE =                 (37754552597.2426380153181895986410367d0,-4436477.61309035798134633873859358838d0)
       ! complex(8), parameter :: ZCOEF10      =             (-0.175541167534927505046985865469493100d0,0.05462104114912648328785153876984631624d0)
       ! complex(8), parameter :: ZCOEF01      =                  (16598.1534018554772046666048067143803d0,272.396494859041276659518202095058496d0)
!
!
      contains
!
!----------------------------------------------------------------------
         function Rfact(N) result(F)
         implicit none
         integer, intent(in) :: N
         real(8) :: F
         integer :: j
!
         if (N.lt.0) then
           write(*,*) 'Rfact: N = ',N
           stop 1
         endif
!
         F=1.d0
         if (N.eq.0) return
         do j=1,N
           F = F*2.d0/j
         enddo
!
         end function
!----------------------------------------------------------------------
!
!
         subroutine real_eval_aaa(Zev,Rev,DRev)
!           Zev: Location of point to evaluate (real number)
!           OUTPUT:
!           Rev :  Images of rational AAA aproximant at points Zev. Vector of Nev reals
!           DRev:  Derivative of rational AAA aproximant at points Zev. Vector of Nev reals
!
            implicit none
            real(8), intent(in)  :: Zev
            real(8), intent(out) :: Rev,DRev
!
!           workspace variables
            real(8) :: pf(NSP),dpf(NSP) ! partial fractions 1/(Zev - ZSP_j) and derivatives -1/(Zev - ZSP_j)**2
            real(8) :: scale,rn,rd,drn,drd ! scale,numerator, denominator and their derivatives
            integer :: j ! loop counters
            integer :: is_zsp ! if Zev coincides with the j-th ZSP, we store j; otherwise 0
! 
#if HP3D_DEBUG
!         ..Set iprint = 0/1 (Non-/VERBOSE)
            integer :: iprint
            iprint = 0
#endif
!           initialize the vector of partial fracions and is_zsp index
            pf = 0.d0
            is_zsp = 0
!           compute scale of support points
            scale = MAXVAL(ABS(ZSP)) 
!           loop over support points
            do j=1,NSP
               if ( ABS(Zev-ZSP(j))/scale.lt.1.d-12) then
                  ! we store j and leave pf(j)=0, so that it does not contribute to numerator or denominator
                  is_zsp = j
               else
                  pf(j) = 1.d0 / (Zev-ZSP(j))
                  dpf(j)=-1.d0 / (Zev-ZSP(j))**2
               endif
            enddo
!           Compute numerator                  rn = \sum_j {(w_j f_j) / (Zev - ZSP_j)}     = WFSP (dot) pf
            rn = DOT_PRODUCT(WFSP,pf)            
!           Compute derivative of numerator    drn= \sum_j {-(w_j f_j) / (Zev - ZSP_j)**2} = WFSP (dot) dpf
            drn= DOT_PRODUCT(WFSP,dpf)            
!           Compute denominator                rd = \sum_j {(w_j) / (Zev - ZSP_j)}         = WSP  (dot) pf
            rd = DOT_PRODUCT(WSP,pf)
!           Compute derivative of denominator  drd= \sum_j {-(w_j) / (Zev - ZSP_j)**2}     = WSP  (dot) dpf
            drd= DOT_PRODUCT(WSP,dpf)            
            if (is_zsp.gt.0) then
               j = is_zsp
!              evaluate using limit expressions for z -> ZSP_j
               Rev = FSP(j)
               DRev= ( rn - FSP(j)*rd )/WSP(j)
            else
               Rev = rn / rd
               DRev= ( drn*rd - rn*drd ) / rd**2
            endif
! 
#if HP3D_DEBUG
            if (iprint.eq.1) then
!           ...Print statements for verification
               write(*,*) 'real_eval_aaa: Zev, is_zsp =',Zev,is_zsp
               write(*,*) 'real_eval_aaa: rn , rd     =',rn,rd
               write(*,*) 'real_eval_aaa: Rev, DRev   =',Rev,DRev
               write(*,*) ' '
            endif
#endif
!
         end subroutine
!----------------------------------------------------------------------
!
!     EVALUATION BY FROBENIUS METHOD
!
!----------------------------------------------------------------------
!
!   subroutine name    - Bessel
!
!-----------------------------------------------------------------------
!
!   latest revision    - Jan 25
!
!   purpose            - evaluate Bessel function of an aribtrary
!                        complex order using Frobenius method
!                        with the expansion at arbitrary real
!                        argument
!   arguments
!     in:
!              Zlambda - complex order of the Bessel fucntion
!              Zc0,Zc1 - the first two complex coefficients in the
!                        Taylor expansion
!              Wavenum - real wavenumber (if c=1, then Wavenum = OMEGA)
!              R0      - point of expansion in r (x_0 = ln R_0)
!              R       - real argument
!     out:
!              Zbess   - complex value of the Bessel function
!              Zdbess  - its first derivative in R (not in x)
!              Zd2bess - its second derivative in R (not in x)
!
!-----------------------------------------------------------------------
!
      subroutine Bessel(Zlambda,Zc0,Zc1,Wavenum,R0,R, &
                        Zbess,Zdbess,Zd2bess)
      ! 
      use iso_fortran_env, wp => real64
      ! 
      implicit none
      integer :: Idec
      complex(wp), intent(in)  :: Zlambda,Zc0,Zc1
      real(8),    intent(in)  :: Wavenum,R0,R
      complex(wp), intent(out) :: Zbess,Zdbess,Zd2bess
      complex(wp) :: zc(0:1000), zsum,zloc,zdloc,zdloc_prev,zb(0:1000)
      integer :: n,j,iprint
      real(8) :: aux,x_0,x,dx
      real(wp) :: eps = 10._wp**(-15)
!
      if (R0.le. 0.d0) then
        write(*,*) 'Bessel: R0 = ',R0
        stop 1
      endif
      if (R.le. 0.d0) then
        write(*,*) 'Bessel: R = ',R
        stop 1
      endif
!
      iprint=0
!
      aux = (Wavenum*R0)**2
      x_0 = log(R0)
      x   = log(R)
      dx = x - x_0
!
      zc = complex(0._wp,0._wp)
      zc(0) = Zc0
      zc(1) = zc1
      zloc = zc(0) + zc(1)*dx
      Zdbess  = zc(1)
      Zd2bess = complex(0._wp,0._wp)
      zdloc_prev = 1._wp
      n=0
      do 
        zsum  = complex(0._wp,0._wp)
        do j=0,n
          zsum = zsum + Rfact(n-j)*zc(j)
        enddo
        zc(n+2) = (Zlambda*zc(n) - aux*zsum)/real( (n+1)*(n+2) ,8)
        zdloc = zc(n+2)*dx**(n+2)
        zloc = zloc + zdloc
        Zdbess  = Zdbess  + zc(n+2)*dx**(n+1) *real(n+2,8)
        Zd2bess = Zd2bess + zc(n+2)*dx**(n)   *real((n+1)*(n+2),8)
        if ((abs(zdloc).lt.eps).and.(abs(zdloc_prev).lt.eps)) exit
        n=n+1
        if (n+2.gt.1000) then
          write(*,*) 'Bessel 1: n = ',n
          stop 1
        endif
        zdloc_prev = zdloc
      enddo
      ! write(*,*) 'Bessel 1: n = ',n
   
      Zd2bess = Zd2bess/R**2 - Zdbess/R**2
      Zdbess  = Zdbess/R
      Zbess   = zloc
      
!
      if (iprint.eq.1) then
        write(*,7010) Zlambda,Zc0,Zc1
 7010   format('Bessel: Zlambda,Zc0,Zc1        = ',3(2e12.5,2x))
        write(*,7020)   Wavenum,R0,R
 7020   format('        Wavenum,R0,R             = ',3(e12.5,2x))
        write(*,7030) n,Zbess,Zdbess
 7030   format('        n,Zbess,Zdbess,Zd2bess = ',i3,3(2x,2e12.5))
!!!!        call pause
      endif
!
      end subroutine Bessel
!
!
!----------------------------------------------------------------------------
!
      subroutine bessel_preset( Wavenum, R0, Reval, Zval, Zdval, Zd2val )
      ! 
      use iso_fortran_env, wp => real64
      ! 
      implicit none
!
      real(8), intent(in)  :: Wavenum,R0,Reval
      complex(wp), intent(out) :: Zval,Zdval,Zd2val
      complex(wp) :: zbess10,zdbess10,zd2bess10,zbess01,zdbess01,zd2bess01, &
                     zone,zero
!
      zone = complex(1._wp,0._wp)
      zero = complex(0._wp,0._wp)
!
!  ...evaluate the solution at Reval
      call Bessel(ZLAMBDA_MODE,zone,zero,Wavenum,R0,Reval, zbess10,zdbess10,zd2bess10)
      call Bessel(ZLAMBDA_MODE,zero,zone,Wavenum,R0,Reval, zbess01,zdbess01,zd2bess01)
!
!  ...value
      Zval   = ZCOEF10*zbess10   + ZCOEF01*zbess01
      Zdval  = ZCOEF10*zdbess10  + ZCOEF01*zdbess01
      Zd2val = ZCOEF10*zd2bess10 + ZCOEF01*zd2bess01
      ! write(*,*) 'bessel_preset: zbess10,zbess01=',zbess10,zbess01
      ! write(*,*) 'bessel_preset: ZCOEF10,ZCOEF01=',ZCOEF10,ZCOEF01
      ! write(*,*) ''
!
      end subroutine bessel_preset
!
!----------------------------------------------------------------------------
!
      subroutine bessel_ivp( Wavenum, R0, Ra, Reval, Zlambda, Zi, Zdi, &
                             Zval, Zdval, Zd2val)
      ! 
      use iso_fortran_env, wp => real64
      ! 
      implicit none
!
      real(8), intent(in)  :: Wavenum,R0,Ra,Reval
      complex(wp), intent(in) :: Zlambda,Zi,Zdi
      complex(wp), intent(out) :: Zval,Zdval,Zd2val
      real(8) :: pi,r,a
      complex(wp) :: z1,z2,zbess10,zdbess10,zd2bess10,zbess01,zdbess01,zd2bess01, &
                     zdet,zone,zero,zv,zv1,zvoid
!
      zone = complex(1._wp,0._wp)
      zero = complex(0._wp,0._wp)
!     
      r = R0 - Ra
      call Bessel(zlambda,zone,zero,Wavenum,R0,r, zbess10,zdbess10,zvoid)
      call Bessel(zlambda,zero,zone,Wavenum,R0,r, zbess01,zdbess01,zvoid)
!
!  ...determine constants z1 and z2 for linear combination
      zdet = zbess10*zdbess01 - zdbess10*zbess01
      z1 = ( Zi*zdbess01 - Zdi*zbess01)/zdet
      z2 = (-Zi*zdbess10 + Zdi*zbess10)/zdet
!
!  ...evaluate the solution at Reval
      call Bessel(zlambda,zone,zero,Wavenum,R0,Reval, zbess10,zdbess10,zd2bess10)
      call Bessel(zlambda,zero,zone,Wavenum,R0,Reval, zbess01,zdbess01,zd2bess01)
!
!  ...value
      Zval   = z1*zbess10   + z2*zbess01
      Zdval  = z1*zdbess10  + z2*zdbess01
      Zd2val = z1*zd2bess10 + z2*zd2bess01
!
      write(*,7010) R0, Zval
 7010 format('bessel_ivp: R0, Zval = ',e12.5,2x,2e12.5)      
!
      end subroutine bessel_ivp
!
!
!-----------------------------------------------------------------------
      end module