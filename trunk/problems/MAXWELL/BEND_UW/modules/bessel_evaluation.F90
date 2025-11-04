!
#include "typedefs.h"
!
      module bessel_evaluation

      use iso_fortran_env, only: real128, quad_real=>real128 
      implicit none
      integer, parameter :: quad_complex=kind((1.0_quad_real,1.0_quad_real))
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
!       real(8), parameter, dimension(NSP) :: FSP = (/ 0.0, 1.6378108732703571, 2.4868093695396d-12, 0.9869115657366958, 1.113411106543641, 0.2577944628424862, 1.3309707657234617, 0.5239226938660501, 0.7296267686829699, 0.20507305782356977, 1.417443896075109 /)
! !
! !     Product of weight times image for all support points
!       real(8), parameter, dimension(NSP) :: WFSP =(/ -0.0, 0.33461941562676195, -3.364587597267269d-13, 0.2464823085119468, 0.5646903849719326, 0.006631114164387208, -0.26476921633224515, -0.23997123482922508, -0.10518686976147264, 0.08100353165133942, -0.6234986989601156 /)
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
!      ! integer, parameter    :: IMODE         =            1
      complex(8), parameter :: ZLAMBDA_MODE_DP =                  (38007107126.7520380677389983893284821d0,0.00000000000000000000000000000000000)
      complex(8), parameter :: ZCOEF10_DP      =                 (1.57184336725612745203169806371674874d0,-0.00000000000000000000000000000000000)
      complex(8), parameter :: ZCOEF01_DP      =                  (2588.85382397807097836830171213751471d0,0.00000000000000000000000000000000000)
! !     mode 2 for k0 = 149.993333460866, r0 = 1300.d0, a=0.5d0       ~~~~~~~~~~ADD PREFIX OF DOUBLE PRECISION NUMBBERS!!!!
!       ! integer, parameter    :: IMODE        =            2
      ! complex(8), parameter :: ZLAMBDA_MODE_DP =                  (37954285745.8141508857519816549451173d0,0.00000000000000000000000000000000000)
      ! complex(8), parameter :: ZCOEF10_DP      =                (0.380494806890728985773127546667222389d0,-0.00000000000000000000000000000000000)
      ! complex(8), parameter :: ZCOEF01_DP      =                 (-8933.70820961992077300315636080164352d0,0.00000000000000000000000000000000000)
! !    mode 3 for k0 = 149.993333460866, r0 = 1300.d0, a=0.5d0       ~~~~~~~~~~ADD PREFIX OF DOUBLE PRECISION NUMBBERS!!!!
!      ! integer, parameter    :: IMODE        =            3
      ! complex(8), parameter :: ZLAMBDA_MODE_DP =                  (37871108378.4249748433402696483647618d0,0.00000000000000000000000000000000000)
      ! complex(8), parameter :: ZCOEF10_DP      =                (-1.02590650828264444059353639246478155d0,-0.00000000000000000000000000000000000)
      ! complex(8), parameter :: ZCOEF01_DP      =                 (-2799.17121012013890047201708963238776d0,0.00000000000000000000000000000000000)
! !      mode 4 for k0 = 149.993333460866, r0 = 1300.d0, a=0.5d0
!      ! integer, parameter    :: IMODE        =            4
      ! complex(8), parameter :: ZLAMBDA_MODE_DP =                  (37754502379.3769398924868581905513854d0,0.00000000000000000000000000000000000)
      ! complex(8), parameter :: ZCOEF10_DP      =               (-0.174962920969328175105637042258464969d0,-0.00000000000000000000000000000000000)
      ! complex(8), parameter :: ZCOEF01_DP      =                  (16577.5092176914679245858580875585293d0,0.00000000000000000000000000000000000)
!

! ! IMPEDANCE MODES WITH k0 = 149.993333460866, r0 = 1300.d0, a=0.5d0, d = -1/c = -8.47252801803306E-01
!       ! integer, parameter    :: IMODE        =            1
       ! complex(8), parameter :: ZLAMBDA_MODE =                 (38007109531.0212535832934389084965127d0,-558800.565253144981038272714515288780d0)
       ! complex(8), parameter :: ZCOEF10      =             (1.57192013800750286812492090876574160d0,-0.009295845543320480213508651202459053779d0)
       ! complex(8), parameter :: ZCOEF01      =                 (2589.47319818030213697897003400320251d0,-136.384737976323732024208621739425211d0)

!       ! integer, parameter    :: IMODE        =            2
       ! complex(8), parameter :: ZLAMBDA_MODE =                 (37954298296.9526978385617700680371866d0,-1263945.84638945095388623562760170783d0)
       ! complex(8), parameter :: ZCOEF10      =             (0.380678894620398692732367563928460725d0,-0.03696625027400351848227028422406469485d0)
       ! complex(8), parameter :: ZCOEF01      =                (-8936.66859328994639309431793099814257d0,-132.634800034923831092032632738758715d0)

!       ! integer, parameter    :: IMODE        =            3
       ! complex(8), parameter :: ZLAMBDA_MODE =                 (37871136582.3208965876604975171195350d0,-2584813.69792156068857807909590339654d0)
       ! complex(8), parameter :: ZCOEF10      =            (-1.02677231006622288890435081016054939d0,-0.008426457738186929433364246017511300460d0)
       ! complex(8), parameter :: ZCOEF01      =                 (-2808.53577462661293111495151099351989d0,515.598170754458927381134060433732170d0)

!       ! integer, parameter    :: IMODE        =            4
       ! complex(8), parameter :: ZLAMBDA_MODE =                 (37754552597.2426380153181895986410367d0,-4436477.61309035798134633873859358838d0)
       ! complex(8), parameter :: ZCOEF10      =             (-0.175541167534927505046985865469493100d0,0.05462104114912648328785153876984631624d0)
       ! complex(8), parameter :: ZCOEF01      =                  (16598.1534018554772046666048067143803d0,272.396494859041276659518202095058496d0)
!
!
!   BENT THREE-LAYER WAVEGUIDE - MODES WITH NEUMANN ON LEFT BDRY, PML AT RIGHT BDRY
!   PARAMETERS: CPML=800, r0=2600
       ! integer, parameter    :: IMODE        =            300
       ! complex(8), parameter :: ZLAMBDA_MODE =                 (3.20269366898994099d11, -3.6352323135487189)
       ! complex(8), parameter :: ZCOEF10      =                 (1.d0,0.d0)
       ! complex(8), parameter :: ZCOEF01      =                 (1.9896129730909189d00, -8.9922785671882991d-08)

       ! integer, parameter    :: IMODE        =            301
       ! complex(8), parameter :: ZLAMBDA_MODE =                 (3.2011601122811722d11, -1.8019164642895208d04)
       ! complex(8), parameter :: ZCOEF10      =                 (-0.16119233143642369d08, 6.4792058794035203d-05)
       ! complex(8), parameter :: ZCOEF01      =                 (0.d0,1.d0)

       ! integer, parameter    :: IMODE        =            302
       ! complex(8), parameter :: ZLAMBDA_MODE =                 (3.1992651581813867d11, -1.0144922279661477d07)
       ! complex(8), parameter :: ZCOEF10      =                 (1.d0,0.d0)
       ! complex(8), parameter :: ZCOEF01      =                 (5.3893057338124903d00, -0.77616216126406644d00)

      ! !!! EVALUATING SOLUTION FOR: even1 -------- bctype = pml -------- r0 = 2600.0
      integer, parameter    :: IMODE        =            300
      complex(quad_complex), parameter :: ZLAMBDA_MODE = (320269366898.994098682030728927095595805868637047582019833380473920493940137629381185178851565885561478_quad_real , -3.63523231354871885693339769698359086091496517864875676271012325202737811452384188023125926024068439506_quad_real)
      complex(quad_complex), parameter :: ZCOEF10      = (1.0_quad_real, 0.0_quad_real)
      complex(quad_complex), parameter :: ZCOEF01      = (1.98961297309091893478652849940301575638240189047012634276987535540896806423028422985006625508207588194_quad_real , -0.0000000899227856718829911309410144836516247383549318625788366469756399938755268589275468289781904093032258958_quad_real)
      complex(quad_complex), parameter :: ZCOEF10_L1   = (0.119069346213287782721447060245441453770055299465827967712409380363447753441896185833673909086937850594_quad_real , 0.00000000331407730027735211080668529470613654136704330110324273610312877373937399258420072243162134253316015265_quad_real)
      complex(quad_complex), parameter :: ZCOEF01_L1   = (1.16149212189883275314379059561332790331444350257830823423202827944442725405652016516140923254957431076_quad_real , 0.0000000290754648274820615391144329511059225457047806290400777014833835973744959836563510172819732680516267853_quad_real)
      complex(quad_complex), parameter :: ZCOEF10_R1   = (0.631711833607984923500230106907207034687534349137211439441151363410289299517303631181001849148770521862_quad_real , -0.0000000826909685271221224525487264552184967499484669623578426894221722190565144245629868476325054651350536624_quad_real)
      complex(quad_complex), parameter :: ZCOEF01_R1   = (-4.65633040977606844030602489771690683278015882869808097595314906158520991300450809670925752732063499274_quad_real , -0.0000000710690213825845800515665876733761145293204391021911258484210166776551889529859062274643972396221016206_quad_real)
      complex(quad_complex), parameter :: ZCOEF10_L2   = (0.000000000000532708985251729840865814582750078673979787642248033029049079036637381066025003859319007651123723693082_quad_real , 0.0000000000000000000426577763757070467948535878587078182322963910949041529830715212677750271723930991141745524668187336234_quad_real)
      complex(quad_complex), parameter :: ZCOEF01_L2   = (0.00000000000708371839050983212587554748363261198287927303960296016307803133881409511876763539023157986377965850762_quad_real , 0.000000000000000000556490180473267813976180814979863550155315354428332403598164921868352660386939332228747225257342962847_quad_real)
      complex(quad_complex), parameter :: ZCOEF10_R2   = (-0.00029528625902482279241329248990027637406912973518034139478948953630434207098914721643016203495732018208_quad_real , -0.0000277665477894270474012233953718870043732376293152907743124647371969974637502104437626943266364696540631_quad_real)
      complex(quad_complex), parameter :: ZCOEF01_R2   = (-0.0000389015756560710820369318352141298450513113776079516982349551037384021194076333155108608541221696486993_quad_real , 0.00150158463579514486712129881798341133031923955629094576087249383545058107734083204902373961834667060547_quad_real)

      ! !!! EVALUATING SOLUTION FOR: odd1 -------- bctype = pml -------- r0 = 2600.0
      ! integer, parameter    :: IMODE        =            301
      ! complex(quad_complex), parameter :: ZLAMBDA_MODE = (320116011228.117221479303269603560379471184030499200846631082044828356745338124538189454119852109146113_quad_real , -18019.1646428952080527734051368663965048143446734561123525288211022906344636048342887718349254395827326_quad_real)
      ! complex(quad_complex), parameter :: ZCOEF10      = (-0.161192331436423687655854082998771310409209297150465006301180487793791276037727653325869561862229084149_quad_real , 0.0000647920587940352025514707908674528374065513381068740961011482923907444280997942534520136311599344097963_quad_real)
      ! complex(quad_complex), parameter :: ZCOEF01      = (1.0_quad_real, 0.0_quad_real)
      ! complex(quad_complex), parameter :: ZCOEF10_L1   = (-0.100324986113793889737611515695405018840458556991978051571200299241896403094078191559932089022252594045_quad_real , 0.00000928061164341713670535709648178324590438023983742859619436809873111505514432658664797749087110726227073_quad_real)
      ! complex(quad_complex), parameter :: ZCOEF01_L1   = (-0.855516207342229560746875675238535699437133792809575381437547602826225559626832823180024494502975922259_quad_real , 0.0000946109329610224922393532905120115250182794411520422600760157894688018641358379590142288240533695771599_quad_real)
      ! complex(quad_complex), parameter :: ZCOEF10_R1   = (0.167805645845308367788491527713816521618838172856320414020135447999108729413212685595712047857235598162_quad_real , -0.0000699568096309775357384767343555227485100971908408487219219224944663390529574432933024744910715841322952_quad_real)
      ! complex(quad_complex), parameter :: ZCOEF01_R1   = (-0.928727344347738654374960882603915217665612372920504277312490183692612766822046119902785721869354488336_quad_real , -0.000187267033225372300145550920810109819647038009075578085742767090801674967493407809157040936856475432335_quad_real)
      ! complex(quad_complex), parameter :: ZCOEF10_L2   = (-0.00000000000450950511138411948544186674625848433520791665224210258488964307368829707903897787703942583256915066866_quad_real , -0.000000000000000865620275798845591801963875065376038296093914714287887968739554968900403772426775077731390880249590576_quad_real)
      ! complex(quad_complex), parameter :: ZCOEF01_L2   = (-0.0000000000559953093628353507708210725641706265792501605174800671342333179006192946826896639904554541970517625243_quad_real , -0.0000000000000102657653277768394261306563096005519570935339330422643343422802219312760170190627641010323898243039868_quad_real)
      ! complex(quad_complex), parameter :: ZCOEF10_R2   = (0.00336145013612462728014735156085174337625242086099971663615686452082292383009053861550742894798588352671_quad_real , -0.00192526376713995316802227842893503625756179099536643517293029180866685072062818098330051823873033988199_quad_real)
      ! complex(quad_complex), parameter :: ZCOEF01_R2   = (-0.0139440782509516575073693959162516209503069101590772018537632436735803428107200174192882107531489800496_quad_real , -0.0228758937019231531041455959414524362705350985496450606309240971990251134656695446256904481443389695841_quad_real)

      ! !!! EVALUATING SOLUTION FOR: even2 -------- bctype = pml -------- r0 = 2600.0
      ! integer, parameter    :: IMODE        =            302
      ! complex(quad_complex), parameter :: ZLAMBDA_MODE = (319926515818.1386733528945469743382071912730506320476612317743459_quad_real , -10144922.279661476989626320316503859991794980759476618689407814845562996703948273620895736789508111965_quad_real)
      ! complex(quad_complex), parameter :: ZCOEF10      = (1.0_quad_real, 0.0_quad_real)
      ! complex(quad_complex), parameter :: ZCOEF01      = (5.389305733812490301546127711645665435051389372658982903074790872_quad_real , -0.776162161264066440622455608793011248631576830920323842468251763546646659823001132899116829802157362217_quad_real)
      ! complex(quad_complex), parameter :: ZCOEF10_L1   = (-0.9031237430153608247102020818191046098917800263357707617444466564_quad_real , 0.0344196007907533224548083002208564954862397354419063234980802244745824898273296115170496643909796943434_quad_real)
      ! complex(quad_complex), parameter :: ZCOEF01_L1   = (-6.069132327032090823852086667626621408752945001033328063565905569_quad_real , 0.329612664409640332775373647269450552928339334374059013693506966106853782223685379438337306358979510412_quad_real)
      ! complex(quad_complex), parameter :: ZCOEF10_R1   = (-1.150279431943048815718936210346622961168026850917599447777724802_quad_real , 0.0880546631921008575473791645731583879321189991480082279896957025449973945023436280068117475461903074655_quad_real)
      ! complex(quad_complex), parameter :: ZCOEF01_R1   = (1.96936763313706616266089571334800365690919430513510029332253496_quad_real , 1.09052446462813677163687309284771711462433154200357833304530650904893844320711070705205175660985991888_quad_real)
      ! complex(quad_complex), parameter :: ZCOEF10_L2   = (-0.0000000009982244922151252141160425076494952383509008036404199441598229820_quad_real , -0.000000000148622882459080130655946285495042339215757323475390263226792136268423989080530393626934833437756841532_quad_real)
      ! complex(quad_complex), parameter :: ZCOEF01_L2   = (-0.00000001122604640888514379006552828139082968208776177174470561483362359_quad_real , -0.00000000160355494002356368864530036953422270605605538565348691966027397109881033052793367639398673142654306957_quad_real)
      ! complex(quad_complex), parameter :: ZCOEF10_R2   = (-0.5981350881691439304360131797939081961590327881913672706669770238_quad_real , 0.0257791607070794002347862155183844157460710336631684662511247372767369692170185414631999500562160719801_quad_real)
      ! complex(quad_complex), parameter :: ZCOEF01_R2   = (0.2444727496912840077736544700427194903362497106851841531774371913_quad_real , 5.20084643035671309247102845217436109233532485096532359250266990497399469521525107083101563980212464972_quad_real)

!   PARAMETERS: CPML=800, r0=1300
       ! integer, parameter    :: IMODE        =            300
       ! complex(8), parameter :: ZLAMBDA_MODE =                 (8.0081126000624088d10, -28701.547450086537d00)
       ! complex(8), parameter :: ZCOEF10      =                 (1.d0,0.d0)
       ! complex(8), parameter :: ZCOEF01      =                 (3.6511740195643883d00, -0.0019136091580319757d00)

       ! integer, parameter    :: IMODE        =            301
       ! complex(8), parameter :: ZLAMBDA_MODE =                 (8.0031153354435613d10, -1.5529950356430494d06)
       ! complex(8), parameter :: ZCOEF10      =                 (-0.58683593789455668d00, 0.10540146192107115d00)
       ! complex(8), parameter :: ZCOEF01      =                 (0.d0,1.d0)

       ! integer, parameter    :: IMODE        =            302
       ! complex(8), parameter :: ZLAMBDA_MODE =                 (7.9978487898740714d10, -9.1995970868529847d06)
       ! complex(8), parameter :: ZCOEF10      =                 (1.d0,0.d0)
       ! complex(8), parameter :: ZCOEF01      =                 (11.124069938102021d00, -7.9270174198735518d00)
!
      contains
!
!----------------------------------------------------------------------------------------------------------------
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
!       compute factorial coefficients for Bessel-Frobenius series
        function Rfact(N) result(F)      
          use iso_fortran_env, wp => real64   
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
          F=1._wp
          if (N.eq.0) return
          do j=1,N
           F = F*2._wp/j
          enddo
!
        end function
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
          zc = (0._wp,0._wp)
          zc(0) = Zc0
          zc(1) = zc1
          zloc = zc(0) + zc(1)*dx
          Zdbess  = zc(1)
          Zd2bess = (0._wp,0._wp)
          zdloc_prev = 1._wp
          n=0
          do 
            zsum  = (0._wp,0._wp)
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
     7020   format('        Wavenum,R0,R           = ',3(e12.5,2x))
            write(*,7030) n,Zbess,Zdbess,zd2bess
     7030   format('        n,Zbess,Zdbess,Zd2bess = ',i3,3(2x,2e12.5))
    !!!!        call call pause
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
          zone = (1._wp,0._wp)
          zero = (0._wp,0._wp)
    !
    !  ...evaluate the solution at Reval
          call Bessel(ZLAMBDA_MODE_DP,zone,zero,Wavenum,R0,Reval, zbess10,zdbess10,zd2bess10)
          call Bessel(ZLAMBDA_MODE_DP,zero,zone,Wavenum,R0,Reval, zbess01,zdbess01,zd2bess01)
    !
    !  ...value
          Zval   = ZCOEF10_DP*zbess10   + ZCOEF01_DP*zbess01
          Zdval  = ZCOEF10_DP*zdbess10  + ZCOEF01_DP*zdbess01
          Zd2val = ZCOEF10_DP*zd2bess10 + ZCOEF01_DP*zd2bess01
          ! write(*,*) 'bessel_preset: zbess10,zbess01=',zbess10,zbess01
          ! write(*,*) 'bessel_preset: ZCOEF10,ZCOEF01=',ZCOEF10,ZCOEF01
          ! write(*,*) ''
  !
        end subroutine bessel_preset
!
!----------------------------------------------------------------------
!     compute factorial coefficients for rescaled Bessel-Frobenius series
        function Rfact_new(N,R0) result(F)
          implicit none
          integer, intent(in) :: N
          real(quad_real), intent(in) :: R0
          real(quad_real) :: F
          integer :: j
    !
          if (N.lt.0) then
            write(*,*) 'Rfact_new: N = ',N
            stop 1
          endif
          if (R0.le.0._quad_real) then
            write(*,*) 'Rfact_new: R0 = ',R0
            stop 1
          endif
    !
          F=1._quad_real
          if (N.eq.0) return
          do j=1,N
            F = F*2._quad_real/real(j,quad_real)/R0
          enddo
  !
        end function Rfact_new
!-----------------------------------------------------------------------
!
!   subroutine name    - Bessel_new
!
!-----------------------------------------------------------------------
!
!   latest revision    - Oct 25
!
!   purpose            - evaluate rescaled Bessel function of an arbitrary
!                        complex order using Frobenius method with the 
!                        expansion at zr=r0 (input zr is complex valued)
!                        
!   arguments
!     in:
!              Zmu     - complex order of the Bessel fucntion
!              Zc0,Zc1 - the first two complex coefficients in the
!                        Taylor expansion
!              Wavenum - real wave number (affected by local refractive index)
!              R0      - rescaling factor
!              Zr      - complex argument
!     out:
!              Zbess   - complex value of the Bessel function
!              Zdbess  - its derivative in x
!              Zd2bess - its second derivative in x
!
!-----------------------------------------------------------------------
!
        subroutine Bessel_new(Zmu,Zc0,Zc1,Wavenum,R0,Zr, &
                              Zbess,Zdbess,Zd2bess)
    !
          implicit none
          integer :: Idec
          complex(quad_complex), intent(in)  :: Zmu,Zc0,Zc1
          real(quad_real),       intent(in)  :: Wavenum,R0
          complex(quad_complex), intent(in)  :: Zr
          complex(quad_complex), intent(out) :: Zbess,Zdbess,Zd2bess
          complex(quad_complex)              :: zc(0:2000), zx, zsum,zloc,zdloc,zdloc_prev,zb(0:2000)
          integer  :: n,j,iprint
          real(quad_real) :: eps,aux
    !
          if (R0.le. 0._quad_real) then
            write(*,*) 'Bessel_new: R0 = ',R0
            stop 1
          endif
          eps = 10._quad_real**(-32)
    !
          iprint=0
       10 continue
          if (iprint.eq.3) then
            write(*,8100) Zmu, Zc0,Zc1, Wavenum,R0,Zr
     8100   format('Bessel_new: Zmu    = ',2e22.15,/,&
                   '            Zc0,Zc1      = ',2(2e22.15,2x),/,&
                   '            Wavenum,R0,Zr = ',e22.15,2x,e22.15,2x,2e22.15)
            call pause
          endif
    !
    !     transform input coordinate into x
          zx = R0*log(Zr/R0)

          aux = Wavenum**2
    !
          zc = (0._quad_real,0._quad_real)
          zc(0) = Zc0
          zc(1) = Zc1
          zloc = zc(0) + zc(1)*zx
          Zdbess = zc(1)
          Zd2bess = (0._quad_real,0._quad_real)
          zdloc_prev = 1._quad_real
          n=0
          do 
            ! zsum  = (0._quad_real,0._quad_real)
            ! do j=0,n
            !   zsum = zsum + rfact_new(n-j,R0)*zc(j)
            ! enddo

            zsum = zc(0)
            do j=0,n-1
              zsum = zsum * 2._quad_real/(R0*(n-j)) + zc(j+1)
            enddo

            zc(n+2) = (Zmu*zc(n) - aux*zsum)/real((n+1)*(n+2),16)
            zdloc = zc(n+2)*zx**(n+2)
            if (iprint.eq.3) then
              write(*,*)'Bessel_new: n, zdloc = ',n, zdloc
            endif
            zloc = zloc + zdloc
            Zdbess = Zdbess + zc(n+2)*zx**(n+1)*(n+2._quad_real)
            Zd2bess = Zd2bess + zc(n+2)*zx**(n)   *real((n+1)*(n+2),16)
            if ((abs(zdloc).lt.eps).and.(abs(zdloc_prev).lt.eps)) exit
            n=n+1
            if (n+2.gt.2000) then
              write(*,*) 'Bessel_new 1: n = ',n
              write(*,*) '             zr = ',zr
              write(*,*) '             zx = ',zx
              iprint=1
              call pause
              go to 10
            endif
            zdloc_prev = zdloc
          enddo
    !!!      write(*,*) 'Bessel_new 1: n = ',n
       
    !     pass to derivatives w.r.t r   
          Zd2bess = Zd2bess*(R0/Zr)**2 - Zdbess*R0/Zr**2
          Zdbess  = Zdbess*R0/Zr
          Zbess   = zloc
    !
          if (iprint.eq.1) then
            write(*,7110) Zmu,Zc0,Zc1
            write(*,7120) Wavenum,R0,Zr
            write(*,7130) n,Zbess,Zdbess,zd2bess
     7110   format('Bessel_new: Zlambda,Zc0,Zc1        = ',3(2e12.5,2x))
     7120   format('            Wavenum,R0,Zr          = ',2(e12.5,2x),2(2x,e12.5))
     7130   format('            n,Zbess,Zdbess,Zd2bess = ',i3,3(2x,2e12.5))
          endif
  !
        end subroutine Bessel_new

!----------------------------------------------------------------------------
!
!
!
!
! 
!
!
!
! !----------------------------------------------------------------------------
! !
!       subroutine bessel_stepindex_preset( Wavenum0,Wavenum1, R0, A, Zr, Zval, Zdval, Zd2val )
!       ! 
!       implicit none
! !
!       real(8), intent(in)  :: Wavenum0,Wavenum1,R0,A
!       complex(wp), intent(in)  :: Zr
!       complex(wp), intent(out) :: Zval,Zdval,Zd2val
!       complex(wp) :: zbess10,zdbess10,zd2bess10,zbess01,zdbess01,zd2bess01, &
!                      zone,zero,zmu,zcoef10_tmp,zcoef01_tmp
! !
!       zone = (1._quad_real,0._quad_real)
!       zero = (0._quad_real,0._quad_real)

      
! !
! !----------------------------------------------------------------------------
!       if (Zr%re.le.R0-A) then

!         ! write(*,*) 'bessel_stepindex_preset:   LEFT    R0,A,Zr=',R0,A,Zr

!         ! first, get the coefficients for this cladding side by evaluating the core at the interface
!         zmu = ZLAMBDA_MODE/R0**2
!         call Bessel_new(zmu,zone,zero,Wavenum0,R0,(R0-A)*zone, zbess10,zdbess10,zd2bess10)
!         call Bessel_new(zmu,zero,zone,Wavenum0,R0,(R0-A)*zone, zbess01,zdbess01,zd2bess01)
!         !  ...save values into the new coefficients
!         zcoef10_tmp  = ZCOEF10*zbess10   + ZCOEF01*zbess01
!         zcoef01_tmp  = ZCOEF10*zdbess10  + ZCOEF01*zdbess01

!         ! now, evaluate with the corresponding base point for this cladding side
!         zmu = ZLAMBDA_MODE/(R0-A)**2  ! adjusted zmu
!         call Bessel_new(zmu,zone,zero,Wavenum1,R0-A,Zr, zbess10,zdbess10,zd2bess10)
!         call Bessel_new(zmu,zero,zone,Wavenum1,R0-A,Zr, zbess01,zdbess01,zd2bess01)
!         !  ...value
!         Zval   = zcoef10_tmp*zbess10   + zcoef01_tmp*zbess01
!         Zdval  = zcoef10_tmp*zdbess10  + zcoef01_tmp*zdbess01
!         Zd2val = zcoef10_tmp*zd2bess10 + zcoef01_tmp*zd2bess01
!       elseif (Zr%re.lt.R0+A) then

!         ! write(*,*) 'bessel_stepindex_preset:   CORE    R0,A,Zr=',R0,A,Zr

!         zmu = ZLAMBDA_MODE/R0**2
!         call Bessel_new(zmu,zone,zero,Wavenum0,R0,Zr, zbess10,zdbess10,zd2bess10)
!         call Bessel_new(zmu,zero,zone,Wavenum0,R0,Zr, zbess01,zdbess01,zd2bess01)
!         !  ...value
!         Zval   = ZCOEF10*zbess10   + ZCOEF01*zbess01
!         Zdval  = ZCOEF10*zdbess10  + ZCOEF01*zdbess01
!         Zd2val = ZCOEF10*zd2bess10 + ZCOEF01*zd2bess01
!       else

!         ! write(*,*) 'bessel_stepindex_preset:   RIGHT   R0,A,Zr=',R0,A,Zr/

!         ! first, get the coefficients for this cladding side by evaluating the core at the interface
!         zmu = ZLAMBDA_MODE/R0**2
!         call Bessel_new(zmu,zone,zero,Wavenum0,R0,(R0+A)*zone, zbess10,zdbess10,zd2bess10)
!         call Bessel_new(zmu,zero,zone,Wavenum0,R0,(R0+A)*zone, zbess01,zdbess01,zd2bess01)
!         !  ...save values into the new coefficients
!         zcoef10_tmp  = ZCOEF10*zbess10   + ZCOEF01*zbess01
!         zcoef01_tmp  = ZCOEF10*zdbess10  + ZCOEF01*zdbess01

!         zmu = ZLAMBDA_MODE/(R0+A)**2  ! adjusted zmu
!         call Bessel_new(zmu,zone,zero,Wavenum1,R0+A,Zr, zbess10,zdbess10,zd2bess10)
!         call Bessel_new(zmu,zero,zone,Wavenum1,R0+A,Zr, zbess01,zdbess01,zd2bess01)
!         !  ...value
!         Zval   = zcoef10_tmp*zbess10   + zcoef01_tmp*zbess01
!         Zdval  = zcoef10_tmp*zdbess10  + zcoef01_tmp*zdbess01
!         Zd2val = zcoef10_tmp*zd2bess10 + zcoef01_tmp*zd2bess01
!       endif
!
!       end subroutine
!
!----------------------------------------------------------------------------
!
!
!     Evaluation based on continuation a principle. The Bessel coefficients in the inner parts of the left and right cladding regions
!     are found from the core region evaluated at the core-interface points (as developed in the Bessel paper), and, starting from 
!     the midpoint of each cladding region, we evaluate the the solution and obtain the coefficients to rescale again the Bessel functions.
!     The eigenvalue ZLAMBDA, and all the coefficients ZCOEF01, ZCOEF10 (core, L1, L2, R1, R2) are precomputed and set as parameters above
!     in this module. Make sure to uncomment the correct set of parameters by checking the corresponding value of the IMODE variable
!     in each block, which must match the ISOL variable given as an input.
!
!
      subroutine bessel_stepindex_preset( Wavenum0_dp,Wavenum1_dp, R0_dp, A_dp, B_dp, Zr_dp, Zval_dp, Zdval_dp, Zd2val_dp )
      ! 
      implicit none
!
      real(8), intent(in)  :: Wavenum0_dp,Wavenum1_dp,R0_dp,A_dp,B_dp
      complex(8), intent(in)  :: Zr_dp
      complex(8), intent(out) :: Zval_dp,Zdval_dp,Zd2val_dp

      ! upgraded I/O variables
      real(quad_real)       :: Wavenum0,Wavenum1,r0,a,b
      complex(quad_complex) :: zr,zval,Zdval,Zd2val
      ! additional variables
      real(quad_real)       :: abmid
      complex(quad_complex) :: zbess10,zdbess10,zd2bess10,zbess01,zdbess01,zd2bess01, &
                               zone,zero,zmu
!
      zone = (1._quad_real,0._quad_real)
      zero = (0._quad_real,0._quad_real)

      wavenum0 = real(Wavenum0_dp, quad_real)
      wavenum1 = real(Wavenum1_dp, quad_real)
      a = real(A_dp, quad_real)
      b = real(B_dp, quad_real)
      r0 = real(R0_dp, quad_real)
      zr = cmplx(Zr_dp%re, Zr_dp%im , quad_complex)



      ! write(*,*) 'bessel_stepindex_preset:'
      ! write(*,*) 'quad_real,quad_complex=',quad_real,quad_complex
      ! write(*,*) 'a_dp, a =',a_dp,a
      ! write(*,*) 'b_dp, b =',b_dp,b
      ! write(*,*) 'r0_dp, r0 =',r0_dp,r0
      ! write(*,*) 'zr_dp, zr =',zr_dp,zr
      ! write(*,*) 'ZLAMBDA_MODE=',ZLAMBDA_MODE
!
!----------------------------------------------------------------------------
      abmid = 0.5_quad_real*(a+b)
      
      if (Zr%re.le.R0-abmid) then

        ! now, evaluate with the corresponding base point for this cladding side, outer part (L2)
        zmu = ZLAMBDA_MODE/(R0-abmid)**2  ! adjusted zmu
        call Bessel_new(zmu,zone,zero,Wavenum1,R0-abmid,Zr, zbess10,zdbess10,zd2bess10)
        call Bessel_new(zmu,zero,zone,Wavenum1,R0-abmid,Zr, zbess01,zdbess01,zd2bess01)
        !  ...value
        Zval   = ZCOEF10_L2*zbess10   + ZCOEF01_L2*zbess01
        Zdval  = ZCOEF10_L2*zdbess10  + ZCOEF01_L2*zdbess01
        Zd2val = ZCOEF10_L2*zd2bess10 + ZCOEF01_L2*zd2bess01
      elseif(Zr%re.le.R0-A) then
        ! now, evaluate with the corresponding base point for this cladding side, inner part (L1)
        zmu = ZLAMBDA_MODE/(R0-A)**2  ! adjusted zmu
        call Bessel_new(zmu,zone,zero,Wavenum1,R0-A,Zr, zbess10,zdbess10,zd2bess10)
        call Bessel_new(zmu,zero,zone,Wavenum1,R0-A,Zr, zbess01,zdbess01,zd2bess01)
        !  ...value
        Zval   = ZCOEF10_L1*zbess10   + ZCOEF01_L1*zbess01
        Zdval  = ZCOEF10_L1*zdbess10  + ZCOEF01_L1*zdbess01
        Zd2val = ZCOEF10_L1*zd2bess10 + ZCOEF01_L1*zd2bess01
      elseif (Zr%re.lt.R0+A) then

        ! write(*,*) 'bessel_stepindex_preset:   CORE    R0,A,Zr=',R0,A,Zr

        zmu = ZLAMBDA_MODE/R0**2
        call Bessel_new(zmu,zone,zero,Wavenum0,R0,Zr, zbess10,zdbess10,zd2bess10)
        call Bessel_new(zmu,zero,zone,Wavenum0,R0,Zr, zbess01,zdbess01,zd2bess01)
        !  ...value
        Zval   = ZCOEF10*zbess10   + ZCOEF01*zbess01
        Zdval  = ZCOEF10*zdbess10  + ZCOEF01*zdbess01
        Zd2val = ZCOEF10*zd2bess10 + ZCOEF01*zd2bess01

      elseif (Zr%re.lt.R0+abmid) then

        ! now, evaluate with the corresponding base point for this cladding side, inner part (R1)
        zmu = ZLAMBDA_MODE/(R0+A)**2  ! adjusted zmu
        call Bessel_new(zmu,zone,zero,Wavenum1,R0+A,Zr, zbess10,zdbess10,zd2bess10)
        call Bessel_new(zmu,zero,zone,Wavenum1,R0+A,Zr, zbess01,zdbess01,zd2bess01)
        !  ...value
        Zval   = ZCOEF10_R1*zbess10   + ZCOEF01_R1*zbess01
        Zdval  = ZCOEF10_R1*zdbess10  + ZCOEF01_R1*zdbess01
        Zd2val = ZCOEF10_R1*zd2bess10 + ZCOEF01_R1*zd2bess01
      else

        ! now, evaluate with the corresponding base point for this cladding side, outer part (R2)
        zmu = ZLAMBDA_MODE/(R0+abmid)**2  ! adjusted zmu
        call Bessel_new(zmu,zone,zero,Wavenum1,R0+abmid,Zr, zbess10,zdbess10,zd2bess10)
        call Bessel_new(zmu,zero,zone,Wavenum1,R0+abmid,Zr, zbess01,zdbess01,zd2bess01)
        !  ...value
        Zval   = ZCOEF10_R2*zbess10   + ZCOEF01_R2*zbess01
        Zdval  = ZCOEF10_R2*zdbess10  + ZCOEF01_R2*zdbess01
        Zd2val = ZCOEF10_R2*zd2bess10 + ZCOEF01_R2*zd2bess01
      endif

      Zval_dp = cmplx(Zval%re, Zval%im , 8)
      Zdval_dp = cmplx(Zdval%re, Zdval%im , 8)
      Zd2val_dp = cmplx(Zd2val%re, zd2val%im , 8)

      ! write(*,*) 'zval_dp, zval =',zval_dp,zval
      ! write(*,*) 'zdval_dp, zdval =',zdval_dp,zdval
      ! write(*,*) 'zd2val_dp, zd2val =',zd2val_dp,zd2val
      ! call pause

      end subroutine


      end module