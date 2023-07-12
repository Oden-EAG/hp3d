!---------------------------------------------------------------------------------------
!> Purpose : display user-defined quantity
!! @param[in] Mdle   - element (middle node) number
!! @param[in] Xi     - master element coordinates
!! @param[in] X      - physical coordinates
!! @param[in] Rn     - outward normal unit vector
!! @param[in] ZsolH  - H1    sol
!! @param[in] ZgradH - H1    grad
!! @param[in] ZsolE  - Hcurl sol
!! @param[in] ZcurlE - Hcurl curl
!! @param[in] ZsolV  - Hdiv  sol
!! @param[in] ZdivV  - Hdiv  div
!! @param[in] ZsolQ  - L2    sol
!!
!! @param[out] val   - quantity to display
!---------------------------------------------------------------------------------------
!
subroutine soldis(Mdle,Xi,X,Rn,ZsolH,ZgradH,ZsolE,ZcurlE,ZsolV,ZdivV,ZsolQ, Val)
  use common_prob_data_UW, only : IEXACT_DISP, ICHOOSE_DISP
  use data_structure3D
  use control
  use parameters
!---------------------------------------------------------------------------------------
  implicit none
  integer,                      intent(in)  :: Mdle
  real*8,dimension(3),          intent(in)  :: Xi,X,Rn
  complex*16,dimension(  MAXEQNH  ),intent(in)  :: ZsolH
  complex*16,dimension(  MAXEQNH,3),intent(in)  :: ZgradH
  complex*16,dimension(3,MAXEQNE  ),intent(in)  :: ZsolE
  complex*16,dimension(3,MAXEQNE  ),intent(in)  :: ZcurlE
  complex*16,dimension(3,MAXEQNV  ),intent(in)  :: ZsolV
  complex*16,dimension(  MAXEQNV  ),intent(in)  :: ZdivV
  complex*16,dimension(  MAXEQNQ  ),intent(in)  :: ZsolQ
  real*8,                           intent(out) :: Val
!---------------------------------------------------------------------------------------
! exact solution workspace
  integer                         :: icase
  complex*16,dimension(  MAXEQNH    ) ::   zvalH
  complex*16,dimension(  MAXEQNH,3  ) ::  zdvalH
  complex*16,dimension(  MAXEQNH,3,3) :: zd2valH
  complex*16,dimension(3,MAXEQNE    ) ::   zvalE
  complex*16,dimension(3,MAXEQNE,3  ) ::  zdvalE
  complex*16,dimension(3,MAXEQNE,3,3) :: zd2valE
  complex*16,dimension(3,MAXEQNV    ) ::   zvalV
  complex*16,dimension(3,MAXEQNV,3  ) ::  zdvalV
  complex*16,dimension(3,MAXEQNV,3,3) :: zd2valV
  complex*16,dimension(  MAXEQNQ    ) ::   zvalQ
  complex*16,dimension(  MAXEQNQ,3  ) ::  zdvalQ
  complex*16,dimension(  MAXEQNQ,3,3) :: zd2valQ
  complex*16                          :: zval
  real*8                              :: dsqrt
  integer                             :: ifl, ndtype, nflag, icomp
!---------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!

  
   call decode(iabs(ICHOOSE_DISP), ifl,ndtype)
   call decode(ifl, nflag,icomp)



   if (NEXACT .ne. 0) then
      icase=0
      call exact(X,icase, zvalH,zdvalH,zd2valH, &
                          zvalE,zdvalE,zd2valE, &
                          zvalV,zdvalV,zd2valV, &
                          zvalQ,zdvalQ,zd2valQ)
   else
      zvalH = ZERO; zdvalH = ZERO; zd2valH = ZERO
      zvalE = ZERO; zdvalE = ZERO; zd2valE = ZERO
      zvalV = ZERO; zdvalV = ZERO; zd2valV = ZERO
      zvalQ = ZERO; zdvalQ = ZERO; zd2valQ = ZERO
   endif   


   select case(ndtype)
!  ...trace of pressure (H1 variable)
      case(1)
        select case(nflag)
!
!     ...FE solution
         case(0); zval = ZsolH(icomp)
!
!     ...exact solution
         case(1); zval = zvalH(icomp)
!
!     ...error function
         case(2); zval = zvalH(icomp)-ZsolH(icomp)
      end select
!      
!  ...trace of velocity (H(div) variable)
      case(3)
         select case(nflag)
!
!     ...FE solution
         case(0); zval = ZsolV(icomp,1)
!
!     ...exact solution
         case(1); zval = zvalV(icomp,1)
!
!     ...error function
         case(2); zval = zvalV(icomp,1)-ZsolV(icomp,1)
      end select
!  ...L2 variables
      case(4)
!
      select case(nflag)
!
!      ...FE solution
         case(0); zval = ZsolQ(icomp)
!
!      ...exact solution
          case(1); zval = zvalQ(icomp)
!
!      ...error function
          case(2); zval = zvalQ(icomp)-ZsolQ(icomp)
      case default
         write(*,7001) ICHOOSE_DISP
 7001    format('soldis: UNKNOWN FLAG, ICHOOSE_DISP = ',i5)
         stop 1
      end select
!      
   case default
      write(*,7001) ICHOOSE_DISP
   end select
! 
   if (ICHOOSE_DISP.gt.0) then
      Val = dreal(zval)
   else
      Val = aimag(zval)
   endif  



!
end subroutine soldis
!
!
!
!---------------------------------------------------------------------------------------
!> Purpose : show the quantities to display
!---------------------------------------------------------------------------------------
   subroutine soldis_select
!      
   use parameters, only : NSTD_OUT
   use common_prob_data_UW, only : ICHOOSE_DISP, IEXACT_DISP
   use control
!   
!---------------------------------------------------------------------------------------
   implicit none
!---------------------------------------------------------------------------------------
!
   100 continue 

   write(*,*) 'SPECIFY THE QUANTITY TO DISPLAY '
   write(*,*) 'USE -/+ FOR IMAGINARY/REAL PART'
   write(*,*) '  '
   write(*,*) ' \hat{p}.....................................11'
   if (NEXACT.ne.0) then
      write(*,*) ' EXACT \hat{p}..............................111'
      write(*,*) ' ERROR \hat{p}..............................211'
      write(*,*) '  '
   endif
!    
   write(*,*) ' \hat{u_x}...................................13'
   if (NEXACT.ne.0) then
      write(*,*) ' EXACT \hat{u_x}............................113'
      write(*,*) ' ERROR \hat{u_x}............................213'
      write(*,*) '  '
   endif
!    
   write(*,*) ' \hat{u_y}...................................23'
   if (NEXACT.ne.0) then
      write(*,*) ' EXACT \hat{u_y}............................123'
      write(*,*) ' ERROR \hat{u_y}............................223'
      write(*,*) '  '
   endif
!
   write(*,*) ' \hat{u_z}...................................33'
   if (NEXACT.ne.0) then
      write(*,*) ' EXACT \hat{u_z}............................133'
      write(*,*) ' ERROR \hat{u_z}............................233'
      write(*,*) '  '
   endif   
!   
   write(*,*) ' p...........................................14'
   if (NEXACT.ne.0) then
      write(*,*) ' EXACT p....................................114'
      write(*,*) ' ERROR in p.................................214'
      write(*,*) '  '
   endif
!    
   write(*,*) ' u_x.........................................24'
   if (NEXACT.ne.0) then
      write(*,*) ' EXACT u_x..................................124'
      write(*,*) ' ERROR in u_x...............................224'
      write(*,*) '  '
   endif
!    
   write(*,*) ' u_y.........................................34'
   if (NEXACT.ne.0) then
      write(*,*) ' EXACT u_y..................................134'
      write(*,*) ' ERROR in u_y...............................234'
      write(*,*) '  '
   endif
!
   write(*,*) ' u_z.........................................44'
   if (NEXACT.ne.0) then
      write(*,*) ' EXACT u_z..................................144'
      write(*,*) ' ERROR in u_z...............................244'
      write(*,*) '  '
   endif   
!    


   read(*,*) ICHOOSE_DISP
      select case(ICHOOSE_DISP)
      case( 11, 111, 211, 13, 113, 213, 23, 123, 223, 33, 133, 233,       &
            14, 114, 214, 24, 124, 224, 34, 134, 234, 334, 44, 144, 244)
      case(-11,-111,-211,-13,-113,-213,-23,-123,-223,-33,-133,-233,       &
           -14,-114,-214,-24,-124,-224,-34,-134,-234,-334,-44,-144,-244)
      case default
        write(*,*) 'soldis_select: WRONG ICHOOSE_DISP = ',ICHOOSE_DISP
        go to 100
      end select
  !
  call disp_soldis(NSTD_OUT)
!
!
!---------------------------------------------------------------------------------------
!> Purpose : Soldis display interface
!---------------------------------------------------------------------------------------
   end subroutine soldis_select

   subroutine disp_soldis(Nstream)
   use common_prob_data_UW, only : ICHOOSE_DISP
!---------------------------------------------------------------------------------------
   implicit none
   integer, intent(in) :: Nstream
   integer             :: ndtype, ifl, nflag, icomp
!---------------------------------------------------------------------------------------

   write(Nstream,100)
   write(Nstream,4000)



   call decode(iabs(ICHOOSE_DISP), ifl,ndtype)
   call decode(ifl, nflag,icomp)

   select case(ndtype)
   case(1)
!  ...trace of pressure (H1 variable)
      write(Nstream,1001)
      select case(nflag)
      case(0);
         write(Nstream,1000)
         if (ICHOOSE_DISP .lt. 0) write(Nstream,2000)
         if (ICHOOSE_DISP .gt. 0) write(Nstream,3000)
      case(1)
         write(Nstream,1010)
         if (ICHOOSE_DISP .lt. 0) write(Nstream,2000)
         if (ICHOOSE_DISP .gt. 0) write(Nstream,3000)
      case(2)
         write(Nstream,1020)
         if (ICHOOSE_DISP .lt. 0) write(Nstream,2000)
         if (ICHOOSE_DISP .gt. 0) write(Nstream,3000)
      end select   
   case(3)
!  ...trace of velocity (Hdiv)
      write(Nstream,1102)
      select case(nflag)
      case(0);
         write(Nstream,1000)
         if (ICHOOSE_DISP .lt. 0) write(Nstream,2000)
         if (ICHOOSE_DISP .gt. 0) write(Nstream,3000)
      case(1)
         write(Nstream,1010)
         if (ICHOOSE_DISP .lt. 0) write(Nstream,2000)
         if (ICHOOSE_DISP .gt. 0) write(Nstream,3000)
      case(2)
         write(Nstream,1020)
         if (ICHOOSE_DISP .lt. 0) write(Nstream,2000)
         if (ICHOOSE_DISP .gt. 0) write(Nstream,3000)
      end select
   case(4)
!  ...L2 variables
      write(Nstream,1202)
      select case(nflag)
      case(0);
         write(Nstream,1000)
         if (ICHOOSE_DISP .lt. 0) write(Nstream,2000)
         if (ICHOOSE_DISP .gt. 0) write(Nstream,3000)
      case(1)
         write(Nstream,1010)
         if (ICHOOSE_DISP .lt. 0) write(Nstream,2000)
         if (ICHOOSE_DISP .gt. 0) write(Nstream,3000)
      case(2)
         write(Nstream,1020)
         if (ICHOOSE_DISP .lt. 0) write(Nstream,2000)
         if (ICHOOSE_DISP .gt. 0) write(Nstream,3000)
      end select
   end select   

      write(Nstream,4000)
!
100 format('DISPLAY SETUP')
1001 format('TRACE OF PRESSURE (H1 VARIABLE)')
1102 format('TRACE OF VELOCITY (H(DIV) VARIABLE)')
1202 format('L2 VARIABLES')
1000 format('APPROXIMATE SOLUTION IS CHOSEN')
1010 format('EXACT SOLUTION IS CHOSEN')
1020 format('ERROR OF THE SOLUTION IS CHOSEN')
2000 format('DISPLAYING THE IMAGINARY PART')
3000 format('DISPLAYING THE REAL PART')
4000 format('----------------------')
!
end subroutine disp_soldis
