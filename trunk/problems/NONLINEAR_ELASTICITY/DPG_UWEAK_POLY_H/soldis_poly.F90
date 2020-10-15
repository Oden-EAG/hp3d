!---------------------------------------------------------------------------------------
!> Purpose : display user-defined quantity - Adjusted to PolyDPG
!! @param[in] Mdle   - element (middle node) number
!! @param[in] Mdlf   - face global number
!! @param[in] Xi     - master element coordinates
!! @param[in] X      - physical coordinates
!! @param[in] Rn     - outward normal unit vector
!! @param[in] ZsolQ  - Face  H1  sol
!! @param[in] ZsolQ  - Face  L2  sol
!! @param[in] ZsolQ  - Field L2  sol
!!
!! @param[out] val   - quantity to display
!---------------------------------------------------------------------------------------
!
subroutine soldis(Mdle,Mdlf,Xi,X,Rn,ZsolU,ZsolF,ZsolQ, Val)
  use common_prob_data, only : IEXACT_DISP, ICHOOSE_DISP
  use data_structure3D_poly
!---------------------------------------------------------------------------------------
  implicit none
  integer,                      intent(in)  :: Mdle,Mdlf
  real*8,dimension(3),          intent(in)  :: Xi,X,Rn
  real*8,dimension(  MAXEQNU  ),intent(in)  :: ZsolU
  real*8,dimension(  MAXEQNF  ),intent(in)  :: ZsolF
  real*8,dimension(  MAXEQNQ  ),intent(in)  :: ZsolQ
  real*8,                       intent(out) :: Val
!---------------------------------------------------------------------------------------
! exact solution workspace
  integer                         :: icase
  real*8,dimension(  MAXEQNH    ) ::   zvalH
  real*8,dimension(  MAXEQNH,3  ) ::  zdvalH
  real*8,dimension(  MAXEQNH,3,3) :: zd2valH
  real*8,dimension(3,MAXEQNE    ) ::   zvalE
  real*8,dimension(3,MAXEQNE,3  ) ::  zdvalE
  real*8,dimension(3,MAXEQNE,3,3) :: zd2valE
  real*8,dimension(3,MAXEQNV    ) ::   zvalV
  real*8,dimension(3,MAXEQNV,3  ) ::  zdvalV
  real*8,dimension(3,MAXEQNV,3,3) :: zd2valV
  real*8,dimension(  MAXEQNU    ) ::   zvalU
  real*8,dimension(  MAXEQNF    ) ::   zvalF
  real*8,dimension(  MAXEQNQ    ) ::   zvalQ
  real*8,dimension(  MAXEQNQ,3  ) ::  zdvalQ
  real*8,dimension(  MAXEQNQ,3,3) :: zd2valQ
!---------------------------------------------------------------------------------------
  real*8 :: tmpVal
!---------------------------------------------------------------------------------------
!
  select case (IEXACT_DISP)
  case(1)
    icase=0
    call exact(X,Rn,icase,zvalH,zdvalH,zd2valH, &
                          zvalE,zdvalE,zd2valE, &
                          zvalV,zdvalV,zd2valV, &
                          zvalU,zvalF,zvalQ,zdvalQ,zd2valQ)
    ! ** Exact  solution
    select case (ICHOOSE_DISP)
    case(1,2,3,4,5,6)
      tmpVal = zvalQ(ICHOOSE_DISP)
    case(7)
      tmpVal = zvalF(1)
    case(8)
      tmpVal = zvalF(2)
    case(9)
      tmpVal = zvalF(3)
    end select 
  case(2)
    ! ** Approximate solution
    select case (ICHOOSE_DISP)
    case(1,2,3,4,5,6)
      tmpVal = ZsolQ(ICHOOSE_DISP)
    case(7)
      tmpVal = ZsolF(1)
    case(8)
      tmpVal = ZsolF(2)
    case(9)
      tmpVal = ZsolF(3)
    end select 
  case(3)
    icase=0
    call exact(X,Rn,icase,zvalH,zdvalH,zd2valH, &
                          zvalE,zdvalE,zd2valE, &
                          zvalV,zdvalV,zd2valV, &
                          zvalU,zvalF,zvalQ,zdvalQ,zd2valQ)
    ! ** Solution error
    select case (ICHOOSE_DISP)
    case(1,2,3,4,5,6)
      tmpVal = zvalQ(ICHOOSE_DISP)-ZsolQ(ICHOOSE_DISP)
    case(7)
      tmpVal = zvalF(1)-ZsolF(1)
    case(8)
      tmpVal = zvalF(2)-ZsolF(2)
    case(9)
      tmpVal = zvalF(3)-ZsolF(3)
    end select 
  end select
! 
  Val = tmpVal
!
end subroutine soldis
!
!
!
!---------------------------------------------------------------------------------------
!> Purpose : show the quantities to display
!---------------------------------------------------------------------------------------
subroutine soldis_select
  use parameters, only : NSTD_OUT
  use common_prob_data, only : ICHOOSE_DISP, IEXACT_DISP
!---------------------------------------------------------------------------------------
  implicit none
!   integer :: iprev
! !---------------------------------------------------------------------------------------
! !
!   write(*,*) 'USE PREVIOUS SELECTION? 1-Y,0-N'
!   read(*,*) iprev
!   if (iprev.eq.1) then
!     call disp_soldis(NSTD_OUT) ; return
!   endif
  !
  write(*,*) 'EXACT SOLUTION (1) or APPROXIMATE SOLUTION (2) or SIGNED ERROR (3)?'
  read(*,*) IEXACT_DISP
    write(*,200)
200 format('SET VARIABLE:  (1) u_1   (2) u_2   (3) u_3   (4) sigma_11    (5) sigma_22',  &
                        '  (6) sigma_33   (7) sigma_n_1    8) sigma_n_2    9) sigma_n_3' )
  read(*,*) ICHOOSE_DISP
  select case(ICHOOSE_DISP)
  case(1,2,3,4,5,6,7,8,9)
  case default
    write(*,200)
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
  use common_prob_data, only : IEXACT_DISP, ICHOOSE_DISP
!---------------------------------------------------------------------------------------
  implicit none
  integer, intent(in) :: Nstream
!---------------------------------------------------------------------------------------

  write(Nstream,100)
  write(Nstream,310)
  select case (IEXACT_DISP)
  case(1); write(Nstream,101)
  case(2); write(Nstream,102)
  case(3); write(Nstream,103)
  end select
  select case (ICHOOSE_DISP)
  case(1); write(Nstream,301)
  case(2); write(Nstream,302)
  case(3); write(Nstream,303)
  case(4); write(Nstream,304)
  case(5); write(Nstream,305)
  case(6); write(Nstream,306)
  case(7); write(Nstream,307)
  case(8); write(Nstream,308)
  case(9); write(Nstream,309)
  end select
  write(Nstream,310)
!
100 format('DISPLAY SETUP')
101 format('EXACT SOLUTION IS CHOSEN')
102 format('APPROX SOLUTION IS CHOSEN')
103 format('(EXACT SOLUTION - APPROX SOLUTION) IS CHOSEN')
301 format('u1 IS CHOSEN')
302 format('u2 IS CHOSEN')
303 format('u3 IS CHOSEN')
304 format('sigma_11 IS CHOSEN')
305 format('sigma_22 IS CHOSEN')
306 format('sigma_33 IS CHOSEN')
307 format('sigma_n_1 IS CHOSEN')
308 format('sigma_n_2 IS CHOSEN')
309 format('sigma_n_3 IS CHOSEN')
310 format('-----------------------')
!
end subroutine disp_soldis