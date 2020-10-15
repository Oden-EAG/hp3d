!---------------------------------------------------------------------------------------
!> Purpose : display user-defined quantity
!! @param[in] Mdle   - element (middle node) number
!! @param[in] Xi     - master element coordinates
!! @param[in] X      - physical coordinates
!! @param[in] Rn     - outward normal unit vector
!! @param[in] SolH  - H1    sol
!! @param[in] GradH - H1    grad
!! @param[in] SolE  - Hcurl sol
!! @param[in] CurlE - Hcurl curl
!! @param[in] SolV  - Hdiv  sol
!! @param[in] DivV  - Hdiv  div
!! @param[in] SolQ  - L2    sol
!!
!! @param[out] val   - quantity to display
!---------------------------------------------------------------------------------------
!
subroutine soldis(Mdle,Xi,X,Rn,SolH,GradH,SolE,CurlE,SolV,DivV,SolQ, Val)
  use common_prob_data, only : IEXACT_DISP, ICHOOSE_DISP
  use data_structure3D
!---------------------------------------------------------------------------------------
  implicit none
  integer,                      intent(in)  :: Mdle
  real*8,dimension(3),          intent(in)  :: Xi,X,Rn
  real*8,dimension(  MAXEQNH  ),intent(in)  :: SolH
  real*8,dimension(  MAXEQNH,3),intent(in)  :: GradH
  real*8,dimension(3,MAXEQNE  ),intent(in)  :: SolE
  real*8,dimension(3,MAXEQNE  ),intent(in)  :: CurlE
  real*8,dimension(3,MAXEQNV  ),intent(in)  :: SolV
  real*8,dimension(  MAXEQNV  ),intent(in)  :: DivV
  real*8,dimension(  MAXEQNQ  ),intent(in)  :: SolQ
  real*8,                       intent(out) :: Val
!---------------------------------------------------------------------------------------
! exact solution workspace
  integer                         :: icase
  real*8,dimension(  MAXEQNH    ) ::   valH
  real*8,dimension(  MAXEQNH,3  ) ::  dvalH
  real*8,dimension(  MAXEQNH,3,3) :: d2valH
  real*8,dimension(3,MAXEQNE    ) ::   valE
  real*8,dimension(3,MAXEQNE,3  ) ::  dvalE
  real*8,dimension(3,MAXEQNE,3,3) :: d2valE
  real*8,dimension(3,MAXEQNV    ) ::   valV
  real*8,dimension(3,MAXEQNV,3  ) ::  dvalV
  real*8,dimension(3,MAXEQNV,3,3) :: d2valV
  real*8,dimension(  MAXEQNQ    ) ::   valQ
  real*8,dimension(  MAXEQNQ,3  ) ::  dvalQ
  real*8,dimension(  MAXEQNQ,3,3) :: d2valQ
!---------------------------------------------------------------------------------------
  real*8 :: tmpH(3),tmpV(3,3),diff
  integer :: m
!---------------------------------------------------------------------------------------
!
  select case (IEXACT_DISP)
  case(1)
    icase=0
    call exact(X,icase, valH,dvalH,d2valH,  &
                        valE,dvalE,d2valE,  &
                        valV,dvalV,d2valV,  &
                        valQ,dvalQ,d2valQ)
    ! ** Exact  solution
    tmpH = valH

 !    write(*,*) 'valH = ', valH
 !    write(*,*) ''
 !    write(*,*) 'solH = ', solH
 !    write(*,*) ''
 !    write(*,*) 'valV = '
 !    write(*,1000) valV
 !    write(*,*) ''
 !    write(*,*) 'solV = '
 !    write(*,1000) solV
 ! 1000 format(3(/,3e14.7))

    ! call pause

    tmpV = valV

  case(2)
    ! ** Approximate solution
    tmpH = solH

    tmpV = solV

  case(3)
    icase=0
    call exact(X,icase, valH,dvalH,d2valH,  &
                        valE,dvalE,d2valE,  &
                        valV,dvalV,d2valV,  &
                        valQ,dvalQ,d2valQ)

    ! ** Solution error
    tmpH = valH-solH

    tmpV = valV-solV

  end select

  select case (ICHOOSE_DISP)
  case(1,2,3)
    Val = tmpH(ICHOOSE_DISP)
  case(4)
    Val = dsqrt(tmpH(1)**2 + tmpH(2)**2 + tmpH(3)**2)
  case(5)
    Val = tmpV(1,1)
  case(6)
    Val = tmpV(1,2)

    diff = tmpV(1,2)-tmpV(2,1)
    if (abs(diff) .lt. 1.d-12) then
      write(*,6000) diff
 6000 format('soldis: sigma_12-sigma_21 = ',e12.5)
    endif
  case(7)
    Val = tmpV(1,3)

    diff = tmpV(1,3)-tmpV(3,1)
    if (abs(diff) .lt. 1.d-12) then
      write(*,6001) diff
 6001 format('soldis: sigma_13-sigma_31 = ',e12.5)
    endif

  case(8)
    Val = tmpV(2,2)
  case(9)
    Val = tmpV(2,3)

    diff = tmpV(2,3)-tmpV(3,2)
    if (abs(diff) .lt. 1.d-12) then
      write(*,6002) diff
 6002 format('soldis: sigma_23-sigma_32 = ',e12.5)
    endif

  case(10)
    Val = tmpV(3,3)
  case(11)
    Val = dsqrt(  tmpV(1,1)**2 + tmpV(1,2)**2 + tmpV(1,3)**2  &
                + tmpV(2,1)**2 + tmpV(2,2)**2 + tmpV(2,3)**2  &
                + tmpV(3,1)**2 + tmpV(3,2)**2 + tmpV(3,3)**2  )
  end select

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
  integer :: iprev
!---------------------------------------------------------------------------------------
!
  write(*,*) 'USE PREVIOUS SELECTION? 1-Y,0-N'
  read(*,*) iprev
  if (iprev.eq.1) then
    call disp_soldis(NSTD_OUT) ; return
  endif
  !
  write(*,*) '1) EXACT SOLUTION  2) APPROXIMATE SOLUTION  3) SIGNED ERROR '
  read(*,*) IEXACT_DISP
10  write(*,*) 'SET VARIABLE     1) u1    2) u2    3) u3    4) |u|'
    write(*,*) '                 5) sigma_11    6) sigma_12    7) sigma_13'
    write(*,*) '                                8) sigma_22    9) sigma_23'
    write(*,*) '                                              10) sigma_33'
    write(*,*) '                11) |sigma| '
  read(*,*) ICHOOSE_DISP
  select case(ICHOOSE_DISP)
  case(1,2,3,4,5,6,7,8,9,10,11)
  case default
    go to 10
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
  write(Nstream,200)
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
  case(10); write(Nstream,310)
  case(11); write(Nstream,311)
  end select
  write(Nstream,200)

100 format('DISPLAY SETUP')
101 format('EXACT SOLUTION IS CHOSEN')
102 format('APPROX SOLUTION IS CHOSEN')
103 format('(EXACT SOLUTION - APPROX SOLUTION) IS CHOSEN')
301 format('u1 IS CHOSEN')
302 format('u2 IS CHOSEN')
303 format('u3 IS CHOSEN')
304 format('|u| IS CHOSEN')
305 format('sigma_11 IS CHOSEN')
306 format('sigma_12 IS CHOSEN')
307 format('sigma_13 IS CHOSEN')
308 format('sigma_22 IS CHOSEN')
309 format('sigma_23 IS CHOSEN')
310 format('sigma_33 IS CHOSEN')
311 format('|sigma| IS CHOSEN')
200 format('-----------------------')

end subroutine disp_soldis