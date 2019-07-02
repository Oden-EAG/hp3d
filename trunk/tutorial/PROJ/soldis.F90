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
subroutine soldis(Mdle,Xi,X,Rn, &
     ZsolH,ZgradH, &
     ZsolE,ZcurlE, &
     ZsolV,ZdivV, &
     ZsolQ, &
     Val)
  use proj
  use data_structure3D
!---------------------------------------------------------------------------------------
  implicit none
  integer,                      intent(in)  :: Mdle
  real*8,dimension(3),          intent(in)  :: Xi,X,Rn
  real*8,dimension(  MAXEQNH  ),intent(in)  :: ZsolH
  real*8,dimension(  MAXEQNH,3),intent(in)  :: ZgradH
  real*8,dimension(3,MAXEQNE  ),intent(in)  :: ZsolE
  real*8,dimension(3,MAXEQNE  ),intent(in)  :: ZcurlE
  real*8,dimension(3,MAXEQNV  ),intent(in)  :: ZsolV
  real*8,dimension(  MAXEQNV  ),intent(in)  :: ZdivV
  real*8,dimension(  MAXEQNQ  ),intent(in)  :: ZsolQ
  real*8,                       intent(out) :: Val

  real*8 :: &
       zvalH(    MAXEQNH    ), &
       zdvalH(   MAXEQNH,3  ), &
       zd2valH(  MAXEQNH,3,3), &
       zvalE(  3,MAXEQNE    ), &
       zdvalE( 3,MAXEQNE,3  ), &
       zd2valE(3,MAXEQNE,3,3), &
       zvalV(  3,MAXEQNV    ), &
       zdvalV( 3,MAXEQNV,3  ), &
       zd2valV(3,MAXEQNV,3,3), &
       zvalQ(    MAXEQNQ    ), &
       zdvalQ(   MAXEQNQ,3  ), &
       zd2valQ(  MAXEQNQ,3,3)
  
!---------------------------------------------------------------------------------------
!
  select case (IEXACT_DISP)
  case(1)
     ! ** exact solution
     call exact(X,NODES(Mdle)%case, &
          zvalH,zdvalH,zd2valH, &
          zvalE,zdvalE,zd2valE, &
          zvalV,zdvalV,zd2valV, &
          zvalQ,zdvalQ,zd2valQ)
     Val = zvalH(1)

  case(2)
     ! ** Approximation solution
     Val = ZsolH(1)
  end select
!  
end subroutine soldis

!
!---------------------------------------------------------------------------------------
!> Purpose : show the quantities to display 
!---------------------------------------------------------------------------------------
subroutine soldis_select
  use proj
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
  write(*,*) 'SET VARIABLE     1) EXACT     2) APPROX'
  read(*,*) iprev
  IEXACT_DISP = iprev
  !
  call disp_soldis(NSTD_OUT)
!  
end subroutine soldis_select
