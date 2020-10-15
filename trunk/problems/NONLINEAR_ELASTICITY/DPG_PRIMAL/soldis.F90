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
  use hyperelasticity, only : get_det_mat, find_material,eval_strain_energy_w_f, &
                              MATERIALS, LINEAR, DEL, get_svd_mat
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
  real*8 :: tmpH(3),tmpV(3),tmpF(3,3),tmpJ,tmpJx,w,dw(3,3),d2w(3,3,3,3), &
            cauchy(3,3), QL(3,3),strn(3,3),stretch(3),tmpLambda(3), &
            QR(3,3),tmpStrs(3),tmpPress
  integer :: m , imat
!---------------------------------------------------------------------------------------
!
  select case (IEXACT_DISP)
  case(1)
    icase=0
    call exact(X,Mdle,icase, valH,dvalH,d2valH,  &
                             valE,dvalE,d2valE,  &
                             valV,dvalV,d2valV,  &
                             valQ,dvalQ,d2valQ)
    ! ** Exact  solution
    tmpH = valH(4:6)

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

 !    call pause

    tmpV(1:3) = valV(1,4:6)*rn(1)  &
              + valV(2,4:6)*rn(2)  &
              + valV(3,4:6)*rn(3)

    tmpF=dvalH(4:6,:)
    tmpF(1,1)=tmpF(1,1)+1.d0
    tmpF(2,2)=tmpF(2,2)+1.d0
    tmpF(3,3)=tmpF(3,3)+1.d0
    call get_det_mat(tmpF,tmpJ)

    call find_material(Mdle,imat)

    call eval_strain_energy_w_f(imat,X,tmpF,w,dw,d2w)

    if (MATERIALS(imat)%CONSTIT.ne.LINEAR) then
      cauchy = tmpJ*matmul(dw,transpose(tmpF))
      strn = tmpF
    else
      cauchy = dw
      strn = tmpF !- DEL
      strn = 0.5d0*(strn+transpose(strn))
    endif

    call get_svd_mat(strn,QL,tmpLambda,QR)
    ! diagonalize cauchy stress tensor
    cauchy = matmul(cauchy,QL)
    cauchy = matmul(transpose(QL),cauchy)
    tmpStrs(1) = cauchy(1,1)
    tmpStrs(2) = cauchy(2,2)
    tmpStrs(3) = cauchy(3,3)

    tmpPress = -sum(tmpStrs)/3.d0

  case(2)
    ! ** Approximate solution
    tmpH = SolH(4:6)

    tmpV(1:3) = solV(1,4:6)*rn(1)  &
              + solV(2,4:6)*rn(2)  &
              + solV(3,4:6)*rn(3)


    tmpF = gradH(4:6,:)
    tmpF(1,1)=tmpF(1,1)+1.d0
    tmpF(2,2)=tmpF(2,2)+1.d0
    tmpF(3,3)=tmpF(3,3)+1.d0
    call get_det_mat(tmpF,tmpJ)


    call find_material(Mdle,imat)

    call eval_strain_energy_w_f(imat,X,tmpF,w,dw,d2w)

    if (MATERIALS(imat)%CONSTIT.ne.LINEAR) then
      cauchy = tmpJ*matmul(dw,transpose(tmpF))
      strn = tmpF
    else
      cauchy = dw
      strn = tmpF !- DEL
      strn = 0.5d0*(strn+transpose(strn))
    endif

    call get_svd_mat(strn,QL,tmpLambda,QR)
    ! diagonalize cauchy stress tensor
    cauchy = matmul(cauchy,QL)
    cauchy = matmul(transpose(QL),cauchy)
    tmpStrs(1) = cauchy(1,1)
    tmpStrs(2) = cauchy(2,2)
    tmpStrs(3) = cauchy(3,3)

    tmpPress = -sum(tmpStrs)/3.d0


  case(3)
    icase=0
    call exact(X,Mdle,icase, valH,dvalH,d2valH,  &
                        valE,dvalE,d2valE,  &
                        valV,dvalV,d2valV,  &
                        valQ,dvalQ,d2valQ)


    call find_material(Mdle,imat)

    ! ** Solution error
    tmpH = valH(4:6)-SolH(4:6)

    tmpV(1:3) = (valV(1,4:6)-SolV(1,4:6))*rn(1)  &
              + (valV(2,4:6)-SolV(2,4:6))*rn(2)  &
              + (valV(3,4:6)-SolV(3,4:6))*rn(3)

    tmpF=dvalH(4:6,:)
    tmpF(1,1)=tmpF(1,1)+1.d0
    tmpF(2,2)=tmpF(2,2)+1.d0
    tmpF(3,3)=tmpF(3,3)+1.d0
    call get_det_mat(tmpF,tmpJx)

    call eval_strain_energy_w_f(imat,X,tmpF,w,dw,d2w)

    if (MATERIALS(imat)%CONSTIT.ne.LINEAR) then
      cauchy = tmpJx*matmul(dw,transpose(tmpF))
      strn = tmpF
    else
      cauchy = dw
      strn = tmpF !- DEL
      strn = 0.5d0*(strn+transpose(strn))
    endif

    call get_svd_mat(strn,QL,tmpLambda,QR)
    ! diagonalize cauchy stress tensor
    cauchy = matmul(cauchy,QL)
    cauchy = matmul(transpose(QL),cauchy)
    tmpStrs(1) = cauchy(1,1)
    tmpStrs(2) = cauchy(2,2)
    tmpStrs(3) = cauchy(3,3)

    tmpPress = -sum(tmpStrs)/3.d0



    tmpF=gradH(4:6,:)
    tmpF(1,1)=tmpF(1,1)+1.d0
    tmpF(2,2)=tmpF(2,2)+1.d0
    tmpF(3,3)=tmpF(3,3)+1.d0
    call get_det_mat(tmpF,tmpJ)

    call eval_strain_energy_w_f(imat,X,tmpF,w,dw,d2w)

    if (MATERIALS(imat)%CONSTIT.ne.LINEAR) then
      cauchy = tmpJ*matmul(dw,transpose(tmpF))
      strn = tmpF
    else
      cauchy = dw
      strn = tmpF !- DEL
      strn = 0.5d0*(strn+transpose(strn))
    endif

    call get_svd_mat(strn,QL,stretch,QR)
    ! diagonalize cauchy stress tensor
    cauchy = matmul(cauchy,QL)
    cauchy = matmul(transpose(QL),cauchy)
    tmpStrs(1) = tmpStrs(1) - cauchy(1,1)
    tmpStrs(2) = tmpStrs(2) - cauchy(2,2) 
    tmpStrs(3) = tmpStrs(3) - cauchy(3,3)

    tmpLambda = tmpLambda - stretch

    tmpPress = tmpPress-(-sum(tmpStrs)/3.d0)

    tmpJ=tmpJx-tmpJ



  end select

  select case (ICHOOSE_DISP)
  case(1,2,3)
    Val = tmpH(ICHOOSE_DISP)
  case(4)
    Val = dsqrt(tmpH(1)**2 + tmpH(2)**2 + tmpH(3)**2)
  case(5,6,7)
    Val = tmpV(ICHOOSE_DISP-4)
  case(8)
    Val = dsqrt(tmpV(1)**2 + tmpV(2)**2 + tmpV(3)**2)
  case(9)
    Val = tmpJ
  case(10)
    Val = tmpPress
  case(11)
    Val = tmpStrs(1)
  case(12)
    Val = tmpStrs(2)
  case(13)
    Val = tmpStrs(3)
  case(14)
    Val = tmpLambda(1)
  case(15)
    Val = tmpLambda(2)
  case(16)
    Val = tmpLambda(3)
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
10  write(*,*) 'SET VARIABLE  1) u1   2) u2   3) u3    4) |u|     5) t1'
    write(*,*) '              6) t2   7) t3   8) |t|   9) det F   10) Mean pressure'
    write(*,*) '             11) Principal Sigma 1    12) Principal Sigma 2    13) Principal Sigma 3'
    write(*,*) '             14) Principal stretch 1  15) Principal stretch 2  16) Principal stretch 3'
  read(*,*) ICHOOSE_DISP
  select case(ICHOOSE_DISP)
  case(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16)
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
  write(Nstream,320)
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
  case(12); write(Nstream,312)
  case(13); write(Nstream,313)
  case(14); write(Nstream,314)
  case(15); write(Nstream,315)
  case(16); write(Nstream,316)
  end select
  write(Nstream,320)

100 format('DISPLAY SETUP')
101 format('EXACT SOLUTION IS CHOSEN')
102 format('APPROX SOLUTION IS CHOSEN')
103 format('(EXACT SOLUTION - APPROX SOLUTION) IS CHOSEN')
301 format('u1 IS CHOSEN')
302 format('u2 IS CHOSEN')
303 format('u3 IS CHOSEN')
304 format('|u| IS CHOSEN')
305 format('t1 IS CHOSEN')
306 format('t2 IS CHOSEN')
307 format('t3 IS CHOSEN')
308 format('|t| IS CHOSEN')
309 format('det F IS CHOSEN')
310 format('Mean pressure IS CHOSEN')
311 format('Principal Sigma 1 IS CHOSEN')
312 format('Principal Sigma 2 IS CHOSEN')
313 format('Principal Sigma 3 IS CHOSEN')
314 format('Principal stretch 1 IS CHOSEN')
315 format('Principal stretch 2 IS CHOSEN')
316 format('Principal stretch 3 IS CHOSEN')

320 format('--------------------------------------------')

end subroutine disp_soldis