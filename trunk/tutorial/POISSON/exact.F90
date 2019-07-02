subroutine exact(Xp, Icase, &
     ZvalH,ZdvalH,Zd2valH, &
     ZvalE,ZdvalE,Zd2valE, &
     ZvalV,ZdvalV,Zd2valV, &
     ZvalQ,ZdvalQ,Zd2valQ)

  use parameters , only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ, ZERO

  implicit none
  !-------------------------------------------
  real*8, dimension(3), intent(in)  :: Xp
  integer,              intent(in)  :: Icase

  real*8,               intent(out) :: &
       ZvalH(MAXEQNH), &
       ZdvalH(MAXEQNH,3),Zd2valH(MAXEQNH,3,3), &
       ZvalE(3,MAXEQNE), &
       ZdvalE(3,MAXEQNE,3),Zd2valE(3,MAXEQNE,3,3), &
       ZvalV(3,MAXEQNV), &
       ZdvalV(3,MAXEQNV,3),Zd2valV(3,MAXEQNV,3,3), &
       ZvalQ(MAXEQNQ), &
       ZdvalQ(MAXEQNQ,3),Zd2valQ(MAXEQNQ,3,3)

  iprint = 0

  ZvalH = ZERO; ZdvalH = ZERO; Zd2valH = ZERO
  ZvalE = ZERO; ZdvalE = ZERO; Zd2valE = ZERO
  ZvalV = ZERO; ZdvalV = ZERO; Zd2valV = ZERO
  ZvalQ = ZERO; ZdvalQ = ZERO; Zd2valQ = ZERO

end subroutine exact
