subroutine exact(Xp, Icase, &
     ZvalH,ZdvalH,Zd2valH, &
     ZvalE,ZdvalE,Zd2valE, &
     ZvalV,ZdvalV,Zd2valV, &
     ZvalQ,ZdvalQ,Zd2valQ)

  use parameters , only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ, ZERO
  use lapl
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
  integer :: iprint
  real*8 &
       fx,   fy,   fz,  &
       dfx,  dfy,  dfz, &
       d2fx, d2fy, d2fz
  real*8 &
       t, x(3), r, r3, dudr, &
       drdx, drdy, drdz, &
       d2rdx2, d2rdy2, d2rdz2, d2rdxy, d2rdxz, d2rdyz, &
       d2udr2
  iprint = 0

  ZvalH = ZERO; ZdvalH = ZERO; Zd2valH = ZERO
  ZvalE = ZERO; ZdvalE = ZERO; Zd2valE = ZERO
  ZvalV = ZERO; ZdvalV = ZERO; Zd2valV = ZERO
  ZvalQ = ZERO; ZdvalQ = ZERO; Zd2valQ = ZERO

  ! set exact solution
  select case(IEXACT_PROB)
  case(IEXACT_HARMONIC_PROB)
     fx   =  sin(Xp(1)); fy   =  cos(Xp(2)); fz  = Xp(3)
     dfx  =  cos(Xp(1)); dfy  = -sin(Xp(2)); dfz = 1.d0
     d2fx = -sin(Xp(1)); d2fy = -cos(Xp(2)); dfz = 0.d0

     ! solution
     ZvalH(1)    = fx *fy *fz

     ! gradient
     ZdvalH(1,1) = dfx*fy *fz
     ZdvalH(1,2) = fx *dfy*fz
     ZdvalH(1,3) = fx *fy *dfz

     ! hessian
     Zd2valH(1,1,1) = d2fx*fy  *fz
     Zd2valH(1,2,2) = fx  *d2fy*fz
     Zd2valH(1,3,3) = fx  *fy  *d2fz

     Zd2valH(1,1,2) = dfx*dfy*fz
     Zd2valH(1,1,3) = dfx*fy *dfz
     Zd2valH(1,2,3) = fx *dfy*dfz

     Zd2valH(1,2,1) = Zd2valH(1,1,2)
     Zd2valH(1,3,1) = Zd2valH(1,1,3)
     Zd2valH(1,3,2) = Zd2valH(1,2,3)

  case(IEXACT_SHOCK_PROB)

     x  = Xp - X0
     r  = sqrt(x(1)**2 + x(2)**2 + x(3)**2)
     r3 = r**3

     drdx = x(1)/r
     drdy = x(2)/r
     drdz = x(3)/r

     d2rdx2 = 1.d0/r - x(1)**2/r3
     d2rdy2 = 1.d0/r - x(2)**2/r3
     d2rdz2 = 1.d0/r - x(3)**2/r3

     d2rdxy = - x(1)*x(2)/r3
     d2rdxz = - x(1)*x(3)/r3
     d2rdyz = - x(2)*x(3)/r3

     ! solution
     t = ALPHA*(r-R0)
     ZvalH(1) = t
     
     ! gradient
     dudr = ALPHA/(1.d0+t**2)
     ZdvalH(1,1) =  dudr*drdx
     ZdvalH(1,2) =  dudr*drdy
     ZdvalH(1,3) =  dudr*drdz

     ! hessian
     d2udr2 = -2.d0*t*dudr**2
     Zd2valH(1,1,1) = d2udr2*drdx**2   + dudr*d2rdx2
     Zd2valH(1,2,2) = d2udr2*drdy**2   + dudr*d2rdy2
     Zd2valH(1,3,3) = d2udr2*drdz**2   + dudr*d2rdz2

     Zd2valH(1,1,2) = d2udr2*drdx*drdy + dudr*d2rdxy
     Zd2valH(1,1,3) = d2udr2*drdx*drdz + dudr*d2rdxz
     Zd2valH(1,2,3) = d2udr2*drdy*drdz + dudr*d2rdyz

     Zd2valH(1,2,1) = Zd2valH(1,1,2)
     Zd2valH(1,3,1) = Zd2valH(1,1,3)
     Zd2valH(1,3,2) = Zd2valH(1,2,3)

  end select

end subroutine exact

