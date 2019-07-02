subroutine exact(Xp, Icase, &
     ZvalH,ZdvalH,Zd2valH, &
     ZvalE,ZdvalE,Zd2valE, &
     ZvalV,ZdvalV,Zd2valV, &
     ZvalQ,ZdvalQ,Zd2valQ)

  use thermo_elast
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

  real*8 :: x,y,z, fx,fy,fz, dfx, dfy, dfz
  
  integer :: iprint, ivar
  !-------------------------------------------
  iprint = 0

  ZvalH = ZERO; ZdvalH = ZERO; Zd2valH = ZERO;
  ZvalE = ZERO; ZdvalE = ZERO; Zd2valE = ZERO;
  ZvalV = ZERO; ZdvalV = ZERO; Zd2valV = ZERO;
  ZvalQ = ZERO; ZdvalQ = ZERO; Zd2valQ = ZERO;

  x = Xp(1); y = Xp(2); z = Xp(3);

  ! a polynomial solution
  fx = x**NPX
  if (NPX.eq.0) then
     dfx = 0.d0
  else
     dfx = NPX*x**(NPX-1)
  endif
  fy = y**NPY
  if (NPY.eq.0) then
     dfy = 0.d0
  else
     dfy = NPY*y**(NPY-1)
  endif
  fz = z**NPZ
  if (NPZ.eq.0) then
     dfz = 0.d0
  else
     dfz = NPZ*z**(NPZ-1)
  endif
  ZvalH( 1  ) = fx*fy*fz
  ZdvalH(1,1) = dfx*fy*fz
  ZdvalH(1,2) = fx*dfy*fz
  ZdvalH(1,3) = fx*fy*dfz
  
  if (iprint.eq.1) then
     write(*,7001)Icase
7001 format(' exact: Icase = ',i2)
     write(*,7002)Xp
7002 format(' exact: Xp = ',3f8.3)
     do ivar=1,1
        write(*,7003) ivar, ZvalH(ivar),ZdvalH(ivar,1:3)
7003    format(' exact: ivar,ZvalH(ivar),ZdvalH(ivar,1:3) = ', &
             i2,2x,e12.5,2x,3e12.5)
     enddo
     call pause
  endif
  
end subroutine exact
