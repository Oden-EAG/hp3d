subroutine exact(Xp, Icase, &
     ZvalH,ZdvalH,Zd2valH, &
     ZvalE,ZdvalE,Zd2valE, &
     ZvalV,ZdvalV,Zd2valV, &
     ZvalQ,ZdvalQ,Zd2valQ)
  use bioem
  use prob_plane
  implicit none
  !-------------------------------------------
  real*8, dimension(3), intent(in)  :: Xp
  integer,              intent(in)  :: Icase

  complex*16,           intent(out) :: &
       ZvalH(MAXEQNH), &
       ZdvalH(MAXEQNH,3),Zd2valH(MAXEQNH,3,3), &
       ZvalE(3,MAXEQNE), &
       ZdvalE(3,MAXEQNE,3),Zd2valE(3,MAXEQNE,3,3), &
       ZvalV(3,MAXEQNV), &
       ZdvalV(3,MAXEQNV,3),Zd2valV(3,MAXEQNV,3,3), &
       ZvalQ(MAXEQNQ), &
       ZdvalQ(MAXEQNQ,3),Zd2valQ(MAXEQNQ,3,3)

  !-------------------------------------------
  integer    :: ivar, itype, iprint
  real*8     :: omega
  complex*16 :: zvalE_loc(3), zdvalE_loc(3,3)
  !-------------------------------------------

  iprint = 0

  ZvalH = Z_0; ZdvalH = Z_0; Zd2valH = Z_0
  ZvalE = Z_0; ZdvalE = Z_0; Zd2valE = Z_0
  ZvalV = Z_0; ZdvalV = Z_0; Zd2valV = Z_0
  ZvalQ = Z_0; ZdvalQ = Z_0; Zd2valQ = Z_0

  call get_prob_type(itype)
  select case(itype)
  case(PROB_PLANE_INCIDENT)
     call get_angular_velocity(omega)
     call get_plane_wave(Xp, omega, zvalE_loc(1:3), zdvalE_loc(1:3,1:3) )

     do ivar=1, MAXEQNE
        ZvalE(1:3,ivar)      = zvalE_loc(1:3)
        ZdvalE(1:3,ivar,1:3) = zdvalE_loc(1:3,1:3)
     enddo
  case(PROB_PLANE_SCATTER)
     
  end select

  if (iprint.eq.1) then
     write(*,*) 'Xp    = ', Xp
     write(*,*) 'zvalE = ', zvalE_loc(1:3)
     write(*,*) 'zdvalE = ', zdvalE_loc(1:3,1:3)
  endif

end subroutine exact
