   subroutine test_exact
   use common_prob_data, ONLY : pi, OMEGA
   use parameters , only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ, ZERO

   integer                          :: Icase
   real*8, dimension(3)             :: Xp, Xr1, Xr2

!..exact solution
!..exact solution
   complex*16,dimension(  MAXEQNH    ) ::   ValH
   complex*16,dimension(  MAXEQNH,3  ) ::  DvalH
   complex*16,dimension(  MAXEQNH,3,3) :: d2valH
   complex*16,dimension(3,MAXEQNE    ) ::   ValE
   complex*16,dimension(3,MAXEQNE,3  ) ::  DvalE
   complex*16,dimension(3,MAXEQNE,3,3) :: d2valE
   complex*16,dimension(3,MAXEQNV    ) ::   ValV
   complex*16,dimension(3,MAXEQNV,3  ) ::  DvalV
   complex*16,dimension(3,MAXEQNV,3,3) :: d2valV
   complex*16,dimension(  MAXEQNQ    ) ::   valQ
   complex*16,dimension(  MAXEQNQ,3  ) ::  dvalQ
   complex*16,dimension(  MAXEQNQ,3,3) :: d2valQ
   complex*16 :: p, gradp(3), grad2p(3,3)

   complex*16 :: temp1, dtemp1,d2temp1  
   complex*16 :: temp2, dtemp2,d2temp2  
   real*8 :: h, h_r
   real*8 :: w0, z_R
!..functions and their derivatives
   real*8     :: w,w_z,w_zz
   real*8     :: f,f_z,f_zz
   real*8     :: r,r_x, r_xx,r_y, r_yy, r_yx
   real*8     :: A,A_x,A_xx,A_xy,A_xz,A_y,A_yx,A_yy,A_yz,A_z,A_zx,A_zy,A_zz
   real*8     :: g,g_x,g_xx,g_xy,g_xz,g_y,g_yx,g_yy,g_yz,g_z,g_zx,g_zy,g_zz
   real*8     :: phi, phi_z, phi_zz
   real*8     :: D, D_z, D_zz
   real*8     :: C,C_x,C_xx,C_xy,C_xz,C_y,C_yx,C_yy,C_yz,C_z,C_zx,C_zy,C_zz
   real*8     :: B,B_x,B_xx,B_xy,B_xz,B_y,B_yx,B_yy,B_yz,B_z,B_zx,B_zy,B_zz
   real*8     :: z,x,y

!
   z = 0.3d0
   x = 0.4d0
   y = 0.6d0

   do i = 1, 10
      h = 2.0d0**(-i)     

!  ...check first derivatives
!  ...derivative with respect to x
      Xp = (/x,y,z/)
      call acoustics_solution(Xp, p, gradp, grad2p)
      dtemp1 = gradp(1)

      Xp = (/x+h,y,z/)
      call acoustics_solution(Xp, p, gradp, grad2p)
      temp1 = p

      Xp = (/x-h,y,z/)
      Xr2 = Xp
      call acoustics_solution(Xp, p, gradp, grad2p)
      temp2 = p
      write(*,*) 'dp_x diff = ',  abs((temp1-temp2)/(2.d0*h) - dtemp1)
!
!  ...derivative with respect to y
      Xp = (/x,y,z/)
      call acoustics_solution(Xp, p, gradp, grad2p)
      dtemp1 = gradp(2)

      Xp = (/x,y+h,z/)
      call acoustics_solution(Xp, p, gradp, grad2p)
      temp1 = p

      Xp = (/x,y-h,z/)
      call acoustics_solution(Xp, p, gradp, grad2p)
      temp2 = p

      write(*,*) 'dp_y diff = ',  abs((temp1-temp2)/(2.0*h) - dtemp1)
!
!
!  ...derivative with respect to z
      Xp = (/x,y,z/)
      call acoustics_solution(Xp, p, gradp, grad2p)
      dtemp1 = gradp(3)

      Xp = (/x,y,z+h/)
      call acoustics_solution(Xp, p, gradp, grad2p)
      temp1 = p

      Xp = (/x,y,z-h/)
      call acoustics_solution(Xp, p, gradp, grad2p)
      temp2 = p

      write(*,*) 'dp_z diff = ',  abs((temp1-temp2)/(2.0*h) - dtemp1)

!  ...check second derivatives
!  ...p_xx
      Xp = (/x,y,z/)
      call exact(Xp,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)

      d2temp1 = d2valQ(1,1,1)

      Xp = (/x+h,y,z/)
      call exact(Xp,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)

      dtemp1 = dvalQ(1,1)

      Xp = (/x-h,y,z/)
      call exact(Xp,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      dtemp2 = dvalQ(1,1)

      write(*,*) 'dp_xx diff = ',  abs((dtemp1-dtemp2)/(2.0*h) - d2temp1)

      Xp = (/x,y,z/)
      call acoustics_solution(Xp, p, gradp, grad2p)
      d2temp1 = grad2p(1,1)

      Xp = (/x+h,y,z/)
      call acoustics_solution(Xp, p, gradp, grad2p)
      dtemp1 = gradp(1)

      Xp = (/x-h,y,z/)
      call acoustics_solution(Xp, p, gradp, grad2p)
      dtemp2 = gradp(1)

      write(*,*) '2:dp_xx diff = ',  abs((dtemp1-dtemp2)/(2.0*h) - d2temp1)

!  ...p_xy
      Xp = (/x,y,z/)
      call exact(Xp,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      d2temp1 = d2valQ(1,1,2)

      Xp = (/x,y+h,z/)
      call exact(Xp,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      dtemp1 = dvalQ(1,1)

      Xp = (/x,y-h,z/)
      call exact(Xp,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      dtemp2 = dvalQ(1,1)

      write(*,*) 'dp_xy diff = ',  abs((dtemp1-dtemp2)/(2.0*h) - d2temp1)

!  ...p_xz
      Xp = (/x,y,z/)
      call exact(Xp,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      d2temp1 = d2valQ(1,1,3)

      Xp = (/x,y,z+h/)
      call exact(Xp,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      dtemp1 = dvalQ(1,1)

      Xp = (/x,y,z-h/)
      call exact(Xp,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      dtemp2 = dvalQ(1,1)

      write(*,*) 'dp_xz diff = ',  abs((dtemp1-dtemp2)/(2.0*h) - d2temp1)

!  ...p_yx
      Xp = (/x,y,z/)
      call exact(Xp,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      d2temp1 = d2valQ(1,2,1)

      Xp = (/x+h,y,z/)
      call exact(Xp,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      dtemp1 = dvalQ(1,2)

      Xp = (/x-h,y,z/)
      call exact(Xp,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      dtemp2 = dvalQ(1,2)

      write(*,*) 'dp_yx diff = ',  abs((dtemp1-dtemp2)/(2.0*h) - d2temp1)

!  ...p_yy
      Xp = (/x,y,z/)
      call exact(Xp,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      d2temp1 = d2valQ(1,2,2)

      Xp = (/x,y+h,z/)
      call exact(Xp,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      dtemp1 = dvalQ(1,2)

      Xp = (/x,y-h,z/)
      call exact(Xp,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      dtemp2 = dvalQ(1,2)

      write(*,*) 'dp_yy diff = ',  abs((dtemp1-dtemp2)/(2.0*h) - d2temp1)

!  ...p_yz
      Xp = (/x,y,z/)
      call exact(Xp,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      d2temp1 = d2valQ(1,2,3)

      Xp = (/x,y,z+h/)
      call exact(Xp,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      dtemp1 = dvalQ(1,2)

      Xp = (/x,y,z-h/)
      call exact(Xp,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      dtemp2 = dvalQ(1,2)

      write(*,*) 'dp_yz diff = ',  abs((dtemp1-dtemp2)/(2.0*h) - d2temp1)

!  ...p_zx
      Xp = (/x,y,z/)
      call exact(Xp,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      d2temp1 = d2valQ(1,3,1)

      Xp = (/x+h,y,z/)
      call exact(Xp,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      dtemp1 = dvalQ(1,3)

      Xp = (/x-h,y,z/)
      call exact(Xp,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      dtemp2 = dvalQ(1,3)

      write(*,*) 'dp_zx diff = ',  abs((dtemp1-dtemp2)/(2.0*h) - d2temp1)

!  ...p_zy
      Xp = (/x,y,z/)
      call exact(Xp,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      d2temp1 = d2valQ(1,3,2)

      Xp = (/x,y+h,z/)
      ! call acoustics_solution(Xp, p, gradp, grad2p)
      ! dtemp1 = gradp(3)
      call exact(Xp,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      dtemp1 = dvalQ(1,3)

      Xp = (/x,y-h,z/)
      call exact(Xp,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      dtemp2 = dvalQ(1,3)

      write(*,*) 'dp_zy diff = ',  abs((dtemp1-dtemp2)/(2.0*h) - d2temp1)

!  ...p_zz
      Xp = (/x,y,z/)
      call exact(Xp,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      d2temp1 = d2valQ(1,3,3)

      Xp = (/x,y,z+h/)
      call exact(Xp,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      dtemp1 = dvalQ(1,3)

      Xp = (/x,y,z-h/)
      call exact(Xp,Icase, ValH,DvalH,d2valH, ValE,DvalE,d2valE,   &
                         ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
      dtemp2 = dvalQ(1,3)

      write(*,*) 'dp_zz diff = ',  abs((dtemp1-dtemp2)/(2.0*h) - d2temp1)      

   enddo

   end subroutine test_exact


