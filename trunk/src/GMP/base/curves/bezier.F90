      module bezier

      contains
!-------------------------------------------------------------------------------------- 
!
      subroutine curve_bezier(Nc,Eta, X,dXdEta,ddXddEta,dddXdddEta,ddddXddddEta)
!
!-------------------------------------------------------------------------------------- 
!     recursive relation for derivative of Berstein polynomials:
!       (B_i^n)' = n(B_(i-1)^(n-1) - B_i^n)
!
!
!-------------------------------------------------------------------------------------- 
      use GMP
!-------------------------------------------------------------------------------------- 
      implicit none
!-------------------------------------------------------------------------------------- 
!     DUMMY ARGUMENTS
      integer,intent(in)                            :: Nc
      real*8,intent(in)                           :: Eta
      real*8,dimension(3),intent(out)             :: X,dXdEta
      real*8,dimension(3), optional, intent(out)  :: ddXddEta
      real*8,dimension(3), optional, intent(out)  :: dddXdddEta
      real*8,dimension(3), optional, intent(out)  :: ddddXddddEta
!-------------------------------------------------------------------------------------- 
      real*8 :: poly0,poly1,poly2,poly3,poly4,dpoly
      integer  :: i,iprint,deg
!-------------------------------------------------------------------------------------- 
!      
      iprint = 0
!!      if (nc.eq.5) iprint=1
!      
      if (iprint.eq.1) then
        write(*,1000)Nc,Eta,CURVES(abs(Nc))%Type
 1000   format(' curve_Bezier: Nc = ',i4,'; Eta = ',e12.5,'; type = ',a12)
      endif
!
      select case(CURVES(abs(Nc))%Type)
!  ...quintic Bezier curve      
      case('5Bezier')
        deg = 5
!  ...sextic Bezier curve        
      case('6Bezier')
        deg = 6
!  ...septic Bezier curve        
      case('7Bezier')
        deg = 7
      case default
        write(*,1005)CURVES(abs(Nc))%Type
 1005   format('curve_Bezier: type not supported, type = ',a12)
        stop
      endselect 
!      
      X = 0.d0;  dXdEta = 0.d0; 
      if (present(ddXddEta))     ddXddEta     = 0.d0
      if (present(dddXdddEta))   dddXdddEta   = 0.d0
      if (present(ddddXddddEta)) ddddXddddEta = 0.d0
!  ...accumulate          
      do i = 0,deg
        call Bernstein_poly(i,deg  ,Eta, poly0,dpoly)
        call Bernstein_poly(i,deg-1,Eta, poly1,dpoly)
        call Bernstein_poly(i,deg-2,Eta, poly2,dpoly)
        call Bernstein_poly(i,deg-3,Eta, poly3,dpoly)
        call Bernstein_poly(i,deg-4,Eta, poly4,dpoly)
!  .....curve    
        X        = X        +             poly0* CURVES(abs(Nc))%Rdata(3*i    :3*i+2) 
!  .....velocity      
        if (i.ge.deg) cycle
        dXdEta   = dXdEta   + deg*        poly1*(CURVES(abs(Nc))%Rdata(3*(i+1):3*i+5) - &
                                                 CURVES(abs(Nc))%Rdata(3*i    :3*i+2))
!  .....acceleration               
        if (present(ddXddEta)) then
        if (i.ge.deg-1) cycle
        ddXddEta = ddXddEta + deg*(deg-1)*poly2*(CURVES(abs(Nc))%Rdata(3*(i+2):3*i+8) - &
                                            2.d0*CURVES(abs(Nc))%Rdata(3*(i+1):3*i+5) + &
                                                 CURVES(abs(Nc))%Rdata(3*i    :3*i+2))
        endif                                         
!  .....3rd derivative                                                
        if (present(dddXdddEta)) then
        if (i.ge.deg-2) cycle
        dddXdddEta = dddXdddEta + deg*(deg-1)*(deg-2)*poly3*                                 &
                                                (CURVES(abs(Nc))%Rdata(3*(i+3):3*(i+3)+2) -  &
                                            3.d0*CURVES(abs(Nc))%Rdata(3*(i+2):3*(i+2)+2) +  &
                                            3.d0*CURVES(abs(Nc))%Rdata(3*(i+1):3*(i+1)+2) -  & 
                                                 CURVES(abs(Nc))%Rdata(3* i   :3* i   +2))
        endif                                         
!  .....4th derivative                                                
        if (present(ddddXddddEta)) then
        if (i.ge.deg-3) cycle
        ddddXddddEta = ddddXddddEta + deg*(deg-1)*(deg-2)*(deg-3)*poly4*                     &
                                                (CURVES(abs(Nc))%Rdata(3*(i+4):3*(i+3)+2) -  &
                                            4.d0*CURVES(abs(Nc))%Rdata(3*(i+3):3*(i+3)+2) +  &
                                            6.d0*CURVES(abs(Nc))%Rdata(3*(i+2):3*(i+2)+2) -  &
                                            4.d0*CURVES(abs(Nc))%Rdata(3*(i+1):3*(i+1)+2) +  & 
                                                 CURVES(abs(Nc))%Rdata(3* i   :3* i   +2))
        endif      
!        
      enddo
!     
      if (iprint.eq.1) then
        write(*,1001)X
 1001   format(' curve_Bezier: X          = ',3(e12.5,2x))
        write(*,1002)dXdEta
 1002   format(' curve_Bezier: dXdEta     = ',3(e12.5,2x))        
        if (present(ddXddEta))    write(*,1003)ddXddEta
 1003   format(' curve_Bezier: ddXddEta   = ',3(e12.5,2x))
        if (present(dddXdddEta))  write(*,1004)dddXdddEta
 1004   format(' curve_Bezier: dddXdddEta = ',3(e12.5,2x))
      endif 
!
!
      end subroutine curve_bezier

      end module bezier
