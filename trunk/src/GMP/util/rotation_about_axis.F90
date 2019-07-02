!----------------------------------------------------------------------
!> @Purpose - routine rotates a point R by an angle Alpha about an axis 
!!            Ax through point Pt
!
!> @param[in]    Ax         - axis
!> @param[in]    Pt         - a point on the axis
!> @param[in]    Alpha      - angle of rotation
!> @param[inout] R          - point to be rotated / rotated point
!> @param[out]   dR_dAlpha  - derivative 
!
!> @revision Feb 13
!----------------------------------------------------------------------
!
subroutine rotation_about_axis(Ax,Pt,Alpha, R,dR_dAlpha)
!
      implicit none
      real*8,dimension(3),intent(in)    :: Ax,Pt
      real*8             ,intent(in)    :: Alpha
      real*8,dimension(3),intent(inout) :: R
      real*8,dimension(3),intent(out)   :: dR_dAlpha
!
      real*8 :: u,v,w,a,b,c,x,y,rz,d,s1,s2,s3,s4
!      
      u=Ax(1)  ;  v=Ax(2)  ;  w=Ax(3)
      a=Pt(1)  ;  b=Pt(2)  ;  c=Pt(3)
      x=R(1)   ;  y=R(2)   ;  rz=R(3)
!
      R(1) = (a*(v**2+w**2)+u*(-b*v-c*w+u*x+v*y+w*rz) +                 &
              ((x-a)*(v**2+w**2)+u*(b*v+c*w-v*y-w*rz))*cos(Alpha) +     &
              sqrt(u**2+v**2+w**2)*(b*w-c*v-w*y+v*rz)* sin(Alpha))/     &
              (u**2+v**2+w**2)
      R(2) = (b*(u**2+w**2)+v*(-a*u-c*w+u*x+v*y+w*rz) +                 &
              ((y-b)*(u**2+w**2)+v*(a*u+c*w-u*x-w*rz))*cos(Alpha) +     &
              sqrt(u**2+v**2+w**2)*(-a*w+c*u+w*x-u*rz)*sin(Alpha))/     &
              (u**2+v**2+w**2)
       R(3) = (c*(u**2+v**2)+w*(-a*u-b*v+u*x+v*y+w*rz) +                &
              ((rz-c)*(u**2+v**2)+w*(a*u+b*v-u*x-v*y))*cos(Alpha) +     &
               sqrt(u**2+v**2+w**2)*(a*v-b*u-v*x+u*y)* sin(Alpha))/     &
               (u**2+v**2+w**2) 
!
      dR_dAlpha(1) = (                                                  &
              ((x-a)*(v**2+w**2)+u*(b*v+c*w-v*y-w*rz))*(-sin(Alpha)) +  &
              sqrt(u**2+v**2+w**2)*(b*w-c*v-w*y+v*rz)* cos(Alpha))/     &
              (u**2+v**2+w**2)
      dR_dAlpha(2) = (                                                  &
              ((y-b)*(u**2+w**2)+v*(a*u+c*w-u*x-w*rz))*(-sin(Alpha)) +  &
              sqrt(u**2+v**2+w**2)*(-a*w+c*u+w*x-u*rz)*cos(Alpha))/     &
              (u**2+v**2+w**2)
      dR_dAlpha(3) = (                                                  &
              ((rz-c)*(u**2+v**2)+w*(a*u+b*v-u*x-v*y))*(-sin(Alpha)) +  &
               sqrt(u**2+v**2+w**2)*(a*v-b*u-v*x+u*y)* cos(Alpha))/     &
               (u**2+v**2+w**2) 
!
!
endsubroutine rotation_about_axis
