!-----------------------------------------------------------------------------------
!
!     routine name      - solution
!
!-----------------------------------------------------------------------------------
!
!     latest revision:  - July 2019
!
!     purpose:          - compute all relevant quantities of the exact solutions
!
!     arguments:
!        in:
!             X         - a point in physical space
!             Icase     - node case (specifies what variables are supported)
!        out:
!             u         - value of the solution at X
!             gradu     - corresponding first derivatives at X
!             gradgradu - corresponding second derivatives at X
!
!-----------------------------------------------------------------------------------
subroutine solution(X, u,gradu,gradgradu)
!
   use data_structure3D
   use common_prob_data, only : PI,IEXACT_PROB
   use parameters      , only : ZERO
!
   implicit none
!
   real*8, dimension(3),   intent(in)  :: X
   real*8,                 intent(out) :: u
   real*8, dimension(3),   intent(out) :: gradu     ! 1st derivative - gradient
   real*8, dimension(3,3), intent(out) :: gradgradu ! 2nd derivative - Hessian
!
   real*8 :: x1,x2,x3,f_x,f_y,f_z,df_x,df_y,df_z,ddf_x,ddf_y,ddf_z
   real*8 :: np_x,np_y,np_z
   integer :: isol,isol_p
!
!--------------------------------------------------------------------------------
!
   x1 = X(1); x2 = X(2); x3 = X(3)
!
   u = ZERO; gradu = ZERO; gradgradu = ZERO
!
!..select exact solution
   select case (IEXACT_PROB)
!
!  ...polynomial solution
      case(0)
!     ...set polynomial order of manufactured solution
         isol_p = 6
         np_x = real(isol_p,8)
         np_y = real(isol_p,8)
         np_z = real(isol_p,8)
!     ...value
         f_x = x1**np_x
         f_y = x2**np_y
         f_z = x3**np_z
!     ...derivatives
         select case(int(np_x))
            case(0); df_x = 0.d0; ddf_x = 0.d0
            case(1); df_x = 1.d0; ddf_x = 0.d0
            case default
               df_x = np_x * x1**(np_x-1.d0)
               ddf_x = np_x * (np_x-1.d0) * x1**(np_x-2.d0)
         end select
         select case(int(np_y))
            case(0); df_y = 0.d0; ddf_y = 0.d0
            case(1); df_y = 1.d0; ddf_y = 0.d0
            case default
               df_y = np_y * x2**(np_y-1.d0)
               ddf_y = np_y * (np_y-1.d0) * x2**(np_y-2.d0)
         end select
         select case(int(np_z))
            case(0); df_z = 0.d0; ddf_z = 0.d0
            case(1); df_z = 1.d0; ddf_z = 0.d0
            case default
               df_z = np_z * x3**(np_z-1.d0)
               ddf_z = np_z * (np_z-1.d0) * x3**(np_z-2.d0)
         end select
!     ...value
         u = f_x*f_y*f_z
!     ...1st order derivatives
         gradu(1)=  df_x *   f_y *   f_z
         gradu(2)=   f_x *  df_y *   f_z
         gradu(3)=   f_x *   f_y *  df_z
!     ...2nd order derivatives
         gradgradu(1,1) = ddf_x *   f_y *   f_z
         gradgradu(1,2) =  df_x *  df_y *   f_z
         gradgradu(1,3) =  df_x *   f_y *  df_z
         gradgradu(2,1) =  gradgradu(1,2)
         gradgradu(2,2) =   f_x * ddf_y *   f_z
         gradgradu(2,3) =   f_x *  df_y *  df_z
         gradgradu(3,1) =  gradgradu(1,3)
         gradgradu(3,2) =  gradgradu(2,3)
         gradgradu(3,3) =   f_x *   f_y * ddf_z
!
!  ...smooth sin sinh z solution on unit cube
!  ...non-vanishing boundary only on x=1,y=1,z=1
!  ...uniform in z
      case(1)
         u = dsin(x1*PI)*dsinh(x2*PI)
         gradu(1) = PI*dcos(x1*PI)*dsinh(x2*PI)
         gradu(2) = PI*dsin(x1*PI)*dcosh(x2*PI)
         gradu(3) = ZERO
         gradgradu(1,1) = -PI*PI*dsin(x1*PI)*dsinh(x2*PI)
         gradgradu(1,2) =  PI*PI*dcos(x1*PI)*dcosh(x2*PI)
         gradgradu(1,3) = 0.d0
         gradgradu(2,1) =  PI*PI*dcos(x1*PI)*dcosh(x2*PI)
         gradgradu(2,2) =  PI*PI*dsin(x1*PI)*dsinh(x2*PI)
         gradgradu(2,3) = 0.d0
         gradgradu(3,1) = 0.d0
         gradgradu(3,2) = 0.d0
         gradgradu(3,3) = 0.d0
!
!  ...smooth sin sinh z solution on unit cube
!  ...non-vanishing boundary only on x=1,y=1,z=1
      case(2)
         u = dsin(x1*PI)*dsinh(x2*PI)*(x3**2)
         gradu(1) = PI*dcos(x1*PI)*dsinh(x2*PI)*(x3**2)
         gradu(2) = PI*dsin(x1*PI)*dcosh(x2*PI)*(x3**2)
         gradu(3) = 2.*dsin(x1*PI)*dsinh(x2*PI)*x3
         gradgradu(1,1) = -PI*PI*dsin(x1*PI)*dsinh(x2*PI)*(x3**2)
         gradgradu(1,2) =  PI*PI*dcos(x1*PI)*dcosh(x2*PI)*(x3**2)
         gradgradu(1,3) =  2.*PI*dcos(x1*PI)*dsinh(x2*PI)*x3
         gradgradu(2,1) =  PI*PI*dcos(x1*PI)*dcosh(x2*PI)*(x3**2)
         gradgradu(2,2) =  PI*PI*dsin(x1*PI)*dsinh(x2*PI)*(x3**2)
         gradgradu(2,3) =  2.*PI*dsin(x1*PI)*dcosh(x2*PI)*x3
         gradgradu(3,1) =  2.*PI*dcos(x1*PI)*dsinh(x2*PI)*x3
         gradgradu(3,2) =  2.*PI*dsin(x1*PI)*dcosh(x2*PI)*x3
         gradgradu(3,3) =  2.*   dsin(x1*PI)*dsinh(x2*PI)
!
      case default
         write(*,*) 'solution: unknown exact solution (isol).'
   end select
!
end subroutine solution
