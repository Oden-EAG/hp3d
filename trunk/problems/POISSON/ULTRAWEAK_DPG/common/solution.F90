!-----------------------------------------------------------------------------------
!
!     routine name      - solution
!
!-----------------------------------------------------------------------------------
!
!     latest revision:  - May 2023
!
!> @brief         - compute all relevant quantities of the exact solutions
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
   use common_prob_data, only : PI,ISOL
   use parameters      , only : ZERO
!
   implicit none
!
   real(8), dimension(3),   intent(in)  :: X
   real(8),                 intent(out) :: u
   real(8), dimension(3),   intent(out) :: gradu     ! 1st derivative - gradient
   real(8), dimension(3,3), intent(out) :: gradgradu ! 2nd derivative - Hessian
!
   real(8) :: x1,x2,x3,f_x,f_y,f_z,df_x,df_y,df_z,ddf_x,ddf_y,ddf_z,eps
   real(8) :: np_x,np_y,np_z
   integer :: isol_p
   real(8) :: x1c,x2c,x3c,alpha,ro
   real(8) :: t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20
   real(8) :: t21,t22,t23,t24,t25,t26
   real(8) :: u1,u2,u3,u1x,u2x,u3x,u1xx,u2xx,u3xx
!
!--------------------------------------------------------------------------------
!
   x1 = X(1); x2 = X(2); x3 = X(3)
!
   u = ZERO; gradu = ZERO; gradgradu = ZERO
!
!..select exact solution
   select case (ISOL)
!
!..polynomial solution
      case(0)
!..set polynomial order of manufactured solution
         isol_p = 3
         np_x = real(isol_p,8)
         np_y = real(isol_p,8)
         np_z = real(isol_p,8)
!..value
         f_x = x1**np_x
         f_y = x2**np_y
         f_z = x3**np_z
!..derivatives
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
!..value
         u = f_x*f_y*f_z
!..1st order derivatives
         gradu(1)=  df_x *   f_y *   f_z
         gradu(2)=   f_x *  df_y *   f_z
         gradu(3)=   f_x *   f_y *  df_z
!..2nd order derivatives
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
!..smooth sin sinh z solution on unit cube
!..non-vanishing boundary only on x=1,y=1,z=1
!..uniform in z
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
!..smooth sin sinh z solution on unit cube
!..non-vanishing boundary only on x=1,y=1,z=1
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
!..arc tan shock like solution
      case(3)
         x1c = -2.5d-1
         x2c = -2.5d-1
         x3c = -2.5d-1
         ro = 3.d0**(5.d-1)
         alpha = 20.d0
         u = -atan(alpha*(ro-sqrt((x1-x1c)**2+(x2-x2c)**2+(x3-x3c)**2)))
         t2 = -x1c
         t3 = -x2c
         t4 = -x3c
         t5 = t2+x1
         t6 = t3+x2
         t7 = t4+x3
         t8 = t5**2
         t9 = t6**2
         t10 = t7**2
         t11 = t8+t9+t10
         gradu(1) = (alpha*1.0D0/sqrt(t11)*(x1*2.0D0-x1c*2.0D0))/(alpha**2*(ro-sqrt(t11))**2*2.0D0+2.0D0)
!
         t2 = -x1c
         t3 = -x2c
         t4 = -x3c
         t5 = t2+x1
         t6 = t3+x2
         t7 = t4+x3
         t8 = t5**2
         t9 = t6**2
         t10 = t7**2
         t11 = t8+t9+t10
         gradu(2) = (alpha*1.0D0/sqrt(t11)*(x2*2.0D0-x2c*2.0D0))/(alpha**2*(ro-sqrt(t11))**2*2.0D0+2.0D0)
!
         t2 = -x1c
         t3 = -x2c
         t4 = -x3c
         t5 = t2+x1
         t6 = t3+x2
         t7 = t4+x3
         t8 = t5**2
         t9 = t6**2
         t10 = t7**2
         t11 = t8+t9+t10
         gradu(3) = (alpha*1.0D0/sqrt(t11)*(x3*2.0D0-x3c*2.0D0))/(alpha**2*(ro-sqrt(t11))**2*2.0D0+2.0D0)
!
         t2 = alpha**2
         t3 = -x1c
         t4 = -x2c
         t5 = -x3c
         t6 = t3+x1
         t7 = t4+x2
         t8 = t5+x3
         t9 = t6**2
         t10 = t7**2
         t11 = t8**2
         t12 = t9+t10+t11
         t13 = sqrt(t12)
         t14 = -t13
         t15 = ro+t14
         t16 = t15**2
         t17 = t2*t16
         t18 = t17+1.0D0
         t19 = 1.0D0/t18
         gradgradu(1,1) = (alpha*t19)/t13-alpha*t9*1.0D0/t13**3*t19+alpha**3*t9*1.0D0/t13**2*t15*t19**2*2.0D0
!
         t2 = alpha**2
         t3 = x1*2.0D0
         t4 = x2*2.0D0
         t5 = x1c*2.0D0
         t6 = x2c*2.0D0
         t7 = -x1c
         t9 = -x2c
         t11 = -x3c
         t8 = -t5
         t10 = -t6
         t12 = t7+x1
         t13 = t9+x2
         t14 = t11+x3
         t15 = t3+t8
         t16 = t4+t10
         t17 = t12**2
         t18 = t13**2
         t19 = t14**2
         t20 = t17+t18+t19
         t21 = sqrt(t20)
         t22 = -t21
         t23 = ro+t22
         t24 = t23**2
         t25 = t2*t24
         t26 = t25+1.0D0
         gradgradu(1,2) = (alpha*t15*t16*1.0D0/t21**3*(-1.0D0/4.0D0))/t26+(alpha**3*t15*t16*1.0D0/t21**2*t23*1.0D0/t26**2)/2.0D0
!
         t2 = alpha**2
         t3 = x1*2.0D0
         t4 = x3*2.0D0
         t5 = x1c*2.0D0
         t6 = x3c*2.0D0
         t7 = -x1c
         t9 = -x2c
         t10 = -x3c
         t8 = -t5
         t11 = -t6
         t12 = t7+x1
         t13 = t9+x2
         t14 = t10+x3
         t15 = t3+t8
         t16 = t4+t11
         t17 = t12**2
         t18 = t13**2
         t19 = t14**2
         t20 = t17+t18+t19
         t21 = sqrt(t20)
         t22 = -t21
         t23 = ro+t22
         t24 = t23**2
         t25 = t2*t24
         t26 = t25+1.0D0
         gradgradu(1,3) = (alpha*t15*t16*1.0D0/t21**3*(-1.0D0/4.0D0))/t26+(alpha**3*t15*t16*1.0D0/t21**2*t23*1.0D0/t26**2)/2.0D0
!
         gradgradu(2,1) = gradgradu(1,2)
!
         t2 = alpha**2
         t3 = -x1c
         t4 = -x2c
         t5 = -x3c
         t6 = t3+x1
         t7 = t4+x2
         t8 = t5+x3
         t9 = t6**2
         t10 = t7**2
         t11 = t8**2
         t12 = t9+t10+t11
         t13 = sqrt(t12)
         t14 = -t13
         t15 = ro+t14
         t16 = t15**2
         t17 = t2*t16
         t18 = t17+1.0D0
         t19 = 1.0D0/t18
         gradgradu(2,2) = (alpha*t19)/t13-alpha*t10*1.0D0/t13**3*t19+alpha**3*t10*1.0D0/t13**2*t15*t19**2*2.0D0
! 
         t2 = alpha**2
         t3 = x2*2.0D0
         t4 = x3*2.0D0
         t5 = x2c*2.0D0
         t6 = x3c*2.0D0
         t7 = -x1c
         t8 = -x2c
         t10 = -x3c
         t9 = -t5
         t11 = -t6
         t12 = t7+x1
         t13 = t8+x2
         t14 = t10+x3
         t15 = t3+t9
         t16 = t4+t11
         t17 = t12**2
         t18 = t13**2
         t19 = t14**2
         t20 = t17+t18+t19
         t21 = sqrt(t20)
         t22 = -t21
         t23 = ro+t22
         t24 = t23**2
         t25 = t2*t24
         t26 = t25+1.0D0
         gradgradu(2,3) = (alpha*t15*t16*1.0D0/t21**3*(-1.0D0/4.0D0))/t26+(alpha**3*t15*t16*1.0D0/t21**2*t23*1.0D0/t26**2)/2.0D0
!
         gradgradu(3,1) = gradgradu(1,3)
         gradgradu(3,2) = gradgradu(2,3)
!
         t2 = alpha**2
         t3 = -x1c
         t4 = -x2c
         t5 = -x3c
         t6 = t3+x1
         t7 = t4+x2
         t8 = t5+x3
         t9 = t6**2
         t10 = t7**2
         t11 = t8**2
         t12 = t9+t10+t11
         t13 = sqrt(t12)
         t14 = -t13
         t15 = ro+t14
         t16 = t15**2
         t17 = t2*t16
         t18 = t17+1.0D0
         t19 = 1.0D0/t18
         gradgradu(3,3) = (alpha*t19)/t13-alpha*t11*1.0D0/t13**3*t19+alpha**3*t11*1.0D0/t13**2*t15*t19**2*2.0D0
!
!..single boundary layer
      case(4)
!
         eps = 5.0d-3
         u = x1 + (exp(x1/eps) - 1.d0)/(1.d0 - exp(1.d0/eps))
!
         gradu = ZERO
         gradu(1) = 1.d0 + (1.d0/eps) * (exp(x1/eps)/(1.d0 - exp(1.d0/eps)))
!
         gradgradu = ZERO
!
         gradgradu(1,1) = (1.d0/eps**2) * (exp(x1/eps)/(1.d0 - exp(1.d0/eps)))
!
!..triple  boundary layer
      case(5)
!
         eps = 5.d-3
         u1 = x1 + (exp(x1/eps) - 1.d0)/(1.d0 - exp(1.d0/eps))
         u2 = x2 + (exp(x2/eps) - 1.d0)/(1.d0 - exp(1.d0/eps))
         u3 = x3 + (exp(x3/eps) - 1.d0)/(1.d0 - exp(1.d0/eps))
!
         u = u1 * u2 * u3
!
         u1x = 1.d0 + (1.d0/eps) * (exp(x1/eps)/(1.d0 - exp(1.d0/eps)))
         u2x = 1.d0 + (1.d0/eps) * (exp(x2/eps)/(1.d0 - exp(1.d0/eps)))
         u3x = 1.d0 + (1.d0/eps) * (exp(x3/eps)/(1.d0 - exp(1.d0/eps)))
!
         u1xx = (1.d0/eps**2) * (exp(x1/eps)/(1.d0 - exp(1.d0/eps)))   
         u2xx = (1.d0/eps**2) * (exp(x2/eps)/(1.d0 - exp(1.d0/eps)))   
         u3xx = (1.d0/eps**2) * (exp(x3/eps)/(1.d0 - exp(1.d0/eps)))   
!
         gradu = ZERO
         gradu(1) = (1.d0 + (1.d0/eps) * (exp(x1/eps)/(1.d0 - exp(1.d0/eps)))) * u2 * u3
         gradu(2) = u1 * (1.d0 + (1.d0/eps) * (exp(x2/eps)/(1.d0 - exp(1.d0/eps)))) * u3
         gradu(3) = u1 * u2 * (1.d0 + (1.d0/eps) * (exp(x3/eps)/(1.d0 - exp(1.d0/eps)))) 
!
         gradgradu = ZERO
!
         gradgradu(1,1) = u1xx * u2 * u3
         gradgradu(1,2) = u1x * u2x * u3
         gradgradu(1,3) = u1x * u2 * u3x
         gradgradu(2,1) = gradgradu(1,2)
         gradgradu(2,2) = u1 * u2xx * u3
         gradgradu(2,3) = u1 * u2x * u3x
         gradgradu(3,1) = gradgradu(1,3)
         gradgradu(3,2) = gradgradu(2,3)
         gradgradu(3,3) = u1 * u2 * u3xx
!
      case default
         write(*,*) 'solution: unknown exact solution (ISOL).'
   end select
!
end subroutine solution
