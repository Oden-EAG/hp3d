!----------------------------------------------------------------------
!
!     routine name      - FicheraCornerDirichlet
!
!----------------------------------------------------------------------
!
!     latest revision:  -  October 2022
!
!     purpose:          - returns Boundary conditions for Fichera Corner Problem
!
!     arguments:
!        in:
!             X         - a point in physical space
!             Icase     - node case (specifies what variables are supported)
!        out:
!             ValH      - value of the H1 solution
!             DvalH     - corresponding first derivatives
!             D2valH    - corresponding second derivatives
!             DvalE     - value of the H(curl) solution
!             DdvalE    - corresponding first derivatives
!             Dd2valE   - corresponding second derivatives
!             DvalV     - value of the H(div) solution
!             DdvalV    - corresponding first derivatives
!             Dd2valV   - corresponding second derivatives
!             DvalQ     - value of the L2 solution
!             DdvalQ    - corresponding first derivatives
!             Dd2valQ   - corresponding second derivatives

subroutine FicheraCornerDirichlet(X,Icase, ValH,DvalH,D2valH, &
    ValE,DvalE,D2valE, &
    ValV,DvalV,D2valV, &
    ValQ,DvalQ,D2valQ)

    use data_structure3D
    !
       implicit none
    !
    !------------------------------------------------------------------------------
        !
        real(8), intent(in)  :: X(3)
        integer, intent(in)  :: Icase
        !
        real(8),dimension(  MAXEQNH    ), intent(out) ::   ValH
        real(8),dimension(  MAXEQNH,3  ), intent(out) ::  DvalH
        real(8),dimension(  MAXEQNH,3,3), intent(out) :: D2valH
        real(8),dimension(3,MAXEQNE    ), intent(out) ::   ValE
        real(8),dimension(3,MAXEQNE,3  ), intent(out) ::  DvalE
        real(8),dimension(3,MAXEQNE,3,3), intent(out) :: D2valE
        real(8),dimension(3,MAXEQNV    ), intent(out) ::   ValV
        real(8),dimension(3,MAXEQNV,3  ), intent(out) ::  DvalV
        real(8),dimension(3,MAXEQNV,3,3), intent(out) :: D2valV
        real(8),dimension(  MAXEQNQ    ), intent(out) ::   ValQ
        real(8),dimension(  MAXEQNQ,3  ), intent(out) ::  DvalQ
        real(8),dimension(  MAXEQNQ,3,3), intent(out) :: D2valQ

        real(8) :: x1,x2,x3,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11
        real(8) :: t12,t13,t14,t15,t16,t17,t18,t19,t20,t21
        real(8) :: u,divq
        real(8), dimension(3) :: q,gradu,zdu
        real(8),dimension(3,MAXEQNV,3) :: Dq

        integer :: isol_p, csn
        real(8) :: np_x,np_y,np_z

        real(8) :: eps,pi,a23,theta,r
        real(8) :: drdx1,drdx2,drdx3,dtdx1,dtdx2,dtdx3

    ! Initialization of arrays
        !..initialize exact solution
        ValH = ZERO ; DvalH = ZERO ; D2valH = ZERO
        ValE = ZERO ; DvalE = ZERO ; D2valE = ZERO
        ValV = ZERO ; DvalV = ZERO ; D2valV = ZERO
        ValQ = ZERO ; DvalQ = ZERO ; D2valQ = ZERO

        q = ZERO; gradu = ZERO; Dq = ZERO
    ! dirichlet boundary conditions on Hdiv trace at y = 0,2 plane, x = 0,2 plane and z = 0,2 plane
    
        x1 = X(1)
        x2 = X(2)
        x3 = X(3)

        csn = 3 !csn stands for case for solution (1 for wrong Rachowicz et al and 2 for analytical sol and 3 for correct rachowicz)

        if(csn .eq. 1) then 
        
            q(1) = (-1.d0/3.d0)*((x1**2 + x2**2)**(-7.d0/6.d0)*x1*x2 + x1*x3 *(x1**2+x3**2)**(-7.d0/6.d0))


            q(2) = (x1**2+x2**2)**(-1.d0/6.d0) - (1.d0/3.d0)*(x1**2 + x2**2)**(-7.d0/6.d0)*x2**2 &
                    - (1.d0/3.d0) * x2*x3 * (x2**2 + x3**2)**(-7.d0/6.d0)


            q(3) = -(1.d0/3.d0)*(((x2**2+x3**2)**(-7.d0/6.d0) + (x1**2+x3**2)**(-7.d0/6.d0))*x3**2) &
                    + (x1**2 + x3**2)**(-1.d0/6.d0) + (x2**2 + x3**2)**(-1.d0/6.d0)


            t2 = x1**2
            t3 = x2**2
            t4 = x3**2
            t5 = t2+t3
            t6 = t2+t4

            Dq(1,1,1) = 1.0d0/t5**(7.0d0/6.0d0)*x2*(-1.0d0/3.0d0)-(1.0d0/t6**(7.0d0/6.0d0 &
            )*x3)/3.0d0+t2*1.0d0/t5**(1.3D+1/6.0d0)*x2*(7.0d0/9.0d0)+t2*1.0d0 &
            /t6**(1.3D+1/6.0d0)*x3*(7.0d0/9.0d0)


            t2 = x1**2
            t3 = x2**2
            t4 = t2+t3
            Dq(1,1,2) = 1.0d0/t4**(7.0d0/6.0d0)*x1*(-1.0d0/3.0d0)+t3*1.0d0/t4**(1.3d+1 &
            /6.0d0)*x1*(7.0d0/9.0d0)
      

            t2 = x1**2
            t3 = x3**2
            t4 = t2+t3
            Dq(1,1,3) = 1.0d0/t4**(7.0d0/6.0d0)*x1*(-1.0d0/3.0d0)+t3*1.0d0/t4**(1.3d+1 &
            /6.0d0)*x1*(7.0d0/9.0D0)
      

  
            Dq(2,1,1) = Dq(1,1,2)


            t2 = x1**2
            t3 = x2**2
            t4 = x3**2
            t5 = t2+t3
            t6 = t3+t4
            Dq(2,1,2) = 1.0d0/t5**(1.3d+1/6.0d0)*x2**3*(7.0d0/9.0d0)-1.0d0/t5**(7.0d0 &
            /6.0d0)*x2-(1.0d0/t6**(7.0d0/6.0d0)*x3)/3.0d0+t3*1.0d0/t6**(1.3d+1  &
            /6.0d0)*x3*(7.0d0/9.0d0)

            t2 = x2**2
            t3 = x3**2
            t4 = t2+t3
            Dq(2,1,3) = 1.0d0/t4**(7.0d0/6.0d0)*x2*(-1.0d0/3.0d0)+t3*1.0d0/t4**(1.3d+1 &
            /6.0d0)*x2*(7.0d0/9.0d0)


            Dq(3,1,1) = Dq(1,1,3)
                


            Dq(3,1,2) = Dq(2,1,3)
                
                
            t2 = x1**2
            t3 = x2**2
            t4 = x3**2
            t5 = x3**3
            t6 = t2+t4
            t7 = t3+t4
            Dq(3,1,3) = t5*1.0d0/t6**(1.3d+1/6.0d0)*(7.0d0/9.0d0)+t5*1.0d0/t7**(1.3d+1 &
            /6.0d0)*(7.0d0/9.0d0)-1.0d0/t6**(7.0d0/6.0d0)*x3-1.0d0/t7**(7.0d0 &
            /6.0d0)*x3
      
            !  u = x2 * (x1**2 + x2**2)**(-1.d0/6.d0) + x3 * (x2**2 + x1**2)**(-1.d0/6.d0) + x3 * (x1**2 + x3**2)**(-1.d0/6.d0)
            !  gradu(1:3) = q 

            u = 0.d0
            gradu(1:3) = ZERO 

            ValV(1:3,1) = q(1:3)  !Hdiv values
            DvalV(1:3,1,1:3) = Dq(1:3,1,1:3) !change it 10^70 so codes breaks
            ValH = ZERO !H1 values
            DvalH(1,1:3) = ZERO !derivative of H1 values
            D2valH(1,1:3,1:3) = ZERO
        
        else if(csn .eq. 2) then



                call exact(X,Icase, ValH,DvalH,D2valH,ValE,DvalE,D2valE, &
                ValV,DvalV,D2valV,ValQ,DvalQ,D2valQ)

        else if(csn .eq. 3) then

            eps = epsilon(1.d0)
            zdu = ZERO
            pi  = atan(1.d0)*4.d0
            a23 = 2.d0/3.d0

            if (r<eps) goto 101
            
            if(x3.ge.abs(x1) .and. x1.le.0.d0)then
                theta = pi/2.d0 + atan(-x1/x3)
            elseif(x1.le.-abs(x3))then
                theta = pi + atan(x3/x1)
            elseif(x3.le.-abs(x1))then
                theta = 1.5d0*pi + atan(-x1/x3)
            elseif(x1.ge.abs(x3) .and. x3.le.0)then
                theta = 2.d0*pi + atan(x3/x1)
            elseif(x1.gt.0.d0 .and. x3 .gt. 0.d0) then
                theta = atan(x3/x1)
            else
                goto 101
            endif

            drdx1 =  x1/r
            drdx3 =  x3/r
            dtdx1 = -x3/r**2
            dtdx3 =  x1/r**2
      
            zdu(1) = zdu(1) + &
                     a23 * r**(a23-1) * drdx1 * sin( (theta-pi/2.d0)*a23 ) &
                     + r**a23 * cos( (theta-pi/2.d0)*a23 ) * a23 * dtdx1

            zdu(3) = zdu(3) + &
                     a23 * r**(a23-1) * drdx3 * sin( (theta-pi/2.d0)*a23 ) &
                     + r**a23 * cos( (theta-pi/2.d0)*a23 ) * a23 * dtdx3 

            101  continue

            x2 = X(2); x3 = X(3); r = sqrt(x2**2+x3**2)
            if (r<eps) goto 102

            if(x3.ge.abs(x2) .and. x2.le.0.d0)then
              theta = pi/2.d0 + atan(-x2/x3)
            elseif(x2.le.-abs(x3))then
              theta = pi + atan(x3/x2)
            elseif(x3.le.-abs(x2))then
              theta = 1.5d0*pi + atan(-x2/x3)
            elseif(x2.ge.abs(x3) .and. x3.le.0)then
              theta = 2.d0*pi + atan(x3/x2)
            elseif(x2.gt.0.d0 .and. x3 .gt. 0.d0) then
                theta = atan(x3/x2)
            else
               goto 102
            endif

            drdx2 =  x2/r
            drdx3 =  x3/r
            dtdx2 = -x3/r**2
            dtdx3 =  x2/r**2
      
            zdu(2) = zdu(2) +  &
                     a23 * r**(a23-1) * drdx2 * sin( (theta-pi/2.d0)*a23 ) &
                   + r**a23 * cos( (theta-pi/2.d0)*a23 ) * a23 * dtdx2  

            zdu(3) = zdu(3) + &
                     a23 * r**(a23-1) * drdx3 * sin( (theta-pi/2.d0)*a23 ) &
                    + r**a23 * cos( (theta-pi/2.d0)*a23 ) * a23 * dtdx3  

            102  continue


            x1 = X(1);x2 = X(2); r = sqrt(x2**2+x1**2)
            if (r < eps) goto 103

            if (x1<=0.d0 .and. x2>=0.d0) then
                theta = atan2(x2,x1) - pi/2.d0
            elseif (x2<0.d0) then
                theta = atan2(x2,x1) + 3.d0*pi/2.d0
            elseif (x2==0.d0) then
                theta = 3.d0*pi/2.d0
            else
                goto 103
            endif

              drdx2 =  x2/r
              drdx1 =  x1/r
              dtdx2 =  x1/r**2
              dtdx1 = -x2/r**2
        
              zdu(2) = zdu(2) + &
                       a23 * r**(a23-1) * drdx2 * sin( theta*a23 )  &
                     + r**a23 * cos( theta*a23 ) * a23 * dtdx2 
              zdu(1) = zdu(1) + &
                       a23 * r**(a23-1) * drdx1 * sin( theta*a23 )  &
                     + r**a23 * cos( theta*a23 ) * a23 * dtdx1 

                !   if (x2 .ge. abs(x1) .and. x1 .le. 0.d0) then
                !     theta = pi/2.d0 + atan(-x1/x2)
                !   elseif (x1 .le. -abs(x2)) then
                !     theta = pi + atan(x2/x1)
                !   elseif (x2 .le. -abs(x1)) then
                !     theta = 1.5d0*pi + atan(-x1/x2)
                !   elseif (x1 .ge. abs(x2) .and. x2 .le. 0.d0)then
                !     theta = 2.d0*pi + atan(x2/x1)
                !   elseif(x1.gt.0.d0 .and. x2 .gt. 0.d0) then
                !     theta = atan(x2/x1)
                !   else
                !      goto 103
                !   endif

                !   drdx2 = x2/r
                !   drdx1 = x1/r
                !   dtdx2 =  x1/r**2
                !   dtdx1 = -x2/r**2
                !   zdu(2) = zdu(2) + &
                !            a23 * r**(a23-1) * drdx2 * sin( (theta-pi/2.d0)*a23 )  &
                !          + r**a23 * cos( (theta-pi/2.d0)*a23 ) * a23 * dtdx2

                !   zdu(1) = zdu(1) + &
                !            a23 * r**(a23-1) * drdx1 * sin( (theta-pi/2.d0)*a23 )  &
                !          + r**a23 * cos( (theta-pi/2.d0)*a23 ) * a23 * dtdx1
        
         103  continue

            ValV(1:3,1) = zdu(1:3)  !Hdiv values
            DvalV(1:3,1,1:3) = 0 !change it 10^70 so codes breaks
            ValH = ZERO !H1 values
            DvalH(1,1:3) = ZERO !derivative of H1 values
            D2valH(1,1:3,1:3) = ZERO

        endif

    end subroutine FicheraCornerDirichlet