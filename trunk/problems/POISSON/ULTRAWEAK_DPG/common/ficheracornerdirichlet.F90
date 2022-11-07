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
        real(8), dimension(3) :: q,gradu
        real(8),dimension(3,MAXEQNV,3) :: Dq

        integer :: isol_p
        real(8) :: np_x,np_y,np_z

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
      


        

        ! else if(x1 .eq. 0.d0) then
        
            u = 0.d0
            gradu(1:3) = ZERO 

        ! else if (x2 .eq. 0.d0) then
        !     u = 0.d0
        !     gradu(1:3) = ZERO

        ! else if (x3 .eq. 0.d0) then
        !     u = 0.d0
        !     gradu(1:3) = ZERO
        
        ! endif

        

        ! if((x1 .eq. 0.d0) .or. (x1 .eq. 2.d0) .or. (x2 .eq. 0.d0) .or. (x2 .eq. 2.d0) &
        !     .or. (x3 .eq. 0.d0) .or. (x3 .eq. 2.d0)) then

        !    isol_p = 3
        ! if((x3 .eq. 1.d0)) then

                ! np_x = real(isol_p,8)
                ! np_y = real(isol_p,8)
                ! np_z = real(isol_p,8)

                ! q(1) = np_x * x1**(np_x-1.d0) * x2**np_y * x3**np_z
                ! q(2) = np_y * x2**(np_y-1.d0) * x1**np_x * x3**np_z
                ! q(3) = np_z * x3**(np_z-1.d0) * x1**np_x * x2**np_y


                ! Dq(1,1,1) =   np_x * (np_x - 1.d0) * x1**(np_x - 2.d0) * x2**np_y * x3**np_z
                ! Dq(1,1,2) =   np_x * x1**(np_x-1.d0) * np_y * x2**(np_y-1.d0) *  x3**np_z
                ! Dq(1,1,3) =   np_x * x1**(np_x-1.d0) * x2**np_y * np_z * x3**(np_z - 1.d0)

                ! Dq(2,1,1) =   Dq(1,1,1)
                ! Dq(2,1,2) =   np_y * (np_y-1.d0) * x2**(np_y-2.d0) * x1**np_x * x3**np_z
                ! Dq(2,1,3) =   np_y * x2**(np_y-1.d0) * x1**np_x * np_z * x3**(np_z - 1.d0) 

                ! Dq(3,1,1) =   Dq(1,1,3)
                ! Dq(3,1,2) =   Dq(2,1,3)
                ! Dq(3,1,3) =   np_z * (np_z - 1.d0) * x3**(np_z-2.d0) * x1**np_x * x2**np_y

    
        ! ! else
                ! np_x = real(isol_p,8)
                ! np_y = real(isol_p,8)
                ! np_z = real(isol_p,8)

                ! u = x1**np_x * x2**np_y * x3**np_z

                ! gradu(1) = np_x * x1**(np_x-1.d0) * x2**np_y * x3**np_z
                ! gradu(2) = np_y * x2**(np_y-1.d0) * x1**np_x * x3**np_z
                ! gradu(3) = np_z * x3**(np_z-1.d0) * x1**np_x * x2**np_y
            
        ! endif

        ValV(1:3,1) = q(1:3)  !Hdiv values
        DvalV(1:3,1,1:3) = Dq(1:3,1,1:3)
        ValH(1) = u !H1 values
        DvalH(1,1:3) = gradu !derivative of H1 values
    end subroutine FicheraCornerDirichlet