!------------------------------------------------------------------------------
!
!     routine name      - FicheraCornerDirichlet
!
!------------------------------------------------------------------------------
!
!     latest revision:  - May 2023
!
!> @brief         - returns Boundary conditions for Fichera Corner Problem
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
subroutine ficheracornerdirichlet(X,Icase, ValH,DvalH,D2valH, &
                                 ValE,DvalE,D2valE, &
                                 ValV,DvalV,D2valV, &
                                 ValQ,DvalQ,D2valQ)
!
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
!
    real(8) :: x1,x2,x3
    real(8) :: zdu(3)
    real(8) :: gradu(3)
    real(8) :: q(3)
    real(8) :: Dq(3,MAXEQNV,3)
!
    integer :: csn
!
    real(8) :: eps,pi,a23,theta,r
    real(8) :: drdx1,drdx2,drdx3,dtdx1,dtdx2,dtdx3
!
!..Initialization of arrays
!
    ValH = ZERO ; DvalH = ZERO ; D2valH = ZERO
    ValE = ZERO ; DvalE = ZERO ; D2valE = ZERO
    ValV = ZERO ; DvalV = ZERO ; D2valV = ZERO
    ValQ = ZERO ; DvalQ = ZERO ; D2valQ = ZERO

    q = ZERO; gradu = ZERO; Dq = ZERO
!..dirichlet boundary conditions on Hdiv trace at
!  y = 0,2 plane, x = 0,2 plane and z = 0,2 plane
!
    x1 = X(1)
    x2 = X(2)
    x3 = X(3)
!    
!..csn stands for case for solution (1 for analytical sol and 2 for Fichera Corner)
    csn = 1
!
    if(csn .eq. 1) then 
!        
        call exact(X,Icase, ValH,DvalH,D2valH,ValE,DvalE,D2valE, &
                   ValV,DvalV,D2valV,ValQ,DvalQ,D2valQ)
!
    else if(csn .eq. 2) then
!
        eps = epsilon(1.d0)
        zdu = ZERO
        pi  = atan(1.d0)*4.d0
        a23 = 2.d0/3.d0
!
        if (r<eps) goto 101
!        
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
!
        drdx1 =  x1/r
        drdx3 =  x3/r
        dtdx1 = -x3/r**2
        dtdx3 =  x1/r**2
!   
        zdu(1) = zdu(1) + &
                 a23 * r**(a23-1) * drdx1 * sin( (theta-pi/2.d0)*a23 ) &
                 +r**a23 * cos( (theta-pi/2.d0)*a23 ) * a23 * dtdx1
!
        zdu(3) = zdu(3) + &
                 a23 * r**(a23-1) * drdx3 * sin( (theta-pi/2.d0)*a23 ) &
                 +r**a23 * cos( (theta-pi/2.d0)*a23 ) * a23 * dtdx3 
!
        101  continue
!
        x2 = X(2); x3 = X(3); r = sqrt(x2**2+x3**2)
        if (r<eps) goto 102
!
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
!
        drdx2 =  x2/r
        drdx3 =  x3/r
        dtdx2 = -x3/r**2
        dtdx3 =  x2/r**2
!   
        zdu(2) = zdu(2) +  &
                    a23 * r**(a23-1) * drdx2 * sin( (theta-pi/2.d0)*a23 ) &
                + r**a23 * cos( (theta-pi/2.d0)*a23 ) * a23 * dtdx2  

        zdu(3) = zdu(3) + &
                    a23 * r**(a23-1) * drdx3 * sin( (theta-pi/2.d0)*a23 ) &
                + r**a23 * cos( (theta-pi/2.d0)*a23 ) * a23 * dtdx3  
!
        102  continue
!
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
!
        drdx2 =  x2/r
        drdx1 =  x1/r
        dtdx2 =  x1/r**2
        dtdx1 = -x2/r**2
!
        zdu(2) = zdu(2) + &
                a23 * r**(a23-1) * drdx2 * sin( theta*a23 )  &
                + r**a23 * cos( theta*a23 ) * a23 * dtdx2 
        zdu(1) = zdu(1) + &
                a23 * r**(a23-1) * drdx1 * sin( theta*a23 )  &
                + r**a23 * cos( theta*a23 ) * a23 * dtdx1 
!   
        103  continue
!
        ValV(1:3,1) = zdu(1:3)  !Hdiv values
        DvalV(1:3,1,1:3) = 0 !change it 10^70 so codes breaks
        ValH = ZERO !H1 values
        DvalH(1,1:3) = ZERO !derivative of H1 values
        D2valH(1,1:3,1:3) = ZERO
!
    endif
!
end subroutine ficheracornerdirichlet
