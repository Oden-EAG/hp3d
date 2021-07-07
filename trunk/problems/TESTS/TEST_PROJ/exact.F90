!-----------------------------------------------------------------------------------
!> Purpose - routine evaluates a manufactured solution
!
!> rev@Oct 14
!-----------------------------------------------------------------------------------
!  REMARKS :
!
!  1. In order to use routine "compute_error" for computing the error, for each
!     energy space variables must be listed in the following order wrt
!     attributes and rhs's :
!
!       comp1,comp2,...
!      +---------------+ +--------
!             attr1           attr2      
!            +--------------------
!                    irhs1
!
!  2. Behavior of routine is controlled by global flags NEDELEC, ICOMP
!
!    NEDELEC = 0 - lower order polynomials for all spaces; for H(curl),H(div)
!                  ICOMP = 1,2,3 selects one component (ICOMP = 0 selects all
!                  components) 
!
!    NEDELEC = 1 - higher order polynomials for H1
!
!    NEDELEC = 2 - higher order polynomials for H(curl); ICOMP = 1,2,3 selects
!                  null component
!
!    NEDELEC = 3 - higher order polynomials for H(div)
!-----------------------------------------------------------------------------------
!
#include"typedefs.h"
!
subroutine exact(Xp,Ivoid, ZvalH,ZdvalH,Zd2valH, &
                           ZvalE,ZdvalE,Zd2valE, &
                           ZvalV,ZdvalV,Zd2valV, &
                           ZvalQ,ZdvalQ,Zd2valQ)
!
      use data_structure3D , only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ,ZERO
      use PROJ             , only : NPX,NPY,NPZ,ICOMP,NEDELEC
!
      implicit none
      real*8, dimension(3),          intent(in ) :: Xp
      integer,                       intent(in ) :: Ivoid
      VTYPE,dimension(  MAXEQNH    ),intent(out) :: ZvalH
      VTYPE,dimension(  MAXEQNH,3  ),intent(out) :: ZdvalH
      VTYPE,dimension(  MAXEQNH,3,3),intent(out) :: Zd2valH
      VTYPE,dimension(3,MAXEQNE    ),intent(out) :: ZvalE
      VTYPE,dimension(3,MAXEQNE,3  ),intent(out) :: ZdvalE
      VTYPE,dimension(3,MAXEQNE,3,3),intent(out) :: Zd2valE
      VTYPE,dimension(3,MAXEQNV    ),intent(out) :: ZvalV
      VTYPE,dimension(3,MAXEQNV,3  ),intent(out) :: ZdvalV
      VTYPE,dimension(3,MAXEQNV,3,3),intent(out) :: Zd2valV
      VTYPE,dimension(  MAXEQNQ    ),intent(out) :: ZvalQ
      VTYPE,dimension(  MAXEQNQ,3  ),intent(out) :: ZdvalQ
      VTYPE,dimension(  MAXEQNQ,3,3),intent(out) :: Zd2valQ
!
      real*8 :: x,y,z,fx,fy,fz,dfx,dfy,dfz,d2fx,d2fy,d2fz,u,v,w,r
      real*8,dimension(3  ) :: du,dv,dw,dr
      real*8,dimension(3,3) :: d2u,d2v,d2w,d2r
      integer :: iprint,i,j
!
!-----------------------------------------------------------------------------------
!
      iprint=0
!      
!     cartesian coordinates
      x=Xp(1) ; y=Xp(2) ; z=Xp(3)
!      
!     initialize
      ZvalH=ZERO ; ZdvalH=ZERO ; Zd2valH=ZERO
      ZvalE=ZERO ; ZdvalE=ZERO ; Zd2valE=ZERO
      ZvalV=ZERO ; ZdvalV=ZERO ; Zd2valV=ZERO
      ZvalQ=ZERO ; ZdvalQ=ZERO ; Zd2valQ=ZERO
!
!     DEFINE MANUFACTURED SOLUTION HERE
!
!     a polynomial solution (good for all seasons!)
      select case(NPX)
      case(0) ; fx = 1.d0 ; dfx = 0.d0 ; d2fx = 0.d0
      case(1) ; fx = x    ; dfx = 1.d0 ; d2fx = 0.d0
      case default
        fx   =             x** NPX
        dfx  = NPX*        x**(NPX-1)
        d2fx = NPX*(NPX-1)*x**(NPX-2)
      endselect
!
      select case(NPY)
      case(0) ; fy = 1.d0 ; dfy = 0.d0 ; d2fy = 0.d0
      case(1) ; fy = y    ; dfy = 1.d0 ; d2fy = 0.d0
      case default
        fy   =             y** NPY
        dfy  = NPY*        y**(NPY-1)
        d2fy = NPY*(NPY-1)*y**(NPY-2)
      endselect
!
      select case(NPZ)
      case(0) ; fz = 1.d0 ; dfz = 0.d0 ; d2fz = 0.d0
      case(1) ; fz = z    ; dfz = 1.d0 ; d2fz = 0.d0
      case default
        fz   =             z** NPZ
        dfz  = NPZ*        z**(NPZ-1)
        d2fz = NPZ*(NPZ-1)*z**(NPZ-2)
      endselect
!
!     set exact solution
      u        =   fx*  fy* fz
      du( 1  ) =  dfx*  fy* fz
      du( 2  ) =   fx* dfy* fz
      du( 3  ) =   fx*  fy*dfz
      d2u(1,1) = d2fx*  fy* fz
      d2u(1,2) =  dfx* dfy* fz
      d2u(1,3) =  dfx*  fy*dfz
      d2u(2,2) =   fx*d2fy* fz
      d2u(2,3) =   fx* dfy*dfz
      d2u(3,3) =   fx*  fy*d2fz
      d2u(2,1) = d2u(1,2)
      d2u(3,1) = d2u(1,3)
      d2u(3,2) = d2u(2,3)
!
      select case(NEDELEC)
!
!-----------------------------------------------------------------------------------
!     L O W E R    O R D E R    P O L Y N O M I A L S                              |
!-----------------------------------------------------------------------------------
      case(0)
!              
!       -- H1 --
        do i=1,MAXEQNH
          ZvalH(  i        ) = u 
          ZdvalH( i,1:3    ) = du(1:3)
          Zd2valH(i,1:3,1:3) = d2u(1:3,1:3)
        enddo
!
!       components
        select case(ICOMP)
!
!       set ALL components
        case(0)
!
          do j=1,3
!          
!           -- H(curl) --
            do i=1,MAXEQNE
              ZvalE(  j,i        ) = u
              ZdvalE( j,i,1:3    ) = du(1:3)
              Zd2valE(j,i,1:3,1:3) = d2u(1:3,1:3)
            enddo
!
!           -- H(div) --
            do i=1,MAXEQNV
              ZvalV(  j,i        ) = u
              ZdvalV( j,i,1:3    ) = du(1:3)
              Zd2valV(j,i,1:3,1:3) = d2u(1:3,1:3)
            enddo
          enddo
!
!       set ONE component
        case(1,2,3)
!
!         -- H(curl) --
          do i=1,MAXEQNE
            ZvalE(  ICOMP,i        ) = u
            ZdvalE( ICOMP,i,1:3    ) = du(1:3)
            Zd2valE(ICOMP,i,1:3,1:3) = d2u(1:3,1:3)
          enddo
!
!         -- H(div) --
          do i=1,MAXEQNV
            ZvalV(  ICOMP,i        ) = u
            ZdvalV( ICOMP,i,1:3    ) = du(1:3)
            Zd2valV(ICOMP,i,1:3,1:3) = d2u(1:3,1:3)
          enddo
!
        case default
          write(*,*)'exact : invalid ICOMP! [1]'      
!          
        endselect  
!
!       -- L2 --
        do i=1,MAXEQNQ
          ZvalQ(  i        ) = u 
          ZdvalQ( i,1:3    ) = du(1:3)
          Zd2valQ(i,1:3,1:3) = d2u(1:3,1:3)
        enddo
!
!-----------------------------------------------------------------------------------
!     H I G H E R    O R D E R    P O L Y N O M I A L S                            |
!-----------------------------------------------------------------------------------
!
!     -- H1 --
      case(1)
!              
        do i=1,MAXEQNH
          ZvalH(  i        ) = u 
          ZdvalH( i,1:3    ) = du(1:3)
          Zd2valH(i,1:3,1:3) = d2u(1:3,1:3)
        enddo
!
!     -- H(curl) --
      case(2)
!
!       components
        select case(ICOMP)
!
!       null 1st component
        case(1)
          call exact_aux( z,(/0.d0, 0.d0,1.d0/),u,du,d2u, v,dv(1:3),d2v(1:3,1:3))
          call exact_aux(-y,(/0.d0,-1.d0,0.d0/),u,du,d2u, w,dw(1:3),d2w(1:3,1:3))
          do i=1,MAXEQNE
            ZvalE(  2,i        ) = v
            ZvalE(  3,i        ) = w
            ZdvalE( 2,i,1:3    ) = dv(1:3)
            ZdvalE( 3,i,1:3    ) = dw(1:3)
            Zd2valE(2,i,1:3,1:3) = d2v(1:3,1:3)
            Zd2valE(3,i,1:3,1:3) = d2w(1:3,1:3)
          enddo
!
!       null 2nd component
        case(2)
          call exact_aux(-z,(/0.d0,0.d0,-1.d0/),u,du,d2u, v,dv(1:3),d2v(1:3,1:3))
          call exact_aux( x,(/1.d0,0.d0, 0.d0/),u,du,d2u, w,dw(1:3),d2w(1:3,1:3))
          do i=1,MAXEQNE
            ZvalE(  1,i        ) = v
            ZvalE(  3,i        ) = w
            ZdvalE( 1,i,1:3    ) = dv(1:3)
            ZdvalE( 3,i,1:3    ) = dw(1:3)
            Zd2valE(1,i,1:3,1:3) = d2v(1:3,1:3)
            Zd2valE(3,i,1:3,1:3) = d2w(1:3,1:3)
          enddo
!
!       null 3rd component
        case(3)
          call exact_aux( y,(/ 0.d0,1.d0,0.d0/),u,du,d2u, v,dv(1:3),d2v(1:3,1:3))
          call exact_aux(-x,(/-1.d0,0.d0,0.d0/),u,du,d2u, w,dw(1:3),d2w(1:3,1:3))
          do i=1,MAXEQNE
            ZvalE(  1,i        ) = v
            ZvalE(  2,i        ) = w
            ZdvalE( 1,i,1:3    ) = dv(1:3)
            ZdvalE( 2,i,1:3    ) = dw(1:3)
            Zd2valE(1,i,1:3,1:3) = d2v(1:3,1:3)
            Zd2valE(2,i,1:3,1:3) = d2w(1:3,1:3)
          enddo
!
        case default
          write(*,*)'exact : invalid ICOMP! [2]'      
!          
        endselect
!
!     -- H(div) --
      case(3)
!              
        call exact_aux(x,(/1.d0,0.d0,0.d0/),u,du,d2u, v,dv(1:3),d2v(1:3,1:3))
        call exact_aux(y,(/0.d0,1.d0,0.d0/),u,du,d2u, w,dw(1:3),d2w(1:3,1:3))
        call exact_aux(z,(/0.d0,0.d0,1.d0/),u,du,d2u, r,dr(1:3),d2r(1:3,1:3))
        do i=1,MAXEQNV
          ZvalV(  1,i        ) = v
          ZvalV(  2,i        ) = w
          ZvalV(  3,i        ) = r
          ZdvalV( 1,i,1:3    ) = dv(1:3)
          ZdvalV( 2,i,1:3    ) = dw(1:3)
          ZdvalV( 3,i,1:3    ) = dr(1:3)
          Zd2valV(1,i,1:3,1:3) = d2v(1:3,1:3)
          Zd2valV(2,i,1:3,1:3) = d2w(1:3,1:3)
          Zd2valV(3,i,1:3,1:3) = d2r(1:3,1:3)
        enddo
! 
      case default
        write(*,*)'exact : invalid NEDELEC!'      
!          
      endselect
!
!
endsubroutine exact
!
!
!
! routine evaluates product of two functions and its derivatives
subroutine exact_aux(U,Ugrad,V,Vgrad,Vhesj, F,Fgrad,Fhesj)
!
      implicit none
      real*8,               intent(in ) :: U,V
      real*8,dimension(3  ),intent(in ) :: Ugrad,Vgrad
      real*8,dimension(3,3),intent(in ) :: Vhesj
      real*8,               intent(out) :: F
      real*8,dimension(3  ),intent(out) :: Fgrad
      real*8,dimension(3,3),intent(out) :: Fhesj
!      
      integer :: i,j
!      
!----------------------------------------------------------------------
!   
!     function
      F = U*V
!      
      do i=1,3
!
!       gradient
        Fgrad(i) = Ugrad(i)*V + U*Vgrad(i)
!
!       Hessian
        do j=1,3
          Fhesj(i,j) = Ugrad(i)*Vgrad(j) + Ugrad(j)*Vgrad(i) + U*Vhesj(i,j)
        enddo
      enddo
!
!
endsubroutine exact_aux
