c----------------------------------------------------------------------
c
c   routine name       - test_rJacobi
c
c----------------------------------------------------------------------
c
      subroutine test_rJacobi
c
#include "syscom.blk"
c      integer :: combination
c
      iprint = 1
c
      eps = 1.d-10
      min_a = 1000
      min_n = 1000
c
      max_a = 10
      max_n = 10
      max_i = 20
c
      write(*,*)'TESTING WEIGHTED JACOBI POLYNOMIALS...'
c
      do ia=0,max_a
        do n=0,max_n
          if (iprint.eq.1) then
            write(*,7002)ia,n
 7002       format(' a,n = ',2(i2,2x))
          endif
          a = real(ia)
c
c  .......values at -1
          r = (-1.d0)**n
          err = abs(rJacobi_W(a,0.d0,-1.d0,n) - r)
          if (err .ge. eps) then
ccc            min_a = min(min_a,ia); min_n = min(min_n,n)
            write(*,7001)ia,n,rJacobi_W(a,0.d0,-1.d0,n),r,err
 7001       format(' a,n = ',2(i2,2x),'; p(-1),r,err = ',3(e12.5,2x))
ccc            call pause
          endif
c
cccc  .......values at 1
ccc          r = combination(int(n+a),n)
ccc          err = abs(rJacobi_W(a,0.d0, 1.d0,n) - r)
ccc          if (err .ge. eps) then
cccccc            min_a = min(min_a,ia); min_n = min(min_n,n)
ccc            write(*,7003)ia,n,rJacobi_W(a,0.d0,1.d0,n),r,err
ccc 7003       format(' a,n = ',2(i2,2x)'; p( 1),r,err = ',3(e12.5,2x))
ccc            call pause
ccc          endif
cccc
cccc  .......derivatives at 1
ccc          r = n*combination(int(n+a+1),n)/2.d0
ccc          err = abs(drJacobi_W(a,0.d0,1.d0,n) - r)
ccc          if (err .ge. eps) then
ccc            min_a = min(min_a,ia); min_n = min(min_n,n)
ccc            write(*,7004)ia,n,drJacobi_W(a,0.d0,1.d0,n),r,err
ccc 7004       format(' a,n = ',2(i2,2x)'; dp(1),r,err = ',3(e12.5,2x))
ccc            call pause
ccc          endif
c
c  .......values in [-1,1]
          if ((ia.ge.1).and.(n.ge.1)) then
            do i=0,max_i
              x = -1.d0 + 2.d0*i/max_i
c  ...........1st test: eq (2.5)
              r = ((a+n)*rJacobi_W(a,0.d0,x,n) -
     .                 n*rJacobi_W(a,0.d0,x,n-1))/(a+2.d0*n)
              err = abs(rJacobi_W(a-1,0.d0,x,n) - r)
              if (err .ge. eps) then
                min_a = min(min_a,ia); min_n = min(min_n,n)
                write(*,7006)ia,n,x,rJacobi_W(a-1,0.d0,x,n),r,err
 7006           format(' I:  a,n = ',2(i2,2x),'; x = ',e12.5,
     .                 '; p(x),r,err = ',3(e12.5,2x))
ccc                call pause
              endif
c  ...........2nd test: eq (2.6)
              r = (2.d0*n+a+1.d0)*((2.d0*n+a+2.d0)*(2.d0*n+a)*x+a**2)/
     .            ((2.d0*n+2.d0)*(n+a+1.d0)*(2.d0*n+a))*
     .            rJacobi_W(a,0.d0,x,n) -
     .            n*(n+a)*(2.d0*n+a+2.d0)/
     .                    ((n+1.d0)*(n+a+1.d0)*(2.d0*n+a))*
     .            rJacobi_W(a,0.d0,x,n-1)
              err = abs(rJacobi_W(a,0.d0,x,n+1) - r)
              if (err .ge. eps) then
                min_a = min(min_a,ia); min_n = min(min_n,n)
                write(*,7007)ia,n,x,rJacobi_W(a,0.d0,x,n+1),r,err
 7007           format(' II: a,n = ',2(i2,2x),'; x = ',e12.5,
     .                 '; p(x),r,err = ',3(e12.5,2x))
ccc                call pause
              endif
            enddo
          endif
        enddo
      enddo
c
      write(*,*)'TESTING WEIGHTED INTEGRATED JACOBI POLYNOMIALS...'
      do ia=1,max_a
        do n=1,max_n
          a = real(ia)
          write(*,7002)ia,n
c
c  .......values at -1
          r = 0.d0
          err = abs(rJacobi_W_I(a,-1.d0,n) - r)
          if (err .ge. eps) then
            min_a = min(min_a,ia); min_n = min(min_n,n)
            write(*,7008)ia,n,rJacobi_W_I(a,-1.d0,n),r,err
 7008       format(' a,n = ',2(i2,2x),'; \hat{p}(-1),r,err = ',
     .              3(e12.5,2x))
ccc            call pause
          endif
c
c  .......values in [-1,1]: eq (2.11)
          do i=0,max_1
            x = -1.d0 + 2.d0*i/max_i
            r=(2.d0*n+a-1.d0)*((2.d0*n+a-2.d0)*(2.d0*n+a)*x+a*(a-2.d0))/
     .          ((2.d0*n+2.d0)*(n+a)*(2.d0*n+a-2.d0))*
     .          rJacobi_W_I(a,x,n) -
     .          (n-1.d0)*(n+a-2.d0)*(2.d0*n+a)/
     .                  ((n+1.d0)*(n+a)*(2.d0*n+a-2.d0))*
     .          rJacobi_W_I(a,x,n-1)
            err = abs(rJacobi_W_I(a,x,n+1) - r)
            if (err .ge. eps) then
              min_a = min(min_a,ia); min_n = min(min_n,n)
              write(*,7009)ia,n,x,rJacobi_W_I(a,x,n+1),r,err
 7009         format(' a,n = ',2(i2,2x),'; x = ',e12.5,
     .               '; \hat{p}(x),r,err = ',3(e12.5,2x))
ccc              call pause
            endif
c
          enddo
        enddo
      enddo
c
      write(*,7010)min_a,min_n
 7010 format('min_a,min_n = ',2(i2))
      call pause

      end
c
c
c
c----------------------------------------------------------------------
c
c   routine name       - build_Jacobi_W_I
c
c----------------------------------------------------------------------
c
c   latest revision    - Oct 10
c
c   purpose            - build block of Weighted Integrated Jacobi
c                        polynomials and their derivatives
c
c   arguments :
c     in:
c                A     - weights
c                X     - variable
c                Nord  - order of polynomial
c     out:
c                Poly  - polynomials
c                Dpoly - derivatives
c
c----------------------------------------------------------------------
c
      subroutine build_Jacobi_W_I(A,X,Nord, Poly,Dpoly)
c
#include "syscom.blk"
c
c  ...set up maximum degree supported
      integer, parameter :: nmax = 20
c  ...output blocks
      dimension Poly(1:nmax),Dpoly(1:nmax)
c
c
      if (Nord.gt.nmax) then
        write(*,*)'build_Jacobi_W_I: Nord = ',Nord
        stop
      endif
c
c  ...build blocks
      do n=1,Nord
        Poly( n) =  rJacobi_W_I(A,X,n)
        Dpoly(n) = drJacobi_W_I(A,X,n)
      enddo
c
      end subroutine
c
c
c
c----------------------------------------------------------------------
c
c   routine name       - rJacobi_W
c
c----------------------------------------------------------------------
c
c   latest revision    - Oct 10
c
c   purpose            - evaluates Weighted Jacobi polynomial
c
c   arguments :
c     in:
c                A,B   - weights
c                X     - variable
c                Nord  - order of polynomial
c
c----------------------------------------------------------------------
c
c   2-terms recurrence relation from Wolfram
c
c   c1 P_{n+1}^(a,b) = (c2 + c3*x) P_n^(a,b) + c3 P_{n-1}^(a,b)
c          P_0^(a,b) = 1
c          P_1^(a,b) = [2(a + 1) + (a + b + 2) (x - 1)]/2
c
c   c1 = 2(n + 1) (n + a + b + 1) (2n + a + b)
c   c2 = (2n + a + b + 1) (a^2 - b^2)
c   c3 = (2n + a + b + 2) (2n + a + b + 1) (2n + a + b)
c   c4 = -2(n + a) (n + b) (2n + a + b + 2)
c
c----------------------------------------------------------------------
c
      double precision function rJacobi_W(A,B,X,Nord)
c
#include "syscom.blk"
c
c  ...set up maximum degree supported
      parameter (nmax = 20)
c
c  ...local variables
      double precision, dimension(0:nmax) :: poly
c
c----------------------------------------------------------------------
c
c  ...test prints:
      iprint=1
c
      if (Nord.gt.nmax) then
        write(*,*)'rJacobi_W: Nord = ',Nord
        stop
      endif
c
c----------------------------------------------------------------------
c
c  ...P_0^(a,b), P_1^(a,b)
      poly(0) = 1.d0
      poly(1) = (2.d0*(A+1.d0) + (A+B+2.d0)*(X-1.d0))/2.d0
c
c  ...P_n^(a,b) , n >=2
      do n=1,Nord-1
c
c  .....coefficients
        c1 = 2.d0*(n+1.d0)*(n+A+B+1.d0)*(2.d0*n+A+B)
        c2 = (2.d0*n+A+B+1.d0)*(A**2-B**2)
        c3 = (2.d0*n+A+B+2.d0)*(2.d0*n+A+B+1.d0)*(2.d0*n+A+B)
        c4 = -2.d0*(n+A)*(n+B)*(2.d0*n+A+B+2.d0)
c  .....recurrence relation
        poly(n+1) = ((c2 + c3*X)*poly(n) + c4*poly(n-1))/c1
      enddo
c
      rJacobi_W = poly(Nord)
c
c----------------------------------------------------------------------
cccc     EXPLICIT FORMULAS FROM MATHEMATICA
ccc      select case(Nord)
ccc      case(0)
ccc        rJacobi_W = 1.d0
ccc      case(1)
ccc        rJacobi_W = (a - b + (2 + a + b)*x)/2.
ccc      case(2)
ccc        rJacobi_W =
ccc     .       ((1 + a)*(2 + a))/2. +
ccc     .    -  ((2 + a)*(3 + a + b)*(-1 + x))/2. +
ccc     .    -  ((3 + a + b)*(4 + a + b)*(-1 + x)**2)/8.
ccc      case(3)
ccc        rJacobi_W =
ccc     .       ((1 + a)*(2 + a)*(3 + a))/6. +
ccc     .    -  ((2 + a)*(3 + a)*(4 + a + b)*(-1 + x))/4. +
ccc     .    -  ((3 + a)*(4 + a + b)*(5 + a + b)*(-1 + x)**2)/8. +
ccc     .    -  ((4 + a + b)*(5 + a + b)*(6 + a + b)*(-1 + x)**3)/48.
ccc       case(4)
ccc        rJacobi_W =
ccc     .       ((1 + a)*(2 + a)*(3 + a)*(4 + a))/24. +
ccc     .    -  ((2 + a)*(3 + a)*(4 + a)*(5 + a + b)*(-1 + x))/12. +
ccc     .    -  ((3 + a)*(4 + a)*(5 + a + b)*(6 + a + b)*(-1 + x)**2)/
ccc     .    -   16. + ((4 + a)*(5 + a + b)*(6 + a + b)*(7 + a + b)*
ccc     .    -     (-1 + x)**3)/48. +
ccc     .    -  ((5 + a + b)*(6 + a + b)*(7 + a + b)*(8 + a + b)*
ccc     .    -     (-1 + x)**4)/384.
ccc      case(5)
ccc        rJacobi_W =
ccc     .       ((1 + a)*(2 + a)*(3 + a)*(4 + a)*(5 + a))/120. +
ccc     .    -  ((2 + a)*(3 + a)*(4 + a)*(5 + a)*(6 + a + b)*(-1 + x))/
ccc     .    -   48. + ((3 + a)*(4 + a)*(5 + a)*(6 + a + b)*(7 + a + b)*
ccc     .    -     (-1 + x)**2)/48. +
ccc     .    -  ((4 + a)*(5 + a)*(6 + a + b)*(7 + a + b)*(8 + a + b)*
ccc     .    -     (-1 + x)**3)/96. +
ccc     .    -  ((5 + a)*(6 + a + b)*(7 + a + b)*(8 + a + b)*
ccc     .    -     (9 + a + b)*(-1 + x)**4)/384. +
ccc     .    -  ((6 + a + b)*(7 + a + b)*(8 + a + b)*(9 + a + b)*
ccc     .    -     (10 + a + b)*(-1 + x)**5)/3840.
ccc      case(6)
ccc        rJacobi_W =
ccc     .       ((1 + a)*(2 + a)*(3 + a)*(4 + a)*(5 + a)*(6 + a))/
ccc     .    -   720. + ((2 + a)*(3 + a)*(4 + a)*(5 + a)*(6 + a)*
ccc     .    -     (7 + a + b)*(-1 + x))/240. +
ccc     .    -  ((3 + a)*(4 + a)*(5 + a)*(6 + a)*(7 + a + b)*
ccc     .    -     (8 + a + b)*(-1 + x)**2)/192. +
ccc     .    -  ((4 + a)*(5 + a)*(6 + a)*(7 + a + b)*(8 + a + b)*
ccc     .    -     (9 + a + b)*(-1 + x)**3)/288. +
ccc     .    -  ((5 + a)*(6 + a)*(7 + a + b)*(8 + a + b)*(9 + a + b)*
ccc     .    -     (10 + a + b)*(-1 + x)**4)/768. +
ccc     .    -  ((6 + a)*(7 + a + b)*(8 + a + b)*(9 + a + b)*
ccc     .    -     (10 + a + b)*(11 + a + b)*(-1 + x)**5)/3840. +
ccc     .    -  ((7 + a + b)*(8 + a + b)*(9 + a + b)*(10 + a + b)*
ccc     .    -     (11 + a + b)*(12 + a + b)*(-1 + x)**6)/46080.
ccc      case default
ccc        write(*,*)'rJacobi_W: Nord = ',Nord
ccc        stop
ccc      end select
c
c
      end function
c
c
c
c----------------------------------------------------------------------
c
c   routine name       - drJacobi_W
c
c----------------------------------------------------------------------
c
c   latest revision    - Oct 10
c
c   purpose            - evaluates derivative of Weighted Jacobi
c                        polynomial
c
c   arguments :
c     in:
c                A,B   - weights
c                X     - variable
c                Nord  - order of polynomial
c
c----------------------------------------------------------------------
c
      double precision function drJacobi_W(A,B,X,Nord)
c
#include "syscom.blk"
c
c----------------------------------------------------------------------
c
c   for Wolfram (and many other sources)
c
c   P_n^(a,b)' = (n + a + b + 1)/2 P_{n-1}^(a+1,b+1)
c
c----------------------------------------------------------------------
c
c
c  ...test prints:
      iprint = 0
c
c----------------------------------------------------------------------
      if (Nord.eq.0) then
        drJacobi_W = 0.d0
      else
        drJacobi_W = rJacobi_W(A+1,B+1,X,Nord-1)*(Nord+A+B+1.d0)/2.d0
      endif
c
c----------------------------------------------------------------------
cccc     EXPLICIT FORMULAS FROM MATHEMATICA
cccc
ccc      select case(Nord)
ccc      case(0)
ccc        drJacobi_W = 0.d0
ccc      case(1)
ccc        drJacobi_W = (2 + a + b)/2.
ccc      case(2)
ccc        drJacobi_W =
ccc     .    ((2 + a)*(3 + a + b))/2. +
ccc     . -  ((3 + a + b)*(4 + a + b)*(-1 + x))/4.
ccc      case(3)
ccc        drJacobi_W =
ccc     .    ((2 + a)*(3 + a)*(4 + a + b))/4. +
ccc     . -  ((3 + a)*(4 + a + b)*(5 + a + b)*(-1 + x))/4. +
ccc     . -  ((4 + a + b)*(5 + a + b)*(6 + a + b)*(-1 + x)**2)/16.
ccc      case(4)
ccc        drJacobi_W =
ccc     .    ((2 + a)*(3 + a)*(4 + a)*(5 + a + b))/12. +
ccc     . -  ((3 + a)*(4 + a)*(5 + a + b)*(6 + a + b)*(-1 + x))/8. +
ccc     . -  ((4 + a)*(5 + a + b)*(6 + a + b)*(7 + a + b)*
ccc     . -     (-1 + x)**2)/16. +
ccc     . -  ((5 + a + b)*(6 + a + b)*(7 + a + b)*(8 + a + b)*
ccc     . -     (-1 + x)**3)/96.
ccc      case(5)
ccc        drJacobi_W =
ccc     .    ((2 + a)*(3 + a)*(4 + a)*(5 + a)*(6 + a + b))/48. +
ccc     . -  ((3 + a)*(4 + a)*(5 + a)*(6 + a + b)*(7 + a + b)*
ccc     . -     (-1 + x))/24. +
ccc     . -  ((4 + a)*(5 + a)*(6 + a + b)*(7 + a + b)*(8 + a + b)*
ccc     . -     (-1 + x)**2)/32. +
ccc     . -  ((5 + a)*(6 + a + b)*(7 + a + b)*(8 + a + b)*
ccc     . -     (9 + a + b)*(-1 + x)**3)/96. +
ccc     . -  ((6 + a + b)*(7 + a + b)*(8 + a + b)*(9 + a + b)*
ccc     . -     (10 + a + b)*(-1 + x)**4)/768.
ccc      case(6)
ccc        drJacobi_W =
ccc     .    ((2 + a)*(3 + a)*(4 + a)*(5 + a)*(6 + a)*(7 + a + b))/
ccc     . -   240. + ((3 + a)*(4 + a)*(5 + a)*(6 + a)*(7 + a + b)*
ccc     . -     (8 + a + b)*(-1 + x))/96. +
ccc     . -  ((4 + a)*(5 + a)*(6 + a)*(7 + a + b)*(8 + a + b)*
ccc     . -     (9 + a + b)*(-1 + x)**2)/96. +
ccc     . -  ((5 + a)*(6 + a)*(7 + a + b)*(8 + a + b)*(9 + a + b)*
ccc     . -     (10 + a + b)*(-1 + x)**3)/192. +
ccc     . -  ((6 + a)*(7 + a + b)*(8 + a + b)*(9 + a + b)*
ccc     . -     (10 + a + b)*(11 + a + b)*(-1 + x)**4)/768. +
ccc     . -  ((7 + a + b)*(8 + a + b)*(9 + a + b)*(10 + a + b)*
ccc     . -     (11 + a + b)*(12 + a + b)*(-1 + x)**5)/7680.
ccc      case default
ccc        write(*,*)'drJacobi: Nord = ',Nord
ccc        stop
ccc      endselect
c
c
      end function
c
c
c
c----------------------------------------------------------------------
c
c   routine name       - rJacobi_W_I
c
c----------------------------------------------------------------------
c
c   latest revision    - Oct 10
c
c   purpose            - evaluates Weighted Integrated Jacobi
c                        polynomial
c
c   arguments :
c     in:
c                A     - weight
c                X     - variable
c                Nord  - order of polynomial
c
c----------------------------------------------------------------------
c
      double precision function rJacobi_W_I(A,X,Nord)
c
#include "syscom.blk"
c
c----------------------------------------------------------------------
c
c   p_n^a := P_n^(a,0)
c
c   \hat{p}_0^a = 1
c   \hat{p}_n^a = \int_{-1}^x p_{n-1}^a , n >= 1
c
c   Recurrence relations from Beuchler's papers:
c
c   n >= 2
c   \hat{p}_n^a = c1 p_n^a + c2 p_{n-1}^a + c3 p_{n-2}^a
c
c   c1 = 2(n + a)/[(2n + a - 1) (2n + a)]
c   c2 = 2a/[(2n + a - 2) (2n + a)]
c   c3 = -2(n - 1)/[(2n + a - 1) (2n + a - 2)]
c
c   n = 1
c   \hat{p}_n^a = c1 ( p_n^{a-1} + p_{n-1}^{a-1} )
c   c1 = 1/(2n + a - 1)
c
c----------------------------------------------------------------------
c
c  ...printing flag
      iprint = 0
c
c----------------------------------------------------------------------
ccc      if (A.eq.0.d0) then
ccc        if (Nord.lt.2) then
ccc          write(*,7001) A,Nord
 7001     format('rJacobi_W_I: A,Nord = ',e25.15,i3)
ccc          stop 1
ccc        endif
ccc        rJacobi_W_I = (1.d0-X**2)/(2.d0*Nord-2.d0)
ccc     .               *rJacobi_W(1.d0,1.d0,X,Nord-2)
ccc      elseif (A.ge.1.d0) then
c
c  ...eq. (7) from Beuchler paper
      if (A.ge.-1.d0) then
        if (Nord.lt.1) then
          write(*,7001) A,Nord
          stop 1
        endif
        rJacobi_W_I = (1.d0+X)/Nord
     .               *rJacobi_W(A-1.d0,1.d0,X,Nord-1)
      else
        write(*,7001) A,Nord
        stop
      endif
c
c----------------------------------------------------------------------
cccc     EXPLICIT FORMULAS FROM MATHEMATICA
ccc      select case(Nord)
ccc      case(0)
ccc        rJacobi_W_I = 1.d0
ccc      case(1)
ccc        rJacobi_W_I =
ccc     .     1 + x
ccc      case(2)
ccc        rJacobi_W_I =
ccc     .   -0.5 + a/4. + (a*x)/2. + x**2/2. + (a*x**2)/4.
ccc      case(3)
ccc        rJacobi_W_I =
ccc     .   (-5*a)/24. + a**2/24. - x/2. - (a*x)/8. + (a**2*x)/8. +
ccc     . -  (3*a*x**2)/8. + (a**2*x**2)/8. + x**3/2. +
ccc     . -  (7*a*x**3)/24. + (a**2*x**3)/24.
ccc      case(4)
ccc        rJacobi_W_I =
ccc     .   (a**3*(1 + x)**4 + 3*a**2*(1 + x)**3*(-3 + 5*x) +
ccc     . -    24*(1 - 6*x**2 + 5*x**4) +
ccc     . -    2*a*(1 + x)**2*(1 + x*(-34 + 37*x)))/192.
ccc      case(5)
ccc        rJacobi_W_I =
ccc     .    (a**4*(1 + x)**5 + 2*a**3*(1 + x)**4*(-7 + 13*x) +
ccc     . -    240*x*(3 - 10*x**2 + 7*x**4) +
ccc     . -    a**2*(1 + x)**3*(11 + x*(-218 + 251*x)) +
ccc     . -    2*a*(1 + x)**2*(73 + x*(-41 + x*(-541 + 533*x))))/
ccc     . -  1920.
ccc      case(6)
ccc        rJacobi_W_I =
ccc     .   (a**5*(1 + x)**6 + 20*a**4*(1 + x)**5*(-1 + 2*x) +
ccc     . -    1440*(-1 + 15*x**2 - 35*x**4 + 21*x**6) +
ccc     . -    100*a**2*(1 + x)**3*(5 + x**2*(-51 + 50*x)) +
ccc     . -    5*a**3*(1 + x)**4*(7 + x*(-106 + 127*x)) +
ccc     . -    12*a*(1 + x)**2*
ccc     . -     (-23 + x*(558 + x*(-388 + x*(-1742 + 1627*x)))))/
ccc     . -  23040.
ccc      case default
ccc        write(*,*)'rJacobi_W_I: Nord = ',Nord
ccc        stop
ccc      end select
c
c----------------------------------------------------------------------
c
      return
c
c

c     O  L  D     I  M  P   L   E   M   E   N   T  A   T   I  O  N
c
      select case (Nord)
      case(0)
        rJacobi_W_I = 1.d0
c  ...use a 2 terms recurrence formula for \hat{p}_1^a
      case(1)
        c1 = 1.d0/(2.d0*Nord + A - 1.d0)
c
        rJacobi_W_I = c1*rJacobi_W(A-1.d0,0.d0,X,Nord)
     .             + c1*rJacobi_W(A-1.d0,0.d0,X,Nord-1)
c  ...use a 3 terms recurrence formula for \hat{p}_n^a , n >=2
      case default
        c1 = 2.d0*(A+Nord)/((2.d0*Nord+A-1.d0)*(2.d0*Nord+A))
        c2 = 2.d0*A/((2.d0*Nord+A-2.d0)*(2.d0*Nord+A))
        c3 = -2.d0*(Nord-1.d0)/((2.d0*Nord+A-1.d0)*(2.d0*Nord+A-2.d0))
c
        rJacobi_W_I = c1*rJacobi_W(A,0.d0,X,Nord) +
     .                c2*rJacobi_W(A,0.d0,X,Nord-1) +
     .                c3*rJacobi_W(A,0.d0,X,Nord-2)
      end select
c
      end function
c
c
c
c----------------------------------------------------------------------
c
c   routine name       - drJacobi_W_I
c
c----------------------------------------------------------------------
c
c   latest revision    - Oct 10
c
c   purpose            - evaluates derivative of Weighted Integrated
c                        Jacobi polynomial
c
c   arguments :
c     in:
c                A     - weight
c                X     - variable
c                Nord  - order of polynomial
c
c----------------------------------------------------------------------
c
      double precision function drJacobi_W_I(A,X,Nord)
c
#include "syscom.blk"
c
c----------------------------------------------------------------------
c
c  ...test prints:
      iprint = 0
c
c----------------------------------------------------------------------
      if (A.ge.-1.d0) then
        if (Nord.lt.1) then
          write(*,7001) A,Nord
 7001     format('drJacobi_W_I: A,Nord = ',e25.15,i3)
          stop 1
        endif
        drJacobi_W_I = 1.d0/Nord
     .               *rJacobi_W(A-1.d0,1.d0,X,Nord-1)
     .               + (1.d0+X)/Nord
     .               *drJacobi_W(A-1.d0,1.d0,X,Nord-1)
      else
        write(*,7001) A,Nord
        stop
      endif
c
c----------------------------------------------------------------------
cccc     EXPLICIT FORMULAS FROM MATHEMATICA
ccc      select case(Nord)
ccc      case(0)
ccc        drJacobi_W_I = 0.d0
ccc      case(1)
ccc        drJacobi_W_I = 1.d0
ccc      case(2)
ccc        drJacobi_W_I =
ccc     . a/2. + x + (a*x)/2.
ccc      case(3)
ccc        drJacobi_W_I =
ccc     .  -0.5 - a/8. + a**2/8. + (3*a*x)/4. + (a**2*x)/4. +
ccc     . -  (3*x**2)/2. + (7*a*x**2)/8. + (a**2*x**2)/8.
ccc     . -  23040.
ccc      case(4)
ccc        drJacobi_W_I =
ccc     .  (15*a**2*(1 + x)**3 + 4*a**3*(1 + x)**3 +
ccc     . -    9*a**2*(1 + x)**2*(-3 + 5*x) +
ccc     . -    2*a*(1 + x)**2*(-34 + 74*x) + 24*(-12*x + 20*x**3) +
ccc     . -    4*a*(1 + x)*(1 + x*(-34 + 37*x)))/192.
ccc      case(5)
ccc        drJacobi_W_I =
ccc     .   (26*a**3*(1 + x)**4 + 5*a**4*(1 + x)**4 +
ccc     . -    8*a**3*(1 + x)**3*(-7 + 13*x) +
ccc     . -    a**2*(1 + x)**3*(-218 + 502*x) +
ccc     . -    240*x*(-20*x + 28*x**3) +
ccc     . -    240*(3 - 10*x**2 + 7*x**4) +
ccc     . -    3*a**2*(1 + x)**2*(11 + x*(-218 + 251*x)) +
ccc     . -    2*a*(1 + x)**2*(-41 + x*(-541 + 533*x) +
ccc     . -       x*(-541 + 1066*x)) +
ccc     . -    4*a*(1 + x)*(73 + x*(-41 + x*(-541 + 533*x))))/1920.
ccc      case(6)
ccc        drJacobi_W_I =
ccc     .   (40*a**4*(1 + x)**5 + 6*a**5*(1 + x)**5 +
ccc     . -    100*a**4*(1 + x)**4*(-1 + 2*x) +
ccc     . -    5*a**3*(1 + x)**4*(-106 + 254*x) +
ccc     . -    1440*(30*x - 140*x**3 + 126*x**5) +
ccc     . -    100*a**2*(1 + x)**3*(50*x**2 + 2*x*(-51 + 50*x)) +
ccc     . -    300*a**2*(1 + x)**2*(5 + x**2*(-51 + 50*x)) +
ccc     . -    20*a**3*(1 + x)**3*(7 + x*(-106 + 127*x)) +
ccc     . -    12*a*(1 + x)**2*
ccc     . -     (558 + x*(-388 + x*(-1742 + 1627*x)) +
ccc     . -       x*(-388 + x*(-1742 + 1627*x) + x*(-1742 + 3254*x)))
ccc     . -      + 24*a*(1 + x)*
ccc     . -     (-23 + x*(558 + x*(-388 + x*(-1742 + 1627*x)))))/
ccc     . -  23040.
ccc      case default
ccc        write(*,*)'drJacobi_W_I: Nord = ',Nord
ccc        stop
ccc      end select
c
      return
c
c
c----------------------------------------------------------------------
c     O  L  D     I  M  P  L  E  M  E  N  T  A  T  I  O  N
      select case (Nord)
      case(0)
        drJacobi_W_I = 0.d0
c  ...use a 2 terms recurrence formula for \hat{p}_1^a
      case(1)
        c1 = 1.d0/(2.d0*Nord+A-1.d0)
c
        drJacobi_W_I = c1*drJacobi_W(A-1.d0,0.d0,X,Nord)
     .              + c1*drJacobi_W(A-1.d0,0.d0,X,Nord-1)
c  ...use a 3 terms recurrence formula for \hat{p}_n^a , n >=2
      case default
        c1 = 2.d0*(A+Nord)/((2.d0*Nord+A-1.d0)*(2.d0*Nord+A))
        c2 = 2.d0*A/((2.d0*Nord+A-2.d0)*(2.d0*Nord+A))
        c3 = -2.d0*(Nord-1.d0)/((2.d0*Nord+A-1.d0)*(2.d0*Nord+A-2.d0))
c
        drJacobi_W_I = c1*drJacobi_W(A,0.d0,X,Nord) +
     .                c2*drJacobi_W(A,0.d0,X,Nord-1) +
     .                c3*drJacobi_W(A,0.d0,X,Nord-2)
      end select
c
      end function
