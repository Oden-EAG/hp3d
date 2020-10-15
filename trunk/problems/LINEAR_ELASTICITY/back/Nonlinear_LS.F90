!
!----------------------------------------------------------------------
!
!   routine name       - fit
!
!----------------------------------------------------------------------
!
!   latest revision    - Sep 13
!
!   purpose            - perform a nonlinear weighted least squares fit
!
!   arguments:
!     in:
!             Nrp      - number of data points to fit
!             Data     - values to fit
!             Weight   - weights
!             N        - number of unknowns
!     in/out:
!             X        - initial/final value of the argument
!
!----------------------------------------------------------------------
!
      subroutine fit(Nrp,Data,Weight,X,N)
!
!
#include "syscom.blk"
#include "cinout.blk"
!
      dimension Data(Nrp),Weight(Nrp),X(N)
!
!  ...gradient, descent direction
      dimension grad(N),dx(N)
      dimension x_n(N),grad_n(N),dx_n(N)
!
      iprint=0
!
!  ...initiate backtracking parameters
      rho = .5d0; c = .1d0; eps = 1.d-15
      call eval(Nrp,Data,Weight,X,N, val,grad,dx)
!
!  ...iterate
      do iter=1,100
        if (iprint.eq.1) then
          write(*,7001) iter,X(1:N)
 7001     format('fit: iter = ',i4,' X = ',10e12.5)
          write(*,7002) val
 7002     format('     value = ',e20.15)
          call pause
        endif
!
!  .....initiate
        alpha = .1d0
        s = 0.d0
        do i=1,N
          s = s + grad(i)*dx(i)
        enddo
        do
          x_n = X + alpha*dx
          call eval(Nrp,Data,Weight,x_n,N, val_n,grad_n,dx_n)
          if (iprint.eq.1) then
            write(*,7003) val_n,alpha
 7003       format('     val_n,alpha = ',2e12.5)
          endif
          if (val_n.le.val+c*alpha*s) exit
          alpha = alpha*rho
        enddo
        diff = abs(val-val_n)
        x = x_n; val = val_n; grad = grad_n; dx = dx_n
        if (diff.lt.eps) exit
      enddo
      if (iprint.ge.1) then
        write(*,7010) X(1:N)
 7010   format('fit: X = ',10e12.5)
      endif
!
      end subroutine
!
!----------------------------------------------------------------------
!
!   routine name       - eval
!
!----------------------------------------------------------------------
!
!   latest revision    - Sep 13
!
!   purpose            - evaluate a weighted least squares misfit
!                        functions, its gradient, and Guass-Newton
!                        direction
!
!   arguments:
!     in:
!             Nrp      - number of data points to fit
!             Data     - values to fit
!             Weight   - weights
!             X        - arguments
!             N        - number of arguments
!     out:
!             Fit      - value of the misfit function at the point
!             Grad     - the corresponding gradient
!             Dx       - direction for minimization
!
!----------------------------------------------------------------------
!
      subroutine eval(Nrp,Data,Weight,X,N, Fit,Grad,Dx)
c
c
#include "syscom.blk"
#include "cinout.blk"
!
      dimension Data(Nrp),Weight(Nrp),X(N),Grad(N),Dx(N)
      dimension fgrad(N),hesj(N,N)
      character*1 uplo
!
      iprint=0
 10   continue
!
!
      Fit = 0.d0; Grad = 0.d0
      hesj = 0.d0
      do np=1,Nrp
        call funct(np,N,X, f,fgrad)
        res = f - Data(np)
        Fit = Fit + .5d0*weight(np)*res**2
        Grad(1:N) = Grad(1:N) + weight(np)*res*fgrad(1:N)
        do m=1,N
          hesj(m,1:N) = hesj(m,1:N) + weight(np)*fgrad(m)*fgrad(1:N)
        enddo
      enddo
      if (iprint.eq.1) then
        write(*,7011) X
 7011   format('eval: X = ',10e12.5)
        write(*,7012) Fit
 7012   format('      Value = ',e20.15)
        write(*,7013) Grad(1:N)
 7013   format('      Grad = ',10e12.5)
        write(*,7014)
 7014   format('      hesj = ')
        do i=1,N
          write(*,7015) hesj(i,1:N)
 7015     format(10e12.5)
        enddo
        call pause
      endif
!
!  ...use Gauss-Newton to determine direction to march
      Dx = -Grad
      uplo = 'U'
      nrhs=1
!
      do i=1,N
        hesj(i,i) = hesj(i,i) + 1.d-5
      enddo
      call DPOSV(uplo,N,nrhs,hesj,N,Dx,N,info)
      if (info.ne.0) then
        write(*,*) 'eval: info = ',info
        iprint=1
        go to 10
        stop1
      endif
!
      if (iprint.eq.1) then
        write(*,7001) X
 7001   format('eval: X   =  ',10e12.5)
        write(*,7002) Fit
 7002   format('      Fit =  ',e12.5)
        write(*,7003) Grad(1:N)
 7003   format('      Grad = ',10e12.5)
        write(*,7004) Dx(1:N)
 7004   format('      Dx   = ',10e12.5)
        call pause
      endif
!
      end subroutine
!
!----------------------------------------------------------------------
!
!   routine name       - funct
!
!----------------------------------------------------------------------
!
!   latest revision    - Sep 13
!
!   purpose            - return value and gradient of a function
!
!   arguments:
!     in:
!             h        - meshsize
!             N        - number of unknowns
!             X        - arguments
!     out:
!             F        - function value
!             Grad     - the corresponding gradient
!
!----------------------------------------------------------------------
!
subroutine funct(h,N,X, F,Grad)
!
  implicit none
!----------------------------------------------------------------------
  real*8,  intent(in)  :: h
  integer, intent(in)  :: N
  real*8,  intent(in)  :: X(3)
  real*8,  intent(out) :: F
  real*8,  intent(out) :: Grad(3)
!----------------------------------------------------------------------
  integer :: iprint
  real*8  :: intErr,c,s,hs
!----------------------------------------------------------------------
!
  dimension X(N),Grad(N)
  iprint=0
!
  if (N.ne.3) then
    write(*,*) 'funct: N = ' N
  endif
!
  intErr = X(1)
  c      = X(2)
  s      = X(3)
!
  hs = h**s
  F = intErr + c*hs
  Grad(1) = 1.d0
  Grad(2) = hs
  Grad(3) = c*s*h**(s-1.d0)
  if (iprint.eq.1) then
    write(*,7001) h,X(1:N)
7001   format('funct: h,X = ',e12.5,2x,3e12.5)
    write(*,7002) F,Grad(1:N)
7002   format('       F = ',e20.5,' Grad = ',3e12.5)
    call pause
  endif
!
end subroutine
