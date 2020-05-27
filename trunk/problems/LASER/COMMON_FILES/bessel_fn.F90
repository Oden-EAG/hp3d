!
#include "implicit_none.h"
!
!---------------------------------------------------------------------------------
!
!!!!!!!!!!!!!!!!!!ROUTINES FOR BESSEL FUNCTIONS FROM JAKE!!!!!!!!!!!!!!!!!
!
! File:          BesselFunctions.f95
!
! Author:        Jacob Grosek
! Start Date:    March 31, 2017
! Last Modified: March 31, 2017
! Air Force Research Laboratory at Kirtland AFB, Albuquerque, NM
! Directed Energy Directorate, Laser Division, Modeling & Simulation Program
!
! Description of the File:
!
! This file contains Fortran Bessel function algorithms, and their
! derivatives.  These Bessel functions only accept real values, and they
! only output real values as well.
!
! This code is modified from "Numerical Recipes in Fortran 77: The Art
! of Scientific Computing," Cambridge University Press (1992)
!
! File Parameters/Variables:
!
! x
!   - (input, real, [ ]) the given point at which the Bessel
!     functions will be calculated
! order
!   - (input, real, [ ]) the order of the Bessel functions
! bessJ
!   - (output, real, [ ]) the value of the Bessel J function
! bessY
!   - (output, real, [ ]) the value of the Bessel Y function
! bessJder
!   - (output, real, [ ]) the value of the derivative of the
!     Bessel J function
! bessYder
!   - (output, real, [ ]) the value of the derivative of the
!     Bessel Y function
!
! File Dependencies:
!
! This subroutine calls upon the subroutine beschbgamma defined in
! this same file.  This subroutine returns the Bessel functions
! bessJ = J_order and bessY = Y_order, and their respective
! derivatives bessJder = J_order′ and bessYder = Y_order′, for
! positive x and of an order greater or equal to 0.  The relative
! accuracy is within one or two significant digits of eps, except
! near a zero of one of the functions, where eps controls its
! absolute accuracy. The parameter "floatingpointmin" is a number
! close to the machine’s smallest floating-point number. All
! internal arithmetic is accomplished with real numbers. In oder to
! convert the entire routine to double precision, use the real(8)
! declaration and decrease eps to 10**−16.  Also convert the
! subroutine beschbgamma.
!
!---------------------------------------------------------------------------------
!
subroutine dbessJY(x, order, bessJ, bessY, bessJder, bessYder)

  implicit none

  integer, parameter :: maxit = 10000
  real(8) :: bessJ,bessJder, bessY, bessYder, x, order
  real(8), parameter :: eps = 1.d-16, floatingpointmin = 1.e-30, &
    pi = 3.141592653589793d0, xmin = 2.d0, &
    emc = 0.577215664901533d0
  integer :: i, signi, l, nl
  real(8) :: a, b, br, bi, c, cr, ci, d, del, del1, den, &
    di, dlr, dli, dr, e, f, fact, fact2, fact3, ff, gam, &
    gam1, gam2, gammi, gampl, h, p, pimu, pimu2, q, r, &
    bessJl, bessJl1, bessJmu, bessJder1, bessJderl, &
    bessJtemp, bessY1, bessYmu, bessYmup, bessYtemp, summ, &
    summ1, temp, w, x2, xi, xi2, xmu, xmu2, gampl1, negate

  if (order < 0.d0) then
    order  = dabs(order)
    negate = (-1.d0)**int(order)
  else
    negate = 1.d0
  endif ! end if statement
  if (dabs(x) <= 1.d-6) then
    if ((dabs(x) <= 1.d-17) .and. (order == 0.d0)) then
      bessJ    = 1.d0 *  negate
      bessJder = 0.d0
      bessY    = -1.d30 * negate
      bessYder = 1.d30 * negate
    elseif ((dabs(x) <= 1.d-17) .and. (order /= 0.d0)) then
      bessJ    = 0.d0
      bessJder = negate * 2.d0 * order * &
        (x / 2.d0)**(order - 1.d0) / gamma(order + 1.d0)
      bessY    = -1.d30 * negate
      bessYder = 1.d30 * negate
    elseif ((dabs(x) > 1.d-17) .and. (order == 0.d0)) then
      bessJder = negate * 2.d0 * order * &
        (x / 2.d0)**(order - 1.d0) / gamma(order + 1.d0)
      bessY    = negate * (2.d0 / pi) * (dlog(x / 2.d0) + emc)
      bessYder = negate * 2.d0 / (pi * x)
    elseif ((dabs(x) > 1.d-17) .and. (order /= 0.d0)) then
      bessJ    = negate * (x / 2.d0)**order / &
        gamma(order + 1.d0)
      bessJder = negate * 2.d0 * order * &
        (x / 2.d0)**(order - 1.d0) / gamma(order + 1.d0)
      bessY    = -negate * (gamma(order) / pi) * &
        (2.d0 / x)**order
      bessYder = negate * (order * gamma(order) / (pi * x)) * &
        (2.d0 / x)**order
    endif ! end if statement
    return
  endif ! end if statement

  ! Here "nl" is the number of downward recurrences of the J’s and
  ! upward recurrences of Y’s. "xmu" lies between −1/2 and 1/2 for
  ! x < xmin, while it is chosen so that x is greater than the
  ! turning point for x >= xmin:
  if (x < xmin) then
    nl = int(order + 0.5d0)
  else
    nl = max(0, int(order - x + 1.5d0))
  endif ! end if statement
  xmu  = order - nl
  xmu2 = xmu * xmu
  if (x /= 0.d0) then
    xi = 1.d0 / x
  else
    xi = 1.d30
  endif ! end if statement
  xi2  = 2.d0 * xi
  w    = xi2 / pi
    ! the Wronskian

  ! Here the first continued fraction is found by the modified
  ! Lentz's method.  The variable "signi" keeps track of sign
  ! changes in the denominator:
  signi = 1
  h     = order * xi
  if (h < floatingpointmin) then
    h = floatingpointmin
  endif ! end if statement
  b = xi2 * order
  d = 0.d0
  c = h
  do i = 1, maxit
    b = b + xi2
    d = b - d
    if (abs(d) < floatingpointmin) then
      d = floatingpointmin
    endif ! end if statement
    c = b - 1.d0 / c
    if (abs(c) < floatingpointmin) then
      c = floatingpointmin
    endif ! end if statement
    d   = 1.d0 / d
    del = c * d
    h   = del * h
    if (d < 0.d0) then
      signi = -signi
    endif ! end if statement
    if (dabs(del - 1.d0) < eps) then
      goto 1
    endif ! end if statement
  enddo
  print *, "x is too large in dbessJY; try asymptotic expansion"
  print *, "x = ", x
  1 continue

  ! Here J_order and J_order' are initialized for downward
  ! recurrence:
  bessJl    = signi * floatingpointmin
  bessJderl = h * bessJl
  bessJl1   = bessJl
    ! store values for future rescaling
  bessJder1 = bessJderl
  fact      = order * xi
  do l = nl, 1, -1
    bessJtemp = fact * bessJl + bessJderl
    fact      = fact - xi
    bessJderl = fact * bessJtemp - bessJl
    bessJl    = bessJtemp
  enddo
  if (bessJl == 0.d0) then
    bessJl = eps
  endif ! end if statement

  ! The subroutine already has calculated unnormalized J_order and
  ! J_order':
  f = bessJderl / bessJl
  if (x < xmin) then ! use series
    x2   = 0.5d0 * x
    pimu = pi * xmu
    if (abs(pimu) < eps) then
      fact = 1.d0
    else
      fact = pimu / dsin(pimu)
    endif ! end if statement
    d = -dlog(x2)
    e = xmu * d
    if (abs(e) < eps) then
      fact2 = 1.d0
    else
      fact2 = dsinh(e) / e
    endif ! end if statement
    call beschbgamma(xmu, gam1, gam2, gampl, gammi)
      ! Chebyshev evaluation of Gamma1 and Gamma2
    ff    = 2.d0 / pi * fact * (gam1 * dcosh(e) + gam2 * fact2 * d)
    e     = dexp(e)
    p     = e  / (gampl * pi)
    q     = 1.d0 / (e * pi * gammi)
    pimu2 = 0.5d0 * pimu
    if (abs(pimu2) < eps) then
      fact3 = 1.d0
    else
      fact3 = dsin(pimu2) / pimu2
    endif ! end if statement
    r     = pi * pimu2 * fact3 * fact3
    c     = 1.d0
    d     = -x2 * x2
    summ  = ff + r * q
    summ1 = p
    do i = 1, maxit
      ff    = (i * ff + p + q) / (i * i - xmu2)
      c     = c * d / i
      p     = p / (i - xmu)
      q     = q / (i + xmu)
      del   = c * (ff + r * q)
      summ  = summ + del
      del1  = c * p - i * del
      summ1 = summ1 + del1
      if (abs(del) < (1.d0 + abs(summ)) * eps) then
        goto 2
      endif ! end if statement
    enddo
    print *, "bessY series failed to converge"
2   continue
    bessYmu  = -summ
    bessY1   = -summ1 * xi2
    bessYmup = xmu * xi * bessYmu - bessY1
    bessJmu  = w / (bessYmup - f * bessYmu)
  else ! evaluate the second continued fraction by Lentz's method
    a    = 0.25d0 - xmu2
    p    = -0.5d0 * xi
    q    = 1.d0
    br   = 2.d0 * x
    bi   = 2.d0
    fact = a * xi / (p * p + q * q)
    cr   = br + q * fact
    ci   = bi + p * fact
    den  = br * br + bi * bi
    dr   = br / den
    di   = -bi / den
    dlr  = cr * dr - ci * di
    dli  = cr * di + ci * dr
    temp = p * dlr - q * dli
    q    = p * dli+q * dlr
    p    = temp
    do i = 2, maxit
      a  = a + 2 * (i - 1)
      bi = bi + 2.d0
      dr = a * dr + br
      di = a * di + bi
      if ((abs(dr) + abs(di)) < floatingpointmin) then
        dr = floatingpointmin
      endif ! end if statement
      fact = a / (cr * cr +  ci * ci)
      cr   = br + cr * fact
      ci   = bi - ci * fact
      if ((abs(cr) + abs(ci)) < floatingpointmin) then
        cr = floatingpointmin
      endif ! end if statement
      den  = dr * dr + di * di
      dr   = dr / den
      di   = -di / den
      dlr  = cr * dr - ci * di
      dli  = cr * di + ci * dr
      temp = p * dlr - q * dli
      q    = p * dli + q * dlr
      p    = temp
      if ((abs(dlr - 1.d0) + abs(dli)) < eps) then
        goto 3
      endif ! end if statement
    enddo
    print *, "the second continued fraction failed in dbessJY"
3   continue
    gam    = (p - f) / q
    bessJmu  = dsqrt(w / ((p - f) * gam + q))
    bessJmu  = sign(bessJmu, bessJl)
    bessYmu  = bessJmu * gam
    bessYmup = bessYmu * (p + q / gam)
    bessY1   = xmu * xi * bessYmu - bessYmup
  endif ! end if statement
  fact     = bessJmu / bessJl
  bessJ    = negate * bessJl1 * fact
    ! scale original J_order and J_order′
  bessJder = bessJder1 * fact

  ! Here is the upward recurrence of Y_order
  do i = 1, nl
    bessYtemp = (xmu + i) * xi2 * bessY1 - bessYmu
    bessYmu   = bessY1
    bessY1    = bessYtemp
  enddo
  bessY    = negate * bessYmu
  bessYder = order * xi * bessYmu - bessY1

  return
end subroutine dbessJY
!
!---------------------------------------------------------------------------------
!
! This subroutine evaluates Gamma1 and Gamma2 by a Chebyshev
! expansion for |x| <= 1/2. Also returns 1 / Gamma(1 + x) and
! 1 / Gamma(1 − x). If converting to real, set nuse1 = 7 and
! nuse2 = 8.  This subroutine calls the chebev function, which is
! contained in this file.
! Subroutine Inputs/Output:
! x     - (input, real, [ ]) the point at which the gamma functions
!         are evaluated
! gam1  - (output, real, [ ]) the function value of Gamma1
! gam2  - (output, real, [ ]) the function value of Gamma2
! gampl - (output, real, [ ]) the function value of 1 / Gamma(1 + x)
! gammi - (output, real, [ ]) the function value of 1 / Gamma(1 − x)
!
subroutine beschbgamma(x, gam1, gam2, gampl, gammi)

  implicit none

  integer, parameter :: nuse1 = 7, nuse2 = 8
  real(8) :: gam1, gam2, gammi, gampl, x
  real(8) :: xx, c1(7), c2(8)
  real(8), external :: chebev
  save c1, c2
  data c1 /-1.142022680371168d0, 6.5165112670737d-3, &
    3.087090173086d-4, -3.4706269649d-6, 6.9437664d-9, &
    3.67795d-11, -1.356d-13/
  data c2 /1.843740587300905d0, -7.68528408447867d-2, &
    1.2719271366546d-3, -4.9717367042d-6, -3.31261198d-8, &
    2.423096d-10, -1.702d-13, -1.49d-15/

  xx    = 8.d0 * x * x - 1.d0
    ! x is multiplied by 2 in order to change its range to
    ! -1 to 1, and then another transformation is applied in order
    ! to evaluate using the even Chebyshev series; since the
    ! function is even it would be wasteful to call the chebev
    ! function with all the odd coefficients being zero; an
    ! approximation of an even function on [-1,1] will only
    ! involve even Chebyshev polynomials.  Once x is in the range
    ! -1 to 1, then the transformation being applied is
    ! T_2n(x) = T_n(2 * x**2 - 1)
  gam1  = chebev(-1.d0, 1.d0, c1, nuse1, xx)
  gam2  = chebev(-1.d0, 1.d0, c2, nuse2, xx)
  gampl = gam2 - x * gam1
  gammi = gam2 + x * gam1

  return
end subroutine beschbgamma
!
!---------------------------------------------------------------------------------
!
! This subroutine calls upon the subroutine beschbgamma defined in
! this same file.  This subroutine returns the modified Bessel
! functions bessI = I_order and bessK = K_order, and their
! respective derivatives bessIder = I_order′ and
! bessKder = K_order′, for positive x and for an order greater or
! equal to 0.  The relative accuracy is within one or two
! significant digits of eps. The parameter "floatingpointmin" is
! a number close to the machine’s smallest floating point number.
! All internal arithmetic is accomplished with real numbers. In
! order to convert the entire routine to double precision, use
! the real(8) declaration and decrease eps to 10**−16.  Also
! convert the subroutine beschbgamma.
! Subroutine Inputs/Output:
! x        - (input, real, [ ]) the given point at which the Bessel
!         functions will be calculated
! order    - (input, real, [ ]) the order of the Bessel functions
! bessI    - (output, real, [ ]) the value of the Bessel I function
! bessK    - (output, real, [ ]) the value of the Bessel K function
! bessIder - (output, real, [ ]) the value of the derivative of the
!         Bessel I function
! bessKder - (output, real, [ ]) the value of the derivative of the
!         Bessel K function
!
subroutine dbessIK(x, order, bessI, bessK, bessIder, bessKder)

  implicit none

  integer, parameter :: maxit = 10000
  real(8) :: bessI, bessIder, bessK, bessKder, x, order
  real(8), parameter :: eps = 1.d-16, floatingpointmin = 1.e-30, &
    pi = 3.141592653589793d0, xmin = 2.d0
  integer :: i, l, nl
  real(8) :: a, a1, b, c, d, del, del1, delh, dels, e, f, &
    fact, fact2, ff, gam1, gam2, gammi, gampl, h, p, pimu, &
    q, q1, q2, qnew, bessIl, bessIl1, bessImu, bessIder1, &
    bessIderl, bessItemp, bessK1, bessKmu, bessKmup, &
    bessKtemp, s, summ, summ1, x2, xi, xi2, xmu, xmu2

  if ((x <= 0.d0) .or. (order < 0.d0)) then
    print *, "Bad arguments were inputed into dbessIK."
    print *, "x = ", x
    print *, "order = ", order
  endif ! end if statement

  ! Here "nl" is the number of downward recurrences of the I’s and
  ! upward recurrences of K’s. "xmu" lies between −1/2 and 1/2:
  nl   = int(order + 0.5d0)
  xmu  = order - nl
  xmu2 = xmu * xmu
  xi   = 1.d0 / x
  xi2  = 2.d0 * xi
  h    = order * xi

  ! Here the first continued fraction is found by the modified
  ! Lentz's method:
  if (h < floatingpointmin) then
    h = floatingpointmin
  endif ! end if statement
  b = xi2 * order
  d = 0.d0
  c = h
  do i = 1, maxit
    b   = b + xi2
    d   = 1.d0 / (b + d)
      ! denominators cannot be zero here, so there is no need for
      ! special precautions
    c   = b + 1.d0 / c
    del = c * d
    h   = del * h
    if (abs(del - 1.d0) < eps) then
      goto 10
    endif ! end if statement
  enddo
  print *, "x is too large in dbessIK; try an asymptotic expansion"
  print *, "x = ", x
10  continue

  ! Here I_order and I_order' are initialized for downward
  ! recurrence:
  bessIl    = floatingpointmin
  bessIderl = h * bessIl
  bessIl1   = bessIl
  bessIder1 = bessIderl
    ! store values for future rescaling
  fact      = order * xi
  do l = nl, 1, -1
    bessItemp = fact * bessIl + bessIderl
    fact      = fact - xi
    bessIderl = fact * bessItemp + bessIl
    bessIl    = bessItemp
  enddo

  ! The subroutine already has calculated unnormalized I_order and
  ! I_order':
  f = bessIderl / bessIl
  if (x <= xmin) then ! use series
    x2   = 0.5d0 * x
    pimu = pi * xmu
    if (abs(pimu) < eps) then
      fact = 1.d0
    else
      fact = pimu / dsin(pimu)
    endif ! end if statement
    d = -dlog(x2)
    e = xmu * d
    if (abs(e) < eps) then
      fact2 = 1.d0
    else
      fact2 = dsinh(e) / e
    endif ! end if statement
    call beschbgamma(xmu, gam1, gam2, gampl, gammi)
      ! Chebyshev evaluation of Gamma1 and Gamma2
    ff    = fact * (gam1 * dcosh(e) + gam2 * fact2 * d)
    summ  = ff
    e     = dexp(e)
    p     = 0.5d0 * e / gampl
    q     = 0.5d0 / (e * gammi)
    c     = 1.d0
    d     = x2 * x2
    summ1 = p
    do i = 1, maxit
      ff    = (i * ff + p + q) / (i * i - xmu2)
      c     = c * d / i
      p     = p / (i - xmu)
      q     = q / (i + xmu)
      del   = c * ff
      summ  = summ + del
      del1  = c * (p - i * ff)
      summ1 = summ1 + del1
      if (abs(del) < (abs(summ) * eps)) then
        goto 12
      endif ! end if statement
    enddo
    print *, "bessK series failed to converge"
12    continue
    bessKmu = summ
    bessK1  = summ1 * xi2
  else ! evaluate the second continued fraction by Steed's algorithm;
     ! note that there can be no zero denominators
    b    = 2.d0 * (1.d0 + x)
    d    = 1.d0 / b
    delh = d
    h    = delh
    q1   = 0.d0       ! initializing for recurrences
    q2   = 1.d0
    a1   = 0.25d0 - xmu2
    c    = a1
    q    = c
    a    = -a1
    s    = 1.d0 + q * delh
    do i = 2, maxit
      a    = a - 2 * (i - 1)
      c    = -a * c / i
      qnew = (q1 - b * q2) / a
      q1   = q2
      q2   = qnew
      q    = q + c * qnew
      b    = b + 2.d0
      d    = 1.d0 / (b + a * d)
      delh = (b * d - 1.d0) * delh
      h    = h + delh
      dels = q * delh
      s    = s + dels
      ! Here because the second continued fraction converges
      ! faster, the convergence of the sum is the only test
      ! needed:
      if (abs(dels / s) < eps) then
        goto 13
      endif ! end if statement
    enddo
    print *, "second continued fraction failed to converge in dbessIK"
13    continue
    h       = a1 * h
    bessKmu = dsqrt(pi / (2.d0 * x)) * dexp(-x) / s
      ! the exp(-x) factor has been omitted in order to rescale
      ! the modified Bessel functions by exp(x) for x >= xmin
    bessK1  = bessKmu * (xmu + x + 0.5d0 - h) * xi
  endif ! end if statement
  bessKmup = xmu * xi * bessKmu - bessK1
  bessImu  = xi  / (f * bessKmu - bessKmup)
    ! I_order is determined from the Wronskian
  bessI    = (bessImu * bessIl1) / bessIl
    ! rescaling of I_order and I_order'
  bessIder = (bessImu * bessIder1) / bessIl

  ! Here are the upward recurrences of K_order
  do i = 1, nl
    bessKtemp = (xmu + i) * xi2 * bessK1 + bessKmu
    bessKmu   = bessK1
    bessK1    = bessKtemp
  enddo
  bessK    = bessKmu
  bessKder = order * xi * bessKmu - bessK1

  return
end subroutine dbessIK
!
!---------------------------------------------------------------------------------
!
! This function evaluates the Chebyshev polynomial
! sum_(k = 1)^(m) c_k * T_(k - 1)(y) - c_1 / 2 at the point
! y = [x - (b + a) / 2] / [(b - a) / 2], where c(m) is an array of
! the Chebyshev coefficients.
! Subroutine Inputs/Output:
! a - (input, real, [ ]) lower endpoint of the interval in which
!     the Chebyshev polynomial is evaluated
! b - (input, real, [ ]) upper endpoint of the interval in which the
!     Chebyshev polynomial is evaluated
! c - (input, real vector (m), [ ]) array of the Chebyshev
!       coefficients
! m - (input, integer, [ ]) number of Chebyshev coefficients to be
!     calculated
! x - (input, real, [ ]) the point that decides where the Chebyshev
!     polynomial is evaluated
!
function chebev(a, b, c, m, x)

  integer, intent(in) :: m
  real(8), intent(in) :: a, b, x, c(m)
  real(8) :: chebev
  real(8) :: d, dd, sv, y, y2
  integer :: j

  if (((x - a) * (x - b)) > 0.d0) then
    print *, "x must be greater than a and b"
    print *, "x = ", x
    print *, "a = ", a, " b = ", b
  endif ! end if statement
  d  = 0.d0
  dd = 0.d0
  y  = (2.d0 * x - a - b) / (b - a)
    ! change of variable
  y2 = 2.d0 * y

  ! Here is Clenshaw's recurrence algorithm:
  do j = m, 2, -1
    sv = d
    d  = y2 * d - dd + c(j)
    dd = sv
  enddo
  chebev = y * d - dd + 0.5d0 * c(1)

  return
end function chebev
!
!---------------------------------------------------------------------------------
!
subroutine d2bessJY(x, order, d2bessJ, d2bessY)
  implicit none
  real(8), intent(in) :: x, order
  real(8), intent(out) :: d2bessJ, d2bessY
  real *8  :: bessJ_1, bessY_1, dbessJ_1, dbessY_1
  real *8  :: bessJ_2, bessY_2, dbessJ_2, dbessY_2
  if(order.lt.1.d0) then
    write(*,*) 'error from d2bessJY: order must be >=1 '
    stop
  endif
  call dbessJY(x, order-1.d0, bessJ_1, bessY_1, dbessJ_1, dbessY_1)
  call dbessJY(x, order+1.d0, bessJ_2, bessY_2, dbessJ_2, dbessY_2)
  d2bessJ = 0.5d0*(dbessJ_1-dbessJ_2)
  d2bessY = 0.5d0*(dbessY_1-dbessY_2)
end subroutine d2bessJY
!
!---------------------------------------------------------------------------------
!
subroutine d2bessIK(x, order, d2bessI, d2bessK)
  implicit none
  real(8), intent(in) :: x, order
  real(8), intent(out) :: d2bessI, d2bessK
  real *8  :: bessI_1, bessK_1, dbessI_1, dbessK_1
  real *8  :: bessI_2, bessK_2, dbessI_2, dbessK_2
  if(order.lt.1.d0) then
    write(*,*) 'error from d2bessIK: order must be >=1 '
    stop
  endif
  call dbessIK(x, order-1.d0, bessI_1, bessK_1, dbessI_1, dbessK_1)
  call dbessIK(x, order+1.d0, bessI_2, bessK_2, dbessI_2, dbessK_2)
  d2bessI =  0.5d0*(dbessI_1+dbessI_2)
  d2bessK = -0.5d0*(dbessK_1+dbessK_2)
  !write(*,*)'from d2bessIK: dbessI_1, dbessI_2 = ',dbessI_1, dbessI_2
end subroutine d2bessIK

