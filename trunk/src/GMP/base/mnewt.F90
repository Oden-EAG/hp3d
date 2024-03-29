!-----------------------------------------------------------------------
!> @brief routine uses Newton-Raphson iterations to solve a nonlinear
!            system of equations defining an implicit parametrization
!            for a GMP point, curve, triangle or rectangle.
!
!> @param[in ]  Nusf    = 1  implicit point
!                       = 2  implicit curve
!                       = 3  implicit triangle
!                       = 4  implicit rectangle
!> @param[in ]  Nsurf   - surface numbers
!> @param[in ]  Eta     - reference coordinates of the point being
!                         determined
!> @param[in ]  Feta    - values of the stretching functions (implicit
!                         triangle or rectangle only)
!> @param[in ]  Xs      - a starting point
!> @param[in ]  Sfact   - renormalization factors for the surfaces
!> @param[out]  X       - physical coordinates of the point
!
!> @date Mar 2023
!-----------------------------------------------------------------------
!
   subroutine mnewt(Nusf,Nsurf,Eta,Feta,Xs,Sfact, X)
!
      implicit none
!
      integer, intent(in)  :: Nusf,Nsurf(6)
      real(8), intent(in)  :: Eta(*),Feta(4),Xs(3),Sfact(4)
      real(8), intent(out) :: X(3)
!
!  ...stiffness matrix and load vector for NR iterations
      real(8) :: alpha(3,3),beta(3)
      real(8) :: aux(3)
      real(8), parameter :: eps = 1.d-12
!
      integer :: i,ifl,j,k,ntrial
      real(8) :: errf,errx
!
      integer :: iprint
      iprint=0
!
!-----------------------------------------------------------------------
!
 10   continue
!
!  ...maximum number of iterations
      ntrial = 10
!
!  ...printing
      if (iprint.eq.1) then
        write(*,1000)Nusf,Nsurf
 1000   format(' Nusf,Nsurf = ',7(i1,2x))
        write(*,1003)Feta,Xs,Sfact
 1003   format(' Feta,Xs,Sfact = ',11(e12.5,2x))
        select case(Nusf)
        case(2)
          write(*,1001) Eta(1)
 1001     format(' Eta = ',e12.5)
        case(3,4)
          write(*,1002) Eta(1:2)
 1002     format(' Eta = ',2(e12.5,2x))
        endselect
      endif
!
!  ...set X equal to initial guess Xs
      X(1:3)=Xs(1:3)
!
!  ...loop through N-R iterations
      do k=1,ntrial
!
!       N-R iteration:  x_n+1 = x_n - DF(x_n)^-T * F(x_n)
!
!  .....construct the N-R system:
        call usrfun(Nusf,Nsurf,Eta,Feta,X,Sfact, alpha,beta)
!       alpha = DF^T (still needs to be inverted!)
!       beta  = -F
!
!  .....printing statements
        if (iprint .eq. 1) then
          select case(Nusf)
!
!  .......implicit curve
          case(2)
            write(*,6001) (Nsurf(i),i=1,4)
 6001       format(' mnewt: IMPLICIT CURVE ', &
                 /,'        SURFACES        = ',4i6)
            write(*,6002) Eta(1),Xs(1:3)
 6002       format('       Eta, Xs = ',f8.3,2x,3f8.3)
!
!  .......implicit rectangle
          case(4)
            write(*,9001) (Nsurf(i),i=1,5),(Feta(i),i=1,4)
 9001       format(' mnewt: IMPLICIT RECTANGLE ', &
                 /,'        SURFACES          = ',5i6, &
                 /,'        STRETCH FNCTS Feta = ',4f8.3)
            write(*,5001) (Eta(i),i=1,2)
          end select
 5001     format(' mnewt: Eta         = ',3f8.3)
          write(*,5002) k,X
 5002     format(' mnewt: k,X         = ',i2,2x,3e12.5)
          write(*,5003)
 5003     format(' mnewt: beta, alpha = ')
          do i = 1, 3
            write(*,5004) beta(i),(alpha(i,j),j=1,3)
 5004       format(e12.5,4x,3(e12.5,2x))
          enddo
          call pause
!  .....end of printing statements
        endif
!
!  .....compute L^1 norm of beta
        errf = 0.d0
        do i = 1, 3
          errf = errf + abs(beta(i))
        enddo
!  .....if beta is small, then return
        if (errf .le. eps) return
!
!  .....solve linear system:  alpha*aux = beta  ==>  aux = -DF^-T * F
        call saruss(alpha,beta, aux,ifl)
!  .....if system is not solved
        if (ifl .ne. 0) then
          iprint = 1
          goto 10
        endif
!
        beta(1) = aux(1)
        beta(2) = aux(2)
        beta(3) = aux(3)
!
!  .....compute L^1 norm of beta and update X
        errx = 0.d0
        do i = 1, 3
          errx = errx + abs(beta(i))
          X(i) = X(i) + beta(i)
        enddo
!  .....if beta is small, then return
        if (errx .le. eps) return
!
!  .....if maximum number of iterations has been reached
        if (k .eq. ntrial) then
          write(*,*)'mnewt: NEWTON-RAPSON METHOD HAS NOT CONVERGED!'
          call pause
          iprint = 1
          goto 10
        endif
!
!  ...end of loop through iterations
      enddo
!
   end subroutine mnewt
!
!
!
!-----------------------------------------------------------------------
!> @brief routine formulates a system of linear equations for the
!            NR solution of a system of algebraic equations defining
!            an implicit point, curve, triangle or rectangle
!
!
!> @param[in ]  Nusf    = 1  implicit point
!                       = 2  implicit curve
!                       = 3  implicit triangle
!                       = 4  implicit rectangle
!> @param[in ]  Nsurf   - surface numbers
!> @param[in ]  Eta     - reference coordinates of the point being
!                         determined
!> @param[in ]  Feta    - values of the stretching functions (implicit
!                         triangle or rectangle only)
!> @param[in ]  X       - current iterate (approximation of the
!                         physical coordinates)
!> @param[in ]  Sfact   - renormalization factors for the surfaces
!> @param[out]  Alpha   - 3 x 3 stiffness matrix
!> @param[out]  Beta    - load vector
!
!> @date Mar 2023
!-----------------------------------------------------------------------
!
   subroutine usrfun(Nusf,Nsurf,Eta,Feta,X,Sfact, Alpha,Beta)
!
      implicit none
!
      real(8) :: ctr50(3)
!      common /sphepoint/ ctr50
!
      integer :: Nusf,Nsurf(6)
      real(8) :: Eta(*),Feta(4),X(3),Sfact(4)
      real(8) :: Alpha(3,3),Beta(3)
!
!  ...surface gradients for different surfaces
      real(8) :: der1(3),der2(3),der3(3),der4(3),der5(3)
      real(8) :: dtt(3),dtt2(3)
!
      integer :: i
      real(8) :: fval,fval1,fval2,fval3,fval4,fval5,tt,tt2,xf1,xf2
!
      integer :: iprint
      iprint=0
!
!  ...select geometry entity
      select case(Nusf)
!
!----------------------------------------------------------------------
!  ...implicit point: intersection of three surfaces
      case(1)
!  .....compute F(x)
        call surf(Nsurf(1),X, Fval1,der1)
        Beta(1) = -Fval1
        call surf(Nsurf(2),X, Fval2,der2)
        Beta(2) = -Fval2
        call surf(Nsurf(3),X, Fval3,der3)
        Beta(3) = -Fval3
!
!  .....compute DF^T(x)
        do i = 1, 3
          Alpha(1,i) = der1(i)
          Alpha(2,i) = der2(i)
          Alpha(3,i) = der3(i)
        enddo
!
!----------------------------------------------------------------------
!  ...implicit curve
      case(2)
      if (iprint.eq.1) then
        write(*,7004) Eta(1), X(1:3)
 7004   format(' usrfun: IMPLICIT CURVE Eta = ',f8.3, &
               ' ITERATE X = ', 3f10.5)
      endif
!
      call surf(Nsurf(1),X, Fval,der1)
      Beta(1) = -Fval
!
      call surf(Nsurf(2),X, Fval,der2)
      Beta(2) = -Fval
!
      call surf(Nsurf(3),X, Fval3,der3)
      call surf(Nsurf(4),X, Fval4,der4)
      fval3 = fval3*Sfact(1); der3(1:3) = der3(1:3)*Sfact(1)
      fval4 = fval4*Sfact(2); der4(1:3) = der4(1:3)*Sfact(2)
!
      Beta(3) = -(1.d0-Eta(1))*Fval3 - Eta(1)*Fval4
      do  i=1,3
        Alpha(1,i) = der1(i)
        Alpha(2,i) = der2(i)
        Alpha(3,i) = (1.d0-Eta(1))*der3(i) + Eta(1)*der4(i)
      enddo
!
!----------------------------------------------------------------------
!  ...implicit triangle
      case(3)
      write(*,*) 'usrfun: NOT TESTED'
!
      call surf(Nsurf(1),X, fval1,der1)
      Beta(1) = -Fval1
!
      call surf(Nsurf(2),X, fval2,der2)
      call surf(Nsurf(4),X, fval4,der4)
      fval2=fval2*Sfact(1) ; der5(1:3)=der2(1:3)*Sfact(1)
      fval4=fval4*Sfact(3) ; der5(1:3)=der4(1:3)*Sfact(3)
      xf2 = 1.d0 - Feta(3)
!
      Beta(2) = -xf2*fval2 - Feta(3)*fval4
!
      do i=1,3
        Alpha(1,i) = der1(i)
        Alpha(2,i) = xf2*der2(i) + Feta(3)*der4(i)
      enddo

      call surf(Nsurf(3),X, fval3,der3)
      fval3=fval3*Sfact(2) ; der3(1:3)=der3(1:3)*Sfact(2)
      call sphere1(X,ctr50,0.d0, fval5,der5)
      xf1 = 1.d0 - Feta(1)
      xf2 = 1.d0 - Feta(4)
!
      Beta(3) = -Eta(1)*(xf1*fval5+Feta(1)*fval3) &
              -  Eta(2)*(xf2*fval5+Feta(4)*fval3)

      do i=1,3
        Alpha(3,i) = Eta(1)*(xf1*der5(i)+Feta(1)*der3(i)) &
                  +  Eta(2)*(xf2*der5(i)+Feta(4)*der3(i))
      enddo
      if (iprint.eq.1) then
         write(*,*)'usrfun:'
         write(*,1000) ctr50
 1000    format(' sph center = ',3(e12.5,2x))
         write(*,*)'der1',der1
         write(*,*)'der2,der4',der2,der4
         write(*,*)'der3,der5',der3,der5
      endif
      return
!
!----------------------------------------------------------------------
!  ...implicit rectangle
      case(4)
!
      call surf(Nsurf(1),X, fval1,der1)
      Beta(1) = -fval1
      do i=1,3
        Alpha(1,i) = der1(i)
      enddo
!
      call surf(Nsurf(5),X, fval5,der5)
      call surf(Nsurf(3),X, fval3,der3)
      fval5 = fval5*Sfact(4); der5(1:3) = der5(1:3)*Sfact(4)
      fval3 = fval3*Sfact(2); der3(1:3) = der3(1:3)*Sfact(2)
!
      xf1 = 1.d0 - Feta(1)
      xf2 = 1.d0 - Feta(3)
!
      tt        =  xf1*fval5     + Feta(1)*fval3
      dtt(1:3)  =  xf1*der5(1:3) + Feta(1)*der3(1:3)
      tt2       =  xf2*fval5     + Feta(3)*fval3
      dtt2(1:3) =  xf2*der5(1:3) + Feta(3)*der3(1:3)
!
      Beta(2) = -(1.d0-Eta(2))*tt -  Eta(2)*tt2
      do i=1,3
        Alpha(2,i) = (1.d0-Eta(2))*dtt(i) &
                   +  Eta(2)*dtt2(i)
      enddo
!
!
      call surf(Nsurf(2),X, fval2,der2)
      call surf(Nsurf(4),X, fval4,der4)
      fval2 = fval2*Sfact(1); der2(1:3) = der2(1:3)*Sfact(1)
      fval4 = fval4*Sfact(3); der4(1:3) = der4(1:3)*Sfact(3)
!
      xf1 = 1.d0 - Feta(4)
      xf2 = 1.d0 - Feta(2)
!
      tt        =  (xf1*fval2     + Feta(4)*fval4)
      dtt(1:3)  =  (xf1*der2(1:3) + Feta(4)*der4(1:3))
!
      tt2       = (xf2*fval2     + Feta(2)*fval4)
      dtt2(1:3) = (xf2*der2(1:3) + Feta(2)*der4(1:3))
!
!
      Beta(3) = -(1.d0-Eta(1))*tt &
              -  Eta(1)*tt2

      do  i=1,3
        Alpha(3,i) = (1.d0-Eta(1))*dtt(i) &
                   +  Eta(1)*dtt2(i)
      enddo
!
      if (iprint.eq.1) then
         write(*,7002) fval1,fval2,fval3,fval4,fval5
 7002    format(' usrfun: fval1,...,5 = ',5(e12.5,2x))
      endif
!
      case default
        write(*,7001) Nusf
 7001   format('usrfun: WRONG INPUT Nusf = ',i5)
        stop
!
!  ...end select geometry entity
      endselect
!
!--------------------------------------------------------------------
!  ...printing statement
      if (iprint .eq. 1) then
        write(*,*) 'usrfun: Beta, Alpha = '
        do i=1,3
          write(*,7003) Beta(i), Alpha(i,1:3)
 7003     format(e12.5,4x,3(e12.5,2x))
        enddo
      endif
!
   end subroutine usrfun
