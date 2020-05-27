!------------------------------------------------------------------------------------
!> Purpose : plane parameterization
!!
!! @param[in]    X    - physical coordinates of a point
!! @param[in]    X0   - a point on the plane
!! @param[inout] Rn   - normal to the plane
!! @param[out]   Fval - function value
!! @param[out]   Dfdx - function gradient
!!
!! @revision Nov 12
!------------------------------------------------------------------------------------
!
subroutine plane1(X,X0,Rn, Fval,Dfdx)
!
      implicit none
      real(8), dimension(3), intent(in )   :: X,X0
      real(8), dimension(3), intent(inout) :: Rn
      real(8)              , intent(out)   :: Fval
      real(8), dimension(3), intent(out)   :: Dfdx
!------------------------------------------------------------------------------------
      integer :: i
      real(8)  :: dRn,dRh,d
!------------------------------------------------------------------------------------
!
      dRn=0.d0
      do 10 i=1,3
        dRh = Rn(i)*Rn(i)
        dRn = dRn+dRh
 10   continue
!
      dRn = dsqrt(dRn)
      do 20 i=1,3
        Rn(i) = Rn(i)/dRn
 20   continue
!
!  ...caculate coefficients of the plane
      d = -( X0(1)*Rn(1) + X0(2)*Rn(2) + X0(3)*Rn(3) )
!
      Fval =  X(1)*Rn(1) + X(2)*Rn(2) + X(3)*Rn(3) + d
!
      Dfdx(1) = Rn(1)
      Dfdx(2) = Rn(2)
      Dfdx(3) = Rn(3)
!
!
end subroutine plane1
!
!
!------------------------------------------------------------------------------------
!> Purpose : plane parameterization
!!
!! @param[in]  X       - physical coordinates of a point
!! @param[in]  X1,2,3  - 3 points defining the plane
!! @param[out] Fval    - function value
!! @param[out] Dfdx    - function gradient
!!
!! @revision Nov 12
!------------------------------------------------------------------------------------
!
subroutine plane2(X,X1,X2,X3, Fval,Dfdx)
!
      implicit none
      real(8), dimension(3), intent(in ) :: X,X1,X2,X3
      real(8)              , intent(out) :: Fval
      real(8), dimension(3), intent(out) :: Dfdx
!------------------------------------------------------------------------------------
!  ...vectors X1X2, X1X3, and the normal versor
      real(8), dimension(3) :: ver1,ver2,rn
      real(8)               :: s,d
      integer              :: iprint
!------------------------------------------------------------------------------------
!
      iprint=0
      if (iprint.eq.1) then
        write(*,7001) X,X1,X2,X3
 7001   format(' plane2: X,X1,X2,X3 = ',4(3f8.3,2x))
        call pause
      endif
!
!  ...caculate the versor orthogonal to the plane
      ver1(1:3) = X2(1:3) - X1(1:3)
      ver2(1:3) = X3(1:3) - X1(1:3)
      call cross_product(ver1,ver2, rn)
      s = sqrt(rn(1)**2+rn(2)**2+rn(3)**2)
      rn(1:3) = rn(1:3)/s
!
      d = -( X1(1)*rn(1) + X1(2)*rn(2) + X1(3)*rn(3) )
      Fval =  X(1)*rn(1) + X(2)*rn(2) + X(3)*rn(3) + d
!
      Dfdx(1) = rn(1)
      Dfdx(2) = rn(2)
      Dfdx(3) = rn(3)
!
!
end subroutine plane2
!
!
!
!------------------------------------------------------------------------------------
!> Purpose : routine determines a unit vector orthogonal to a plane that passes
!!           through a point X and is at a minimum distance from points Y1,2,3
!!
!! @param[in ] X,Y  - coordinates of points
!! @param[out] Rn   - unit normal to the plane
!!
!! @revision Nov 12
!------------------------------------------------------------------------------------
subroutine determine_plane1(X,Y, Rn)
!
      implicit none
      real(8),dimension(3  ),intent(in ) :: X
      real(8),dimension(3,3),intent(in ) :: Y
      real(8),dimension(3  ),intent(out) :: Rn
!
!  ...auxiliary matrices
      real(8),dimension(3,3) :: a,b
!
!  ...eigenvalues and work space
      real(8),dimension(3)   :: w
      real(8),dimension(102) :: work
!
      real(8) :: s,rlambda
      integer :: i,j,k,iprint,info
!------------------------------------------------------------------------------------
!
      iprint=0
!
!  ...determine the auxiliary matrices
      do i=1,3
        do k=1,3
          a(i,k) = X(i) - Y(i,k)
        enddo
      enddo
      do i=1,3
        do j=1,3
          s = 0.d0
          do k=1,3
            s = s + a(i,k)*a(j,k)
          enddo
          b(i,j) = s
        enddo
      enddo
!
      if (iprint.eq.1) then
        write(*,7001)
 7001   format('determine_plane1: b = ')
        do i=1,3
          write(*,7002) (b(i,j),j=1,3)
 7002     format(3e12.5)
        enddo
        call pause
      endif
!
!  ...solve the eigenvalue problem and determine the eigenvector
!     corresponding to the smallest eigenvalue
      call dsyev('V','L',3,b,3,w,work,102,info)
      rlambda = w(1)
      Rn(1:3) = b(1:3,1)
!
      if (iprint.eq.1) then
        write(*,7003) rlambda,RN(1:3)
 7003   format('determine_plane1: rlambda, Rn = ',e12.5,3x,3e12.5)
        call pause
      endif
!
!
end subroutine determine_plane1
!
!
!------------------------------------------------------------------------------------
!> Purpose : routine determines a unit vector orthogonal to a plane that passes
!!           through points X1,2 and is at a minimum distance from points Y1,2
!!
!! @param[in ] X,Y  - coordinates of points
!! @param[out] Rn   - unit normal to the plane
!!
!! @revision Nov 12
!------------------------------------------------------------------------------------
subroutine determine_plane2(X,Y, Rn)
!
      implicit none
      real(8),dimension(3,2),intent(in ) :: X
      real(8),dimension(3,2),intent(in ) :: Y
      real(8),dimension(3  ),intent(out) :: Rn
!
!  ...basis of normal vectors
      real(8),dimension(3,2) :: rna

!  ...auxiliary matrices
      real(8),dimension(2)   :: x0,x1
      real(8),dimension(3)   :: x12,y12,c,d
      real(8),dimension(2,2) :: a,b
!
!  ...eigenvalues and work space
      real(8),dimension(2)  :: w
      real(8),dimension(68) :: work
!
      real(8) :: s,rlambda
      integer :: i,j,k,iprint,info,l
!------------------------------------------------------------------------------------
!
      iprint=0
      if (iprint.eq.1) then
        write(*,7007) (X(1:3,i),i=1,2),(Y(1:3,i),i=1,2)
 7007   format('determine_plane2: ',/,'POINTS X = ',2(3e12.5,2x), &
               /,'POINTS Y = ',2(3e12.5,2x))
        call pause
      endif
!
!  ...midpoint of line X_1,X_2
      x12(1:3) = (X(1:3,1) + X(1:3,2))/2.d0
!
!  ...midpoint of line Y_1,Y_2
      y12(1:3) = (Y(1:3,1) + Y(1:3,2))/2.d0
!
!  ...vector X_1 X_2
      c(1:3) = X(1:3,2) - X(1:3,1)
!
!  ...vector X_12 Y_12
      d(1:3) = y12(1:3) - x12(1:3)
!
!  ...first basis vector
      rna(1,1) = c(2)*d(3) - c(3)*d(2)
      rna(2,1) = c(3)*d(1) - c(1)*d(3)
      rna(3,1) = c(1)*d(2) - c(2)*d(1)
!
!  ...second basis vector
      rna(1,2) = c(2)*rna(3,1) - c(3)*rna(2,1)
      rna(2,2) = c(3)*rna(1,1) - c(1)*rna(3,1)
      rna(3,2) = c(1)*rna(2,1) - c(2)*rna(1,1)
!
!  ...normalize the basis vectors
      do j=1,2
        s = 0.d0
        do i=1,3
          s = s + rna(i,j)**2
        enddo
        s = sqrt(s)
        rna(1:3,j) = rna(1:3,j)/s
      enddo
!
      if (iprint.eq.1) then
        write(*,7005) rna(1:3,1),rna(1:3,2)
 7005   format('determine_plane2: UNIT VECTORS = ',2(3e12.5,2x))
      endif
!
!  ...determine the auxiliary matrices
      do i=1,2
        do k=1,2
          s = 0.d0
          do l=1,3
            s = s + rna(l,i)*(x12(l) - Y(l,k))
          enddo
          a(i,k) = s
        enddo
      enddo
!
      if (iprint.eq.1) then
        write(*,7006)
 7006   format('determine_plane2: a = ')
        do i=1,2
          write(*,7002) (a(i,j),j=1,2)
        enddo
      endif
!
      do i=1,2
        do j=1,2
          s = 0.d0
          do k=1,2
            s = s + a(i,k)*a(j,k)
          enddo
          b(i,j) = s
        enddo
      enddo
!
      if (iprint.eq.1) then
        write(*,7001)
 7001   format('determine_plane2: b = ')
        do i=1,2
          write(*,7002) (b(i,j),j=1,2)
 7002     format(2e12.5)
        enddo
      endif
!
!  ...solve the eigenvalue problem and determine the eigenvector
!     corresponding to the smallest eigenvalue
      call dsyev('V','L',2,b,2,w,work,68,info)
      if (info.ne.0) then
        write(*,*) 'determine_plane2: ERROR IN EIGENSOLVER'
        stop
      endif
      rlambda = w(1)
!wr10.05.05
!cccc x1(1:3) = b(1:3,1)
      x1(1:2) = b(1:2,1)
!
      if (iprint.eq.1) then
        write(*,7003) rlambda,x1(1:2)
 7003   format('determine_plane2: rlambda, x1 = ',e12.5,3x,2e12.5)
      endif
!
      do i=1,3
        Rn(i) = x1(1)*rna(i,1) + x1(2)*rna(i,2)
      enddo
!
      if (iprint.eq.1) then
        write(*,7004) Rn(1:3)
 7004   format('determine_plane2: Rn = ',3e12.5)
        call pause
      endif
!
!
end subroutine determine_plane2
