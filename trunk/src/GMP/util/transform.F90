!-----------------------------------------------------------------------
!> @Purpose - routine transforms components of a vector from a Cartesian
!!            to a Cartesian system of coordinates
!
!> @param[in]    Idec   = 12  transform from system 1 to system 2
!>                        21  transform from system 2 to system 1
!> @param[in]    Transf - transformation matrix from system 1 to
!>                        system 2, i.e. columns of the matrix are
!>                        versors of system 2 according to system 1
!> @param[inout] X1(X2) - components of a vector in system 1(2)
!> @param[inout] X2(X1) - components of a vector in system 2(1)
!
!> @revision Feb 13
!-----------------------------------------------------------------------
!
subroutine transform(Idec,Transf,X1,X2)
!
      implicit none
      integer, intent(in   ) :: Idec
      real(8), intent(in   ) :: Transf(3,3)
      real(8), intent(inout) :: X1(3)
      real(8), intent(inout) :: X2(3)
!
      real(8) :: s
      integer :: i,j
!-----------------------------------------------------------------------
!
      select case(Idec)
!
      case(12)
        do i=1,3
          s=0.d0
          do j=1,3 ; s=s+Transf(j,i)*X1(j) ; enddo
          X2(i)=s
        enddo
!
      case(21)
        do i=1,3
          s=0.d0
          do j=1,3 ; s=s+Transf(i,j)*X2(j) ; enddo
          X1(i)=s
        enddo
!
      case default
        write(*,*) 'transform: WRONG Idec = ',Idec
        stop
      endselect
!
!
end subroutine transform
