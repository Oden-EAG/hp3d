#if HP3D_USE_X11

!----------------------------------------------------------------------
!
!   routine name       - cartobs
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - routine computes matrix of transformation
!                        from cartesian coordinates to coordinates
!                        in an observer's system
!
!----------------------------------------------------------------------
!
   subroutine cartobs
!
      use graphmod
!
      implicit none
!
!  ...ortogonal projections of Cartesian unit vectors on the projection
!     plane
      real(8) :: eproj(3,3)
!
!  ...unit vectors of the 3D Cartesian system with the first two axes
!     on the projection plane and the third in the dircetion of the
!     projection vector
      real(8) :: uvect(3,3)
!
!  ...versor of the projection vector
      real(8) :: rnu(3)
!
      integer :: i,j
      real(8) :: s
!
!-----------------------------------------------------------------------
!
!  ...project Cartesian unit vectors on the projection plane
      rnu(1:3) = RN(1:3); call normalize(rnu)
      eproj(1:3,1:3) = 0.d0
      do i=1,3
        eproj(i,i) = 1.d0
        eproj(1:3,i) = eproj(1:3,i) - rnu(i)*rnu(1:3)
      enddo
!
!  ...select a 2D Cartesian on the projection plane
!
!  ...if not a top view...
      s = sqrt(rnu(1)**2+rnu(2)**2)
      if (s.gt..1d0) then
!
!  .....choose the projection of e_3 vector for the vertical unit vector
        uvect(1:3,2) = eproj(1:3,3)
        call cross_product(uvect(1:3,2),rnu, uvect(1:3,1))
      else
!
!  .....choose the projection of e_1 vector for the horizontal unit vector
        uvect(1:3,1) = eproj(1:3,1)
        call cross_product(rnu,uvect(1:3,1), uvect(1:3,2))
      endif
      call normalize(uvect(1:3,1))
      call normalize(uvect(1:3,2))
!
!  ...set the third unit vector
      uvect(1:3,3) = rnu(1:3)
!
!  ...compute the transformation matrix
!
!  ...loop through the observer's coordinates
      do i=1,3
!
!  .....loop through the 3D Cartesian coordinates
        do j=1,3
          RMTR(i,j) = uvect(j,i)
        enddo
      enddo
!
!
   end subroutine cartobs

#endif
