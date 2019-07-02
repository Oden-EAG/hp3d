c----------------------------------------------------------------------
c
c   routine name       - cartobs
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Sep 07
c
c   purpose            - routine computes matrix of transformation
c                        from cartesian coordinates to coordinates 
c                        in an observer's system
c
c   arguments          - none
c
c   required  routines - none
c
c----------------------------------------------------------------------
c
      subroutine cartobs
c
      use graphmod
#include "syscom.blk"
c
c  ...ortogonal projections of Cartesian unit vectors on the projection 
c     plane
      dimension eproj(3,3)
c
c  ...unit vectors of the 3D Cartesian system with the first two axes
c     on the projection plane and the third in the dircetion of the
c     projection vector
      dimension uvect(3,3)
c
c  ...versor of the projection vector
      dimension rnu(3)
c
c-----------------------------------------------------------------------
c
c  ...project Cartesian unit vectors on the projection plane
      rnu(1:3) = RN(1:3); call normalize(rnu)
      eproj(1:3,1:3) = 0.d0
      do i=1,3
        eproj(i,i) = 1.d0
        eproj(1:3,i) = eproj(1:3,i) - rnu(i)*rnu(1:3)
      enddo
c
c  ...select a 2D Cartesian on the projection plane
c
c  ...if not a top view...
      s = sqrt(rnu(1)**2+rnu(2)**2)
      if (s.gt..1d0) then
c
c  .....choose the projection of e_3 vector for the vertical unit vector
        uvect(1:3,2) = eproj(1:3,3)
        call cross_product(uvect(1:3,2),rnu, uvect(1:3,1))
      else
c
c  .....choose the projection of e_1 vector for the horizontal unit vector
        uvect(1:3,1) = eproj(1:3,1)
        call cross_product(rnu,uvect(1:3,1), uvect(1:3,2))
      endif
      call normalize(uvect(1:3,1))
      call normalize(uvect(1:3,2))
c
c  ...set the third unit vector
      uvect(1:3,3) = rnu(1:3)
c
c  ...compute the transformation matrix
c
c  ...loop through the observer's coordinates
      do i=1,3
c
c  .....loop through the 3D Cartesian coordinates
        do j=1,3
          RMTR(i,j) = uvect(j,i)
        enddo
      enddo
c
c
      end subroutine cartobs
