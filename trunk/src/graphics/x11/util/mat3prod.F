#if HP3D_USE_X11

c
c----------------------------------------------------------------------
c
c   routine name       - mat3prod
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 2024
c
c   purpose            - routine computes the product of two given
c                        matrices 3x3
c
c   usage              - call mat3prod(Ma,Mb,Mc)
c
c   arguments
c             in       - Ma, Mb - matrices 3x3
c             out      - Mc - their product
c
c----------------------------------------------------------------------
c
      subroutine mat3prod(Ma,Mb,Mc)
c
      implicit none
c
      real(8), intent(in)  :: Ma(3,3), Mb(3,3)
      real(8), intent(out) :: Mc(3,3)
c
      integer :: i,j,k
      real(8) :: pr
c
      do i=1,3
         do j=1,3
            pr = 0.d0
            do k=1,3
               pr = pr + Ma(i,k)*Mb(k,j)
            enddo
            Mc(i,j) = pr
         enddo
      enddo
c
      end subroutine

#endif
