#if HP3D_USE_X11

!
!----------------------------------------------------------------------
!
!   routine name       - mat3prod
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - routine computes the product of two given
!                        matrices 3x3
!
!   usage              - call mat3prod(Ma,Mb,Mc)
!
!   arguments
!             in       - Ma, Mb - matrices 3x3
!             out      - Mc - their product
!
!----------------------------------------------------------------------
!
      subroutine mat3prod(Ma,Mb,Mc)
!
      implicit none
!
      real(8), intent(in)  :: Ma(3,3), Mb(3,3)
      real(8), intent(out) :: Mc(3,3)
!
      integer :: i,j,k
      real(8) :: pr
!
      do i=1,3
         do j=1,3
            pr = 0.d0
            do k=1,3
               pr = pr + Ma(i,k)*Mb(k,j)
            enddo
            Mc(i,j) = pr
         enddo
      enddo
!
      end subroutine

#endif
