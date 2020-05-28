c
c----------------------------------------------------------------------
c
c   routine name       - mat3prod
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - May 2020
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
c   required  routines - none
c
c----------------------------------------------------------------------
c
      subroutine mat3prod(Ma,Mb,Mc)
c
      implicit real(8) (a-h,m-z)
c
      dimension Ma(3,3), Mb(3,3), Mc(3,3)
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

