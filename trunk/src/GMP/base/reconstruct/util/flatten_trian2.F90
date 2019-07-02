!--------------------------------------------------------------------------
subroutine flatten_trian2(Nt)
!
      use GMP , only : TRIANGLES
!      
      implicit none
      integer,intent(in)  :: Nt
      real*8,dimension(3) :: temp,v1,v2,v3
      integer             :: i,j
      real*8              :: x,y
      integer, parameter  :: deg=7
      integer, external   :: bijec
!
!--------------------------------------------------------------------------
!  STEP 1 : modify interior control points of triangle                    |
!--------------------------------------------------------------------------
!
      call trian2verts(Nt, v1,v2,v3)
      do j=1,deg-1
        do i=1,(deg-1-j)
          x=float(i)/float(deg)
          y=float(j)/float(deg)
          temp(1:3)=(1.d0-x-y)*v1(1:3) + x*v2(1:3) + y*v3(1:3)
          TRIANGLES(Nt)%Rdata(bijec(i,j):bijec(i,j)+2)=temp(1:3)
        enddo
      enddo
!
!
end subroutine flatten_trian2
!--------------------------------------------------------------------------
