!> @brief Performs H(curl) Piola transform
!> @date  Apr 2024
subroutine piola_hcurl(Dxidx,Nrdof,ShapE, pShapE)
!
   implicit none
!
   real(8), intent(in)  :: Dxidx(3,3)
   integer, intent(in)  :: Nrdof
   real(8), intent(in)  :: ShapE(3,Nrdof)
   real(8), intent(out) :: pShapE(3,Nrdof)
!
#if HP3D_DEBUG
   real(8) :: xShapE(3,Nrdof)
   integer :: i,k
   do k = 1,Nrdof
      xShapE(1:3,k) = ShapE(1,k)*Dxidx(1,1:3) &
                    + ShapE(2,k)*Dxidx(2,1:3) &
                    + ShapE(3,k)*Dxidx(3,1:3)
   enddo
#endif
!
   call DGEMM('T','N',3,Nrdof,3,1.0d0, Dxidx (1:3,1:3    ),3,  &
                                       ShapE (1:3,1:Nrdof),3,  &
                                0.0d0,pShapE (1:3,1:Nrdof),3)
!
!..TODO: add piola unit test instead
#if HP3D_DEBUG
   do k = 1,Nrdof
      do i = 1,3
         if (abs(pShapE(i,k)-xShapE(i,k)) > 1.0d-12) then
            write(*,*) 'soleval: pShapE,xShapE = ',pShapE(i,k),xShapE(i,k)
            stop
         endif
      enddo
   enddo
#endif
!
end subroutine piola_hcurl

!> @brief Performs H(curl) Piola transform
!!        with transposed output matrix
!> @date  Apr 2024
subroutine piola_hcurl_trans(Dxidx,Nrdof,ShapE, pShapE)
!
   implicit none
!
   real(8), intent(in)  :: Dxidx(3,3)
   integer, intent(in)  :: Nrdof
   real(8), intent(in)  :: ShapE(3,Nrdof)
   real(8), intent(out) :: pShapE(Nrdof,3)
!
#if HP3D_DEBUG
   real(8) :: xShapE(3,Nrdof)
   integer :: i,k
   do k = 1,Nrdof
      xShapE(1:3,k) = ShapE(1,k)*Dxidx(1,1:3) &
                    + ShapE(2,k)*Dxidx(2,1:3) &
                    + ShapE(3,k)*Dxidx(3,1:3)
   enddo
#endif
!
   call DGEMM('T','N',Nrdof,3,3,1.0d0, ShapE (1:3,1:Nrdof),3,  &
                                       Dxidx (1:3,1:3    ),3,  &
                                0.0d0,pShapE (1:Nrdof,1:3),Nrdof)
!
#if HP3D_DEBUG
   do k = 1,Nrdof
      do i = 1,3
         if (abs(pShapE(k,i)-xShapE(i,k)) > 1.0d-12) then
            write(*,*) 'soleval: pShapE,xShapE = ',pShapE(k,i),xShapE(i,k)
            stop
         endif
      enddo
   enddo
#endif
!
end subroutine piola_hcurl_trans


!> @brief Performs H(div) Piola transform
!> @date  Apr 2024
subroutine piola_hdiv(Dxdxi,Rjac,Nrdof,ShapV, pShapV)
!
   implicit none
!
   real(8), intent(in)  :: Dxdxi(3,3)
   real(8), intent(in)  :: Rjac
   integer, intent(in)  :: Nrdof
   real(8), intent(in)  :: ShapV(3,Nrdof)
   real(8), intent(out) :: pShapV(3,Nrdof)
!
   real(8) :: iRjac
!
#if HP3D_DEBUG
   real(8) :: xShapV(3,Nrdof)
   integer :: i,k
   do k = 1,Nrdof
      xShapV(1:3,k) = Dxdxi(1:3,1)*ShapV(1,k) &
                    + Dxdxi(1:3,2)*ShapV(2,k) &
                    + Dxdxi(1:3,3)*ShapV(3,k)
      xShapV(1:3,k) = xShapV(1:3,k) / Rjac
   enddo
#endif
!
   iRjac = 1.d0 / Rjac
   call DGEMM('N','N',3,Nrdof,3,iRjac,Dxdxi (1:3,1:3    ),3,  &
                                      ShapV (1:3,1:Nrdof),3,  &
                                0.0d0,pShapV(1:3,1:Nrdof),3)
!
#if HP3D_DEBUG
   do k = 1,Nrdof
      do i = 1,3
         if (abs(pShapV(i,k)-xShapV(i,k)) > 1.0d-12) then
            write(*,*) 'soleval: pShapV,xShapV = ',pShapV(i,k),xShapV(i,k)
            stop
         endif
      enddo
   enddo
#endif
!
end subroutine piola_hdiv

!> @brief Performs H(div) Piola transform
!!        with transposed output matrix
!> @date  Apr 2024
subroutine piola_hdiv_trans(Dxdxi,Rjac,Nrdof,ShapV, pShapV)
!
   implicit none
!
   real(8), intent(in)  :: Dxdxi(3,3)
   real(8), intent(in)  :: Rjac
   integer, intent(in)  :: Nrdof
   real(8), intent(in)  :: ShapV(3,Nrdof)
   real(8), intent(out) :: pShapV(Nrdof,3)
!
   real(8) :: iRjac
!
#if HP3D_DEBUG
   real(8) :: xShapV(3,Nrdof)
   integer :: i,k
   do k = 1,Nrdof
      xShapV(1:3,k) = Dxdxi(1:3,1)*ShapV(1,k) &
                    + Dxdxi(1:3,2)*ShapV(2,k) &
                    + Dxdxi(1:3,3)*ShapV(3,k)
      xShapV(1:3,k) = xShapV(1:3,k) / Rjac
   enddo
#endif
!
   iRjac = 1.d0 / Rjac
   call DGEMM('T','T',3,Nrdof,3,iRjac,ShapV (1:3,1:Nrdof),3,  &
                                      Dxdxi (1:3,1:3    ),3,  &
                                0.0d0,pShapV(1:Nrdof,1:3),3)
!
#if HP3D_DEBUG
   do k = 1,Nrdof
      do i = 1,3
         if (abs(pShapV(k,i)-xShapV(i,k)) > 1.0d-12) then
            write(*,*) 'soleval: pShapV,xShapV = ',pShapV(k,i),xShapV(i,k)
            stop
         endif
      enddo
   enddo
#endif
!
end subroutine piola_hdiv_trans
