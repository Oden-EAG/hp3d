!-----------------------------------------------------------------------
!
!   routine name       - geom
!
!-----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - routine determines derivatives of the master
!                        element coordinates wrt to the physical
!                        coordinates and the jacobian for a 3D element
!
!   arguments :
!     in:
!               Dxdxi  - derivatives of the physical coordinates wrt
!                        master element coordinates
!                             first index is the physical coord
!                             second index is the master coord
!     out:
!               Dxidx  - derivatives of the master coord wrt to the
!                        physical ones
!                             first index is the master coordinate
!                             second index is the physical coordinate
!               Rjac   - the jacobian
!               Iflag  = 0 if the jacobian is positive
!                      = 1 otherwise
!
!-----------------------------------------------------------------------
!
   subroutine geom(Dxdxi, Dxidx,Rjac,Iflag)
!
      implicit none
!
      real(8), intent(in)  :: Dxdxi(3,3)
      real(8), intent(out) :: Dxidx(3,3)
      real(8), intent(out) :: Rjac
      integer, intent(out) :: Iflag
!
      real(8) :: det1,det2,det3
!
#if HP3D_DEBUG
      real(8) :: a(3,3),a1,a2,a3,s
      integer :: i,j,k
      integer :: iprint
      iprint = 0
      if (iprint.eq.1) then
        do i=1,3
          write(*,9999) i,Dxdxi(i,1:3)
 9999     format(' geom: i,Dxdxiz(i,1:3) = ',i2,2x,3(e12.5,2x))
        enddo
      endif
#endif
!
      Dxidx=0; Iflag=0
!
!  ...jacobian
      Rjac = Dxdxi(1,1)*Dxdxi(2,2)*Dxdxi(3,3)   &
           + Dxdxi(2,1)*Dxdxi(3,2)*Dxdxi(1,3)   &
           + Dxdxi(3,1)*Dxdxi(1,2)*Dxdxi(2,3)   &
           - Dxdxi(3,1)*Dxdxi(2,2)*Dxdxi(1,3)   &
           - Dxdxi(1,1)*Dxdxi(3,2)*Dxdxi(2,3)   &
           - Dxdxi(2,1)*Dxdxi(1,2)*Dxdxi(3,3)
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,*) 'geom: Rjac = ',Rjac
      endif
#endif
!
      if (Rjac.lt.0.d0) then
        Iflag=1
!
#if HP3D_DEBUG
        if (iprint.eq.1) then
!
!  .......compute vector product of first two derivatives
          a1 = Dxdxi(2,1)*Dxdxi(3,2) - Dxdxi(3,1)*Dxdxi(2,2)
          a2 = Dxdxi(3,1)*Dxdxi(1,2) - Dxdxi(1,1)*Dxdxi(3,2)
          a3 = Dxdxi(1,1)*Dxdxi(2,2) - Dxdxi(2,1)*Dxdxi(1,2)
          write(*,7001) Dxdxi(1:3,1), Dxdxi(1:3,2)
 7001     format('geom: Dxdxi1 = ',3e12.5,' Dxdxi2 = ',3e12.5)
          write(*,7002) a1,a2,a3,  Dxdxi(1:3,3)
 7002     format('geom: Dxdxi1 x Dxdxi2 = ',3e12.5,' Dxdxi3 = ',3e12.5)
          write(*,7003) Rjac,a1*Dxdxi(1,3)+a2*Dxdxi(2,3)+a3*Dxdxi(3,3)
 7003     format('geom: Rjac = ',e12.5,'(Dxdxi1 x Dxdxi2)oDxdxi3 = ',e12.5)
          call pause
        endif
#endif
      endif
!
      det1 =   Dxdxi(2,2)*Dxdxi(3,3) - Dxdxi(3,2)*Dxdxi(2,3)
      det2 = - Dxdxi(2,1)*Dxdxi(3,3) + Dxdxi(3,1)*Dxdxi(2,3)
      det3 =   Dxdxi(2,1)*Dxdxi(3,2) - Dxdxi(3,1)*Dxdxi(2,2)
!
      Dxidx(1,1) = det1/Rjac
      Dxidx(2,1) = det2/Rjac
      Dxidx(3,1) = det3/Rjac
!
      det1 =   Dxdxi(3,2)*Dxdxi(1,3) - Dxdxi(1,2)*Dxdxi(3,3)
      det2 =   Dxdxi(1,1)*Dxdxi(3,3) - Dxdxi(3,1)*Dxdxi(1,3)
      det3 = - Dxdxi(1,1)*Dxdxi(3,2) + Dxdxi(3,1)*Dxdxi(1,2)
!
      Dxidx(1,2) = det1/Rjac
      Dxidx(2,2) = det2/Rjac
      Dxidx(3,2) = det3/Rjac
!
      det1 =   Dxdxi(1,2)*Dxdxi(2,3) - Dxdxi(2,2)*Dxdxi(1,3)
      det2 = - Dxdxi(1,1)*Dxdxi(2,3) + Dxdxi(2,1)*Dxdxi(1,3)
      det3 =   Dxdxi(1,1)*Dxdxi(2,2) - Dxdxi(2,1)*Dxdxi(1,2)
!
      Dxidx(1,3) = det1/Rjac
      Dxidx(2,3) = det2/Rjac
      Dxidx(3,3) = det3/Rjac
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
!
!  .....check the inversion
        a(1:3,1:3) = 0.d0
        do i=1,3
          do j=1,3
            s=0.d0
            do k=1,3
              s = s + Dxdxi(i,k)*Dxidx(k,j)
            enddo
            a(i,j) = s
          enddo
        enddo
        do i=1,3
          write(*,7004) i, a(i,1:3)
 7004     format('geom: i, Dxdxi(i,*)*Dxidx(*,j),j=1,2,3 = ', &
          i1,2x,3e12.5)
        enddo
        call pause
      endif
#endif
!
   end subroutine geom
