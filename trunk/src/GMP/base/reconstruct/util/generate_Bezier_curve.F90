!--------------------------------------------------------------------------------
      subroutine generate_Bezier_curve(Nc,X,V,A)
!--------------------------------------------------------------------------------
      use GMP
!--------------------------------------------------------------------------------     
      implicit none
!--------------------------------------------------------------------------------      
!     DUMMY ARGUMENTS
!  ...curve number
      integer,intent(in) :: Nc
!  ...positions, velocities, accelerations      
      real(8),dimension(3,2),intent(in) :: X,V,A
!--------------------------------------------------------------------------------
      integer :: is,iprint,i
!
      iprint=0
      if (iprint.eq.1) then
        write(*,*)'generate_Bezier_curve: Nc = ',Nc
      endif
!
!  ...rename curve
      CURVES(Nc)%Type = '5Bezier'
!  ...allocate control points
      deallocate (CURVES(Nc)%Rdata, STAT=is)
      if (is.ne.0) then
        write(*,*)'generate_Bezier_curve: Rdata not deallocted for Nc = ',Nc
        stop
      endif
      allocate (CURVES(Nc)%Rdata(-3:20), STAT=is)
      if (is.ne.0) then
        write(*,*)'generate_Bezier_curve: Rdata not allocted for Nc = ',Nc
        stop
      endif
      CURVES(Nc)%Rdata = 0.d0
!
!  ...endpoints:                                        b_0, b_5
      CURVES(Nc)%Rdata(0:2)   = X(1:3,1)
      CURVES(Nc)%Rdata(15:17) = X(1:3,2)
!  ...velocity control points:                          b_1, b_4
      CURVES(Nc)%Rdata(3:5)   = X(1:3,1) + V(1:3,1)/5.d0
      CURVES(Nc)%Rdata(12:14) = X(1:3,2) - V(1:3,2)/5.d0
!  ...acceleration control points:                      b_2, b_3      
      CURVES(Nc)%Rdata(6:8)   = X(1:3,1) + 2.d0*V(1:3,1)/5.d0 + A(1:3,1)/20.d0
      CURVES(Nc)%Rdata(9:11) = X(1:3,2) - 2.d0*V(1:3,2)/5.d0 + A(1:3,2)/20.d0

      if (iprint.eq.1) then
        write(*,*)'generate_Bezier_curve: control points'
        do i = 0,5
          write(*,7000)i,CURVES(Nc)%Rdata(3*i:3*i+2)
 7000     format(' ',i1,':   ',3(e12.5,2x))
        enddo
        write(*,*)'----------------------------------------------------'
      endif
!            
!
      end subroutine
