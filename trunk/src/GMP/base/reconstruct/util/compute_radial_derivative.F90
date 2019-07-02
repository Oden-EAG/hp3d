!-------------------------------------------------------------------------------
!
      subroutine compute_radial_derivative(Nt,Ie,S, Der,DDer,DDDer)
!      
!-------------------------------------------------------------------------------
      use GMP
!-------------------------------------------------------------------------------
      implicit none
!-------------------------------------------------------------------------------
!     DUMMY ARGUMENTS
      integer,intent(in)                :: Nt
      integer,intent(in)                :: Ie
      real*8,intent(in)               :: S
      real*8,dimension(3),intent(out) :: Der
      real*8,dimension(3),intent(out) :: DDer
      real*8,dimension(3),intent(out) :: DDDer
!-------------------------------------------------------------------------------
!     LOCAL VARIABLES
      real*8 :: c0,c1,c2,poly00,poly10,poly11,poly20,poly21,poly22,dpoly
      integer :: i,j,k,l,iprint
      integer, parameter :: deg = 7
!-------------------------------------------------------------------------------
!     EXTERNAL PROCEDURES
      integer, external :: bijec
!-------------------------------------------------------------------------------
!
      iprint=0
!
      if (iprint.eq.1) then
!        write(*,1000)Nt,Ie,S
 1000   format(' compute_radial_derivative: Nt = ',i4,'; Ie = ',i1,'; S = ',e12.5)
      endif
!
!  ...accumulate
      Der = 0.d0;  DDer = 0.d0;  DDDer = 0.d0
!  ...loop over control points      
      do j = 0,1
        do i = 0,deg-j
          k = deg-i-j
!  .......compute coefficients            
          call Bernstein_poly(i,deg-j,S, poly00,dpoly)
          c0 = deg*poly00*(-1.d0)**(j+1)
!          
          call Bernstein_poly(i-1,deg-j-1,S, poly10,dpoly)
          call Bernstein_poly(i  ,deg-j-1,S, poly11,dpoly)
          c1 = deg*(deg-j)*(poly10 - poly11)*(-1.d0)**(j+1)
!          
          call Bernstein_poly(i-2,deg-j-2,S, poly20,dpoly)
          call Bernstein_poly(i-1,deg-j-2,S, poly21,dpoly)
          call Bernstein_poly(i  ,deg-j-2,S, poly22,dpoly)
          c2 = deg*(deg-j)*(deg-j-1)*(poly20 - 2.d0*poly21 + poly22)*(-1.d0)**(j+1)
!          
          if (iprint.eq.1) then
!            write(*,1002)i,j,c0
 1002       format('   i = ',i1,'; j = ',i1,' --> c0 = ',e12.5)
!            write(*,1005)c1
 1005       format('                --> c1 = ',e12.5)
!            write(*,1006)c2
 1006       format('                --> c2 = ',e12.5)
          endif
!  .......account for edge 
          select case(ie)
          case(1)
            l = bijec(i,j)      
          case(2)
            l = bijec(k,i)      
          case(3)        
            l = bijec(j,k)      
          endselect
          if (iprint.eq.1) then
!            write(*,1003)l,l+2
 1003       format('   control point = Rdata(',i2,':',i2,')')           
          endif
          Der   = Der   + c0*TRIANGLES(Nt)%Rdata(l:l+2)
          DDer  = DDer  + c1*TRIANGLES(Nt)%Rdata(l:l+2)
          DDDer = DDDer + c2*TRIANGLES(Nt)%Rdata(l:l+2)
          if (iprint.eq.1)then
!            write(*,1009)Der
 1009       format('                --> Der   = ',3(e12.5,2x))      
!            write(*,1010)DDer
 1010       format('                --> DDer  = ',3(e12.5,2x))      
!            write(*,1011)DDDer
 1011       format('                --> DDDer = ',3(e12.5,2x))      
          endif
        enddo
!  ...end of loop over control points involved          
      enddo 
!
      if (iprint.eq.1) then
!        write(*,1001)Der
 1001   format('   Der   = ',3(e12.5,2x))      
!        write(*,1007)DDer
 1007   format('   DDer  = ',3(e12.5,2x))      
!        write(*,1008)DDDer
 1008   format('   DDDer = ',3(e12.5,2x))      
      endif
!
      end
