!------------------------------------------------------------------------------------      
!
      subroutine check_derivative(Nt,Ie)
!      
!------------------------------------------------------------------------------------      
      use GMP
!------------------------------------------------------------------------------------      
      implicit none
!------------------------------------------------------------------------------------      
!     DUMMY ARGUMENTS
      integer, intent(in)  :: Nt,Ie
      real(8) :: Der(3)
!------------------------------------------------------------------------------------      
!     LOCAL VARIABLES      
      real(8) :: c
      integer :: bijec,fact
      integer :: i,j,irow
!------------------------------------------------------------------------------------      
!
      der = 0.d0
      do irow = 0,1
        do i = 0,(5-irow)
          j = 5 - i - irow
!  .......compute coefficient          
          if (irow.eq.0) then
            c = fact(5)/(fact(i)*fact(j))*5.d0*0.5d0**5
          else
            c = -fact(5)/(fact(i)*fact(j))*0.5d0**4              
          endif
!  .......account for edge number          
          select case(Ie)
          case(1)
            Der = Der + c*TRIANGLES(Nt)%Rdata(bijec(i,irow):bijec(i,irow)+2)
          case(2)
            Der = Der + c*TRIANGLES(Nt)%Rdata(bijec(i,j):bijec(i,j)+2)
          case(3)  
            Der = Der + c*TRIANGLES(Nt)%Rdata(bijec(irow,j):bijec(irow,j)+2)
          endselect
        enddo
      enddo
      write(*,1000)Nt,Ie,Der(1:3)
 1000 format(' nt = ',i3,', ie = ',i1,'; Der = ',3(e12.5,2x))

!
      end
!
!
!
!
!---------------------------------------------------------------------------
      integer function bijec(i,j)
!---------------------------------------------------------------------------
      implicit none
!---------------------------------------------------------------------------
      integer, intent(in) :: i,j
      integer :: deg,n,m
!---------------------------------------------------------------------------
!
      deg = 7
      bijec = 0
!  ...accumulate      
      do n = 0,(j - 1)
        do m = 0,(deg - n)
          bijec = bijec + 3
        enddo
      enddo
      do n = 0,(i - 1)
        bijec = bijec + 3
      enddo
!
!  ...debugging
!      write(*,1000)i,j,bijec
 1000 format(' bijec: i = ',i1,', j = ',i1,' --> k = ',i2)
!
      end
!
!
!
!!!!---------------------------------------------------------------------------
!!!      integer function bijec7(i,j)
!!!!---------------------------------------------------------------------------
!!!      implicit none
!!!!---------------------------------------------------------------------------
!!!      integer, intent(in) :: i,j
!!!      integer :: deg,n,m
!!!!---------------------------------------------------------------------------
!!!!
!!!      deg = 7
!!!      bijec = 0
!!!!  ...accumulate      
!!!      do n = 0,(j - 1)
!!!        do m = 0,(deg - n)
!!!          bijec = bijec + 3
!!!        enddo
!!!      enddo
!!!      do n = 0,(i - 1)
!!!        bijec = bijec + 3
!!!      enddo
!!!!
!!!!  ...debugging
!!!!      write(*,1000)i,j,bijec
!!! 1000 format(' bijec: i = ',i1,', j = ',i1,' --> k = ',i2)
!!!!
!!!      end
