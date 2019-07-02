      subroutine check_normal(Iwork,iprint, error)
      use GMP
      implicit none
      integer, dimension(2,2),intent(in) :: Iwork
      real*8,dimension(2) :: eta
      real*8,dimension(3) :: x,vec
      real*8,dimension(3,2) :: dx,vec_aux
      real*8 :: nsub = 12.d0
      integer:: isub = 12
      real*8 :: error,s
      integer :: nt,i,it,iprint


      if (iprint.eq.1) then  
      write(*,*)'check_normal...'
      endif
      error = 0.d0
!  ...loop over subdivisions      
      do i=0, isub
        do it = 1,2
          nt = Iwork(it,1)
          select case(Iwork(it,2))
          case(1)
            eta(2) = 0.d0
            eta(1) = i/nsub
          case(-1)        
            eta(2) = 0.d0
            eta(1) = 1.d0 - i/nsub     
          case(2)
            eta(1) = 1.d0 - i/nsub
            eta(2) = i/nsub
          case(-2)
            eta(1) = i/nsub
            eta(2) = 1.d0 - i/nsub
          case(3)
            eta(1) = 0.d0
            eta(2) = i/nsub
          case(-3)        
            eta(1) = 0.d0
            eta(2) = 1.d0 - i/nsub      
          endselect
          call trian(nt,eta,x, dx)
!  .......compute normal          
          call cross_product(dx(:,1),dx(:,2), vec)
          call normalize(vec)
          vec_aux(1:3,it) = vec
          if (iprint.eq.1) then
            write(*,8000)i,nt,Iwork(it,2),eta(1:2),vec(1:3)
 8000       format(' isub = ',i2,',   nt = ',i3,', ie = ',i2,' eta = (',2(e12.5,2x),'), &
                     norm = (',3(e12.5,2x),')')
          endif
!  .....end of loop over attached triangles
        enddo
!  .....compute diffence b/w normals        
        vec_aux(1:3,1) = vec_aux(1:3,1) - vec_aux(1:3,2) 
        call norm(vec_aux(1:3,1), s)
!  .....update error        
        error = error + s
!  ...end of loop over subdivisions        
      enddo
      if (iprint.eq.1) then
        write(*,8001)error
 8001   format(' --> error = ',e12.5)
        write(*,*)'-----------------------------------------------------'
      endif

      end
