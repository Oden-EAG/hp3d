      subroutine test_normal

      use GMP
      implicit none
      real*8 :: x,y,rz,error,error_save,x_save,rz_save,x_in,y_in,rz_in
      integer, dimension(2,2) :: iwork
      integer :: idec,i,j,k,iidec,iflag,irow,bijec,iprint,it,ie,nt,i1,j1,i2,j2,i3,j3,i4,j4
      integer :: nsub
      real*8 :: x1,x2,rz1,rz2,temp,r_c,y1,y2,y_save
      real*8,dimension(3) :: r_der,b_lo_save,b_up_save
      real*8,external :: Bern_poly
      integer, external :: fact,imod
      real*8, dimension(3) :: XX
      real*8, dimension(3,2) :: Dxdeta

      iprint=0

      nsub = 100
!  ...triangles data
      iwork(1,1) =  3; iwork(1,2) =  1
      iwork(2,1) = 14; iwork(2,2) = -1
!

 10   continue
!
      write(*,*)'====================================================='
      write(*,*)'QUIT ...............................................0'
      write(*,*)''
      write(*,*)'move b_12 ..........................................1'
      write(*,*)'move b_22 ..........................................2'
      write(*,*)'move b_21 ..........................................3'
      write(*,*)'optimal location of b21 ............................4'
      write(*,*)'move b_11 and b_31 .................................5'
      write(*,*)'test Bernstein polynomial...........................6'
      write(*,*)'test factorial .          ..........................7'
      write(*,*)'test normal derivative .  ..........................8'
      write(*,*)'test imod .               ..........................9'
      write(*,*)'test trian G1               .......................10'
      write(*,*)'compute radial derivative   .......................11'
      write(*,*)'test triangle 3   .................................12'
      write(*,*)'====================================================='
      read(*,*)iidec
!
      select case(iidec)
      case(12)
 76     continue
        write(*,*)'x = '
        read(*,*)x
        write(*,*)'y = '
        read(*,*)y
        call trian(3,(/x,y/), XX,Dxdeta)
!        write(*,*)'Dxdeta_1 = ',Dxdeta(1:3,1)
!        write(*,*)'Dxdeta_2 = ',Dxdeta(1:3,2)
        goto 76
      case(11)
 15     continue
        write(*,*)'ie = '
        read(*,*)i
        write(*,*)'s = '
        read(*,*)x
        call compute_radial_derivative(3,i,x, XX,r_der,b_lo_save)
        goto 15

      case(10)
        call generate_Bezier_triangle(3)

      case(9)
 14     continue
        write(*,*)'i = '
        read(*,*)i
        write(*,*)'mod = '
        read(*,*)k
        write(*,*)'imod = ',imod(i,k)
        goto 14

      case(8)
 13     continue
        write(*,*)'s = '
        read(*,*)r_c
        call normal_derivative(r_c)
        goto 13

      case(7)
 12     continue
        write(*,*)'n = '
        read(*,*)i
        write(*,*)'fact = ',fact(i)
        goto 12

      case(6)
 11     continue
        write(*,*)'i = '
        read(*,*)i
        write(*,*)'n = '
        read(*,*)k
        write(*,*)'t = '
        read(*,*)r_c
        write(*,*)' Bern_poly = ',Bern_poly(i,k,r_c)
        goto 11
!
      case(0)
        return
!
      case(1)
!  .....save initial values
        x_in  = TRIANGLES(3)%Rdata(36)
        y_in  = TRIANGLES(3)%Rdata(37)
        rz_in = TRIANGLES(3)%Rdata(38)
!
 20     continue
        write(*,1000)TRIANGLES(3)%Rdata(36:38)
 1000   format(' b_12 = ',3(e12.5,2x))
        call check_normal(iwork,1, error)
        write(*,*)'move control point? 1-Yes'
        read(*,*)idec
        if (idec.eq.1) then
          write(*,*)'set x:'
          read(*,*)x
          write(*,*)'set y:'
          read(*,*)y
          TRIANGLES(3)%Rdata(36) = x
          TRIANGLES(3)%Rdata(37) = y
          TRIANGLES(3)%Rdata(38) = x
          goto 20
        endif
!
!  .....reset initial values
        TRIANGLES(3)%Rdata(36) = x_in
        TRIANGLES(3)%Rdata(37) = y_in
        TRIANGLES(3)%Rdata(38) = rz_in
!
!
      case(2)
!  .....save initial values
        x_in  = TRIANGLES(3)%Rdata(39)
        y_in  = TRIANGLES(3)%Rdata(40)
        rz_in = TRIANGLES(3)%Rdata(41)
!
 30     continue
        write(*,1001)TRIANGLES(3)%Rdata(39:41)
 1001   format(' b_22 = ',3(e12.5,2x))
        call check_normal(iwork,1, error)
        write(*,*)'move control point? 1-Yes'
        read(*,*)idec
        if (idec.eq.1) then
          write(*,*)'set x:'
          read(*,*)x
          write(*,*)'set y:'
          read(*,*)y
          TRIANGLES(3)%Rdata(39) = x
          TRIANGLES(3)%Rdata(40) = y
          TRIANGLES(3)%Rdata(41) = y
          goto 30
        endif
!
!  .....reset initial values
        TRIANGLES(3)%Rdata(39) = x_in
        TRIANGLES(3)%Rdata(40) = y_in
        TRIANGLES(3)%Rdata(41) = rz_in
!
!
      case(3)
!  .....save initial values
        x_in  = TRIANGLES(3)%Rdata(24)
        y_in  = TRIANGLES(3)%Rdata(25)
        rz_in = TRIANGLES(3)%Rdata(26)
!
 40     continue
        write(*,1002)TRIANGLES(3)%Rdata(24:26)
 1002   format(' b_21 = ',3(e12.5,2x))
        call check_normal(iwork,1, error)
        iflag = 0
        write(*,*)'reset x? 1-Yes'
        read(*,*)idec
        if (idec.eq.1) then
          write(*,*)'set x:'
          read(*,*)x
          TRIANGLES(3)%Rdata(24) = x
          TRIANGLES(3)%Rdata(25) = x
          iflag = iflag + 1
        endif
        write(*,*)'reset z? 1-Yes'
        read(*,*)idec
        if (idec.eq.1) then
          write(*,*)'set z:'
          read(*,*)rz
          TRIANGLES(3)%Rdata(26) = rz
          iflag = iflag + 1
        endif
        if (iflag.gt.0) goto 40
!
!  .....reset initial values
        TRIANGLES(3)%Rdata(24) = x_in
        TRIANGLES(3)%Rdata(25) = y_in
        TRIANGLES(3)%Rdata(26) = rz_in
!
!
      case(4)
!  .....save initial values
        x_in  = TRIANGLES(3)%Rdata(24)
        y_in  = TRIANGLES(3)%Rdata(25)
        rz_in = TRIANGLES(3)%Rdata(26)
!
 50     continue
        write(*,1001)TRIANGLES(3)%Rdata(24:26)
        call check_normal(iwork,1, error)
        x1 = TRIANGLES(3)%Rdata(24); x2 = x1
        rz1 = TRIANGLES(3)%Rdata(26); rz2 = rz1
        write(*,*)'reset range for x? 1-Yes'
        read(*,*)idec
        if (idec.eq.1) then
          write(*,*)'set x1,x2:'
          read(*,*)x1,x2
          if (x2.lt.x1) then
            temp = x1
            x1 = x2
            x2 = temp
          endif
        endif
        write(*,*)'reset range for z? 1-Yes'
        read(*,*)idec
        if (idec.eq.1) then
          write(*,*)'set z1,z2:'
          read(*,*)rz1,rz2
          if (rz2.lt.rz1) then
            temp = rz1
            rz1 = rz2
            rz2 = temp
          endif
        endif
!
 60     continue
        error_save = 10000.d0
!
!  .....loop over subdivisions
        do i = 0,nsub
          do j = 0,nsub
            x  = x1 + i*(x2 - x1)/nsub
            rz = rz1 + j*(rz2 - rz1)/nsub
!           write(*,*)'x = ',x
!           write(*,*)'z = ',rz
!  .........update control point for triangles
            TRIANGLES(3)%Rdata(24) = x
            TRIANGLES(3)%Rdata(25) = x
            TRIANGLES(14)%Rdata(24) = x
            TRIANGLES(14)%Rdata(25) = x
            TRIANGLES(3)%Rdata(26) = rz
            TRIANGLES(14)%Rdata(26) = -rz
!  .........compute error
            call check_normal(iwork,0, error)
!!!         write(*,*)'error = ',error
!  .........save error and control point if needed
            if (error.lt.error_save) then
              error_save = error
              x_save = x
              rz_save = rz
!!           write(*,*)'error_save = ',error_save
!!           write(*,*)'z_save = ',rz_save
            endif
          enddo
!  .....end of loop over subdivisions
        enddo
!
!  .....set control point to best one
        TRIANGLES(3)%Rdata(24) = x_save
        TRIANGLES(3)%Rdata(25) = x_save
        TRIANGLES(14)%Rdata(24) = x_save
        TRIANGLES(14)%Rdata(25) = x_save
        TRIANGLES(3)%Rdata(26) = rz_save
        TRIANGLES(14)%Rdata(26) = -rz_save
!
!  .....check normal for new location of control point
        call check_normal(iwork,1, error)
!
        write(*,*)'---------------------------------------------'
        write(*,*)'error_save = ',error_save
        write(*,*)'x_save     = ',x_save
        write(*,*)'z_save     = ',rz_save
        write(*,*)'---------------------------------------------'
!
        write(*,*)'display graphics? 1-Yes'
        read(*,*)idec
        if (idec.eq.1) call graphg
!
        write(*,*)'iterate? 1-Yes'
        read(*,*)idec
        if (idec.eq.1) then
!  .......reset range for x and z around control point
          x1  = x_save  - 2.d0*(x2  - x1) /nsub
          x2  = x_save  + 2.d0*(x2  - x1) /nsub
          rz1 = rz_save - 2.d0*(rz2 - rz1)/nsub
          rz2 = rz_save + 2.d0*(rz2 - rz1)/nsub
          write(*,2000) x1,x2
 2000     format(' new range: ',e12.5,' <= x <= ',e12.5)
          write(*,2001) rz1,rz2
 2001     format(' new range: ',e12.5,' <= z <= ',e12.5)
          goto 60
        endif
!
!  .....reset initial values
        TRIANGLES(3)%Rdata(24) = x_in
        TRIANGLES(3)%Rdata(25) = y_in
        TRIANGLES(3)%Rdata(26) = rz_in
        TRIANGLES(14)%Rdata(24) = x_in
        TRIANGLES(14)%Rdata(25) = y_in
        TRIANGLES(14)%Rdata(26) = -rz_in
        goto 50
!
!
!
!
      case (5)
!  .....save initial values
!!        x_in  = TRIANGLES(3)%Rdata(21)
!!        y_in  = TRIANGLES(3)%Rdata(22)
!!        rz_in = TRIANGLES(3)%Rdata(23)
!---------------------------------------------------------------------------
!
!  STEP 1: compute normal derivative
!
        r_der = 0.d0
!  .....loop over attached triangles
        do i = 1,2
          select case(iwork(i,2))
          case(1,-1)
            i1 = 0; j1 = 0
            i2 = 0; j2 = 1
            i3 = 5; j3 = 0
            i4 = 4; j4 = 1
          case(2,-2)
            i1 = 5; j1 = 0
            i2 = 4; j2 = 0
            i3 = 0; j3 = 5
            i4 = 0; j4 = 4
          case(3,-3)
            i1 = 0; j1 = 5
            i2 = 1; j2 = 4
            i3 = 0; j3 = 0
            i4 = 1; j4 = 0
          endselect
          r_der = r_der + 1.25d0*(-1.d0)**(i+1)* &
           (TRIANGLES(iwork(i,1))%Rdata(bijec(i1,j1):bijec(i1,j1)+2) - &
            TRIANGLES(iwork(i,1))%Rdata(bijec(i2,j2):bijec(i2,j2)+2) +  &
           TRIANGLES(iwork(i,1))%Rdata(bijec(i3,j3):bijec(i3,j3)+2) - &
           TRIANGLES(iwork(i,1))%Rdata(bijec(i4,j4):bijec(i4,j4)+2))
!  .....end of loop over attached triangles
        enddo
!!        write(*,*)'desired normal derivative at midpoint:'
!!        write(*,4002)r_der
!! 4002   format(' r_der = ',3(e12.5,2x))
!
 70     continue
        write(*,*)'-------------------------------------------------'
        write(*,*)'control points:'
        write(*,3000)TRIANGLES(3)%Rdata(21:23)
 3000   format(' b_11 = ',3(e12.5,2x))
        write(*,3001)TRIANGLES(3)%Rdata(24:26)
 3001   format(' b_21 = ',3(e12.5,2x))
        write(*,3002)TRIANGLES(3)%Rdata(27:29)
 3002   format(' b_31 = ',3(e12.5,2x))
        call check_normal(iwork,1, error)
        write(*,*)'dislay graphics? 1-Yes'
        read(*,*)idec
        if (idec.eq.1) call graphg
!!        write(*,*)'checking derivative of upper triangle:'
!!        call check_derivative(3,1)
!!        write(*,*)'checking derivative of lower triangle:'
!!        call check_derivative(14,1)

!
!--------------------------------------------------------------------------------
!
!  STEP 2: set range for y and loop over subdivisions
!
        write(*,*)'set range for y'
        read(*,*)y1,y2

        nsub = 100000
        error_save = 100000.d0
!  .....loop over subdivisions
        do k = 0,nsub
          y  = y1 + k*(y2 - y1)/nsub
!  .......reset control points of 1st triangle
          TRIANGLES(3)%Rdata(22) = y
          TRIANGLES(3)%Rdata(23) = y
          TRIANGLES(3)%Rdata(27) = y
          TRIANGLES(3)%Rdata(29) = y
!  .......reset control points of 2nd triangle
          TRIANGLES(14)%Rdata(21) = y
          TRIANGLES(14)%Rdata(23) = -y
          TRIANGLES(14)%Rdata(28) = y
          TRIANGLES(14)%Rdata(29) = -y
!
!!          write(*,*)'-----------------------------------------------'
!!          write(*,*)'control points of upper triangle'
!!          write(*,3000)TRIANGLES(3)%Rdata(21:23)
!!          write(*,3002)TRIANGLES(3)%Rdata(27:29)
!!          write(*,*)'-----------------------------------------------'
!!          write(*,*)'control points of lower triangle'
!!          write(*,3000)TRIANGLES(14)%Rdata(21:23)
!!          write(*,3002)TRIANGLES(14)%Rdata(27:29)
!
!-------------------------------------------------------------------------------
!
!  STEP 2a: recompute interior control points
!
!  .......reinitialize to 0
          TRIANGLES(3)%Rdata(24:26) = 0.d0
          TRIANGLES(14)%Rdata(24:26) = 0.d0
!  .......loop over attached triangles
          do it = 1,2
            nt = iwork(it,1);  ie = iwork(it,2)
            r_der = (-1.d0)**(it+1)*r_der
            if (iprint.eq.1) then
              write(*,*)'------------------------------------------------'
              write(*,1013)it,r_der
 1013         format(' it = ',i1,'; derivative = ',3(e12.5,2x))
            endif
            do irow = 0,1
              do i = 0,(5-irow)
                j = 5 - i - irow
!  .............cycle if on control point
                if (i.eq.j) cycle
!!                write(*,*)'i,j = ',i,j
                if (irow.eq.0) then
                  r_c = -fact(5)/(fact(i)*fact(j))*5.d0*0.5d0**5
                else
                  r_c = fact(5)/(fact(i)*fact(j))*0.5d0**4
                endif
!  .............account for edge number
                select case(ie)
                case(1,-1)
                  TRIANGLES(nt)%Rdata(bijec(2,1):bijec(2,1)+2) =  &
                   TRIANGLES(nt)%Rdata(bijec(2,1):bijec(2,1)+2) +  &
                   r_c*TRIANGLES(nt)%Rdata(bijec(i,irow):bijec(i,irow)+2)
                  if (iprint.eq.1) then
                    write(*,1010)i,irow,  &
                   TRIANGLES(nt)%Rdata(bijec(2,1):bijec(2,1)+2)
 1010               format(' + c*b_',i1,i1,' = ',3(e12.5,2x))
                  endif
                case(2,-2)
                  TRIANGLES(nt)%Rdata(bijec(2,2):bijec(2,2)+2) =  &
                   TRIANGLES(nt)%Rdata(bijec(2,2):bijec(2,2)+2) +  &
                   r_c*TRIANGLES(nt)%Rdata(bijec(i,j):bijec(i,j)+2)
                  if (iprint.eq.1) then
                    write(*,1010)i,j,   &
                    TRIANGLES(nt)%Rdata(bijec(2,2):bijec(2,2)+2)
                  endif
                case(3,-3)
                  TRIANGLES(nt)%Rdata(bijec(1,2):bijec(1,2)+2) =  &
                   TRIANGLES(nt)%Rdata(bijec(1,2):bijec(1,2)+2) +   &
                   r_c*TRIANGLES(nt)%Rdata(bijec(irow,j):bijec(irow,j)+2)
                  if (iprint.eq.1) then
                    write(*,1010)irow,j,  &
                   TRIANGLES(nt)%Rdata(bijec(1,2):bijec(1,2)+2)
                  endif
                endselect
              enddo
            enddo
!  .........divide control point by coefficient
            select case(ie)
            case(1,-1)
              TRIANGLES(nt)%Rdata(bijec(2,1):bijec(2,1)+2) = -8.d0*  &
               (TRIANGLES(nt)%Rdata(bijec(2,1):bijec(2,1)+2) + r_der)/15.d0
              if (iprint.eq.1) then
                write(*,1011)TRIANGLES(nt)%Rdata(bijec(2,1):bijec(2,1)+2)
 1011           format(' --> b_21 = ',3(e12.5,2x))
              endif
            case(2,-2)
              TRIANGLES(nt)%Rdata(bijec(2,2):bijec(2,2)+2) = -8.d0*  &
             (TRIANGLES(nt)%Rdata(bijec(2,2):bijec(2,2)+2) + r_der)/15.d0
              if (iprint.eq.1) then
                write(*,1012)TRIANGLES(nt)%Rdata(bijec(2,2):bijec(2,2)+2)
 1012           format(' --> b_22 = ',3(e12.5,2x))
              endif
            case(3,-3)
              TRIANGLES(nt)%Rdata(bijec(1,2):bijec(1,2)+2) = -8.d0*  &
             (TRIANGLES(nt)%Rdata(bijec(1,2):bijec(1,2)+2) + r_der)/15.d0
              if (iprint.eq.1) then
                write(*,1014)TRIANGLES(nt)%Rdata(bijec(1,2):bijec(1,2)+2)
 1014           format(' --> b_12 = ',3(e12.5,2x))
              endif
            endselect
!  .......end of loop over attached triangles
          enddo
!
!  .......compute and update error if necessary
          call check_normal(iwork,0, error)
!  .......save error and control point if needed
          if (error.lt.error_save) then
            error_save = error
            y_save = y
            b_up_save = TRIANGLES(3)%Rdata(24:26)
            b_lo_save = TRIANGLES(14)%Rdata(24:26)
          endif
!
!  .......end of loop over subdivisions
          enddo
!
!-------------------------------------------------------------------------
!
!  STEP 2c: update values
!
!  .......midedge control points
          TRIANGLES(3)%Rdata(24:26) = b_up_save
          TRIANGLES(14)%Rdata(24:26) = b_lo_save
!  .......1st triangle corner points
          TRIANGLES(3)%Rdata(22) = y_save
          TRIANGLES(3)%Rdata(23) = y_save
          TRIANGLES(3)%Rdata(27) = y_save
          TRIANGLES(3)%Rdata(29) = y_save
!  .......2nd triangle corner points
          TRIANGLES(14)%Rdata(21) = y_save
          TRIANGLES(14)%Rdata(23) = -y_save
          TRIANGLES(14)%Rdata(28) = y_save
          TRIANGLES(14)%Rdata(29) = -y_save
!
!         write(*,*)'final normal derivative at midpoint:'
!         call check_derivative(3,1)
!         call check_derivative(14,1)
!

!         write(*,*)'---------------------------------------------------'
!         write(*,*)'new control points:'
!         write(*,4000)TRIANGLES(3)%Rdata(24:26)
! 4000    format(' upper b_21 = ',3(e12.5,2x))
!         write(*,4001)TRIANGLES(14)%Rdata(24:26)
! 4001    format(' lower b_21 = ',3(e12.5,2x))



          goto 70
!!        endif
!
!  .....reset initial values
!        TRIANGLES(3)%Rdata(36) = x_in
!        TRIANGLES(3)%Rdata(37) = y_in
!        TRIANGLES(3)%Rdata(38) = rz_in
!

!

!
      case default
        goto 10
!
      end select
!
      goto 10
!
!
      end


