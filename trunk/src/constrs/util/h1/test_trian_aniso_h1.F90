#if HP3D_DEBUG

!----------------------------------------------------------------------
!
!   routine name       - test_trian_aniso_h1
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - routine evaluates constraint coefficients for
!                        the triangular master element
!
!----------------------------------------------------------------------
!
   subroutine test_trian_aniso_h1(Iref)
!
!--------------------------------------------------------------------
!
      use parameters
      use constraints
!
!     NOTE : define implicit variables
      implicit none
!
      real(8) :: shapsma(MAXquadH,MAXquadH), &
                 shapbig(MAXquadH,MAXquadH), &
                 void(2,MAXquadH),xibig(2),xisma(2)
!
!  ...order and edge orientations
      integer :: norder(5),norient(4)
!
      real(8) :: val
      integer :: i,k,nrb,nrs,nrdof
!
      integer :: iprint
      iprint=0
!
 777  continue
!
      write(*,*) 'test_trian_aniso_h1: SET nord '
      read(*,*) nord
!
!  ...set uniform order MAXP:
      norder(1:5) = nord
      norder(5)   = nord*10 + nord
!
!  ...set orientations
      norient = 0
!
!*********************************************************************
!
!  ...input from small quad and get xibig with bilinear map
      write(*,*) 'test_trian_aniso_h1: SET xisma '
      read(*,*) xisma(1:2)
      call map_quad(Iref, xisma, xibig)
      write(*,*) 'xisma =', xisma(1:2)
      write(*,*) 'xibig =', xibig(1:2)
!
!  ...shape functions of big element
!      call shapeHt(xibig,norder,norient, &
!                   nrdof,shapbig(1:MAXtriaH,1),void)
      call shape2DH(TRIA,xibig,norder,norient,  &
                    nrdof,shapbig(:,1),void)
      write(*,*) 'shapbig = '
      do k=1,nrdof
        write(*,*) k,shapbig(k,1)
      enddo
!
!  ...find small shape functions at this point:
!      call shapeHq(xisma,norder,norient, &
!                   nrdof,shapsma(1:MAXquadrH,1),void)
      call shape2DH(QUAD,xisma,norder,norient,  &
                    nrdof,shapsma(:,1),void)
!
!*********************************************************************
!
!  ...verify if big shape functions are right combinations of small
!     ones:
!
!  ...loop through big bubble functions - central node:
      nrb = 3*nord
      do i=1,(nord-1)*(nord-2)/2
!
!  .....initiate value of the linear combination:
        val = 0.d0
!
!  .......small node 2 ( quad ):
          nrs = 4*nord
cc          do j2=1,nord-1
cc          do j1=1,nord-1
cc            jj = (j2-1)*(MAXP-1)+j1
cc            j = (j2-1)*(nord-1)+j1
cc            val = val + RRQH(Iref,1,i,jj)*shapsma(nrs+j,1)
cc          enddo
cc          enddo
          do j=1,(nord-1)**2
            val = val + get_rrqh(Iref, nord, 1, i, j)*shapsma(nrs+j,1)
          enddo
!
!  .......small node 3 ( medg ):
          nrs = 4 + (nord-1)*(mod(Iref+3,3)+1)
          do j=1,nord-1
            val = val + RRTH(1,i, Iref+3,j)*shapsma(nrs+j,1)
          enddo
!
!
         write(*,7002) i,shapbig(nrb+i,1),abs(val-shapbig(nrb+i,1))
 7002    format('test_trian_aniso_h1: i,shapebig,difference = ',i3,2e12.5)
!
!  ...end of loop through shape functions of big element
      enddo
!
      write(*,*) 'test_trian_aniso_h1: CONTINUE ?(1/0)'
      read(*,*) ians
      if (ians.eq.1) go to 777
!
      return
!
   end subroutine test_trian_aniso_h1

#endif
