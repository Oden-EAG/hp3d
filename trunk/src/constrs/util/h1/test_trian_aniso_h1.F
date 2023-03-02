#if DEBUG_MODE

c----------------------------------------------------------------------
c
c   routine name       - test_trian_aniso_h1
c
c----------------------------------------------------------------------
c
c   latest revision    - May 2020
c
c   purpose            - routine evaluates constraint coefficients for
c                        the triangular master element
c
c----------------------------------------------------------------------
c
      subroutine test_trian_aniso_h1(Iref)
c
c--------------------------------------------------------------------
c
      use parameters
      use constraints
c
c     NOTE : define implicit variables
      implicit none
c
      real(8) :: shapsma(MAXquadH,MAXquadH),
     .           shapbig(MAXquadH,MAXquadH),
     .           void(2,MAXquadH),xibig(2),xisma(2)
c
c  ...order and edge orientations
      integer :: norder(5),norient(4)
c
      real(8) :: val
      integer :: i,k,nrb,nrs,nrdof
c
      integer :: iprint = 0
c
 777  continue
c
      write(*,*) 'test_trian_aniso_h1: SET nord '
      read(*,*) nord
c
c  ...set uniform order MAXP:
      norder(1:5) = nord
      norder(5)   = nord*10 + nord
c
c  ...set orientations
      norient = 0
c
c*********************************************************************
c
c  ...input from small quad and get xibig with bilinear map
      write(*,*) 'test_trian_aniso_h1: SET xisma '
      read(*,*) xisma(1:2)
      call map_quad(Iref, xisma, xibig)
      write(*,*) 'xisma =', xisma(1:2)
      write(*,*) 'xibig =', xibig(1:2)
c
c  ...shape functions of big element
c      call shapeHt(xibig,norder,norient,
c     .             nrdof,shapbig(1:MAXtriaH,1),void)
      call shape2DH(TRIA,xibig,norder,norient,
     .              nrdof,shapbig(:,1),void)
      write(*,*) 'shapbig = '
      do k=1,nrdof
        write(*,*) k,shapbig(k,1)
      enddo
c
c  ...find small shape functions at this point:
c      call shapeHq(xisma,norder,norient,
c     .             nrdof,shapsma(1:MAXquadrH,1),void)
      call shape2DH(QUAD,xisma,norder,norient,
     .              nrdof,shapsma(:,1),void)
c
c*********************************************************************
c
c  ...verify if big shape functions are right combinations of small
c     ones:
c
c  ...loop through big bubble functions - central node:
      nrb = 3*nord
      do i=1,(nord-1)*(nord-2)/2
c
c  .....initiate value of the linear combination:
        val = 0.d0
c
c  .......small node 2 ( quad ):
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
c
c  .......small node 3 ( medg ):
          nrs = 4 + (nord-1)*(mod(Iref+3,3)+1)
          do j=1,nord-1
            val = val + RRTH(1,i, Iref+3,j)*shapsma(nrs+j,1)
          enddo
c
c
         write(*,7002) i,shapbig(nrb+i,1),abs(val-shapbig(nrb+i,1))
 7002    format('test_trian_aniso_h1: i,shapebig,difference = ',i3,2e12.5)
c
c  ...end of loop through shape functions of big element
      enddo
c
      write(*,*) 'test_trian_aniso_h1: CONTINUE ?(1/0)'
      read(*,*) ians
      if (ians.eq.1) go to 777
c
      return
c
      end subroutine test_trian_aniso_h1

#endif
