c
      subroutine plotgeomVTK
c
      use GMP
      use vtk_interf
c
#include "syscom.blk"
c
      allocatable :: x(:,:)
      dimension eta(3),dxdeta(3,3)
c
      dimension zdummy(1,1)
c
c-----------------------------------------------------------------------
c
      iprint = 1
c
      isol = 0; zdummy = 0.d0
c
      call init_vtk_geom
c
      call set_nrfig(isol,NRHEXAS+NRTETRA+NRPRISM+NRPYRAM)
c
c  ...set or read resolution and set number of subdivisions
      nres = -1
      do while(nres<0)
        write(*,*)'plotgeoVTK: ENTER NONNEGATIVE RESOLUTION:'
        read(*,*) nres
      enddo
      nsub = 2**nres
c
c  ...allocate array of points (physical coordinates)
      allocate(x(3,(nsub+1)**3), STAT=ic)
      if (ic/=0) then
        write(*,*)'plotgeoVTK: FAILED ALLOCATION [x]!'
        stop
      endif
c
      deta = 1.d0/dble(nsub)
c
      if (iprint>=1) write(*,100) NRHEXAS
 100  format('plotgeomVTK: NRHEXAS =',i6)
c
c  ...loop through hexahedra
      do no=1,NRHEXAS
c
c  .....generate array of points ordered lexicographically with respect
c       to reference coordinates
        l = 0
        do l3=0,nsub
        do l2=0,nsub
        do l1=0,nsub
          l = l+1
c
          eta(1) = l1*deta
          eta(2) = l2*deta
          eta(3) = l3*deta
c
c  .......evaluate physical coordinates
          call hexa(no,eta, x(1,l),dxdeta)
c
        enddo
        enddo
        enddo
c
c  .....pass data to VTK
        call add_hexa(isol,x,zdummy,nsub,no)
c
      enddo
c
      if (iprint>=1) write(*,200) NRTETRA
 200  format('plotgeomVTK: NRTETRA =',i6)
c
c  ...loop through tetrahedra
      do no=1,NRTETRA
c
c  .....generate array of points ordered lexicographically with respect
c       to reference coordinates
        l = 0
        do l3=0,nsub
        do l2=0,nsub-l3
        do l1=0,nsub-l3-l2
          l = l+1
c
          eta(1) = l1*deta
          eta(2) = l2*deta
          eta(3) = l3*deta
c
c  .......evaluate physical coordinates
          call tetra(no,eta, x(1,l),dxdeta)
c
        enddo
        enddo
        enddo
c
c  .....pass data to VTK
        call add_tetra(isol,x,zdummy,nsub,no)
c
      enddo
c
      if (iprint>=1) write(*,300) NRPRISM
 300  format('plotgeomVTK: NRPRISM =',i6)
c
c  ...loop through prisms
      do no=1,NRPRISM
c
c  .....generate array of points ordered lexicographically with respect
c       to reference coordinates
        l = 0
        do l3=0,nsub
        do l2=0,nsub
        do l1=0,nsub-l2
          l = l+1
c
          eta(1) = l1*deta
          eta(2) = l2*deta
          eta(3) = l3*deta
c
c  .......evaluate physical coordinates
          call prism(no,eta, x(1,l),dxdeta)
c
        enddo
        enddo
        enddo
c
c  .....pass data to VTK
        call add_prism(isol,x,zdummy,nsub,no)
c
      enddo
c
      if (iprint>=1) write(*,400) NRPYRAM
 400  format('plotgeomVTK: NRPYRAM =',i6)
c
c  ...loop through pyramids
      do no=1,NRPYRAM
c
c  .....generate array of points ordered lexicographically with respect
c       to reference coordinates
        l = 0
        do l3=0,nsub
        do l2=0,nsub-l3
        do l1=0,nsub-l3
          l = l+1
c
          eta(1) = l1*deta
          eta(2) = l2*deta
          eta(3) = l3*deta
c
c  .......evaluate physical coordinates
          call pyram(no,eta, x(1,l),dxdeta)
c
        enddo
        enddo
        enddo
c
c  .....pass data to VTK
        call add_pyramid(isol,x,zdummy,nsub,no)
c
      enddo
c
      deallocate(x, STAT=ic)
      if (ic/=0) then
        write(*,*)'plotgeomVTK: FAILED DEALLOCATION [x]!'
        stop
      endif
c
      call delete_vtk_geom
c
      end subroutine plotgeomVTK
