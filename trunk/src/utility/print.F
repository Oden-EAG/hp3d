c----------------------------------------------------------------------
c
c     routine name       - print_shapeE
c
c----------------------------------------------------------------------
c
c     latest revision    - Jan 11
c
c     purpose            - print element shape functions
c
c     arguments :
c     in:
c       Type      - element type
c       Norder    - element order of approximation
c       ShapE     - values of shape functions
c       CurlE     - the corresponding curl
c
c----------------------------------------------------------------------
c
      subroutine print_shapeE(Type,Norder,ShapE,CurlE)
c
      use element_data
      use parameters
#include "syscom.blk"
c
      character(len=4) :: Type
      dimension Norder(19),
     .          ShapE(3,MAXbrickE),CurlE(3,MAXbrickE)
c
      k=0
      do i=1,nedge(Type)+nface(Type)+1
        j=nvert(Type)+i
        call ndof_nod(type_nod(Type, j),Norder(i),
     .                ndofH,ndofE,ndofV,ndofQ)
        write(*,7001) j,type_nod(Type, j), Norder(i), ndofE
 7001   format('NODE = ',i2,' TYPE = ',a4,' NORD = ',i2, ' NDOF = ',i3)

        do j=1,ndofE
          k=k+1
          write(*,7002) k, ShapE(1:3,k), CurlE(1:3,k)
 7002     format('k= ',i3,' ShapE = ',3(e12.5,2x),
     .           '  CurlE = ',3(e12.5,2x))
        enddo
      enddo
c
c
      end subroutine print_shapeE
c----------------------------------------------------------------------
c
c     routine name       - print_shapeH
c
c----------------------------------------------------------------------
c
c     latest revision    - Jan 11
c
c     purpose            - print element shape functions
c
c     arguments :
c     in:
c       Type      - element type
c       Norder    - element order of approximation
c       ShapH     - values of shape functions
c       DshapH    - the corresponding curl
c
c----------------------------------------------------------------------
c
      subroutine print_shapeH(Type,Norder,ShapH,DshapH)
c
      use element_data
      use parameters
#include "syscom.blk"
c
      character(len=4) :: Type
      dimension Norder(19),
     .          ShapH(MAXbrickE),DshapH(3,MAXbrickE)
c
      k=0
      do i=1,nvert(Type)
        write(*,7001) i,type_nod(Type, i), 0, ndofH
        do j=1,ndofH
          k=k+1
          write(*,7002) k, ShapH(k), DshapH(1:3,k)
        enddo
      enddo

      do i=1,nedge(Type)+nface(Type)+1
        j=nvert(Type)+i
        call ndof_nod(type_nod(Type, j),Norder(i),
     .                ndofH,ndofE,ndofV,ndofQ)
        write(*,7001) j,type_nod(Type, j), Norder(i), ndofH
 7001   format('NODE = ',i2,' TYPE = ',a4,' NORD = ',i2, ' NDOF = ',i3)

        do j=1,ndofH
          k=k+1
          write(*,7002) k, shapH(k), dshapH(1:3,k)
 7002     format('k= ',i3,' ShapH = ',(e12.5,2x),
     .           '  DshapH = ',3(e12.5,2x))
        enddo
      enddo
c
c
      end subroutine print_shapeH
