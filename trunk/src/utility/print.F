c----------------------------------------------------------------------
c
c     routine name       - print_shapeE
c
c----------------------------------------------------------------------
c
c     latest revision    - Feb 2023
c
c     purpose            - print element shape functions
c
c     arguments :
c     in:
c       Ntype     - element type
c       Norder    - element order of approximation
c       ShapE     - values of shape functions
c       CurlE     - the corresponding curl
c
c----------------------------------------------------------------------
c
      subroutine print_shapeE(Ntype,Norder,ShapE,CurlE)
c
      use element_data
      use parameters
c
      implicit none
c
      integer, intent(in) :: Ntype,Norder(19)
      real(8), intent(in) :: ShapE(3,MAXbrickE),CurlE(3,MAXbrickE)
c
      integer :: i,j,k,ndofH,ndofE,ndofV,ndofQ
c
      k=0
      do i=1,nedge(Ntype)+nface(Ntype)+1
        j=nvert(Ntype)+i
        call ndof_nod(type_nod(Ntype,j),Norder(i),
     .                ndofH,ndofE,ndofV,ndofQ)
        write(*,7001) j,S_Type(type_nod(Ntype,j)), Norder(i), ndofE
 7001   format('NODE = ',i2,' TYPE = ',a4,' NORD = ',i2, ' NDOF = ',i3)

        do j=1,ndofE
          k=k+1
          write(*,7002) k, ShapE(1:3,k), CurlE(1:3,k)
 7002     format('k= ',i3,' ShapE = ',3(e12.5,2x),
     .           '  CurlE = ',3(e12.5,2x))
        enddo
      enddo
c
      end subroutine print_shapeE
c
c----------------------------------------------------------------------
c
c     routine name       - print_shapeH
c
c----------------------------------------------------------------------
c
c     latest revision    - Feb 2023
c
c     purpose            - print element shape functions
c
c     arguments :
c     in:
c       Ntype     - element type
c       Norder    - element order of approximation
c       ShapH     - values of shape functions
c       DshapH    - the corresponding curl
c
c----------------------------------------------------------------------
c
      subroutine print_shapeH(Ntype,Norder,ShapH,DshapH)
c
      use element_data
      use parameters
c
      implicit none
c
      integer, intent(in) :: Ntype,Norder(19)
      real(8), intent(in) :: ShapH(MAXbrickE),DshapH(3,MAXbrickE)
c
      integer :: i,j,k,ndofH,ndofE,ndofV,ndofQ
c
      k=0
      do i=1,nvert(Ntype)
        write(*,7001) i, S_Type(type_nod(Ntype,i)), 0, ndofH
        do j=1,ndofH
          k=k+1
          write(*,7002) k, ShapH(k), DshapH(1:3,k)
        enddo
      enddo

      do i=1,nedge(Ntype)+nface(Ntype)+1
        j=nvert(Ntype)+i
        call ndof_nod(type_nod(Ntype,j),Norder(i),
     .                ndofH,ndofE,ndofV,ndofQ)
        write(*,7001) j, S_Type(type_nod(Ntype,j)), Norder(i), ndofH
 7001   format('NODE = ',i2,' TYPE = ',a4,' NORD = ',i2, ' NDOF = ',i3)

        do j=1,ndofH
          k=k+1
          write(*,7002) k, shapH(k), dshapH(1:3,k)
 7002     format('k= ',i3,' ShapH = ',(e12.5,2x),
     .           '  DshapH = ',3(e12.5,2x))
        enddo
      enddo
c
      end subroutine print_shapeH
