c----------------------------------------------------------------------
c
c   routine name       - shape3E_0
c
c----------------------------------------------------------------------
c
c   latest revision    - Nov 10
c
c   purpose            - routine evaluates H(curl) shape functions for
c                        3D elements of various types (version developed
c                        by K. Kim)
c
c   arguments :
c     in:
c         Type         - element type
c         Xi           - master element coordinates
c         Norder       - order of approximation for the nodes
c         Nedge_orient - edge orientations
c         Nface_orient - face orientations
c     out:
c         Nrdof        - number of the element shape functions
c         ShapE        - shape functions
c         CurlE        - curl of shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape3E_0(Type,Xi,Norder,Nedge_orient,Nface_orient,
     .                     Nrdof,ShapE,CurlE)
c
      use parameters
#include "syscom.blk"
c
      character (len=*),intent(in) :: Type
      dimension Nedge_orient(12),Nface_orient(6),Norder(19)
      dimension Xi(3),ShapE(3,MAXbrickE),CurlE(3,MAXbrickE)
c
c----------------------------------------------------------------------
c
      iprint=0
c
      if (iprint.eq.1) then
        write(*,7001)Type,Xi(1:3)
7001    format(' shape3E: Type,Xi = ',a5,2x,3(e12.5,2x))
      end if
c
      select case(Type)
      case('pris','mdlp')
        call shape3pE(Xi,Norder,Nedge_orient,Nface_orient,
     .                Nrdof,ShapE,CurlE)
      case('tetr','mdln')
        call shape3nE(Xi,Norder,Nedge_orient,Nface_orient,
     .                Nrdof,ShapE,CurlE)
      case('bric','mdlb')
        call shape3bE(Xi,Norder,Nedge_orient,Nface_orient,
     .                Nrdof,ShapE,CurlE)
      case('pyra','mdld')
        write(*,*) 'shape3E: Hcurl pyramid functions not implemented.'
        write(*,*) '         Switch to SHAPE=1.'
        stop
      case default
        write(*,*) 'shape3E_0: error, Type = ', Type
        stop
      endselect
c
c
      endsubroutine shape3E_0
