c----------------------------------------------------------------------
c
c   routine name       - hpvert
c
c----------------------------------------------------------------------
c
c   latest revision    - Dec 07
c
c   purpose            - routine determines geometry dof for a vertex
c                        node      
c
c   arguments :
c     in:
c            Iflag     - a flag specifying which of the objects the
c                        vertex node is generated on
c                        = 5  for a prism
c                        = 6  for a hexahedron
c                        = 7  for a tetrahedron
c                        = 8  for a pyramid
c            No        - number of a specific object (point,curve,triangle,
c                        rectangle...)
c            Xi        - reference coordinates of the point
c     out:
c            Xnod      - geometry dof      
c
c----------------------------------------------------------------------
c
      subroutine hpvert(Iflag,No,Xi, Xnod)
c     
      use parameters , only : NDIMEN
#include "syscom.blk"
c
      dimension Xi(3),Xnod(3)
      dimension void(3,3)
c
c----------------------------------------------------------------------
c
      iprint=0
c
c  ...determine coordinates of the point
      select case(Iflag)
      case(5);            call prism(No,Xi, Xnod,void)
      case(6);            call hexa( No,Xi, Xnod,void)
      case(7);            call tetra(No,Xi, Xnod,void)
      case(8);            call pyram(No,Xi, Xnod,void)
      case default
         print *, "hpvert: unknown node type =", iflag
         stop 1
      end select
c     
c  ...printing
      if (iprint.eq.1) then
        write(*,7005) Iflag,No,Xi,Xnod
 7005   format('hpvert: Iflag,No = ',2i8,' Xi = ',3f8.3,
     .          '  Xnod = ',3f8.3)
        call pause
      endif
c     
c
      end
