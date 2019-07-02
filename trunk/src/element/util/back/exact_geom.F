c----------------------------------------------------------------------
c
c   routine name       - exact_geom
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 09
c
c   purpose            - routine computes physical coordinates for
c                        a point in the master element using the
c                        exact geometry maps
c
c   arguments :
c     in:
c             Mdle     - an element middle node, identified with
c                        the element
c             Xi       - a point within the corresponding master
c                        element
c     out:
c             X        - coordinates of the corresponding physical
c                        point
c             Dxdxi    - derivatives wrt master element coordinates
c
c   required  routines -
c
c----------------------------------------------------------------------
c
      subroutine exact_geom(Mdle,Xi, X,Dxdxi)
c
      use GMP
      use data_structure3D
#include "syscom.blk"
c
      dimension Xi(3),X(3),Dxdxi(3,3)
      dimension eta(3),dxdeta(3,3)
c
      iprint=0
c
      nfath = NODES(Mdle)%father
      if (nfath.gt.0) then
        write(*,*) 'exact_geom: UNFINISHED '
        stop 1
      endif
      nel = -nfath
      eta(1:3) = Xi(1:3)
c
c  ...the corresponding GMP block
      nbl = ELEMS(nel)%GMPblock
      call decode(nbl, no,lab)
      select case(lab)
      case(1)
        call prism(no,eta, X,dxdeta)
      case(2)
        call hexa(no,eta, X,dxdeta)
      case(3)
        call tetra(no,eta, X,dxdeta)
      case(4)
        call pyram(no,eta, X,dxdeta)
      case default
        write(*,*) 'exact_geom: lab = ',lab
        stop 1
      end select
c
      Dxdxi(1:3,1:3) = dxdeta(1:3,1:3)
c
c
      end



