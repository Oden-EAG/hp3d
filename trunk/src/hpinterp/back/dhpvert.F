c----------------------------------------------------------------------
c
c   routine name       - dhpvert
c
c----------------------------------------------------------------------
c
c   latest revision    - Dec 10
c
c   purpose            - routine determines Dirichlet dof for a vertex
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
c            Icase     - node case
c     out:
c            Xnod      - geometry dof
c
c----------------------------------------------------------------------
c
      subroutine dhpvert(Iflag,No,Xi,Icase, ZdofH)
c
      use data_structure3D
#include "syscom.blk"
c
c
c  ...ATTENTION: output array is dimensioned with the number
c     of components depending upon "Icase", this facilitates
c     the use of the routine in update_dirichlet
      dimension Xi(NDIMEN),
     .          ZdofH(NRCOMS*NREQNH(Icase))
c
c  ...Dirichlet BC data
      dimension zvalH(MAXEQNH),zdvalH(MAXEQNH,3),
     .          zvalE(3,MAXEQNE),
     .          zvalV(3,MAXEQNV)
c
c  ...point coordinates
      dimension x(3),void(3,3)
c
c  ...physics flags decoded
      dimension ncase(NR_PHYSA)
c
c----------------------------------------------------------------------
c
      iprint=0
c
c  ...determine coordinates of the point
      select case(Iflag)
      case(5);        call prism(No,Xi, x,void)
      case(6);        call  hexa(No,Xi, x,void)
      case(7);        call tetra(No,Xi, x,void)
      case(8);        call pyram(No,Xi, x,void)
      end select
c
      if (iprint.eq.1) then
        write(*,8000)Icase,Xi(1:3),x(1:3)
 8000   format(' dhpvert: Icase,Xi,x = ',i2,2x,
     .          3(e12.5,2x),2x,3(e12.5,2x))
      endif
c
      call dirichlet(x,Icase, zvalH,zdvalH,zvalE,zvalV)
c
c  ...shift the data skipping irrelevant entries
      call decod(Icase,2,NR_PHYSA, ncase)
      if (iprint.eq.1) then
        write(*,7002) ncase
 7002   format('dhpvert: ncase = ',5i1)
      endif
c
      ivarH=0;  nvarH=0
c
      if (iprint.eq.1) then
        write(*,*) 'dhpvert: NRCOMS = ',NRCOMS
      endif
c
c  ...loop through multiple copies of variables
      do j=1,NRCOMS
c
c  .....loop through physical attributes
        do i=1,NR_PHYSA
c
c  .......loop through components
          do k=1,NR_COMP(i)
            select case(DTYPE(i))
            case('contin')
              ivarH = ivarH+1
              if (ncase(i).eq.1) then
                nvarH = nvarH+1
                ZdofH(nvarH) = zvalH(ivarH)
              endif
            end select
          enddo
        enddo
      enddo
c
c
      end
