c----------------------------------------------------------------------
c
c   routine name       - select
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - dec, 1993
c
c   purpose            - routine establishes a flag indicating
c                        a quanntity to be displaed in the graphics
c                        routines
c
c   usage              - call select
c
c   arguments          - none
c
c   required  routines -
c
c----------------------------------------------------------------------
c
      subroutine select
c
      use graphmod
      use physics
c
#include "syscom.blk"
c
   10 continue
      write(*,7000)
 7000 format('SELECT THE QUANTITY TO DISPLAY:')
c
c  ...loop through physics attributes
      do iphys=1,NR_PHYSA
        write(*,7001) PHYSA(iphys)
 7001   format('ATTRIBUTE = ',a5)
        select case(DTYPE(iphys))
        case('contin','tangen')
          do i=1,NR_COMP(iphys)
            write(*,7002) i,i*100+1+iphys
 7002 format('APPROXIMATE SOLUTION ',i1,' COMPONENT ....',i3)
            write(*,7003) i,i*100+2+iphys
 7003 format('      EXACT SOLUTION ',i1,' COMPONENT ....',i3)
          enddo
        end select
      enddo
c
      read(*,*) ISELECT
c
c
      end



