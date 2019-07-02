c----------------------------------------------------------------------
c
c   routine name       - find_case
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Dec 07
c
c   purpose            - routine identifies the node case number
c                        based on physics info
c
c   arguments :
c     in:
c          Number      - number of physics attributes for a specific
c                        node
c          Phys        - the corresponding list of the attributes 
c     out:
c          Icase       - node case number
c
c   required  routines - none
c
c----------------------------------------------------------------------
c 
      subroutine find_case(Number,Phys, Icase)
c
      use physics
c
#include "syscom.blk"
c
      character(len=5) :: Phys(Number)
      dimension narray(10)
c
      iprint=0
c
      narray=0
c
      if (iprint.eq.1) then
        write(*,7001) Number, Phys(1:Number)
 7001   format(' find_case: Phys(1:',i1,') = ',5(a5,' ; '))
      endif
c
c  ...loop through physics attributes
      do j=1,NR_PHYSA
        call locate_char(PHYSA(j),Phys,Number, ii)
        if (ii.ne.0) narray(j)=1
      enddo
c
c  ...encode the case number
      call encod(narray,2,NR_PHYSA, Icase)
      if (iprint.eq.1) then
        write(*,7002) Icase
 7002   format(' find_case: Icase = ',i2)
        call pause
      endif
c
c
      end
