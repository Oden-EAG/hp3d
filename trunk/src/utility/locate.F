c----------------------------------------------------------------------
c
c   routine name       - locate
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Jun 07
c
c   purpose            - routine locates an integer on a list
c
c   arguments :
c     in:
c              Nel     - an integer to be located
c              List    - a list of integers
c              Nlist   - length of the list
c     out:
c              Number  - the element location
c
c   required  routines - 
c
c----------------------------------------------------------------------
c
      subroutine locate(Nel,List,Nlist, Number)
c
      dimension List(Nlist)
c
      do i=1,Nlist
        if (List(i).eq.Nel) then 
          Number=i
          return
        endif
      enddo
      Number=0
c
      end



