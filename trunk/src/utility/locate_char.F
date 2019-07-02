c----------------------------------------------------------------------
c
c   routine name       - locate_char
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Dec 06
c
c   purpose            - routine locates a character unknown
c                        on a list 
c
cc   arguments :
c     in:
c              
c              List    - a list
c              Nlist   - length of the list
c     out:
c              Number  - the element location
c
c   required  routines - 
c
c----------------------------------------------------------------------
c
      subroutine locate_char(Nel,List,Nlist, Number)
c
      character (len = *) Nel,List(Nlist)
c
      if (Nlist.eq.0) then
        Number=0
        return
      else
        do i=1,Nlist
          if (List(i).eq.Nel) then 
            Number=i
            return
          endif
        enddo
        Number=0
        return
      endif
c
      end






