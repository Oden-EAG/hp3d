c----------------------------------------------------------------------
c
c   routine name       - findmin
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Jul 98
c
c   purpose            - find a minimum integer in an array and its
c                        location
c
c   usage              - call findmin(Int_array,N, Min_element,Loc)
c
c   arguments :
c     in:
c          Int_array   - an integer array
c          N           - dimension of the array
c     out:
c          Min_element - minimum element
c          Loc         - location of the element
c
c   required  routines - none
c
c----------------------------------------------------------------------
c
      subroutine findmin(Int_array,N, Min_element,Loc)
c
      dimension Int_array(N)
c
cwr03.15.02
ccccc Min_element = 10**10
      Min_element = 10**8
      Loc=0
      do i=1,N
        if (Min_element.gt.Int_array(i)) then
          Min_element = Int_array(i)
          Loc = i
        endif
      enddo
c
      if (Loc.eq.0) then
        write(*,*) 'findmin: ERROR !!'
        stop 1
      endif
c
c
      return
      end



