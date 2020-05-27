      subroutine decode2(Nick, Nick1,Nick2)
c
c-----------------------------------------------------------------------c
c   routine name       - decode2
c
c-----------------------------------------------------------------------c
c   computer           - machine independent
c
c   latest revision    - dec, 1991
c
c   purpose            - routine decodes a nickname
c
c   usage              - call decode2(Nick, Nick1,Nick2)
c
c   arguments :
c     in:       Nick   - a nickname
c     out:
c               Nick1  - the first part of the nickname
c               Nick2  - the second two-digit part of the nickname
c
c   required  routines - none
c
c----------------------------------------------------------------------
c
      Nick1 = Nick/100
      Nick2 = Nick - Nick1*100
c
      return
      end
