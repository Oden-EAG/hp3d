c-----------------------------------------------------------------------
c
c   routine name       - decode
c
c-----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - May 20
c
c   purpose            - routine decodes a nickname consisting of
c                        2 digits
c
c   arguments :
c     in:       Nick   - a nickname
c     out:
c               J1     - the first digit  of the nickname
c               J2     - the second digit of the nickname
c
c   required  routines - none
c
c-----------------------------------------------------------------------
c
      subroutine decode(Nick, J1,J2)
      implicit none
      integer :: Nick, J1, J2
c
      J1 = Nick/10
      J2 = Nick - J1*10
c
      end subroutine decode
c
c-----------------------------------------------------------------------
c
c   routine name       - ddecode
c
c-----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - May 20
c
c   purpose            - routine decodes a nickname consisting of
c                        3 digits
c
c   arguments :
c     in:       Nick   - a nickname
c     out:
c               J1     - the first digit  of the nickname
c               J2     - the second digit of the nickname
c               J3     - the third digit of the nickname
c
c   required  routines - none
c
c-----------------------------------------------------------------------
c
      subroutine ddecode(Nick, J1,J2,J3)
      implicit none
      integer :: Nick, J1, J2,J3
c
      integer :: nick2
c
      nick2 = Nick/10
      J3 = Nick - nick2*10
c
      J1 = nick2/10
      J2 = nick2 - J1*10
c
      end subroutine ddecode
