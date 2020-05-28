c----------------------------------------------------------------------
c
c   routine name       - decod
c
c----------------------------------------------------------------------
c
c   latest revision    - Aug 07
c
c   purpose            - routine decodes an integer into a sequence
c                        of digits organized from left to right
c
c   arguments :
c     in:
c               Nick   - a nickname to be decoded
c               Mod    - system number
c               N      - number of the digits stored in the nickname
c     out:
c               Narray - the decoded digits
c
c----------------------------------------------------------------------
c
      subroutine decod(Nick,Mod,N, Narray)
c
      implicit none
c
      integer,             intent(in)    :: Nick,Mod,N
      integer,dimension(N),intent(inout) :: Narray
c
c  ...local variables
      integer :: i,nick1,nick2
c
      nick1 = Nick
      do i=1,N-1
        nick2 = nick1/Mod
        Narray(N+1-i) = nick1 - nick2*Mod
        nick1 = nick2
      enddo
      Narray(1)=nick1
c
c
      end subroutine decod
c
c
c----------------------------------------------------------------------
c
c   routine name       - decodLong
c
c----------------------------------------------------------------------
c
c   latest revision    - Jun 15
c
c   purpose            - routine decodes a LONG integer into a sequence
c                        of digits organized from left to right
c
c   arguments :
c     in:
c               Nick   - a nickname to be decoded (LONG integer)
c               Mod    - system number
c               N      - number of the digits stored in the nickname
c     out:
c               Narray - the decoded digits
c
c----------------------------------------------------------------------
c
      subroutine decodLong(Nick,Mod,N, Narray)
c
      implicit none
c
      integer(8),          intent(in)  :: Nick
      integer,             intent(in)  :: Mod,N
      integer,dimension(N),intent(out) :: Narray
c
c  ...local variables
      integer    :: i
      integer(8) :: nick1,nick2,longmod
c
c  ...be careful with integer type - need to convert!
      longmod = int(Mod,8)
      nick1 = Nick
      do i=1,N
        nick2 = nick1/longmod
c    ...be careful with integer type - need to convert!
        Narray(N+1-i) = int(nick1 - nick2*longmod,4)
        nick1 = nick2
      enddo
c
c
      end subroutine decodLong
