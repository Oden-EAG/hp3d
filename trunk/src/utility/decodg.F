c----------------------------------------------------------------------
c
c   routine name       - decodg
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Nov 98
c
c   purpose            - routine decodes an integer into a sequence
c                        of digits
c
c   usage              - call decodg(Nick,Mod,N, Narray)
c
c   arguments :
c     in:
c               Nick   - a nickname to be decoded
c               Mod    - system number
c               N      - number of the digits stored in the nickname
c     out:
c               Narray - the decoded digits
c
c   required  routines - 
c
c----------------------------------------------------------------------
c
      subroutine decodg(Nick,Mod,N, Narray)
c
      dimension Narray(N)
c
      nick1 = Nick    
      do i=1,N
        nick2 = nick1/Mod
        Narray(i) = nick1 - nick2*Mod
        nick1 = nick2
      enddo
c
c
      return
      end


     
        
      

