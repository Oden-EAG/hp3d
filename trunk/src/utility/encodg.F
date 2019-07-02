c----------------------------------------------------------------------
c
c   routine name       - encodg
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Nov 98
c
c   purpose            - routine encodes a sequence of digits
c                        into a nickname
c
c   usage              - call encodg(Narray,Mod,N, Nick)
c
c   arguments :
c     in:
c               Narray - digits to be encoded
c               Mod    - system number
c               N      - number of the digits stored in the nickname
c     out:
c               Nick   - the nickname 
c
c   required  routines - 
c
c----------------------------------------------------------------------
c
      subroutine encodg(Narray,Mod,N, Nick)
c
      dimension Narray(N)
c
      Nick=Narray(N)
      do i=N-1,1,-1
        Nick = Nick*Mod + Narray(i)
      enddo
c
c
      end subroutine encodg
