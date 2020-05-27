      subroutine clfil
c              
c-----------------------------------------------------------------------c
c   routine name       - clfil
c
c-----------------------------------------------------------------------c
c   computer           - machine independent
c
c   latest revision    - mar, 91
c
c   purpose            - routine closes files
c
c   usage              - call clfil
c
c   arguments          - none
c
c   required  routines - none
c
c-----------------------------------------------------------------------c
#include "syscom.blk"
#include "cinout.blk"
c      
      close(NIN)
c
      close(NOUT)
c      
      close(N_GEO)
c
      close(NPROFB)
c
c      
      return
      end
