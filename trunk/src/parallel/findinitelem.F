c-----------------------------------------------------------------------
c
c   routine name       - findinitelem
c
c-----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Jan 08
c
c   purpose            - routine finds initial mesh element 
c                        for given middle node
c
c   usage              - findinitelem
c
c   arguments          - Mdle - middle node (IN)
c                        Iel - initial mesh element (OUT)
c
c   required  routines - 
c
c-----------------------------------------------------------------------
c
      subroutine findinitelem(Mdle,Iel)
      
      use data_structure3D
      
#include "syscom.blk"
c
c  ...temporary version for initial mesh elements only
      Iel = abs(NODES(Mdle)%father)

      end subroutine findinitelem
