      subroutine findan(nel, nelanc, ndummy)

      use data_structure3D

#include "syscom.blk"

      nelanc = nel


      return
c==========================================================
      nfather = NODES(nel)%father

      ndummy = 0
   
      do while(nfather.gt.0)
         nelanc = nfather
         nfather = NODES(nfather)%father
         ndummy = ndummy + 1
         enddo

c      write(*,*)'FINDAN: NEL,NELANC=',nel,nelanc
c      call pause

      end
