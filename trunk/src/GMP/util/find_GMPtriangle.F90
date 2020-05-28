!----------------------------------------------------------------------
!   routine name       - find_GMPtriangle
!----------------------------------------------------------------------
!   latest revision    - Dec 07
!
!   purpose            - For given three points, routine finds
!                        the corresponding GMP triangle
!                        connected to the nodes, a brute force
!                        search routine
!
!   arguments
!     in:
!             Nvert    - vertex points
!     out:
!             Nt       = triangle number if the search has been
!                        succesful
!                      = 0 otherwise
!---------------------------------------------------------------------
subroutine find_GMPtriangle(Nvert, Ntfound)
!
      use GMP
#include "syscom.blk"
#include "cinout.blk"
!
      dimension Nvert(3)
!
      iprint = 0
 10   continue
      if (iprint .eq. 1) then
        write(*,7001) Nvert(1:3)
 7001   format('find_GMPtriangle: Nvert = ',3i8)
      endif
      Ntfound = 0
!
!  ...loop through triangles
      do nt = 1, NRTRIAN
        iflag = 0
        do iv = 1, 3
          call locate(TRIANGLES(nt)%VertNo(iv),Nvert,3, ii)
          if (ii .ne. 0)  iflag = iflag + 1
        enddo
        if (iflag .eq. 3) then
          Ntfound = nt
          return
        endif
      enddo
      write(*,7002)
 7002 format('find_GMPtriangle: HAVE NOT FOUND THE TRIANGLE, INITIATING ', &
                               'A BRUTE FORCE SEARCH THROUGH TETS')
      do ntet = 1, NRTETRA
        iflag = 0
        do iv = 1, 4
          call locate(TETRAS(ntet)%VertNo(iv),Nvert,3, ii)
          if (ii .ne. 0)  iflag = iflag + 1
        enddo
        if (iflag .eq. 3) then
          write(*,7003) ntet
 7003     format('find_GMPtriangle: ntet = ',i8)
          call pause
          go to 10
        endif
      enddo
!
end subroutine find_GMPtriangle
