      subroutine dest (Nd, Aldest ,Nee, Amdest, Andest)
c
#include "syscom.blk"
c***===***===***===***===***===***===***===***===***===***===***===***==
c FUNCTION:
c converts nodal dest. vectors to dof dest. vectors
c equations to be eliminated are written to Andest
c    giving current location in front
c if dof i  is making its last appearance then
c    Amdest(i) is lt zero                 c
c**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
c ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
c
c Typ Name      Function
c I   Nd      number of nodal destination vectors
c I   Aldest  nodal destination vectors
c O   Nee     number of dof to eliminate from the front
c O   Amdest  dof destination vectors
c O   Andest  destination in the front where the remaining values
c              in the front will transfer to as each exiting dof is elim
c*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
c LATEST REVISION: Mar 2023
c++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++=
c NAMING CONVENTIONS:
c     AAAAAAAA    Variables in COMMON & PARAMETERS
c     Aaaaaaaa    Variables as ARGUMENTS
c     aaaaaaaa    LOCAL Variables
c         7xxx    FORMAT Statements
c         9xxx    ERROR Handling
c+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
c
#include "surfsc1.blk"
#include "surfsc2.blk"
c
      real(8) :: Aldest(*), Amdest(*), Andest(*)
c
#if DEBUG_MODE
      integer :: iprint=0
#endif
c
c internal function : modr
c ------------------------
c
      modr(i,j) = i - (i/j)*j
c
#if DEBUG_MODE
      if (iprint.eq.1) then
        write(*,*) 'DEST: Nd, Aldest = ',Nd
        write(*,7001) (Aldest(iii),iii=1,Nd)
 7001   format(1x,6f10.1)
        call pause
      endif
#endif
c
c initialize
c
      NFW = 0
      km = 1
      kn = 1
      NDOFM = 0
      ne = 0
c
c create the element dof destination vectors   (Amdest)
c =======================================================
c loop over the number of nodal destination vectors
c
      do 50 i = 1,Nd
c     --------------
c set the nodal dest vec to an integer
c note: destination vecs are packed as:
c
cld,aug,1991 wrong!!!  desvec = (destflag + NICMUL*ndof + 100*frontdest)
c  desvec = (destflag + 10*ndof + NICMUL*10*frontdest)
c
         ldest = Aldest(i)
c
c pull the destination flag for this node
c  0 = first occurance of this node
c  1 = intermediate occurance of this node
c  2 = final and only occurance of this node
c  3 = final occurance of this node
c
         m = modr(ldest,10)
c
c pull the number of dof associatted with this node
cwb >
cwb          n = modr(ldest,100)/10
         n = modr(ldest,10*NICMUL)/10
cwb <
c increment the local total of dof for this elem
c
         NDOFM = NDOFM + n
c
c increment the 'number of dof to eliminate' counter
c
         if (m .ge. 2) ne = ne + n
c
c pull the destination in the front (and set to a zero pointer)
cwb >
cwb          l = ldest/100 - 1
         l = ldest/(10*NICMUL) - 1
cwb <
c loop over dof
c
         do 10 j = 1,n
c
c set the front location (ie: the dof destination vector array)
c
            Amdest(km) = l + j
c
            if (m .ge. 2) then
c
c flag dof to be eliminated with a negative
c
               Amdest(km) = -Amdest(km)
c
c set into elimination dof destination vector array (Andest)
c  note: this is changed to its final form below
c
               Andest(kn) = l + j
               kn = kn+1
c
            endif
c
            km = km+1
c
   10    continue
c
c set the current frontwidth number
cwb >
         l = abs(nint(Amdest(km-1)))
cwb           l = abs(idnint(Amdest(km-1)))
cwb <
         if (l .gt. NFW) NFW = l
c
   50 continue
c  -----------
c
c create the elimination destination vector array (Andest)
c   note: Andest points to the destination in the front where
c         the remaining values in the front will transfer to as
c         each exiting dof is eliminated
c
c ie:  Amdest = [-1,-2,3,4,5,6,-7,-8,9,10]
c then Andest = [1,1,5,5]  (ne=4)
c so we end up with : [3,4,5,6,9,10]
c========================================================
c loop over the number of dof to eliminate
c
      do 70 i = 1,ne
c     --------------
         j = i + 1
c
c ALLIANT directives
cvd$ select (vector)
c ARDENT directives
c$doit VBEST
c
c loop over the remaining number of dof to eliminate
c
         do 60 l = j,ne
c
            ndi = Andest(i)
            ndl = Andest(l)
c
            if (ndi .lt. ndl) ndl = ndl - 1
            Andest(l) = ndl
c
60       continue
c
   70 continue
c  -----------
      Nee = ne
c      write(*,*) 'in dest Nee = ',Nee
c      pause
c
c debug print
c
      if(IPRDES .ne. 1) return
cwb >
cwb       write(NFSOUT,7000) (idnint(Aldest(i)),i=1,nd)
cwb       write(NFSOUT,7010) (idnint(Amdest(i)),i=1,NDOFM)
cwb       if (ne .ne. 0) write(NFSOUT,7020) (idnint(Andest(i)),i=1,ne)
      write(NFSOUT,7000) (nint(Aldest(i)),i=1,nd)
      write(NFSOUT,7010) (nint(Amdest(i)),i=1,NDOFM)
      if (ne .ne. 0) write(NFSOUT,7020) (nint(Andest(i)),i=1,ne)
cwb <
7000    format(/,1x,'IN DEST: NODAL DESTINATION VECTORS',10i7,
     .        10(/,34x,10i7))
7010    format(11x,'DOF DESTINATION VECTORS',10i7,10(/,35x,10i7))
7020    format(9x,'ELIM. DESTINATION VECTORS',10i7,10(/,35x,10i7))
c
c
      end subroutine dest
