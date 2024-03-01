!***===***===***===***===***===***===***===***===***===***===***===***==
! FUNCTION:
! converts nodal dest. vectors to dof dest. vectors
! equations to be eliminated are written to Andest
!    giving current location in front
! if dof i  is making its last appearance then
!    Amdest(i) is lt zero                 !
!**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
! ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
!
! Typ Name      Function
! I   Nd      number of nodal destination vectors
! I   Aldest  nodal destination vectors
! O   Nee     number of dof to eliminate from the front
! O   Amdest  dof destination vectors
! O   Andest  destination in the front where the remaining values
!              in the front will transfer to as each exiting dof is elim
!*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*=*
! LATEST REVISION: Mar 2023
!++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++==++=
! NAMING CONVENTIONS:
!     AAAAAAAA    Variables in COMMON & PARAMETERS
!     Aaaaaaaa    Variables as ARGUMENTS
!     aaaaaaaa    LOCAL Variables
!         7xxx    FORMAT Statements
!         9xxx    ERROR Handling
!+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+
!
      subroutine dest (Nd, Aldest, Nee, Amdest, Andest)
!
      use surfsc1
      use surfsc2
!
      implicit none
!
      integer :: Nd, Nee
      real(8) :: Aldest(*), Amdest(*), Andest(*)
!
      integer :: i,j,km,kn,l,ldest,m,n,ndi,ndl,ne
!
#if HP3D_DEBUG
      integer :: iprint
#endif
!
#if HP3D_DEBUG
      iprint = 0
      if (iprint.eq.1) then
        write(*,*) 'DEST: Nd, Aldest = ',Nd
        write(*,7001) (Aldest(iii),iii=1,Nd)
 7001   format(1x,6f10.1)
        call pause
      endif
#endif
!
! initialize
!
      NFW = 0
      km = 1
      kn = 1
      NDOFM = 0
      ne = 0
!
! create the element dof destination vectors   (Amdest)
! =======================================================
! loop over the number of nodal destination vectors
!
      do 50 i = 1,Nd
!     --------------
! set the nodal dest vec to an integer
! note: destination vecs are packed as:
!
cld,aug,1991 wrong!!!  desvec = (destflag + NICMUL*ndof + 100*frontdest)
!  desvec = (destflag + 10*ndof + NICMUL*10*frontdest)
!
         ldest = int(Aldest(i))
!
! pull the destination flag for this node
!  0 = first occurance of this node
!  1 = intermediate occurance of this node
!  2 = final and only occurance of this node
!  3 = final occurance of this node
!
         m = modr(ldest,10)
!
! pull the number of dof associatted with this node
!cwb >
!cwb          n = modr(ldest,100)/10
         n = modr(ldest,10*NICMUL)/10
!cwb <
! increment the local total of dof for this elem
!
         NDOFM = NDOFM + n
!
! increment the 'number of dof to eliminate' counter
!
         if (m .ge. 2) ne = ne + n
!
! pull the destination in the front (and set to a zero pointer)
!cwb >
!cwb          l = ldest/100 - 1
         l = ldest/(10*NICMUL) - 1
!cwb <
! loop over dof
!
         do 10 j = 1,n
!
! set the front location (ie: the dof destination vector array)
!
            Amdest(km) = l + j
!
            if (m .ge. 2) then
!
! flag dof to be eliminated with a negative
!
               Amdest(km) = -Amdest(km)
!
! set into elimination dof destination vector array (Andest)
!  note: this is changed to its final form below
!
               Andest(kn) = l + j
               kn = kn+1
!
            endif
!
            km = km+1
!
   10    continue
!
! set the current frontwidth number
!cwb >
         l = abs(nint(Amdest(km-1)))
!cwb           l = abs(idnint(Amdest(km-1)))
!cwb <
         if (l .gt. NFW) NFW = l
!
   50 continue
!  -----------
!
! create the elimination destination vector array (Andest)
!   note: Andest points to the destination in the front where
!         the remaining values in the front will transfer to as
!         each exiting dof is eliminated
!
! ie:  Amdest = [-1,-2,3,4,5,6,-7,-8,9,10]
! then Andest = [1,1,5,5]  (ne=4)
! so we end up with : [3,4,5,6,9,10]
!========================================================
! loop over the number of dof to eliminate
!
      do 70 i = 1,ne
!     --------------
         j = i + 1
!
! ALLIANT directives
!cvd$ select (vector)
! ARDENT directives
!$doit VBEST
!
! loop over the remaining number of dof to eliminate
!
         do 60 l = j,ne
!
            ndi = int(Andest(i))
            ndl = int(Andest(l))
!
            if (ndi .lt. ndl) ndl = ndl - 1
            Andest(l) = ndl
!
60       continue
!
   70 continue
!  -----------
      Nee = ne
!      write(*,*) 'in dest Nee = ',Nee
!      pause
!
! debug print
!
      if(IPRDES .ne. 1) return
!cwb >
!cwb       write(NFSOUT,7000) (idnint(Aldest(i)),i=1,nd)
!cwb       write(NFSOUT,7010) (idnint(Amdest(i)),i=1,NDOFM)
!cwb       if (ne .ne. 0) write(NFSOUT,7020) (idnint(Andest(i)),i=1,ne)
      write(NFSOUT,7000) (nint(Aldest(i)),i=1,nd)
      write(NFSOUT,7010) (nint(Amdest(i)),i=1,NDOFM)
      if (ne .ne. 0) write(NFSOUT,7020) (nint(Andest(i)),i=1,ne)
!cwb <
7000    format(/,1x,'IN DEST: NODAL DESTINATION VECTORS',10i7, &
     .        10(/,34x,10i7))
7010    format(11x,'DOF DESTINATION VECTORS',10i7,10(/,35x,10i7))
7020    format(9x,'ELIM. DESTINATION VECTORS',10i7,10(/,35x,10i7))
!
!
      end subroutine dest
