      subroutine desvec (In, Ia, Ibb, Ic, Ntot)
c
c***===***===***===***===***===***===***===***===***===***===***===***==
c FUNCTION:
c calculation of nodal destination vectors from nicknames
c------------------------------------------------------------------
c
c note: destination vecs are packed as:
c
c        desvec = (destflag + 10*ndof + NICMUL*10*frontdest)
c
c where : destflag  = 0 : first occurance of this node
c                   = 1 : intermediate occurance of this node
c                   = 2 : final and only occurance of this node
c                   = 3 : final occurance of this node
c
c         ndof      = number of dof associated w/ this nickname
c
c         frontdest = destination into the current front
c
c
c**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
c ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
c
c Typ Name      Function
c I   In    array (length numel) containing the number of nodes+multipli
c             essentially the number of nicknames/element
c             **note: this must be pre-built by the calling program
c
c I   Ia    array that contains element nick names
c             nicknames are packed as:
c
c            nick = (ndof + NICMUL*nodenum)
c
c            where :   ndof    = number of dof associated w/ this nickna
c                      nodenum = node number or some other unique identi
c                                 reflecting the connectivity or relatio
c                                 to the system of equations
c
c           **note: this must be pre-built by the calling program
c
c W   Ibb    array for destination vectors
c I   Ntot   total number of degrees of freedom
c
c W   Ic     temporary workspace for calculating the front elimination
c..........................................................
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
c
      use surfsc1
      use surfsc2
#include "syscomz.blk"
      dimension In(*), Ia(*), Ibb(*), Ic(*)
c
#if DEBUG_MODE
      integer :: iprint=0
#endif
c
c internal function: modr
c
      modr(i,j) = i - (i/j)*j
c
c initialize
c
      MDOF = 0
      MFW = 0
      ides = 1
      ip = 0
      jdn = 0
      Ntot=0
c
c initialize the destination vector storage
c
      do 10 i = 1,NFN
         Ibb(i) = 0
   10 continue
#if DEBUG_MODE
      if (iprint.eq.1) write(NFSOUT,*) 'in desvec.f, numelm=',numelm
#endif
c
c loop thru the elements
c
      do 100 iel = 1,NUMELM
c     **********************
c
c if the element contains no nodes or multipliers, then skip it
c
#if DEBUG_MODE
         if (iprint.eq.1) then
           write(*,*) 'desvec: iel = ',iel
         endif
#endif
         n = In(iel)
         if (n .eq. 0) go to 100
c
c initialize
c
         nt = 0
         ips = ip
         ipc = 1
         ne = 0
         ntt = 0
c
c loop thru the nicknames
c
         do 60 id = 1,n
c        ***************
c
#if DEBUG_MODE
            if (iprint.eq.1) then
              write(*,*) 'desvec: id = ',id
            endif
#endif
            ifrst = 0
c
c pull up the nickname
c
            ip = ip + 1
            inic = Ia(ip)
c
c unpack the number of dof associated w/ this nickname
cwb >
#if DEBUG_MODE
            if (iprint.eq.1) then
              write(*,*) 'desvec: inic,NICMUL = ',inic,NICMUL
            endif
#endif
            NDOFM = modr(inic,NICMUL)
#if DEBUG_MODE
            if (iprint.eq.1) then
              write(*,*) 'desvec: NDOFM = ',NDOFM
            endif
#endif
cwb             NDOFM = modr(inic,10)
cwb <
            nt = nt + NDOFM
c
c set the front destination
c---------------------------
c if there is a destvec already associated with this nickname
c   then the front dest is what it was previously
c
            if (Ibb(ip) .gt. 0) then
               jdes = Ibb(ip)
               ifrst = 1
c
c if there is not a destvec associated with this nickname
c  then pick up the next location in the front
c  also if this is greater than the current max frontwidth
c  then update the current max frontwidth
c
            else
               jdes = ides
               ides = ides + NDOFM
               if(ides-1 .gt. mfw) mfw = ides-1
            endif
c
c pack the destination vector (preliminary; destflag is not set)
cwb >
            Ibb(ip) = jdes*10*NICMUL + NDOFM*10
cwb             Ibb(ip) = jdes*100 + NDOFM*10
cwb <
            jp = ip + 1
            if (jp .le. nfn) then
c
c loop thru all the remaining nicknames and search for the next matching
c ----------------------------------------------------------------------
               do 40 jjd = jp,nfn
c              ******************
                  jd = jjd
c
c check for a matching nickname
c
                  if (inic .eq. Ia(jd)) then
c
c a match was found,
c   so set the next occurance to the current destination
c
                     Ibb(jd) = jdes
c
c set the destflag for this dest vec
c  (either the first or intermediate occurance)    (ie: ifrst=0 or 1)
c
                     Ibb(ip) = Ibb(ip) + ifrst
                     if(jd .gt. jdn) jdn = jd
c
c skip to the next nickname
c
                     go to 60
                  endif
c
   40          continue
c  ********************
            endif
c
c no matching nickname was found, so this one is to exit the front
c-----------------------------------------------------------------
c  set the final occurance flag in destvec
c  (ifrst=0 or 1) so we are setting destflag=2 or 3
c
            Ibb(ip) = Ibb(ip) + 2 + ifrst
c
c stick the front dest and number of dof for each elimination node
c  in this element into temp storage
c
            Ic(ipc) = jdes
            Ic(ipc+1) = NDOFM
            ipc = ipc + 2
c
c increment the number of nodes to eliminate
c
            ne = ne + 1
c
c increment the total number of dof to eliminate for this element
c
            ntt = ntt + NDOFM
c
   60    continue
c  ***************
c
c increment the total number of dof for the model
c
         Ntot = Ntot + ntt
c
c check if there is a new max dof/node
c
         if (nt .gt. MDOF) MDOF = nt
c
c if this isnt the last elem and we need to eliminate nodes from the fro
c
         if (iel.lt.NUMELM .and. ne.ne.0) then
c
c set the new total number of number of dof in the front
c
            ides = ides - ntt
            jp = ips + n + 1
c
c for nodes in later elements which we given front destinations
c  we must now change their front dests to reflect that some nodes are e
c-----------------------------------------------------------------------
            if (jp .le. jdn) then
c
c search thru the remaining nicknames
c
               do 80 jd = jp,jdn
c
c check if this nickname was given a front dest already
c
                  if(Ibb(jd) .ne. 0) then
c
                     ipc = 1
                     nt = 0
c
c look thru those front dests which are being eliminated
c
                     do 70 i = 1,ne
c
c if the front dest of this nickname is ge than the one we are eliminati
c  then increment the number of dof by that being eliminated
c
                        if (Ibb(jd) .ge. Ic(ipc))  nt = nt + Ic(ipc+1)
                        ipc = ipc + 2
   70                continue
c
c reduce the front dest of this nickname by the number of dof eliminated
c  below it in the front
c
                     Ibb(jd) = Ibb(jd) - nt
                  endif
c
   80          continue
            endif
         endif
c
c call the user-supplied routine which stores the element destination ve
c  into some global storage scheme
c
         call preout (iel, n, Ia(ips+1), Ibb(ips+1))
c        --------------------------------------------

czka000811         write(30,*)iel,'  =IEL'
czka000811         write(30,*)n
czka000811         write(30,*)(Ibb(ii),ii=ips+1,ips+n)
 1000    format(1x,10i12)

c debug print
c
         if (IPRDES .eq. 1) then
           write(NFSOUT,7701) iel,n,(Ibb(ips+i),i=1,n)
 7701      format(1x,'DESVEC: iel,nnod',2i8,'  (Ibb(ips+i),i=1,n)',
     .            (/,10x,10i8))
           call pause
         endif
c
  100 continue
c
c
      end subroutine desvec
