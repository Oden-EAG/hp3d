!***===***===***===***===***===***===***===***===***===***===***===***==
! FUNCTION:
! calculation of nodal destination vectors from nicknames
!------------------------------------------------------------------
!
! note: destination vecs are packed as:
!
!        desvec = (destflag + 10*ndof + NICMUL*10*frontdest)
!
! where : destflag  = 0 : first occurance of this node
!                   = 1 : intermediate occurance of this node
!                   = 2 : final and only occurance of this node
!                   = 3 : final occurance of this node
!
!         ndof      = number of dof associated w/ this nickname
!
!         frontdest = destination into the current front
!
!
!**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**==**=
! ARGUMENTS:  (I : input, O : output, IO : input & output, W : workspace
!
! Typ Name      Function
! I   In    array (length NUMELM) containing the number of nodes+multipli
!             essentially the number of nicknames/element
!             **note: this must be pre-built by the calling program
!
! I   Ia    array that contains element nick names
!             nicknames are packed as:
!
!            nick = (ndof + NICMUL*nodenum)
!
!            where :   ndof    = number of dof associated w/ this nickna
!                      nodenum = node number or some other unique identi
!                                 reflecting the connectivity or relatio
!                                 to the system of equations
!
!           **note: this must be pre-built by the calling program
!
! W   Ibb    array for destination vectors
! I   Ntot   total number of degrees of freedom
!
! W   Ic     temporary workspace for calculating the front elimination
!..........................................................
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
   subroutine desvec (In, Ia, Ibb, Ic, Ntot)
!
      use surfsc1
      use surfsc2
!
      implicit none
!
      integer :: In(*), Ia(*), Ibb(*), Ic(*)
      integer :: Ntot
!
      integer :: id,ides,i,iel,ifrst,inic,ip,ipc,ips
      integer :: jd,jdes,jdn,jjd,jp,n,ne,nt,ntt
!
#if HP3D_DEBUG
      integer :: iprint
#endif
!
! initialize
!
      MDOF = 0
      MFW = 0
      ides = 1
      ip = 0
      jdn = 0
      Ntot=0
!
! initialize the destination vector storage
!
      do i = 1,NFN
         Ibb(i) = 0
      enddo
#if HP3D_DEBUG
      iprint = 0
      if (iprint.eq.1) write(NFSOUT,*) 'in desvec: numelm=',numelm
#endif
!
! loop thru the elements
!
      do iel = 1,NUMELM
!     **********************
!
! if the element contains no nodes or multipliers, then skip it
!
#if HP3D_DEBUG
         if (iprint.eq.1) then
           write(*,*) 'desvec: iel = ',iel
         endif
#endif
         n = In(iel)
         if (n .eq. 0) cycle
!
! initialize
!
         nt = 0
         ips = ip
         ipc = 1
         ne = 0
         ntt = 0
!
! loop thru the nicknames
!
         do id = 1,n
!        ***************
!
#if HP3D_DEBUG
            if (iprint.eq.1) then
              write(*,*) 'desvec: id = ',id
            endif
#endif
            ifrst = 0
!
! pull up the nickname
!
            ip = ip + 1
            inic = Ia(ip)
!
! unpack the number of dof associated w/ this nickname
!wb >
#if HP3D_DEBUG
            if (iprint.eq.1) then
              write(*,*) 'desvec: inic,NICMUL = ',inic,NICMUL
            endif
#endif
            NDOFM = modr(inic,NICMUL)
#if HP3D_DEBUG
            if (iprint.eq.1) then
              write(*,*) 'desvec: NDOFM = ',NDOFM
            endif
#endif
!wb             NDOFM = modr(inic,10)
!wb <
            nt = nt + NDOFM
!
! set the front destination
!---------------------------
! if there is a destvec already associated with this nickname
!   then the front dest is what it was previously
!
            if (Ibb(ip) .gt. 0) then
               jdes = Ibb(ip)
               ifrst = 1
!
! if there is not a destvec associated with this nickname
!  then pick up the next location in the front
!  also if this is greater than the current max frontwidth
!  then update the current max frontwidth
!
            else
               jdes = ides
               ides = ides + NDOFM
               if(ides-1 .gt. mfw) mfw = ides-1
            endif
!
! pack the destination vector (preliminary; destflag is not set)
!wb >
            Ibb(ip) = jdes*10*NICMUL + NDOFM*10
!wb             Ibb(ip) = jdes*100 + NDOFM*10
!wb <
            jp = ip + 1
            if (jp .le. nfn) then
!
! loop thru all the remaining nicknames and search for the next matching
! ----------------------------------------------------------------------
               do jjd = jp,nfn
!              ******************
                  jd = jjd
!
! check for a matching nickname
!
                  if (inic .eq. Ia(jd)) then
!
! a match was found,
!   so set the next occurance to the current destination
!
                     Ibb(jd) = jdes
!
! set the destflag for this dest vec
!  (either the first or intermediate occurance)    (ie: ifrst=0 or 1)
!
                     Ibb(ip) = Ibb(ip) + ifrst
                     if(jd .gt. jdn) jdn = jd
!
! skip to the next nickname
!
                     goto 60
                  endif
!
               enddo
!  ********************
            endif
!
! no matching nickname was found, so this one is to exit the front
!-----------------------------------------------------------------
!  set the final occurance flag in destvec
!  (ifrst=0 or 1) so we are setting destflag=2 or 3
!
            Ibb(ip) = Ibb(ip) + 2 + ifrst
!
! stick the front dest and number of dof for each elimination node
!  in this element into temp storage
!
            Ic(ipc) = jdes
            Ic(ipc+1) = NDOFM
            ipc = ipc + 2
!
! increment the number of nodes to eliminate
!
            ne = ne + 1
!
! increment the total number of dof to eliminate for this element
!
            ntt = ntt + NDOFM
!
   60    enddo
!  ***************
!
! increment the total number of dof for the model
!
         Ntot = Ntot + ntt
!
! check if there is a new max dof/node
!
         if (nt .gt. MDOF) MDOF = nt
!
! if this isnt the last elem and we need to eliminate nodes from the fro
!
         if (iel.lt.NUMELM .and. ne.ne.0) then
!
! set the new total number of number of dof in the front
!
            ides = ides - ntt
            jp = ips + n + 1
!
! for nodes in later elements which we given front destinations
!  we must now change their front dests to reflect that some nodes are e
!-----------------------------------------------------------------------
            if (jp .le. jdn) then
!
! search thru the remaining nicknames
!
               do jd = jp,jdn
!
! check if this nickname was given a front dest already
!
                  if(Ibb(jd) .ne. 0) then
!
                     ipc = 1
                     nt = 0
!
! look thru those front dests which are being eliminated
!
                     do i = 1,ne
!
! if the front dest of this nickname is ge than the one we are eliminati
!  then increment the number of dof by that being eliminated
!
                        if (Ibb(jd) .ge. Ic(ipc))  nt = nt + Ic(ipc+1)
                        ipc = ipc + 2
                     enddo
!
! reduce the front dest of this nickname by the number of dof eliminated
!  below it in the front
!
                     Ibb(jd) = Ibb(jd) - nt
                  endif
!
               enddo
            endif
         endif
!
! call the user-supplied routine which stores the element destination ve
!  into some global storage scheme
!
         call preout (iel, n, Ia(ips+1), Ibb(ips+1))
!        --------------------------------------------

!zka000811         write(30,*)iel,'  =IEL'
!zka000811         write(30,*)n
!zka000811         write(30,*)(Ibb(ii),ii=ips+1,ips+n)
 1000    format(1x,10i12)

! debug print
!
         if (IPRDES .eq. 1) then
           write(NFSOUT,7701) iel,n,(Ibb(ips+i),i=1,n)
 7701      format(1x,'DESVEC: iel,nnod',2i8,'  (Ibb(ips+i),i=1,n)', &
                  (/,10x,10i8))
           call pause
         endif
!
      enddo
!
!
   end subroutine desvec
