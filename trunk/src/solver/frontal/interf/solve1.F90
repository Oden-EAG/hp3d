!------------------------------------------------------------------
!
!   routine name       - solve1
!
!------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine solves a user defined system
!                        of equations interfacing with the
!                        frontal solver routines
!
!   arguments
!   in:
!        Number_of_RHS  - number of right-hand sides (load vectors)
!
!------------------------------------------------------------------
#include "typedefs.h"
   subroutine solve1(Number_of_RHS)
!
      use data_structure3D
      use element_data
      use assembly
      use frsolmod
      use control
      use assembly_sc, only: NRDOF_TOT, NRDOF_CON
      use stc,         only: stc_alloc, stc_dealloc, stc_get_nrdof
      use par_mesh,    only: DISTRIBUTED,HOST_MESH
      use mpi_param,   only: RANK,ROOT
!
!  ...frontal solver common blocks
      use surfsc1
      use surfsc2
!
      implicit none
!
      integer, intent(in) :: Number_of_RHS
!
!  ...nodes for a modified element and the corresponding number
!     of H1,H(curl),H(div) and L2 dof
      integer :: nodm(MAXNODM),ndofmH(MAXNODM),ndofmE(MAXNODM), &
                               ndofmV(MAXNODM),ndofmQ(MAXNODM)
!
!  ...number of variables for each physics attribute for an element
      integer :: nrdofs(NR_PHYSA)
!
!  ...element geometry dof, direction vector
      real(8) :: xnod(3,MAXbrickH),xc(3),direction(3)
!
!  ...number of local element dof for each physics variable
      integer :: nrdofi(NR_PHYSA),nrdofb(NR_PHYSA)
!
      integer :: inick,iel,iel1,i,j,iv,kk,kel
      integer :: max_active_node_no,mdle,mr,ms,mu,nod,nrnodm
      integer :: nrnick,nrdof,nrdofm,nrdofc,nrnod1,nrdof_mdl
      VTYPE   :: zvoid
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!
!-----------------------------------------------------------------------
!
      if (RANK .ne. ROOT) return
      if (DISTRIBUTED .and. (.not. HOST_MESH)) then
        write(*,*) 'solve1: mesh is distributed. returning...'
        return
      endif
!
#if HP3D_COMPLEX
      IDUMPWR = 0
#endif
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,*) 'solve1: NRELES,Number_of_RHS = ', &
                            NRELES,Number_of_RHS
        call pause
      endif
#endif
!
!  ...save the number of right-hand sides (load vectors)
      NR_RHS = Number_of_RHS
      if (NR_RHS .ne. NRRHS) then
         write(*,*) 'NR_RHS, NRRHS = ',NR_RHS,NRRHS
         stop
      endif
!
!  ...preliminary loop through elements to collect the necessary data
!     for allocating memory
      nrnick = 0
      allocate(MAXDOFS(NR_PHYSA)); MAXDOFS = 0
!
!  ...use the potential maximum number of dofs for extraction and
!     Dirichlet dof vectors
!
!  ...case with static condensation
      if (ISTC_FLAG) then
        MAXDOFM = (MAXbrickH-MAXmdlbH)*NRHVAR &
                + (MAXbrickE-MAXmdlbE)*NREVAR &
                + (MAXbrickV-MAXmdlbV)*NRVVAR
!
!  ...no static condensation
      else
        MAXDOFM = MAXbrickH*NRHVAR &
                + MAXbrickE*NREVAR &
                + MAXbrickV*NRVVAR &
                + MAXbrickQ*NRQVAR
      endif
!
      allocate(NEXTRACT(MAXDOFM))
      allocate(IDBC(MAXDOFM))
      allocate(ZDOFD(MAXDOFM,NR_RHS))
!
!  ...initiate maximum active node number
      max_active_node_no = 0
!
!  ...determine the actual maximum number of element dof to allocate
!     element matrices
      MAXDOFM = 0; MAXDOFC=0
      nrdof=0; nrdof_mdl=0
      call reset_visit
      do iel=1,NRELES
        mdle = ELEM_ORDER(iel)
!
!  .....determine nodes of the modified element
#if HP3D_DEBUG
        if (iprint.eq.1) then
          write(*,*) 'solve1: iel', iel, 'mdle', mdle
          call pause
        endif
#endif
        if (ISTC_FLAG) then
          call celem_systemI(iel,mdle,1, &
                     nrdofs,nrdofm,nrdofc, &
                     nodm,ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm, &
                     zvoid,zvoid)
        else
          call celem(mdle,1, &
                     nrdofs,nrdofm,nrdofc, &
                     nodm,ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm, &
                     zvoid,zvoid)
        endif
!
!  .....compute the total number of mdle node dof
        if (ISTC_FLAG) then
           call stc_get_nrdof(mdle, nrdofi,nrdofb)
           nrdof_mdl = nrdof_mdl + sum(nrdofb)
        endif
!
!  .....compute the number of node nicknames for the element
        nrnod1=0
        do i=1,nrnodm
          if (ndofmH(i).ne.0) nrnod1 = nrnod1+1
          if (ndofmE(i).ne.0) nrnod1 = nrnod1+1
          if (ndofmV(i).ne.0) nrnod1 = nrnod1+1
          if (ndofmQ(i).ne.0) nrnod1 = nrnod1+1
!
          nod = nodm(i)
          if (NODES(nod)%visit .eq. 0) then
            nrdof = nrdof+ndofmH(i)+ndofmE(i)+ndofmV(i)+ndofmQ(i)
            NODES(nod)%visit = 1
          endif
!
!  .......update the maximum active node number
          max_active_node_no = max(max_active_node_no,nodm(i))
        enddo
        nrnick = nrnick + nrnod1
!
!  .....update the maximum number of local dof
        do i=1,NR_PHYSA
          MAXDOFS(i) = max0(MAXDOFS(i),nrdofs(i))
        enddo
!
!  .....update the maximum number of modified element dof
        MAXDOFM = max0(MAXDOFM,nrdofm)
!
!  .....update the maximum number of modified element dof after
!       compression
        MAXDOFC = max0(MAXDOFC,nrdofc)
!
!  ...end of loop through elements
      enddo
!
      NRDOF_CON = nrdof
      NRDOF_TOT = nrdof + nrdof_mdl
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7006) NR_PHYSA
 7006   format("number of physical attributes, NR_PHYSA = ",i3)
        write(*,7005) MAXDOFS(1:NR_PHYSA)
 7005   format('solve1: MAXDOFS = ',8i10)
        write(*,7001) MAXDOFM,MAXDOFC
 7001   format('solve1: MAXDOFM,MAXDOFC = ',2i10)
        write(*,7007) max_active_node_no
 7007   format('solve1: max_active_node_no = ',i8)
      endif
#endif
!
!-----------------------------------------------------------------------
!
      if (REORDER) then
        if (.not.allocated(NEW_ELEM_ORDER)) then
!
!  .......reorder the elements in order to minimize the bandwidth...
          direction(1:3) = 0.d0; direction(1) = 1.d0;
          write(*,*) 'solve1: REORDERING ELEMENTS ALONG x DIRECTION'
          allocate(NEW_ELEM_ORDER(NRELES))
          allocate(ELEM_CENTER(NRELES))
          do iel=1,NRELES
            mdle = ELEM_ORDER(iel)
            call nodcor(mdle, xnod)
            xc(1:3) = 0.d0
            do iv=1,nvert(NODES(mdle)%ntype)
              xc(1:3) = xc(1:3) + xnod(1:3,iv)
            enddo
            xc(1:3) = xc(1:3)/nvert(NODES(mdle)%ntype)
            call scalar_product(xc,direction, ELEM_CENTER(iel))
          enddo
          call sortm(NRELES,NEW_ELEM_ORDER,ELEM_CENTER)
          do iel=1,NRELES
            iel1=NEW_ELEM_ORDER(iel)
            NEW_ELEM_ORDER(iel) = ELEM_ORDER(iel1)
            write(*,*) 'iel,NEW_ELEM_ORDER(iel) = ', &
                        iel,NEW_ELEM_ORDER(iel)
            if (iel/20*20.eq.iel) call pause
          enddo
          deallocate(ELEM_CENTER)
        endif
      else
        if (.not.allocated(NEW_ELEM_ORDER)) then
          allocate(NEW_ELEM_ORDER(NRELES))
          NEW_ELEM_ORDER(1:NRELES) = ELEM_ORDER(1:NRELES)
        endif
      endif
!
!-----------------------------------------------------------------------
!
!  ...define nicknames.....
      allocate(IN(NRELES))
      allocate(IAWORK(2*nrnick+3*MAXDOFC))
!
!  ...initiate the sequential element number
      kel = 0
!
!  ...initiate the node nickname counter
      kk=0
!
!  ...loop through elements
      do iel=1,NRELES
        mdle = NEW_ELEM_ORDER(iel)
!
!  .....determine nodes of the modified element
        if (ISTC_FLAG) then
          call celem_systemI(iel,mdle,1, &
                     nrdofs,nrdofm,nrdofc, &
                     nodm,ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm, &
                     zvoid,zvoid)

        else
          call celem(mdle,1, &
                     nrdofs,nrdofm,nrdofc, &
                     nodm,ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm, &
                     zvoid,zvoid)
        endif
!
!  .....update the sequential element number
        kel = kel+1
!
!  .....initiate the element counter for the nodes
        inick=0
!
!  .....H1 dof ........
!
!  .....loop through nodes in the reversed order
        do i=nrnodm,1,-1
          if (ndofmH(i).gt.0) then
            inick = inick+1
            kk = kk+1
            IAWORK(kk) = nodm(i)*1000 + ndofmH(i)
            if (IAWORK(kk).le.0) then
              write(*,*) 'solve1: NICKNAMES INTEGER OVERFLOW'; stop 1
            endif
          endif
        enddo
!
!  .....H(curl) dof ........
!
!  .....loop through nodes in the reversed order
        do i=nrnodm,1,-1
          if (ndofmE(i).gt.0) then
            inick = inick+1
            kk = kk+1
            IAWORK(kk) = (nodm(i)+max_active_node_no)*1000 + ndofmE(i)
            if (IAWORK(kk).le.0) then
              write(*,*) 'solve1: NICKNAMES INTEGER OVERFLOW'; stop 1
            endif
          endif
        enddo
!
!  .....H(div) dof ........
!
!  .....loop through nodes in the reversed order
        do i=nrnodm,1,-1
          if (ndofmV(i).gt.0) then
            inick = inick+1
            kk = kk+1
            IAWORK(kk) = (nodm(i)+2*max_active_node_no)*1000 + ndofmV(i)
            if (IAWORK(kk).le.0) then
              write(*,*) 'solve1: NICKNAMES INTEGER OVERFLOW'; stop 1
            endif
          endif
        enddo
!
!  .....L2 dof ........
        if (.not.ISTC_FLAG) then
!
!  .......middle node only
          i=nrnodm
          if (ndofmQ(i).gt.0) then
            inick = inick+1
            kk = kk+1
            IAWORK(kk) = (nodm(i)+3*max_active_node_no)*1000 + ndofmQ(i)
            if (IAWORK(kk).le.0) then
              write(*,*) 'solve1: NICKNAMES INTEGER OVERFLOW'; stop 1
            endif
          endif
        endif
!
        IN(kel) = inick
!
!  ...end of the loop through elements
      enddo
!
!  ...exit if no free dof
      write(*,*) 'solve1: TOTAL NUMBER OF NICKNAMES = ',kk
      if (kk.eq.0) then
        deallocate(MAXDOFS,NEXTRACT,IDBC,ZDOFD,IN,IAWORK,NEW_ELEM_ORDER)
        return
      endif
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
   10   write(*,7002) NRELES
 7002   format(' solve1: SET ELEMENT NUMBER TO REVIEW NICKNAMES,', &
               ' NRELES = ',i8)
        read(*,*) mdle
        if (mdle.eq.0) then
          goto 20
        elseif (mdle.gt.0) then
          call locate(mdle,NEW_ELEM_ORDER,NRELES, iel1)
          write(*,*) 'iel1 = ',iel1
        else
          iel1=-mdle
          mdle = NEW_ELEM_ORDER(iel1)
        endif
        kk=0
        do iel=1,iel1-1
          kk = kk + IN(iel)
        enddo
        write(*,7003) iel1,mdle,S_Type(NODES(mdle)%ntype)
 7003   format('        NICKNAMES FOR iel1,mdle',i8,i10,' TYPE = ',a5)
        do i=1,IN(iel1)
          kk=kk+1
          write(*,*) 'i,nick(i) = ',i,IAWORK(kk)
        enddo
        goto 10
      endif
#endif
!
 20   continue
!
!----------------------------------------------------------------------
!
!  filling up common surfs1 (input for frontal solver)
!
      NFSOUT = 0
!  ...hardwire ISYM to symmetry flag from control file
      ISYM   = ISYM_FLAG
      IRESOL = 0
      NRHS   = NR_RHS
      IWRT   = 0
      IASSEM = 1
      IPRSTR = 0
      IPRPIV = 1
      IPRDES = 0
      IERR   = 0
      IPFSLF = 0
      IPFSLE = 0
      IPFSXX = 0
      IPFSBK = 0
!
!  ...set these parameters to negative sequential element number
!     for debugging forward elimination
      IPFSST = 0
      IPFSLH = 0
      IPFSRH = 0
      IPFSLE = 0
      IPFSLF = 0
      IPFSRF = 0
!
      NICMUL = 1000
!
      NUMELM = NRELES
!
      NPDESV = 1
      allocate(IDESVE(nrnick))
      allocate(NDESVE(2,NRELES))
!
!  ...call prefront
!!!      write(*,*) 'solve1: nrnick,MAXDOFC = ',nrnick,MAXDOFC
      MA = 2*nrnick+3*MAXDOFC
#if HP3D_DEBUG
      if (iprint.eq.1) write(*,*) 'solve1: CALLING PREFRONT'
#endif
      call surfsp(IN,IAWORK, ms,mu,mr)
      if (IERR.eq.1) stop 1
      deallocate(IN, IAWORK)
!
!  ...allocate element matrices
      allocate(BLOC(NR_PHYSA))
      allocate(AAUX(NR_PHYSA))
      allocate(ALOC(NR_PHYSA,NR_PHYSA))
      do i=1,NR_PHYSA
        BLOC(i)%nrow = MAXDOFS(i)
        BLOC(i)%ncol = NR_RHS
        allocate(BLOC(i)%array(MAXDOFS(i),NR_RHS))
        do j=1,NR_PHYSA
          ALOC(i,j)%nrow = MAXDOFS(i)
          ALOC(i,j)%ncol = MAXDOFS(j)
          allocate(ALOC(i,j)%array(MAXDOFS(i),MAXDOFS(j)))
        enddo
        AAUX(i)%nrow = MAXDOFM
        AAUX(i)%ncol = MAXDOFS(i)
        allocate(AAUX(i)%array(MAXDOFM,MAXDOFS(i)))
      enddo
      allocate(ZBMOD(MAXDOFM,NR_RHS))
      allocate(ZAMOD(MAXDOFM,MAXDOFM))
!
!  ...check the size of workspace for surfss:
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (ISYM.eq.1 .or. ISYM.eq.4) then
!!!        write(*,7004) ms+mr
!!! 7004   format('solve1: MINIMUM MEMORY  = ',i10)
        if (ms+mr.gt.MFRSOL) then
          write(*,*)'solve1: MFRSOL TOO SMALL',ms,mr,MFRSOL
          stop 1
        endif
      else
!!!        write(*,7004) mu+mr
        if (mu+mr.gt.MFRSOL) then
          write(*,*)'solve1: MFRSOL TOO SMALL',mu,mr,MFRSOL
          stop 1
        endif
      endif
!
!  ...allocate work space
#if HP3D_DEBUG
      if (iprint.eq.1) write(*,*) 'solve1: MFRSOL = ',MFRSOL
#endif
      allocate(ZWORKFRS(MFRSOL))
!
!  ...allocate static condensation data structures
      if (ISTC_FLAG) call stc_alloc
!
!  ...call frontal solver
      MA = MFRSOL
      ZWORKFRS(1:MFRSOL) = ZERO
#if HP3D_DEBUG
      if (iprint.eq.1) write(*,*) 'solve1: CALLING FRONTAL SOLVER'
#endif
!
      call surfss(ZWORKFRS)
      if (IERR.ne.0) then
        write(*,*) 'ERROR IN SURFSS'
        stop 1
      endif
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,*) 'solve1: SURFSS COMPLETED'
      endif
#endif
!
!  ...deallocate static condensation data structures
      if (ISTC_FLAG) call stc_dealloc
!
!  ...deallocate ALL arrays used by frontal solver
      do i=1,NR_PHYSA
        if (allocated(BLOC(i)%array)) deallocate(BLOC(i)%array)
        do j=1,NR_PHYSA
          if (allocated(ALOC(i,j)%array)) deallocate(ALOC(i,j)%array)
        enddo
        if (allocated(AAUX(i)%array)) deallocate(AAUX(i)%array)
      enddo
      deallocate(BLOC,AAUX,ALOC,ZBMOD,ZAMOD,NEXTRACT,IDBC,ZDOFD,MAXDOFS)
      deallocate(IDESVE,NDESVE)
      deallocate(ZWORKFRS)
      deallocate(NEW_ELEM_ORDER)
!
#if HP3D_DEBUG
      if (iprint.eq.1) write(*,*) 'solve1: EXIT ROUTINE SOLVE1'
#endif
!
   end subroutine solve1
