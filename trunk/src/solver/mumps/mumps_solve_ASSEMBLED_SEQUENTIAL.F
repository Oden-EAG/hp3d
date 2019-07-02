c-----------------------------------------------------------------------
c
c   routine name       - mumps_solve_seq
c
c-----------------------------------------------------------------------
c
c   latest revision    - Mar 11
c
c   purpose            - routine assembles the stiffness matrix,
c                        generates the interface with MUMPS, solves the
c                        resulting matrix, and stores the solution in
c                        module data_structure3D
c
c   arguments          - none
c
c-----------------------------------------------------------------------
c
      subroutine mumps_solve_seq(Number_of_RHS)
c
      use data_structure3D
      use element_data
      use assembly
      use frsolmod
      use mumps_data
      use physics

      use control , only : ISYM_FLAG, SYMMETRIC, UNSYMMETRIC
      use parallelism
#include "syscom.blk"
#include "surfsc1.blk"
#include "mpif.h"
c
      common /ccelem/ iprint_celem
c
c  ...nodes for a modified element and the corresponding number
c     of H1,H(curl),H(div) and L2 dof
      dimension nodm(MAXNODM),ndofmH(MAXNODM),ndofmE(MAXNODM),
     .                        ndofmV(MAXNODM),ndofmQ(MAXNODM)
c
c  ...number of variables for each physics attribute for an element
      dimension nrdofs(NR_PHYSA)
c
c  ...decoded index for a node
      dimension index(NRINDEX)
c
c  ...element geometry dof, direction vector
      dimension xnod(3,MAXbrickH),xc(3),direction(3)
c
c  ...modified symmetric (or not) element stiffness matrix and load vector
      allocatable :: zstiff(:),zxload(:)
c
c  ...multiple right hand sides
      allocatable :: zrhs_for_MUMPS(:,:)
c
c  ...bijection between local and global dof
      allocatable :: iglobal(:)
c
c  ...component counters for the nodes (use in case of multiple loads)
      dimension mvarH(MAXNODM),mvarE(MAXNODM),
     .          mvarV(MAXNODM),mvarQ(MAXNODM)
c
      integer(kind=8) :: t1,t2,clock_rate,clock_max
c
c
c-----------------------------------------------------------------------
c     STEP I: INITIALIZE DATA AND DEFINE MUMPS_PAR%N
c-----------------------------------------------------------------------
c
      iprint = 0
      idebug = 0
c  ...monitors calls to mumps
      iprint_celem = 0
c  ...prints system of equations, row by row
      iprint_system = 0
c  ...compares result with actual result stored in NODES
c       comming e.g. from the frontal solver
      icompare = 0
c  ...test rhs
      itestrhs=0
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call start_clock(iclock_init)
c
c  ...it is symmetric (or maybe not? don't hard-wire, for Christ sake!)
ccc      ISYM=1
      ISYM = ISYM_FLAG
c
c  ...save the number of right-hand sides (load vectors)
      NR_RHS = Number_of_RHS
      NRHS = Number_of_RHS
c
c  ...start mumps and allocate memory for MUMPS
c     DO NOT remove start_mumps....
      call start_mumps
      call alloc_mumps
c
c-----------------------------------------------------------------------
c
c  ...use the potential maximum number of dofs for extraction and
c     Dirichlet dof vectors
      MAXDOFM = MAXbrickH*NRHVAR + MAXbrickE*NREVAR
     .        + MAXbrickV*NRVVAR + MAXbrickQ*NRQVAR
      allocate(NEXTRACT(MAXDOFM))
      allocate(IDBC(MAXDOFM))
      allocate(ZDOFD(MAXDOFM,NR_RHS))
c
c set:
c MUMPS_PAR%N = total number of dof
c
c  ...only the host process assembles data for MUMPS
      if (MUMPS_PAR%MYID == 0) then
c
c  ...set values of N and NELT
c TO BE UPDATED LATER (THIS IS NOT COMPRESSED)
      MUMPS_PAR%N    = NRDOFSH + NRDOFSE + NRDOFSV + NRDOFSQ
c
c----------------------------------------------------------------------
c     STEP II: COMPUTE AN UPPER BOUND OF THE NUMBER OF NONZERO
c              ENTRIES AND THE MAXIMUM NUMBER OF DOF PER ELEMENT
c----------------------------------------------------------------------
c compute:
c   nr_nonzero = total number of non-zero entries
c     (over all elements, with repetitions)
c   MAXDOFM/MAXDOFC = maximum number of dof per element
c     (max. over elements) before and after compression

c  ...initialize the number of nonzero entries
      nr_nonzero=0 ; MAXDOFC=0 ; MAXDOFM=0

c  ...maximum (over all elems) nr of dof for each physics
      allocate(MAXDOFS(NR_PHYSA)); MAXDOFS = 0

      MAXDOFM = 0; MAXDOFC=0
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
c
c  .....determine nodes of the modified element
        call celem(mdle,1,
     .             nrdofs,nrdofm,nrdofc,
     .             nodm,ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,
     .             zvoid,zvoid)
c
c  .....update the maximum number of local dof
        do i=1,NR_PHYSA
          MAXDOFS(i) = max0(MAXDOFS(i),nrdofs(i))
        enddo
c
c  .....update the maximum number of modified element dof
        MAXDOFM = max0(MAXDOFM,nrdofm)
c
c  .....update the maximum number of modified element dof after
c       compression
        MAXDOFC = max0(MAXDOFC,nrdofc)
c
        nr_nonzero = nr_nonzero + nrdofc**2
c
c  ...end of loop through elements
      enddo
c
c
c-----------------------------------------------------------------------
c     STEP III: ALLOCATE MEMORY
c-----------------------------------------------------------------------
c
c IRN,JCN,A are allocated with nr_nonzero
c RHS with %N
c zstiff, zxload with MAXDOFM
c ZALOC, ZBLOC with MAXDOFC
c
c  ...allocate memory
      allocate(MUMPS_PAR%IRN(nr_nonzero),        STAT=i1)
      allocate(MUMPS_PAR%JCN(nr_nonzero),        STAT=i2)
      allocate(MUMPS_PAR%A(nr_nonzero),          STAT=i3)
      allocate(MUMPS_PAR%RHS(MUMPS_PAR%N),       STAT=i4)
      if ((ISYM_FLAG.eq.SYMMETRIC)) then
        allocate(zstiff(MAXDOFM*(MAXDOFM+1)/2),      STAT=i5)
      else
        allocate(zstiff(MAXDOFM*MAXDOFM),            STAT=i5)
      endif
      allocate(zxload(MAXDOFM*NRHS),               STAT=i6)
      if (i1+i2+i3+i4+i5+i6 /= 0) then
        write(*,*)'mumps_solve: FAILED SECOND ALLOCATION!'
        write(*,*)'i1-8 = ',i1,i2,i3,i4,i5,i6,i7,i8
        write(*,*)'MAXDOFM,MAXDOFC=',MAXDOFM,MAXDOFC
        stop
      endif
c
c  ...allocate element matrices
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
c
      write(*,*)'mumps_solve_seq: ALOCs, BLOCs allocated'
c
c  ...only the host process assembles data for MUMPS
      endif
c
#if PARALLEL_MODE
      nout = 0;
      if(MYRANK.eq.0)n_out = MUMPS_PAR%N
      call mpi_bcast(n_out,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(idebug.eq.1)then
        write(*,*)MYRANK,'allocate(zrhs_for_MUMPS(',n_out,NRHS,')'
      endif
      allocate(zrhs_for_MUMPS(n_out,NRHS), STAT=i9)
#else
      allocate(zrhs_for_MUMPS(MUMPS_PAR%N,NRHS), STAT=i9)
#endif
      if (i9 /= 0) then
        write(*,*)'mumps_solve_seq: FAILED SECOND ALLOCATION!'
        write(*,*)'i9',i9
        stop
      endif
c
c  ...only the host process assembles data for MUMPS
      if (MUMPS_PAR%MYID == 0) then
c
c  ...initialize the right hand side
      zrhs_for_MUMPS=ZERO
c
c  ...save the maximum number of dof
      MAXDOF = MAXDOFC
c
c  ...only the host process assembles data for MUMPS
      endif
c
      call mpi_bcast(MAXDOF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
c
      allocate(iglobal(MAXDOF),STAT=i1)
      if(i1/=0)then
        write(*,*)'mumps_solve_seq: Error allocation iglobal'
        stop
      endif
c
c-----------------------------------------------------------------------
c
c  ...reorder the elements in order to minimize the bandwidth...
c     OLD/NEW ELEM_ORDER are defined in frsolmod
      if(.not.allocated(NEW_ELEM_ORDER))then
        direction(1:3) = 0.d0; direction(2) = 1.d0;
        allocate(OLD_ELEM_ORDER(NRELES))
        allocate(NEW_ELEM_ORDER(NRELES))
        allocate(ELEM_CENTER(NRELES))
        mdle=0
        do iel=1,NRELES
          call nelcon(mdle, mdle)
          OLD_ELEM_ORDER(iel)=mdle
          call nodcor(mdle, xnod)
          xc(1:3) = 0.d0
          do iv=1,nvert(NODES(mdle)%type)
            xc(1:3) = xc(1:3) + xnod(1:3,iv)
          enddo
          xc(1:3) = xc(1:3)/nvert(NODES(mdle)%type)
          call scalar_product(xc,direction, ELEM_CENTER(iel))
        enddo
        call sortm(NRELES,NEW_ELEM_ORDER,ELEM_CENTER)
        do iel=1,NRELES
          iel1=NEW_ELEM_ORDER(iel)
          NEW_ELEM_ORDER(iel) = OLD_ELEM_ORDER(iel1)
ccc        write(*,*) 'iel,NEW_ELEM_ORDER(iel) = ',iel,NEW_ELEM_ORDER(iel)
ccc        if (NEW_ELEM_ORDER(iel).eq.1200) stop 1
ccc        if (iel/20*20.eq.iel) call pause
        enddo
        deallocate(OLD_ELEM_ORDER,ELEM_CENTER)
        write(*,*)'mumps_solve_seq: ELEMS sorted'
      else
        write(*,*)'mumps_solve_seq: ELEMS sorted previously'
      endif
c
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c     STEP IV: ASSEMBLE THE STIFFNESS MATRIX AND THE LOAD VECTOR
c-----------------------------------------------------------------------
c
c  ...allocate data structure for storing inverse maps
c-> MP 21.4.2014
c     allocate(NODE_MUMPS_DATA(NRNODS),STAT=i1)
      allocate(NODE_MUMPS_DATA_H1(NRNODS),STAT=i1)
      allocate(NODE_MUMPS_DATA_Hcurl(NRNODS),STAT=i2)
      allocate(NODE_MUMPS_DATA_Hdiv(NRNODS),STAT=i3)
      allocate(NODE_MUMPS_DATA_L2(NRNODS),STAT=i4)
c<- MP 21.4.2014
      if(i1+i2+i3+i4/=0)then
        write(*,*)'NODE_MUMPS_DATA_H1 allocation error',i1
        write(*,*)'NODE_MUMPS_DATA_Hcurl allocation error',i2
        write(*,*)'NODE_MUMPS_DATA_Hdiv allocation error',i3
        write(*,*)'NODE_MUMPS_DATA_L2 allocation error',i4
        stop
      endif
c
c  ...the map must be found on all processors
c  ...initialize visited flags
      do inode=1,NRNODS
c     -> MP 21.4.2014
c        NODE_MUMPS_DATA(inode)%visited=0
        NODE_MUMPS_DATA_H1(inode)%visited=0
        NODE_MUMPS_DATA_Hcurl(inode)%visited=0
        NODE_MUMPS_DATA_Hdiv(inode)%visited=0
        NODE_MUMPS_DATA_L2(inode)%visited=0
c<- MP 21.4.2014
      enddo
c
      write(*,*)'mumps_solve_seq: NRELES = ',NRELES
c
c  ...initialize global dof counter
      iglobaldof=0 ; inzero = 0
      iprint_celem = 0
c
      do iel=1,NRELES
        mdle = NEW_ELEM_ORDER(iel)
c
c  .....only the host process assembles data for MUMPS
        if (MUMPS_PAR%MYID == 0) then
c
c  .......evaluate the element load vector and stiffness matrix
          if (iprint_celem==1) then
            write(*,*)MYRANK,'call celem(',mdle,',2,...)'
          endif
          call celem(mdle,2,
     .               nrdofs,nrdofm,nrdofc,
     .               nodm,ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,
     .               zxload,zstiff)
          if (iprint_celem==1) then
            write(*,*)MYRANK,'nrdofs=',nrdofs
            write(*,*)MYRANK,'nrdofm=',nrdofm
            write(*,*)MYRANK,'nrdofc=',nrdofc
            write(*,*)MYRANK,'ndofmH(1:',nrnodm,')=',ndofmH
            write(*,*)MYRANK,'ndofmE(1:',nrnodm,')=',ndofmE
            write(*,*)MYRANK,'ndofmV(1:',nrnodm,')=',ndofmV
            write(*,*)MYRANK,'ndofmQ(1:',nrnodm,')=',ndofmQ
          endif
c
c  .....only the host process assembles data for MUMPS
        else
c
c  .......determine nodes of the modified element
          if (iprint_celem==1) then
            write(*,*)MYRANK,'call celem(',mdle,',1,...)'
          endif
          call celem(mdle,1,
     .               nrdofs,nrdofm,nrdofc,
     .               nodm,ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,
     .               zvoid,zvoid)
          if (iprint_celem==1) then
            write(*,*)MYRANK,'after celem(',mdle,',1,...)'
            write(*,*)MYRANK,'nrdofs=',nrdofs
            write(*,*)MYRANK,'nrdofm=',nrdofm
            write(*,*)MYRANK,'nrdofc=',nrdofc
            write(*,*)MYRANK,'ndofmH(1:',nrnodm,')=',ndofmH
            write(*,*)MYRANK,'ndofmE(1:',nrnodm,')=',ndofmE
            write(*,*)MYRANK,'ndofmV(1:',nrnodm,')=',ndofmV
            write(*,*)MYRANK,'ndofmQ(1:',nrnodm,')=',ndofmQ
          endif
c
c  .....only the host process assembles data for MUMPS
        endif
c
c-----------------------------------------------------------------------
c       STEP IV.1: FIND BIJECTION BETWEEN LOCAL AND GLOBAL DOF
c-----------------------------------------------------------------------
c
c  .....initialize local dof counter
        ilocaldof=0

c  .....initiate the element counter for the nodes
        inick=0
c
c  .....H1 dof ........
c
c  .....loop through nodes in the reversed order
        do i=nrnodm,1,-1
c
c  .......ndofmH(i) = nrdof for i-th node
          if (ndofmH(i).gt.0) then
c MP 21.4.2014 add _H1
            if (NODE_MUMPS_DATA_H1(Nodm(i))%visited.eq.0) then
              NODE_MUMPS_DATA_H1(Nodm(i))%visited=iglobaldof+1
              do idof=1,ndofmH(i)
                ilocaldof=ilocaldof+1
                iglobaldof=iglobaldof+1
                if (ilocaldof.gt.MAXDOF) then
                  write(*,*)'increase size of iglobal'
                  write(*,*) '1,ilocaldof,MAXDOF = ',ilocaldof,MAXDOF
                  stop
                endif
                iglobal(ilocaldof)=iglobaldof
              enddo
            else
              do idof=1,ndofmH(i)
                ilocaldof=ilocaldof+1
                if (ilocaldof.gt.MAXDOF) then
                  write(*,*)'increase size of iglobal'
                  write(*,*) '2,ilocaldof,MAXDOF = ',ilocaldof,MAXDOF
                  stop
                endif
                iglobal(ilocaldof) =
     .            NODE_MUMPS_DATA_H1(Nodm(i))%visited-1+idof
              enddo
            endif
          endif
        enddo
c
c
c  .....H(curl) dof ........
c
c  .....loop through nodes in the reversed order
        do i=nrnodm,1,-1
          if (ndofmE(i).gt.0) then
c MP 21.4.2014 add _Hcurl
            if (NODE_MUMPS_DATA_Hcurl(Nodm(i))%visited.eq.0) then
              NODE_MUMPS_DATA_Hcurl(Nodm(i))%visited=iglobaldof+1
              do idof=1,ndofmE(i)
                ilocaldof=ilocaldof+1
                iglobaldof=iglobaldof+1
                if(ilocaldof.gt.MAXDOF)then
                  write(*,*)'increase size of iglobal'
                  write(*,*) '3,ilocaldof,MAXDOF = ',ilocaldof,MAXDOF
                  stop
                endif
                iglobal(ilocaldof)=iglobaldof
              enddo
            else
              do idof=1,ndofmE(i)
                ilocaldof=ilocaldof+1
                if(ilocaldof.gt.MAXDOF)then
                  write(*,*)'increase size of iglobal'
                  write(*,*) '4,ilocaldof,MAXDOF = ',ilocaldof,MAXDOF
                  stop
                endif
                iglobal(ilocaldof) =
     .            NODE_MUMPS_DATA_Hcurl(Nodm(i))%visited-1+idof
              enddo
            endif
          endif
        enddo
c
c  .....H(div) dof ........
c
c  .....loop through nodes in the reversed order
        do i=nrnodm,1,-1
          if (ndofmV(i).gt.0) then
c MP 21.4.2014 add _Hdiv
            if(NODE_MUMPS_DATA_Hdiv(Nodm(i))%visited.eq.0)then
              NODE_MUMPS_DATA_Hdiv(Nodm(i))%visited=iglobaldof+1
              do idof=1,ndofmV(i)
                ilocaldof=ilocaldof+1
                iglobaldof=iglobaldof+1
                if(ilocaldof.gt.MAXDOF)then
                  write(*,*)'increase size of iglobal'
                  write(*,*) '5,ilocaldof,MAXDOF = ',ilocaldof,MAXDOF
                  stop
                endif
                iglobal(ilocaldof)=iglobaldof
              enddo
            else
              do idof=1,ndofmV(i)
                ilocaldof=ilocaldof+1
                if(ilocaldof.gt.MAXDOF)then
                  write(*,*)'increase size of iglobal'
                  write(*,*) '6,ilocaldof,MAXDOF = ',ilocaldof,MAXDOF
                  stop
                endif
                iglobal(ilocaldof) =
     .            NODE_MUMPS_DATA_Hdiv(Nodm(i))%visited-1+idof
              enddo
            endif
          endif
        enddo
c
c  .....L2 dof ........
c
c  .....middle node only
        i=nrnodm
        if (ndofmQ(i).gt.0) then
c MP 21.4.2014 add _L2
          if(NODE_MUMPS_DATA_L2(Nodm(i))%visited.eq.0)then
            NODE_MUMPS_DATA_L2(Nodm(i))%visited=iglobaldof+1
            do idof=1,ndofmQ(i)
              ilocaldof=ilocaldof+1
              iglobaldof=iglobaldof+1
              if(ilocaldof.gt.MAXDOF)then
                write(*,*)'increase size of iglobal'
                write(*,*) '7,ilocaldof,MAXDOF = ',ilocaldof,MAXDOF
                stop
              endif
              iglobal(ilocaldof)=iglobaldof
            enddo
          else
            do idof=1,ndofmQ(i)
              ilocaldof=ilocaldof+1
              if(ilocaldof.gt.MAXDOF)then
                write(*,*)'increase size of iglobal'
                write(*,*) '8,ilocaldof,MAXDOF = ',ilocaldof,MAXDOF
                stop
              endif
              iglobal(ilocaldof) =
     .          NODE_MUMPS_DATA_L2(Nodm(i))%visited-1+idof
            enddo
          endif
        endif
c
c  .....save the total number of local dof
        itotaldof=ilocaldof
c
c
c  ...only the host process assembles data for MUMPS
      if (MUMPS_PAR%MYID == 0) then
c
c  .....loop through columns of the local stiffness matrix
        do icol=1,itotaldof
c
c  .......symmetric
          if (ISYM_FLAG.eq.SYMMETRIC) then
c
c  .........loop through rows of the local stiffness matrix
            do irow=icol,itotaldof
              ind = (irow-1)*irow/2 + icol
c
c  ...........generate a new nonzero entry in the stiffness matrix
              if (abs(zstiff(ind)).gt.ZERO_MUMPS) then
                inzero=inzero+1
                if(inzero.gt.nr_nonzero)then
                  write(*,*)'increase size of %A'
                  stop
                endif
                MUMPS_PAR%A(inzero)=zstiff(ind)
                MUMPS_PAR%JCN(inzero)=iglobal(icol)
                MUMPS_PAR%IRN(inzero)=iglobal(irow)
              endif
            enddo
          else
c
c  .........loop through rows of the local stiffness matrix
            do irow=1,itotaldof
c
c  ...........generate a new nonzero entry on the stiffness matrix
c     ind=(icol-1)*itotaldof+irow
               ind=(irow-1)*itotaldof+icol
              if (abs(zstiff(ind)).gt.ZERO_MUMPS) then
                inzero=inzero+1
                if(inzero.gt.nr_nonzero)then
                  write(*,*)'increase size of %A'
                  stop
                endif
                MUMPS_PAR%A(inzero)=zstiff(ind)
                MUMPS_PAR%JCN(inzero)=iglobal(icol)
                MUMPS_PAR%IRN(inzero)=iglobal(irow)
              endif
            enddo
          endif
c
c  .......load vector (assembling)
          do irhs=1,NRHS
            if(iglobal(icol).gt.MUMPS_PAR%N)then
              write(*,*)'iglobal(icol).gt.MUMPS_PAR%N',
     .          iglobal(icol),MUMPS_PAR%N
              stop
            endif
            zrhs_for_MUMPS(iglobal(icol),irhs) =
     .        zrhs_for_MUMPS(iglobal(icol),irhs) +
     .          zxload(icol+nrdofc*(irhs-1))
          enddo
        enddo
c
c  ...only the host process assembles data for MUMPS
      endif
c
c
c  ...end of loop through elements
      enddo
c
c  ...only the host process sets %N
      if (MUMPS_PAR%MYID == 0) then
c ......here, iglobaldof is total number of utilized d.o.f.
        if(MUMPS_PAR%N.ne.iglobaldof)then
          write(*,*)'MUMPS_PAR%N=iglobaldof',MUMPS_PAR%N,iglobaldof
        endif
        MUMPS_PAR%N=iglobaldof
      endif
c
      if(itestrhs.eq.1)then
        do irhs=1,NRHS
          write(*,*)'Rhs:',irhs
          do i=1,MUMPS_PAR%N
            if(abs(zrhs_for_MUMPS(i,irhs)).gt.1.d-14)then
              write(*,*)i,zrhs_for_MUMPS(i,irhs)
            endif
          enddo
          call pause
        enddo
      endif
c
      write(*,*)'mumps_solve_seq: DATA assembled'
c
c  ...only the host process assembles data for MUMPS
      if (MUMPS_PAR%MYID == 0) then
c
c  ...write the number of nonzero entrie
      MUMPS_PAR%NZ=inzero
c      write(*,*)'MUMPS_PAR%NZ=',MUMPS_PAR%NZ

c     deallocate local memory
      deallocate(zstiff, STAT=i1)
      deallocate(zxload, STAT=i2)
      if (i1+i2 /= 0) then
        write(*,*)'mumps_solve_seq: FAILED FIRST DEALLOCATION!'
        write(*,*)'i1-4 = ',i1,i2
        stop
      endif

      write(*,*)'mumps_solve_seq: arrays deallocated'

c-----------------------------------------------------------------------
      iprint_system=0
      if(iprint_system==1.and.MYRANK.eq.0)then
        write(*,*)'MUMPS_PAR%N',MUMPS_PAR%N
        do irow=1,MUMPS_PAR%N
          write(*,*)'row',irow
          do icol=1,MUMPS_PAR%N
             do k=1,MUMPS_PAR%NZ
               if(MUMPS_PAR%IRN(k).eq.irow
     .       .and.MUMPS_PAR%JCN(k).eq.icol)then
                 write(*,*)'A(',irow,icol,')=',MUMPS_PAR%A(k)
                 exit
               endif
             enddo
          enddo
c          call pause
        enddo
      endif
c
c-----------------------------------------------------------------------
c
c     end of host process' assembly code
      endif
c
c-----------------------------------------------------------------------
c              STEP IV: SOLVE THE PROBLEM USING MUMPS
c-----------------------------------------------------------------------
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)
c
      call stop_clock(dtime,iclock_init)
      ASSEMBLY_TIME = dtime
      write(*,*)MYRANK,'DATA PREPARATION',dtime
c
c  ...run analysis/factorization phase of MUMPS
      call system_clock(t1,clock_rate,clock_max)
      call run_mumps_solve
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call start_clock(iclock_rhs)
c
c  ...loop through right hand sides
      do irhs=1,NRHS
c
c     right hand side solved on host
      if (MUMPS_PAR%MYID == 0) then
c
c  .....copy rhs
        MUMPS_PAR%RHS(1:MUMPS_PAR%N)=zrhs_for_MUMPS(1:MUMPS_PAR%N,irhs)
c
c     right hand side solved on host
      endif
c
c  .....run solve (backward elimination) phase of MUMPS
        call run_mumps_rhs
        call system_clock(t2,clock_rate,clock_max)
        write(*,4001) real(t2-t1,8)/real(clock_rate,8)
 4001   format('mumps_solve_ASSEMBLED_SEQUENTIAL = ', f13.4)
c
c     right hand side solved on host
      if (MUMPS_PAR%MYID == 0) then
c
c  .....copy solution
        zrhs_for_MUMPS(1:MUMPS_PAR%N,irhs)=MUMPS_PAR%RHS(1:MUMPS_PAR%N)
c
        if(itestrhs.eq.1)then
          write(*,*)'solution:',irhs
          do i=1,MUMPS_PAR%N
            if(abs(zrhs_for_MUMPS(i,irhs)).gt.1.d-14)then
              write(*,*)i,zrhs_for_MUMPS(i,irhs)
             endif
          enddo
          call pause
        endif
c
c     right hand side solved on host
      endif
c
      enddo
c
c  ...DEBUG.............................................................
      if (iprint>=1) then
        write(*,3000) MYRANK,MUMPS_PAR%INFO(16)
 3000   format(i4,':MUMPS TOTAL MEMORY:       ',i4,' MB')
      endif
c  ...record memory used for factorization:
      FACTOR_MEM = MUMPS_PAR%INFO(16)
c
c-----------------------------------------------------------------------
c            STEP V: STORE SOLUTION IN MODULE DATA_STRUCTURE
c-----------------------------------------------------------------------
#if PARALLEL_MODE
c  ...solution only at host processor
#if C_MODE
      call mpi_bcast(zrhs_for_MUMPS,n_out*NRHS,
     .  MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierr)
#else
      call mpi_bcast(zrhs_for_MUMPS,n_out*NRHS,
     .  MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif
#endif

c
      diff=0.d0
c
c  ...reset global node and dof counters
      nodg  = 0 ; ndofg = 0
c
c  ...loop over active elements
      do iel=1,NRELES
        mdle = NEW_ELEM_ORDER(iel)
c
c  .....get element connectivities
        call celem(mdle,1,
     .             nrdofs,nrdofm,nrdofc,
     .             nodm,ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,
     .             zvoid,zvoid)
c
c  .....initiate component counters
        mvarH = 0; mvarE = 0; mvarV = 0; mvarQ = 0
c
c  .....loop through loads
        do irhs=1,NR_RHS
c
c----------------------------------------------------------------------
c
c  .......H1 dofs......................................................
c
          if (NRHVAR.eq.0) go to 200
c
c  .......loop over element nodes in the reversed order
          do i=nrnodm,1,-1
            nod = nodm(i)
c
c MP 21.4.2014 add _H1
            nn = NODE_MUMPS_DATA_H1(Nodm(i))%visited - 1
c
c  .........compute the number of active H1 variables for the node
            call get_index(nod, index)
            if (iprint.eq.1) then
              write(*,7020) nod,index(1:NRINDEX)
 7020         format('mumps_solve: nod,index = ',i8,2x,10i2)
            endif
            nvarH=0
            do k=1,NRINDEX
              if (index(k).eq.2) nvarH=nvarH+1
            enddo
            if (nvarH.eq.0) cycle
c
c  .........loop through the nodal dof
            do j=1,ndofmH(i)/nvarH
c
c  ...........loop through the components
              ivar=mvarH(i)
              do k=1,NRINDEX
                select case(index(k))
                case(1)
                  ivar=ivar+1
                case(2)
                  ivar=ivar+1
                  nn=nn+1
c
c ................in solout, nn is the number of local dof for the node,
c                 here, we refer to global numbering
c                 (since zrhs_for_MUMPS is global)
c
                  if (icompare==1) then
ccc                 if(iel.eq.1)then
ccc                   write(*,*)'iel,nod,ivar,load,val',
ccc     .             iel,nod,ivar,irhs,zrhs_for_MUMPS(nn,irhs)
ccc                 endif
c
c  .................compare with currently stored data
                    diff = diff + abs(NODES(nod)%zdofH(ivar,j)
     .                               -zrhs_for_MUMPS(nn,irhs))
                  else
c
c  .................copy the dof into the data structure
                    if (iprint.eq.1) then
                      write(*,*) 'mumps_solve: irhs,i,ivar,nn = ',
     .                                         irhs,i,ivar,nn
                    endif
                    NODES(nod)%zdofH(ivar,j) = zrhs_for_MUMPS(nn,irhs)
                  endif
                  if (iprint.eq.1) then
                    write(*,7006) nn,zrhs_for_MUMPS(nn,irhs)
 7006               format('mumps_solve: nn, zrhs_for_MUMPS(nn) = ',
     .                      i4,2e12.5)
                    write(*,7007) nod,ivar,j,NODES(nod)%zdofH(ivar,j)
 7007               format('mumps_solve: nod,ivar,j,',
     .                     'NODES(nod)%zdofH(ivar,j) = ',
     .                      i5,i3,i3,2x,2e12.5)
                  endif
                end select
c
c  ...........end of loop through the components
              enddo
c
c  .........end of loop through the nodal dof
            enddo
c
c  .........update the number of components stored so far
            mvarH(i) = ivar
c
c  .......loop over element nodes in reversed order
          enddo
c
c-----------------------------------------------------------------------
c
c  ...H(curl) dof ......................................................
  200 if (NREVAR.eq.0) go to 300
c
c  ...loop over element nodes in reversed order
      do i=nrnodm,1,-1
        nod = nodm(i)
c
c MP 21.4.2014 add _Hcurl
        nn = NODE_MUMPS_DATA_Hcurl(Nodm(i))%visited - 1
c
c
c  .....compute the number of active H(curl) variables for the node
        call get_index(nod, index)
c
        nvarE=0
        do k=1,NRINDEX
          if (index(k).eq.4) nvarE=nvarE+1
        enddo
        if (nvarE.eq.0) cycle
c
c  .....loop through the nodal dof
        do j=1,ndofmE(i)/nvarE
c
c  .......loop through the components
          ivar=mvarE(i)
          do k=1,NRINDEX
            select case(index(k))
            case(3)
              ivar=ivar+1
            case(4)
              ivar=ivar+1
              nn=nn+1
c
c  ...........copy the dof
              NODES(nod)%zdofE(ivar,j) = zrhs_for_MUMPS(nn,irhs)
              if (iprint.eq.1) then
                write(*,7006) nn,zrhs_for_MUMPS(nn,irhs)
                write(*,7009) nod,ivar,j,NODES(nod)%zdofE(ivar,j)
 7009           format('mumps_solve: nod,ivar,j,'
     .                 'NODES(nod)%zdofE(ivar,j) = ',
     .                  i5,i3,i3,2x,2e12.5)
              endif
            end select
c  .......loop over compnents
          enddo
c  .....loop over node H(curl) dofs
        enddo
c
c  .....update the number of components stored so far
        mvarE(i) = ivar
c
c  ...loop over element nodes in reversed order
      enddo
c
c-----------------------------------------------------------------------
c
c  ...H(div) dof .................................
  300 if (NRVVAR.eq.0) go to 400
c
c  ...loop through nodes in the reversed order
      do i=nrnodm,1,-1
        nod = nodm(i)
c
c MP 21.4.2014 add _Hdiv
        nn = NODE_MUMPS_DATA_Hdiv(Nodm(i))%visited - 1
c
c
c  .....compute the number of active H(div) variables for the node
        call get_index(nod, index)
        nvarV=0
        do k=1,NRINDEX
          if (index(k).eq.6) nvarV=nvarV+1
        enddo
        if (nvarV.eq.0) cycle
c
c  .....loop through the nodal dof
        do j=1,ndofmV(i)/nvarV
c
c  .......loop through the components
          ivar=mvarV(i)
          do k=1,NRINDEX
            select case(index(k))
            case(5)
              ivar=ivar+1
            case(6)
              ivar=ivar+1
              nn=nn+1
c
c  ...........copy the dof
              NODES(nod)%zdofV(ivar,j) = zrhs_for_MUMPS(nn,irhs)
              if (iprint.eq.1) then
                write(*,7006) nn,zrhs_for_MUMPS(nn,irhs)
                write(*,7010) nod,j,ivar,NODES(nod)%zdofV(ivar,j)
 7010           format('mumps_solve: nod,j,ivar,',
     .                 'NODES(nod)%zdofV(ivar,1) = ',
     .                  i5,i3,i3,2x,2e12.5)
              endif
            end select
          enddo
        enddo
c
c  .....update the number of components stored so far
        mvarV(i) = ivar
c
      enddo
c-----------------------------------------------------------------------
c
c  ...L2 dof .................................
  400 if (NRQVAR.eq.0) go to 999
c
c  ...middle node only
      i=nrnodm
      nod = nodm(i)
c
c MP 21.4.2014 add _L2
        nn = NODE_MUMPS_DATA_L2(Nodm(i))%visited - 1
c
c
c  ...compute the number of active L2 variables for the node
      call get_index(nod, index)
      nvarQ=0
      do k=1,NRINDEX
        if (index(k).eq.8) nvarQ=nvarQ+1
      enddo
      if (nvarQ.eq.0) go to 999
c
c  ...loop through the nodal dof
      do j=1,ndofmQ(i)/nvarQ
c
c  .....loop through the components
        ivar=mvarQ(i)
        do k=1,NRINDEX
          select case(index(k))
          case(7)
            ivar=ivar+1
          case(8)
            ivar=ivar+1
            nn=nn+1
c
c  .........copy the dof
            NODES(nod)%zdofQ(ivar,j) = zrhs_for_MUMPS(nn,irhs)
            if (iprint.eq.1) then
              write(*,7006) nn,zrhs_for_MUMPS(nn,irhs)
              write(*,7011) nod,j,ivar,NODES(nod)%zdofQ(ivar,j)
 7011         format('mumps_solve: nod,j,ivar,',
     .               'NODES(nod)%zdofQ(ivar,1) = ',
     .                i5,i3,i3,2x,2e12.5)
            endif
          end select
        enddo
      enddo
c
c  ...update the number of components stored so far
      mvarQ(i) = ivar
c
c-------------------------------------------------------------------------

 999  continue

c
c  ...end of loop through rhs
      enddo
c
c  ...end loop through active elements
      enddo
c
      if(icompare==1)then
        write(*,*)'TOTAL DIFFERENCE',diff
        write(*,*)'DIFFERENCE PER ELEMENT',diff/NRELES
      endif
c-------------------------------------------------------------------------

      call stop_clock(dtime,iclock_rhs)
      write(*,*)MYRANK,'PROCESSING RHS (BS+STORE)',dtime

      call start_clock(iclock_dealloc)

c      write(*,*)MYRANK,'zrhs_for_MUMPS',zrhs_for_MUMPS(1:10,1)
c      call mpi_barrier(MPI_COMM_WORLD,ierr)
c      stop

c
c  ...only the host process assembles data for MUMPS
      if (MUMPS_PAR%MYID == 0) then
c
c  ...deallocate mumps data
      deallocate(MUMPS_PAR%JCN,    STAT=i1)
      deallocate(MUMPS_PAR%IRN,    STAT=i2)
      deallocate(MUMPS_PAR%A,      STAT=i3)
      deallocate(MUMPS_PAR%RHS,    STAT=i4)
      if (i1+i2+i3+i4 /= 0) then
        write(*,*)'mumps_solve: FAILED SECOND DEALLOCATION!'
        write(*,*)'i1-4 = ',i1,i2,i3,i4
        stop
      endif
c
c  ...only the host process assembles data for MUMPS
      endif
c
      deallocate(zrhs_for_MUMPS,    STAT=i6)
      if (i6 /= 0) then
        write(*,*)'mumps_solve: FAILED SECOND DEALLOCATION!'
        write(*,*)'i6 = ',i6
        stop
      endif
c
c  ...deallocate MUMPS memory
      call dealloc_mumps
c
c MP 21.4.2014 deallocate H1/Hdiv/Hcurl/L2
      deallocate(NODE_MUMPS_DATA_H1, STAT=i1)
      deallocate(NODE_MUMPS_DATA_Hcurl, STAT=i1)
      deallocate(NODE_MUMPS_DATA_Hdiv, STAT=i1)
      deallocate(NODE_MUMPS_DATA_L2, STAT=i1)
      if (i1 /= 0) then
        write(*,*)'mumps_solve: FAILED DEALLOCATION of NODE_MUMPS_DATA'
        write(*,*)'i6 = ',i6
        stop
      endif

c
c  ...deallocate ALL arrays used
      do i=1,NR_PHYSA
        deallocate(BLOC(i)%array)
        do j=1,1,NR_PHYSA
          deallocate(ALOC(i,j)%array)
        enddo
        deallocate(AAUX(i)%array)
      enddo
      deallocate(BLOC,AAUX,ALOC,ZBMOD,ZAMOD,NEXTRACT,IDBC,ZDOFD,MAXDOFS)
      deallocate(NEW_ELEM_ORDER)
c
      call stop_clock(dtime,iclock_dealloc)
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      write(*,*)MYRANK,'MUMPS DEALLOC',dtime

      end subroutine mumps_solve_seq
