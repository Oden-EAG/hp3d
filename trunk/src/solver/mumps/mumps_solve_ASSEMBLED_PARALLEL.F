c-----------------------------------------------------------------------
c
c   routine name       - mumps_solve_par
c
c-----------------------------------------------------------------------
c
c   latest revision    - May 07
c
c   purpose            - routine assembles the stiffness matrix,
c                        generates the interface with MUMPS, solves the
c                        resulting matrix, and stores the solution in
c                        module data_structure2D
c
c   arguments          - none
c
c-----------------------------------------------------------------------
c
      subroutine mumps_solve_par(Number_of_RHS)
c
      use data_structure3D
      use element_data
      use assembly
      use frsolmod
      use mumps_data
      use physics

      use control
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
c  ...bijection between local and global dof
      allocatable :: iglobal(:)


c  ...multiple right hand sides
      allocatable :: zrhs_for_MUMPS(:,:)
      allocatable :: zrhs_for_MUMPS_out(:,:)
c
c  ...component counters for the nodes (use in case of multiple loads)
      dimension mvarH(MAXNODM),mvarE(MAXNODM),
     .          mvarV(MAXNODM),mvarQ(MAXNODM)
c
c-----------------------------------------------------------------------
c     STEP I: INITIALIZE DATA AND DEFINE MUMPS_PAR%N
c-----------------------------------------------------------------------
c
      iprint = 0
      iprint_celem=0
      idebug = 0
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call start_clock(iclock_init)
c
c
c  ...save the number of right-hand sides (load vectors)
      NR_RHS = Number_of_RHS
      NRHS = Number_of_RHS
c
c  ...allocate memory for MUMPS
      call start_mumps
      call alloc_mumps
c
c set:
c MUMPS_PAR%N = total number of dof
c
c  ...only the host process sets %N
      if (MUMPS_PAR%MYID == 0) then
c
c  ...set values of N and NELT
c TO BE UPDATED LATER (THIS IS NOT COMPRESSED)
      MUMPS_PAR%N    = NRDOFSH + NRDOFSE + NRDOFSV + NRDOFSQ
c
c  ...only the host process sets %N
      endif
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
c-----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c     STEP II: COMPUTE AN UPPER BOUND OF THE NUMBER OF NONZERO
c              ENTRIES AND THE MAXIMUM NUMBER OF DOF PER ELEMENT
c----------------------------------------------------------------------
c
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

c        call findinitelem(mdle,iel_fath)
c        if(ELEMS_PARALLEL_DATA(iel_fath)%processor_owner.ne.
c     .    MYRANK)cycle
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
c-----------------------------------------------------------------------
c     STEP III: ALLOCATE MEMORY
c-----------------------------------------------------------------------
c
c  ...allocate memory
#if PARALLEL_MODE
      allocate(MUMPS_PAR%IRN_loc(nr_nonzero),        STAT=i1)
      allocate(MUMPS_PAR%JCN_loc(nr_nonzero),        STAT=i2)
      allocate(MUMPS_PAR%A_loc(nr_nonzero),          STAT=i3)
#else
      allocate(MUMPS_PAR%IRN_loc(nr_nonzero),        STAT=i1)
      allocate(MUMPS_PAR%JCN_loc(nr_nonzero),        STAT=i2)
      allocate(MUMPS_PAR%A_loc(nr_nonzero),          STAT=i3)
#endif
c  ...only host processor allocates RHS
      if (MUMPS_PAR%MYID == 0) then
        allocate(MUMPS_PAR%RHS(MUMPS_PAR%N),       STAT=i4)
      else
        i4 = 0
      endif
      if (ISYM_FLAG.eq.SYMMETRIC) then
        allocate(zstiff(MAXDOFM*(MAXDOFM+1)/2),      STAT=i5)
      else
        allocate(zstiff(MAXDOFM*MAXDOFM),            STAT=i5)
      endif
      allocate(zxload(MAXDOFM*NRHS),               STAT=i6)
      if (i1+i2+i3+i4+i5+i6 /= 0) then
        write(*,*)'mumps_solve_assembled_host: ',
     .    'FAILED SECOND ALLOCATION!'
        write(*,*)'i1-8 = ',i1,i2,i3,i4,i5,i6
        write(*,*)'MAXDOFM,MAXDOFC=',MAXDOFM,MAXDOFC
        stop
      endif
c
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
#if PARALLEL_MODE
      nout = 0;
      if(MYRANK.eq.0)n_out = MUMPS_PAR%N
      call mpi_bcast(n_out,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      if(idebug.eq.1)then
        write(*,*)MYRANK,'allocate(zrhs_for_MUMPS(',n_out,NRHS,')'
      endif
      allocate(zrhs_for_MUMPS(n_out,NRHS), STAT=i9)
      allocate(zrhs_for_MUMPS_out(n_out,NRHS), STAT=i9)
#else
      allocate(zrhs_for_MUMPS(MUMPS_PAR%N,NRHS), STAT=i9)
      allocate(zrhs_for_MUMPS_out(MUMPS_PAR%N,NRHS), STAT=i9)
#endif
      if (i9 /= 0) then
        write(*,*)'mumps_solve_assembled_host: ',
     .    'FAILED SECOND ALLOCATION!'
        write(*,*)'i9',i9
        stop
      endif
c
c  ...initialize the right hand side
      zrhs_for_MUMPS=ZERO
      zrhs_for_MUMPS_out=ZERO
c
c  ...save the maximum number of dof
      MAXDOF = MAXDOFC
c
      allocate(iglobal(MAXDOF),STAT=i1)
      if(i1/=0)then
        write(*,*)'mumps_solve_AS:Error allocation iglobal'
        stop
      endif
c
c-----------------------------------------------------------------------
c
      if(.not.allocated(NEW_ELEM_ORDER))then
c
c  ...reorder the elements in order to minimize the bandwidth...
c     OLD/NEW ELEM_ORDER are defined in frsolmod
      direction(1:3) = 0.d0; direction(1) = 1.d0;
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
c
      endif
c
c
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c     STEP IV: ASSEMBLE THE STIFFNESS MATRIX AND THE LOAD VECTOR
c-----------------------------------------------------------------------
c
c  ...allocate data structure for storing inverse maps
      allocate(NODE_MUMPS_DATA(NRNODS),STAT=i1)
      if(i1/=0)then
        write(*,*)'NODE_MUMPS_DATA allocation error',i1
        stop
      endif

c  ...the map must be found on all processors
c  ...initialize visited flags
      do inode=1,NRNODS
        NODE_MUMPS_DATA(inode)%visited=0
      enddo

c  ...initialize global dof counter
      iglobaldof=0 ; inzero = 0
c
c  ...loop through elements
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
c
c  ...only processor owner performs integration
      call findinitelem(mdle,iel_fath)
      if(ELEMS_PARALLEL_DATA(iel_fath)%processor_owner.eq.MYRANK)then
c
c  .....evaluate the element load vector and stiffness matrix
        if(iprint_celem==1)then
          write(*,*)MYRANK,'call celem(',mdle,',2,...)'
        endif
        call celem(mdle,2,
     .             nrdofs,nrdofm,nrdofc,
     .             nodm,ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,
     .             zxload,zstiff)
        if(iprint_celem==1)then
          write(*,*)'nrdofs=',nrdofs
          write(*,*)'nrdofm=',nrdofm
          write(*,*)'nrdofc=',nrdofc
          write(*,*)'ndofmH(1:',nrnodm,')=',ndofmH
          write(*,*)'ndofmE(1:',nrnodm,')=',ndofmE
          write(*,*)'ndofmV(1:',nrnodm,')=',ndofmV
          write(*,*)'ndofmQ(1:',nrnodm,')=',ndofmQ
        endif
c
c  ...only processor owner performs integration
      else
c
c  .....determine nodes of the modified elemen
        if(iprint_celem.eq.1)then
          write(*,*)MYRANK,'call celem(',mdle,',1,...)'
        endif
        call celem(mdle,1,
     .             nrdofs,nrdofm,nrdofc,
     .             nodm,ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,
     .             zvoid,zvoid)
c
c  ...only processor owner performs integration
      endif
c
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
c  .......ndofmH(i) = nrdof for i-th node
          if (ndofmH(i).gt.0) then
            if(NODE_MUMPS_DATA(Nodm(i))%visited.eq.0)then
              NODE_MUMPS_DATA(Nodm(i))%visited=iglobaldof+1
              do idof=1,ndofmH(i)
                ilocaldof=ilocaldof+1
                iglobaldof=iglobaldof+1
                if(ilocaldof.gt.MAXDOF)then
                  write(*,*)'increase size of iglobal'
                  stop
                endif
                iglobal(ilocaldof)=iglobaldof
              enddo
            else
              do idof=1,ndofmH(i)
                ilocaldof=ilocaldof+1
                if(ilocaldof.gt.MAXDOF)then
                  write(*,*)'increase size of iglobal'
                  stop
                endif
                iglobal(ilocaldof) =
     .            NODE_MUMPS_DATA(Nodm(i))%visited-1+idof
              enddo
            endif
          endif
        enddo
c
c  .....H(curl) dof ........
c
c  .....loop through nodes in the reversed order
        do i=nrnodm,1,-1
          if (ndofmE(i).gt.0) then
            if(NODE_MUMPS_DATA(Nodm(i))%visited.eq.0)then
              NODE_MUMPS_DATA(Nodm(i))%visited=iglobaldof+1
              do idof=1,ndofmE(i)
                ilocaldof=ilocaldof+1
                iglobaldof=iglobaldof+1
                if(ilocaldof.gt.MAXDOF)then
                  write(*,*)'increase size of iglobal'
                  stop
                endif
                iglobal(ilocaldof)=iglobaldof
              enddo
            else
              do idof=1,ndofmE(i)
                ilocaldof=ilocaldof+1
                if(ilocaldof.gt.MAXDOF)then
                  write(*,*)'increase size of iglobal'
                  stop
                endif
                iglobal(ilocaldof) =
     .            NODE_MUMPS_DATA(Nodm(i))%visited-1+idof
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
            if(NODE_MUMPS_DATA(Nodm(i))%visited.eq.0)then
              NODE_MUMPS_DATA(Nodm(i))%visited=iglobaldof+1
              do idof=1,ndofmV(i)
                ilocaldof=ilocaldof+1
                iglobaldof=iglobaldof+1
                if(ilocaldof.gt.MAXDOF)then
                  write(*,*)'increase size of iglobal'
                  stop
                endif
                iglobal(ilocaldof)=iglobaldof
              enddo
            else
              do idof=1,ndofmV(i)
                ilocaldof=ilocaldof+1
                if(ilocaldof.gt.MAXDOF)then
                  write(*,*)'increase size of iglobal'
                  stop
                endif
                iglobal(ilocaldof) =
     .            NODE_MUMPS_DATA(Nodm(i))%visited-1+idof
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
          if(NODE_MUMPS_DATA(Nodm(i))%visited.eq.0)then
            NODE_MUMPS_DATA(Nodm(i))%visited=iglobaldof+1
            do idof=1,ndofmQ(i)
              ilocaldof=ilocaldof+1
              iglobaldof=iglobaldof+1
              if(ilocaldof.gt.MAXDOF)then
                write(*,*)'increase size of iglobal'
                stop
              endif
              iglobal(ilocaldof)=iglobaldof
            enddo
          else
            do idof=1,ndofmQ(i)
              ilocaldof=ilocaldof+1
              if(ilocaldof.gt.MAXDOF)then
                write(*,*)'increase size of iglobal'
                stop
              endif
              iglobal(ilocaldof) =
     .          NODE_MUMPS_DATA(Nodm(i))%visited-1+idof
            enddo
          endif
        endif
c
c  .....save the total number of local dof
        itotaldof=ilocaldof
c
c  ...only the processor owner assembled data for mumps
      call findinitelem(mdle,iel_fath)
      if(ELEMS_PARALLEL_DATA(iel_fath)%processor_owner.eq.MYRANK)then
c
c  .....loop through columns of the local stiffness matrix
        do icol=1,itotaldof
c
c  .......symmetric
          if (ISYM_FLAG.eq.SYMMETRIC) then
c
c  .........loop through rows of the local stiffness matrix
            do irow=icol,itotaldof
              int_index = (irow-1)*irow/2 + icol
c
c  ...........generate a new nonzero entry on the stiffness matrix
c              if (abs(zstiff(index)).gt.ZERO_MUMPS) then
                inzero=inzero+1
                if(inzero.gt.nr_nonzero)then
                  write(*,*)'inzero.gt.nr_nonzero',inzero,nr_nonzero
                  stop
                endif
                MUMPS_PAR%A_loc(inzero)=zstiff(int_index)
                MUMPS_PAR%JCN_loc(inzero)=iglobal(icol)
                MUMPS_PAR%IRN_loc(inzero)=iglobal(irow)
c              endif
            enddo
          else
c
c  .........loop through rows of the local stiffness matrix
            do irow=1,itotaldof
c
c  ...........generate a new nonzero entry on the stiffness matrix
              int_index=(icol-1)*itotaldof+irow
c              if (abs(zstiff(index)).gt.ZERO_MUMPS) then
                inzero=inzero+1
                if(inzero.gt.nr_nonzero)then
                  write(*,*)'inzero.gt.nr_nonzero',inzero,nr_nonzero
                  stop
                endif
                MUMPS_PAR%A_loc(inzero)=zstiff(int_index)
                MUMPS_PAR%JCN_loc(inzero)=iglobal(icol)
                MUMPS_PAR%IRN_loc(inzero)=iglobal(irow)
c              endif
            enddo
          endif
c
c  .......load vector (assembling)
          do irhs=1,NRHS
            if(iglobal(icol).gt.n_out)then
              write(*,*)'iglobal(',icol,').gt.n_out',iglobal(icol),n_out
              stop
            endif
            zrhs_for_MUMPS(iglobal(icol),irhs) =
     .        zrhs_for_MUMPS(iglobal(icol),irhs) +
     .          zxload(icol+(irhs-1)*nrdofc)
         enddo
        enddo
c
c  ...only the processor owner does assembling
      endif
c
c  ...end of loop through elements
      enddo
c
c  ...only the host process sets %N
      if (MUMPS_PAR%MYID == 0) then
c ......here, iglobaldof is total number of utilized d.o.f.
        MUMPS_PAR%N=iglobaldof
      endif
c
c  ...write the number of nonzero entries
      MUMPS_PAR%NZ_loc=inzero
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)
c
c     deallocate local memory
      deallocate(zstiff, STAT=i1)
      deallocate(zxload, STAT=i2)
      if (i1+i2 /= 0) then
        write(*,*)'mumps_solve_assembled_host: ',
     .    'FAILED FIRST DEALLOCATION!'
        write(*,*)'i1-4 = ',i1,i2
        stop
      endif

#if PARALLEL_MODE
c  ...contributions to RHSs should be collected (only on processor 0)
#if C_MODE
      call mpi_allreduce(zrhs_for_MUMPS,zrhs_for_MUMPS_out,
     .  n_out*NRHS,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
#else
      call mpi_allreduce(zrhs_for_MUMPS,zrhs_for_MUMPS_out,
     .  n_out*NRHS,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif
      zrhs_for_MUMPS(1:n_out,1:NRHS)=zrhs_for_MUMPS_out(1:n_out,1:NRHS)
#endif
c
c-----------------------------------------------------------------------
c              STEP IV.3: SOLVE THE PROBLEM USING MUMPS
c-----------------------------------------------------------------------
c
      call stop_clock(dtime,iclock_init)
      write(*,*)MYRANK,'DATA PREPARATION',dtime

c  ...run analysis/factorization phase of MUMPS
      call run_mumps_solve
c
c  ...DEBUG.............................................................
c      if (iprint>=1) then
        write(*,3000) MYRANK,MUMPS_PAR%INFO(16)
 3000   format(i4,':run_mumps_solve MUMPS TOTAL MEMORY:',i4,' MB')
c      endif
c
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call start_clock(iclock_rhs)

c  ...loop through right hand sides
      do irhs=1,NRHS
c
c     right hand side stored on host
      if (MUMPS_PAR%MYID == 0) then
c
c  .....copy rhs
        MUMPS_PAR%RHS(1:MUMPS_PAR%N)=zrhs_for_MUMPS(1:MUMPS_PAR%N,irhs)
c
c     right hand side stored on host
      endif
c
c  .....run solve (backward elimination) phase of MUMPS
        call run_mumps_rhs
c
c     right hand side solved on host
      if (MUMPS_PAR%MYID == 0) then
c
c  .....copy solution
        zrhs_for_MUMPS(1:MUMPS_PAR%N,irhs)=MUMPS_PAR%RHS(1:MUMPS_PAR%N)
c
c     right hand side solved on host
      endif
c
      enddo
c  ...DEBUG.............................................................
c      if (iprint>=1) then
        write(*,3017) MYRANK,MUMPS_PAR%INFO(16)
 3017   format(i4,':run_mumps_rhs MUMPS TOTAL MEMORY:       ',i4,' MB')
c      endif
c
c-----------------------------------------------------------------------
c            STEP V: STORE SOLUTION IN MODULE DATA_STRUCTURE2
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
c  ...reset global node and dof counters
      nodg  = 0 ; ndofg = 0
c
      do iel=1,NRELES
        mdle = NEW_ELEM_ORDER(iel)
c
c  .....get element connectivities
        call celem(mdle,1,
     .             nrdofs,nrdofm,nrdofc,
     .             nodm,ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,
     .             zvoid,zvoid)

c
c  ...initiate component counters
      mvarH = 0; mvarE = 0; mvarV = 0; mvarQ = 0
c
      do irhs=1,NRHS
c
c  ...H1 dof .................................
      if (NRHVAR.eq.0) go to 200
c
c  ...loop through nodes in the reversed order
      do i=nrnodm,1,-1
        nod = nodm(i)
c->
        nn = NODE_MUMPS_DATA(Nodm(i))%visited - 1
c<-
c
c  .....compute the number of active H1 variables for the node
        call get_index(nod, index)
        if (iprint.eq.1) then
          write(*,7020) nod,index(1:NRINDEX)
 7020     format('solout: nod,index = ',i8,2x,10i2)
        endif
        nvarH=0
        do k=1,NRINDEX
          if (index(k).eq.2) nvarH=nvarH+1
        enddo
        if (nvarH.eq.0) cycle
c
c  .....loop through the nodal dof
        do j=1,ndofmH(i)/nvarH
c
c  .......loop through the components
          ivar=mvarH(i)
          do k=1,NRINDEX
            select case(index(k))
            case(1)
              ivar=ivar+1
            case(2)
              ivar=ivar+1
              nn=nn+1
c
c  ...........copy the dof
              NODES(nod)%zdofH(ivar,j) = zrhs_for_MUMPS(nn,irhs)
              if (iprint.eq.1) then
                write(*,7006) nn,zrhs_for_MUMPS(nn,irhs)
 7006           format('solout: nn, zrhs_for_MUMPS(nn) = ',i4,2e12.5)
                write(*,7007) nod,j,ivar,NODES(nod)%zdofH(ivar,j)
 7007           format('solout: nod,j,ivar,NODES(nod)%zdofH(ivar,1) = ',
     .                          i5,i3,i3,2x,2e12.5)
              endif
            end select
          enddo
        enddo
c
c  .....update the number of components stored so far
        mvarH(i) = ivar
c
      enddo
c
c-----------------------------------------------------------------------
c
c  ...H(curl) dof .................................
  200 if (NREVAR.eq.0) go to 300
c
c  ...loop through nodes in the reversed order
      do ii=nrnodm,1,-1
        nod = nodm(i)
c->
        nn = NODE_MUMPS_DATA(Nodm(i))%visited - 1
c<-
c
c  .....compute the number of active H(curl) variables for the node
        call get_index(nod, index)
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
                write(*,7009) nod,j,ivar,NODES(nod)%zdofE(ivar,j)
 7009           format('solout: nod,j,ivar,NODES(nod)%zdofE(ivar,1) = ',
     .                          i5,i3,i3,2x,2e12.5)
              endif
            end select
          enddo
        enddo
c
c  .....update the number of components stored so far
        mvarE(i) = ivar
c
      enddo
c-----------------------------------------------------------------------
c
c  ...H(div) dof .................................
  300 if (NRVVAR.eq.0) go to 400
c
c  ...loop through nodes in the reversed order
      do ii=nrnodm,1,-1
        nod = nodm(i)
c->
        nn = NODE_MUMPS_DATA(Nodm(i))%visited - 1
c<-
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
 7010           format('solout: nod,j,ivar,NODES(nod)%zdofV(ivar,1) = ',
     .                          i5,i3,i3,2x,2e12.5)
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
c->
        nn = NODE_MUMPS_DATA(Nodm(i))%visited - 1
c<-
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
 7011         format('solout: nod,j,ivar,NODES(nod)%zdofQ(ivar,1) = ',
     .                      i5,i3,i3,2x,2e12.5)
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
c  ...end loop through active elements
      enddo
c
c  ...end of loop through RHS
      enddo
c
c-------------------------------------------------------------------------
c
      call stop_clock(dtime,iclock_rhs)
      write(*,*)MYRANK,'PROCESSING RHS (BS+STORE)',dtime

      call start_clock(iclock_dealloc)
c  ...deallocate mumps data
c----->
c#if PARALLEL_MODE
c      i1=0; i2=0; i3=0
c#else
      deallocate(MUMPS_PAR%JCN_loc,    STAT=i1)
      deallocate(MUMPS_PAR%IRN_loc,    STAT=i2)
      deallocate(MUMPS_PAR%A_loc,      STAT=i3)
c#endif
c<-----
c  ...only the host process assembled RHS
      if (MUMPS_PAR%MYID == 0) then
        deallocate(MUMPS_PAR%RHS,    STAT=i4)
      else
        i4=0
      endif
      if (i1+i2+i3+i4 /= 0) then
        write(*,*)'mumps_solve_assembled_host: ',
     .    'FAILED SECOND DEALLOCATION!'
        write(*,*)'i1-4 = ',i1,i2,i3,i4
        stop
      endif
c
      deallocate(zrhs_for_MUMPS,    STAT=i6)
      deallocate(zrhs_for_MUMPS_out,STAT=i7)
      if (i6+i7 /= 0) then
        write(*,*)'mumps_solve_assembled_host: ',
     .    'FAILED SECOND DEALLOCATION!'
        write(*,*)'i6 = ',i6,i7
        stop
      endif
c
#if PARALLEL_MODE
#else
c  ...deallocate MUMPS memory
      call dealloc_mumps
#endif
c
      deallocate(NODE_MUMPS_DATA, STAT=i1)
      if (i1 /= 0) then
        write(*,*)'mumps_solve: FAILED DEALLOCATION of NODE_MUMPS_DATA'
        write(*,*)'i6 = ',i6
        stop
      endif

c
c  ...deallocate ALL arrays used by frontal solver
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
c
      end subroutine mumps_solve_par
