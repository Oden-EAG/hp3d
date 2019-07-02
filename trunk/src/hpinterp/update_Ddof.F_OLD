c-----------------------------------------------------------------------
c
c   routine name       - update_Ddof
c
c-----------------------------------------------------------------------
c
c   latest revision    - Jan 11
c
c   purpose            - routine updates values of solution degrees
c                        of freedom for Dirichlet nodes
c
c-----------------------------------------------------------------------
c
      subroutine update_Ddof
c
      use data_structure3D
      use environment, only : QUIET_MODE
#include "syscom.blk"
c
c  ...orientation of element nodes
      dimension nedge_orient(12),nface_orient(6)
c
c  ...element nodes
      dimension nodesl(27), norientl(27)
c
c  ...order of approximation for an element
      dimension norder(19)
c
c  ...reference coordinates for an element
      dimension xsub(3,8)
c
c  ...solution dof for an element
      dimension zdofH(MAXEQNH,MAXbrickH)
      dimension zdofE(MAXEQNE,MAXbrickE)
      dimension zdofV(MAXEQNV,MAXbrickV)
      dimension zdofQ(MAXEQNQ,MAXbrickQ)
c
c  ...auxiliary variables for timing
      real*8 :: time1, time2, OMP_get_wtime
c
c-----------------------------------------------------------------------
c
      time1 = OMP_get_wtime()
      iprint=0
c
c  ...lower the GMP interface flag for all nodes
      do nod=1,NRNODS
        NODES(nod)%geom_interf=0
      enddo
c
c-----------------------------------------------------------------------
c
c  Step 1: Update   V E R T   dof for Dirichlet nodes
c
      if (iprint.ge.1) write(*,*) 'update_Ddof: VERTEX NODES...'
c
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        call refel(mdle, iflag,no,xsub)
        call elem_nodes(mdle, nodesl, norientl)
        do iv=1,nvert(NODES(mdle)%type)
          nod = nodesl(iv)
          if (nod.lt.0) cycle
          if ((NODES(nod)%geom_interf.eq.0).and.is_dirichlet(nod)) then
            if (iprint.eq.2) then
              write(*,7002) mdle,iflag,no,iv,nod
 7002         format(' update_Ddof: mdle,iflag,no,iv,nod = ',
     .                i8,2x,i1,2x,i8,2x,i2,2x,i8)
            endif
c  .........cycle if H1 dofs are not supported by the node
            if (.not.associated(NODES(nod)%zdofH)) cycle

            if (is_dirichlet_homogeneous(nod)) then
              NODES(nod)%zdofH = ZERO
            else
              call dhpvert(mdle,iflag,no,xsub(1:3,iv),NODES(nod)%case,
     .          NODES(nod)%zdofH)
            endif
            if (iprint.eq.2) then
              do ivar=1,NREQNH(NODES(nod)%case)
                write(*,7003) ivar,NODES(nod)%zdofH(ivar,1)
 7003           format(' update_Ddof: ivar,zdofH(ivar) = ',i2,2x,2e12.5)
              enddo
              call pause
            endif
            NODES(nod)%geom_interf=1
          endif
        enddo
      enddo
c
c-----------------------------------------------------------------------
c
c  Step 2: Update   E D G E   dof for Dirichlet nodes
c
      if (iprint.ge.1) write(*,*) 'update_Ddof: MID-EDGE NODES...'
c
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        iprint_parameter = mdle
        call refel(mdle, iflag,no,xsub)
        call elem_nodes(mdle, nodesl, norientl)
        call find_orient(mdle, nedge_orient,nface_orient)
        call find_order(mdle, norder)
c
c  .....compute solution dofs (need for H1 update)
        call solelm(mdle, zdofH,zdofE,zdofV,zdofQ)
        do ie=1,nedge(NODES(mdle)%type)
          in = nvert(NODES(mdle)%type)+ie
          nod = nodesl(in)
c
          if (nod.lt.0) cycle
          if ((NODES(nod)%geom_interf.eq.0).and.is_dirichlet(nod)) then

c
c  .........update H1 Dirichlet dofs
            if (associated(NODES(nod)%zdofH)) then
              if (is_dirichlet_homogeneous(nod)) then
                NODES(nod)%zdofH = ZERO
              else
                call dhpedgeH(mdle,iflag,no,xsub,
     .                        NODES(mdle)%type,NODES(nod)%case,
     .                        nedge_orient,nface_orient,norder,ie,
     .                        zdofH, NODES(nod)%zdofH)
              endif
            endif
c
c  .........update H(curl) Dirichlet dofs
            if (associated(NODES(nod)%zdofE)) then
              if (is_dirichlet_homogeneous(nod)) then
                NODES(nod)%zdofE = ZERO
              else
                call dhpedgeE(mdle,iflag,no,xsub,
     .                        NODES(mdle)%type,NODES(nod)%case,
     .                        nedge_orient,nface_orient,norder,ie,
     .                        NODES(nod)%zdofE)
              endif
            endif
c
            if (iprint.eq.3) then
              call find_ndof(nod, ndofH,ndofE,ndofV,ndofQ)
              write(*,7011) nod,norder(ie)
 7011         format(' update_Ddof: nod,nord = ',i8,i3,
     .               ' NODES(nod)%zdofH = ')
              nvarH = NREQNH(NODES(nod)%case)*NRCOMS
              do k=1,ndofH
                write(*,7012) k,NODES(nod)%zdofH(1:nvarH,k)
#if C_MODE
 7012           format(i2,2x,30(/,5(2e12.5,2x)))
#else
 7012           format(i2,2x,30(/,10(e12.5,2x)))
#endif
              enddo
              write(*,7013)
 7013         format(' update_Ddof: NODES(nod)%zdofE = ')
              nvarE = NREQNE(NODES(nod)%case)*NRCOMS
              do k=1,ndofE
                write(*,7012) k,NODES(nod)%zdofE(1:nvarE,k)
              enddo
              call pause
            endif
            NODES(nod)%geom_interf=1
          endif
        enddo
      enddo
c
c-----------------------------------------------------------------------
c
c  Step 3: Update   F A C E   dof for Dirichlet nodes
c
      if (iprint.ge.1) write(*,*) 'update_Ddof: MID-FACE NODES...'
c
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        call refel(mdle, iflag,no,xsub)
        call elem_nodes(mdle, nodesl, norientl)
        call find_orient(mdle, nedge_orient,nface_orient)
        call find_order(mdle, norder)
c
c  .....compute solution dofs (needed for H1 and Hcurl update)
        call solelm(mdle, zdofH,zdofE,zdofV,zdofQ)
c
c  .....loop over faces
        do if=1,nface(NODES(mdle)%type)
c
c  .......get local node number
          in = nvert(NODES(mdle)%type)+nedge(NODES(mdle)%type)+if
c
c  .......get global node number
          nod = nodesl(in)
          if (nod.lt.0) cycle
          if (NODES(nod)%geom_interf.eq.0.and.is_dirichlet(nod)) then
c
c           -- H1 --
            if (associated(NODES(nod)%zdofH)) then
              if (is_dirichlet_homogeneous(nod)) then
                NODES(nod)%zdofH = ZERO
              else
                 call dhpfaceH(mdle,iflag,no,xsub,
     .                         NODES(mdle)%type,NODES(nod)%case,
     .                         nedge_orient,nface_orient,norder,if,
     .                         zdofH, NODES(nod)%zdofH)
               endif
             endif
c
c           -- H(curl) --
            if (associated(NODES(nod)%zdofE)) then
              if (is_dirichlet_homogeneous(nod)) then
                NODES(nod)%zdofE = ZERO
              else
                call dhpfaceE(mdle,iflag,no,xsub,
     .                        NODES(mdle)%type,NODES(nod)%case,
     .                        nedge_orient,nface_orient,norder,if,
     .                        zdofE, NODES(nod)%zdofE)
               endif
             endif
c
c           -- H(div) --
            if (associated(NODES(nod)%zdofV)) then
              if (is_dirichlet_homogeneous(nod)) then
                NODES(nod)%zdofV = ZERO
              else
                call dhpfaceV(mdle,iflag,no,xsub,
     .                        NODES(mdle)%type,NODES(nod)%case,
     .                        nedge_orient,nface_orient,norder,if,
cccc
cccc LD, Mar 17, 2017: have eliminated zdofV form the list of arguments
cccc     .                        zdofV, NODES(nod)%zdofV)
     .                        NODES(nod)%zdofV)
               endif
             endif
c
            if (iprint.eq.4) then
              call find_ndof(nod, ndofH,ndofE,ndofV,ndofQ)
              write(*,7011) nod,norder(ie)
              nvarH = NREQNH(NODES(nod)%case)
              do k=1,ndofH
                write(*,7012) k,NODES(nod)%zdofH(1:nvarH,k)
              enddo
            endif
            NODES(nod)%geom_interf=1
          endif
c  .....end of loop over faces
        enddo
      enddo
c
      if (.not. QUIET_MODE) then
         time2 = OMP_get_wtime()
         write(*,8010) time2-time1
 8010    format(' update_Ddof: ',f12.5,'  seconds',/)
      endif
c
      if (iprint.ge.1) write(*,*) 'update_Ddof: DONE...'
c
      end subroutine
