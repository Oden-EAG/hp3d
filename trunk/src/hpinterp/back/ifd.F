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
c
c-----------------------------------------------------------------------
c
      iprint=1
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
          if ((NODES(nod)%geom_interf.eq.0).and.(ifD(nod).eq.1)) then
            if (iprint.eq.2) then
              write(*,7002) mdle,iflag,no,iv,nod
 7002         format(' update_Ddof: mdle,iflag,no,iv,nod = ',
     .                i8,2x,i1,2x,i8,2x,i2,2x,i8)
            endif
c  .........cycle if H1 dofs are not supported by the node
            if (.not.associated(NODES(nod)%zdofH)) cycle
            call dhpvert(iflag,no,xsub(1:3,iv),NODES(nod)%case,
     .                   NODES(nod)%zdofH)
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
c  .....compute solution dofs (need for H1 update)
        call solelm(mdle, zdofH,zdofE,zdofV,zdofQ)
        do ie=1,nedge(NODES(mdle)%type)
          in = nvert(NODES(mdle)%type)+ie
          nod = nodesl(in)
c
          if (nod.lt.0) cycle
          if ((NODES(nod)%geom_interf.eq.0).and.(ifD(nod).eq.1)) then
cc            write(*,*) 'update_Ddof : MAXEQNV ', MAXEQNV

c
c  .........update H^1 Dirichlet dofs
            if (associated(NODES(nod)%zdofH)) then
              call dhpedgeH(iflag,no,xsub,
     .                      NODES(mdle)%type,NODES(nod)%case,
     .                      nedge_orient,nface_orient,norder,ie,
     .                      zdofH, NODES(nod)%zdofH)
            endif
c
c  .........update H(curl) Dirichlet dofs
            if (associated(NODES(nod)%zdofE)) then
              call dhpedgeE(iflag,no,xsub,
     .                      NODES(mdle)%type,NODES(nod)%case,
     .                      nedge_orient,nface_orient,norder,ie,
     .                      NODES(nod)%zdofE)
            endif
c
c  .........update H(div) Dirichlet dofs
            if (associated(NODES(nod)%zdofV)) then
              call dhpedgeV(iflag,no,xsub,
     .                      NODES(mdle)%type,NODES(nod)%case,
     .                      nedge_orient,nface_orient,norder,ie,
     .                      NODES(nod)%zdofV)
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
 7012           format(i2,2x,30(/,4(2e12.5,2x)))
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
c  .....compute solution dofs (needed for H1 and Hcurl update)
        call solelm(mdle, zdofH,zdofE,zdofV,zdofQ)
c  .....loop over faces
        do if=1,nface(NODES(mdle)%type)
c  .......get local node number
          in = nvert(NODES(mdle)%type)+nedge(NODES(mdle)%type)+if
c  .......get global node number
          nod = nodesl(in)
          if (nod.lt.0) cycle
          if (NODES(nod)%geom_interf.eq.0.and.(ifD(nod).eq.1))then
c
c  .........update H^1 Dirichlet dofs
            if (associated(NODES(nod)%zdofH)) then
              call dhpfaceH(iflag,no,xsub,
     .                      NODES(mdle)%type,NODES(nod)%case,
     .                      nedge_orient,nface_orient,norder,if,
     .                      zdofH, NODES(nod)%zdofH)
            endif
c
c  .........update H(curl) Dirichlet dofs
            if (associated(NODES(nod)%zdofE)) then
              call dhpfaceE(iflag,no,xsub,
     .                      NODES(mdle)%type,NODES(nod)%case,
     .                      nedge_orient,nface_orient,norder,if,
     .                      zdofE, NODES(nod)%zdofE)
            endif
c
c  .........update H(div) Dirichlet dofs
            if (associated(NODES(nod)%zdofV)) then
              call dhpfaceV(iflag,no,xsub,
     .                      NODES(mdle)%type,NODES(nod)%case,
     .                      nedge_orient,nface_orient,norder,if,
     .                      zdofV, NODES(nod)%zdofV)
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
c
      end
c
c-----------------------------------------------------------------------
c
c   routine name       - IfD
c
c-----------------------------------------------------------------------
c
c
c   latest revision    - Jan 07
c
c   purpose            - routine checks if the node has a Dirichlet BC
c                        flag
c
c   arguments
c     in:
c                 Nod  - node number
c     out:
c                 IfD  = 'true' if a Dirichlet node
c                      = 'false' otherwise
c
c-----------------------------------------------------------------------
c
      integer function IfD(Nod)
c
      use data_structure3D
c
c  ...node BC flags
      integer ibc(NR_PHYSA)
c
c  ...decode the node BC flag
      call decod(NODES(Nod)%bcond,10,NR_PHYSA, ibc)
      IfD = 0
      do iphys=1,NR_PHYSA
        select case(ibc(iphys))
        case(1,6,7,8)
          IfD = 1
        end select
      enddo
c
c
      end
