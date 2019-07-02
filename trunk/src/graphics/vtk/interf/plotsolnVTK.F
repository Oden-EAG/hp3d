c----------------------------------------------------------------------
c
c   routine name       - plotsolnVTK
c
c----------------------------------------------------------------------
c
c   latest revision    - Nov 07
c
c   purpose            - interface to VTK graphics library
c
c   arguments          - none
c
c----------------------------------------------------------------------
c
      subroutine plotsolnVTK
c
      use data_structure3D
      use vtk_interf
c
#include "syscom.blk"
c
      common /cgeom/ iprint_geom
c
c     element orders of approximation
      dimension norder(19)
c
c     element geometry and solution dof
      dimension xnod(3,MAXbrickH)
      dimension zdofH(NRHVAR,MAXbrickH),zdofH2(MAXEQNH,MAXbrickH)
      dimension zdofE(NREVAR,MAXbrickE),zdofE2(MAXEQNE,MAXbrickE)
      dimension zdofV(NRVVAR,MAXbrickV),zdofV2(MAXEQNV,MAXbrickV)
      dimension zdofQ(NRQVAR,MAXbrickQ),zdofQ2(MAXEQNQ,MAXbrickQ)
c
c     element shape functions
      dimension shapH(MAXbrickH),gradH(3,MAXbrickH)
c
      dimension xi(3),dxdxi(3,3),dxidx(3,3)
      allocatable :: x(:,:),z(:,:)
c
c     solution and derivatives
      dimension zH(NRHVAR),  zgradH(NRHVAR,3),
     .          zE(NREVAR,3),zcurlE(NREVAR),
     .          zV(NRVVAR,3), zdivV(NRVVAR),
     .          zQ(NRQVAR)
c
      dimension nedg_orient(12),nface_orient(6)
c
c     decoded middle node case
      dimension ncase(NR_PHYSA)
c
      character(len=1) name
c
c-----------------------------------------------------------------------
c
      do i=1,NR_PHYSA
        write(*,1000) PHYSA(i),i
      enddo
      write(*,1001) NR_PHYSA+1
 1000 format(a5,  ' GRAPHICS..........',i1)
 1001 format('GLOBAL PRESSURE.........',i1)
      read(*,*) iphys
c
c-----------------------------------------------------------------------
c
      iprint      = 1
      iprint_geom = 0
c
      isol = 1
c
      if (NRCOMS>1) then
        write(*,1010) NRCOMS
 1010   format('ENTER DATA STRUCTURE COMPONENT (NRCOMS=',i2,'):')
        read(*,*) icom
      else
        icom = 1
      endif
c
c     set number of complex-valued solution components to be displayed
      if (iphys<=NR_PHYSA) then
        select case(PHYSA(iphys))
        case('elast')
          nrcomp = 4
          name   = 'u'
        case('acous')
          nrcomp = 1
          name   = 'p'
        end select
      elseif (iphys==NR_PHYSA+1) then
        nrcomp = 1
        name = 'p'
      endif
      if (iprint>=1) write(*,*)'plotsolnVTK: nrcomp =',nrcomp
c
      call init_vtk_soln(nrcomp)
c
c     count number of figures
      nrfig = 0
      mdle = 0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
c
        if (iphys<=NR_PHYSA) then
          select case(PHYSA(iphys))
          case('elast')
            if (NODES(mdle)%case/=2) cycle
          case('acous')
            if (NODES(mdle)%case/=1) cycle
          end select
        endif
        nrfig = nrfig+1
c
      enddo
c
      call set_nrfig(isol,nrfig)
c
      nres = -1
      do while(nres<0)
        write(*,*)'ENTER NONNEGATIVE RESOLUTION:'
        read(*,*) nres
      enddo
      nsub = 2**nres
c
      allocate(x(3,(nsub+1)**3),z(nrcomp,(nsub+1)**3), STAT=ic)
      if (ic/=0) then
        write(*,*)'plotsolnVTK: FAILED ALLOCATION!'
        stop
      endif
c
c     loop through elements
      mdle = 0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
c
        if (iphys<=NR_PHYSA) then
          select case(PHYSA(iphys))
          case('elast')
            if (NODES(mdle)%case/=2) cycle
          case('acous')
            if (NODES(mdle)%case/=1) cycle
          end select
        endif
c
c       get order of approximation
        call find_order(mdle, norder)
        call celndof(NODES(mdle)%type,norder,
     .               nrdofH,nrdofE,nrdofV,nrdofQ)
c
c       get node orientations
        call find_orient(mdle, nedg_orient,nface_orient)
c
c       get geometry and solution dof
        call nodcor(mdle, xnod)
        call solelm(mdle, zdofH2,zdofE2,zdofV2,zdofQ2)
c
c       copy selected component
        do k=1,nrdofH
        do i=1,NRHVAR
          zdofH(i,k) = zdofH2((icom-1)*NRHVAR+i,k)
        enddo
        enddo
c        do k=1,nrdofE
c        do i=1,NREVAR
c          zdofE(i,k) = zdofE2((icom-1)*NREVAR+i,k)
c        enddo
c        enddo
c        do k=1,nrdofV
c        do i=1,NRVVAR
c          zdofV(i,k) = zdofV2((icom-1)*NRVVAR+i,k)
c        enddo
c        enddo
c        do k=1,nrdofQ
c        do i=1,NRQVAR
c          zdofQ(i,k) = zdofQ2((icom-1)*NRQVAR+i,k)
c        enddo
c        enddo
c
        dxi = 1.d0/dble(nsub)
c
        select case(NODES(mdle)%type)
c
c-----------------------------------------------------------------------
c       BRICK
c
        case('mdlb')
c
c         loop through points
          l = 0
          do l3=0,nsub
          do l2=0,nsub
          do l1=0,nsub
            l = l+1
c
c           set reference coordinates
            xi(1) = dxi*dble(l1)
            xi(2) = dxi*dble(l2)
            xi(3) = dxi*dble(l3)
c
c           compute physical coordinates and derivatives
            call shape3H(NODES(mdle)%type,xi,norder,
     .                   nedge_orient,nface_orient,
     .                   nrdofH,shapH,gradH)
ccc            call shape3_new(NODES(mdle)%type,nedge_orient,nface_orient,
ccc     .                      xi,norder, nrdofH,shapH,gradH)
c
            x(:,l) = 0.d0
            dxdxi  = 0.d0
            do k=1,nrdofH
              x(:,l)     = x(:,l)     + xnod(:,k)*shapH(k)
              dxdxi(:,1) = dxdxi(:,1) + xnod(:,k)*gradH(1,k)
              dxdxi(:,2) = dxdxi(:,2) + xnod(:,k)*gradH(2,k)
              dxdxi(:,3) = dxdxi(:,3) + xnod(:,k)*gradH(3,k)
            enddo
            call geom(dxdxi, dxidx,rjac,igeom)
c
c           compute solution
            call eval_vtk(x(1,l),xi,dxdxi,dxidx,rjac,
     .                    nrdofH,shapH,gradH,zdofH, zH,zgradH)
c
c           select quantity to be displayed
            call disp_vtk(mdle,x(1,l),iphys,nrcomp,
     .                    zH,zgradH,zE,zcurlE,zV,zdivV,zQ,
     .                    z(1,l))
c
          enddo
          enddo
          enddo
c
          call add_hexa(isol,x,z,nsub,mdle)
c
c-----------------------------------------------------------------------
c       TETRAHEDRON
c
        case('mdln')
c
c         loop through points
          l = 0
          do l3=0,nsub
          do l2=0,nsub-l3
          do l1=0,nsub-l3-l2
            l = l+1
c
c           set reference coordinates
            xi(1) = dxi*dble(l1)
            xi(2) = dxi*dble(l2)
            xi(3) = dxi*dble(l3)
c
c           compute physical coordinates and derivatives
            call shape3H(NODES(mdle)%type,xi,norder,
     .                   nedge_orient,nface_orient,
     .                   nrdofH,shapH,gradH)

ccc            call shape3_new(NODES(mdle)%type,nedge_orient,nface_orient,
ccc     .                      xi,norder, nrdofH,shapH,gradH)
c
            x(:,l) = 0.d0
            dxdxi  = 0.d0
            do k=1,nrdofH
              x(:,l)     = x(:,l)     + xnod(:,k)*shapH(k)
              dxdxi(:,1) = dxdxi(:,1) + xnod(:,k)*gradH(1,k)
              dxdxi(:,2) = dxdxi(:,2) + xnod(:,k)*gradH(2,k)
              dxdxi(:,3) = dxdxi(:,3) + xnod(:,k)*gradH(3,k)
            enddo
            call geom(dxdxi, dxidx,rjac,igeom)
c
c           compute solution
            call eval_vtk(x(1,l),xi,dxdxi,dxidx,rjac,
     .                    nrdofH,shapH,gradH,zdofH, zH,zgradH)
c
c           select quantity to be displayed
            call disp_vtk(mdle,x(1,l),iphys,nrcomp,
     .                    zH,zgradH,zE,zcurlE,zV,zdivV,zQ,
     .                    z(1,l))
c
          enddo
          enddo
          enddo
c
          call add_tetra(isol,x,z,nsub,mdle)
c
c-----------------------------------------------------------------------
c       PRISM
c
        case('mdlp')
c
c         loop through points
          l = 0
          do l3=0,nsub
          do l2=0,nsub
          do l1=0,nsub-l2
            l = l+1
c
c           set reference coordinates
            xi(1) = dxi*dble(l1)
            xi(2) = dxi*dble(l2)
            xi(3) = dxi*dble(l3)
c
c           compute physical coordinates and derivatives
            call shape3H(NODES(mdle)%type,xi,norder,
     .                   nedge_orient,nface_orient,
     .                   nrdofH,shapH,gradH)
ccc            call shape3_new(NODES(mdle)%type,nedge_orient,nface_orient,
ccc     .                      xi,norder, nrdofH,shapH,gradH)
c
            x(:,l) = 0.d0
            dxdxi  = 0.d0
            do k=1,nrdofH
              x(:,l)     = x(:,l)     + xnod(:,k)*shapH(k)
              dxdxi(:,1) = dxdxi(:,1) + xnod(:,k)*gradH(1,k)
              dxdxi(:,2) = dxdxi(:,2) + xnod(:,k)*gradH(2,k)
              dxdxi(:,3) = dxdxi(:,3) + xnod(:,k)*gradH(3,k)
            enddo
            call geom(dxdxi, dxidx,rjac,igeom)
c
c           compute solution
            call eval_vtk(x(1,l),xi,dxdxi,dxidx,rjac,
     .                    nrdofH,shapH,gradH,zdofH, zH,zgradH)
c
c           select quantity to be displayed
            call disp_vtk(mdle,x(1,l),iphys,nrcomp,
     .                    zH,zgradH,zE,zcurlE,zV,zdivV,zQ,
     .                    z(1,l))
c
          enddo
          enddo
          enddo
c
          call add_prism(isol,x,z,nsub,mdle)
c
c-----------------------------------------------------------------------
c       PYRAMID
c
        case('mdld')
c
c         loop through points
          l = 0
          do l3=0,nsub
          do l2=0,nsub-l3
          do l1=0,nsub-l3
            l = l+1
c
c           set reference coordinates
            xi(1) = dxi*dble(l1)
            xi(2) = dxi*dble(l2)
            xi(3) = dxi*dble(l3)
c
c           compute physical coordinates and derivatives
            call shape3H(NODES(mdle)%type,xi,norder,
     .                   nedge_orient,nface_orient,
     .                   nrdofH,shapH,gradH)
ccc            call shape3_new(NODES(mdle)%type,nedge_orient,nface_orient,
ccc     .                      xi,norder, nrdofH,shapH,gradH)
c
            x(:,l) = 0.d0
            dxdxi  = 0.d0
            do k=1,nrdofH
              x(:,l)     = x(:,l)     + xnod(:,k)*shapH(k)
              dxdxi(:,1) = dxdxi(:,1) + xnod(:,k)*gradH(1,k)
              dxdxi(:,2) = dxdxi(:,2) + xnod(:,k)*gradH(2,k)
              dxdxi(:,3) = dxdxi(:,3) + xnod(:,k)*gradH(3,k)
            enddo
            call geom(dxdxi, dxidx,rjac,igeom)
c
c           compute solution
            call eval_vtk(x(1,l),xi,dxdxi,dxidx,rjac,
     .                    nrdofH,shapH,gradH,zdofH, zH,zgradH)
c
c           select quantity to be displayed
            call disp_vtk(mdle,x(1,l),iphys,nrcomp,
     .                    zH,zgradH,zE,zcurlE,zV,zdivV,zQ,
     .                    z(1,l))
c
          enddo
          enddo
          enddo
c
          call add_pyramid(isol,x,z,nsub,mdle)
c
        case default
          write(*,*)'plotsolnVTK: UNKNOWN TYPE =',NODES(mdle)%type
          stop
        end select
c
      enddo
c
      deallocate(x,z, STAT=ic)
      if (ic/=0) then
        write(*,*)'plotsolnVTK: FAILED DEALLOCATION!'
        stop
      endif
c
c-----------------------------------------------------------------------
c
      call delete_vtk_soln
c
      end subroutine plotsolnVTK
