c----------------------------------------------------------------------
c
c   routine name       - prism
c
c---------------------------------------------------------------------
c
c   latest revision    - Dec 08
c
c   purpose            - controls parametrization of prisms
c
c   arguments :
c     in:
c               No     - prism number
c               Eta    - reference coordinates of a point
c                        in the reference prism
c     out:
c               X      - physical coordinates of the point
c               Dxdeta - derivatives of the physical coordinates wrt
c                        to the reference coordinates
c
c   required  routines -
c
c---------------------------------------------------------------------
c
      subroutine prism_OLD(No,Eta, X,Dxdeta)
c
      use GMP
#include "syscom.blk"
      common /cprism/ iprint
      common /cprism_TI/ iprint_prism_T
c
      dimension Eta(3),X(3),Dxdeta(3,3)
c
c  ...vertex coordinates
      dimension xv(3,6)
c
c  ...order,linear shape functions
      dimension norder(19),vshap(6),dvshap(3,6)
c
c  ...edge and face orientations
      dimension nedge_orient(12),nface_orient(6)
      data nedge_orient /0,0,0,0,0,0,0,0,0,0,0,0/
      data nface_orient /0,0,0,0,0,0/
c
c-----------------------------------------------------------------------
c
ccc      iprint=0
      iprint_prism_T = iprint
c
      if (iprint.eq.1) then
        write(*,7001) No,PRISMS(No)%Type
 7001   format('prism: No,PRISMS(No)%Type = ',i8,2x,a10)
        call pause
      endif
c
      select case(PRISMS(No)%Type)
      case('Linear')
c
c  .....get the vertex coordinates
        do iv=1,6
          np = PRISMS(No)%VertNo(iv)
          call pointr(np, xv(1:3,iv))
        enddo
c
c  .....evaluate vertex shape functions
        norder(1:11)=1; norder(12:15) = 11
        call shapeHp(Eta,norder,nedge_orient,nface_orient,
     .               nrdof,vshap,dvshap)
c
        X(1:3)=0.d0; Dxdeta(1:3,1:3) = 0.d0
        do k=1,6
          X(1:3) = X(1:3) + xv(1:3,k)*vshap(k)
          do i=1,3
            Dxdeta(1:3,i) = Dxdeta(1:3,i) + xv(1:3,k)*dvshap(i,k)
          enddo
        enddo
      case('TIprism')
        call prism_TI(No,Eta, X,Dxdeta)
      case default
        write(*,7002) PRISMS(No)%Type
 7002   format('prism: unknown prism type = ',a10)
        stop 1
      end select
c
c
      end subroutine prism_OLD
c
c----------------------------------------------------------------------
c
c   routine name       - prism_TI
c
c---------------------------------------------------------------------
c
c   latest revision    - Feb 09
c
c   purpose            - routine defines the parameterization
c                        for a prism based on a transfinite
c                        interpolation (with linear and singular
c                        blending functions) of the paramaterization
c                        of its faces
c
c   arguments :
c     in:
c               No     - prism number
c               Eta    - reference coordinates of a point
c                        in the reference prism
c     out:
c               X      - physical coordinates of the point
c               Dxdeta - derivatives of the physical coordinates wrt
c                        to the reference coordinates
c
c---------------------------------------------------------------------
c
      subroutine prism_TI_OLD(No,Eta, X,Dxdeta)
c
      use control
      use GMP
      use element_data
#include "syscom.blk"
      common /cprism_TI/ iprint
      common /ctrianB/ iprint_trianB
c
      dimension Eta(3), X(3),Dxdeta(3,3)
c
c  ...derivatives of edge coordinate
      dimension dtedeta(3)
c
c  ...affine coordinates for triangular faces
      dimension vshapt(3),dvshapt(2,3)
c
c  ...1D vertex shape functions in the xi_3 direction
      dimension vshap(2),dvshap(2)
c
c  ...blending function
      dimension dblend(2)
c
c  ...edge kernels
      dimension xe(3),dxedt(3)
c
c  ...face kernels
      dimension xf(3),dxfdtf(3,2)
c
c-----------------------------------------------------------------------
c
ccc      iprint=0
      iprint_trianB = iprint
      if (iprint.eq.1) then
        write(*,7001) No,Eta(1:3)
 7001   format('prism_TI: No,Eta = ',i5,2x,3e12.5)
      endif
c
c  ...affine coordinates for the triangular faces
      vshapt(1) = 1.d0 - Eta(1) - Eta(2)
      dvshapt(1,1) = -1.d0; dvshapt(2,1) = -1.d0
      vshapt(2) = Eta(1)
      dvshapt(1,2) = 1.d0; dvshapt(2,2) = 0.d0
      vshapt(3) = Eta(2)
      dvshapt(1,3) = 0.d0; dvshapt(2,3) = 1.d0
c
c  ...1D shape functions in the xi_3 direction
      vshap(1) = 1.d0 - Eta(3); dvshap(1) = -1.d0
      vshap(2) = Eta(3); dvshap(2) = 1.d0
c
      X(1:3) = 0.d0; Dxdeta(1:3,1:3) = 0.d0
c
c  ...vertex interpolant
      iv=0
      do j=1,2
      do i=1,3
        iv=iv+1
        np = PRISMS(No)%VertNo(iv)
        X(1:3) = X(1:3) + POINTS(np)%Rdata(1:3)*vshapt(i)*vshap(j)
        do k=1,2
          Dxdeta(1:3,k) = Dxdeta(1:3,k)
     .                  + POINTS(np)%Rdata(1:3)*dvshapt(k,i)*vshap(j)
        enddo
        Dxdeta(1:3,3) = Dxdeta(1:3,3)
     .                + POINTS(np)%Rdata(1:3)*vshapt(i)*dvshap(j)
      enddo
      enddo
      if (iprint.eq.1) then
        write(*,*) 'prism_TI: VERTEX INTERPOLANT = '
        do ivar=1,3
          write(*,7011) X(ivar),Dxdeta(ivar,1:3)
 7011     format(e12.5,3x,3e12.5)
        enddo
      endif
c
c-----------------------------------------------------------------------
c
c  ...loop through horizontal edges
      ie=0
      do j=1,2
      do i=1,3
        ie=ie+1
c
        nc = PRISMS(No)%EdgeNo(ie); norient=0
        if (nc.lt.0) then
          nc = -nc; norient=1
        endif
        if (CURVES(nc)%Type.eq.'Seglin') cycle
c
c  .....get the edge vertices specifying the local edge orientation
        iv1 =  TRIAN_EDGE_TO_VERT(1,i); iv2 =  TRIAN_EDGE_TO_VERT(2,i)
c
c  .....project Eta(1:2) onto the edge
        call proj_t2e(Eta(1:2),iv1,iv2,vshapt,dvshapt, te,dtedeta(1:2))
        if ((abs(te).lt.GEOM_TOL).or.(abs(1.d0-te).lt.GEOM_TOL)) cycle
c
        if (iprint.eq.1) then
          write(*,7012) ie,nc,CURVES(nc)%Type
 7012     format('prism_TI: ie,nc,Type = ',i2,i5,2x,a5)
        endif
c
c  .....evaluate edge kernel function
        call curveK(nc,te,norient, xe,dxedt)
c
c  .....2D blending function
        blend = vshapt(iv1)*vshapt(iv2)
        dblend(1:2) = dvshapt(1:2,iv1)*vshapt(iv2)
     .              + vshapt(iv1)*dvshapt(1:2,iv2)
c
c  .....add edge contribution
        X(1:3) = X(1:3) + xe(1:3)*blend*vshap(j)
        do k=1,2
          Dxdeta(1:3,k) = Dxdeta(1:3,k)
     .                  + dxedt(1:3)*dtedeta(k)*blend*vshap(j)
     .                  + xe(1:3)*dblend(k)*vshap(j)
        enddo
        Dxdeta(1:3,3) = Dxdeta(1:3,3)
     .                + xe(1:3)*blend*dvshap(j)
      if (iprint.eq.1) then
        write(*,*) 'prism_TI: AFTER HORIZONTAL EDGE i,j',i,j
        do ivar=1,3
          write(*,7011) X(ivar),Dxdeta(ivar,1:3)
        enddo
      endif
      enddo
      enddo
c
c-----------------------------------------------------------------------
c
c  ...loop through vertical edges
      do i=1,3
        ie=ie+1
c
        nc = PRISMS(No)%EdgeNo(ie); norient=0
        if (nc.lt.0) then
          nc = -nc; norient=1
        endif
        if (CURVES(nc)%Type.eq.'Seglin') cycle
c
c  .....evaluate edge bubble function
        call curveB(nc,Eta(3),norient, xe,dxedt)
c
c  .....2D blending function
        blend = vshapt(i)
        dblend(1:2) = dvshapt(1:2,i)
c
c  .....add edge contribution
        X(1:3) = X(1:3) + xe(1:3)*blend
        do k=1,2
          Dxdeta(1:3,k) = Dxdeta(1:3,k) + xe(1:3)*dblend(k)
        enddo
        Dxdeta(1:3,3) = Dxdeta(1:3,3) + dxedt(1:3)*blend
      enddo
      if (iprint.eq.1) then
        write(*,*) 'prism_TI: AFTER VERTICAL EDGES = '
        do ivar=1,3
          write(*,7011) X(ivar),Dxdeta(ivar,1:3)
        enddo
      endif
c
c-----------------------------------------------------------------------
c
c  ...loop through horizontal faces
      if=0
      do j=1,2
        if=if+1
        call decode(PRISMS(No)%FigNo(if), nt,norient)
        if ((TRIANGLES(nt)%Type.eq.'TransTri').or.
     .      (TRIANGLES(nt)%Type.eq.'PlaneTri')) cycle
        if (iprint.eq.1) then
          write(*,7013) if,nt,TRIANGLES(nt)%Type
 7013     format('prism_TI: if,nt,Type = ',i2,i5,2x,a5)
        endif
c
c  .....compute the face bubble
        call trianB(nt,Eta(1:2),norient, xf,dxfdtf)
        if (iprint.eq.1) then
          do ivar=1,3
            write(*,7033) ivar,xf(ivar),dxfdtf(ivar,1:2)
 7033       format('prism_TI: ivar,xf,dxdtf = ',i2,2x,e12.5,2x,2e12.5)
          enddo
        endif
c
c  .....add face contribution
        X(1:3) = X(1:3) + xf(1:3)*vshap(j)
        do k=1,2
          Dxdeta(1:3,k) = Dxdeta(1:3,k)
     .                  + dxfdtf(1:3,k)*vshap(j)
        enddo
        Dxdeta(1:3,3) = Dxdeta(1:3,3)
     .                + xf(1:3)*dvshap(j)
      enddo
      if (iprint.eq.1) then
        write(*,*) 'prism_TI: AFTER HORIZONTAL FACES = '
        do ivar=1,3
          write(*,7011) X(ivar),Dxdeta(ivar,1:3)
        enddo
        call pause
      endif
c
c-----------------------------------------------------------------------
c
c  ...loop through vertical faces
      do ii=1,3
        if=if+1
        call decode(PRISMS(No)%FigNo(if), nr,norient)
        if ((RECTANGLES(nr)%Type.eq.'BilQua').or.
     .      (RECTANGLES(nr)%Type.eq.'TraQua').or.
     .      (RECTANGLES(nr)%Type.eq.'PTIRec')) cycle
        if (iprint.eq.1) then
          write(*,7014) if,nr,RECTANGLES(nr)%Type
 7014     format('prism_TI: if,nr,Type = ',i2,i5,2x,a5)
        endif
c
c  .....get the edge vertices specifying the local edge orientation
        iv1 =  TRIAN_EDGE_TO_VERT(1,ii); iv2 =  TRIAN_EDGE_TO_VERT(2,ii)
c
c  .....project Eta(1:2) onto the edge
        call proj_t2e(Eta(1:2),iv1,iv2,vshapt,dvshapt, te,dtedeta(1:2))
c
c  .....to be finished....
        write(*,*) 'prism_TI: UNFINISHED'
        stop2
      enddo
c
c
      end subroutine prism_TI_OLD

