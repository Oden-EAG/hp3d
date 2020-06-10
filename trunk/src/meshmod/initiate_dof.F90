#include "typedefs.h"
!----------------------------------------------------------------------------------------
!> Purpose : initiate H1, H(curl) and H(div), and L2 dof for the sons of 
!            the middle node of an element
!
!> @date Jun 20
!
!> @param[in]  Mdle - mdle node number
!
!> @param[out] dof stored directly in the data structure arrays
!
!----------------------------------------------------------------------------------------
!
      subroutine initiate_dof(Mdle,Xnod_fath,ZdofH_fath,ZdofE_fath,ZdofV_fath,ZdofQ_fath, &
                                   NH,NE,NV,NQ)
!
      use refinements_geometry
      use refinements
      use data_structure3D
      implicit none
!
      integer, intent(in)  :: Mdle,NH,NE,NV,NQ
      real*8, dimension(NDIMEN,  NH),  intent(in)  :: Xnod_fath
      VTYPE,  dimension(MAXEQNH, NH),  intent(in)  :: ZdofH_fath
      VTYPE,  dimension(MAXEQNE, NE),  intent(in)  :: ZdofE_fath
      VTYPE,  dimension(MAXEQNV, NV),  intent(in)  :: ZdofV_fath
      VTYPE,  dimension(MAXEQNV, NQ),  intent(in)  :: ZdofQ_fath
!
!  ...locals
      real*8, dimension(NDIMEN,  MAXbrickH)  :: xnod
      VTYPE,  dimension(MAXEQNH, MAXbrickH)  :: zdofH
      VTYPE,  dimension(MAXEQNE, MAXbrickE)  :: zdofE
      VTYPE,  dimension(MAXEQNV, MAXbrickV)  :: zdofV
      VTYPE,  dimension(MAXEQNV, MAXbrickQ)  :: zdofQ
!
      real*8, allocatable, dimension(:)   :: aloc
      real*8, allocatable, dimension(:,:) :: bloc
!
      character(len=4) :: type_fath, type
!
!  ...element order, orientation for edges and faces
      integer :: kref_fath, nrs,is,mdle_son, nint, l, iflag, k1, k2, nk, k, info, ivar, iprint, &
                 norder_fath(19), norient_edge_fath(12), norient_face_fath(6), &
                 norder(19), norient_edge(12), norient_face(6), norder1(19)
!
!  ...a son's vertex coordinates wrt the father
      real*8 :: xvert(NDIMEN,8), wa, weight
!
!  ...geometry
      real*8 :: xi(3), x(3), dxdxi(3,3), dxidx(3,3), rjac
!
!  ...H1 shape functions
      integer :: nrdofH_fath, nrdofH
      real*8  :: shapH_fath(MAXbrickH), gradH_fath(3,MAXbrickH)
      real*8  :: shapH(MAXbrickH), gradH(3,MAXbrickH)
!
!  ...H(curl) shape functions
      integer :: nrdofE_fath, nrdofE  
      real*8  :: shapE_fath(3,MAXbrickE), curlE_fath(3,MAXbrickE)
      real*8  :: shapE(3,MAXbrickE), curlE(3,MAXbrickE), E1(3),E2(3)
!
!  ...H(div) shape functions
      integer :: nrdofV_fath, nrdofV
      real*8  :: shapV_fath(3,MAXbrickV), divV_fath(MAXbrickV)
      real*8  :: shapV(3,MAXbrickV), divV(MAXbrickV),  V1(3),V2(3)
!
!  ...L2 shape functions
      integer :: nrdofQ_fath, nrdofQ
      real*8  :: shapQ_fath(MAXbrickQ)
      real*8  :: shapQ(MAXbrickQ)
!
!  ...3D quadrature data
      real*8 :: xiloc(3,MAX_NINT3), waloc(MAX_NINT3)
!
!  ...solutions to be projected
      real*8  :: xp(NDIMEN)
      VTYPE   :: zsolH(MAXEQNH), zsolE(3,MAXEQNE), zsolV(3,MAXEQNV), zsolQ(MAXEQNQ) 
!
!  ...LAPACK solver
      character uplo
!
!  ...element nodes and orientations
      integer :: nodesl(27), norientl(27), i,ibeg,iend,jbeg,jend
!
      integer :: nrv,nre,nrf,nrn,j,nod, &
                 ndofH,ndofE,ndofV,ndofQ, nvarH,nvarE,nvarV,nvarQ
!
!  ...decoded index for a node
      integer :: index(NRINDEX)
!
!  ...work space for printing
      real*8 :: aux(100)
!
!  ...function for vector storage of Hermitian Gram matrix
      nk(k1,k2) = (k2-1)*k2/2+k1
!
!----------------------------------------------------------------------
!
      iprint=0
      nvarH = NRHVAR*NRCOMS; nvarE = NREVAR*NRCOMS
      nvarV = NRVVAR*NRCOMS; nvarQ = NRQVAR*NRCOMS
!
!  ...element type
      type_fath = NODES(Mdle)%type
      kref_fath = NODES(Mdle)%ref_kind
      call nr_mdle_sons(type_fath,kref_fath, nrs)
!
      call find_order(Mdle, norder_fath)
      call find_orient(Mdle, norient_edge_fath,norient_face_fath)
!
      if (iprint.eq.1) then
        call celndof(type_fath,norder_fath, nrdofH_fath,nrdofE_fath,nrdofV_fath,nrdofQ_fath)
        write(*,*) 'initiate_dof: Xnod_fath = '
        do ivar=1,NDIMEN
          write(*,7010) Xnod_fath(ivar,1:nrdofH_fath)
 7010     format(30f7.2)
        enddo
        write(*,*) 'initiate_dof: ZdofH_fath = '
        do ivar=1,nvarH
          write(*,7010) ZdofH_fath(ivar,1:nrdofH_fath)
        enddo
        write(*,*) 'initiate_dof: ZdofE_fath = '
        do ivar=1,nvarE
          write(*,7010) ZdofE_fath(ivar,1:nrdofE_fath)
        enddo
        write(*,*) 'initiate_dof: ZdofV_fath = '
        do ivar=1,nvarV
          write(*,7010) ZdofV_fath(ivar,1:nrdofV_fath)
        enddo
        write(*,*) 'initiate_dof: ZdofQ_fath = '
        do ivar=1,nvarQ
          write(*,7010) ZdofQ_fath(ivar,1:nrdofQ_fath)
        enddo
        call pause
      endif
!
!  ...loop through the element middle node sons
      do is=1,nrs
        mdle_son = son(mdle,is)
        type = NODES(mdle_son)%type
        call find_order(mdle_son, norder)
        call find_orient(mdle_son, norient_edge,norient_face)
        call set_3Dint(type,norder, nint,xiloc,waloc)
        call get_son_coord(type_fath,kref_fath,is, nrv,xvert)
        if (iprint.eq.1) write(*,*) 'initiate_dof: is,mdle_son = ',is,mdle_son
        call celndof(type,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
        call set_linear_order(type, norder1)
!
!  *********************************** H1 variables ***********************************
!
        allocate(aloc(nrdofH*(nrdofH+1)))
        aloc = 0.d0
#if C_MODE
        allocate(bloc(nrdofH,NDIMEN+2*nvarH))
#else
        allocate(bloc(nrdofH,NDIMEN+nvarH))
#endif
        bloc = 0.d0
!
!  .....loop over integration points
        do l=1,nint
!
!  .......coordinates and weight of this integration point
          xi(1:3) = xiloc(1:3,l)
          wa=waloc(l)
          call shape3DH(type,xi,norder,norient_edge,norient_face, nrdofH,shapH,gradH)
          call geom3D(mdle_son,xi,xvert,shapH,gradH,nrv, x,dxdxi,dxidx,rjac,iflag)
!
!  .......integration weight
          weight = rjac*wa
!
!  .......compute the father solution to be projected
          call shape3DH(type_fath,x,norder_fath,norient_edge_fath,norient_face_fath, &
                        nrdofH_fath,shapH_fath,gradH_fath)
          xp  = 0.d0; zsolH = ZERO
          do k1=1,nrdofH_fath
            xp(1:NDIMEN)  = xp(1:NDIMEN) + Xnod_fath(1:NDIMEN,k1)*shapH_fath(k1)
            zsolH(1:MAXEQNH) = zsolH(1:MAXEQNH) + ZdofH_fath(1:MAXEQNH,k1)*shapH_fath(k1)
          enddo
!
!  .......first loop through dof
          do k1=1,nrdofH
!
!  .........accumulate for the load vector
            bloc(k1,1:NDIMEN) = bloc(k1,1:NDIMEN) + xp(1:NDIMEN)*shapH(k1)*weight
#if C_MODE
            do ivar=1,nvarH
              bloc(k1,NDIMEN+ivar) = bloc(k1,NDIMEN+ivar) &
                                    + dreal(zsolH(ivar))*shapH(k1)*weight
              bloc(k1,NDIMEN+nvarH+ivar) = bloc(k1,NDIMEN+nvarH+ivar) &
                                    + dimag(zsolH(ivar))*shapH(k1)*weight
            enddo
#else
            do ivar=1,nvarH
              bloc(k1,NDIMEN+ivar) = bloc(k1,NDIMEN+ivar) + zsolH(ivar)*shapH(k1)*weight
            enddo
#endif
!
!  .........second loop through dof
            do k2=k1,nrdofH
!
!  ...........accumulate for the stiffness matrix
              k = nk(k1,k2)
              aloc(k) = aloc(k) + shapH(k1)*shapH(k2)*weight
            enddo
          enddo
       enddo
!
      if (iprint.eq.2) then
  100   write(*,7011)  nrdofH
 7011   format('initiate_dof: nrdofH = ',i2,/,'SET ibeg,iend,jbeg,jend for aloc')
        read(*,*) ibeg,iend,jbeg,jend
        if (ibeg.ne.0) then
          do j=jbeg,jend
            write(*,7110) j
 7110       format('j = ',i3)
            k=0
            do i=ibeg,iend
              k=k+1
              if (i.ge.j) then
                l=nk(j,i)
              else
                l=nk(i,j)
              endif
              aux(k) =  aloc(l)
            enddo
            write(*,7015) aux(1:k)
 7015       format(10e12.5)
          enddo
          go to 100
        endif
      endif
!
!  ....factorize
       uplo = 'U'
       call DPPTRF(uplo, nrdofH, aloc, info) 
#if C_MODE
       call DPPTRS(uplo, nrdofH, NDIMEN+2*nvarH, aloc, bloc, nrdofH, info)
       do k=1,nrdofH
         xnod(1:NDIMEN,k) = bloc(k,1:NDIMEN)
         do j=1,nvarH
           zdofH(j,k) = cmplx(bloc(k,NDIMEN+j),bloc(k,NDIMEN+nvarH+j))
         enddo
       enddo
#else
       call DPPTRS(uplo, nrdofH, NDIMEN+nvarH, aloc, bloc, nrdofH, info)
       do k=1,nrdofH
         xnod(1:NDIMEN,k) = bloc(k,1:NDIMEN)
         zdofH(1:nvarH,k) = bloc(k,NDIMEN+1:NDIMEN+nvarH)
       enddo
#endif
       deallocate(aloc,bloc)       
       if (iprint.eq.1) write(*,*) 'initiate_dof: DONE WITH H1 DOF'
!
!  *********************************** H(curl) variables ***********************************
!
        allocate(aloc(nrdofE*(nrdofE+1)))
        aloc = 0.d0
#if C_MODE
        allocate(bloc(nrdofE,2*nvarE))
#else
        allocate(bloc(nrdofE,nvarE))
#endif
        bloc = 0.d0
!
!  .....loop over integration points
        do l=1,nint
!
!  .......coordinates and weight of this integration point
          xi(1:3) = xiloc(1:3,l)
          wa=waloc(l)
          call shape3DH(type,xi,norder1,norient_edge,norient_face, nrdofH,shapH,gradH)
          call shape3DE(type,xi,norder, norient_edge,norient_face, nrdofE,shapE,curlE)
          call geom3D(mdle_son,xi,xvert,shapH,gradH,nrv, x,dxdxi,dxidx,rjac,iflag)
!
!  .......integration weight
          weight = rjac*wa
!
!  .......compute the father solution to be projected
          call shape3DE(type_fath,x,norder_fath,norient_edge_fath,norient_face_fath, &
                        nrdofE_fath,shapE_fath,curlE_fath)
          zsolE = ZERO
          do k1=1,nrdofE_fath
            do ivar=1,MAXEQNE
              zsolE(1:3,ivar) = zsolE(1:3,ivar) + ZdofE_fath(ivar,k1)*shapE_fath(1:3,k1)
            enddo
          enddo
!
!  .......first loop through dof
          do k1=1,nrdofE
!
! ..........Piola transform
            E1(1:3) = shapE(1,k1)*dxidx(1,1:3)  &
                    + shapE(2,k1)*dxidx(2,1:3)  &
                    + shapE(3,k1)*dxidx(3,1:3)
!
!  .........accumulate for the load vector
#if C_MODE
            do ivar=1,nvarE
              bloc(k1,ivar) = bloc(k1,ivar) &
              + dreal(zsolE(1,ivar)*E1(1)+zsolE(2,ivar)*E1(2)+zsolE(3,ivar)*E1(3))*weight
              bloc(k1,nvarE+ivar) = bloc(k1,nvarE+ivar) &
              + dimag(zsolE(1,ivar)*E1(1)+zsolE(2,ivar)*E1(2)+zsolE(3,ivar)*E1(3))*weight
            enddo
#else
            do ivar=1,nvarE
              bloc(k1,ivar) = bloc(k1,ivar) &
              + (zsolE(1,ivar)*E1(1)+zsolE(2,ivar)*E1(2)+zsolE(3,ivar)*E1(3))*weight
            enddo
!!            write(*,*)'l,k1,bloc(k1,1:nvarE) = ',l,k1,bloc(k1,1:nvarE)
#endif
!
!  .........second loop through dof
            do k2=k1,nrdofE
!
! ............Piola transform
              E2(1:3) = shapE(1,k2)*dxidx(1,1:3)  &
                      + shapE(2,k2)*dxidx(2,1:3)  &
                      + shapE(3,k2)*dxidx(3,1:3)
!
!  ...........accumulate for the stiffness matrix
              k = nk(k1,k2)
              aloc(k) = aloc(k) + (E1(1)*E2(1) + E1(2)*E2(2) + E1(3)*E2(3))*weight
            enddo
          enddo
       enddo
!
      if (iprint.eq.3) then
  200   write(*,7021)  nrdofE
 7021   format('initiate_dof: nrdofE = ',i2,/,'SET ibeg,iend,jbeg,jend for aloc')
        read(*,*) ibeg,iend,jbeg,jend
        if (ibeg.ne.0) then
          do j=jbeg,jend
            write(*,7110) j
            k=0
            do i=ibeg,iend
              k=k+1
              if (i.ge.j) then
                l=nk(j,i)
              else
                l=nk(i,j)
              endif
              aux(k) =  aloc(l)
            enddo
            write(*,7015) aux(1:k)
          enddo
          go to 200
        endif
      endif
!
!  ....factorize
       uplo = 'U'
       call DPPTRF(uplo, nrdofE, aloc, info) 
#if C_MODE
       call DPPTRS(uplo, nrdofE, 2*nvarE, aloc, bloc, nrdofE, info)
       do k=1,nrdofE
         do j=1,nvarE
           zdofE(j,k) = cmplx(bloc(k,j),bloc(k,nvarE+j))
         enddo
       enddo
#else
       call DPPTRS(uplo, nrdofE, nvarE, aloc, bloc, nrdofE, info)
       do k=1,nrdofE
         zdofE(1:nvarE,k) = bloc(k,1:nvarE)
       enddo
#endif
       deallocate(aloc,bloc)       
       if (iprint.eq.1) write(*,*) 'initiate_dof: DONE WITH H(curl) DOF'
!
!
!  *********************************** H(div) variables ***********************************
!
        allocate(aloc(nrdofV*(nrdofV+1)))
        aloc = 0.d0
#if C_MODE
        allocate(bloc(nrdofV,2*nvarV))
#else
        allocate(bloc(nrdofV,nvarV))
#endif
        bloc = 0.d0
!
!  .....loop over integration points
        do l=1,nint
!
!  .......coordinates and weight of this integration point
          xi(1:3) = xiloc(1:3,l)
          wa=waloc(l)
          call shape3DH(type,xi,norder1,norient_edge,norient_face, nrdofH,shapH,gradH)
          call shape3DV(type,xi,norder,norient_face, nrdofV,shapV,divV)
          call geom3D(mdle_son,xi,xvert,shapH,gradH,nrv, x,dxdxi,dxidx,rjac,iflag)
!
!  .......integration weight
          weight = rjac*wa
!
!  .......compute the father solution to be projected
          call shape3DV(type_fath,x,norder_fath,norient_face_fath, &
                        nrdofV_fath,shapV_fath,divV_fath)
          zsolV = ZERO
          do k1=1,nrdofV_fath
            do ivar=1,MAXEQNV
              zsolV(1:3,ivar) = zsolV(1:3,ivar) + ZdofV_fath(ivar,k1)*shapV_fath(1:3,k1)
            enddo
          enddo
!
!  .......first loop through dof
          do k1=1,nrdofV
!
! ..........Piola transform
            V1(1:3) = (dxdxi(1:3,1)*shapV(1,k1)   &
                    +  dxdxi(1:3,2)*shapV(2,k1)   &
                    +  dxdxi(1:3,3)*shapV(3,k1))/rjac
!
!  .........accumulate for the load vector
#if C_MODE
            do ivar=1,nvarV
              bloc(k1,ivar) = bloc(k1,ivar) &
              + dreal(zsolV(1,ivar)*V1(1)+zsolV(2,ivar)*V1(2)+zsolV(3,ivar)*V1(3))*weight
              bloc(k1,nvarV+ivar) = bloc(k1,nvarV+ivar) &
              + dimag(zsolV(1,ivar)*V1(1)+zsolV(2,ivar)*V1(2)+zsolV(3,ivar)*V1(3))*weight
            enddo
#else
            do ivar=1,nvarV
              bloc(k1,ivar) = bloc(k1,ivar) &
              + (zsolV(1,ivar)*V1(1)+zsolV(2,ivar)*V1(2)+zsolV(3,ivar)*V1(3))*weight
            enddo
#endif
!
!  .........second loop through dof
            do k2=k1,nrdofV
!
! ............Piola transform
              V2(1:3) = (dxdxi(1:3,1)*shapV(1,k2)   &
                      +  dxdxi(1:3,2)*shapV(2,k2)   &
                      +  dxdxi(1:3,3)*shapV(3,k2))/rjac
!
!  ...........accumulate for the stiffness matrix
              k = nk(k1,k2)
              aloc(k) = aloc(k) + (V1(1)*V2(1) + V1(2)*V2(2) + V1(3)*V2(3))*weight
            enddo
          enddo
       enddo
!
!  ....factorize
       uplo = 'U'
       call DPPTRF(uplo, nrdofV, aloc, info) 
#if C_MODE
       call DPPTRS(uplo, nrdofV, 2*nvarV, aloc, bloc, nrdofV, info)
       do k=1,nrdofV
         do j=1,nvarV
           zdofV(j,k) = cmplx(bloc(k,j),bloc(k,nvarV+j))
         enddo
       enddo
#else
       call DPPTRS(uplo, nrdofV, nvarV, aloc, bloc, nrdofV, info)
       do k=1,nrdofV
         zdofV(1:nvarV,k) = bloc(k,1:nvarV)
       enddo
#endif
       deallocate(aloc,bloc)       
       if (iprint.eq.1) write(*,*) 'initiate_dof: DONE WITH H(div) DOF'
!
!  *********************************** L2 variables ***********************************
!
        allocate(aloc(nrdofQ*(nrdofQ+1)))
        aloc = 0.d0
#if C_MODE
        allocate(bloc(nrdofQ,2*nvarQ))
#else
        allocate(bloc(nrdofQ,nvarQ))
#endif
        bloc = 0.d0
!
!  .....loop over integration points
        do l=1,nint
!
!  .......coordinates and weight of this integration point
          xi(1:3) = xiloc(1:3,l)
          wa=waloc(l)
          call shape3DH(type,xi,norder1,norient_edge,norient_face, nrdofH,shapH,gradH)
          call shape3DQ(type,xi,norder, nrdofQ,shapQ)
          call geom3D(mdle_son,xi,xvert,shapH,gradH,nrv, x,dxdxi,dxidx,rjac,iflag)
!
!  .......integration weight
          weight = rjac*wa
!
!  .......compute the father solution to be projected
          call shape3DQ(type_fath,x,norder_fath, nrdofQ_fath,shapQ_fath)
          zsolQ = ZERO
          do k1=1,nrdofQ_fath
            zsolQ(1:MAXEQNQ) = zsolQ(1:MAXEQNQ) + ZdofQ_fath(1:MAXEQNQ,k1)*shapQ_fath(k1)
          enddo
!
!  .......first loop through dof
          do k1=1,nrdofQ
!
!  .........accumulate for the load vector
#if C_MODE
            do ivar=1,nvarQ
              bloc(k1,ivar) = bloc(k1,ivar) &
                            + dreal(zsolQ(ivar))*shapQ(k1)/rjac*weight
              bloc(k1,nvarQ+ivar) = bloc(k1,nvarQ+ivar) &
                                  + dimag(zsolH(ivar))*shapQ(k1)/rjac*weight
            enddo
#else
            do ivar=1,nvarQ
              bloc(k1,ivar) = bloc(k1,ivar) + zsolQ(ivar)*shapQ(k1)/rjac*weight
            enddo
#endif
!
!  .........second loop through dof
            do k2=k1,nrdofQ
!
!  ...........accumulate for the stiffness matrix
              k = nk(k1,k2)
              aloc(k) = aloc(k) + shapQ(k1)*shapQ(k2)/rjac**2*weight
            enddo
          enddo
       enddo
!
!  ....factorize
       uplo = 'U'
       call DPPTRF(uplo, nrdofQ, aloc, info) 
#if C_MODE
       call DPPTRS(uplo, nrdofQ, 2*nvarQ, aloc, bloc, nrdofQ, info)
       do k=1,nrdofQ
         do j=1,nvarQ
           zdofQ(j,k) = cmplx(bloc(k,j),bloc(k,nvarQ+j))
         enddo
       enddo
#else
       call DPPTRS(uplo, nrdofQ, nvarQ, aloc, bloc, nrdofQ, info)
       do k=1,nrdofQ
         zdofQ(1:nvarQ,k) = bloc(k,1:nvarQ)
       enddo
#endif
       deallocate(aloc,bloc)       
       if (iprint.eq.1) write(*,*) 'initiate_dof: DONE WITH L2 DOF'
!
!*************************************************************************************
!
!  ....determine element nodes
       call elem_nodes(mdle_son, nodesl,norientl)
       nrn  = nvert(type)+nedge(type)+nface(type)+1
       if (iprint.eq.1) then
         write(*,7020) nodesl(1:nrn)
 7020    format(' initiate_dof: nodesl = ',8i5,3x,12i5,3x,6i5,3x,i5)
       endif
!
!  ....loop through the element nodes
       nrdofH = 0; nrdofE = 0; nrdofV = 0; nrdofQ = 0; 
       do j=1,nrn
         nod = nodesl(j)
         call ndof_nod(NODES(nod)%type,NODES(nod)%order, ndofH,ndofE,ndofV,ndofQ)
         if ((iprint.eq.1).and.(j.le.nrv)) then
           write(*,7030) j, xnod(1:NDIMEN,nrdofH+1:nrdofH+ndofH),XVERT_mdlb111(1:3,j,is)
 7030      format(' VERTEX = ',i1,' COORD = ',3f8.3,' XVERT = ',3f8.3)
         endif
         if (iprint.eq.1) then
           write(*,*) 'initiate_dof: DOF for j,nod = ',j,nod
           write(*,7040) zdofH(1,nrdofH+1:nrdofH+ndofH)
 7040      format(15e12.5)
           write(*,7040) zdofE(1,nrdofE+1:nrdofE+ndofE)
           write(*,7040) zdofV(1,nrdofV+1:nrdofV+ndofV)
           write(*,7040) zdofQ(1,nrdofQ+1:nrdofQ+ndofQ)
           call pause
         endif

!
!  ......save the new dof if the node is active
         if (NODES(nod)%act) then
            NODES(nod)%dof%coord(1:NDIMEN,1:ndofH) = xnod(1:NDIMEN,nrdofH+1:nrdofH+ndofH)
            call dof_in(nod,zdofH(1,nrdofH+1),zdofE(1,nrdofE+1),zdofV(1,nrdofV+1),zdofQ(1,nrdofQ+1))
         endif
         nrdofH = nrdofH + ndofH; nrdofE = nrdofE + ndofE
         nrdofV = nrdofV + ndofV; nrdofQ = nrdofQ + ndofQ 
!
!  ....end of loop through the element nodes
       enddo
!
!  ...end of loop through the element sons
      enddo
      if (iprint.eq.1) then
        write(*,*) 'initiate_dof: EXITING'
        call pause
      endif
!
!
      end subroutine initiate_dof



      subroutine set_linear_order(Type,Norder)
      character(len=4),  intent(in) :: Type
      integer,           intent(out) :: Norder(19)
      select case(Type)
      case('medg')
        Norder(1)=1
      case('mdlq')
        Norder(1:4)=1; Norder(5)=11
      case('mdlt')
        Norder(1:4)=1
      case('mdlb')
        Norder(1:12)=1; Norder(13:18)=11; Norder(19)=111
      case('mdlp')
        Norder(1:9)=1; Norder(10:11)=1; Norder(12:15)=11
      case('mdln')
        Norder(1:11)=1
      case('mdld')
        Norder(1:8)=1; Norder(9)=11; Norder(10:14)=1
      end select
      end subroutine set_linear_order
