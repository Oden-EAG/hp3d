c----------------------------------------------------------------------
c
c   routine name       - refel
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 2023
c
c   purpose            - for an element, routine returns reference
c                        coordinates of its vertices wrt particular
c                        GMP block, the element is lying on
c
c   arguments :
c     in:
c             Mdle     - an element middle node, identified with
c                        the element
c     out:
c             Iflag    = 5 GMP prism
c                      = 6 GMP hexa
c                      = 7 GMP tetra
c                      = 8 GMP pyramid
c             No       - GMP block number
c             Xsubel   - GMP block reference coordinates
c                        for the element vertices
c
c----------------------------------------------------------------------
c
      subroutine refel(Mdle, Iflag,No,Xsubel)
c
      use parameters
      use refinements
      use element_data
      use data_structure3D
c
      implicit none
c
      integer, intent(in)  :: Mdle
      integer, intent(out) :: Iflag,No
      real(8), intent(out) :: Xsubel(3,8)
c
      real(8) :: x(3,8),x_new(3,8)
c
c  ...history information: father, its type, refinement kind, son number
      integer :: ftype
      integer :: nfather(MAXGEN),father_type(MAXGEN),son_type(MAXGEN),
     .           nfather_ref_kind(MAXGEN),noson(MAXGEN)
c
      integer :: nson,nfath,nrgen,nel
      integer :: i,j,ie,is,iv,ifc,igen,jp,kref1,kref2,kref3
      integer :: nrsons,nv1,nv2,nv3,nv4,lab
c
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
c
c----------------------------------------------------------------------
c
      Xsubel = 0.d0
c
#if HP3D_DEBUG
      if (iprint.ge.1) then
        write(*,7000)mdle
 7000   format(' refel: mdle = ',i10)
      endif
#endif
c
c  ...go up the tree to find the initial mesh ancestor
      nson = Mdle; nfath = NODES(nson)%father; igen=0
      do while(nfath.gt.0)
        igen = igen + 1
        nfather(igen) = nfath
        father_type(igen) = NODES(nfath)%ntype
        son_type(igen) = NODES(nson)%ntype
        nfather_ref_kind(igen) = NODES(nfath)%ref_kind
        call nr_mdle_sons(NODES(nfath)%ntype, NODES(nfath)%ref_kind,
     .                    nrsons)
c        call locate(nson,NODES(nfath)%sons,nrsons, noson(igen))
        noson(igen) = nson - NODES(nfath)%first_son + 1
c        if (noson(igen)<0 .or. noson(igen)>nrsons) call pause
c
c  .....extensive printing
#if HP3D_DEBUG
        if (iprint.ge.2) then
          write(*,7002)igen,nfath,nson,S_Type(father_type(igen)),
     .                 S_Type(son_type(igen)),
     .                 nfather_ref_kind(igen),noson(igen)
7002      format(' igen,nfath,nson,nfath_type,nson_type,nfath_refk,
     .ison = ',i2,4x,i8,2x,i8,2x,a4,2x,a4,2x,i2,2x,i2)
        endif
#endif
c
c  .....update
        nson = nfath
        nfath = NODES(nson)%father
      enddo
      nrgen = igen
c
c  ...initial mesh ancestor is
      nel = abs(nfath)
c
c  ...reference coordinates of ancestor
      select case(NODES(nel)%ntype)
      case(MDLB); x(1:3,1:8) = BRICK_COORD(1:3,1:8)
      case(MDLP); x(1:3,1:6) = PRISM_COORD(1:3,1:6)
      case(MDLN); x(1:3,1:4) = TETRA_COORD(1:3,1:4)
      case(MDLD); x(1:3,1:5) = PYRAM_COORD(1:3,1:5)
      end select
c
#if HP3D_DEBUG
c  ...printing
      if (iprint.ge.2) then
        do i=1,Nvert(NODES(nel)%ntype)
          write(*,8005)i,x(1:3,i)
        enddo
8005    format( 'i,x(1:3,i) = ',i2,2x,3(e12.5,2x))
      endif
#endif
c
c  ...corresponding GMP block
      call decode(ELEMS(nel)%GMPblock, No,lab)
      Iflag=4+lab
c
c  ...go down the tree reconstructing reference coordinates
      do igen=nrgen,1,-1
        ftype = father_type(igen)
c
c  .....loop through son's vertex nodes
        do j=1,nvert(son_type(igen))
c
c  .......decode father refinement flag
          call decode_ref(ftype,nfather_ref_kind(igen),
     .                    kref1,kref2,kref3)
c
c  .......find the parent node of the vertex
          jp = npar_ref(ftype,j,noson(igen),kref1,kref2,kref3)
c
c  .......find the son number of the vertex
          is = nson_ref(ftype,j,noson(igen),kref1,kref2,kref3)
c
c  .......node shared with the father
          if (is.eq.0) then
            if (jp.eq.0) jp=j
            x_new(1:3,j) = x(1:3,jp)
c
#if HP3D_DEBUG
            if (iprint.ge.2) then
              write(*,*)'node shared with the father: jp = ',jp
              write(*,8000)j,x_new(1:3,j)
            endif
#endif
c
c  .......node generated through refinement
          else
c
c  .........parent edge node
            if (jp.le.nvert(ftype)+nedge(ftype)) then
              ie = jp-nvert(ftype)
              call edge_to_vert(ftype,ie, nv1,nv2)
              x_new(1:3,j) = (x(1:3,nv1) + x(1:3,nv2))/2.d0
c
#if HP3D_DEBUG
              if (iprint.ge.2) then
                write(*,*)'parent edge node: nv1,nv2 = ',nv1,nv2
                write(*,8000)j,x_new(1:3,j)
              endif
#endif
c
c  .........parent rectangular face node
            elseif (jp.le.nvert(ftype)+nedge(ftype)+nface(ftype))
     .        then
              ifc = jp-nvert(ftype)-nedge(ftype)
              call face_to_vert(ftype,ifc, nv1,nv2,nv3,nv4)
              x_new(1:3,j) = (x(1:3,nv1) + x(1:3,nv2)
     .                      + x(1:3,nv3) + x(1:3,nv4))/4.d0
c
#if HP3D_DEBUG
              if (iprint.ge.2) then
                write(*,*)'parent rectangular face node: nv1,nv2,nv3,
     . nv4 = ',nv1,nv2,nv3,nv4
                write(*,8000)j,x_new(1:3,j)
              endif
#endif
c
c  .........parent middle node (brick element only)
            else
              x_new(1:3,j) = 0.d0
              do iv=1,8
                x_new(1:3,j) = x_new(1:3,j) + x(1:3,iv)
              enddo
              x_new(1:3,j) = x_new(1:3,j)/8.d0
c
#if HP3D_DEBUG
              if (iprint.ge.2) then
                write(*,*)'parent middle node (brick only!)'
                write(*,8000)j,x_new(1:3,j)
              endif
#endif
c
            endif
          endif
c
 8000   format(' j,x_new(1:3,j) = ',i2,2x,3(e12.5,2x))
c
c  .....end of loop through vertex nodes
        enddo
c
#if HP3D_DEBUG
c  .....printing
        if (iprint.ge.2) then
          write(*,7003)igen,Nvert(son_type(igen))
 7003     format(' igen,nvert_son = ',i2,2x,i2)
          do i=1,Nvert(son_type(igen))
            write(*,7004)i,x_new(1:3,i)
 7004       format(' i,x_new(1:3,i) = ',i2,2x,3(e12.5,2x))
          enddo
        endif
#endif
c
c  .....update
        x = x_new
c
c  ...end of loop through generations
      enddo
c
      Xsubel(1:3,1:nvert(NODES(Mdle)%ntype)) =
     .     x(1:3,1:nvert(NODES(Mdle)%ntype))
c
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7001) Mdle,Xsubel(1:3,1:nvert(NODES(Mdle)%ntype))
 7001   format('refel: Mdle = ',i6,' Xsubel = ',8(/,3(f8.3,2x)))
        call pause
      endif
#endif
c
c
      end subroutine refel
