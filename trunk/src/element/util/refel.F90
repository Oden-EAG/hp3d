!----------------------------------------------------------------------
!
!   routine name       - refel
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2024
!
!   purpose            - for an element, routine returns reference
!                        coordinates of its vertices wrt particular
!                        GMP block, the element is lying on
!
!   arguments :
!     in:
!             Mdle     - an element middle node, identified with
!                        the element
!     out:
!             Iflag    = 5 GMP prism
!                      = 6 GMP hexa
!                      = 7 GMP tetra
!                      = 8 GMP pyramid
!             No       - GMP block number
!             Xsubel   - GMP block reference coordinates
!                        for the element vertices
!
!----------------------------------------------------------------------
!
   subroutine refel(Mdle, Iflag,No,Xsubel)
!
      use parameters
      use refinements
      use element_data
      use data_structure3D
!
      implicit none
!
      integer, intent(in)  :: Mdle
      integer, intent(out) :: Iflag,No
      real(8), intent(out) :: Xsubel(3,8)
!
      real(8) :: x(3,8),x_new(3,8)
!
!  ...history information: father, its type, refinement kind, son number
      integer :: ftype
      integer :: nfather(MAXGEN),father_type(MAXGEN),son_type(MAXGEN), &
                 nfather_ref_kind(MAXGEN),noson(MAXGEN)
!
      integer :: nson,nfath,nrgen,nel
      integer :: j,ie,is,iv,ifc,igen,jp,kref1,kref2,kref3
      integer :: nrsons,nv1,nv2,nv3,nv4,lab
!
#if HP3D_DEBUG
      integer :: i
      integer :: iprint
      iprint=0
#endif
!
!----------------------------------------------------------------------
!
      Xsubel = 0.d0
!
#if HP3D_DEBUG
      if (iprint.ge.1) then
        write(*,7000)mdle
 7000   format(' refel: mdle = ',i10)
      endif
#endif
!
!  ...go up the tree to find the initial mesh ancestor
      nson = Mdle; nfath = NODES(nson)%father; igen=0
      do while(nfath.gt.0)
        igen = igen + 1
        nfather(igen) = nfath
        father_type(igen) = NODES(nfath)%ntype
        son_type(igen) = NODES(nson)%ntype
        nfather_ref_kind(igen) = NODES(nfath)%ref_kind
        call nr_mdle_sons(NODES(nfath)%ntype, NODES(nfath)%ref_kind,nrsons)
!        call locate(nson,NODES(nfath)%sons,nrsons, noson(igen))
        noson(igen) = nson - NODES(nfath)%first_son + 1
!        if (noson(igen)<0 .or. noson(igen)>nrsons) call pause
!
!  .....extensive printing
#if HP3D_DEBUG
        if (iprint.ge.2) then
          write(*,7002)igen,nfath,nson,S_Type(father_type(igen)), &
                       S_Type(son_type(igen)),                    &
                       nfather_ref_kind(igen),noson(igen)
7002      format(' igen,nfath,nson,nfath_type,nson_type,nfath_refk, &
                   ison = ',i2,4x,i8,2x,i8,2x,a4,2x,a4,2x,i2,2x,i2)
        endif
#endif
!
!  .....update
        nson = nfath
        nfath = NODES(nson)%father
      enddo
      nrgen = igen
!
!  ...initial mesh ancestor is
      nel = abs(nfath)
!
!  ...reference coordinates of ancestor
      select case(NODES(nel)%ntype)
      case(MDLB); x(1:3,1:8) = BRICK_COORD(1:3,1:8)
      case(MDLP); x(1:3,1:6) = PRISM_COORD(1:3,1:6)
      case(MDLN); x(1:3,1:4) = TETRA_COORD(1:3,1:4)
      case(MDLD); x(1:3,1:5) = PYRAM_COORD(1:3,1:5)
      end select
!
#if HP3D_DEBUG
!  ...printing
      if (iprint.ge.2) then
        do i=1,Nvert(NODES(nel)%ntype)
          write(*,8005)i,x(1:3,i)
        enddo
8005    format( 'i,x(1:3,i) = ',i2,2x,3(e12.5,2x))
      endif
#endif
!
!  ...corresponding GMP block
      call decode(ELEMS(nel)%GMPblock, No,lab)
      Iflag=4+lab
!
!  ...go down the tree reconstructing reference coordinates
      do igen=nrgen,1,-1
        ftype = father_type(igen)
!
!  .....loop through son's vertex nodes
        do j=1,nvert(son_type(igen))
!
!  .......decode father refinement flag
          call decode_ref(ftype,nfather_ref_kind(igen), &
                          kref1,kref2,kref3)
!
!  .......find the parent node of the vertex
          jp = npar_ref(ftype,j,noson(igen),kref1,kref2,kref3)
!
!  .......find the son number of the vertex
          is = nson_ref(ftype,j,noson(igen),kref1,kref2,kref3)
!
!  .......node shared with the father
          if (is.eq.0) then
            if (jp.eq.0) jp=j
            x_new(1:3,j) = x(1:3,jp)
!
#if HP3D_DEBUG
            if (iprint.ge.2) then
              write(*,*)'node shared with the father: jp = ',jp
              write(*,8000)j,x_new(1:3,j)
            endif
#endif
!
!  .......node generated through refinement
          else
!
!  .........parent edge node
            if (jp.le.nvert(ftype)+nedge(ftype)) then
              ie = jp-nvert(ftype)
              call edge_to_vert(ftype,ie, nv1,nv2)
              x_new(1:3,j) = (x(1:3,nv1) + x(1:3,nv2))/2.d0
!
#if HP3D_DEBUG
              if (iprint.ge.2) then
                write(*,*)'parent edge node: nv1,nv2 = ',nv1,nv2
                write(*,8000)j,x_new(1:3,j)
              endif
#endif
!
!  .........parent rectangular face node
            elseif (jp.le.nvert(ftype)+nedge(ftype)+nface(ftype)) then
              ifc = jp-nvert(ftype)-nedge(ftype)
              call face_to_vert(ftype,ifc, nv1,nv2,nv3,nv4)
              x_new(1:3,j) = (x(1:3,nv1) + x(1:3,nv2) &
                            + x(1:3,nv3) + x(1:3,nv4))/4.d0
!
#if HP3D_DEBUG
              if (iprint.ge.2) then
                write(*,*)'parent rectangular face node: ', &
                           'nv1,nv2,nv3,nv4 = ',            &
                            nv1,nv2,nv3,nv4
                write(*,8000)j,x_new(1:3,j)
              endif
#endif
!
!  .........parent middle node (brick element only)
            else
              x_new(1:3,j) = 0.d0
              do iv=1,8
                x_new(1:3,j) = x_new(1:3,j) + x(1:3,iv)
              enddo
              x_new(1:3,j) = x_new(1:3,j)/8.d0
!
#if HP3D_DEBUG
              if (iprint.ge.2) then
                write(*,*)'parent middle node (brick only!)'
                write(*,8000)j,x_new(1:3,j)
              endif
#endif
!
            endif
          endif
!
 8000   format(' j,x_new(1:3,j) = ',i2,2x,3(e12.5,2x))
!
!  .....end of loop through vertex nodes
        enddo
!
#if HP3D_DEBUG
!  .....printing
        if (iprint.ge.2) then
          write(*,7003)igen,Nvert(son_type(igen))
 7003     format(' igen,nvert_son = ',i2,2x,i2)
          do i=1,Nvert(son_type(igen))
            write(*,7004)i,x_new(1:3,i)
 7004       format(' i,x_new(1:3,i) = ',i2,2x,3(e12.5,2x))
          enddo
        endif
#endif
!
!  .....update
        x = x_new
!
!  ...end of loop through generations
      enddo
!
      Xsubel(1:3,1:nvert(NODES(Mdle)%ntype)) = &
           x(1:3,1:nvert(NODES(Mdle)%ntype))
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7001) Mdle,Xsubel(1:3,1:nvert(NODES(Mdle)%ntype))
 7001   format('refel: Mdle = ',i6,' Xsubel = ',8(/,3(f8.3,2x)))
        call pause
      endif
#endif
!
!
   end subroutine refel
