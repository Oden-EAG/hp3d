#if HP3D_USE_X11

!----------------------------------------------------------------------
!
!   routine name       - lsvistr
!
!----------------------------------------------------------------------
!
!   latest revision    - Mar 2023
!
!   purpose            - routine makes the list of all visible
!                        triangles and stores rescaled coordinates
!                        of their end-points (vertices) in an array
!                        RGTR
!
!   required  routines - trian,recta,sortz,clip2,trobs
!
!----------------------------------------------------------------------
!
   subroutine lsvistr
!
      use GMP
      use graphmod
!
      implicit none
!
      real(8) :: xmax(3),xmin(3)
      real(8) :: xi(2),xcoord(3,3),dxdxi(NDIM,2)
      real(8) :: pm(2),pn(2),xtro(3,3)
!
!  ...integer coordinates of points on a rectangular face
      integer :: ixi(2,3)
!
!  ...face rectangles numbers and their orientation for a block
      integer :: nsid_rect(6),nsid_orient(6),nsid_tria(4)
!
!  ...visibility flags for two sides of a figure
      integer :: nvis_flag(2)
!
!  ...binary flags for plotting boundary of a triangle
      integer :: nplotb(3)
!
      real(8) :: cl,sign,vis
      integer :: i,iclip,iflag,ii,iloc,inv2,is,isub,ivar,j,k,l,lab
      integer :: nb,nbb,nbl,ndom,nhalf,nick,nickb,nnn,nr,nsub1,nsub2,nt
!
#if HP3D_DEBUG
      integer :: iprint
#endif
!
      real(8) :: bigp,bign,small,one
      data bigp,bign,small,one /1.d30,-1.d30,1.d-6,1.d0/
!
!-----------------------------------------------------------------------
#if HP3D_DEBUG
      iprint=0
#endif
!
!  ...set default bounds for the picture
      do ivar=1,3
        XEX(2*ivar-1) = bigp
        XEX(2*ivar)   = bign
        xmin(ivar)    = bigp
        xmax(ivar)    = bign
      enddo
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7011) NRSUB
 7011   format('lsvistr: NRSUB = ', i2)
        write(*,7012) NDIM,MANDIM
 7012   format('lsvistr: NDIM,MANDIM = ', 2i2)
        write(*,7013) NRRECTA,NRTRIAN
 7013   format('lsvistr: NRRECTA,NRTRIAN = ', 2i8)
        call pause
      endif
#endif
!
!-----------------------------------------------------------------------
!
!  ...loop through GMP triangles
      NRVISTR = 0
      do nt=1,NRTRIAN
!
        if (MANDIM.eq.3) then
!
!  .......initiate number of neighbors
          nnn=0
!
!  .......initiate both sides visibility flags
          nvis_flag(1:2) = 0
!
!  .......loop through the adjacent blocks
          do is=1,2
!
!  .........the adjacent block's nickname  is
            nb = abs(TRIANGLES(nt)%BlockNo(is))
!
            if (nb.eq.0) then
!
!  ...........there is no block, this side is on the boundary,
!             reset the visibility flag
              nvis_flag(is)=1
            else
!
              nnn=nnn+1
!
!  ...........the block exists, check if not on the list
!             of invisible blocks
              call locate(nb,IGINV,NRINVBL, iloc)
!
!  ...........check if the block belongs to an invisible subdomain
              call decode(nb, nbb,lab)
              select case(lab)
              case(1)
                ndom = PRISMS(nbb)%domain
                if (NDOMAIN(ndom).eq.0) iloc=1
              case(3)
                ndom = TETRAS(nbb)%domain
                if (NDOMAIN(ndom).eq.0) iloc=1
              case(4)
                ndom = PYRAMIDS(nbb)%domain
                if (NDOMAIN(ndom).eq.0) iloc=1
              case default
                write(*,*) 'lsvistr: nb = ',nb
                stop 1
              end select
!
              if (iloc.ne.0) then
!
!  .............the block is invisible, this side should be visible
!               then, reset the visibility flag
                nvis_flag(is)=1
              endif
            endif
!
!  .......end of loop through sides of the triangle
          enddo
!
#if HP3D_DEBUG
          if (iprint.ge.1) then
            write(*,7003) nt,nvis_flag
          endif
#endif
!
!  .......the triangle will be visible if only one of the flags
!         is raised
          select case(nvis_flag(1)*10+nvis_flag(2))
          case(00)
!
!  .........skip the triangle, there is no visible block on
!           either side
            goto 60
          case(11)
!
!  .........skip the triangle, it is in between two visible
!           blocks (an interior face) and, therefore, it is
!           invisible
            goto 60
          case(01)
!
!  .........the block on the first side is visible, on the second
!           side there is no visible block, display the triangle
!           as a face of the first block
            nb = TRIANGLES(nt)%BlockNo(1)
          case(10)
!
!  .........the block on the second side is visible, on the first
!           side there is no visible block, display the triangle
!           as a face of the second block
            nb = TRIANGLES(nt)%BlockNo(2)
          end select
!
!  .......check if on the list of curvilinear blocks
          call locate(nb,NLINBLOCKS,NRCURVBL, iloc)
!
#if HP3D_DEBUG
          if (iprint.ge.1) then
            write(*,7010) nt,nb,iloc
 7010       format('lsvistr: DISPLAYING TRIANGLE nt,nb,iloc = ',3i8)
            call pause
          endif
#endif
!
!  .......find the orientation of the triangle wrt the adjacent block
          call decode(abs(nb), nbl,lab)
          select case(lab)
!
!  .......prism
          case(1)
            do is=1,2
              call decode(PRISMS(nbl)%FigNo(is), &
                          nsid_tria(is),nsid_orient(is))
            enddo
!
#if HP3D_DEBUG
            if (iprint.eq.2) then
              write(*,7067) nbl,nsid_tria(1:2),nsid_orient(1:2)
 7067         format('lsvistr: PRISM = ',i6,' SIDES WITH ORIENT = ', &
                      2i6,2x,2i2)
            endif
#endif
!
            call locate(nt,nsid_tria(1:2),2, ii)
            if (ii.eq.0) then
              write(*,*) 'lsvistr: INCONSISTENCY 1 !'
              write(*,*) 'TRIANGLE, PRISM = ',nt,nbl
              call print_GMP
              stop 1
            endif
            if (nsid_orient(ii).le.2) then
              sign =   one
            else
              sign = - one
            endif
!
!  .......tetrahedron
          case(3)
            do is=1,4
              call decode(TETRAS(nbl)%FigNo(is), &
                          nsid_tria(is),nsid_orient(is))
            enddo
!
!  .........locate the triangle
            call locate(nt,nsid_tria,4, ii)
            if (ii.eq.0) then
              write(*,*) 'lsvistr: INCONSISTENCY 2 !'
              stop 1
            endif
            if (nsid_orient(ii).le.2) then
              sign =   one
            else
              sign = - one
            endif
!
!  .......pyramid
          case(4)
            do is=1,4
              call decode(PYRAMIDS(nbl)%FigNo(1+is), &
                          nsid_tria(is),nsid_orient(is))
            enddo
            call locate(nt,nsid_tria(1:4),4, ii)
            if (ii.eq.0) then
              write(*,*) 'lsvistr: INCONSISTENCY 3 !'
              stop 1
            endif
            if (nsid_orient(ii).le.2) then
              sign =   one
            else
              sign = - one
            endif
!
          case default
            write(*,*) 'lsvistr: lab = ',lab
            stop 1
          end select
!
#if HP3D_DEBUG
          if (iprint.ge.1) then
            write(*,7002) nbl,ii,nsid_orient(ii),sign
          endif
#endif
!
!  .....the case of a 2D manifold
        elseif (MANDIM.eq.2) then
!
!  .......loop through the list of invisible triangles and check
!         whether the current rectangle is on the list
          do inv2 = 1,NRINVBL
!
!  .........skip if on the list
            if ((nt*10+1).eq.IGINV(inv2)) goto 60
          enddo
          sign = one
!
        else
          write(*,*) 'lsvistr: MANDIM = ',MANDIM
          stop 1
        endif
!
!
!  .....loop through subtriangles...
        do j=1,NRSUB
          do i=1,NRSUB-j+1
            do l=1,2
!
#if HP3D_DEBUG
              if (iprint.eq.2) then
                write(*,7039) j,i,l
 7039           format('lsvistr: LOOPING THROUGH SUBTRIANGLES j,i,l = ', &
                        3i3)
              endif
#endif
!
!  ...........find coordinates of vertices and transform them
              select case(l)
!
!  ..........."lower" triangle
              case(1)
              xi(1) = i*DX
              xi(2) = (j-1)*DX
              call trian(nt,xi, xcoord(1,1),dxdxi)
              call trobs(xcoord(1,1), xtro(1,1))
!
              xi(1) = (i-1)*DX
              xi(2) = j*DX
              call trian(nt,xi, xcoord(1,2),dxdxi)
              call trobs(xcoord(1,2), xtro(1,2))
!
              xi(1) = (i-1)*DX
              xi(2) = (j-1)*DX
              call trian(nt,xi, xcoord(1,3),dxdxi)
              call trobs(xcoord(1,3), xtro(1,3))
!
!  ..........."upper" triangle
              case(2)
              if (j.eq.NRSUB-i+1) then
                goto 30
              else
                xi(1) = i*DX
                xi(2) = (j-1)*DX
                call trian(nt,xi, xcoord(1,1),dxdxi)
                call trobs(xcoord(1,1), xtro(1,1))
!
                xi(1) = i*DX
                xi(2) = j*DX
                call trian(nt,xi, xcoord(1,2),dxdxi)
                call trobs(xcoord(1,2), xtro(1,2))
!
                xi(1) = (i-1)*DX
                xi(2) = j*DX
                call trian(nt,xi, xcoord(1,3),dxdxi)
                call trobs(xcoord(1,3), xtro(1,3))
              endif
              end select
!
!  ...........update bounds for the object
              do k=1,3
                do ivar=1,3
                  XEX(2*ivar-1)=min(XEX(2*ivar-1),xcoord(ivar,k))
                  XEX(2*ivar)=max(XEX(2*ivar),xcoord(ivar,k))
                enddo
              enddo
!
!  ...........determine whether the triangle is to be displayed or not
!        .....first perform section plane clipping
              do k = 1,3
                cl = CLPL(1)*xcoord(1,k) + CLPL(2)*xcoord(2,k) &
                   + CLPL(3)*xcoord(3,k) + CLPL(4)
                if (cl.gt.0.d0) goto 30
              enddo
!
!      .......next 2D clipping wrt the current window
              call clip2(xtro,iclip)
              if (iclip.eq.0) goto 30
!
!  ...........and finally check the orientation in space
              do ivar=1,2
                pm(ivar) = xtro(ivar,2) - xtro(ivar,1)
                pn(ivar) = xtro(ivar,3) - xtro(ivar,1)
              enddo
              vis = pm(1)*pn(2)-pm(2)*pn(1)
              vis = vis*sign
!!!              if ((vis.lt.small).and.(MANDIM.eq.3)) goto 30
!
!
!  ...........update list of visible facets and bounds for the picture
              NRVISTR = NRVISTR+1
              if (NRVISTR.gt.MXIGTR) then
                write(*,*) 'lsvistr: INCREASE MXIGTR'
                stop 1
              endif
              do k=1,3
                do ivar=1,3
                  isub = (NRVISTR-1)*9+(k-1)*3+ivar
                  RGTR(isub) = xtro(ivar,k)
                  xmin(ivar) = min(xmin(ivar),xtro(ivar,k))
                  xmax(ivar) = max(xmax(ivar),xtro(ivar,k))
                enddo
              enddo
!
!  ...........update list of triangles' midpoints
              RGTRZ(NRVISTR) = (xtro(3,1)+xtro(3,2)+xtro(3,3))/3
!
              IGTRCU(NRVISTR) = 0
!
!
!  ...........encode block's number
              select case(MANDIM)
              case(2)
                nick=10*nt+1
              case(3)
                nick = nbl
              end select
              if ((i.eq.NRSUB/2).and.(j.eq.NRSUB/2).and. &
                  (l.eq.1)) then
                IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+10000*nick
              endif
!
!  ...........encode the color
              select case(nnn)
              case(0)
                write(*,*) 'ERROR !!'; stop 1
              case(1)
                IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+1000
              case(2)
                IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+2000
              end select
!
!  ...........update list of visible curves when on boundary
              if (l.eq.1) then
                nick=0; nplotb(1:3)=0
                if (i.eq.1) then
                  nplotb(2)=1
                  if (j.eq.(NRSUB)/2) nick=abs(TRIANGLES(nt)%EdgeNo(3))
                endif
                if (j.eq.1) then
                  nplotb(3)=1
                  if (i.eq.(NRSUB)/2) nick=abs(TRIANGLES(nt)%EdgeNo(1))
                endif
                if (j.eq.NRSUB-i+1) then
                  nplotb(1)=1
                  if (i.eq.(NRSUB)/2) nick=abs(TRIANGLES(nt)%EdgeNo(2))
                endif
!======================================================================
! Activate to make parameterization visible
              nplotb(1:3)=1
!======================================================================
                call encodg(nplotb,2,3, nickb)
!
!ldem, encoding of curve numbers disabled
!!!                IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+10000*nick
                IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+nickb
              endif
   30       enddo
          enddo
        enddo
   60 enddo
!
!-----------------------------------------------------------------------
!
!  ...loop through GMP rectangles
      do nr=1,NRRECTA
!
        if (MANDIM.eq.3) then
!
!  .......initiate number of neighbors
          nnn=0
!
!  .......initiate both sides visibility flags
          nvis_flag(1:2) = 0
!
!  .......loop through the adjacent blocks
          do is=1,2
!
!  .........the adjacent block is
            nb = RECTANGLES(nr)%BlockNo(is)
!
#if HP3D_DEBUG
            if (iprint.ge.1) then
              write(*,7035) nr,is,nb
 7035         format('lsvistr: nr,is,nb = ',i6,i3,i7)
            endif
#endif
!
            if (nb.eq.0) then
!
!  ...........there is no block, this side is on the boundary,
!             reset the visibility flag
              nvis_flag(is)=1
            else
              nnn=nnn+1
!
!  ...........the block exists, check if not on the list
!             of invisible blocks
              call locate(nb,IGINV,NRINVBL, iloc)
!
#if HP3D_DEBUG
              if (iprint.ge.1) then
                write(*,7029) nr,is,iloc
              endif
#endif
!
!  ...........check if the block belongs to an invisible subdomain
              call decode(nb, nbb,lab)
              select case(lab)
              case(1)
                ndom = PRISMS(nbb)%domain
                if (NDOMAIN(ndom).eq.0) iloc=1
#if HP3D_DEBUG
                if (iprint.ge.1) write(*,*) 'is = ',is
#endif
              case(2)
                ndom = HEXAS(nbb)%domain
                if (NDOMAIN(ndom).eq.0) iloc=1
              case(4)
                ndom = PYRAMIDS(nbb)%domain
                if (NDOMAIN(ndom).eq.0) iloc=1
              case default
                write(*,*) 'lsvistr: nb = ',nb
                stop 1
              end select
!
#if HP3D_DEBUG
              if (iprint.ge.1) then
                write(*,7029) nr,is,iloc
 7029           format('lsvistr: nr,is,iloc = ',i6,i2,i4)
              endif
#endif
!
              if (iloc.ne.0) then
!
!  .............the block is invisible, this side should be visible
!               then, reset the visibility flag
                nvis_flag(is)=1
              endif
            endif
!
#if HP3D_DEBUG
            if (iprint.ge.1) write(*,*) 'is = ',is
#endif
!  .......end of loop through sides
          enddo
!
#if HP3D_DEBUG
          if (iprint.ge.1) then
            write(*,7003) nr,nvis_flag
 7003       format('lsvistr: nr,nvis_flag = ',i8,2x,2i2)
          endif
#endif
!
!  .......the rectangle will be visible if only one of the flags
!         is raised
          select case(nvis_flag(1)*10+nvis_flag(2))
          case(00)
!
!  .........skip the rectangle, there is no visible block on
!           either side
            goto 160
          case(11)
!
!  .........skip the rectangle, it is in between two visible
!           blocks (an interior face) and, therefore, it is
!           invisible
            goto 160
          case(01)
!
!  .........the block on the first side is visible, on the second
!           side there is no visible block, display the rectangle
!           as a face of the first block
            nb = abs(RECTANGLES(nr)%BlockNo(1))
          case(10)
!
!  .........the block on the second side is visible, on the first
!           side there is no visible block, display the rectangle
!           as a face of the second block
            nb = abs(RECTANGLES(nr)%BlockNo(2))
          end select
!
!  .......check if on the list of curvilinear blocks
          call locate(nb,NLINBLOCKS,NRCURVBL, iloc)
!
#if HP3D_DEBUG
          if (iprint.ge.1) then
            write(*,7001) nr,nb,iloc
 7001       format('lsvistr: DISPLAYING RECTANGLE nr,nb,iloc = ',2i6,i2)
            call pause
          endif
#endif
!
!  .......find the orientation of the rectangle wrt the adjacent block
          call decode(abs(nb), nbl,lab)
          select case(lab)
!
!  .......prism
          case(1)
            do is=1,3
              call decode(PRISMS(nbl)%FigNo(2+is), &
                          nsid_rect(is),nsid_orient(is))
            enddo
            call locate(nr,nsid_rect,3, ii)
            if (ii.eq.0) then
              write(*,*) 'lsvistr: INCONSISTENCY 4 !'
              stop 1
            endif
            if (nsid_orient(ii).le.3) then
              sign =   one
            else
              sign = - one
            endif
            is=2+ii
!
!  .......hexahedron
          case(2)
            do is=1,6
              call decode(HEXAS(nbl)%FigNo(is), &
                          nsid_rect(is),nsid_orient(is))
            enddo
!
!  .........locate the rectangle
            call locate(nr,nsid_rect,6, is)
            select case(is)
              case(2,3,4)
                if (nsid_orient(is).le.3) then
                  sign =   one
                else
                  sign = - one
                endif
              case(1,5,6)
                if (nsid_orient(is).le.3) then
                  sign = - one
                else
                  sign =   one
                endif
              case default
                write(*,*) 'lsvistr: ERROR '
                stop 1
            end select
!
!  .......pyramid
          case(4)
            do is=1,1
              call decode(PYRAMIDS(nbl)%FigNo(is), &
                          nsid_rect(is),nsid_orient(is))
            enddo
            call locate(nr,nsid_rect(1),1, ii)
            if (ii.eq.0) then
              write(*,*) 'lsvistr: INCONSISTENCY 5 !'
              stop 1
            endif
            if (nsid_orient(ii).le.3) then
              sign =   one
            else
              sign = - one
            endif
!
          case default
            write(*,*) 'lsvistr: lab = ',lab
            stop 1
          end select
!
#if HP3D_DEBUG
          if (iprint.ge.1) then
            write(*,7002) nbl,is,nsid_orient(is),sign
 7002       format('lsvistr: nbl,is,nsid_orient(is),sign = ',3i4,f7.2)
          endif
#endif
!
!  .....the case of a 2D manifold
        elseif (MANDIM.eq.2) then
!
!  .......loop through the list of invisible rectangles and check
!         whether the current rectangle is on the list
          do inv2 = 1,NRINVBL
!
!  .........skip if on the list
            if ((nr*10+abs(nr)).eq.IGINV(inv2)) goto 160
          enddo
          sign = one
!
        else
          write(*,*) 'lsvistr: MANDIM = ',MANDIM
          stop 1
        endif
!
!        DX = 1
        nhalf=NRSUB/2; nsub1 = NRSUB; nsub2 = NRSUB
!
!  .....loop through subtriangles...
        do i=1,NRSUB
          do j=1,NRSUB
            do l=1,2
!
!  ...........set integer coordinates for the triangle vertices
              if ((i.gt.nhalf  .and. j.le.nhalf) &
                .or. (i.le.nhalf .and. j.gt.nhalf)) then
!
!  .............flag for drawing the boundary
                iflag = 1
                select case(l)
                case(1)
                  ixi(1,1)=i;   ixi(1,2)=i-1; ixi(1,3)=i-1
                  ixi(2,1)=j-1; ixi(2,2)=j;   ixi(2,3)=j-1
                case(2)
                  ixi(1,1)=i;   ixi(1,2)=i;   ixi(1,3)=i-1
                  ixi(2,1)=j-1; ixi(2,2)=j;   ixi(2,3)=j
                end select
              else
!
!  .............flag for drawing the boundary
                iflag = 2
                select case(l)
                case(1)
                  ixi(1,1)=i-1; ixi(1,2)=i;   ixi(1,3)=i
                  ixi(2,1)=j-1; ixi(2,2)=j-1; ixi(2,3)=j
                case(2)
                  ixi(1,1)=i-1; ixi(1,2)=i;   ixi(1,3)=i-1
                  ixi(2,1)=j-1; ixi(2,2)=j;   ixi(2,3)=j
                end select
              endif
!
!  ...........find coordinates of vertices and transform them
              do k=1,3
                xi(1) = ixi(1,k)*DX
                xi(2) = ixi(2,k)*DX
                if (iloc.gt.0) then
                  call recta(nr,xi, xcoord(1,k),dxdxi)
                else
!!!                  call recta_linear(nr,xi, xcoord(1,k),dxdxi)
                  call recta(nr,xi, xcoord(1,k),dxdxi)
                endif
                call trobs(xcoord(1,k), xtro(1,k))
              enddo
!
!  ...........update bounds for the object
              do k=1,3
                do ivar=1,3
                  XEX(2*ivar-1) = min(XEX(2*ivar-1),xcoord(ivar,k))
                  XEX(2*ivar)   = max(XEX(2*ivar),  xcoord(ivar,k))
                enddo
              enddo
!
!  ...........determine whether the triangle is to be displayed or not
!          ...first perform section plane clipping
              do k = 1,3
                cl = CLPL(1)*xcoord(1,k) + CLPL(2)*xcoord(2,k) &
                   + CLPL(3)*xcoord(3,k) + CLPL(4)
                if (cl.gt.0.d0) cycle
              enddo
!
!      .......next 2D clipping wrt the current window
              call clip2(xtro,iclip)
              if (iclip.eq.0) cycle
!
!      .......and finally check the orientation in space
              do ivar=1,2
                pm(ivar) = xtro(ivar,2) - xtro(ivar,1)
                pn(ivar) = xtro(ivar,3) - xtro(ivar,1)
              enddo
              vis = pm(1)*pn(2)-pm(2)*pn(1)
!
!  ...........update list of visible facets and bounds for the picture
              NRVISTR = NRVISTR+1
              if (NRVISTR.gt.MXIGTR) then
                write(*,*) 'lsvistr: INCREASE MXIGTR = ',MXIGTR
                stop 1
              endif
              do k=1,3
                do ivar=1,3
                  isub = (NRVISTR-1)*9+(k-1)*3+ivar
                  RGTR(isub) = xtro(ivar,k)
                  xmin(ivar) = min(xmin(ivar),xtro(ivar,k))
                  xmax(ivar) = max(xmax(ivar),xtro(ivar,k))
                enddo
              enddo
!
!  ...........update list of triangles' midpoints
              RGTRZ(NRVISTR) = (xtro(3,1)+xtro(3,2)+xtro(3,3))/3
!
              IGTRCU(NRVISTR) = 0
!
!  ...........encode block's number
              select case(MANDIM)
              case(2)
                nick=10*nr+2
              case(3)
                nick = nbl
              end select
              if ((i.eq.(1+NRSUB/2)).and.(j.eq.(1+NRSUB/2)).and. &
                  (l.eq.1)) then
                IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+10000*nick
              endif
!
!  ...........encode the color
              select case(nnn)
              case(0)
                write(*,*) 'ERROR !!'; stop 1
              case(1)
                IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+1000
              case(2)
                IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+2000
              end select
!!              select case(MANDIM)
!!              case(2)
!!                if (vis.lt.small) &
!!                  IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+1000
!!              case(3)
!!                if (iloc.eq.0) then
!!                  if (sign.gt.0.d0) then
!!                    IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+1000
!!                  else
!!                    IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+2000
!!                  endif
!!                elseif (iloc.gt.0) then
!!                  if (sign.gt.0.d0) then
!!                    IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+3000
!!                  else
!!                    IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+4000
!!                  endif
!!                endif
!!              end select
!
!  ...........update list of visible curves when on boundary
              nick=0; nplotb(1:3)=0
              select case(iflag)
              case(1)
                select case(l)
                case(1)
                  if (j.eq.1) then
                    nplotb(3)=1
!!!                    IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+3
                  endif
                  if (i.eq.1) then
                    nplotb(2)=1
!!!                    IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+2
                  endif
                case(2)
                  if (i.eq.nsub1) then
                    nplotb(1)=1
!!!                    IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+1
                  endif
                  if (j.eq.nsub2) then
                    nplotb(2)=1
!!!                    IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+2
                  endif
                end select
              case(2)
                select case(l)
                case(1)
                  if (j.eq.1) then
                    nplotb(1)=1
!!!                    IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+1
                  endif
                  if (i.eq.nsub1) then
                    nplotb(2)=1
!!!                    IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+2
                  endif
                case(2)
                  if (i.eq.1) then
                    nplotb(3)=1
!!!                    IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+3
                  endif
                  if (j.eq.nsub2) then
                    nplotb(2)=1
!!!                    IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+2
                  endif
                end select
              end select
              nick=0
              if (l.eq.1) then
                if (j.eq.1) then
                 if (i.eq.(NRSUB)/2) nick=abs(RECTANGLES(nr)%EdgeNo(1))
                endif
                if (i.eq.NRSUB) then
                 if (j.eq.(NRSUB)/2) nick=abs(RECTANGLES(nr)%EdgeNo(2))
                endif
              else
                if (i.eq.1) then
                 if (j.eq.(NRSUB)/2) nick=abs(RECTANGLES(nr)%EdgeNo(4))
                endif
                if (j.eq.NRSUB) then
                 if (i.eq.(NRSUB)/2) nick=abs(RECTANGLES(nr)%EdgeNo(3))
                endif
              endif
              IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+10000*nick
!================================================================
!  Activate to make the parameterization visible
              nplotb(1:3)=1
!================================================================
              call encodg(nplotb,2,3, nickb)
              IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+nickb


!
!  .........end of loop through l (lower/upper triangle)
            enddo
!
!  .......end of loop through j
          enddo
!
!  .....end of loop through i
        enddo
  160 enddo
!
!
!  ...compute object's dimensions and coordinates of its
!     central point
      DIMOB(1) = (xmax(1)-xmin(1))/2.d0
      DIMOB(2) = (xmax(2)-xmin(2))/2.d0
      DIMOB(3) = (xmax(3)-xmin(3))/2.d0
      XCENTR(1) = (xmax(1)+xmin(1))/2.d0
      XCENTR(2) = (xmax(2)+xmin(2))/2.d0
      XCENTR(3) = (xmax(3)+xmin(3))/2.d0
!
!  ...sort a!!ording to z-coordinate in order back-to-front
      call sortz
!
   end subroutine lsvistr

#endif
