#if HP3D_USE_X11

!----------------------------------------------------------------------
!
!   routine name       - display_face
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine makes a list of all visible
!                        triangles and stores rescaled coordinates
!                        of their end-points (vertices) in an array
!                        RGTR for an element face
!
!   arguments :
!     in:
!         Numlev       - number of levels to plot solution values
!                        (0-encode p-approximation colors)
!         Mdle         - element number, same as the middle node
!         Ifc          - face number
!         Nodesl       - element nodes
!         Nedge_orient - edge orientation
!         Nface_orient - face orientation
!         Norder       - order of approximation for the element
!         Iflagn       = 1 encode element (middle node) numbers
!                      = 2 encode face nodes numbers
!                      = 3 color code BC's flags
!         Xnod         - geometry dof for the element
!         ZdofH        - H1 dof for the element
!         ZdofE        - H(curl) dof for the element
!         ZdofV        - H(div) dof for the element
!         ZdofQ        - L2 dof for the element
!         Solev        - scale of solution levels
!         Xmin,Xmax    - bounds for the picture
!
!----------------------------------------------------------------------
#include "typedefs.h"
   subroutine display_face(Numlev,Mdle,Ifc,Nodesl,                   &
                              Nedge_orient,Nface_orient,Norder,Iflagn,  &
                              Xnod,ZdofH,ZdofE,ZdofV,ZdofQ,             &
                              Solev,Xmin,Xmax,                          &
                              Nrintertot,Nrinter,Xlocinter)
!
      use element_data , only : NFAXES
      use data_structure3D
      use graphmod
      use physics
!
      implicit none
!
      integer :: Numlev,Mdle,Ifc,Iflagn,Nrintertot,Nrinter
      integer :: Nedge_orient(12),Nface_orient(6),Norder(19),Nodesl(27)
      real(8) :: Xnod(3,MAXbrickH),Solev(*),Xmin(3),Xmax(3)
      VTYPE   :: ZdofH(MAXEQNH,MAXbrickH),ZdofE(MAXEQNE,MAXbrickE),  &
                 ZdofV(MAXEQNV,MAXbrickV),ZdofQ(MAXEQNQ,MAXbrickQ)
      real(8) :: Xlocinter(3,500)
!
!  ...local face nodes numbers
      integer :: nface_nodes(9)
!
!  ...global face nodes numbers
      integer :: nface_nodes_global(9)
!
!  ...local order of approximation (4 edges + 1 middle)
      integer :: norder_face(5)
!
!  ...master face and element coordinates
      integer :: it(2,3,2),id_flag(2)
      real(8) :: t(2),xi(3,3),pm(2),pn(2)
      real(8) :: xtro(3,3),xcoord(3,3),val(3),costr(2,3,3)
      real(8) :: utemp(3),cotriancut(2,3,3)
      
!
!  ...color code for subtriangles, nicknames encoding which edges
!     of a subtriangle should be drawn
      integer :: mpcol(2),nsid(2)
!
!  ...BC's flags
      integer :: ibc(6,NRINDEX)
!
!  ...misc
      real(8) :: dt,sdw,sup,vis
      integer :: i,iaux,iclip,idec,iflag,ifree,il,ilev
      integer :: isav,isflag,istr,isub,itriancut,ivar,j,k,l
      integer :: ntype,nfstr,nhalf,nick,nod,nordh,nordv,nrfn
      integer :: nrtriancut,nstr,nstrl,nsub1,nsub2
!
      real(8) :: small = 1.d-6
!
#if HP3D_DEBUG
      integer :: iprint,nrdofH,nrdofE,nrdofV,nrdofQ
#endif
!
!---------------------------------------------------------------------
!
      ntype = NODES(Mdle)%ntype
!
#if HP3D_DEBUG
      iprint = 0
      if (iprint.eq.1) then
        write(*,*) 'display_face: type,Mdle,Ifc = ', &
                                  S_Type(ntype),Mdle,Ifc
        write(*,*) '              face_type(ntype,Ifc) = ', &
                                  face_type(ntype,Ifc)
        write(*,*) '              Xnod = '
        call celndof(ntype,Norder, nrdofH,nrdofE,nrdofV,nrdofQ)
        do k=1,nrdofH
          write(*,7003) k,Xnod(1:3,k)
 7003     format('k = ',i4,3e12.5)
        enddo
        call pause
      endif
#endif
!
!  ...get the element boundary conditions flags
      call find_bc(Mdle, ibc)
!
!  ...determine the local number of the nodes on the face
!     [vertices ; edges ; middle ]
      call face_nodes(ntype,Ifc, nface_nodes,nrfn)
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7010) Ifc, nface_nodes(1:nrfn)
 7010   format('display_face: FACE = ',i2,' FACE NODES = ',9i3)
        call pause
      endif
#endif
!
!  ...loop over face nodes and record:
!      - edge nodes and mid-face node
!      - edge nodes orders of approximation
!      - mid-face node order of approximation a!!ording to face
!        LOCAL system of coordinates
      j=0
      do i=1,nrfn
        k=nface_nodes(i)-nvert(ntype)
!
!  .....record global node number
        l=nface_nodes(i)
        nface_nodes_global(i)=Nodesl(l)
!
!  .....if node is not a vertex
        if (k.gt.0) then
          j=j+1
!
!  .......record order of approximation a!!ounting for orientation
!!!!!          nod=nface_nodes(j)
          nod = nface_nodes_global(i)
          select case(NODES(nod)%ntype)
!
!  .......rectangular face node
          case(MDLQ,RECT)
            select case(NFAXES(3,Nface_orient(Ifc)))
!  .........the face axes have NOT been reversed
            case(0) ; norder_face(j)=Norder(k)
!
!  .........the face axes HAVE been reversed
            case(1)
              call decode(Norder(k), nordh,nordv)
              norder_face(j)=nordv*10+nordh
            endselect
!
!  .......all other nodes
          case default
            norder_face(j)=Norder(k)
          endselect
        endif
      enddo
!
      Nrinter=0
!
!  ...set number of subdivisions
      nsub1=NRSUB
      nhalf=NRSUB/2
      dt=1.d0/NRSUB
!
!  ...loop through subtriangles...
      do i=1,nsub1
!
        select case(face_type(ntype,Ifc))
          case(TRIA); nsub2=NRSUB+1-i
          case(RECT); nsub2=NRSUB
        end select
!
        do j=1,nsub2
!
!  .......determine integer coordinates for the subtriangles
!         and color flags for the mesh visualization

          call int_coord(Numlev,NRSUB,face_type(ntype,Ifc),norder_face, &
                         i,j,it(1:2,1:3,1),it(1:2,1:3,2),               &
                         mpcol(1),mpcol(2),nsid(1),nsid(2))
          if (Iflagn.eq.3) then
!...TODO: this is not correct after updating BC implementation: NRPHY_DISP
!         now has to refer to a particular component (not physics attr)
            mpcol(1:2) = ibc(Ifc,NRPHY_DISP)+1
          endif
#if HP3D_DEBUG
          if (iprint.eq.1) then
            write(*,7002) i,j,mpcol(1),mpcol(2)
 7002       format('display_face: i,j,mpcol = ',4i3)
          endif
#endif
!
!  .......loop through the lower/upper triangle
          do l=1,2
!
!  .........check if the triangle is to be displayed
            if (mpcol(l).eq.0) goto 20
!
!  .........find coordinates of vertices
            do k=1,3
              t(1:2) = it(1:2,k,l)*dt
              call compute_face(Numlev,Mdle,Ifc,                  &
                                Nedge_orient,Nface_orient,Norder, &
                                Xnod,ZdofH,ZdofE,ZdofV,ZdofQ,t,   &
                                xcoord(1:3,k),val(k))
!
!  ...........transform the coordinates to observer's system
              call trobs(xcoord(1,k), xtro(1,k))
!
!  .........end of loop through the vertices
            enddo
#if HP3D_DEBUG
            if (iprint.eq.1) then
              write(*,8002) i,j,l
 8002         format('display_face: i,j,l = ',3i3)
              write(*,8001) (xi(1:3,k),k=1,3)
 8001         format('display_face: xi = ',3(3f7.2,2x))
              write(*,8003) (xcoord(1:3,k),k=1,3)
 8003         format('display_face: xcoord = ',3(3f7.2,2x))
              write(*,8004) (xtro(1:3,k),k=1,3)
 8004         format('display_face: xtro = ',3(3f7.2,2x))
              write(*,8006) val(1:3)
 8006         format('display_face: val = ',3e12.5)
            endif
#endif
!
!  .........perform 2D clipping wrt the current window
            call clip2(xtro,iclip)
!
!  .........and check the orientation in space
            do ivar=1,2
              pm(ivar) = xtro(ivar,2) - xtro(ivar,1)
              pn(ivar) = xtro(ivar,3) - xtro(ivar,1)
            enddo
            vis = pm(1)*pn(2)-pm(2)*pn(1)
!
!  .........update bounds for the object
            do k=1,3
              do ivar=1,3
                XEX(2*ivar-1)=min(XEX(2*ivar-1),xcoord(ivar,k))
                XEX(2*ivar)=max(XEX(2*ivar),xcoord(ivar,k))
              enddo
            enddo
!
!  .........create one or two visible triangles depending upon
!           intersection with the clipping plane
            utemp(1)=CLPL(1)*xcoord(1,1)+CLPL(2)*xcoord(2,1)+  &
                     CLPL(3)*xcoord(3,1)+CLPL(4)
            utemp(2)=CLPL(1)*xcoord(1,2)+CLPL(2)*xcoord(2,2)+  &
                     CLPL(3)*xcoord(3,2)+CLPL(4)
            utemp(3)=CLPL(1)*xcoord(1,3)+CLPL(2)*xcoord(2,3)+  &
                     CLPL(3)*xcoord(3,3)+CLPL(4)
!
            call fincut(xi,utemp, nrtriancut,cotriancut,isflag)
#if HP3D_DEBUG
            if (iprint.eq.1) then
              write(*,*) 'display_face: isflag = ',isflag
            endif
#endif
!
!  .........loop over visible cut triangles
            do itriancut=1,nrtriancut
!
!  ...........if we have  new triangles
              if (isflag.eq.1) then
!
!  .............collect new vertices
                Nrinter=Nrinter+1
                Nrintertot=Nrintertot+1
                do il=1,3
                  Xlocinter(il,Nrintertot)=cotriancut(1,il,1)
                enddo
                Nrinter=Nrinter+1
                Nrintertot=Nrintertot+1
                do il=1,3
                  Xlocinter(il,Nrintertot)=cotriancut(1,il,2)
                enddo
!
!  .............skip not visible small triangles
                if (iclip.eq.0) goto 20
                select case (Ifc)
!
!  ...............orientation of the face consistent with the external
!                 normal
                  case(2,3,4)
                    if (vis.lt.small) goto 20
!
!  ...............orientation of the face opposite to the external
!                 normal
                  case(1,5,6)
                    if (-vis.lt.small) goto 20
                end select
!
!
                iflag=0
!
!  .............substitute new local coordinates
                do k=1,3
                  do il=1,3
                    xi(il,k) = cotriancut(itriancut,il,k)
                  enddo
                  call compute_face(Numlev,Mdle,Ifc,                    &
                                    Nedge_orient,Nface_orient,Norder,   &
                                    Xnod,ZdofH,ZdofE,ZdofV,ZdofQ,t,     &
                                    xcoord(1:3,k),val(k))
!
!  ...............transform the coordinates to observer's system
                  call trobs(xcoord(1,k), xtro(1,k))
!
                enddo
!
              endif
!
!  ...........skip not visible small triangles
#if HP3D_DEBUG
              if (iprint.eq.1) then
                write(*,*) 'display_face: iclip,vis = ',iclip,vis
              endif
#endif
              if (iclip.eq.0) goto 20
!!!              if (vis.lt.small) goto 20
!
!  ...........store display information
              NRVISTR = NRVISTR+1
              if (NRVISTR.gt.MXIGTR) then
                write(*,9997)MXIGTR
 9997           format(' Increase MXIGTR = ',i15)
                write(*,*)'Dump out? 1-Yes/0-No'
                read(*,*)idec
                if (idec.ne.0) then
                  call dumpout_physics
                  call dumpout_hp3d('files/dumpc3Dhp')
                endif
                stop
              endif
              IGTRCU(NRVISTR) = 0
              IGTRNO(NRVISTR) = 0
!
!  ...........for p-approximation
              if (Numlev.eq.0) then
                do  k=1,3
                  do  ivar=1,3
                    isub = (NRVISTR-1)*9+(k-1)*3+ivar
                    if (isub.gt.MXRGTRZ) then
                      write(*,9996)MXRGTRZ
 9996                 format(' Increase MXRGTRZ = ',i15)
                      write(*,*)'Dump out? 1-Yes/0-No'
                      read(*,*)idec
                      if (idec.ne.0) then
                        call dumpout_physics
                        call dumpout_hp3d('files/dumpc3Dhp')
                      endif
                      stop
                    endif
                    RGTR(isub) = xtro(ivar,k)
                    Xmin(ivar) = min(Xmin(ivar),xtro(ivar,k))
                    Xmax(ivar) = max(Xmax(ivar),xtro(ivar,k))
                  enddo
                enddo
!!!                write(*,8020) Xmin,Xmax
 8020           format('Xmin,Xmax = ',2(3f8.3,2x))
!
!  ...........for contour plot
              else
!
!  .............find entries in matrix of small triangles
                if (NRVISTR.eq.1) then
                  ifree=0
                else
                  nfstr=IGTRNO(NRVISTR-1)/100
                  nstr=IGTRNO(NRVISTR-1)-100*nfstr
                  ifree=nfstr+nstr
                endif
                nstr=0
                do iaux=1,Numlev+2
                  ilev=iaux-1
                  if (ilev.eq.0) then
                    sup=Solev(1)
                    sdw=-1.d10
                  else if(ilev.le.Numlev) then
                    sup=Solev(ilev+1)
                    sdw=Solev(ilev)
                  else
                    sup=1.d10
                    sdw=Solev(Numlev+1)
                  endif
                  call finstr(xtro,val,sdw,sup, nstrl,costr,id_flag)
!
!  ...............for each small triangle
                  do istr=1,nstrl
                    if (ifree+nstr+istr.gt.MXIGSTR) then
                      write(*,9999)MXIGSTR
 9999                 format(' Increase MXIGSTR = ',i15)
                      write(*,*)'Dump out? 1-Yes/0-No'
                      read(*,*)idec
                      if (idec.ne.0) then
                        call dumpout_physics
                        call dumpout_hp3d('files/dumpc3Dhp')
                      endif
                      stop
                    endif
!
!  .................encode the color and coordinates
                    if (IDISPLAY_TYPE.eq.1) then
                      IGSTR(ifree+nstr+istr)=id_flag(istr)
                    else
                      IGSTR(ifree+nstr+istr)=ilev
                    endif
                    do k=1,3
                      do ivar=1,3
                        isub = (ifree+nstr+istr-1)*9+(k-1)*3+ivar
                        if (isub.gt.MXRGTRZ) then
                          write(*,9998)MXRGTRZ
 9998                     format(' Increase MXRGTRZ = ',i15)
                          write(*,*)'Dump out? 1-Yes/0-No'
                          read(*,*)idec
                          if (idec.ne.0) then
                            call dumpout_physics
                            call dumpout_hp3d('files/dumpc3Dhp')
                          endif
                          stop
                        endif
                        RGTR(isub) = costr(istr,ivar,k)
                        Xmin(ivar) = min(Xmin(ivar),costr(istr,ivar,k))
                        Xmax(ivar) = max(Xmax(ivar),costr(istr,ivar,k))
                      enddo
                    enddo
                  enddo
                  nstr=nstr+nstrl
!
!  .............end of loop through additional triangles
                enddo
!
!  .............save storage information
                isav=nstr+100*ifree
                IGTRNO(NRVISTR)=isav
              endif
!
!  ...........update list of triangles midpoints
#if HP3D_DEBUG
              if (iprint.eq.1) then
                write(*,7001) NRVISTR
 7001           format(' display_face: NRVISTR = ',i8)
              endif
#endif
              RGTRZ(NRVISTR) = (xtro(3,1)+xtro(3,2)+xtro(3,3))/3
!
!  ...........update list of visible curves when on boundary
              IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+nsid(l)
!
!  ...........encode the color
              IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+10*mpcol(l)
!
!  ...........encode nodes numbers
              if (Numlev.eq.0) then
!
!  .............encode element (middle node) numbers
                if (Iflagn.eq.1) then
                  select case(face_type(ntype,Ifc))
                  case(TRIA)
!
!  .................encode mid-face node number
                    nick = Mdle
                    if (l.eq.2) then
                      if((i.eq.(nsub1/3+1)).and.(j.eq.(nsub1/3+1))) then
                        IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+1000*nick
                      endif
                    endif
                 case(RECT)
!
!  .................encode mid-face node number
                    nick = Mdle
                    if (l.eq.2) then
                      if ((i.eq.(nsub1/2)).and.(j.eq.(nsub1/2))) then
                        IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+1000*nick
                      endif
                    endif
                  end select
!
!  .............encode face nodes numbers
                elseif (Iflagn.eq.2) then
                  select case(face_type(ntype,Ifc))
                  case(TRIA)
!
!  .................encode mid-face node number
                    nick = nface_nodes_global(7)
                    if (l.eq.2) then
                      if((i.eq.(nsub1/3+1)).and.(j.eq.(nsub1/3+1))) then
                        IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+1000*nick
                      endif
                    endif
!
!  .................encode nodes numbers
                    if (l.eq.1) then
                      if ((i.eq.1).and.(j.eq.1)) then
                        nick=nface_nodes_global(1)
                      elseif ((i.eq.1).and.(j.eq.nsub2)) then
                        nick=nface_nodes_global(3)
                      elseif ((i.eq.nsub1).and.(j.eq.1)) then
                        nick=nface_nodes_global(2)
                      elseif((i.eq.(nsub1/3+1)).and.(j.eq.(nsub1/3+1))) then
                        nick=nface_nodes_global(7)
                      elseif ((i.eq.1).and.(j.eq.(nsub2/2))) then
                        nick=nface_nodes_global(6)
                      elseif((j.eq.1).and.(i.eq.(nsub1/2))) then
                        nick=nface_nodes_global(4)
                      elseif((j.eq.nsub2).and.(i.eq.(nsub1/2))) then
                        nick=nface_nodes_global(5)
                      else
                        nick=0
                      endif
                      IGTRNO(NRVISTR) = IGTRNO(NRVISTR)+nick
                    endif
                 case(RECT)
!
!  .................encode mid-face node number
                    nick = nface_nodes_global(9)
                    if (l.eq.2) then
                      if ((i.eq.(nsub1/2)).and.(j.eq.(nsub1/2))) then
                        IGTRCU(NRVISTR) = IGTRCU(NRVISTR)+1000*nick
                      endif
                    endif
!
!  .................encode nodes numbers
                    if (l.eq.1) then
                      if ((i.eq.1).and.(j.eq.1)) then
                        nick=nface_nodes_global(1)
                      elseif ((i.eq.1).and.(j.eq.nsub2)) then
                        nick=nface_nodes_global(4)
                      elseif ((i.eq.nsub1).and.(j.eq.1)) then
                        nick=nface_nodes_global(2)
                      elseif ((i.eq.nsub1).and.(j.eq.nsub2)) then
                        nick=nface_nodes_global(3)
                      elseif((i.eq.(nsub1/2)).and.(j.eq.(nsub1/2))) then
                        nick=nface_nodes_global(9)
                      elseif ((i.eq.1).and.(j.eq.(nsub2/2))) then
                        nick=nface_nodes_global(8)
                      elseif ((i.eq.nsub1).and.(j.eq.(nsub2/2))) then
                        nick=nface_nodes_global(6)
                      elseif((j.eq.1).and.(i.eq.(nsub1/2))) then
                        nick=nface_nodes_global(5)
                      elseif((j.eq.nsub2).and.(i.eq.(nsub1/2))) then
                        nick=nface_nodes_global(7)
                      else
                        nick=0
                      endif
                      IGTRNO(NRVISTR) = IGTRNO(NRVISTR)+nick
                    endif
                  end select
                endif
              endif
#if HP3D_DEBUG
              if (iprint.eq.1) then
                write(*,8005) NRVISTR,IGTRNO(NRVISTR)
 8005           format('display_face: NRVISTR,IGTRNO(NRVISTR) = ',2i5)
              endif
#endif
!
            enddo
!
!  .......end of loop through lower/upper triangle
   20     enddo
!
!  .....end of loop through the vertical integer coordinate
        enddo
!
!  ...end of loop through the horizontal integer coordinate
      enddo
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,*) 'display_face: EXITING...'
        call pause
      endif
#endif
!
!
   end subroutine display_face

#endif
