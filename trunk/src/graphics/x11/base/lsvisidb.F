c----------------------------------------------------------------------
c
c   routine name       - lsvisidb
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 10
c
c   purpose            - routine makes a list of all visible
c                        triangles and stores rescaled coordinates
c                        of their end-points (vertices) in an array
c                        RGTR
c
c   arguments :
c     in:
c               Numlev = 0 display graphical image of the mesh
c                      > 0 - number of levels to plot solution values
c               Iflagn = 1 encode element (middle node) numbers
c                      = 2 encode face nodes numbers
c
c
c----------------------------------------------------------------------
c
      subroutine lsvisidb(Numlev,Iflagn)
c
      use element_data
      use data_structure3D
      use graphmod
c
#include "syscom.blk"
c
c
c  ...geometry and solution dof
      dimension xnod(3,MAXbrickH),
     .          zdofH(MAXEQNH,MAXbrickH),zdofE(MAXEQNE,MAXbrickE),
     .          zdofV(MAXEQNV,MAXbrickV),zdofQ(MAXEQNQ,MAXbrickQ)
c
      dimension xmax(3),xmin(3)
      dimension nedge_orient(12),nface_orient(6),norder(19)
c
c  ...element neighbors through faces
      dimension neig(4,6)
c
c  ...element nodes and their orientation
      dimension nodesl(27),norientl(27)
      dimension solev(NR_COLORS-10)

      dimension nrinter(8)
      dimension xlocinter(3,500)
      dimension xcenter(3)
c
c---------------------------------------------------------------
c
      data bigp,bign,small,one /1.d30,-1.d30,1.d-6,1.d0/
c
      iprint=0
c
c  ...prints vertex coordinates
      iprint1=0
c
c  ...find limiting values for contour plot
      if (Numlev.gt.0) call finlimb(numlev,solev)
c
c  ...set default bounds for the picture
      do ivar=1,3
        XEX(2*ivar-1)=bigp
        XEX(2*ivar)=bign
        xmin(ivar) = bigp
        xmax(ivar) = bign
      enddo
c
c  ...loop through elements
      NRVISTR = 0
      mdle=0
      do 20 iel=1,NRELES
        call nelcon(mdle, mdle)
        if (iprint.eq.1) then
          write(*,*) 'lsvisidb: mdle = ',mdle
        endif
c
c  .....check visibility (elements)

        if (IBOX_CUT.ne.0) then
          call find_center(mdle, xcenter)
          if (IBOX_CUT.gt.0) then
            if (xcenter(IBOX_CUT).gt.BOX_CUT(IBOX_CUT,1)) then
cc              write(*,*) 'box mdle ', mdle
              go to 20
            endif
          else
            if (xcenter(-IBOX_CUT).lt.BOX_CUT(-IBOX_CUT,1)) then
cc              write(*,*) 'box mdle ', mdle
              go to 20
            endif
          endif
        endif
c
        call locate(mdle,IGINV,NRINVBL, loc)
        if (loc.gt.0) go to 20
        call find_domain(mdle, ndom)
        if (iprint.eq.1) then
          write(*,7056) mdle,ndom,NDOMAIN(ndom)
 7056     format('lsvisidb: mdle,ndom,NDOMAIN(ndom) = ',i6,i3,i2)
        endif
        if (NDOMAIN(ndom).eq.0) go to 20
c
        call find_orient(mdle, nedge_orient,nface_orient)
        call find_order(mdle, norder)
c
        call nodcor(mdle,xnod)
        if (iprint.eq.1) then
          write(*,*)'lsvisidb: VERTEX COORDINATES = '
          do ivar=1,3
            write(*,8001) xnod(ivar,1:nvert(NODES(mdle)%type))
 8001       format(8(f8.5,2x))
          enddo
          call pause
        endif
c
c  .....get neighbors
        call find_neig(mdle, neig)
        if (iprint.eq.1) then
          write(*,8002) neig(1,1:nface(NODES(mdle)%type))
 8002     format('lsvisidb: neig = ',6(i8,2x))
          call pause
        endif
c
c  .....get node numbers
        call elem_nodes(mdle, nodesl,norientl)
c
c  .....determine element dof
        call solelm(mdle, zdofH,zdofE,zdofV,zdofQ)
c
c  .....check whether the element is sliced by the cutting plane
        icut=0
        daux=1
        do iv=1,nvert(NODES(mdle)%type)
          cl = CLPL(1)*xnod(1,iv) + CLPL(2)*xnod(2,iv) +
     .         CLPL(3)*xnod(3,iv) + CLPL(4)
          if (iv.eq.1)then
            daux1=cl
          else
            daux=cl*daux1
          endif
          if (daux.le.0) then
            icut=1
            exit
          endif
        enddo
c
c  .....skip elements totally in front of the cutting plane
        if ((icut.eq.0).and.(cl.gt.0)) go to 20
c
        nrintertot=0
c
c  .....loop through element faces
        do 10 if=1,nface(NODES(mdle)%type)
          if (iprint.eq.1) then
            write(*,7003) if
 7003       format('lsvisidb: if = ',i2)
          endif
c
c  .......for not sliced elements
          if (icut.eq.0) then
c
c  .........check if adjacent to the boundary or to the invisible
c           element, continue only if it is
c
c  .........flag : invisible
ccc            ivis=1
            ivis=0
            if (neig(1,if).lt.0) then
              call locate(-neig(1,if),IGINV,NRINVBL, iloc)
              if (iloc.gt.0) ivis=0
              call find_domain(-neig(1,if), ndom)
              if (NDOMAIN(ndom).eq.0) ivis=0
            elseif (neig(1,if).eq.0) then
              ivis=0
            else
ccc              do ineig=1,4
              do ineig=1,1
                call locate(neig(ineig,if),IGINV,NRINVBL, iloc)
                if (iloc.gt.0) ivis=0
                call find_domain(neig(ineig,if), ndom)
                if (NDOMAIN(ndom).eq.0) ivis=0
              enddo
            endif
            if (ivis.eq.1) go to 10
c
c  .......end for uncut elements
          endif
c
c  .......prepare the list of triangles for the face
          call display_face(Numlev,mdle,if,nodesl,
     .                      nedge_orient,nface_orient,norder,Iflagn,
     .                      xnod,zdofH,zdofE,zdofV,zdofQ,
     .                      solev,xmin,xmax,
     .                      nrintertot,nrinter(if),xlocinter)
c
   10   continue
c
        if ((icut.eq.1).and.(nrintertot.gt.0)) then
c
          if (nrintertot.ge.500) then
            write(*,7001)
 7001       format('lsvisidb: TOO MANY INTERSECTION POINTS')
            stop 1
          endif
c
cc          call display_slice(Numlev,eltype,mdle
cc     .                 ,nodesl,norder,xnod,zdofH,solev,xmin,xmax
cc     .                 ,nrintertot,nrinter,xlocinter)

        endif
   20 continue
c
c  ...compute objects dimensions and coordinates of its
c     central point
      DIMOB(1:3)  = (xmax(1:3)-xmin(1:3))/2.d0
      XCENTR(1:3) = (xmax(1:3)+xmin(1:3))/2.d0
      write(*,7035) xmin,xmax
 7035 format('lsvisidb: xmin,xmax = ',2(3f8.3,2x))
c
c  ...sort according to z-coordinate in order back-to-front
      if (NRVISTR.gt.0) then
        call sortz
      else
        write(*,7002)
 7002   format('lsvisidb: NO TRIANGLES TO BE DRAWN...')
      endif
c
      if (iprint.eq.1) then
        write(*,*) 'lsvisidb: EXITING...'
        call pause
      endif
c
c
      return
      end
