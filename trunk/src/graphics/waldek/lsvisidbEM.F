c----------------------------------------------------------------------
c
c   routine name       - lsvisidb
c
c----------------------------------------------------------------------
c
c   latest revision    - Jan 08
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
      subroutine lsvisidbEM(Numlev,Iflagn)
c
      use element_data
      use data_structure3D
      use graphmod
c
#include "syscom.blk"
#include "coltab.blk"
#include "gscale.blk"
#include "graphsp.blk"
#include "gparams.blk"
#include "csets.blk"
c
c
c  ...geometry and solution dof
      dimension xnod(3,MAXbrickH),
     .          zdofH(MAXEQNH,MAXbrickH),zdofE(MAXEQNE,MAXbrickE),
     .          zdofV(MAXEQNV,MAXbrickV),zdofQ(MAXEQNQ,MAXbrickQ)
c
      dimension xmax(3),xmin(3)
      dimension nedge_orient(12),nface_orient(6),norder(19)
      dimension neig(4,6),nodesl(27)
      dimension solev(NR_COLORS-10)

      dimension nrinter(8)
      dimension xlocinter(3,500)
c
c---------------------------------------------------------------
c  ...zooming
      common /cizoo1/ izoo
      common /cizu2/ xmi,xma,ymi,yma,wmi,wma
      real*8, allocatable :: xnod0(:,:)
      dimension ibf(6)

      common /cnocuts/ NOCUTS
c
c
      data bigp,bign,small,one /1.d30,-1.d30,1.d-6,1.d0/
c
      iprint=0
c
c  ...prints vertex coordinates
      iprint1=0
c
c  ...find limiting values for contour plot
      if (Numlev.gt.0) call finlimbEM(numlev,solev)
c

      if(izoo.eq.1)then
        allocate(xnod0(3,MAXbrickH), STAT=i1)
      endif

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
          write(*,*) 'lsvisidb: iel,mdle = ',iel,mdle
        endif

c          write(*,1000)iel,mdle,NODES(mdle)%type
c          write(*,1002)iel,mdle,ELEMS(iel)%nodes(27)
c 1000     format('lsvisidb: iel,mdle = ',2i6,a5)
c 1002     format('lsvisidb: iel,mdle = ',3i6)
c




c    ...if elsets are activated, verify if mdle belongs to elset ISETACT
c       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        if(ISETACT.gt.0)
     .    then
c          ...find initial element ancestor:
              call findan(mdle, nelanc,ndummy)
c
              do i=1,NRELSET(1,ISETACT)
                 iels = NRELSET(2,ISETACT) + i-1

c                if(nelanc.eq.ELEMS(NELISET(iels))%nodes(27))then
                 if(nelanc.eq.NELISET(iels))then
                   go to 1001
                 endif

                 enddo
              go to 20
c
          elseif(ISETACT.lt.0)then
c          ...find initial element ancestor:
              call findan(mdle, nelanc,ndummy)
c
              do i=1,NRELSET(1,-ISETACT)
                 iels = NRELSET(2,-ISETACT) + i-1

                 if(nelanc.eq.ELEMS(NELISET(iels))%nodes(27))then
                 endif

                 if(nelanc.eq.ELEMS(NELISET(iels))%nodes(27)) go to 20
                 enddo
              go to 1001
c
          endif
 1001     continue
c------------------------------------




c  .....check visibility (elements)
        call locate(mdle,IGINV,NRINVBL, loc)
        if (loc.gt.0) go to 20
        call find_domain(mdle, ndom)
ccccc   if (NDOMAIN(ndom).eq.0) go to 20


c        if (NDOMAIN(ndom).eq.0)then
c          write(*,*)'avoid mdle,iel=',mdle,iel ; call pause
c          go to 20
c        endif
c
        call find_orient(mdle, nedge_orient,nface_orient)
        call find_order(mdle, norder)
        call nodcor(mdle, xnod)


              small=1.d-5
              if(izoo.eq.1) then
               xm0 = 1.e+10
               xm1 =-1.e+10
               ym0 = 1.e+10
               ym1 =-1.e+10
               wm0 = 1.e+10
               wm1 =-1.e+10
               if(ELEMS(mdle)%type.eq.'tetr')naux=4
               if(ELEMS(mdle)%type.eq.'pris')naux=6
               do    i=1,naux
                     xm0 = min(xm0,xnod(1,i))
                     xm1 = max(xm1,xnod(1,i))
                     ym0 = min(ym0,xnod(2,i))
                     ym1 = max(ym1,xnod(2,i))
                     wm0 = min(wm0,xnod(3,i))
                     wm1 = max(wm1,xnod(3,i))
                     enddo
               if(xm0 .lt. xmi-small) go to 20
               if(xm1 .gt. xma+small) go to 20
               if(ym0 .lt. ymi-small) go to 20
               if(ym1 .gt. yma+small) go to 20
               if(wm0 .lt. wmi-small) go to 20
               if(wm1 .gt. wma+small) go to 20
               endif




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
 8002     format(6(i8,2x))
          call pause
        endif


c  ...verify if 'mdle' is on boundary of group of elments being drawn
c     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (izoo.eq.1) then
c
c    ...loop through 6 faces
        if(NODES(mdle)%type.eq.'mdlp')ns0=4
        if(NODES(mdle)%type.eq.'mdln')ns0=5
        do is=1,ns0
           if(neig(1,is).eq.0)then
             ibf(is)=1
             go to 66
           endif
c
c       ...loop though 4 possible neighbors through face
           do in=1,4
              call nodcor(neig(in,is), xnod0)
               xm0 = 1.e+10
               xm1 =-1.e+10
               ym0 = 1.e+10
               ym1 =-1.e+10
               wm0 = 1.e+10
               wm1 =-1.e+10

               if(ELEMS(neig(in,is))%type.eq.'tetr')naux=4
               if(ELEMS(neig(in,is))%type.eq.'pris')naux=6
               do    i=1,naux
                     xm0 = min(xm0,xnod0(1,i))
                     xm1 = max(xm1,xnod0(1,i))
                     ym0 = min(ym0,xnod0(2,i))
                     ym1 = max(ym1,xnod0(2,i))
                     wm0 = min(wm0,xnod0(3,i))
                     wm1 = max(wm1,xnod0(3,i))
                     enddo
c           ...if at least one vertex of at least one neighboris outside
               if(xm0 .lt. xmi-small .or. xm1 .gt. xma+small .or.
     .            ym0 .lt. ymi-small .or. ym1 .gt. yma+small .or.
     .            wm0 .lt. wmi-small .or. wm1 .gt. wma+small)then

                 if(nocuts.eq.0) ibf(is) = 1
                 go to 66
               endif

              enddo
c          ...end of loop through neighbors across face

 66        continue
           enddo
c       ...end of loop through faces
      endif




c
c  .....get node numbers
        call find_nodes(mdle, nodesl)
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
ccccc   if ((icut.eq.0).and.(cl.gt.0)) go to 20
        if ((icut.eq.0).and.(cl.gt.0))then
ccc          write(*,*)'avoid at 2' ; call pause
          go to 20
        endif
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
            ivis=1
            if (neig(1,if).lt.0) then
              call locate(-neig(1,if),IGINV,NRINVBL, iloc)
              if (iloc.gt.0) ivis=0
              call find_domain(-neig(1,if), ndom)
              if (NDOMAIN(ndom).eq.0) ivis=0
            elseif (neig(1,if).eq.0) then
              ivis=0
            else
              do ineig=1,4
                call locate(neig(ineig,if),IGINV,NRINVBL, iloc)
                if (iloc.gt.0) ivis=0
                call find_domain(neig(ineig,if), ndom)
                if (NDOMAIN(ndom).eq.0) ivis=0
              enddo
            endif

c        ...face is visible if on boundary
            if(ibf(if).ne.0)ivis=0

c        ...brute force!!!
            ivis=0

ccccc       if (ivis.eq.1) go to 10
            if (ivis.eq.1)then
ccc              write(*,*)'avoid at 3' ; call pause
              go to 10
            endif
c
c  .......end for uncut elements
          endif
c
c  .......prepare the list of triangles for the face
          call display_face(Numlev,NODES(mdle)%type,mdle,if,nodesl,
     .                      nedge_orient,nface_orient,norder,Iflagn,
     .                      xnod,zdofH,solev,xmin,xmax,
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

      if(izoo.eq.1)then
        deallocate(xnod0, STAT=i1)
      endif


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
