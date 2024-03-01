#if HP3D_USE_X11

!----------------------------------------------------------------------
!
!   routine name       - lsvisidb
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine makes a list of all visible
!                        triangles and stores rescaled coordinates
!                        of their end-points (vertices) in an array
!                        RGTR
!
!   arguments :
!     in:
!               Numlev = 0 display graphical image of the mesh
!                      > 0 number of levels to plot solution values
!               Iflagn = 1 encode element (middle node) numbers
!                      = 2 encode face nodes numbers
!
!----------------------------------------------------------------------
#include "typedefs.h"
   subroutine lsvisidb(Numlev,Iflagn)
!
      use element_data
      use data_structure3D
      use graphmod
!
      implicit none
!
      integer, intent(in) :: Numlev,Iflagn
!
!  ...geometry and solution dof
      real(8) :: xnod(3,MAXbrickH)
      VTYPE   :: zdofH(MAXEQNH,MAXbrickH),zdofE(MAXEQNE,MAXbrickE),  &
                 zdofV(MAXEQNV,MAXbrickV),zdofQ(MAXEQNQ,MAXbrickQ)
!
      real(8) :: xmax(3),xmin(3)
      integer :: nedge_orient(12),nface_orient(6),norder(19)
!
!  ...element neighbors through faces
      integer :: neig(4,6)
!
!  ...element nodes and their orientation
      integer :: nodesl(27),norientl(27)
      real(8) :: solev(NR_COLORS-10)

      integer :: nrinter(8)
      real(8) :: xlocinter(3,500)
      real(8) :: xcenter(3)
!
      real(8) :: cl,daux,daux1
      integer :: icut,iel,ifc,iloc,ineig,iv,ivar,ivis
      integer :: loc,mdle,ndom,nrintertot
!
      real(8), parameter :: bigp  =  1.d30
      real(8), parameter :: bign  = -1.d30
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!
!---------------------------------------------------------------
!
!  ...find limiting values for contour plot
      if (Numlev.gt.0) call finlimb(Numlev,solev)
!
!  ...set default bounds for the picture
      do ivar=1,3
        XEX(2*ivar-1) = bigp
        XEX(2*ivar)   = bign
        xmin(ivar)    = bigp
        xmax(ivar)    = bign
      enddo
!
!  ...loop through elements
      NRVISTR = 0
      mdle=0
      do 20 iel=1,NRELES
        call nelcon(mdle, mdle)
#if HP3D_DEBUG
        if (iprint.eq.1) then
          write(*,*) 'lsvisidb: mdle = ',mdle
        endif
#endif
!
!  .....check visibility (elements)

        if (IBOX_CUT.ne.0) then
          call find_center(mdle, xcenter)
          if (IBOX_CUT.gt.0) then
            if (xcenter(IBOX_CUT).gt.BOX_CUT(IBOX_CUT,1)) then
!!              write(*,*) 'box mdle ', mdle
              go to 20
            endif
          else
            if (xcenter(-IBOX_CUT).lt.BOX_CUT(-IBOX_CUT,1)) then
!!              write(*,*) 'box mdle ', mdle
              go to 20
            endif
          endif
        endif
!
        call locate(mdle,IGINV,NRINVBL, loc)
        if (loc.gt.0) go to 20
        call find_domain(mdle, ndom)
#if HP3D_DEBUG
        if (iprint.eq.1) then
          write(*,7056) mdle,ndom,NDOMAIN(ndom)
 7056     format('lsvisidb: mdle,ndom,NDOMAIN(ndom) = ',i6,i3,i2)
        endif
#endif
        if (NDOMAIN(ndom).eq.0) go to 20
!
        call find_orient(mdle, nedge_orient,nface_orient)
        call find_order(mdle, norder)
!
        call nodcor(mdle,xnod)
#if HP3D_DEBUG
        if (iprint.eq.1) then
          write(*,*)'lsvisidb: VERTEX COORDINATES = '
          do ivar=1,3
            write(*,8001) xnod(ivar,1:nvert(NODES(mdle)%ntype))
 8001       format(8(f8.5,2x))
          enddo
          call pause
        endif
#endif
!
!  .....get neighbors
        call find_neig(mdle, neig)
#if HP3D_DEBUG
        if (iprint.eq.1) then
          write(*,8002) neig(1,1:nface(NODES(mdle)%ntype))
 8002     format('lsvisidb: neig = ',6(i8,2x))
          call pause
        endif
#endif
!
!  .....get node numbers
        call elem_nodes(mdle, nodesl,norientl)
!
!  .....determine element dof
        call solelm(mdle, zdofH,zdofE,zdofV,zdofQ)
!
!  .....check whether the element is sliced by the cutting plane
        icut=0
        daux=1
        do iv=1,nvert(NODES(mdle)%ntype)
          cl = CLPL(1)*xnod(1,iv) + CLPL(2)*xnod(2,iv) + &
               CLPL(3)*xnod(3,iv) + CLPL(4)
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
!
!  .....skip elements totally in front of the cutting plane
        if ((icut.eq.0).and.(cl.gt.0)) go to 20
!
        nrintertot=0
!
!  .....loop through element faces
        do 10 ifc=1,nface(NODES(mdle)%ntype)
#if HP3D_DEBUG
          if (iprint.eq.1) then
            write(*,7003) ifc
 7003       format('lsvisidb: ifc = ',i2)
          endif
#endif
!
!  .......for not sliced elements
          if (icut.eq.0) then
!
!  .........check if adjacent to the boundary or to the invisible
!           element, continue only if it is
!
!  .........flag : invisible
!!!            ivis=1
            ivis=0
            if (neig(1,ifc).lt.0) then
              call locate(-neig(1,ifc),IGINV,NRINVBL, iloc)
              if (iloc.gt.0) ivis=0
              call find_domain(-neig(1,ifc), ndom)
              if (NDOMAIN(ndom).eq.0) ivis=0
            elseif (neig(1,ifc).eq.0) then
              ivis=0
            else
!!!              do ineig=1,4
              do ineig=1,1
                call locate(neig(ineig,ifc),IGINV,NRINVBL, iloc)
                if (iloc.gt.0) ivis=0
                call find_domain(neig(ineig,ifc), ndom)
                if (NDOMAIN(ndom).eq.0) ivis=0
              enddo
            endif
            if (ivis.eq.1) go to 10
!
!  .......end for uncut elements
          endif
!
!  .......prepare the list of triangles for the face
          call display_face(Numlev,mdle,ifc,nodesl,                  &
                            nedge_orient,nface_orient,norder,Iflagn, &
                            xnod,zdofH,zdofE,zdofV,zdofQ,            &
                            solev,xmin,xmax,                         &
                            nrintertot,nrinter(ifc),xlocinter)
!
   10   continue
!
        if ((icut.eq.1).and.(nrintertot.gt.0)) then
!
          if (nrintertot.ge.500) then
            write(*,7001)
 7001       format('lsvisidb: TOO MANY INTERSECTION POINTS')
            stop 1
          endif
!
!!          call display_slice(Numlev,eltype,mdle                    &
!!                       ,nodesl,norder,xnod,zdofH,solev,xmin,xmax   &
!!                       ,nrintertot,nrinter,xlocinter)

        endif
   20 continue
!
!  ...compute objects dimensions and coordinates of its
!     central point
      DIMOB(1:3)  = (xmax(1:3)-xmin(1:3))/2.d0
      XCENTR(1:3) = (xmax(1:3)+xmin(1:3))/2.d0
      write(*,7035) xmin,xmax
 7035 format('lsvisidb: xmin,xmax = ',2(3f8.3,2x))
!
!  ...sort a!!ording to z-coordinate in order back-to-front
      if (NRVISTR.gt.0) then
        call sortz
      else
        write(*,7002)
 7002   format('lsvisidb: NO TRIANGLES TO BE DRAWN...')
      endif
!
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,*) 'lsvisidb: EXITING...'
        call pause
      endif
#endif
!
!
   end subroutine lsvisidb

#endif
