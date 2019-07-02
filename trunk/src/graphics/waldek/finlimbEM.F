c----------------------------------------------------------------------
c
c   routine name       - finlimb
c
c----------------------------------------------------------------------
c
c   latest revision    - Jan 08
c
c   purpose            - routine finds 'numlev'- levels of solution values
c
c   arguments :
c     in:
c               Numlev - number of levels to plot solution values
c     out:
c               Solev  - limiting values for each level 
c
c----------------------------------------------------------------------
c
      subroutine finlimbEM(Numlev, Solev)
c
      use data_structure3D
#include "syscom.blk"
#include "coltab.blk"
#include "gparams.blk"
#include "gselect.blk"
#include "gscale.blk"
#include "graphsp.blk"
#include "csets.blk"
c
      dimension Solev(NR_COLORS-10)
      dimension nedge_orient(12),nface_orient(6),norder(19)
c
c  ...geometry and actual dof, BC flags
      dimension xnod(3,MAXbrickH),ibf(6)
      dimension zdofH(MAXEQNH,MAXbrickH)
      dimension zdofE(MAXEQNE,MAXbrickE)
      dimension zdofV(MAXEQNV,MAXbrickV)
      dimension zdofQ(MAXEQNQ,MAXbrickQ)
c
c  ...master face coordinates of a point
      dimension t(2)
c
c  ...master and physical coordinates of a point, solution value
      dimension xi(3),xpoint(3),
     .          zsolH(MAXEQNH)

c  ...zooming
      common /cizoo1/ izoo
      common /cizu2/ xmi,xma,ymi,yma,wmi,wma

c
c
      iprint=0
      if (iprint.eq.1) then
        write(*,*) 'finlimb: Numlev = ',Numlev
        call pause
      endif
c
c  ...increment in master element coordinates
      dxi = DX
c
c  ...in search for biggest and smallest value
      solmax =-1.d10
      solmin = 1.d10
c
c  ...loop through elements      
      mdle=0
      do iel=1,NRELES
c
        call nelcon(mdle, mdle)
c       write(*,*)'iel=',iel


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
c                   write(*,*)'finlimb, draw mdle=',mdle ; call pause
                   go to 1001
                 endif

                 enddo

c              write(*,*)'finlimb, avoid mdle=',mdle ; call pause
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


        call locate(mdle,IGINV,NRINVBL, loc)
        if (loc.gt.0) cycle
c
        call find_domain(mdle, ndom)
cwr07.15.09
ccccc   if (NDOMAIN(ndom).eq.0) cycle
c
        call find_bc(mdle, ibf)
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






        call solelm(mdle, zdofH,zdofE,zdofV,zdofQ)
        call celndof(NODES(mdle)%type,norder, 
     .               nrdofH,nrdofE,nrdofV,nrdofQ)
c
c  .....loop through element faces
        do if=1,nface(NODES(mdle)%type)
c
c  .......loop through points on the face
          nsub1=NRSUB
          do j=0,nsub1

            select case(face_type(NODES(mdle)%type,if))
            case('tria'); nsub2=NRSUB-j
            case('rect'); nsub2=NRSUB
            end select



            do i=0,nsub2
              t(1) = i*dxi; t(2) = j*dxi
              call compute_face(Numlev,NODES(mdle)%type,mdle,if,
     .                          nedge_orient,nface_orient,norder,
     .                          xnod,zdofH,t,
     .                          xpoint,val)
c
c  ...........update extremes  
              solmax = max(solmax,val)
              solmin = min(solmin,val)
            enddo
          enddo

c
c  .....end of loop through faces
        enddo
c        write(*,*)'end iel=',iel
c
 20   continue
c  ...end of loop through elements
      enddo
c
c
      write(*,7010) solmin,solmax
 7010 format('finlimb: EXTREME VALUES = ',2e16.6)
c
c      write(*,*) 'DO YOU WANT TO CHANGE THE RANGE OF COLORS ?'
c      write(*,*) '0....NO'
c      write(*,*) '1....USE COMMON BOUND FOR NEG AND POS VALUES'
c      write(*,*) '2....SET UP YOUR OWN BOUNDS'
c      read(*,*) idec
       idec=0
      select case(idec)
      case(1)
        solmax = max(abs(solmax),abs(solmin))
        solmin = -solmax
        write(*,7010) solmin,solmax
      case(2)
        write(*,*) 'finlinmb: GIVE LOWER AND UPPER BOUND'
        read(*,*) solmin,solmax
      end select
c
c  ...divide the whole range into levels
c
c  ...use the volume values....
      dsol = (solmax-solmin)/float(Numlev)
      Solev(1) = solmin
c
      do i=1,Numlev
        Solev(i+1) = Solev(1) +float(i)*dsol
      enddo
c
        write(*,*) 'finlimb: exit'
c        call pause
c
      end

