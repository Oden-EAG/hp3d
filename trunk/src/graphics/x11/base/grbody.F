#if HP3D_USE_X11

c----------------------------------------------------------------------
c
c   routine name       - grbody
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 2023
c
c   purpose            - routine displays 3D FEM grid or solution
c
c   arguments
c        in:
c                Iwind - type of window (for colormap)
c
c----------------------------------------------------------------------
c
      subroutine grbody(Iwind)
c
      use GMP , only: NRDOMAIN
      use physics , only : NR_PHYSA,PHYSA
      use data_structure3D , ONLY: NRELES
      use graphmod
c
      implicit none
c
      integer :: Iwind
c
      real(8) :: xcl(3),xmod(2)
      integer :: igv(10),im(2)
c
      integer :: i,ibl,iblno,ichoice,idefcol,idom,iel,iflagn,inickbl
      integer :: input_mode,ivis,mdle,mdlep,nrelem,nrsub1,numlev
      real(8) :: psi,q,q1,rnsc,theta
c
      real(8), parameter :: pi = 3.14159265358979312d0
      real(8), parameter :: bigp  =  1.d30
c
c  ...set initial coordinates of section plane
      CLPL(1)=1.d0
      CLPL(2)=0.d0
      CLPL(3)=0.d0
      CLPL(4)=-1.d10
c
c  ...initialize array of invisible elements
      NRINVBL = 0
      IGINV(1) = 0
c
c  ...set initial scaling constants
      DIMIM = bigp
      XCIM(1) = 0.d0
      XCIM(2) = 0.d0
c
c  ...set up number of line segments for drawing a curve
      nrsub1=2
      NRSUB = 2**nrsub1
      DX=1.d0/NRSUB
c
c  ...determine data for projection...

      theta = pi/5.d0+pi
      psi   = pi/6.d0
C       theta = pi/5.d0
C       psi = pi/4.d0
c
      RN(1) = cos(psi)*cos(theta)
      RN(2) = cos(psi)*sin(theta)
      RN(3) = sin(psi)
c
c  ...set all domains to be displayed
      if (NRDOMAIN.gt.MAXNRDOMAIN) then
        write(*,*) 'grbody: INCREASE MAXNRDOMAIN'
        stop 1
      endif
      NDOMAIN(1:NRDOMAIN)=1
ccc      NDOMAIN(1) = 1
c
c  ...choose the meaning of colours
      write(*,*) 'DEFINE COLORS :'
      write(*,*) '0-ORDER OF P-APPROXIMATION WITH ELEMENTS NOs'
      write(*,*) '1-ORDER OF P-APPROXIMATION WITH NODES NOs'
      write(*,*) '2-BOUNDARY CONDITIONS FLAGS'
      write(*,*) '3-SOLUTION VALUES'
      read(*,*) idefcol
c
      select case(idefcol)
      case(0)
        numlev=0; iflagn=1
      case(1)
        numlev=0; iflagn=2
      case(2)
        numlev=0; iflagn=3
        write(*,*) 'SET THE ATTRIBUTE NUMBER'
        do i=1,NR_PHYSA
          write(*,9999)i,PHYSA(i)
9999      format(1x,i1,' - ',a5)
        enddo
        read(*,*) NRPHY_DISP
      case(3)
        write(*,*) 'GIVE THE NUMBER OF LEVELS FOR SOLUTION PLOT '
        read(*,*) numlev
c  .....select quantity to be displayed
        call soldis_select
      end select
c
      write(*,*) '...PLEASE WAIT, PREPARING IMAGE'
c
c  ...create global image
      call cartobs
c
      call lsvisidb(numlev,iflagn)
      XCIM(1) = XCENTR(1)
      XCIM(2) = XCENTR(2)
c
c  ...open the window
      call selwin(Iwind)
c
c  ...decide whether to scale by x or y axes
      q = DIMOB(1)/DIMOB(2)
      q1 = (xlength-2*rmargin)/(ylength-2*rmargin)
      if (q.gt.q1) then
        DIMIM = DIMOB(1)
        SIZE  = (xlength-2*rmargin)/2.d0
      else
        DIMIM = DIMOB(2)
        SIZE  = (ylength-2*rmargin)/2.d0
      endif
      XCWIN(1) = rmargin + (xlength-2*rmargin)/2.d0
      XCWIN(2) = rmargin + (ylength-2*rmargin)/2.d0
c
c  ...display image
      call dpvisid(numlev*10+0)
c
c  ...display menu in infinite loop
   10 continue
c
c  ...close window
      write(*,*)'grbody: call dpborder with numlev=',numlev
      call dpborder(numlev)
      write(*,*) 'PLEASE CLICK THE MOUSE INSIDE THE GRAPHICS'
      write(*,*) 'WINDOW TO CONTINUE ...'
      call closwind(1)
c
      write(*,*) 'SELECT OPTION :'
      write(*,*) '0 - EXIT'
      write(*,*) '1 - CHANGE THE POINT OF VIEW (WITH RESCALING)'
      write(*,*) '2 - CHANGE THE POINT OF VIEW (WITHOUT RESCALING)'
      write(*,*) '3 - CHANGE THE CENTRAL POINT OF THE IMAGE AND',
     .           ' RESCALE THE IMAGE'
ccc      write(*,*) '3 - CHANGE THE CENTRAL POINT OF THE IMAGE'
      write(*,*) '4 - RESCALE THE IMAGE'
      write(*,*) '5 - DRAW A CROSS-SECTION'
      write(*,*) '6 - DISPLAY THE WHOLE OBJECT'
      write(*,*) '7 - MAKE SOME ELEMENTS INVISIBLE'
      write(*,*) '8 - LEAVE ONLY FEW ELEMENTS VISIBLE'
      write(*,*) '9 - MAKE ALL ELEMENTS VISIBLE'
      write(*,*) '10 - CHANGE THE DEGREE OF COARSENESS'
      if (numlev.gt.0) then
        write(*,*) '11 - CHANGE THE NUMBER OF LEVELS FOR SOLUTION PLOT'
      endif
      write(*,*) '12 - SELECT DOMAINS TO DISPLAY'
      write(*,*) '13 - ROUGH CUT'
      read(*,*) ichoice
c
      select case(ichoice)
      case(0)
        return
c
c  ...change the point of view...
      case(1,2)
        write(*,6024) psi/pi*180,theta/pi*180,RN(1:3)
 6024   format('CURRENT psi,theta,RN = ',2f6.1, 2x, 3f8.3)
        write(*,*) 'SET NEW psi,theta IN DEGREES'
        read(*,*)  psi,theta
        psi = psi/180*pi
        theta = theta/180*pi
        RN(1) = cos(psi)*cos(theta)
        RN(2) = cos(psi)*sin(theta)
        RN(3) = sin(psi)
        write(*,6025) psi/pi*180,theta/pi*180,RN(1:3)
 6025   format('NEW psi,theta,RN = ',2f6.1, 2x, 3f8.3)
        call cartobs
        if (ichoice.eq.1) DIMIM = bigp
c
c  ...change the central point of the image
      case(3)
        write(*,*) 'CLICK AT THE NEW CENTRAL POINT'
        call xmousepos(im(1),im(2))
        XCIM(1) = (im(1) - XCWIN(1))/SIZE*DIMIM + XCIM(1)
        XCIM(2) = (im(2) - XCWIN(2))/SIZE*DIMIM + XCIM(2)
        write(*,*) 'XCIM = ',XCIM
        write(*,*) 'SET MAGNIFICATION FACTOR'
        read(*,*) rnsc
        rnsc=dabs(rnsc)
        DIMIM=DIMIM/rnsc
ccc        write(*,*) 'SET NEW COORDINATES OF THE CENTRAL POINT.'
ccc        write(*,*) 'MUST BE BETWEEN -1 AND 1 (WRT ACTUAL WINDOW).'
ccc        read(*,*) rncp1,rncp2
ccc        write(*,*) 'XCIM OLD = ',XCIM
ccc        XCIM(1)=XCIM(1)+rncp1*DIMIM
ccc        XCIM(2)=XCIM(2)+rncp2*DIMIM
ccc        write(*,*) 'XCIM NEW = ',XCIM
c
c  ...rescale the image...
      case(4)
        write(*,*) 'SET MAGNIFICATION FACTOR'
        read(*,*) rnsc
        rnsc=dabs(rnsc)
        DIMIM=DIMIM/rnsc
c
c  ...draw a cross-section...
      case(5)
        write(*,*) 'THE EXTREME VALUES OF COORDINATES '
        write(*,*) 'ARE (IN ORIGINAL SYSTEM OXYZ) :'
        write(*,*) 'X - ',XEX(1),XEX(2)
        write(*,*) 'Y - ',XEX(3),XEX(4)
        write(*,*) 'Z - ',XEX(5),XEX(6)
        write(*,*) 'SET COORDINATES OF A VECTOR'
        write(*,*) 'NORMAL TO THE SECTION PLANE'
        read(*,*) CLPL(1),CLPL(2),CLPL(3)
        write(*,*) 'SET COORDINATES OF A POINT'
        write(*,*) 'BELONGING TO THE SECTION PLANE'
        read(*,*) xcl(1),xcl(2),xcl(3)
        CLPL(4)=-CLPL(1)*xcl(1)-CLPL(2)*xcl(2)-CLPL(3)*xcl(3)
c
c  ...display the whole object...
      case(6)
        DIMIM = bigp
        CLPL(1) = 1.d0
        CLPL(2) = 0.d0
        CLPL(3) = 0.d0
        CLPL(4) = -1.d10
        IBOX_CUT = 0
c
c  ...make some elements invisible...
      case(7)
        write(*,*) 'DISPLAY ELEMENT NUMBERS? (0-No/1-Yes)'
        read(*,*) iblno
        if (iblno.eq.1) then
          call selwin(iwind)
          call dpvisid(numlev*10+iblno)
          write(*,*) 'PLEASE CLICK THE MOUSE INSIDE THE GRAPHICS'
          write(*,*) 'WINDOW TO CONTINUE ...'
          call closwind(1)
        endif

cldem 12.99
        write(*,*) 'GIVE THE NUMBER OF ELEMENTS TO MAKE INVISIBLE'
        read(*,*) inickbl
        write(*,*) 'GIVE THE ELEMENT NUMBERS TO MAKE INVISIBLE'
        do 15 ibl=1,inickbl
          read(*,*) nrelem
          NRINVBL = NRINVBL + 1
          IGINV(NRINVBL) = nrelem
   15   continue
        go to 16
        write(*,*) 'CLICK ON ELEMENTS TO BE INVISIBLE'
        write(*,*) 'CLICK OUTSIDE THE MESH TO TERMINATE INPUT'
        mdlep=-1
        do while (mdlep.ne.0)
c
c  .......read in the element number and ref kind
          call xmousepos(im(1),im(2))
          xmod(1) = (im(1) - XCWIN(1))/SIZE*DIMIM + XCIM(1)
          xmod(2) = (im(2) - XCWIN(2))/SIZE*DIMIM + XCIM(2)
          call elem_no(xmod, mdle)
          write(*,*) 'SELECTED xmod,mdle = ',xmod(1:2),mdle
          if (mdle.ne.0) then
            NRINVBL = NRINVBL + 1
            IGINV(NRINVBL) = mdle
          endif
          mdlep=mdle
        enddo
 16     continue
c
c  ...leave only a few elements visible...
      case(8)
        write(*,*) 'DISPLAY ELEMENT NUMBERS? (0-NO/1-YES)'
        read(*,*) iblno
        if (iblno.eq.1) then
          call selwin(iwind)
          call dpvisid(numlev*10+iblno)
          write(*,*) 'PLEASE CLICK THE MOUSE INSIDE THE GRAPHICS'
          write(*,*) 'WINDOW TO CONTINUE ...'
          call closwind(1)
        endif
c
        write(*,*) 'INPUT BY MOUSE OR NUMBERS (1/2)'
        read(*,*) input_mode
        select case(input_mode)
        case(1)
          write(*,*) 'CLICK ON ELEMENTS TO STAY VISIBLE ( .le.10 )'
          write(*,*) 'CLICK OUTSIDE THE MESH TO TERMINATE INPUT'
          inickbl=0
          mdlep=-1
          do while (mdlep.ne.0)
c
c  .........read in the element number
            call xmousepos(im(1),im(2))
            xmod(1) = (im(1) - XCWIN(1))/SIZE*DIMIM + XCIM(1)
            xmod(2) = (im(2) - XCWIN(2))/SIZE*DIMIM + XCIM(2)
            call elem_no(xmod, mdle)
            if (mdle.ne.0) then
              inickbl = inickbl + 1
              igv(inickbl) = mdle
            endif
            mdlep=mdle
          enddo
        case(2)
          write(*,*) 'GIVE THE NUMBER OF ELEMENTS TO STAY VISIBLE'
          read(*,*) inickbl
          write(*,*) 'GIVE THE ELEMENT NUMBERS TO STAY VISIBLE'
          do ibl=1,inickbl
            read(*,*) igv(ibl)
          enddo
        end select
c
        NRINVBL = 0
        mdle=0
        do iel=1,NRELES
          call nelcon(mdle, mdle)
          call locate(mdle,igv,inickbl, ivis)
          if (ivis.eq.0) then
            NRINVBL = NRINVBL+1
            if (NRINVBL.gt.MAXNRINVBL) then
              write(*,*) 'grbody: INCREASE MAXNRINVBL '
              stop 1
            endif
            IGINV(NRINVBL) = mdle
          endif
        enddo
c
c  ...make all elements visible
      case(9)
        NRINVBL = 0
        IGINV(1) = 0
c
c  ...change the degree of coarseness...
      case(10)
        write(*,*)'SET NEW DEGREE OF COARSENESS FOR OVALS (2,3,4)'
        write(*,*)'CURRENT VALUE IS ', nrsub1
        read(*,*) nrsub1
        if (nrsub1.gt.4) then
          write(*,*)'CHANGED TO 4'
          nrsub1 = 4
        else if (NRSUB.lt.2) then
          write(*,*)'CHANGED TO 2'
          nrsub1 = 2
        endif
        NRSUB=2**nrsub1
        DX=1.d0/NRSUB
c
c  ...change the number of levels for the solution plot...
      case(11)
        write(*,*) 'GIVE NEW NUMBER OF LEVELS FOR SOLUTION PLOT '
        read(*,*) numlev
c
c  ...select domains to display...
      case(12)
ccc     ...FF Jul 15: commented this call as it is useless and
ccc        inconvenient (it requires creating an empty routine with
ccc        the name display_domains EVERY time a problem is set).
ccc        call display_domains
        write(*,*)'INPUT DOMAIN NUMBERS (0-Exit)'
        idom=1 ; NDOMAIN(1:NRDOMAIN)=0
        do while (idom.gt.0)
          write(*,*) 'DOMAIN ='
          read(*,*)idom
          if (idom.gt.0) NDOMAIN(idom)=1
        enddo
        write(*,7038) NDOMAIN(1:NRDOMAIN)
 7038   format('object3: NDOMAIN = ',20i2)

      case(13)
        write(*,*) 'SET POINT '
        read(*,*) BOX_CUT(1:3,1)
        write(*,*) 'SET CUTTING DIRECTION ( 0:no-cut, 1:x, 2:y, 3:z )'
        read(*,*) IBOX_CUT
      end select
      write(*,*) '...PLEASE WAIT, PREPARING IMAGE'
      call lsvisidb(numlev,iflagn)
c
c  ...rescale when necessary
      if((ichoice.eq.1).or.(ichoice.eq.6)) then
        XCIM(1) = XCENTR(1)
        XCIM(2) = XCENTR(2)
c
c  .....decide whether to scale by x or y axes
        q = DIMOB(1)/DIMOB(2)
        q1 = xlength/ylength
        if (q.gt.q1) then
          DIMIM = DIMOB(1)
          SIZE = xlength/2.d0
        else
          DIMIM = DIMOB(2)
          SIZE = ylength/2.d0
        endif
c
      endif
c
c  ...open the window again
      call selwin(iwind)
c
      select case(idefcol)
      case(0)
        write(*,*) 'DISPLAY ELEMENT NUMBERS? (0-No/1-Yes) '
        read(*,*) iblno
      case(1)
        write(*,*) 'DISPLAY NODE NUMBERS? (0-No/2-Yes) '
        read(*,*) iblno
      case(2)
        write(*,*) 'DISPLAY ELEMENT NUMBERS? (0-No/1-Yes) '
        read(*,*) iblno
      case(3)
        iblno=0
      end select
      call dpvisid(numlev*10+iblno)
c
      goto 10
c
c
      end subroutine grbody

#endif
