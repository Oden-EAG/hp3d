c----------------------------------------------------------------------
c
c   routine name       - grbody
c
c----------------------------------------------------------------------
c
c   latest revision    - Jan 08
c
c   purpose            - routine displays 3D FEM grid
c
c   arguments
c        in:
c                Iwind - type of window (for colormap)
c
c----------------------------------------------------------------------
c
      subroutine grbodyEM(Iwind)
c
      use GMP , only: NRDOMAIN
      use data_structure3D , ONLY: NRELES
      use graphmod
c
#include "syscom.blk"
c
#include "coltab.blk"
#include "gscale.blk"
#include "graphsp.blk"
#include "gtransf.blk"
#include "gwind.blk"
#include "gparams.blk"

#include "gselect.blk"
#include "csets.blk"

      common /phys_attribute/ no_attribute
c
      dimension xcl(3),igv(10),im(2),xmod(2)
c
c----------------------------------------------------------------------
c  ...zooming
      common /cizoo1/ izoo
      common /cizu2/ xmi,xma,ymi,yma,wmi,wma

c  ...communication with routine soldis (use spherical coordinates)
      common /csoldis_wr/ nsph


      common /cnocuts/ NOCUTS
c
      save irout5, numlev0, icomp, xcl, ixcl, icart, nrsub1, ielevat,
     .     RN1, REDUCE, NOSCALE, MONO_OLD, iblno, iedge
      save numlev

      save no_comp,ichoose
c----------------------------------------------------
      data irou5/0/
      data pi /3.14159265358979312d0/
      data bigp,bign,small /1.d30,-1.d30,1.d-14/

      if(irout5.ne.100)then
        irout5=100

c    ...no spherical coordinates for displaying
        nsph=0

        no_comp=1
        ichoose=1
        ISELECT = no_comp*1000+ichoose


c    ...set initial coordinates of section plane
        CLPL(1)=1.d0
        CLPL(2)=0.d0
        CLPL(3)=0.d0
        CLPL(4)=-1.d10
c
c    ...initialize array of invisible elements
        NRINVBL = 0
        IGINV(1) = 0
c
c    ...set initial scaling constants
        DIMIM = bigp
        XCIM(1) = 0.d0
        XCIM(2) = 0.d0
c
c    ...set up number of line segments for drawing a curve
        nrsub1=2
        NRSUB = 2**nrsub1
        DX=1.d0/NRSUB
c
c    ...determine data for projection...
        theta = pi/6
        psi = pi/6
        RN(1) = cos(psi)*cos(theta)
        RN(2) = cos(psi)*sin(theta)
        RN(3) = sin(psi)

c    ...display solution
        idefcol=3

        NUMLEV=14
        numlev0=14

        NOCUTS = 0
      endif

      ielevat=0
c
c  ...set all domains to be displayed
      if (NRDOMAIN.gt.MAXNRDOMAIN) then
        write(*,*) 'grbody: INCREASE MAXNRDOMAIN'
        stop 1
      endif
      NDOMAIN(1:NRDOMAIN)=0; NDOMAIN(1)=1



c
cc  ...choose the meaning of colours
c      write(*,*) 'DEFINE COLORS :'
c      write(*,*) '0-ORDER OF P-APPROXIMATION WITH ELEMENTS NOs'
c      write(*,*) '1-ORDER OF P-APPROXIMATION WITH NODES NOs'
c      write(*,*) '2-BOUNDARY CONDITIONS FLAGS'
c      write(*,*) '3-SOLUTION VALUES'
c      read(*,*) idefcol
c      select case(idefcol)
c      case(0)
c        numlev=0; iflagn=1
c      case(1)
c        numlev=0; iflagn=2
c      case(2)
c        numlev=0; iflagn=3
c        write(*,*) 'SET THE ATTRIBUTE NUMBER'
c        read(*,*) no_attribute
c      case(3)
c        write(*,*) 'GIVE THE NUMBER OF LEVELS FOR SOLUTION PLOT '
c        read(*,*) numlev
cc  .....select quantity to be displayed
c        call select
c      end select
c


   12 write(*,*)'*****************************************************'
      write(*,*)' '
      write(*,*)'           R O U T I N E   G R B O D Y E M'
      write(*,*)' '
      write(*,*)'*****************************************************'
      write(*,*)' '
      write(*,*)'MESH...................................DISPLAY MESH'
      write(*,*)' '
      write(*,*)'TWODIM.............................DISPLAY SOLUTION'
      write(*,*)'      EYE=?.......................OBSERVATION POINT'
      write(*,*)'      LEVEL=?......................NUMBER of LEVELS'
      write(*,*)'      DIVISION=?................................2-4'
      write(*,*)'      COMP=?.....................SOLUTION COMPONENT'
      write(*,*)'      ELEMENT.......................ELEMENT NUMBERS'
      write(*,*)'      NODE.............................NODE NUMBERS'
      write(*,*)' '
      write(*,*)'SLICE=?,?,?.................VECTOR OF SLICING PLANE'
      write(*,*)'    X=?,?,?..................POINT OF SLICING PLANE'
      write(*,*)' '
      write(*,*)'SLICE3D=?,?,?...............VECTOR OF SLICING PLANE'
      write(*,*)'CUTONLY=0/1..................IF DRAW ONLY ON SLICES'
      write(*,*)'END............................................EXIT'


c
      call ffldin(5,6,25,ier)
      if (ier.eq.1) go to 12
c
c---------------------------------------------------------------------
      call findop('ELSET' ,noc,loc)
       if (noc.gt.0)
     .  then
            ISETACT = int( rfldata(loc) )
        endif
c---------------------------------------------------------------------
      call findop('IFIRST' ,noc,loc)
       if (noc.gt.0)
     .  then
            IFIRST = int( rfldata(loc) )
        endif
c---------------------------------------------------------------------
      call findop('MESH' ,noc,loc)
      if (noc.gt.0)
     .  then
c        ...to display mesh:
            numlev = 0
            ielevat = 0
        endif
c---------------------------------------------------------------------
      call findop('BOUND' ,noc,loc)
      if (noc.gt.0)
     .  then
            ibound  = int( rfldata(loc) )
            numlev = 0
            ielevat = 0
        endif
c---------------------------------------------------------------------
      call findop('TWODIM' ,noc,loc)
      if (noc.gt.0)
     .  then
            iians = 2
            DIMIM = bigp
            CLPL(1) = 1.d0
            CLPL(2) = 0.d0
            CLPL(3) = 0.d0
            CLPL(4) = -1.d10
            numlev = numlev0
            ielevat = 0
        endif
c---------------------------------------------------------------------
      call findop('LEVEL' ,noc,loc)
      if (noc.gt.0)
     .  then
            numlev0= int( rfldata(loc) )
            numlev = numlev0
        endif
c---------------------------------------------------------------------
      call findop('DIVISION' ,noc,loc)
      if (noc.gt.0)
     .  then
            nrsub1 = int( rfldata(loc) )
ccccc       NRSUB = 2**nrsub1
            NRSUB = 2**(nrsub1+1)
            DX=1.d0/NRSUB
        endif
c---------------------------------------------------------------------
      call findop('COMP' ,noc,loc)
      if (noc.gt.0)
     .  then
            ichoose = int( rfldata(loc) )
        endif
c---------------------------------------------------------------------
      call findop('NRVAR' ,noc,loc)
      if (noc.gt.0)
     .  then
            no_comp = int( rfldata(loc) )
        endif
c---------------------------------------------------------------------
      call findop('EYE' ,noc,loc)
      if (noc.gt.0)
     .  then
            if(ielevat.eq.0)then
              RN(1) =  rfldata(loc)
              RN(2) =  rfldata(loc+1)
              RN(3) =  rfldata(loc+2)
              write(*,*)'rn=',rn
              DIMIM = bigp
            else
              write(*,*)'NO EYE OPTION IN DISPLAY3D: use EYE1 !!!'
            endif

c        ...necessary call routine cartobs
            icart = 0
        endif
c---------------------------------------------------------------------
      call findop('X' ,noc,loc)
      if (noc.gt.0)
     .  then
            XCL(1) =  rfldata(loc)
            XCL(2) =  rfldata(loc+1)
            XCL(3) =  rfldata(loc+2)
            ixcl = 1
        endif
c---------------------------------------------------------------------
      call findop('SLICE' ,noc,loc)
      if (noc.gt.0)
     .  then
            CLPL(1) =  rfldata(loc)
            CLPL(2) =  rfldata(loc+1)
            CLPL(3) =  rfldata(loc+2)
            CLPL(4)=-CLPL(1)*xcl(1)-CLPL(2)*xcl(2)-CLPL(3)*xcl(3)
        endif
c---------------------------------------------------------------------
      call findop('SLICE3D' ,noc,loc)
      if (noc.gt.0)
     .  then
            ielevat = 1
            numlev = 0
            CLPL(1) =  rfldata(loc)
            CLPL(2) =  rfldata(loc+1)
            CLPL(3) =  rfldata(loc+2)
            CLPL(4)=-CLPL(1)*xcl(1)-CLPL(2)*xcl(2)-CLPL(3)*xcl(3)
c        ...make RN coincide with CLPL:
            RN(1) = CLPL(1)
            RN(2) = CLPL(2)
            RN(3) = CLPL(3)
        endif
c---------------------------------------------------------------------
      call findop('SOLSCALE' ,noc,loc)
      if (noc.gt.0)
     .  then
            SOLSCALE =  rfldata(loc)
        endif
c---------------------------------------------------------------------
      call findop('ELEMENT' ,noc,loc)
      if (noc.gt.0)
     .  then
c        ...to display element numbers
            iblno = int( rfldata(loc) )
            write(*,*)'iblno=',iblno
        endif
c---------------------------------------------------------------------
      call findop('NODE' ,noc,loc)
      if (noc.gt.0)
     .  then
c        ...to display node numbers
            iblno = 2
        endif
c---------------------------------------------------------------------
      call findop('REDUCE' ,noc,loc)
      if (noc.gt.0)
     .  then
            REDUCE =  rfldata(loc)
        endif
c---------------------------------------------------------------------
      call findop('NOSCALE' ,noc,loc)
      if (noc.gt.0)
     .  then
            NOSCALE = int( rfldata(loc) )
        endif
c---------------------------------------------------------------------
      call findop('CUTONLY' ,noc,loc)
      if (noc.gt.0)
     .  then
            ICUTONLY = int( rfldata(loc) )
            write(*,*)'ICUTONLY=',icutonly
        endif
c---------------------------------------------------------------------
      call findop('MONO' ,noc,loc)
      if (noc.gt.0)
     .  then
            MONO = int( rfldata(loc) )
            write(*,*)'MONO=',MONO
        endif
c---------------------------------------------------------------------
      call findop('SCLOG' ,noc,loc)
      if (noc.gt.0)
     .  then
            SCLOG = rfldata(loc)
        endif
c---------------------------------------------------------------------
      call findop('INVSCALE' ,noc,loc)
      if (noc.gt.0)
     .  then
            INVSCALE = int(rfldata(loc))
        endif
c--------------------------------------------------------------------
      call findop('DLAP' ,noc,loc)
      if (noc.gt.0)
     .  then
            aux =  rfldata(loc)
            if(aux.gt.1.d-10)
     .        then
                  DLAP = aux
                  ILAP = 1
                  write(*,*)'SHIFTING VERTICES...',DLAP
              else
                  DLAP = 0.d0
                  ILAP = 0
              endif
        endif
c--------------------------------------------------------------------
      call findop('ZOO',noc,iloc)
      if (noc.gt.0) then
         xmi=rfldata(iloc)
         xma=rfldata(iloc+1)
         ymi=rfldata(iloc+2)
         yma=rfldata(iloc+3)
         wmi=rfldata(iloc+4)
         wma=rfldata(iloc+5)
         write(*,*)'ZOOMING WITH: ',xmi,xma,ymi,yma,wmi,wma
         izoo = 1

c     ...give up zooming:
         if((xma-xmi).le.1.d-10 .or. (yma-ymi).le.1.d-10 .or.
     .      (wma-wmi).le.1.d-10) then
           write(*,*)'GIVE UP ZOOMING'
           izoo=0
         endif
      endif
c---------------------------------------------------------------------
c     ...give up zooming:
      call findop('NOZOO',noc,iloc)
      if (noc.gt.0) then
         write(*,*)'QUIT ZOOMING'
         izoo = 0
      endif
c---------------------------------------------------------------------
      call findop('NRSOL' ,noc,loc)
      if (noc.gt.0)
     .  then
            NRSOLUT = int(rfldata(loc))
        endif
c--------------------------------------------------------------------
      call findop('NOCUTS' ,noc,loc)
      if (noc.gt.0)
     .  then
            NOCUTS = int(rfldata(loc))
        endif
c--------------------------------------------------------------------
      call findop('TOTAL' ,noc,loc)
      if (noc.gt.0)
     .  then
            ITOTAL = int(rfldata(loc))
        endif
c--------------------------------------------------------------------
      call findop('SPH' ,noc,loc)
      if (noc.gt.0)
     .  then
            nsph = int(rfldata(loc))
        endif
c--------------------------------------------------------------------
      call findop('END' ,noc,loc)
      if (noc.gt.0) return


      ISELECT = no_comp*1000+ichoose

      write(*,*) '...PLEASE WAIT, PREPARING IMAGE'
c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


c
c  ...create global image
      call cartobs
c
      write(*,*)'call lsvisidb...'
      call lsvisidbEM(numlev,iflagn)
      XCIM(1) = XCENTR(1)
      XCIM(2) = XCENTR(2)
c
c  ...open the window
      call selwin(Iwind)
c
c  ...decide whether to scale by x or y axes
      q = DIMOB(1)/DIMOB(2)
      q1 = xlength/ylength
      if (q.gt.q1) then
        DIMIM = DIMOB(1)
        SIZE = xlength/2.d0
      else
        DIMIM = DIMOB(2)
        SIZE = ylength/2.d0
      endif
      XCWIN(1) = rmargin + xlength/2.d0
      XCWIN(2) = rmargin + ylength/2.d0
c
c  ...display image
      write(*,*)'call dpvisid...iblno=',iblno
c     call dpvisid(numlev*10+0)
      call dpvisid(numlev*10+iblno)



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

      go to 12


c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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
      if (numlev.gt.0) write(*,*) '11 - CHANGE THE NUMBER OF LEVELS
     . FOR SOLUTION PLOT'
      write(*,*) '12 - SELECT DOMAINS TO DISPLAY'
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
c
c  ...make some elements invisible...
      case(7)
        write(*,*) 'DISPLAY ELEMENT NUMBERS ? (0-NO/1-YES)'
        read(*,*) iblno
        if (iblno.eq.1) then
          call selwin(iwind)
          call dpvisid(numlev*10+iblno)
          write(*,*) 'PLEASE CLICK THE MOUSE INSIDE THE GRAPHICS'
          write(*,*) 'WINDOW TO CONTINUE ...'
          call closwind(1)
        endif


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
          if (mdle.ne.0) then
            NRINVBL = NRINVBL + 1
            IGINV(NRINVBL) = mdle
          endif
          mdlep=mdle
        enddo
c
c  ...leave only a few elements visible...
      case(8)
        write(*,*) 'DISPLAY ELEMENT NUMBERS ? (0-NO/1-YES)'
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
        write(*,*)'SET NEW DEGREE OF COARSENESS FOR OVALS.'
        write(*,*)'MUST BE BETWEEN 1 AND 4'
        write(*,*)'CURRENT VALUE IS ', nrsub1
        read(*,*) nrsub1
        if (nrsub1.gt.4) then
          write(*,*)'HAVE CHANGED TO 4'
          nrsub1 = 4
        else if (NRSUB.lt.1) then
          write(*,*)'HAVE CHANGED TO 1'
          nrsub1 = 1
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
        write(*,*) 'SELECT WHICH DOMAINS TO DISPLAY'
        write(*,*) 'NRDOMAIN = ',NRDOMAIN
        read(*,*) NDOMAIN(1:NRDOMAIN)
        write(*,7038) NDOMAIN(1:NRDOMAIN)
 7038   format('object3: NDOMAIN = ',20i2)
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
        write(*,*) 'DISPLAY ELEMENT NUMBERS ? (0/1) '
        read(*,*) iblno
      case(1)
        write(*,*) 'DISPLAY NODE NUMBERS ? (0/2) '
        read(*,*) iblno
      case(2)
        write(*,*) 'DISPLAY ELEMENT NUMBERS ? (0/1) '
        read(*,*) iblno
      case(3)
        iblno=0
      end select
      call dpvisid(numlev*10+iblno)
c
      go to 10
c
c
      end

