c//////////////////////////////////////////////////////////////////////
c P. Gatto, Nov 12 : this is a LEGACY routine. Do not use it, unless
c                    you revise it!
c//////////////////////////////////////////////////////////////////////
c
c
c----------------------------------------------------------------------
c   latest revision    - Jul 04
c
c   purpose            - routine reads in geometry data for the
c                        GMP package in the OLD FORMAT, allocates
c                        dynamically the necessary data structure
c                        arrays, and produces a copy of the input
c                        file in the NEW FORMAT
c
c   arguments          - none
c---------------------------------------------------------------------
c
      subroutine input_geometry
c
      use GMP
#include "syscom.blk"
#include "cinout.blk"
c
      character*9 type
      integer label, Nrcurv
      integer  tem(8),ttem(12),tttem(6)
c
c----------------------------------------------------------------------
c
      iprint=0
c
c---------------------------------------------------------------------
c
      iprint=1
c
c  ...read in the dimension of the problem and the manifold
      read(NIN,*) NDIM,MANDIM
c
      if ((NDIM.ne.2).and.(NDIM.ne.3)) then
        write(*,*) 'input_geometry: ERROR !!'
        write(*,*) 'WRONG DIMENSION OF THE PROBLEM. MUST BE 2 OR 3 !'
        write(*,*) 'NDIM = ',NDIM
        stop 1
      endif
c
      if ((MANDIM.ne.2).and.(MANDIM.ne.3)) then
        write(*,*) 'input_geometry: ERROR !!'
        write(*,*) 'WRONG DIMENSION OF THE MANIFOLD. MUST BE 2 OR 3 !'
        write(*,*) 'MANDIM = ',MANDIM
        stop 1
      endif
c
      write(NOUT,*) 'INPUT ECHO.............................'
      write(NOUT,*) 'NDIM,MANDIM = ',NDIM,MANDIM
      write(NWE,*)NDIM,MANDIM,'     ...NDIM,MANDIM'
      write(NWE,*)''
      write(NWE,*)''
c
c
c ....................................................................
c
c           dynamically allocate the memory for surfaces
c           and read in geometry data for surfaces
c
c.....................................................................
c
c  ...read in number of surfaces
      read(NIN,*) NRSURFS
      write(NOUT,*) 'NRSURFS = ',NRSURFS
      write(NWE,*)NRSURFS, '                 ...NUMBER OF SURFACES '
      write(NWE,*)''
c
c  ...allocate the memory dynamically
         allocate( SURFACES(NRSURFS), STAT=ic )
         if (ic.ne.0) then
            print*,'SURFACES ARE NOT ALLOCATED '
            stop 1
         endif
c
c  ...read in and store the data on surfaces...
      do isur=1,NRSURFS
c
c  .....transfer label of the points  to type of the points
        read(NIN,*) label
        select case(label)
c
c  .......plane with normal vector and one point
          case(1)
            Type = 'VecPt'
c
c  .......plane with three points
          case(2)
            Type = 'ThrPt'
c
c  .......sphere with radius and center
          case(3)
            Type = 'Sphere'
c
c  .......cylinder with vector of direction ,  point of begining and
c         radius
          case(4)
            Type = 'Cylinder'
c
c  .......ellipsoid with semi-axes a,b,c along x,y,z and given central
c         point
          case(7)
            Type = 'Ellipsoid'
c
c  .......other cases
          case default
            write(*,*)'input_geometry.f: WRONG SURFACE LABEL !!'
            stop 1
c
        end select
        SURFACES(isur)%Type = Type
        write(NWE,*) Type, '                    ...TYPE OF SURFACE',
     .               isur
c
c
        select case(Type)
c
c  .......plane which is normal to a given vector and passing through
c         a point, orientation specified by the normal
          case('VecPt')
            allocate( SURFACES(isur)%Rdata(6), STAT = ic)
            if (ic.ne.0) then
              write(*,*) 'Rdate IS NOT ALLOCATED',
     .                   'SURFACE NUMBER = ',isur
              stop 1
            endif
            read(NIN,*)(SURFACES(isur)%Rdata(j),j=1,3)
            read(NIN,*)(SURFACES(isur)%Rdata(3+j),j=1,3)
            write(NWE,300)(SURFACES(isur)%Rdata(j),j=1,3),
     .           '   ...COORDINATES OF THE POINT'
            write(NWE,300)(SURFACES(isur)%Rdata(3+j),j=1,3),
     .           '   ...COORDINATES OF THE NORMAL VECTOR'
c
c  .......plane passing through three given points, A,B,C
c         orientation specified by AB x AC
          case('ThrPt')
            allocate( SURFACES(isur)%Rdata(9), STAT = ic)
            if(ic.ne.0)then
              write(*,*) 'Rdate IS  NOT ALLOCATED',
     .                   'SURFACE NUMBER = ',isur
              stop 1
            endif
            read(NIN,*)( SURFACES(isur)%Rdata(j),j=1,3)
            read(NIN,*)(SURFACES(isur)%Rdata(3+j),j=1,3)
            read(NIN,*)(SURFACES(isur)%Rdata(6+j),j=1,3)
            write(NWE,300)( SURFACES(isur)%Rdata(j),j=1,3),
     .           '   ...COORDINATES OF THE FIRST POINT'
            write(NWE,300)(SURFACES(isur)%Rdata(3+j),j=1,3),
     .           '   ...COORDINATES OF THE SECOND POINT'
            write(NWE,300)(SURFACES(isur)%Rdata(6+j),j=1,3),
     .           '   ...COORDINATES OF THE THIRD POINT'
c
c  .......sphere
          case('Sphere')
            allocate( SURFACES(isur)%Rdata(4), STAT = ic)
            if(ic.ne.0)then
              write(*,*) 'SURFACES COORDINATES ARE NOT ALLOCATED',
     .                   'SURFACE NUMBER = ',isur
              stop 1
            endif
            read(NIN,*)(SURFACES(isur)%Rdata(j),j=1,3)
            read(NIN,*)SURFACES(isur)%Rdata(4)
            write(NWE,300)(SURFACES(isur)%Rdata(j),j=1,3),
     .           '   ...COORDINATES OF THE CENTER'
            write(NWE,100)SURFACES(isur)%Rdata(4),
     .           '                     ...RADIUS OF THE SPHERE'
c
c  .......cylinder
          case('Cylinder')
            allocate( SURFACES(isur)%Rdata(7), STAT = ic)
            if(ic.ne.0)then
              write(*,*) 'Rdata IS  NOT ALLOCATED',
     .                   'SURFACE NUMBER = ',isur
              stop 1
            endif
            read(NIN,*)(SURFACES(isur)%Rdata(j),j=1,3)
            read(NIN,*)(SURFACES(isur)%Rdata(3+j),j=1,3)
            read(NIN,*) SURFACES(isur)%Rdata(7)
            write(NWE,300)(SURFACES(isur)%Rdata(j),j=1,3)
            write(NWE,300)(SURFACES(isur)%Rdata(3+j),j=1,3)
            write(NWE,*) SURFACES(isur)%Rdata(7)
c
c  .......ellipsoid with semi-axes a,b,c along x,y,z and given central
c         point
          case('Ellipsoid')
            allocate( SURFACES(isur)%Rdata(6), STAT = ic)
            if(ic.ne.0)then
              write(*,*) 'Rdata IS  NOT ALLOCATED',
     .                   'SURFACE NUMBER = ',isur
              stop 1
            endif
c
c  .........read coordinates of the central point
            read(NIN,*)(SURFACES(isur)%Rdata(j),j=1,3)
c
c  .........read semiaxes lengths
            read(NIN,*)(SURFACES(isur)%Rdata(3+j),j=1,3)
c
c  .........write out the data to the new format input file
            write(NWE,300)(SURFACES(isur)%Rdata(j),j=1,3)
            write(NWE,300)(SURFACES(isur)%Rdata(3+j),j=1,3)
c
c  .......other cases
          case default
            write(*,*)'The wrong suface type'
        end select
c
        write(NWE,*)''
      enddo
c
c
c ......................................................................
c
c           dynamically allocate the memory for points
c           and read in geometry data for points
c
c.......................................................................
c
c  ...read in number of points
      read(NIN,*) NRPOINT
      write(NOUT,*) 'NRPOINT = ',NRPOINT
      write(NWE,*)''
      write(NWE,*)''
      write(NWE,*)NRPOINT, '                 ...NUMBER OF POINTS'
      write(NWE,*)''
c
c  ...allocate the memory dynamically
      allocate( POINTS(NRPOINT), STAT=ic )
      if( ic.ne.0 )then
        print*,'POINTS ARE NOT ALLOCATED'
        stop
      endif
c
c  ...read in and store the points data
c
c
      do ip=1,NRPOINT
c
c  .....read in the label and transfer to point type
        read(NIN,*) label
        select case(label)
        case(1)
          Type = 'Regular'
        case(50)
          Type = 'Implicit'
        case(30)
          Type = 'CoorNrm'
        case default
          write(*,*)'THE WRONG LABEL'
        end select
        POINTS(ip)%Type  = Type
        write(NWE,*)Type, '                    ...TYPE OF POINT  ',ip
c
c  .....read in and store the number of curves
        read(NIN,*) NrCurv
        POINTS(ip)%NrCurv = NrCurv
        write(NOUT,*)'POINT NO IS', IP
        write(NOUT,*) 'NUMBER OF CURVES MEETING',
     .               ' AT THE POINT = ',NrCurv
c
c  .....read in and store the curves' numbers
        allocate( POINTS(ip)%CurvNo(NrCurv), STAT = ic)
        if (ic.ne.0) then
           write(*,*) 'CurvNo IS  NOT ALLOCATED',
     .          'POINT NUMBER IS', ip
           stop 1
        endif
        if(nrcurv.gt.0)then
           read(NIN,*) (POINTS(ip)%CurvNo(j),j=1,NrCurv)
           write(NOUT,*) 'THE CURVES NUMBERS = ',
     .          (POINTS(ip)%CurvNo(j),j=1,NrCurv)
        endif
c
c  .....read in and store the coordinates for regular points
        select case(Type)
        case('Regular')
          allocate( POINTS(ip)%Rdata(NDIM), STAT = ic)
          if (ic.ne.0) then
            write(*,*) 'Rdate IS  NOT ALLOCATED',
     .                 'POINT NUMBER IS', ip
            stop 1
          endif
          read(NIN,*) (POINTS(ip)%Rdata(k),k=1,NDIM)
          write(NOUT,*) 'COORDINATES OF THE POINT IS ',
     .                  (POINTS(ip)%Rdata(k),k=1,NDIM)
          select case(NDIM)
          case(2)
          write(NWE,3000)(POINTS(ip)%Rdata(k),k=1,NDIM),
     .         '   ...COORDINATES OF THE POINT'
          case(3)
          write(NWE,300)(POINTS(ip)%Rdata(k),k=1,NDIM),
     .         '   ...COORDINATES OF THE POINT'
          end select
c
        case('CoorNrm')
          allocate( POINTS(ip)%Rdata(6), STAT = ic)
          if (ic.ne.0) then
            write(*,*) 'Rdate IS  NOT ALLOCATED',
     .                 'POINT NUMBER IS', ip
            stop 1
          endif
          read(NIN,*) (POINTS(ip)%Rdata(k),k=1,6)
          write(NOUT,*) 'COORDINATES OF THE POINT IS ',
     .         (POINTS(ip)%Rdata(k),k=1,NDIM)
          write(NWE,600)(POINTS(ip)%Rdata(k),k=1,6),
     .          '   ...COORDINATES AND NORMAL VECTORS OF THE POINT'
c
c  .......read in and store the intersecting surface and
c         starting point for implicit  points
        case('Implicit')
          if (NDIM.ne.3) then
            write(*,*) 'input_geometry: IMPLICIT POINT NOT ALLOWED',
     .                 ' FOR 2D PROBLEMS'
            stop 1
          endif
          allocate( POINTS(ip)%Idata(3), STAT = ic)
          if (ic.ne.0) then
            write(*,*) 'Idata IS  NOT ALLOCATED',
     .                 'POINT NUMBER IS', ip
            stop 1
          endif
          read(NIN,*) (POINTS(ip)%Idata(j),j=1,3)
          write(NWE,30) 'SURFACES: ',(POINTS(ip)%Idata(j),j=1,3)
          write(NOUT,*) 'SURFACES: ',(POINTS(ip)%Idata(j),j=1,3)
          allocate( POINTS(ip)%Rdata(3), STAT = ic)
          if (ic.ne.0) then
            write(*,*) 'Rdata IS  NOT ALLOCATED',
     .                 'POINT NUMBER IS', ip
            stop 1
          endif
          read(NIN,*) (POINTS(ip)%Rdata(j),j=1,3)
          select case(NDIM)
          case(2)
          write(NWE,3000) (POINTS(ip)%Rdata(j),j=1,3),
     .         '   ...INITIAL POINT'
          case(3)
          write(NWE,300) (POINTS(ip)%Rdata(j),j=1,3),
     .         '   ...INITIAL POINT'
          end select
        case default
          write(*,*)'WRONG TYPE FOR THE POINT'
          stop 1
        end select
        write(NOUT,*)''
        write(NWE,*)''
c
      enddo
c
c
c .........................................................
c
c           dynamically allocate the memory for curves
c           and read in geometry data for curves
c
c..........................................................

c  ...Read in number of curves
      read(NIN,*) NRCURVE
      write(NOUT,*) 'NRCURVE = ',NRCURVE
      write(NWE,*)''
      write(NWE,*)NRCURVE, '                 ...NUMBER OF CURVES'
      write(NWE,*)''
c
c  ...Allocate the memory dynamically
      allocate( CURVES(NRCURVE), STAT=ic )
      if( ic.ne.0 )then
        print*,'CURVES ARE NOT ALLOCATED'
        stop
      endif
c
c  ...Read in and store the curves data
c
c
      do icur=1,NRCURVE
c
c  .....read in and store the curve type
        read(NIN,*) label
        select case(label)
          case(1)
            Type = 'Seglin'
          case(-1)
            Type = 'QuaCir'
          case(-2)
            Type = 'SegCir'
          case(-3)
            Type = 'QuaEl1'
          case(-4)
            Type = 'QuaEl2'
          case(-5)
            Type = 'QuaSEl'
          case(30)
             Type = 'HermCur'
          case(50)
            Type = 'ImpCir'
          case default
            write(*,*) 'input_geometry: WRONG CURVE TYPE !!'
            stop 1
        end select
        CURVES(icur)%Type  = Type
        write(NWE,*)Type, '                    ...TYPE OF CURVE',icur
c
c  .....read in and store endpoints for the curve
        read(NIN,*) (CURVES(icur)%EndPoNo(j),j=1,2)
        write(NOUT,*)'CURVE NUMBER IS',ICUR
        write(NOUT,*) 'ENDPOINTS OF THE CURVE IS= ',
     .                (CURVES(icur)%EndPoNo(j),j=1,2)
        write(NWE,20)(CURVES(icur)%EndPoNo(j),j=1,2),
     .       '                      ...ENDPOINTS OF THE CURVE'
c  ...read and store number of figures meeting along the curve
        read(NIN,*) NrFig
        CURVES(icur)%NrFig = NrFig
        write(NOUT,*) 'NUMBER OF FIGURES MEETING ALONG THE CURVE = ',
     .                 NrFig
c
c  .....read in and store the figures numbers
        if (NrFig.gt.0) then
           allocate( CURVES(icur)%FigNo(NrFig), STAT = ic)
           if(ic.ne.0)then
              write(*,*) 'FigNo IS NOT ALLOCATED',
     .             'THE CURVE NUMBER IS ', icur
              stop 1
           endif
           read(NIN,*) (CURVES(icur)%FigNo(j),j=1,NrFig)
           write(NOUT,*) 'THE FIGURES NICKNAMES =',
     .          (CURVES(icur)%FigNo(j),j=1,nrfig)
        else
           read(NIN,*)
        endif
c
c  .....read in and store the coordinates of the center of a circle
        if ((Type.eq.'QuaCir').or.(Type.eq.'SegCir')) then
          allocate( CURVES(icur)%Rdata(NDIM), STAT = ic)
          if(ic.ne.0)then
            write(*,*) 'Rdata IS NOT ALLOCATED',
     .                 'THE CURVE NUMBER IS ', icur
            stop 1
          endif
          read(NIN,*) (CURVES(icur)%Rdata(k),k=1,NDIM)
          write(NOUT,*) 'COORDINATES: ',
     .                  (CURVES(icur)%Rdata(k),k=1,NDIM)
          select case(NDIM)
          case(2)
          write(NWE,3000) (CURVES(icur)%Rdata(k),k=1,NDIM),
     .         '   ...CENTER OF THE CURVE'
          case(3)
          write(NWE,300) (CURVES(icur)%Rdata(k),k=1,NDIM),
     .         '   ...CENTER OF THE CURVE'
          end select
c
c  .....read in and store the coordinates of the center of an ellipse
        else if ((Type.eq.'QuaEl1').or.(Type.eq.'QuaEl2')) then
          allocate( CURVES(icur)%Rdata(NDIM), STAT = ic)
          if(ic.ne.0)then
            write(*,*) 'Rdata IS NOT ALLOCATED',
     .                 'THE CURVE NUMBER IS ', icur
            stop 1
          endif
          read(NIN,*) (CURVES(icur)%Rdata(k),k=1,NDIM)
          write(NOUT,*) 'COORDINATES: ',
     .                  (CURVES(icur)%Rdata(k),k=1,NDIM)
          select case(NDIM)
          case(2)
          write(NWE,3000) (CURVES(icur)%Rdata(k),k=1,NDIM),
     .         '   ...CENTER OF THE CURVE'
          case(3)
          write(NWE,300) (CURVES(icur)%Rdata(k),k=1,NDIM),
     .         '   ...CENTER OF THE CURVE'
          end select
c
c  .....read in and store the coordinates of the center of a
c  .....superellipse and its associated powers          
        else if (Type.eq.'QuaSEl') then
          allocate( CURVES(icur)%Rdata(NDIM+2), STAT = ic)
          if(ic.ne.0)then
            write(*,*) 'Rdata IS NOT ALLOCATED',
     .                 'THE CURVE NUMBER IS ', icur
            stop 1
          endif
          read(NIN,*) (CURVES(icur)%Rdata(k),k=1,NDIM+2)
          write(NOUT,*) 'COORDINATES: ',
     .                  (CURVES(icur)%Rdata(k),k=1,NDIM)
          write(NOUT,*) 'POWERS OF SUPERELLIPSE: ',
     .                  (CURVES(icur)%Rdata(k),k=NDIM+1,NDIM+2)
          select case(NDIM)
          case(2)
          write(NWE,3000) (CURVES(icur)%Rdata(k),k=1,NDIM),
     .         '   ...CENTER OF THE CURVE'
          case(3)
          write(NWE,300) (CURVES(icur)%Rdata(k),k=1,NDIM),
     .         '   ...CENTER OF THE CURVE'
          end select
          write(NWE,3000) (CURVES(icur)%Rdata(k),k=NDIM+1,NDIM+2),
     .         '   ...POWERS OF SUPERELLIPSE'
c
c  .....read in and store the intersecting surface for implicit curves
        else if (Type.eq.'ImpCir')then
          allocate( CURVES(icur)%Idata(4), STAT = ic)
          if (ic.ne.0)then
            write(*,*) 'Idata IS NOT ALLOCATED',
     .                 'THE CURVE NUMBER IS ', icur
            stop 1
          endif
          read(NIN,*) (CURVES(icur)%Idata(k),k=1,4)
c
c  .......ignore the signs that might be present in some old
c         input files
          do k=1,4
            CURVES(icur)%Idata(k) = iabs(CURVES(icur)%Idata(k))
          enddo
          write(NWE,40) (CURVES(icur)%Idata(k),k=1,4),
     .         '              ...SURFACE NUMBERS CONSTITUTING THE CURVE'
c
       else if (Type.eq.'HermCur')then
          allocate( CURVES(icur)%Rdata(12), STAT = ic)
          if(ic.ne.0)then
             write(*,*) 'Rdata IS NOT ALLOCATED',
     .            'THE CURVE NUMBER IS ', icur
             stop 1
          endif
          read(NIN,*)CURVES(icur)%Rdata(1:12)
          write(NWE,*)CURVES(icur)%Rdata(1:12)
          write(NOUT,*)CURVES(icur)%Rdata(1:12)
c
        endif
c
        write(NOUT,*)''
        write(NWE,*)''
      enddo
      write(NWE,*)''
c
c
c
c .........................................................
c
c           dynamically allocate the memory for triangles
c           and read in geometry data for triangles
c
c..........................................................

c  ...Read in number of triangles
      read(NIN,*) NRTRIAN
      write(NOUT,*) 'NRTRIAN = ',NRTRIAN
      write(NWE,*)NRTRIAN, '                 ...NUMBER OF TRIANGLES'
c
c  ...Allocate the memory dynamically
      allocate( TRIANGLES(NRTRIAN), STAT=ic )
      if( ic.ne.0 )then
        print*,'TRIANGLES ARE NOT ALLOCATED'
        stop
      endif
c
c  ...Read in and store the triangles data
c
c
      do itri=1,NRTRIAN
c  ...read the labers and transfer to the triangle types
        read(NIN,*) label
        select case(label)
          case(1)
            Type = 'PlaneTri'
          case(2)
            Type = 'TransTri'
          case(-1)
            Type = 'SpherTri'
          case(-2)
            Type = 'QtCirTri'
          case(-3)
            Type = 'PaSphTri'
          case(50)
            Type = 'ImpliTri'
          case(-4)
            Type = 'ImSphTri'
          case default
            write(*,*) 'input_geometry: WRONG TRIANGLE LABEL!'
          end select
        TRIANGLES(itri)%Type  = Type
         write(NWE,*)Type, '        ...TYPE OF THE TRIANGLE',itri
c
c  ...read in and store the edges' numbers that consistute
c     the triangle
        read(NIN,*) (TRIANGLES(itri)%EdgeNo(j),j=1,3)
        write(NOUT,*) 'THE EDGE NICKNAMES = ',
     .                 (TRIANGLES(itri)%EdgeNo(j),j=1,3)
c
         do i =1,3
            nc = TRIANGLES(itri)%EdgeNo(i)
            if(nc.gt.0)then
               tem(i) = CURVES(nc)%EndPoNo(1)
            else
               tem(i) = CURVES(-nc)%EndPoNo(2)
            endif
         enddo
         write(NWE,30)tem(1:3),'   ...THREE VERTEXES OF THE TRIANGLE'
c
c  ...read and store nicknames of blocks adjacent to the triangle
        if (MANDIM.eq.3) then
          read(NIN,*) (TRIANGLES(itri)%BlockNo(k),k=1,2)
          write(NOUT,*) 'ADJACENT BLOCKS NUMBERS = ',
     .                   (TRIANGLES(itri)%BlockNo(k),k=1,2)
        endif
c
c
c  ...read in and store the surface that define the implicit  triangles
        if (type.eq.'ImpliTri')then
          allocate(TRIANGLES(itri)%Idata(4), STAT = ic)
          if(ic.ne.0)then
             write(*,*) 'Idata IS NOT ALLOCATED',
     .                 'THE TRIANGLE NUMBER IS ', itri
            stop 1
          endif
          read(NIN,*) (TRIANGLES(itri)%Idata(j),j=1,4)
          write(NOUT,*) 'TRIANGLES%Idata: ',
     .                (TRIANGLES(itri)%Idata(j),j=1,4)
          write(NWE,40)(TRIANGLES(itri)%Idata(j),j=1,4)
        endif
c
      enddo
c
c
c .........................................................
c
c           dynamically allocate the memory for rectangles
c           and read in geometry data for rectangles
c
c..........................................................

c  ...Read in number of rectangles
      read(NIN,*) NRRECTA
      write(NOUT,*) 'NRRECTA = ',NRRECTA
      write(NWE,*)''
      write(NWE,*)''
      write(NWE,*)NRRECTA,'                 ...NUMBER OF RECTANGLES'
      write(NWE,*)''
c
c  ...Allocate the memory dynamically
      allocate( RECTANGLES(NRRECTA), STAT=ic )
      if( ic.ne.0 )then
         print*,'RECTANGLES ARE NOT ALLOCATED'
         stop 1
      endif
c
c  ...Read in and store the rectangles data
c
c
      do irec=1,NRRECTA
        if (irec.eq.14) then
          iprint=1
        else
          iprint=0
        endif
c
c  .....read in and store the type
        read(NIN,*) label
        select case(label)
          case(1)
            Type = 'BilQua'
          case(2)
            Type = 'TraQua'
          case(3)
            Type = 'TraQua'
          case(-1)
            Type = 'CylRec'
          case(40)
            Type = 'HermRec'
          case(50)
            Type = 'ImpRec'
          case default
            write(*,*)'A WRONG RECTANGLE TYPE'
        end select
        RECTANGLES(irec)%Type  = Type
        write(NWE,*)Type, '                    ...TYPE OF RECTANGLE'
     .             ,irec
c
c  .....read in and store the edge curves numbers that constitute
c       the rectangle
        read(NIN,*) (RECTANGLES(irec)%EdgeNo(j),j=1,4)
        write(NOUT,*)'RECTANGULE NUMBER IS',IREC
        write(NOUT,*) 'THE EDGE NICKNAMES = ',
     .               (RECTANGLES(irec)%EdgeNo(j),j=1,4)
        tem = 0
        do i =1,4
          nc = RECTANGLES(irec)%EdgeNo(i)
          if (nc.ne.0) then
            if (nc.gt.0) then
              tem(i) = CURVES(nc)%EndPoNo(1)
            else
              tem(i) = CURVES(-nc)%EndPoNo(2)
            endif
          endif
        enddo
        write(NWE,40) tem(1:4),
     .        '              ...FOUR VERTICES OF THE RECTANGLE'
c
c  .....read and store nicknames of blocks adjacent to the rectangle
        if (MANDIM.eq.3) then
          read(NIN,*) (RECTANGLES(irec)%BlockNo(k), k=1,2)
          write(NOUT,*) 'ADJACENT BLOCKS NUMBERS = ',
     .                  (RECTANGLES(irec)%BlockNo(k), k=1,2)
        endif
c
c
c  .....read in and store surfaces that define the implicit rectangle
        if (Type.eq.'ImpRec') then
c
          allocate(RECTANGLES(irec)%Idata(5), STAT = ic)
          if (ic.ne.0) then
            write(*,*) 'Idata IS NOT ALLOCATED',
     .                 'RECTANGLE NUMBER IS', irec
            stop 1
          endif
          call read_imprec(irec)
c
          write(NOUT,*) 'IDATA OF THE RECTANGLE IS',
     .                 (RECTANGLES(irec)%Idata(j),j=1,5)
          write(NWE,50) (RECTANGLES(irec)%Idata(j),j=1,5),
     .        '          ...SURFACE NUMBERS CONSTITUTING THE RECTANGLE'

        elseif (type.eq.'HermRec')then
           allocate(RECTANGLES(irec)%Rdata(24), STAT = ic)
           if (ic.ne.0)then
              write(*,*) 'Idata IS NOT ALLOCATED',
     .                   'RECTANGLE NUMBER IS', irec
              stop 1
           endif
           read(NIN,*) (RECTANGLES(irec)%Rdata(j),j=1,24)
           write(NOUT,*) 'RDATA OF THE RECTANGLE IS',
     .                  (RECTANGLES(irec)%Rdata(j),j=1,24)
           write(NWE,2400) (RECTANGLES(irec)%Rdata(j),j=1,24),
     .          '          ...COEFFICIENTS OF THE RECTANGLE'

        endif

c
        write(NOUT,*)''
        write(NWE,*)''
      enddo
c
c
c .........................................................
c
c           dynamically allocate the memory for prisms
c           and read in geometry data for prisms
c
c..........................................................
      if (MANDIM.eq.3) then
c
c  .....read in number of prisms
        read(NIN,*) NRPRISM
        write(NOUT,*) 'NRPRISM = ',NRPRISM
       write(NWE,*)''
       write(NWE,*) NRPRISM,   '                 ...NUMBER OF PRISMS'
       write(NWE,*)''
c
c  .....allocate the memory dynamically
        allocate( PRISMS(NRPRISM), STAT=ic )
        if( ic.ne.0 )then
           print*,'PRISMS ARE NOT ALLOCATED'
           stop 1
        endif
c
c  .....read in and store the prisms data
c
        do ipri=1,NRPRISM
c
c  .......read in and store the type
          read(NIN,*) label
          select case(label)
          case(1)
            Type = 'LiLiPrism'
          case(-1)
            Type = 'TriaShell'
          case default
            write(*,*)'THE WRONG LABEL FOR PRISMS'
          end select
          PRISMS(ipri)%Type  = Type
c
          write(NWE,*) Type, '                       ...TYPE OF THE ',
     .         'PRISM',ipri
c  .......read in and store the figures' numbers that consistute
c         the prism
          read(NIN,*) (PRISMS(ipri)%FigNo(j), j = 1,5)
        enddo
c
c
c.......................................................................
c
c           dynamically allocate the memory for hexas
c           and read in geometry data for hexas
c
c......................................................................
c
c  .....read in the number of hexas
        read(NIN,*) NRHEXAS
        write(NOUT,*) 'NRHEXAS = ',NRHEXAS
        write(NWE,*)NRHEXAS,'                 ...NUMBER OF HEXAHEDRONS'
        write(NWE,*)''
c
c  .....allocate the memory dynamically
        allocate( HEXAS(NRHEXAS), STAT=ic )
        if( ic.ne.0 )then
          print*,'HEXAHEDRAL ARE NOT ALLOCATED'
          stop 1
        endif
c
c  .....read in and store the hexas data
        do ihex=1,NRHEXAS
c
c  .......read in and store the hexa type
          read(NIN,*) label
          select case(label)
          case(1)
            Type = 'TriLiHex'
          case(-1)
            Type = 'RecShell'
          case default
            write(*,*)'THE WRONG LABEL FOR HEXAS'
          end select
          HEXAS(ihex)%Type  = Type
c
c  .......read in and store the figures' numbers that constitute
c         the hexahedron faces
          read(NIN,*) (HEXAS(ihex)%FigNo(j),j=1,6)
          allocate( HEXAS(ihex)%Idata(1), STAT=ic )
          if( ic.ne.0 )then
             print*,'HEXAHEDRAL ARE NOT ALLOCATED'
             stop 1
          endif
          write(NOUT,*)'HEXA NUMBER IS',IHEX
          write(NOUT,*)'FACES NUMBER ARE',(HEXAS(ihex)%FigNo(j),j=1,6)
          write(NOUT,*)''
c
        enddo
c
c----------------------------------------------------------------------
c
c  .....determine all orientations and extra connectivities
        call orientations
c
c----------------------------------------------------------------------
c
c  .....output the data on hexahedra
c
c  .....read in and store the hexas data
        do ihex=1,NRHEXAS
c
          Type = HEXAS(ihex)%Type
          write(NWE,*) Type, '                    ...TYPE OF THE ',
     .         'HEXAHEDRON',ihex
c
c  .......output the vertex points connectivities defining the
c         hexa in the new format
          write(NWE,80) (HEXAS(ihex)%VertNo(j),j=1,8),
     .         '      ...EIGHT VERTICES OF THE HEXAHEDRON'
          write(NWE,*)''
c
        enddo
c
c  ...if a 3D manifold
      endif



c
c
 100  format(1f16.10,5a)
 300  format(3f16.10,5a)
 3000 format(2f16.10,5a)
 600  format(6f16.10,5a)
 2400  format(24f16.10,3a)
 20   format(2I6,5a)
 30   format(3I6,5a)
 40   format(4I6,5a)
 50   format(5I6,5a)
 80   format(8I6,5a)
c
      return
      end
c
c
c
c----------------------------------------------------------------------
c
c   routine name       - orientations
c
c---------------------------------------------------------------------
c
c
c   computer           - machine independent
c
c   latest revision    - Jul 04
c
c   purpose            - routine determines orientations:
c                          figure to edges
c                          hexa to faces, hexa to edges
c                        and stores them in the GMP data structure
c                        arrays
c
c   arguments          - none
c
c   required  routines -
c
c---------------------------------------------------------------------
c
      subroutine orientations
c
      use GMP
#include "syscom.blk"
c
c  ...face rectangles numbers
      dimension nof(6)
c
c  ...edge curves numbers
      dimension noedg(12)
c
c  ...for each edge iedg=1,12
c        nedg_to_face(1,1:2,iedg) contains numbers of faces that
c        share the edge
c        nedg_to_face(2,1:2,iedg) contains the corresponding edge
c        number in hexa local face coordinates
      dimension nedg_to_face(2,2,12)
c
c  ...for each face nsid=1,6 and its edge j=1,4 in face loacl numbering
c        naxes(j,nsid) = edge number in face RECTANGLE coordinates
      dimension naxes(4,6)
c
c  ...face orientations
      dimension norient(6)
c
c  ...for each vertex point i=1,4 of a rectangle,
c        nvert_to_edge(1:2,i) contains numbesr of adjacent edges
      dimension nvert_to_edge(2,4)
c
c  ...sign factors
      dimension nfactor(0:7)
c
      data nedg_to_face /1,1,3,1, 1,2,4,1, 1,3,5,1, 1,4,6,1,
     .                   2,1,3,3, 2,2,4,3, 2,3,5,3, 2,4,6,3,
     .                   3,4,6,4, 3,2,4,4, 4,2,5,2, 5,4,6,2/
c
      data nvert_to_edge /4,1, 1,2, 2,3, 3,4/
c
      data nfactor/1,1,1,1, -1,-1,-1,-1/
c
c----------------------------------------------------------------------
c
      iprint=0
c
c
c  ...loop through rectangles in the GMP model
      do irec=1,NRRECTA
c
c  .....determine the rectangle edges
        do i=1,4
          noedg(i) = iabs(RECTANGLES(irec)%EdgeNo(i))
        enddo
        if (iprint.eq.1) then
          write(*,7004) irec,(noedg(i),i=1,4)
 7004     format('orientations: irec,noedg = ',i4,2x,4i6)
        endif
        if (noedg(1).eq.0) stop 1
c
c  .....use the first vertex to determine orientations of the
c       first and fourth edge
        nc1 = noedg(1); nc2 = noedg(4)
        do i1=1,2
          nv1 = CURVES(nc1)%EndPoNo(i1)
          do i2=1,2
            nv2 = CURVES(nc2)%EndPoNo(i2)
            if (nv1.eq.nv2) then
              select case(i1*10+i2)
              case(11)
                RECTANGLES(irec)%EdgeNo(1) =   noedg(1)
                RECTANGLES(irec)%EdgeNo(4) = - noedg(4)
              case(12)
                RECTANGLES(irec)%EdgeNo(1) =   noedg(1)
                RECTANGLES(irec)%EdgeNo(4) =   noedg(4)
              case(21)
                RECTANGLES(irec)%EdgeNo(1) = - noedg(1)
                RECTANGLES(irec)%EdgeNo(4) = - noedg(4)
              case(22)
                RECTANGLES(irec)%EdgeNo(1) = - noedg(1)
                RECTANGLES(irec)%EdgeNo(4) =   noedg(4)
              end select
            endif
          enddo
        enddo
c
c  .....use the third vertex to determine orientations of the
c       second and third edge
        nc1 = noedg(2); nc2 = noedg(3)
        do i1=1,2
          nv1 = CURVES(nc1)%EndPoNo(i1)
          do i2=1,2
            nv2 = CURVES(nc2)%EndPoNo(i2)
            if (nv1.eq.nv2) then
              select case(i1*10+i2)
              case(11)
                RECTANGLES(irec)%EdgeNo(2) = - noedg(2)
                RECTANGLES(irec)%EdgeNo(3) =   noedg(3)
              case(12)
                RECTANGLES(irec)%EdgeNo(2) = - noedg(2)
                RECTANGLES(irec)%EdgeNo(3) = - noedg(3)
              case(21)
                RECTANGLES(irec)%EdgeNo(2) =   noedg(2)
                RECTANGLES(irec)%EdgeNo(3) =   noedg(3)
              case(22)
                RECTANGLES(irec)%EdgeNo(2) =   noedg(2)
                RECTANGLES(irec)%EdgeNo(3) = - noedg(3)
              end select
            endif
          enddo
        enddo
c
c  .....determine the vertices
        do i=1,4
          nc = iabs(RECTANGLES(irec)%EdgeNo(i))
          if (RECTANGLES(irec)%EdgeNo(i).gt.0) then
            RECTANGLES(irec)%VertNo(i) = CURVES(nc)%EndPoNo(1)
          else
            RECTANGLES(irec)%VertNo(i) = CURVES(nc)%EndPoNo(2)
          endif
        enddo
c
        if (iprint.eq.1) then
          write(*,7001) irec
 7001     format('orientations: EDGES AND VERTICES FOR RECTANGLE = ',
     .            i6)
          write(*,7002) (RECTANGLES(irec)%EdgeNo(j),j=1,4)
 7002     format(12i6)
          write(*,7002) (RECTANGLES(irec)%VertNo(j),j=1,4)
        endif
c
c  ...end of loop through rectangles
      enddo
c
c----------------------------------------------------------------------
c
c  ...loop through hexahedra in the GMP model
      do ihex=1,NRHEXAS
c
c  .....determine the hexahedron face rectangles numbers
        do is=1,6
          call decode(iabs(HEXAS(ihex)%FigNo(is)), nof(is),ivoid)
        enddo
c
c  .....determine naxes (see above)
        do iedg=1,12
          isid1 = nedg_to_face(1,1,iedg); n1 = nedg_to_face(2,1,iedg)
          isid2 = nedg_to_face(1,2,iedg); n2 = nedg_to_face(2,2,iedg)
          nrect1 = nof(isid1)
          nrect2 = nof(isid2)
c
c  .......find the edge common to both face rectangles
          do i1=1,4
            nc1 = iabs(RECTANGLES(nrect1)%EdgeNo(i1))
            do i2=1,4
              nc2 = iabs(RECTANGLES(nrect2)%EdgeNo(i2))
              if (nc1.eq.nc2) go to 10
            enddo
          enddo
 10       naxes(n1,isid1) = i1; naxes(n2,isid2) = i2
          noedg(iedg) = nc1
c
c  .....end of loop through edges
        enddo
c
c  .....determine the face orientations
        do is=1,6
          nick = naxes(1,is)*10 + naxes(4,is)
          select case(nick)
          case(14)
            HEXAS(ihex)%FigNo(is) = nof(is)*10 + 0
          case(43)
            HEXAS(ihex)%FigNo(is) = nof(is)*10 + 1
          case(32)
            HEXAS(ihex)%FigNo(is) = nof(is)*10 + 2
          case(21)
            HEXAS(ihex)%FigNo(is) = nof(is)*10 + 3
          case(41)
            HEXAS(ihex)%FigNo(is) = nof(is)*10 + 4
          case(12)
            HEXAS(ihex)%FigNo(is) = nof(is)*10 + 5
          case(23)
            HEXAS(ihex)%FigNo(is) = nof(is)*10 + 6
          case(34)
            HEXAS(ihex)%FigNo(is) = nof(is)*10 + 7
          case default
            write(*,*) 'find_hexa_to_face_orient: ERROR !'
            stop 1
          end select
        enddo
c
c  .....determine the edge orientations
c
c  .....use the bottom face to fix the orientations of the first four
c       edges
        call decode(HEXAS(ihex)%FigNo(1), nrect, io)
        nfact = nfactor(io)
c
c  .....loop through the edges
        do i=1,4
c
c  .......loop through the sides of the rectangle
          do j=1,4
            if (iabs(RECTANGLES(nrect)%EdgeNo(j)).eq.Noedg(i)) then
              nsign = RECTANGLES(nrect)%EdgeNo(j)/Noedg(i)
              if (i.gt.2) nsign = -nsign
              HEXAS(ihex)%EdgeNo(i) = Noedg(i)*nsign*nfact
            endif
          enddo
        enddo
c
c  .....use the top face to fix the orientations of the next four
c       edges
        call decode(HEXAS(ihex)%FigNo(2), nrect, io)
        nfact = nfactor(io)
c
c  .....loop through the edges
        do i=5,8
c
c  .......loop through the sides of the rectangle
          do j=1,4
            if (iabs(RECTANGLES(nrect)%EdgeNo(j)).eq.Noedg(i)) then
              nsign = RECTANGLES(nrect)%EdgeNo(j)/Noedg(i)
              if (i.gt.6) nsign = -nsign
              HEXAS(ihex)%EdgeNo(i) = Noedg(i)*nsign*nfact
            endif
          enddo
        enddo
c
c  .....use the front face to fix the orientations of the next two
c       edges
        call decode(HEXAS(ihex)%FigNo(3), nrect, io)
        nfact = nfactor(io)
c
c  .....loop through the edges
        do i=9,10
c
c  .......loop through the sides of the rectangle
          do j=1,4
            if (iabs(RECTANGLES(nrect)%EdgeNo(j)).eq.Noedg(i)) then
              nsign = RECTANGLES(nrect)%EdgeNo(j)/Noedg(i)
              if (i.eq.9) nsign = -nsign
              HEXAS(ihex)%EdgeNo(i) = Noedg(i)*nsign*nfact
            endif
          enddo
        enddo
c
c  .....use the back face to fix the orientations of the last two
c       edges
        call decode(HEXAS(ihex)%FigNo(5), nrect, io)
        nfact = nfactor(io)
c
c  .....loop through the edges
        do i=11,12
c
c  .......loop through the sides of the rectangle
          do j=1,4
            if (iabs(RECTANGLES(nrect)%EdgeNo(j)).eq.Noedg(i)) then
              nsign = RECTANGLES(nrect)%EdgeNo(j)/Noedg(i)
              if (i.eq.12) nsign = -nsign
              HEXAS(ihex)%EdgeNo(i) = Noedg(i)*nsign*nfact
            endif
          enddo
        enddo
c
c  .....determine the vertex points
        nc = iabs(HEXAS(ihex)%EdgeNo(1))
        if (HEXAS(ihex)%EdgeNo(1).gt.0) then
          HEXAS(ihex)%VertNo(1) = CURVES(nc)%EndPoNo(1)
          HEXAS(ihex)%VertNo(2) = CURVES(nc)%EndPoNo(2)
        else
          HEXAS(ihex)%VertNo(2) = CURVES(nc)%EndPoNo(1)
          HEXAS(ihex)%VertNo(1) = CURVES(nc)%EndPoNo(2)
        endif
        nc = iabs(HEXAS(ihex)%EdgeNo(3))
        if (HEXAS(ihex)%EdgeNo(3).gt.0) then
          HEXAS(ihex)%VertNo(4) = CURVES(nc)%EndPoNo(1)
          HEXAS(ihex)%VertNo(3) = CURVES(nc)%EndPoNo(2)
        else
          HEXAS(ihex)%VertNo(3) = CURVES(nc)%EndPoNo(1)
          HEXAS(ihex)%VertNo(4) = CURVES(nc)%EndPoNo(2)
        endif
        nc = iabs(HEXAS(ihex)%EdgeNo(5))
        if (HEXAS(ihex)%EdgeNo(5).gt.0) then
          HEXAS(ihex)%VertNo(5) = CURVES(nc)%EndPoNo(1)
          HEXAS(ihex)%VertNo(6) = CURVES(nc)%EndPoNo(2)
        else
          HEXAS(ihex)%VertNo(6) = CURVES(nc)%EndPoNo(1)
          HEXAS(ihex)%VertNo(5) = CURVES(nc)%EndPoNo(2)
        endif
        nc = iabs(HEXAS(ihex)%EdgeNo(7))
        if (HEXAS(ihex)%EdgeNo(7).gt.0) then
          HEXAS(ihex)%VertNo(8) = CURVES(nc)%EndPoNo(1)
          HEXAS(ihex)%VertNo(7) = CURVES(nc)%EndPoNo(2)
        else
          HEXAS(ihex)%VertNo(7) = CURVES(nc)%EndPoNo(1)
          HEXAS(ihex)%VertNo(8) = CURVES(nc)%EndPoNo(2)
        endif
c
c
        if (iprint.eq.1) then
          write(*,7003) ihex
 7003     format('orientations: FACES, EDGES AND VERTS FOR HEXA = ',i6)
          write(*,7002) (HEXAS(ihex)%FigNo(j),j=1,6)
          write(*,7002) (HEXAS(ihex)%EdgeNo(j),j=1,12)
          write(*,7002) (HEXAS(ihex)%VertNo(j),j=1,8)
        endif
c
c  ...end of loop through the hexahedra
      enddo
c
      return
      end
