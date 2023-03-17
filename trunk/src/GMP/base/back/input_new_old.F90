!-----------------------------------------------------------------------
!
!   routine name       - new_old
!
!-----------------------------------------------------------------------
!
!   latest revision    - Aug 04
!
!   purpose            - routine reads in geometry data for GMP
!                        in the NEW FORMAT, allocates dynamically
!                        the necessary data structure arrays, and
!                        produces automatically a copy of the input
!                        file in the OLD FORMAT
!
!-----------------------------------------------------------------------
!
      subroutine new_old
!
      use GMP
#include "syscom.blk"
#include "cinout.blk"
!
      character(10) type
!
!-----------------------------------------------------------------------
!
      iprint=1
!
!  ...read in the dimension of the problem and the manifold
      read(NWE,*)  NDIM,MANDIM
      write(NIN,*) NDIM,MANDIM,'             ...NDIM,MANDIM'
      write(NIN,*) ''
      write(NIN,*) ''
!
!  ...read in surfaces data
      read(NWE,*)  NRSURFS
      write(NIN,*) NRSURFS,   &
           '                         ...NUMBER OF SURFACES '
      write(NIN,*)''
!
!
!  ...allocate dynamically the memory
      allocate( SURFACES(NRSURFS), STAT=ic )
      if (ic.ne.0) then
        print*,'SURFACES ARE NOT ALLOCATED '
        stop
      endif
      if (iprint.eq.1) then
        write(*,*) 'new_old: READING SURFACES...'
      endif
!
!  ...read in and store the data on surfaces...
      do isur=1,NRSURFS
!
!  .....transfer label of the points  to type of the points
        read(NWE,*) Type
        select case(Type)
!
!    ...plane with normal vector and one point
        case('VecPt')
          label= 1
!
!    ...plane with three points
        case('ThrPt')
          label = 2
!
!    ...sphere with radius and central point
        case('Sphere')
          label = 3
!
!    ...cylinder with vector of direction ,  point of begining and
!       radius
        case('Cylinder')
          label = 4
!
!  .....ellipsoid with a given center and axes aligned with
!       the global system of coordinates
        case('Ellipsoid')
          label = 7
!
!    ...other cases
        case default
          write(*,*)'input_geometry.f: WRONG TYPE FOR THE SURFACE!!'
        end select
        SURFACES(isur)%Type = Type
        write(NIN,*) label, &
        '                         ...LABEL OF SURFACE', isur
!
!
        select case(label)
!
!  .......plane which is normal to a given vector and passing through
!         a point, orientation specified by the normal
          case(1)
            allocate( SURFACES(isur)%Rdata(6), STAT = ic)
            if (ic.ne.0) then
              write(*,*) 'Rdate IS NOT ALLOCATED', &
                         'SURFACE NUMBER = ',isur
              stop
            endif
            read(NWE,*)(SURFACES(isur)%Rdata(j),j=1,3)
            write(NIN,300)SURFACES(isur)%Rdata(1:3), &
                 '          ...COORDINATES OF THE POINTS'
            read(NWE,*)(SURFACES(isur)%Rdata(3+j),j=1,3)
            write(NIN,300)(SURFACES(isur)%Rdata(3+j),j=1,3), &
                 '          ...COORDINATES OF THE NORMAL VECTOR'
!
!  .......plane passing through three given points, A,B,C
!         orientation specified by AB x AC
          case(2)
             allocate( SURFACES(isur)%Rdata(9), STAT = ic)
             if(ic.ne.0)then
                write(*,*) 'Rdate IS  NOT ALLOCATED', &
                     'SURFACE NUMBER = ',isur
                stop
            endif
            read(NWE,*)( SURFACES(isur)%Rdata(j),j=1,3)
            write(NIN,300)(SURFACES(isur)%Rdata(j),j=1,3), &
                 '          ...COORDINATES OF THE FIRST POINT'
            read(NWE,*)(SURFACES(isur)%Rdata(3+j),j=1,3)
            write(NIN,300)(SURFACES(isur)%Rdata(3+j),j=1,3), &
                 '          ...COORDINATES OF THE SECOND POINT'
            read(NWE,*)(SURFACES(isur)%Rdata(6+j),j=1,3)
            write(NIN,300)(SURFACES(isur)%Rdata(6+j),j=1,3), &
                 '          ...COORDINATES OF THE THIRD POINT'
!
!  .......sphere
          case(3)
            allocate( SURFACES(isur)%Rdata(4), STAT = ic)
            if(ic.ne.0)then
              write(*,*) 'SURFACES COORDINATES ARE NOT ALLOCATED', &
                         'SURFACE NUMBER = ',isur
              stop
            endif
            read(NWE,*)(SURFACES(isur)%Rdata(j),j=1,3)
            write(NIN,300)(SURFACES(isur)%Rdata(j),j=1,3),  &
                 '          ...COORDINATES OF THE CENTER'
            read(NWE,*)SURFACES(isur)%Rdata(4)
            write(NIN,100)SURFACES(isur)%Rdata(4),  &
                 '                            ...RADIUS OF THE SPHERE'
!
!  .......cylinder
          case(4)
            allocate( SURFACES(isur)%Rdata(7), STAT = ic)
            if(ic.ne.0)then
              write(*,*) 'Rdate IS  NOT ALLOCATED', &
                         'SURFACE NUMBER = ',isur
              stop
            endif
            read(NWE,*)(SURFACES(isur)%Rdata(j),j=1,3)
            write(NIN,300)(SURFACES(isur)%Rdata(j),j=1,3)
            read(NWE,*)(SURFACES(isur)%Rdata(3+j),j=1,3)
            write(NIN,300)(SURFACES(isur)%Rdata(3+j),j=1,3)
            read(NWE,*) SURFACES(isur)%Rdata(7)
            write(NIN,*) SURFACES(isur)%Rdata(7)
!
!  .......ellipsoid with semi-axes a,b,c along x,y,z and given central
!         point
          case(7)
            allocate( SURFACES(isur)%Rdata(6), STAT = ic)
            if(ic.ne.0)then
              write(*,*) 'Rdata IS  NOT ALLOCATED', &
                         'SURFACE NUMBER = ',isur
              stop
            endif
!
!  .........read coordinates of the central point
            read(NWE,*)(SURFACES(isur)%Rdata(j),j=1,3)
!
!  .........read semiaxes lengths
            read(NWE,*)(SURFACES(isur)%Rdata(3+j),j=1,3)
!
!  .........write out the data to the new format input file
            write(NIN,300)(SURFACES(isur)%Rdata(j),j=1,3)
            write(NIN,300)(SURFACES(isur)%Rdata(3+j),j=1,3)
!
!
!  .......other cases
          case default
            write(*,*)'new_old: WRONG SURFACE TYPE'
            stop
        end select
!
        write(NIN,*)''
      enddo
!
!----------------------------------------------------------------------
!
!                            P O I N T S
!
!  ...read in number of points
      read(NWE,*) NRPOINT
      write(NIN,*)''
      write(NIN,*)NRPOINT, &
      '                         ...NUMBER OF POINTS'
      write(NIN,*)''
!
!  ...allocate the memory dynamically
      allocate( POINTS(NRPOINT), STAT=ic )
      if (ic.ne.0) then
        write(*,*) 'new_old: POINTS ARE NOT ALLOCATED'
        stop
      endif
!
      if (iprint.eq.1) then
        write(*,*) 'new_old: READING POINTS...'
      endif
!
!  ...read in and store the points data
      do ip=1,NRPOINT
        read(NWE,*)Type
        POINTS(ip)%Type  = Type
        POINTS(ip)%NrCurv = 0
!
        select case(Type)
!
!  .....read in and store the coordinates for regular points
        case('Regular')
          allocate( POINTS(ip)%Rdata(3), STAT = ic)
          if (ic.ne.0) then
            write(*,*) 'new_old: Rdata IS  NOT ALLOCATED', &
                       'POINT NUMBER IS', ip
            stop
          endif
          read(NWE,*) (POINTS(ip)%Rdata(k),k=1,3)
!
        case('CoorNrm')
          allocate( POINTS(ip)%Rdata(6), STAT = ic)
          if (ic.ne.0) then
            write(*,*) 'new_old: Rdata IS  NOT ALLOCATED', &
                       'POINT NUMBER IS', ip
            stop
          endif
          read(NWE,*) (POINTS(ip)%Rdata(k),k=1,6)
!
!  .....read in and store the intersecting surfaces data and
!       starting point for implicit  points
        case('Implicit')
          allocate (POINTS(ip)%Idata(3), STAT = ic)
          if (ic.ne.0) then
            write(*,*) 'new_old: Idata IS  NOT ALLOCATED', &
                       'POINT NUMBER IS', ip
            stop
          endif
          read(NWE,*) (POINTS(ip)%Idata(j),j=1,3)
          allocate( POINTS(ip)%Rdata(3), STAT = ic)
          if (ic.ne.0) then
            write(*,*) 'new_old: Rdata IS  NOT ALLOCATED', &
                       'POINT NUMBER IS', ip
            stop
          endif
          read(NWE,*) (POINTS(ip)%Rdata(j),j=1,3)
!
!
        case default
           write(*,*)'new_old: WRONG POINT TYPE'
        end select
!
      enddo
!
!----------------------------------------------------------------------
!
!                        C U R V E S
!
!
!  ...read in number of curves
      read(NWE,*) NRCURVE
!
      allocate( CURVES(NRCURVE), STAT=ic )
      if( ic.ne.0 )then
        print*,'CURVES ARE NOT ALLOCATED'
        stop
      endif
!
      if (iprint.eq.1) then
        write(*,*) 'new_old: READING CURVES...'
      endif
!
      do ic = 1, NRCURVE
        read(NWE,*) type
        CURVES(ic)%NrFig = 0
        CURVES(ic)%Type  = type
!
!  .....read in and store endpoints for the curve
        read(NWE,*) (CURVES(ic)%EndPoNo(j),j=1,2)
!
        select case (type)
        case('QuaCir','SegCir','QuaEl1','QuaEl2')
!
!  .......read in and store the coordinates of the center of a circle
          allocate(CURVES(ic)%Rdata(NDIM))
          read(NWE,*) (CURVES(ic)%Rdata(k),k=1,NDIM)
!
        case('QuaSEl')
!
!  .......read in and store the coordinates of the center of a
!         superellipse and powers
          allocate(CURVES(ic)%Rdata(NDIM+2))
          read(NWE,*) (CURVES(ic)%Rdata(k),k=1,NDIM+2)
!
        case('ImpCir')
!
!  .......read in and store the intersecting surfaces
          allocate( CURVES(ic)%Idata(4), STAT = icc)
          if (icc.ne.0) then
            write(*,*) 'Idata IS NOT ALLOCATED', &
                       'THE CURVE NUMBER IS ', icur
            stop
          endif
          read(NWE,*) (CURVES(ic)%Idata(k),k=1,4)
        end select
!
      enddo
!
!---------------------------------------------------------------------
!
!                       T R I A N G L E S
!
!  ...read in number of triangles
      read(NWE,*) NRTRIAN
      if (NRTRIAN.ne.0) then
        write(*,*) 'new_old: UNFINISHED '
        stop
      endif
!
!
!  ...allocate dynamically the memory
      allocate( TRIANGLES(NRTRIAN), STAT=ic )
      if( ic.ne.0 )then
        print*,'TRIANGLES ARE NOT ALLOCATED'
        stop
      endif
!
!---------------------------------------------------------------------
!
!                        R E C T A N G L E S
!
!  ...read in and store the rectangles data
      read(NWE,*) NRRECTA
!
!  ...allocate dynamically the memory
      allocate( RECTANGLES(NRRECTA), STAT=ic )
      if (ic.ne.0) then
        print*,'RECTANGLES ARE NOT ALLOCATED'
        stop
      endif
!
      if (iprint.eq.1) then
        write(*,*) 'new_old: READING RECTANGLES...'
      endif
!
      do irec=1,NRRECTA
        read(NWE,*) Type
        RECTANGLES(irec)%EdgeNo(1:4) = 0
        RECTANGLES(irec)%BlockNo(1:2) = 0
        RECTANGLES(irec)%Type  = Type
!
        read(NWE,*) RECTANGLES(irec)%VertNo(1:4)
!
!  .....read in and store surfaces that define the implicit rectangle
        if (type.eq.'ImpRec')then
          allocate(RECTANGLES(irec)%Idata(5), STAT = ic)
          if (ic.ne.0) then
            write(*,*) 'Idata IS NOT ALLOCATED', &
                       'RECTANGLE NUMBER IS', irec
            stop
          endif
          read(NWE,*) (RECTANGLES(irec)%Idata(j),j=1,5)
        endif
!
      enddo
!
      read(NWE,*) NRPRISM
!
!
!  ...read in number of hexas
      read(NWE,*) NRHEXAS
!
!  ...allocate dynamically the memory
      allocate(HEXAS(NRHEXAS), STAT=ic )
      if (ic.ne.0) then
        print*,'new_old: HEXAS ARE NOT ALLOCATED'
        stop
      endif
!
      if (iprint.eq.1) then
        write(*,*) 'new_old: READING HEXAS...'
      endif
!
      do ih=1,NRHEXAS
        read(NWE,*) HEXAS(ih)%Type
        read(NWE,*) HEXAS(ih)%VertNo(1:8)
!
!  ...end of loop through hexahedra
      enddo
!
!  ...determine missing connectivities and orientations
      call connect
!
!
!
!*************************************************************
!...write down the result into input file
!*************************************************************
!
      do i=1,NRPOINT
        select case(POINTS(i)%Type)
        case('Regular')
          write(NIN,*)1,  &
               '                         ...LABEL OF POINT',i
          write(NIN,*)POINTS(i)%NrCurv, &
               '                         ...NUMBER OF CURVES', &
               ' MEETING AT THE POINT '
          NrCurv = POINTS(i)%NrCurv
          select case(nrcurv)
            case(4)
               write(NIN,40)POINTS(i)%CurvNo(1:NrCurv), &
                    '             ...CURVES NUMBERS'
            case(5)
               write(NIN,50)POINTS(i)%CurvNo(1:NrCurv), &
                    '       ...CURVES NUMBERS'
            case(6)
               write(NIN,60)POINTS(i)%CurvNo(1:NrCurv), &
                    '       ...CURVES NUMBERS'
          end select

          write(NIN,300) (POINTS(i)%Rdata(k),k=1,3),'   ' &
               ,'       ...COORDINATES OF THE POINT'
          write(NIN,*)''
!
        case('CoorNrm')
          write(NOUT,*)30, &
               '                         ...LABEL OF POINT',i
          write(NIN,*)POINTS(i)%NrCurv, &
               '                   ...NUMBER OF CURVES', &
               ' MEETING AT THE POINT '
          NrCurv = POINTS(i)%NrCurv
          select case(nrcurv)
            case(4)
               write(NIN,40)POINTS(i)%CurvNo(1:NrCurv), &
                    '             ...CURVES NUMBERS'
            case(5)
               write(NIN,50)POINTS(i)%CurvNo(1:NrCurv), &
                    '       ...CURVES NUMBERS'
            case(6)
               write(NIN,60)POINTS(i)%CurvNo(1:NrCurv), &
                    '       ...CURVES NUMBERS'
          end select
          write(NIN,600) (POINTS(i)%Rdata(k),k=1,6),'   ' &
               ,'      ...COORDINATES AND NORMAL VECTORS'
          write(NIN,*)''
!
        case('Implicit')
          write(NIN,*) 50, &
               '                         ...LABEL OF POINT',i
          write(NIN,*) POINTS(i)%NrCurv, &
               '                         ...NUMBER OF CURVES', &
               ' MEETING AT THE POINT '
          NrCurv = POINTS(i)%NrCurv
          select case(nrcurv)
            case(4)
               write(NIN,40) POINTS(i)%CurvNo(1:NrCurv), &
                    '             ...CURVES NUMBERS'
            case(5)
               write(NIN,50) POINTS(i)%CurvNo(1:NrCurv), &
                    '       ...CURVES NUMBERS'
            case(6)
               write(NIN,60) POINTS(i)%CurvNo(1:NrCurv), &
                    '       ...CURVES NUMBERS'
          end select
            write(NIN,30)(POINTS(i)%Idata(j),j=1,3), &
                 '                   ...SURFACES THAT CONSTITUTE THE', &
                 ' POINT'
          write(NIN,300)(POINTS(i)%Rdata(k),k=1,3), &
                 '          ...INITIAL APPROXIMATION POINT'
!
        case default
          write(*,*) 'new_old: WRONG POINT LABEL', POINTS(i)%Type
        end select
      enddo
!
!
      write(NIN,*)''
      write(NIN,*)NRCURVE  ,'                         ...NUMBER', &
           ' OF CURVES'
      write(NIN,*)''
      do i=1,NRCURVE
        if(CURVES(i)%Type.eq.'NormCoord')then
            write(NIN,*)30,'                         ...TYPE OF CURVE',i
            write(NIN,*) (CURVES(i)%EndPoNo(j),j=1,2), &
                 '             ...END POINTS OF THE CURVE'
            NrFig = CURVES(i)%NrFig
            write(NIN,*)NrFig, &
                 '                         ...NUMBER OF FIGURES', &
                 ' MEETING AT THE CURVE'
            select case(nrfig)
            case(2)
               write(NIN,20)CURVES(i)%FigNo(1:nrfig), &
                    '                 ...FIGURES NUMBERS'
            case(3)
               write(NIN,30)CURVES(i)%FigNo(1:nrfig), &
                    '                 ...FIGURES NUMBERS'
            case(4)
               write(NIN,40)CURVES(i)%FigNo(1:nrfig), &
                    '               ...FIGURES NUMBERS'
            case(5)
               write(NIN,50)CURVES(i)%FigNo(1:nrfig), &
                    '               ...FIGURES NUMBERS'
            case(6)
               write(NIN,60)CURVES(i)%FigNo(1:nrfig), &
                    '              ...FIGURES NUMBERS'
            end select
            write(NIN,*)''
!
        elseif(CURVES(i)%Type.eq.'ImpCir')then
            write(NIN,*)50,'                         ...TYPE OF CURVE',i
            write(NIN,*) (CURVES(i)%EndPoNo(j),j=1,2), &
                 '             ...END POINTS OF THE CURVE'
            NrFig = CURVES(i)%NrFig
            write(NIN,*)NrFig, &
                 '                         ...NUMBER OF FIGURES', &
                 ' MEETING AT THE CURVES'
            select case(nrfig)
            case(2)
               write(NIN,20)CURVES(i)%FigNo(1:nrfig), &
                    '                   ...FIGURES NUMBERS'
            case(3)
               write(NIN,30)CURVES(i)%FigNo(1:nrfig), &
                    '                   ...FIGURES NUMBERS'
            case(4)
               write(NIN,40)CURVES(i)%FigNo(1:nrfig), &
                    '               ...FIGURES NUMBERS'
            case(5)
               write(NIN,50)CURVES(i)%FigNo(1:nrfig), &
                    '               ...FIGURES NUMBERS'
            case(6)
               write(NIN,60)CURVES(i)%FigNo(1:nrfig), &
                    '              ...FIGURES NUMBERS'
            end select
            write(NIN,40)(CURVES(i)%Idata(k),k=1,4), &
                 '             ...SURFACES CONSTITUTING THE CURVE'
            write(NIN,*)''
!
        elseif(CURVES(i)%Type.eq.'Seglin')then
            write(NIN,*)1,'                         ...TYPE OF CURVE',i
            write(NIN,*) (CURVES(i)%EndPoNo(j),j=1,2), &
                 '             ...END POINTS OF THE CURVE'
            NrFig = CURVES(i)%NrFig
            write(NIN,*)NrFig, &
                 '                         ...NUMBER OF FIGURES', &
                 ' MEETING AT THE CURVE'
            select case(nrfig)
            case(2)
               write(NIN,20)CURVES(i)%FigNo(1:nrfig), &
                    '                     ...FIGURES NUMBERS'
            case(3)
               write(NIN,30)CURVES(i)%FigNo(1:nrfig), &
                    '                   ...FIGURES NUMBERS'
            case(4)
               write(NIN,40)CURVES(i)%FigNo(1:nrfig), &
                    '             ...FIGURES NUMBERS'
            case(5)
               write(NIN,50)CURVES(i)%FigNo(1:nrfig), &
                    '                 ...FIGURES NUMBERS'
            case(6)
               write(NIN,60)CURVES(i)%FigNo(1:nrfig), &
                    '                ...FIGURES NUMBERS'
            end select
            write(NIN,*)''

        elseif(CURVES(i)%Type.eq.'QuaCir')then
            write(NIN,*)-1,'                         ...TYPE OF CURVE',i
            write(NIN,*) (CURVES(i)%EndPoNo(j),j=1,2), &
                 '             ...END POINTS OF THE CURVE'
            NrFig = CURVES(i)%NrFig
            write(NIN,*)NrFig, &
                 '                         ...NUMBER OF FIGURES', &
                 ' MEETING AT THE CURVE'
            select case(nrfig)
            case(2)
               write(NIN,20)CURVES(i)%FigNo(1:nrfig), &
                    '                     ...FIGURES NUMBERS'
            case(3)
               write(NIN,30)CURVES(i)%FigNo(1:nrfig), &
                    '                   ...FIGURES NUMBERS'
            case(4)
               write(NIN,40)CURVES(i)%FigNo(1:nrfig), &
                    '                 ...FIGURES NUMBERS'
            case(5)
               write(NIN,50)CURVES(i)%FigNo(1:nrfig), &
                    '                 ...FIGURES NUMBERS'
            case(6)
               write(NIN,60)CURVES(i)%FigNo(1:nrfig), &
                    '                ...FIGURES NUMBERS'
            end select
            write(NIN,*)''

         elseif(CURVES(i)%Type.eq.'SegCir')then
            write(NIN,*)-2,'                         ...TYPE OF CURVE',i
            write(NIN,*) (CURVES(i)%EndPoNo(j),j=1,2), &
                 '             ...END POINTS OF THE CURVE'
            NrFig = CURVES(i)%NrFig
            write(NIN,*)NrFig, &
                 '                         ...NUMBER OF FIGURES', &
                 ' MEETING AT THE CURVE'
            select case(nrfig)
            case(2)
               write(NIN,20)CURVES(i)%FigNo(1:nrfig), &
                    '                     ...FIGURES NUMBERS'
            case(3)
               write(NIN,30)CURVES(i)%FigNo(1:nrfig), &
                    '                   ...FIGURES NUMBERS'
            case(4)
               write(NIN,40)CURVES(i)%FigNo(1:nrfig), &
                    '                 ...FIGURES NUMBERS'
            case(5)
               write(NIN,50)CURVES(i)%FigNo(1:nrfig), &
                    '                 ...FIGURES NUMBERS'
            case(6)
               write(NIN,60)CURVES(i)%FigNo(1:nrfig), &
                    '                ...FIGURES NUMBERS'
            end select
            write(NIN,*)''

        elseif(CURVES(i)%Type.eq.'QuaEl1')then
            write(NIN,*)-3,'                         ...TYPE OF CURVE',i
            write(NIN,*) (CURVES(i)%EndPoNo(j),j=1,2), &
                 '             ...END POINTS OF THE CURVE'
            NrFig = CURVES(i)%NrFig
            write(NIN,*)NrFig, &
                 '                         ...NUMBER OF FIGURES', &
                 ' MEETING AT THE CURVE'
            select case(nrfig)
            case(2)
               write(NIN,20)CURVES(i)%FigNo(1:nrfig), &
                    '                     ...FIGURES NUMBERS'
            case(3)
               write(NIN,30)CURVES(i)%FigNo(1:nrfig), &
                    '                   ...FIGURES NUMBERS'
            case(4)
               write(NIN,40)CURVES(i)%FigNo(1:nrfig), &
                    '                 ...FIGURES NUMBERS'
            case(5)
               write(NIN,50)CURVES(i)%FigNo(1:nrfig), &
                    '                 ...FIGURES NUMBERS'
            case(6)
               write(NIN,60)CURVES(i)%FigNo(1:nrfig), &
                    '                ...FIGURES NUMBERS'
            end select
            write(NIN,*)''

        elseif(CURVES(i)%Type.eq.'QuaEl2')then
            write(NIN,*)-4,'                         ...TYPE OF CURVE',i
            write(NIN,*) (CURVES(i)%EndPoNo(j),j=1,2), &
                 '             ...END POINTS OF THE CURVE'
            NrFig = CURVES(i)%NrFig
            write(NIN,*)NrFig, &
                 '                         ...NUMBER OF FIGURES', &
                 ' MEETING AT THE CURVE'
            select case(nrfig)
            case(2)
               write(NIN,20)CURVES(i)%FigNo(1:nrfig), &
                    '                     ...FIGURES NUMBERS'
            case(3)
               write(NIN,30)CURVES(i)%FigNo(1:nrfig), &
                    '                   ...FIGURES NUMBERS'
            case(4)
               write(NIN,40)CURVES(i)%FigNo(1:nrfig), &
                    '                 ...FIGURES NUMBERS'
            case(5)
               write(NIN,50)CURVES(i)%FigNo(1:nrfig), &
                    '                 ...FIGURES NUMBERS'
            case(6)
               write(NIN,60)CURVES(i)%FigNo(1:nrfig), &
                    '                ...FIGURES NUMBERS'
            end select
            write(NIN,*)''

        elseif(CURVES(i)%Type.eq.'QuaSEl')then
            write(NIN,*)-5,'                         ...TYPE OF CURVE',i
            write(NIN,*) (CURVES(i)%EndPoNo(j),j=1,2), &
                 '             ...END POINTS OF THE CURVE'
            NrFig = CURVES(i)%NrFig
            write(NIN,*)NrFig, &
                 '                         ...NUMBER OF FIGURES', &
                 ' MEETING AT THE CURVE'
            select case(nrfig)
            case(2)
               write(NIN,20)CURVES(i)%FigNo(1:nrfig), &
                    '                     ...FIGURES NUMBERS'
            case(3)
               write(NIN,30)CURVES(i)%FigNo(1:nrfig), &
                    '                   ...FIGURES NUMBERS'
            case(4)
               write(NIN,40)CURVES(i)%FigNo(1:nrfig), &
                    '                 ...FIGURES NUMBERS'
            case(5)
               write(NIN,50)CURVES(i)%FigNo(1:nrfig), &
                    '                 ...FIGURES NUMBERS'
            case(6)
               write(NIN,60)CURVES(i)%FigNo(1:nrfig), &
                    '                ...FIGURES NUMBERS'
            end select
            write(NIN,*)''
        else
            write(*,*)'INPUT_HEAD: WRONG CURVE LABEL', CURVES(i)%Type
        endif
      enddo
!
!
      write(NIN,*)''
      write(NIN,*) NRTRIAN,  &
           '                         ...NUMBER OF TRIANGLES'
!
!
      write(NIN,*)''
      write(NIN,*)''
      write(NIN,*)NRRECTA ,  &
           '                         ...NUMBER OF RECTANGLES'
      do i = 1,NRRECTA
         if(RECTANGLES(i)%Type.eq.'BilQua')then
            write(NIN,*)''
            write(NIN,*)1, &
                 '                         ...TYPE OF RECTANGLE',i
            write(NIN,40)(RECTANGLES(i)%EdgeNo(j),j=1,4), &
                 '             ...FOUR CURVES NUMBERS'
            write(NIN,20)RECTANGLES(i)%BlockNo(1:2), &
                 '                         ...ADJACENT', &
                 'HEXAHEDRONS NUMBERS'
!
         else if(RECTANGLES(i)%Type.eq.'TraQua')then
            write(NIN,*)''
            write(NIN,*)3, &
                 '                         ...TYPE OF RECTANGLE',i
            write(NIN,40)(RECTANGLES(i)%EdgeNo(j),j=1,4), &
                 '             ...FOUR CURVES NUMBERS'
            write(NIN,20)RECTANGLES(i)%BlockNo(1:2), &
            '                         ...ADJACENT HEXAHEDRONS NUMBERS'
!
          else if(RECTANGLES(i)%Type.eq.'RatioRec')then
            write(NIN,*)''
            write(NIN,*)30, &
                 '                         ...TYPE OF RECTANGLE',i
            write(NIN,40)(RECTANGLES(i)%EdgeNo(j),j=1,4), &
                 '             ...FOUR CURVES NUMBERS'
            write(NIN,20)RECTANGLES(i)%BlockNo(1:2), &
            '                         ...ADJACENT HEXAHEDRON NUMBERS'
!
         else if(RECTANGLES(i)%Type.eq.'ImpRec')then
            write(NIN,*)''
            write(NIN,*)50, &
                 '                         ...TYPE OF RECTANGLE',i
            write(NIN,40)(RECTANGLES(i)%EdgeNo(j),j=1,4), &
                 '             ...FOUR CURVES NUMBERS'
            write(NIN,20)RECTANGLES(i)%BlockNo(1:2), &
            '                         ...ADJACENT HEXAHEDRONS NUMBERS'
            write(NIN,50)RECTANGLES(i)%Idata(1:5), &
                 '       ...SURFACES THAT CONSISTUTE THE RECTANGLE'
!
         elseif(RECTANGLES(i)%Type.eq.'CylRec')then
            write(NIN,*)''
            write(NIN,*)-1, &
                 '                        ...TYPE OF RECTANGLE',i
            write(NIN,40)(RECTANGLES(i)%EdgeNo(j),j=1,4), &
                 '             ...FOUR CURVES NUMBERS'
            write(NIN,20)RECTANGLES(i)%BlockNo(1:2), &
            '                         ...ADJACENT HEXAHEDRONS NUMBERS'
         else
            write(*,*)'INPUT_HEAD: WRONG RECTANGLE LABEL',  &
                 RECTANGLES(i)%Type
         endif
      enddo
!
!
      write(NIN,*)''
      write(NIN,*) NRPRISM  , &
           '                         ...NUMBER OF PRISMS'
!
!
      write(NIN,*)''
      write(NIN,*)''
      write(NIN,*) NRHEXAS  ,  &
           '                         ...NUMBER OF HEXAHEDRONS'
      do iii = 1,  NRHEXAS
         write(NIN,*)''
         write(NIN,*)1, &
              '                         ...TYPE OF HEXAHEDRON',III
         write(NIN,60) (HEXAS(iii)%FigNo(j), j = 1,6), &
              ' ...RECTANGULAR FACES NUMBERS'
       enddo
!
 100   format(1x,1f9.3,5a)
 300   format(1x,3f9.3,5a)
 600   format(1x,6f9.3,5a)
!
 20    format(1x,2i6,9a)
 30    format(1x,3i6,9a)
 40    format(1x,4i6,9a)
 50    format(1x,5i6,9a)
 60    format(1x,6i6,9a)
!
!
     end subroutine new_old
