!-----------------------------------------------------------------------
!> @brief interactive routine to print GMP data structure
!!
!!
!> @date Jul 2023
!-----------------------------------------------------------------------
!
      subroutine GMP_result
!
      use GMP
      implicit none
!
!  ...misc
      integer :: idec,i
!
!-----------------------------------------------------------------------
!
   10 write(*,*)
      write(*,*) 'GMP_result: SELECT OBJECT AND NUMBER'
      write(*,*) '            EXIT...................0'
      write(*,*) '            SURFACE  ..............1'
      write(*,*) '            POINT  ................2'
      write(*,*) '            CURVE  ................3'
      write(*,*) '            TRIANGLE  .............4'
      write(*,*) '            RECTANGLE  ............5'
      write(*,*) '            PRISM  ................6'
      write(*,*) '            HEXAHEDRON ............7'
      write(*,*) '            TETRAHEDRON ...........8'
      write(*,*) '            PYRAMID ...............9'
!
      write(*,7005) NRSURFS,NRPOINT,NRCURVE,NRTRIAN,NRRECTA,NRPRISM,NRHEXAS,NRTETRA,NRPYRAM

 7005 format(' NRSURFS,NRPOINT,NRCURVE,NRTRIAN,NRRECTA,NRPRISM,NRHEXAS,NRTETRA,NRPYRAM = ',&
             /,9i8)
!
      read(*,*) idec,i
      if (idec.eq.0) return
      select case(idec)
!
!  ...surface data
      case(1)
        if ((i.lt.1).or.(i.gt.NRSURFS)) go to 10
        write(*,7010) SURFACES(i)%Type,SURFACES(i)%Idata
 7010   format('SURFACES(i)%Type  = ',a10,/,&
               'SURFACES(i)%Idata = ',10i6)
        write(*,7011) SURFACES(i)%Rdata
 7011   format('SURFACES(i)%Rdata = ',10f8.3)
!
!  ...point data
      case(2)
        if ((i.lt.1).or.(i.gt.NRPOINT)) go to 10
        write(*,7020) POINTS(i)%Type,POINTS(i)%NrCurv,POINTS(i)%CurvNo
 7020   format('POINTS(i)%Type   = ',a10,/,&
               'POINTS(i)%NrCurv = ',i3,/,&
               'POINTS(i)%CurvNo = ',30i6)
        if (associated(POINTS(i)%Idata)) write(*,7021) POINTS(i)%Idata
 7021   format('POINTS(i)%Idata  = ',10i6)
        if (associated(POINTS(i)%Rdata)) write(*,7022) POINTS(i)%Rdata
 7022   format('POINTS(i)%Rdata  = ',10f8.3)
!
!  ...curve data
      case(3)
        if ((i.lt.1).or.(i.gt.NRCURVE)) go to 10
        write(*,7030) CURVES(i)%Type,CURVES(i)%EndPoNo,&
                      CURVES(i)%NrFig,CURVES(i)%FigNo
 7030   format('CURVES(i)%Type    = ',a10,/,&
               'CURVES(i)%EndPoNo = ',2i6,/,&
               'CURVES(i)%NrFig   = ',i4,/,&
               'CURVES(i)%FigNo   = '10i6)
        if (associated(CURVES(i)%Idata)) write(*,7031) CURVES(i)%Idata
 7031   format('CURVES(i)%Idata = ',10i6)
        if (associated(CURVES(i)%Rdata)) write(*,7032) CURVES(i)%Rdata
 7032   format('CURVES(i)%Rdata = ',10f8.3)
!
!  ...triangle data
      case(4)
        if ((i.lt.1).or.(i.gt.NRTRIAN)) go to 10
        write(*,7040) TRIANGLES(i)%Type,TRIANGLES(i)%VertNo,&
                      TRIANGLES(i)%EdgeNo,TRIANGLES(i)%BlockNo
 7040   format('TRIANGLES(i)%Type   = ',a10,/,&
               'TRIANGLES(i)%VertNo = ',3i6,/,&
               'TRIANGLES(i)%EdgeNo = ',3i6,/,&
               'TRIANGLES(i)%BlockNo = ',2i6)
        if (associated(TRIANGLES(i)%Idata)) write(*,7041) TRIANGLES(i)%Idata
 7041   format('TRIANGLES(i)%Idata = ',30i10)
!
!  ...rectangle data
      case(5)
        if ((i.lt.1).or.(i.gt.NRRECTA)) go to 10
        write(*,7050) RECTANGLES(i)%Type,RECTANGLES(i)%VertNo, &
                      RECTANGLES(i)%EdgeNo,RECTANGLES(i)%BlockNo 
 7050   format('RECTANGLES(i)%Type    = ',a10,/,&
               'RECTANGLES(i)%VertNo  = ',4i6,/,&
               'RECTANGLES(i)%EdgeNo  = ',4i6,/,&
               'RECTANGLES(i)%BlockNo = ',2i6)
        if (associated(RECTANGLES(i)%Idata)) write(*,7051) RECTANGLES(i)%Idata
 7051   format('RECTANGLES(i)%Idata = ',30i10)
        if (associated(RECTANGLES(i)%Rdata)) write(*,7052) RECTANGLES(i)%Rdata
 7052   format('RECTANGLES(i)%Rdata = ',30f8.3)
!
!  ...prism data
      case(6)
        if ((i.lt.1).or.(i.gt.NRPRISM)) go to 10
        write(*,7060) PRISMS(i)%Type, PRISMS(i)%VertNo,PRISMS(i)%EdgeNo, &
                      PRISMS(i)%FigNo,PRISMS(i)%Domain
 7060   format('PRISMS(i)%Type   = ',a10,/,&
               'PRISMS(i)%VertNo = ',6i6,/,&
               'PRISMS(i)%EdgeNo = ',9i6,/,&
               'PRISMS(i)%FigNo  = ',5i6,/,&
               'PRISMS(i)%Domain = ',i3)
!
!  ...hexa data
      case(7)
        if ((i.lt.1).or.(i.gt.NRHEXAS)) go to 10
        write(*,7070) HEXAS(i)%Type,HEXAS(i)%VertNo,HEXAS(i)%EdgeNo,&
                      HEXAS(i)%FigNo,HEXAS(i)%Domain
 7070   format('HEXAS(i)%Type   = ',a10,/,&
               'HEXAS(i)%VertNo = ',8i6,/,&
               'HEXAS(i)%EdgeNo = ',12i6,/,&
               'HEXAS(i)%FigNo  = ',6i6,/,&
               'HEXAS(i)%Domain = ',i3)
        if (associated(HEXAS(i)%Idata)) write(*,7071) HEXAS(i)%Idata
 7071   format('HEXAS(i)%Idata = ',10i10)
!
!  ...tet data
      case(8)
        if ((i.lt.1).or.(i.gt.NRTETRA)) go to 10
        write(*,7080) TETRAS(i)%Type,TETRAS(i)%VertNo,&
                      TETRAS(i)%EdgeNo,TETRAS(i)%FigNo,TETRAS(i)%Domain
 7080   format('TETRAS(i)%Type   = ',a10,/,&
               'TETRAS(i)%VertNo = ',4i6,/,&
               'TETRAS(i)%EdgeNo = ',6i6,/,&
               'TETRAS(i)%FigNo  = ',4i6,/,&
               'TETRAS(i)%Domain = ',i4)
        if (associated(TETRAS(i)%Idata)) write(*,7081) TETRAS(i)%Idata
 7081   format('TETRAS(i)%Idata = ',30i10)
!
!  ...pyramid data
      case(9)
        if ((i.lt.1).or.(i.gt.NRPYRAM)) go to 10
        write(*,7090) PYRAMIDS(i)%Type,PYRAMIDS(i)%VertNo,PYRAMIDS(i)%EdgeNo, &
                      PYRAMIDS(i)%FigNo,PYRAMIDS(i)%Domain 
 7090   format('PYRAMIDS(i)%Type   = ',a10,/, &
               'PYRAMIDS(i)%VertNo = ',5i6,/,&
               'PYRAMIDS(i)%EdgeNo = ',8i6,/,&
               'PYRAMIDS(i)%FigNo  = ',5i6,/,&
               'PYRAMIDS(i)%Domain = ',i3)
        if (associated(PYRAMIDS(i)%Idata)) write(*,7091) PYRAMIDS(i)%Idata
 7091   format('PYRAMIDS(i)%Idata = ',30i10)
      case default
        write(*,*) 'GMP_result: UNDEFINED OPTION'
        go to 10
      end select
      go to 10
!
!    
      end subroutine GMP_result

