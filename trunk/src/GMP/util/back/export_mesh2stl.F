c---------------------------------------------------------------------
c
c   routine name       - export mesh
c
c---------------------------------------------------------------------
c
c   latest revision    - Feb 11
c
c   purpose            - given a mesh in RECONSTRUCT format, routine
c                        exports 'G1RecTri' triangles is STL format
c
c   arguments :
c                 in:  - Fp: file name
c
c-----------------------------------------------------------------------
c
      subroutine export_mesh2stl(Fname)
c
      use GMP
#include "syscom.blk"
c
c  ...file name
      character(len=64) :: Fname
c  ...vertices of triangle ordered in counter clockwise fashion
      dimension :: nvert(3),rnorm(3),aux1(3),aux2(3)
c
      nfile=33
c
c  ...open file for writing
      open(unit=nfile, file=Fname, form='formatted')
c
      write(*,*)'exporting mesh to stl format...'
c
      write(nfile,*)'solid '
      ntri = 0
c  ...loop over triangles
      do i = 1,NRTRIAN
        if (TRIANGLES(i)%Type.ne.'G1RecTri') cycle
        ntri = ntri + 1
c  .....get adjacent tet
        nick = TRIANGLES(i)%BlockNo(1)
        call decode(nick, ntet,lab)
c  .....loop over tet faces
        do ifig = 1,4
          call decode(TETRAS(ntet)%FigNo(ifig), nt,no)
c  .......if face is found
          if (nt.eq.i) then
c  .........select local face number
            select case(ifig)
c  .........clockwise faces (from element_data module!!!)
            case(1,4)
              if (no.le.2) then
                norient = 1
              elseif (no.le.5) then
                norient = 0
              else
                write(*,*)'export_mesh2STL: invalid no = ',no
                stop
              endif
c  .........counterclockwise faces (from element_data module!!!)
            case(2,3)
              if (no.le.2) then
                norient = 0
              elseif (no.le.5) then
                norient = 1
              else
                write(*,*)'export_mesh2STL: invalid no = ',no
                stop
              endif
            end select
            exit
          endif
c  .....end of loop over faces
        enddo
c  .....save vertices in counter clockwise fashion
        nvert(1:3) = TRIANGLES(i)%VertNo(1:3)
        if (norient.eq.1) then
          ntemp = nvert(2)
          nvert(2) = nvert(3)
          nvert(3) = ntemp
        endif
c  .....compute normal
        aux1(1:3) = POINTS(nvert(2))%Rdata(1:3) -
     .              POINTS(nvert(1))%Rdata(1:3)
        aux2(1:3) = POINTS(nvert(3))%Rdata(1:3) -
     .              POINTS(nvert(1))%Rdata(1:3)
        call cross_product(aux1(1:3),aux2(1:3), rnorm(1:3))
        call normalize(rnorm(1:3))
c
c  .....P  R  I  N  T  I  N  G
        write(nfile,1001)rnorm(1:3)
 1001   format(' facet normal ',3(e20.13,3x))
        write(nfile,*)'  outer loop'
c  .....loop over vertices
        do j = 1,3
          np = nvert(j)
          write(nfile,1000)POINTS(np)%Rdata(1:3)
 1000     format('     vertex  ',3(e20.13,3x))
c  .....end of loop over vertices
        enddo
        write(nfile,*)'  endloop'
        write(nfile,*)'endfacet'
c  ...end of loop over triangles
      enddo
c
      write(nfile,*)'endsolid '
c
c  ...close file
      close(nfile)
c
      write(*,*)'mesh exported to stl format!'
c
c
      end
