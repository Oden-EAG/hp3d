c----------------------------------------------------------------------
c
c   routine name       - export_gmps2vtk
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 11
c
c   purpose            - export GMP outer surface to vtk legacy format
c
c   Remarks: routine is written for tets only!
c
c----------------------------------------------------------------------
c
      subroutine export_gmps2vtk(Fname)
c
      use GMP
      use error
#include "syscom.blk"
c
c  ...element nodes and their orientations
      dimension nvert(3)
      character(len=64) :: Fname
c
      nfile = 31
c
c  ...open file for writing
      open(unit=nfile, file=Fname, form='formatted')
c
c  ...deallocate Idata for all points
      do i=1,NRPOINT
        if(associated(POINTS(i)%Idata)) then
          write(*,*)'Idata associated for point = ',i
          call pause
          deallocate(POINTS(nv)%Idata, stat=istat)
          if (istat.ne.SUCCESS) then
            call logic_error(ERR_ALLOC_FAILURE,
     .        __FILE__,__LINE__)
          endif
        endif
      enddo
c
c  ...Step 1: loop over tets, count and rename surface triangles
      ntria=0
      do i=1,NRTETRA
        do j=1,4
          call decode(TETRAS(i)%FigNo(j), nt,lab)
          if (TRIANGLES(nt)%BlockNo(2).eq.0) then
            ntria=ntria+1
            TRIANGLES(nt)%Type = 'PTITri'
          endif
        enddo
      enddo
c
c  ...Step 2: loop over surface triangles and renumber vertices
      n=0
      do i=1,NRTRIAN
        if (TRIANGLES(i)%Type.ne.'PTITri') cycle
        do j=1,3
          nv = TRIANGLES(i)%VertNo(j)
c  .......if 1st visit to point
          if (.not.associated(POINTS(nv)%Idata)) then
            n = n + 1
            allocate(POINTS(nv)%Idata(1), stat=istat)
            if (istat.ne.SUCCESS) then
              call logic_error(ERR_ALLOC_FAILURE,
     .          __FILE__,__LINE__)
            endif
            POINTS(nv)%Idata(1) = n
          endif
        enddo
      enddo
c
c  ...Step 3: loop over surface triangles and write point's coordinates
      write(nfile,1000)
 1000 format('# vtk DataFile Version 2.0')
      write(nfile,1001)
 1001 format('hp3d export mesh')
      write(nfile,1002)
 1002 format('ASCII')
      write(nfile,1003)
 1003 format('DATASET UNSTRUCTURED_GRID')
      write(nfile,1004)n
 1004 format('POINTS ',i10,' double')
      do i=1,NRTRIAN
        if (TRIANGLES(i)%Type.ne.'PTITri') cycle
        do j=1,3
          nv = TRIANGLES(i)%VertNo(j)
c  .......if 1st visit to point
          if (POINTS(nv)%Idata(1).gt.0) then
c  .........convert coordinates from [m] to [cm]
            POINTS(nv)%Rdata(1:3) = POINTS(nv)%Rdata(1:3)*100.d0
            write(nfile,7001) POINTS(nv)%Rdata(1:3)
 7001       format(3(f20.12,2x))
c  .........flip sign
            POINTS(nv)%Idata(1) = -POINTS(nv)%Idata(1)
          endif
        enddo
      enddo
c
c  ...Step 4: write elements
      write(nfile,1005)ntria,ntria*4
 1005 format('CELLS ',2(i10,2x))
      do i=1,NRTRIAN
        if (TRIANGLES(i)%Type.ne.'PTITri') cycle
c  .....collect vertices
        do j=1,3
          nv = TRIANGLES(i)%VertNo(j)
          nvert(j) = -POINTS(nv)%Idata(1) - 1
        enddo
        write(nfile,7003) nvert(1:3)
 7003   format('3 ',3(i10,2x))
      enddo
      write(nfile,1006)ntria
 1006 format('CELL_TYPES ',i10)
      do i=1,ntria
        write(nfile,1007)
 1007  format('5')
      enddo
c
c  ...Step 5: write surface number
      write(nfile,1008)ntria
 1008 format('CELL_DATA ',i10)
      write(nfile,1009)
 1009 format('SCALARS  int  1')
      write(nfile,1010)
 1010 format('LOOKUP_TABLE default')
 1011 format('14')
      do i=1,ntria
        write(nfile,1011)
      enddo
c
c  ...Step 5: deallocate Idata
      do i=1,NRPOINT
        if (.not.associated(POINTS(i)%Idata)) cycle
        deallocate(POINTS(i)%Idata, stat=istat)
        if (istat.ne.SUCCESS) then
          call logic_error(ERR_ALLOC_FAILURE,
     .      __FILE__,__LINE__)
        endif
      enddo
c
      close(nfile)
c
      end subroutine
