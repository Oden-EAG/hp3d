c---------------------------------------------------------------------
c
c   routine name       - export_gmp2vtk
c
c---------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - export phantom to vtk legacy format
c
c-----------------------------------------------------------------------
c
      subroutine export_gmp2vtk(Fname)
c
      use GMP
      use error
#include "syscom.blk"
c
c  ...element nodes and their orientations
      dimension nvert(4)
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
c  ...Step 1: loop over tetra and renumber vertices
      nver = 0
      do i=1,NRTETRA
        do j=1,4
          nv = TETRAS(i)%VertNo(j)
c  .......if 1st visit to point
          if (.not.associated(POINTS(nv)%Idata)) then
            nver = nver + 1
            allocate(POINTS(nv)%Idata(1), stat=istat)
            if (istat.ne.SUCCESS) then
              call logic_error(ERR_ALLOC_FAILURE,
     .          __FILE__,__LINE__)
            endif
            POINTS(nv)%Idata(1) = nver
          endif
        enddo
      enddo
c
c  ...Step 2: loop over tetras and write point's coordinates
      write(nfile,1000)
 1000 format('# vtk DataFile Version 2.0')
      write(nfile,1001)
 1001 format('hp3d export mesh')
      write(nfile,1002)
 1002 format('ASCII')
      write(nfile,1003)
 1003 format('DATASET UNSTRUCTURED_GRID')
      write(nfile,1004)nver
 1004 format('POINTS ',i10,' double')
      do i=1,NRTETRA
        do j=1,4
          nv = TETRAS(i)%VertNo(j)
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
c  ...Step 3: write elements
      write(nfile,1005)NRTETRA,NRTETRA*5
 1005 format('CELLS ',2(i10,2x))
      do i=1,NRTETRA
c  .....collect vertices
        do j=1,4
          nv = TETRAS(i)%VertNo(j)
          nvert(j) = -POINTS(nv)%Idata(1) - 1
        enddo
        write(nfile,7003) nvert(1:4)
 7003   format('4 ',4(i10,2x))
      enddo
      write(nfile,1006)NRTETRA
 1006 format('CELL_TYPES ',i10)
      do i=1,NRTETRA
        write(nfile,1007)
 1007  format('10')
      enddo
c
c  ...Step 4: deallocate Idata
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
