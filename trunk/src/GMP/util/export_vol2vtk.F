c----------------------------------------------------------------------
c
c   routine name       - export_vol2vtk
c
c----------------------------------------------------------------------
c
c   latest revision    - Jan 11
c
c   purpose            - export domains to vtk format
c
c   arguments:
c                in:   - Fname: file name      
c                      - Nrd  : number of surfaces to be exported
c                      - Nlist: list of surfaces to be exported
c
c-----------------------------------------------------------------------
c
      subroutine export_vol2vtk(Fname,Nrd,Nlist)
c
      use error
      use GMP
#include "syscom.blk"
c
c  ...vertices of surface triangle in the new enumeration
      dimension nvert(8)
      character(len=4),dimension(100000) :: ctype
      character(len=64) :: Fname
      dimension Nlist(Nrd)
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
c  ...Step 1: loop over elements and renumber vertices
c      
      ncell=0 ; npoin=0 ; ncoun=0
      do n=1,NRPRISM
        ndom = PRISMS(n)%Domain
        call locate(ndom,Nlist,Nrd, nloc)
        if (nloc.eq.0) cycle
        ncell=ncell+1
        ncoun=ncoun+7
        ctype(ncell)='pris'
        do i=1,6
          nv=PRISMS(n)%VertNo(i)
c  .......if 1st visit to point          
          if (.not.associated(POINTS(nv)%Idata)) then
            npoin=npoin+1
            allocate(POINTS(nv)%Idata(1), stat=istat)
            if (istat.ne.SUCCESS) then
              call logic_error(ERR_ALLOC_FAILURE,
     .          __FILE__,__LINE__)
            endif
            POINTS(nv)%Idata(1) = npoin
          endif
        enddo
      enddo
c
      do n=1,NRHEXAS
        ndom = HEXAS(n)%Domain
        call locate(ndom,Nlist,Nrd, nloc)
        if (nloc.eq.0) cycle
        ncell=ncell+1
        ncoun=ncoun+9
        ctype(ncell)='hexa'
        do i=1,8
          nv=HEXAS(n)%VertNo(i)
c  .......if 1st visit to point          
          if (.not.associated(POINTS(nv)%Idata)) then
            npoin=npoin+1
            allocate(POINTS(nv)%Idata(1), stat=istat)
            if (istat.ne.SUCCESS) then
              call logic_error(ERR_ALLOC_FAILURE,
     .          __FILE__,__LINE__)
            endif
            POINTS(nv)%Idata(1) = npoin
          endif
        enddo
      enddo
c
      do n=1,NRTETRA
        ndom = TETRAS(n)%Domain
        call locate(ndom,Nlist,Nrd, nloc)
        if (nloc.eq.0) cycle

        call locate_tetra(n, iflag)
        if (iflag.eq.0) cycle

        ncell=ncell+1
        ncoun=ncoun+5
        ctype(ncell)='tetr'
        do i=1,4
          nv=TETRAS(n)%VertNo(i)
c  .......if 1st visit to point          
          if (.not.associated(POINTS(nv)%Idata)) then
            npoin=npoin+1
            allocate(POINTS(nv)%Idata(1), stat=istat)
            if (istat.ne.SUCCESS) then
              call logic_error(ERR_ALLOC_FAILURE,
     .          __FILE__,__LINE__)
            endif
            POINTS(nv)%Idata(1) = npoin
          endif
        enddo
      enddo
c
      do n=1,NRPYRAM
        ndom = PYRAMIDS(n)%Domain
        call locate(ndom,Nlist,Nrd, nloc)
        if (nloc.eq.0) cycle
        ncell=ncell+1
        ncoun=ncoun+6
        ctype(ncell)='pyra'
        do i=1,5
          nv=PYRAMIDS(n)%VertNo(i)
c  .......if 1st visit to point          
          if (.not.associated(POINTS(nv)%Idata)) then
            npoin=npoin+1
            allocate(POINTS(nv)%Idata(1), stat=istat)
            if (istat.ne.SUCCESS) then
              call logic_error(ERR_ALLOC_FAILURE,
     .          __FILE__,__LINE__)
            endif
            POINTS(nv)%Idata(1) = npoin
          endif
        enddo
      enddo
c
c  ...Step 2
c
      write(nfile,1000) 
 1000 format('# vtk DataFile Version 2.0')
      write(nfile,1001)
 1001 format('hp3d export mesh')
      write(nfile,1002)
 1002 format('ASCII')
      write(nfile,1003)
 1003 format('DATASET UNSTRUCTURED_GRID')
      write(nfile,1004)npoin
 1004 format('POINTS ',i10,' double')
c
      do n=1,NRPRISM
        ndom = PRISMS(n)%Domain
        call locate(ndom,Nlist,Nrd, nloc)
        if (nloc.eq.0) cycle
        do i=1,6
          nv=PRISMS(n)%VertNo(i)
c  .......if 1st visit to point          
          if (POINTS(nv)%Idata(1).gt.0) then
            write(nfile,7001) POINTS(nv)%Rdata(1:3)
 7001       format(3(f20.12,2x))     
c  .........flip sign            
            POINTS(nv)%Idata(1) = -POINTS(nv)%Idata(1)
          endif
        enddo
      enddo
c
      do n=1,NRHEXAS
        ndom = HEXAS(n)%Domain
        call locate(ndom,Nlist,Nrd, nloc)
        if (nloc.eq.0) cycle
        do i=1,8
          nv=HEXAS(n)%VertNo(i)
c  .......if 1st visit to point          
          if (POINTS(nv)%Idata(1).gt.0) then
            write(nfile,7001) POINTS(nv)%Rdata(1:3)
c  .........flip sign            
            POINTS(nv)%Idata(1) = -POINTS(nv)%Idata(1)
          endif
        enddo
      enddo
c
      do n=1,NRTETRA
        ndom = TETRAS(n)%Domain
        call locate(ndom,Nlist,Nrd, nloc)
        if (nloc.eq.0) cycle

        call locate_tetra(n, iflag)
        if (iflag.eq.0) cycle

        do i=1,4
          nv=TETRAS(n)%VertNo(i)
c  .......if 1st visit to point          
          if (POINTS(nv)%Idata(1).gt.0) then
            write(nfile,7001) POINTS(nv)%Rdata(1:3)
c  .........flip sign            
            POINTS(nv)%Idata(1) = -POINTS(nv)%Idata(1)
          endif
        enddo
      enddo
c
      do n=1,NRPYRAM
        ndom = PYRAMIDS(n)%Domain
        call locate(ndom,Nlist,Nrd, nloc)
        if (nloc.eq.0) cycle
        do i=1,5
          nv=PYRAMIDS(n)%VertNo(i)
c  .......if 1st visit to point          
          if (POINTS(nv)%Idata(1).gt.0) then
            write(nfile,7001) POINTS(nv)%Rdata(1:3)
c  .........flip sign            
            POINTS(nv)%Idata(1) = -POINTS(nv)%Idata(1)
          endif
        enddo
      enddo
c
c  ...Step 3: write elements
      write(nfile,1005)ncell,ncoun
 1005 format('CELLS ',2(i10,2x))
c
      do n=1,NRPRISM
        ndom = PRISMS(n)%Domain
        call locate(ndom,Nlist,Nrd, nloc)
        if (nloc.eq.0) cycle
        do i=1,6
          nv=PRISMS(n)%VertNo(i)
          nvert(i) = -POINTS(nv)%Idata(1) - 1
        enddo
        write(nfile,7003) nvert(1:6)
 7003   format('6 ',4x,6(i6,2x)) 
      enddo
c
      do n=1,NRHEXAS
        ndom = HEXAS(n)%Domain
        call locate(ndom,Nlist,Nrd, nloc)
        if (nloc.eq.0) cycle
        do i=1,8
          nv=HEXAS(n)%VertNo(i)
          nvert(i) = -POINTS(nv)%Idata(1) - 1
        enddo
        write(nfile,7004) nvert(1:8)
 7004   format('8 ',4x,8(i6,2x)) 
      enddo
c
      do n=1,NRTETRA
        ndom = TETRAS(n)%Domain
        call locate(ndom,Nlist,Nrd, nloc)
        if (nloc.eq.0) cycle

        call locate_tetra(n, iflag)
        if (iflag.eq.0) cycle

        do i=1,4
          nv=TETRAS(n)%VertNo(i)
          nvert(i) = -POINTS(nv)%Idata(1) - 1
        enddo
        write(nfile,7005) nvert(1:4)
 7005   format('4 ',4x,4(i6,2x)) 
      enddo
c
      do n=1,NRPYRAM
        ndom = PYRAMIDS(n)%Domain
        call locate(ndom,Nlist,Nrd, nloc)
        if (nloc.eq.0) cycle
        do i=1,5
          nv=PYRAMIDS(n)%VertNo(i)
          nvert(i) = -POINTS(nv)%Idata(1) - 1
        enddo
        write(nfile,7006) nvert(1:5)
 7006   format('5 ',4x,5(i6,2x)) 
      enddo
c
      write(nfile,1006)ncell
 1006 format('CELL_TYPES ',i10)
 1007 format('13')
 1008 format('11')
 1009 format('10')
 1010 format('14')
      do i=1,ncell
        select case(ctype(i))
        case('pris');  write(nfile,1007)
        case('hexa');  write(nfile,1008)
        case('tetr');  write(nfile,1009)
        case('pyra');  write(nfile,1010)
        endselect
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
c
      end subroutine export_vol2vtk
c
c
c      
c----------------------------------------------------------------------
c
      subroutine locate_tetra(Nt, Iflag)
c     Iflag=1   iniside
c     Iflag=0   outside
c
      use GMP
#include "syscom.blk"
c
      dimension temp(3)      
c
      Iflag=0 ; temp=0.d0
      do i=1,4
        nv = TETRAS(Nt)%VertNo(i)
        temp = temp + POINTS(nv)%Rdata(1:3)
      enddo
      temp = temp/4.d0
      r = sqrt(temp(1)**2+temp(2)**2+temp(3)**2)
      if (r.lt.8.d0) Iflag=1
      write(*,*)'iflag = ',Iflag
c
c
      endsubroutine
