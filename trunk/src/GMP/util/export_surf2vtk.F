c----------------------------------------------------------------------
c
c   routine name       - export_surf2vtk
c
c----------------------------------------------------------------------
c
c   latest revision    - Jan 11
c
c   purpose            - export surface mesh to vtk format
c
c   arguments:
c                in:   - Fname: file name      
c                      - Nrs  : number of surfaces to be exported
c                      - Nlist: list of surfaces to be exported
c
c-----------------------------------------------------------------------
c
      subroutine export_surf2vtk(Fname,Nrs,Nlist)
c
      use error
      use GMP
#include "syscom.blk"
c
c  ...file name      
      character(len=64) :: Fname
c  ...list of surfaces to be exported      
      dimension Nlist(Nrs)
c  ...vertices of surface triangle in the new enumeration
      dimension nvert(3)
c  ...shift factor      
      dimension rshift(3)
c  ...triangle surface number      
      dimension nsurf(NRTRIAN)
c
      nfile = 31
c
c     customiziation flag
c     0 - no customization
c     1 - inner ear
c     2 - phantom
c     3 - truncating sphere
      icust=0
c
c  ...set up stretch and shift factors      
      select case(icust)
c  ...inner ear      
      case(1)
        rstret = 1.d0
        rshift = (/-1.d0,0.d0,0.d0/)
c  ...truncating sphere        
      case(3)
        rstret = 16.5d0
        rshift = (/0.d0,-3.3d0,0.d0/)
      case default        
      endselect      
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
c  ...Step 1: loop over surface triangles, count and re-number vertices
      n=0
      do i=1,NRTRIAN
c  .....cycle if not on a surface      
        if (TRIANGLES(i)%Type.ne.'PTITri') cycle
        ns = TRIANGLES(i)%Idata(1)
        call locate(ns,Nlist,Nrs, nloc)
c  .....cycle if not on a desired surface        
        if (nloc.eq.0) cycle
c
c......................................................................        
c       C U S T O M I Z A T I O N        
          if (icust.eq.1) then      
            if (ns.eq.9) then
              call locate_trian(i, iflag)
              write(*,*)'i,iflag = ',i,iflag
              if (iflag.eq.1) cycle
            endif
          endif
c......................................................................        
c
        do j=1,3
          nv=TRIANGLES(i)%VertNo(j)
c  .......if 1st visit to point          
          if (.not.associated(POINTS(nv)%Idata)) then
            n=n+1
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
c  ...Step 2: loop over surface triangles, count surface triangles
c             and write point's coordinates
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
      n=0 ; nsurf=0
      do i=1,NRTRIAN
c  .....cycle if not on a surface      
        if (TRIANGLES(i)%Type.ne.'PTITri') cycle
        ns = TRIANGLES(i)%Idata(1)
        call locate(ns,Nlist,Nrs, nloc)
c  .....cycle if not on a desired surface        
        if (nloc.eq.0) cycle
c
c......................................................................        
c       C U S T O M I Z A T I O N         
        if (icust.eq.1) then
          if (ns.eq.9) then
            call locate_trian(i, iflag)
            write(*,*)'i,iflag = ',i,iflag
            if (iflag.eq.1) cycle
          endif
        endif
c......................................................................        
c
        n=n+1
        nsurf(n)=ns
        do j=1,3
          nv=TRIANGLES(i)%VertNo(j)
c  .......if 1st visit to point          
          if (POINTS(nv)%Idata(1).gt.0) then
            write(nfile,7001) rstret*POINTS(nv)%Rdata(1:3)+rshift(1:3)
 7001       format(3(f20.12,2x))     
c  .........flip sign            
            POINTS(nv)%Idata(1) = -POINTS(nv)%Idata(1)
          endif
        enddo
      enddo
c
c  ...Step 3: write elements
      write(nfile,1005)n,n*4
 1005 format('CELLS ',2(i10,2x))
      do i=1,NRTRIAN
c  .....cycle if not on a surface      
        if (TRIANGLES(i)%Type.ne.'PTITri') cycle
        ns = TRIANGLES(i)%Idata(1)
        call locate(ns,Nlist,Nrs, nloc)
c  .....cycle if not on a desired surface        
        if (nloc.eq.0) cycle
c
c......................................................................        
c       C U S T O M I Z A T I O N         
        if (icust.eq.1) then
          if (ns.eq.9) then   
            call locate_trian(i, iflag)
            if (iflag.eq.1) cycle
          endif
        endif
c......................................................................        
c        
c  .....collect vertices        
        do j=1,3
          nv=TRIANGLES(i)%VertNo(j)
          nvert(j) = -POINTS(nv)%Idata(1) - 1
        enddo
        write(nfile,7003) nvert(1:3)
 7003   format('3 ',4x,3(i6,2x))       
      enddo
c
c  ...Step 4: write elements type      
      write(nfile,1006)n
 1006 format('CELL_TYPES ',i10)
      do i=1,n
        write(nfile,1007)
 1007  format('5')
      enddo
c
c  ...Step 5: write elemenets data (suraface number)      
      write(nfile,1008)n
 1008 format('CELL_DATA ',i10)     
      write(nfile,1009)
 1009 format('SCALARS  int  1')
      write(nfile,1010)
 1010 format('LOOKUP_TABLE default')
 1011 format(i1)
 1012 format(i2)
 1013 format('14')
 1014 format('15')
      do i=1,n
        select case(icust) 
c  .....phantom        
        case(2) ; write(nfile,1013)
c  .....truncating sphere          
        case(3) ; write(nfile,1014)
        case default
          if (nsurf(i).lt.10) then
            write(nfile,1011)nsurf(i)
          else
            write(nfile,1012)nsurf(i)
          endif
        endselect
      enddo
c
c  ...Step 6: deallocate Idata
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
      end subroutine
c
c      
c----------------------------------------------------------------------
c
      subroutine locate_trian(Nt, Iflag)
c
      use GMP
#include "syscom.blk"
c
      Iflag=0
      do i=1,3
        nv=TRIANGLES(Nt)%VertNo(i)
        r=sqrt(POINTS(nv)%Rdata(2)**2+POINTS(nv)%Rdata(3)**2)
        write(*,*)'r = ',r
        if (r.gt.0.851d0) Iflag = 1
      enddo
c
c
      endsubroutine
