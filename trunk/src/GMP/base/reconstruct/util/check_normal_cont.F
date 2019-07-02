c----------------------------------------------------------------------
c
c   routine name       - check_normal_cont
c
c----------------------------------------------------------------------
c
c   latest revision    - Apr 10
c
c   purpose            - routine performs G^1 reconstruction of
c                        geometry for a triangulated surface
c
c   arguments          - none   
c
c----------------------------------------------------------------------
c
      subroutine check_normal_cont
c
      use GMP
      use control
c
#include "syscom.blk"

      integer, parameter :: nsub = 3
c      
      dimension eta(2),void(3),dxdeta(3,2)
      dimension rnorm1(1:3,0:nsub),rnorm2(1:3,0:nsub)
      dimension nt_save(2),nedg_save(2),nblock_save(2)
      dimension ncurve_fail(100)
c
c
      iprint=0
c      
      nfail = 0
c      
c  ...loop through all curves
      do nc = 1,NRCURVE
c  .....skip if not a quintic Bezier curve
        if ((CURVES(nc)%Type.ne.'5Bezier').or.
     .      (CURVES(nc)%Type.ne.'7Bezier')) cycle
c  .....skip if curve is a sharpe edge
        call check_edge(nc, idec); if (idec.eq.2) cycle
ccc        write(*,8062) nc
ccc 8062   format(' reconstruct: checking normal along nc = ',i5)
c  .....get the beginning point
        np = CURVES(nc)%EndPoNo(1)
c  .....loop through connected triangles
        ifound = 0
        do it = 1,CURVES(nc)%NrFig
          call decode(abs(CURVES(nc)%FigNo(it)), nt,lab)
c  .......cycle if not a triangle          
          if (lab.ne.1) cycle
c  .......cycle if not a G1RecTri          
          if (TRIANGLES(nt)%Type.ne.'G1RecTri') cycle
c  .......loop through edges of the triangle
          do ie = 1,3
            nc1 = iabs(TRIANGLES(nt)%EdgeNo(ie))
            if (nc1.ne.nc) cycle
            ifound = ifound + 1
            nt_save(ifound) = nt; nedg_save(ifound) = ie
c  .........loop through adjacent blocks            
            do is = 1,2
              if (TRIANGLES(nt)%BlockNo(is).ne.0) then
                nblock_save(ifound) = TRIANGLES(nt)%BlockNo(is)
              endif
            enddo
c
c  for each subdivision point collect normals   
c            
c  .........loop over subdivision points            
            do i = 0,nsub
              xi = i*1.d0/nsub
              select case(ie)
              case(1)
                if (TRIANGLES(nt)%EdgeNo(ie).lt.0) xi = 1.d0 - xi
                eta(1) = xi; eta(2) = 0.d0
              case(2)
                if (TRIANGLES(nt)%EdgeNo(ie).lt.0) xi = 1.d0 - xi
                eta(1) = 1.d0 - xi; eta(2) = xi
              case(3)
                if (TRIANGLES(nt)%EdgeNo(ie).lt.0) xi = 1.d0 - xi
                eta(1) = 0.d0; eta(2) = xi
              end select
              call trian(nt,eta, void,dxdeta)
              call cross_product(dxdeta(1:3,1),dxdeta(1:3,2), void)
              call normalize(void)
c  ...........use 1st subdivision point to determine triangle orientation              
              if (i.eq.0) then
                select case (POINTS(np)%Type)
                case ('CoorNrm')
                  call scalar_product(void,POINTS(np)%Rdata(4:6), s)
                case ('SharpPt')
                  n = size(POINTS(np)%Idata)
                  iflag = 0
                  do j = 1,n
                    call decode(POINTS(np)%Idata(j), it,lab)
                    if (it.ne.nt) cycle
                    select case(lab)
                    case (1)
                      call scalar_product(void,POINTS(np)%Rdata(4:6), s)
                      iflag = 1
                    case (2)
                      call scalar_product(void,POINTS(np)%Rdata(7:9), s)
                      iflag = 1
                    case default
                      write(*,*)'check_normal_cont: inconsistent label,
     . lab = ',lab
                      stop
                    endselect
                  enddo
                  if (iflag.eq.0) then
                    write(*,*)'check_normal_cont: topology 
     . inconsistency!'
                    write(*,*)'np,it,Idata = ',np,it,POINTS(np)%Idata
                    call print_GMP
                    stop
                  endif        
                case default
                  write(*,*)'check_normal_cont: inconsistent point 
     . type!'
                  write(*,*)'np,type = ',np,POINTS(np)%Type
                  stop
                endselect  
                factor = 1.d0
c  .............sign function returns 1 if s > 0, -1 otherwise                
                rsign = sign(factor,s)
              endif
              if (iprint.eq.1) then
         write(*,2001)ifound,nt_save(ifound),nedg_save(ifound),xi,rsign,
     . void
 2001    format(' ifound = ',i1,'; nt = ',i5,'; ie = ',i1,'; xi = ',
     . e12.5,' --> sign = ',e12.5,'; nor = ',3(e12.5,2x))
              endif
              select case(ifound)
              case(1)
                rnorm1(1:3,i) = rsign*void(1:3)
              case(2)
                rnorm2(1:3,i) = rsign*void(1:3)
              case default
                write(*,*)'check_normal_cont: INCONSISTENCY '
                stop
              end select
c  .........end of loop over subdivision points              
            enddo
c  .......end of loop over triangle edges            
          enddo
c  .....end of loop over connected triangles          
        enddo
c  .....printing
        if (iprint.eq.1) then
          do i = 0,nsub
            write(*,2045)rnorm1(1:3,i)
 2045       format('   rnorm1 = ',3(e12.5,2x))
            write(*,2046)rnorm2(1:3,i)
 2046       format('   rnorm2 = ',3(e12.5,2x))
            enddo
            call pause
        endif 
c
        iflag = 0 
c  for each subdivision point compare normals        
c  .....loop over subdivisions        
        do i = 0,nsub
          if (iprint.eq.1) then
            write(*,2045)rnorm1(1:3,i)
            write(*,2046)rnorm2(1:3,i)
            call pause
          endif
c  .......loop over components        
          do j = 1,3
c  .........if geometrical tollerance is exceeded          
            if (abs(rnorm2(j,i)-rnorm1(j,i)).gt.GEOM_TOL) then
c  ...........raise flag
              iflag = 1                    
              write(*,3334)nc,i,rnorm1(1:3,i)     
 3334         format(' nc = ',i8,'; sub = ',i2,'; norm1 = ',3(e12.5,2x))
              write(*,3335)rnorm2(1:3,i)     
 3335         format('                          norm2 = ',3(e12.5,2x))
cc              write(*,3333)nc,i,j,abs(rnorm2(j,i)-rnorm1(j,i))
cc 3333         format(' nc = ',i3,'; sub = ',i2,'; comp = ',i1,'; error
cc     .  = ',e12.5)
cc              write(*,*) 'ERROR IN COMP j = ',j,
cc     .                    abs(rnorm2(j,i)-rnorm1(j,i))
cc              write(*,8056) np,POINTS(np)%Rdata(4:6)
cc 8056         format('np,rn = ',i5,2x,3f8.3)
cc              write(*,8019) nt_save(1:2),nedg_save(1:2),nblock_save(1:2)
cc 8019         format('rectangles,edges,blocks = ',3(2i7,2x))
cc              write(*,8020) nc,i,rnorm1(1:3,i),rnorm2(1:3,i)
cc 8020         format('reconstruct: nc,i,rnorm1(1:3,i),
cc     . rnorm2(1:3,i) = ',i4,i2,2x,3f8.3,2x,3f8.3)
cc              write(*,8021) CURVES(nc)%EndPoNo(1:2)
cc 8021         format('             endpoints = ',2i5)
cc              call pause
cc              go to 50
            endif
c  .......end of loop over components            
          enddo
c  .....end of loop over subdivision points          
        enddo
c  .....increment counter and save curve number
        if (iflag.eq.1) then
          nfail = nfail + 1; ncurve_fail(nfail) = nc
        endif
c 50     continue
cc      call pause  
c  ...end of loop through curves
      enddo
c      
      if (nfail.eq.0) then
        write(*,*)'G1 OBJECT!!!'
      else        
        write(*,*)'=================================================='
        write(*,555)NRCURVE,nfail
 555    format(' NRCURVE = ',i10,'; nfail = ',i10)       
        do i = 1,nfail
        write(*,557)i,ncurve_fail(i),CURVES(ncurve_fail(i))%EndPoNo(1:2)
 557      format(' i = ',i2,'; nc = ',i7,', endpoints = ',i8,', ',i8)
        enddo     
        write(*,*)'=================================================='
      endif        
c
c
      end
