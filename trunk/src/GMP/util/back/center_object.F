c---------------------------------------------------------------------
c      
c   routine name       - center_object
c
c---------------------------------------------------------------------
c
c   latest revision    - Jan 11
c
c   purpose            - routine centers an object around the origin
c
c---------------------------------------------------------------------
c
      subroutine center_object
c
      use GMP
#include "syscom.blk"
c
      dimension void(3)
c
c-----------------------------------------------------------------------
c
      iprint=1
c
c  ...Step 1: loop over points and determine dimensions of object      
      xmin  = 1.d30; xmax = -1.d30
      ymin  = 1.d30; ymax = -1.d30
      rzmin = 1.d30; rzmax = -1.d30
      do i = 1,NRPOINT 
        xmin  = min(POINTS(i)%Rdata(1),xmin )
        xmax  = max(POINTS(i)%Rdata(1),xmax )
        ymin  = min(POINTS(i)%Rdata(2),ymin )
        ymax  = max(POINTS(i)%Rdata(2),ymax )
        rzmin = min(POINTS(i)%Rdata(3),rzmin)
        rzmax = max(POINTS(i)%Rdata(3),rzmax)
      enddo
c  ...printing      
      if (iprint.eq.1) then
        write(*,*)'xmin,xmax = ',xmin,xmax
        write(*,*)'ymin,ymax = ',ymin,ymax
        write(*,*)'zmin,zmax = ',rzmin,rzmax
        call pause
      endif
c
c  ...Step 2: compute shift vector
      void(1) = (xmin  + xmax )/2.d0      
      void(2) = (ymin  + ymax )/2.d0      
      void(3) = (rzmin + rzmax)/2.d0     
c      
c  ...Step 3: shift object
c
c     P  O  I  N  T  S
      do i = 1,NRPOINT
        POINTS(i)%Rdata(1:3) = POINTS(i)%Rdata(1:3) - void(1:3)
      enddo
c      
c     C U R V E S
      do i = 1,NRCURVE
        if (CURVES(i)%Type.ne.'5Bezier') cycle     
        do j = 0,5
          CURVES(i)%Rdata(3*j:3*j+2) = CURVES(i)%Rdata(3*j:3*j+2) -
     .                                 void(1:3)
        enddo
      enddo 
c      
c     T R I A N G L E S
      do i = 1,NRTRIAN
        if (TRIANGLES(i)%Type.ne.'G1RecTri') cycle
        do j = 0,35
        TRIANGLES(i)%Rdata(3*j:3*j+2) = TRIANGLES(i)%Rdata(3*j:3*j+2) -
     .                                  void(1:3)
        enddo
      enddo 
c
c
      end
