      module vtk_interf
c
      integer, private :: NGEOM,NSOLN
c
      integer, private :: NCOMP
c
      contains
c
c=======================================================================
c
      subroutine init_vtk_geom
      NGEOM = 25
      open(unit=NGEOM,
     .     file='VTK/geom',
     .     form='formatted',
     .     access='sequential',
     .     status='unknown')
      end subroutine init_vtk_geom
c
c=======================================================================
c
      subroutine delete_vtk_geom
      close(NGEOM)
      end subroutine delete_vtk_geom
c
c=======================================================================
c
      subroutine init_vtk_soln(Nrcomp)
      NSOLN = 26
      open(unit=NSOLN,
     .     file='VTK/soln',
     .     form='formatted',
     .     access='sequential',
     .     status='unknown')
      write(NSOLN,*) C_MODE
      write(NSOLN,*) Nrcomp
      NCOMP = Nrcomp
      end subroutine init_vtk_soln
c
c=======================================================================
c
      subroutine delete_vtk_soln
      close(NSOLN)
      end subroutine delete_vtk_soln
c
c=======================================================================
c
      subroutine set_nrfig(Isol,Nrfig)
c
c     output unit
      if (Isol==0) then
        nunit = NGEOM
      else
        nunit = NSOLN
      endif
c
      write(nunit,*) Nrfig
c
      end subroutine set_nrfig
c
c=======================================================================
c
      subroutine add_hexa(Isol,X,Z,Nsub,Num)
#include "syscom.blk"
c
      dimension X(:,:),Z(:,:)
c
c     output unit
      if (Isol==0) then
        nunit = NGEOM
      else
        nunit = NSOLN
      endif
c
c     type (hexa)
      write(nunit,*) 1
c
c     number
      write(nunit,*) Num
c
c     number of subdivisions
      write(nunit,*) Nsub
c
c     number of points
      nrpt = (Nsub+1)**3
c
c     point coordinates
      do j=1,nrpt
      do i=1,3
        write(nunit,*) X(i,j)
      enddo
      enddo
c
c     solution values
      if (Isol==1) then
        do j=1,nrpt
        do i=1,NCOMP
          write(nunit,*) dreal_part(Z(i,j)),dimag_part(Z(i,j))
        enddo
        enddo
      endif
c
      end subroutine add_hexa
c
c=======================================================================
c
      subroutine add_tetra(Isol,X,Z,Nsub,Num)
#include "syscom.blk"
c
      dimension X(:,:),Z(:,:)
c
c     output unit
      if (Isol==0) then
        nunit = NGEOM
      else
        nunit = NSOLN
      endif
c
c     type (tetra)
      write(nunit,*) 2
c
c     number
      write(nunit,*) Num
c
c     number of subdivisions
      write(nunit,*) Nsub
c
c     number of points
      nrpt = (Nsub+1)*(Nsub+2)*(Nsub+3)/6
c
c     point coordinates
      do j=1,nrpt
      do i=1,3
        write(nunit,*) X(i,j)
      enddo
      enddo
c
c     solution values
      if (Isol==1) then
        do j=1,nrpt
        do i=1,NCOMP
          write(nunit,*) dreal_part(Z(i,j)),dimag_part(Z(i,j))
        enddo
        enddo
      endif
c
      end subroutine add_tetra
c
c=======================================================================
c
      subroutine add_prism(Isol,X,Z,Nsub,Num)
#include "syscom.blk"
c
      dimension X(:,:),Z(:,:)
c
c     output unit
      if (Isol==0) then
        nunit = NGEOM
      else
        nunit = NSOLN
      endif
c
c     type (prism)
      write(nunit,*) 3
c
c     number
      write(nunit,*) Num
c
c     number of subdivisions
      write(nunit,*) Nsub
c
c     number of points
      nrpt = (Nsub+1)**2*(Nsub+2)/2
c
c     point coordinates
      do j=1,nrpt
      do i=1,3
        write(nunit,*) X(i,j)
      enddo
      enddo
c
c     solution values
      if (Isol==1) then
        do j=1,nrpt
        do i=1,NCOMP
          write(nunit,*) dreal_part(Z(i,j)),dimag_part(Z(i,j))
        enddo
        enddo
      endif
c
      end subroutine add_prism
c
c=======================================================================
c
      subroutine add_pyramid(Isol,X,Z,Nsub,Num)
#include "syscom.blk"
c
      dimension X(:,:),Z(:,:)
c
c     output unit
      if (Isol==0) then
        nunit = NGEOM
      else
        nunit = NSOLN
      endif
c
c     type (pyramid)
      write(nunit,*) 4
c
c     number
      write(nunit,*) Num
c
c     number of subdivisions
      write(nunit,*) Nsub
c
c     number of points
      nrpt = (Nsub+1)*(Nsub+2)*(2*Nsub+3)/6
c
c     point coordinates
      do j=1,nrpt
      do i=1,3
        write(nunit,*) X(i,j)
      enddo
      enddo
c
c     solution values
      if (Isol==1) then
        do j=1,nrpt
        do i=1,NCOMP
          write(nunit,*) dreal_part(Z(i,j)),dimag_part(Z(i,j))
        enddo
        enddo
      endif
c
      end subroutine add_pyramid
c
c=======================================================================
c
      end module vtk_interf
