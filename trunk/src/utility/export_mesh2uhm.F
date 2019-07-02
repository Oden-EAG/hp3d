c---------------------------------------------------------------------
c
c   routine name       - export_mesh2uhm
c
c---------------------------------------------------------------------
c
c   latest revision    - Feb 2010
c
c   purpose            - export mesh to vtk legacy format
c
c   arguments :
c
c-----------------------------------------------------------------------
c
      subroutine export_mesh2uhm(Fname)
c
      use data_structure3D
      use element_data
#include "syscom.blk"
c
c  ...element nodes and their orientations
      dimension nverl(4),nodesl(27), norientl(27) 
      character(len=64) :: Fname
c
      nfile = 31
c
c  ...open file for writing
      open(unit=nfile, file=Fname, form='formatted')

c  ...initialize visit
      call reset_visit
c
c  ...Step 1 : visit all nodes 
      n_element = 0
      mdle = 0
      do iel=1, NRELES
         call nelcon(mdle, mdle)
         select case (NODES(mdle)%type)
         case ('mdln', 'tetr')
            n_element = n_element + 1
            call elem_nodes(mdle, nodesl, norientl)
            do ino=1,15
               NODES(nodesl(ino))%visit = 101
            end do
         end select
      end do
      
c     ...node counting
      n_node = 0
      do ino=1, NRNODS
         if (NODES(ino)%visit.eq.101) then
            n_node = n_node + 1
         endif
      enddo
c     
c     ...Step 2 : write nodes to file
      write(nfile,*) '# UHM format'
      write(nfile,*) 'UHM'
      write(nfile,*) '# number of nodes'
      write(nfile,*) n_node
c
      do ino=1, NRNODS
         if (NODES(ino)%visit.eq.101) then
            call find_ndof(ino, ndofH,ndofE,ndofV,ndofQ)
            write(nfile,*) ino, ' 0 ', ndofH, NODES(ino)%order, ' 0 '
         endif
      enddo
c
c     ...Step 3 : write element to file
      write(nfile,*) '# number of elements'
      write(nfile,*)  n_element
      
      mdle = 0
      do iel=1,NRELES
         call nelcon(mdle, mdle)
         select case (NODES(mdle)%type)
         case ('mdln','tetr')
            call elem_nodes(mdle, nodesl,norientl)
            write(nfile,*) ' 15  # ', mdle
            do ino=1,15
               write(nfile,*) nodesl(ino), ' 0 '
            end do
         end select
      enddo
c     
      close(nfile)         
      call reset_visit
c
      end subroutine
