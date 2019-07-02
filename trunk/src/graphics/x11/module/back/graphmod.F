c----------------------------------------------------------------------
c
c   module name        - graphmod
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Oct 98
c
c   purpose            - module defines space for graphics
c
c
c----------------------------------------------------------------------
c
      module graphmod
c
      integer, allocatable, dimension(:) :: IGTRCU, IGTRNO, IGTR, IGSTR
c
      double precision, allocatable, dimension(:) :: RGTR, RGTRZ
c
c  ...dimensions
      integer :: MXIGTR,MXIGSTR,MXRGTRZ
c
c  ...explanation of variables:
c     RGTR   - real storage place for points coordinates
c              for graphics (in observer's system) for
c              triangles vertices
c     IGTRCU - information on triangle sides, color etc
c     IGTRNO - information about nodes numbers or small triangles
c     IGSTR  - information about colour of small triangles
c     RGTRZ  - depth-coordinates of triangles' mid-points
c     IGTR   - ordered list of triangles (back-to-front)
c     IGINV  - invisible blocks
c
      end module graphmod
