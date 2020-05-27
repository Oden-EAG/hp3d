#include "implicit_none.h"
!
!
subroutine mesh2vtk_geom_upscale(Nout, Isel, Is)
!
      use data_structure3D
      use element_data
      use upscale
!
      implicit none
      integer, intent(in) :: Nout, Isel(10), Is
!
!     object containing all data for visualization
!
!
      type(vis) :: vis_obj
!
      character(len=4) :: type
!
       integer :: &
       i, iv, iel, ioffs, mdle, ndom, &
       nodesl(27), norder(19), &
       norientl(27), nedge_orient(12), nface_orient(6), nverl(8), &
       nr_vert, nr_elem, nr_size

       real(8) :: xi(3), x(3), xnod(3,MAXbrickH), dxdxi(3,3)

       VTYPE :: zdofH(MAXEQNH,MAXbrickH),zdofE(MAXEQNE,MAXbrickE), &
                zdofV(MAXEQNV,MAXbrickV),zdofQ(MAXEQNQ,MAXbrickQ)
!
       VTYPE :: zsolH(  MAXEQNH),zgradH(  MAXEQNH,3), &
                zsolE(3,MAXEQNE),zcurlE(3,MAXEQNE  ), &
                zsolV(3,MAXEQNV),zdivV(   MAXEQNV  ), &
                zsolQ(  MAXEQNQ)
!
!-------------------------------------------------------------------------------------
!
!     Step 0 : Compute nr_vert and nr_cell
!
      nr_vert=0 ; nr_elem=0 ; nr_size=0
!
      mdle=0
      do i=1,NRELES
        call nelcon(mdle, mdle)
!
        call find_domain(mdle, ndom)
        if (Isel(ndom+Is) == 0)  cycle
!
!       select appropriate visualization object (TETR_VIS, PRIS_VIS, HEXA_VIS)
        type = NODES(mdle)%type
        vis_obj = vis_on_type(type)
!
        nr_vert = nr_vert + vis_obj%nr_vert
        nr_elem = nr_elem + vis_obj%nr_elem
        nr_size = nr_size + vis_obj%nr_elem*(nvert(type)+1)
      enddo
!
!     Step 1 : Print Header
      write(Nout,6000)
      write(Nout,6001)
      write(Nout,6002)
      write(Nout,6003)
      write(Nout,6004) nr_vert
 6000 format('# vtk DataFile Version 2.0')
 6001 format('hp3d export mesh ')
 6002 format('ASCII')
 6003 format('DATASET UNSTRUCTURED_GRID')
 6004 format('POINTS ',i12, '  double')
!
!     Step 2 : Points
      mdle=0
      do i=1,NRELES
        call nelcon(mdle, mdle)
!
        call find_domain(mdle, ndom)
        if (Isel(ndom+Is) == 0)  cycle
!
        type = NODES(mdle)%type
        call elem_nodes(mdle, nodesl, norientl)
        call find_order_from_list(type, nodesl, norder)
        call find_orient_from_list(type, norientl, nedge_orient,nface_orient)
!
        xnod = 0.d0
        call nodcor(mdle, xnod)
        call solelm(mdle, zdofH,zdofE,zdofV,zdofQ)
!
!       select appropriate visualization object (TETR_VIS, PRIS_VIS, HEXA_VIS)
        type = NODES(mdle)%type
        vis_obj = vis_on_type(type)
!
!       loop over vertices of visualization object
        do iv=1,vis_obj%nr_vert
!
!         vertex coordinates (indexing starts at 0)
          call get_vis_point(vis_obj,iv-1, xi)
!
!         evaluate geometry and solution
          call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod, &
               zdofH,zdofE,zdofV,zdofQ,0,x,dxdxi,zsolH,zgradH,zsolE,zcurlE,zsolv,zdivV,zsolQ)
!
!         print coordinates in physical space
          write(Nout,6005) x(1:3)
 6005     format('   ', 3(f20.15,2x))
        enddo
      enddo
!
!     Step 4 : Cell header
      write(Nout,6006) nr_elem, nr_size
 6006 format('CELLS ', i12, i12)
!
!     Step 5 : Write vis elem
      mdle=0 ; ioffs=0
      do i=1,NRELES
         call nelcon(mdle, mdle)
!
         call find_domain(mdle, ndom)
         if (Isel(ndom+Is) == 0)  cycle
!
!        select appropriate visualization object (TETR_VIS, PRIS_VIS, HEXA_VIS)
         type = NODES(mdle)%type
         vis_obj = vis_on_type(type)
!
!        loop over elements of visualization object
         do iel=1,vis_obj%nr_elem
!
!          element vertices
           call get_vis_elem(vis_obj,iel,ioffs, nverl)
!
!          print vertices
           write(Nout,6007) nvert(type), nverl(1:nvert(type))
 6007      format(i10, 8i10)
         enddo
!
!        update offset for vertex enumeration
         ioffs = ioffs + vis_obj%nr_vert
      enddo
!
!     Step 6 : VTK cell type
      write(Nout,6008) nr_elem
 6008 format('CELL_TYPES ', i12)

      mdle=0
      do i=1,NRELES
        call nelcon(mdle, mdle)
!
        call find_domain(mdle, ndom)
        if (Isel(ndom+Is) == 0)  cycle
!
!       select appropriate visualization object (TETR_VIS, PRIS_VIS, HEXA_VIS)
        type = NODES(mdle)%type
        vis_obj = vis_on_type(type)
!
!       loop over elements of visualization object
        do iel=1,vis_obj%nr_elem
!
!         print cell type
          write(Nout,6009) ivtk_type(type)
 6009     format(i10)
        enddo
      enddo

      write(Nout,5000) nr_vert ! point data
 5000 format('POINT_DATA ', i12)
!5001 format('CELL_DATA ', i12)

end subroutine mesh2vtk_geom_upscale
