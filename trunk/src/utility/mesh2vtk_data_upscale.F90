#include "implicit_none.h"
subroutine mesh2vtk_data_upscale(Nout, Title, Isel, Z_harmonic, Show_Data )
  !
  use data_structure3D
  use element_data
  use upscale
  implicit none 
  !
  integer, intent(in)  :: Nout, Isel(10)
  VTYPE, intent(in) :: Z_harmonic
  character(len=*), intent(in) :: Title
  !
  type(vis) :: vis_obj
  character(len=4) :: type
  integer :: &
       i, mdle, &
       nr_vert, nr_elem

  interface
     subroutine Show_Data(Nout, Isel, Mdle, Z_harmonic)
        integer,                intent(in) :: Nout, Mdle
        integer, dimension(10), intent(in) :: Isel
        VTYPE,                  intent(in) :: Z_harmonic
     end subroutine Show_Data
  end interface
  
  ! Step 0: Compute nr_vert and nr_cell
  write(*,*) 'mesh2vtk:: Compupting nr verts and nr cells '
  mdle = 0; nr_vert = 0; nr_elem = 0;
  do i=1,NRELES
     call nelcon(mdle, mdle)
     type = NODES(mdle)%type
     vis_obj = vis_on_type(type)
     nr_vert = nr_vert + vis_obj%nr_vert
     nr_elem = nr_elem + vis_obj%nr_elem
  enddo

  ! Step 1: Header

!   select case(Isel(1))
!   case(0); write(Nout,5000) nr_vert ! point data
!   case(1); write(Nout,5001) nr_elem ! cell data
!   end select

! 5000 format('POINT_DATA ', i12)
! 5001 format('CELL_DATA ', i12)

  select case(Isel(2))
  case(0); 
     write(Nout,6000) title
     write(Nout,8000) 
  case(1); 
     write(Nout,6001) title
!     write(Nout,8000) 
  end select

6000 format('SCALARS  ',a16,' double 1')
6001 format('VECTORS  ',a16,' double ')

8000 format('LOOKUP_TABLE default')

  write(*,*) 'mesh2vtk:: Writing data '

  ! Step 2: Points Data
  mdle = 0
  do i=1,NRELES
     call nelcon(mdle, mdle)
     call Show_Data(Nout, Isel, Mdle, Z_harmonic)
  end do

  write(*,*) 'mesh2vtk:: Done '
end subroutine mesh2vtk_data_upscale

! Soleval some times behave very differently...
! Definitely data corruption exist. Fortran mystery...
!
#include "implicit_none.h"
subroutine dummy_soleval(Mdle)
  use data_structure3D
  use upscale
  implicit none
  ! ** Arguements                                                                       
  !----------------------------------------------------------------------               
  integer, intent(in) :: Mdle
  character(len=4) :: type
  !                                                                                     
  integer :: &
       i, iv,  &
       nodesl(27), norder(19), &
       norientl(27), nedge_orient(12), nface_orient(6)
  real*8 :: &
       xi(3), x(3), xnod(3,MAXbrickH), dxdxi(3,3)
  !                                                                                     
  VTYPE :: &
       zdofH(MAXEQNH,MAXbrickH),zdofE(MAXEQNE,MAXbrickE), &
       zdofV(MAXEQNV,MAXbrickV),zdofQ(MAXEQNQ,MAXbrickQ)
  !                                                                                     
  VTYPE :: &
       zsolH(MAXEQNH),  zgradH(MAXEQNH,3), &
       zsolE(3,MAXEQNE),zcurlE(3,MAXEQNE), &
       zsolV(3,MAXEQNV),zdivV(MAXEQNV), &
       zsolQ(MAXEQNQ)
  !                                                                                     
  VTYPE:: &
       zvalH(MAXEQNH), &
       zdvalH(MAXEQNH,3),zd2valH(MAXEQNH,3,3), &
       zvalE(3,MAXEQNE), &
       zdvalE(3,MAXEQNE,3),zd2valE(3,MAXEQNE,3,3), &
       zvalV(3,MAXEQNV), &
       zdvalV(3,MAXEQNV,3),zd2valV(3,MAXEQNV,3,3), &
       zvalQ(MAXEQNQ), &
       zdvalQ(MAXEQNQ,3),zd2valQ(MAXEQNQ,3,3)
    type(vis) :: vis_obj
    nodesl = 0; norientl = 0; norder = 0;
    nedge_orient = 0; nface_orient = 0;
    call elem_nodes(Mdle, nodesl, norientl)
    call find_orient_from_list(type, norientl, nedge_orient,nface_orient)
    call find_order_from_list(type, nodesl, norder)
    xnod = 0.d0; zdofH = ZERO; zdofE = ZERO;
    call nodcor(Mdle, xnod)
    call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
    type = NODES(Mdle)%type
    vis_obj = vis_on_type(type)
    do iv=1, vis_obj%nr_vert
       call get_vis_point(vis_obj, iv-1, xi)
       call soleval( &
            mdle, xi, &
            nedge_orient, nface_orient, norder, &
            xnod, zdofH,zdofE,zdofV,zdofQ, &
            1, x, dxdxi, &
            zsolH,zgradH,zsolE,zcurlE,zsolv,zdivV,zsolQ)
    end do
end subroutine dummy_soleval
