!----------------------------------------------------------------------------------------
!> Purpose : write geometry to .h5 file
!!
!> @param[in ] Sname -  
!> @param[in ] Sfile - geometry file
!> @param[out] IcE   -
!> @param[out] IcN   -
!> @param[out] IcP   -
!!
!> @date Apr 2019
!----------------------------------------------------------------------------------------
!
#include "implicit_none.h"
!
subroutine geom2vtk(Sname,Sfile, IcE,IcN,IcP)
!
      use data_structure3D
      use upscale
      use paraview , only : PARAVIEW_GEOM,PARAVIEW_EXGEOM,PARAVIEW_ISOGEOM,PARAVIEW_DOMAIN
!
      implicit none
!
      character(len=*),  intent(in ) :: Sname, Sfile
      integer,           intent(out) :: IcE, IcN, IcP
!
      type(vis) :: vis_obj
!
      real*8,dimension(  MAXbrickH) :: shapH
      real*8,dimension(3,MAXbrickH) :: gradH
      character(len=4) :: etype
      integer :: i, iv, nV, iel, ioffs, mdle, nrdofH, j,k, ndom
      integer :: nodesl(27), norder(19), nverl(8)
      integer :: norientl(27), nedge_orient(12), nface_orient(6)
!
      real*8 :: xi(3), x(3), xnod(3,MAXbrickH), dxdxi(3,3)
!
!   ..OpenMP parallelization: auxiliary variables
      integer, dimension(NRELES) :: n_elem, n_vert_offset, n_elem_vert
      real*8, allocatable :: val_array(:,:)
!
!----------------------------------------------------------------------------------------
!
!     Step 0 : Clear vis geometry
      call vis_geometry_clear
!
!     Step 1 : Points

!   ..create list of mdle node numbers,
!     and count number of visualization points
      mdle=0; IcP=0
      do i=1,NRELES
        call nelcon(mdle, mdle)
        n_elem(i) = mdle
        if (PARAVIEW_DOMAIN.ne.0) then
          call find_domain(mdle, ndom)
          if (ndom.ne.PARAVIEW_DOMAIN) cycle
        endif
        etype = NODES(mdle)%type
        vis_obj = vis_on_type(etype)
        n_vert_offset(i) = IcP
        n_elem_vert(i) = vis_obj%nr_vert
        IcP = IcP + vis_obj%nr_vert
      enddo
      allocate(val_array(3,IcP))
!
!$OMP PARALLEL DO                               &
!$OMP PRIVATE(mdle,ndom,etype,iv,xi,nV,         &
!$OMP         norder,nedge_orient,nface_orient, &
!$OMP         xnod,x,dxdxi,nrdofH,shapH,gradH)  &
!$OMP SCHEDULE(DYNAMIC)
      do i=1,NRELES
        mdle = n_elem(i)
!
!       SKIP IF WRONG DOMAIN
        if (PARAVIEW_DOMAIN.ne.0) then
          call find_domain(mdle, ndom)
          if (ndom.ne.PARAVIEW_DOMAIN) cycle
        endif
!
!       gdofs, solution dofs, order, orientations
        call nodcor(mdle, xnod)
        call find_elem_nodes(mdle, norder,nedge_orient,nface_orient)
!
!       select appropriate visualization object
        etype = NODES(mdle)%type
        nV = n_elem_vert(i)
!
        do iv=1,nV
!
!         visualization point accounting for VLEVEL
          call get_vis_point(vis_on_type(etype),iv-1, xi)
!
!         H1 shape functions at visualization point
          call shape3H(etype,xi,norder,nedge_orient,nface_orient, nrdofH,shapH,gradH)
!
!         geometry map
          select case(PARAVIEW_GEOM)
!
!         -- ISOPARAMETRIC GEOMETRY MAP --
          case(PARAVIEW_ISOGEOM)
            x(1:3)=0.d0
            do k=1,nrdofH
              x(1:3) = x(1:3) + xnod(1:3,k)*shapH(k)
            enddo
!
!         -- EXACT GEOMETRY MAP --
          case(PARAVIEW_EXGEOM)
            call exact_geom(mdle,xi, x,dxdxi)
          endselect
!
!         add point to visualization
          val_array(1:3,n_vert_offset(i)+iv) = x
        enddo
      enddo
!$OMP END PARALLEL DO
!
      do i=1,IcP
        call vis_geometry_add_point(val_array(1:3,i));
      enddo
      deallocate(val_array)
!
!     Step 2 : Elements
      IcE=0; IcN=0
!
      do i=1,NRELES
        mdle = n_elem(i)
!
!       SKIP IF WRONG DOMAIN
        if (PARAVIEW_DOMAIN.ne.0) then
          call find_domain(mdle, ndom)
          if (ndom.ne.PARAVIEW_DOMAIN) cycle
        endif
!        
        call elem_nodes(mdle, nodesl,norientl)
!    
!       select appropriate visualization object  
        etype = NODES(mdle)%type
        vis_obj = vis_on_type(etype)
    
        do iel=1, vis_obj%nr_elem
!
!         visualization element accounting for VLEVEL
          call get_vis_elem(vis_obj,iel,n_vert_offset(i), nverl)
!
!         add element to visualization
          call vis_geometry_add_object(ivis_type(etype), nverl)
          IcE = IcE + 1
          IcN = IcN + nvert(etype)
        enddo
      enddo
!
!     Step 3 : write to file
      call vis_geometry_write(Sname,len(Sname),Sfile,len(Sfile))
!    
!     Step 4 : Clear vis geometry
      call vis_geometry_clear
!
!
end subroutine geom2vtk
