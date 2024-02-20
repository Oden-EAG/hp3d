#include "typedefs.h"
!----------------------------------------------------------------------------------------
!> @brief write geometry to .h5 file
!!
!> @param[in ] Sname - name/description (e.g., 'Geometry')
!> @param[in ] Sfile - geometry file name (e.g., ../outputs/paraview/geom_00000.h5)
!> @param[out] IcE   - Number of all vis object subelements in the active mesh
!> @param[out] IcN   - Number of all vis object subelement vertices in the active mesh
!> @param[out] IcP   - Number of points (coords) of all vis objects in the active mesh
!!
!> @date Sep 2023
!----------------------------------------------------------------------------------------
subroutine geom2vtk(Sname,Sfile, IcE,IcN,IcP)
!
   use data_structure3D
   use upscale
   use paraview
   use mpi_wrapper
   use par_mesh      , only: DISTRIBUTED,HOST_MESH
!
   implicit none
!
   character(len=*),  intent(in ) :: Sname, Sfile
   integer,           intent(out) :: IcE, IcN, IcP
!
   type(vis) :: vis_obj
!
   real(8) :: shapH(  MAXbrickH)
   real(8) :: gradH(3,MAXbrickH)
   integer :: ntype
   integer :: iv, nV, ico, iel, mdle, nrdofH, i, j, k, l, ivis, ndom
   integer :: nodesl(27), norder(19), nverl(27)
   integer :: norientl(27), nedge_orient(12), nface_orient(6)
!
   real(8) :: xi(3), x(3), xnod(3,MAXbrickH), dxdxi(3,3)
!
!..auxiliary variables for VTU output
   integer, allocatable       :: VTU_element_type_offset(:)
   integer                    :: m
!
!..OpenMP parallelization: auxiliary variables
   integer, dimension(NRELES) :: n_vert_offset,n_obj_offset,n_elem_vert
!
!..MPI
   integer :: ierr,count,subd,ice_subd,icn_subd
!
!----------------------------------------------------------------------------------------
!
!..Step 1 : Preliminary calculations (offsets, etc.)
!
!..create list of mdle node numbers,
!  and count number of visualization points
   IcP=0; ico=0; m = 0;
!
   if (VIS_VTU) then
      allocate(VTU_element_type_offset(NRELES))
   endif
!
   do iel=1,NRELES
      mdle = ELEM_ORDER(iel)
      if (PARAVIEW_DOMAIN.ne.0) then
         call find_domain(mdle, ndom)
         if (ndom.ne.PARAVIEW_DOMAIN) cycle
      endif
      ntype = NODES(mdle)%ntype
      vis_obj = vis_on_type(ntype)
      n_vert_offset(iel) = IcP
      n_obj_offset(iel) = ico
      n_elem_vert(iel) = vis_obj%NR_VERT
!
      if (VIS_VTU) VTU_element_type_offset(iel) = m
!
      IcP = IcP + vis_obj%NR_VERT
!
      if (SECOND_ORDER_VIS) then
         ico = ico + vis_obj%NR_ELEM*(nobj_conf(ntype))
         m = m + 1
      else
         ico = ico + vis_obj%NR_ELEM*(nobj_conf(ntype)+1)
         m = m + vis_obj%NR_ELEM
      endif
   enddo
!
   call geometry_init(ico,Icp)
   GEOM_PTS(1:3,1:Icp) = 0.d0; GEOM_OBJ(1:ico) = 0
!
   ice_subd = 0
!
!$OMP PARALLEL
!
!..Step 2 : Points
!
!$OMP DO                                        &
!$OMP PRIVATE(mdle,ndom,ntype,iv,xi,nV,subd,    &
!$OMP         norder,nedge_orient,nface_orient, &
!$OMP         xnod,x,dxdxi,nrdofH,shapH,gradH)  &
!$OMP SCHEDULE(DYNAMIC)
   do iel=1,NRELES
      mdle = ELEM_ORDER(iel)
!
!  ...SKIP IF WRONG DOMAIN
      if (PARAVIEW_DOMAIN.ne.0) then
         call find_domain(mdle, ndom)
         if (ndom.ne.PARAVIEW_DOMAIN) cycle
      endif
!
      if (DISTRIBUTED) then
         call get_subd(mdle, subd)
         if (RANK .ne. subd) cycle
      endif
!
!  ...gdofs, solution dofs, order, orientations
      call nodcor(mdle, xnod)
      call find_elem_nodes(mdle, norder,nedge_orient,nface_orient)
!
!  ...select appropriate visualization object
      ntype = NODES(mdle)%ntype
      nV = n_elem_vert(iel)
!
      do iv=1,nV
!
!     ...visualization point accounting for VLEVEL
         call get_vis_point(vis_on_type(ntype),iv-1, xi)
!
!     ...H1 shape functions at visualization point
         call shape3DH(ntype,xi,norder,nedge_orient,nface_orient, nrdofH,shapH,gradH)
!
!     ...geometry map
         select case(PARAVIEW_GEOM)
!
!     -- ISOPARAMETRIC GEOMETRY MAP --
         case(PARAVIEW_ISOGEOM)
            x(1:3)=0.d0
            do k=1,nrdofH
               x(1:3) = x(1:3) + xnod(1:3,k)*shapH(k)
            enddo
!
!     -- EXACT GEOMETRY MAP --
         case(PARAVIEW_EXGEOM)
            call exact_geom(mdle,xi, x,dxdxi)
         end select
!
!     ...add point to visualization
         GEOM_PTS(1:3,n_vert_offset(iel)+iv) = x
      enddo
   enddo
!$OMP END DO
!
!..Step 3 : Elements
!
!..Computing the total number of elements (incl. subelements if VLEVEL > 0)
!..VTU Format needs this information a-priori as VTU needs headers
!  with this information before appending the data.
   if (VIS_VTU) then
   !$OMP DO                            &
   !$OMP PRIVATE(mdle,ndom,ntype,ivis) &
   !$OMP REDUCTION(+:ice_subd)
      do iel = 1,NRELES
         mdle = ELEM_ORDER(iel)
         if (PARAVIEW_DOMAIN.ne.0) then
            call find_domain(mdle, ndom)
            if (ndom.ne.PARAVIEW_DOMAIN) cycle
         endif
         if (SECOND_ORDER_VIS) then
            ice_subd = ice_subd + 1
         else
            ntype = NODES(mdle)%ntype
            call get_vis_nrelem(ntype, ivis)
            ice_subd = ice_subd + ivis
         endif
      enddo
   !$OMP END DO
   endif
!$OMP END PARALLEL
!
   if (VIS_VTU) then
      allocate(VTU_ELEM_TYPES(ice_subd))
      VTU_ELEM_TYPES = 0
   endif
!
   ice_subd=0; icn_subd=0
!
!$OMP PARALLEL
!$OMP DO                                     &
!$OMP PRIVATE(mdle,ntype,i,j,k,l,ivis,ndom,  &
!$OMP         nodesl,norientl,nverl,subd)    &
!$OMP REDUCTION(+:ice_subd,icn_subd)         &
!$OMP SCHEDULE(DYNAMIC)
   do iel=1,NRELES
      mdle = ELEM_ORDER(iel)
!
!  ...SKIP IF WRONG DOMAIN
      if (PARAVIEW_DOMAIN.ne.0) then
         call find_domain(mdle, ndom)
         if (ndom.ne.PARAVIEW_DOMAIN) cycle
      endif
!
      if (DISTRIBUTED) then
         call get_subd(mdle, subd)
         if (RANK .ne. subd) cycle
      endif
!
      call elem_nodes(mdle, nodesl,norientl)
!
!  ...select appropriate visualization object
      ntype = NODES(mdle)%ntype
!
      j = nobj_conf(ntype)
      l = ivis_type(ntype)
      k = n_obj_offset(iel)
!
      call get_vis_nrelem(ntype, ivis)
!
      do i=1,ivis
!
         if (VIS_VTU) VTU_ELEM_TYPES(VTU_element_type_offset(iel) + i) = l
!
!     ...visualization element accounting for VLEVEL
         call get_vis_elem(vis_on_type(ntype),i,n_vert_offset(iel), nverl)
!
         if (SECOND_ORDER_VIS) then
            GEOM_OBJ(k+1:k+j) = nverl(1:j)
         else
            GEOM_OBJ(k+1) = l
            GEOM_OBJ(k+2:k+1+j) = nverl(1:j)
         endif
!
         k = k + (1+j)
         ice_subd = ice_subd + 1
         icn_subd = icn_subd + j
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL
!
!..Step 4 : Collect on host
!
   if (DISTRIBUTED .and. .not. HOST_MESH) then
      count = 1; IcE=0; IcN=0;
      call MPI_REDUCE(ice_subd,IcE,count,MPI_INTEGER,MPI_SUM,ROOT,MPI_COMM_WORLD, ierr)
      call MPI_REDUCE(icn_subd,IcN,count,MPI_INTEGER,MPI_SUM,ROOT,MPI_COMM_WORLD, ierr)
      if (RANK .eq. ROOT) then
         count = 3*IcP
         call MPI_REDUCE(MPI_IN_PLACE,GEOM_PTS,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD, ierr)
         count = ico
         call MPI_REDUCE(MPI_IN_PLACE,GEOM_OBJ,count,MPI_INTEGER,MPI_SUM,ROOT,MPI_COMM_WORLD, ierr)
         if (VIS_VTU) then
            count = size(VTU_ELEM_TYPES)
            call MPI_REDUCE(MPI_IN_PLACE,VTU_ELEM_TYPES,count,MPI_INTEGER,MPI_SUM,ROOT,MPI_COMM_WORLD, ierr)
         endif
      else
         count = 3*IcP
         call MPI_REDUCE(GEOM_PTS,GEOM_PTS,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD, ierr)
         count = ico
         call MPI_REDUCE(GEOM_OBJ,GEOM_OBJ,count,MPI_INTEGER,MPI_SUM,ROOT,MPI_COMM_WORLD, ierr)
         if (VIS_VTU) then
            count = size(VTU_ELEM_TYPES)
            call MPI_REDUCE(VTU_ELEM_TYPES,VTU_ELEM_TYPES,count,MPI_INTEGER,MPI_SUM,ROOT,MPI_COMM_WORLD, ierr)
         endif
      endif
   else
      IcE = ice_subd
      IcN = icn_subd
   endif
!
!..Step 5 : Write to file
!
   if (RANK .eq. ROOT) then
      if (VIS_VTU) then
!     ...with standard I/O
         call write_VTU_geom(IcE)
      else
!     ...with HDF5
         call geometry_write(Sname,len(Sname),Sfile,len(Sfile))
      endif
   endif
!
!..Step 6 : Deallocate
!
   if (VIS_VTU) then
      deallocate(VTU_ELEM_TYPES)
      deallocate(VTU_element_type_offset)
   endif
!
   call geometry_close()
!
end subroutine geom2vtk
