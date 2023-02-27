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
!> @date Feb 2023
!----------------------------------------------------------------------------------------
subroutine geom2vtk(Sname,Sfile, IcE,IcN,IcP)
!
   use data_structure3D
   use upscale
   use paraview
   use MPI           , only: MPI_COMM_WORLD,MPI_SUM,MPI_INTEGER,MPI_Wtime
   use mpi_param     , only: RANK,ROOT,NUM_PROCS
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
   integer :: nodesl(27), norder(19), nverl(8)
   integer :: norientl(27), nedge_orient(12), nface_orient(6)
!
   real(8) :: xi(3), x(3), xnod(3,MAXbrickH), dxdxi(3,3)
!
!..OpenMP parallelization: auxiliary variables
   integer, dimension(NRELES) :: n_vert_offset,n_obj_offset,n_elem_vert
!
!..timer
   !real(8) :: start_time,end_time
!
!..MPI
   integer :: ierr,count,subd,ice_subd,icn_subd
!
!----------------------------------------------------------------------------------------
!
!   if (RANK .eq. ROOT) write(*,*) 'geom2vtk:'
!   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!
!..Step 1 : Preliminary calculations (offsets, etc.)
!
!..create list of mdle node numbers,
!  and count number of visualization points
   IcP=0; ico=0;
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
      IcP = IcP + vis_obj%NR_VERT
      ico = ico + vis_obj%NR_ELEM*(nobj_conf(ntype)+1)
   enddo
!
!   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
!   if (RANK .eq. ROOT) write(*,300) end_time - start_time
!  300 format(' timer: ',f12.5,' seconds')
!
   call geometry_init(ico,Icp)
   GEOM_PTS(1:3,1:Icp) = 0.d0; GEOM_OBJ(1:ico) = 0
!
!   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
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
   ice_subd=0; icn_subd=0
!
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
      call get_vis_nrelem(ntype, ivis)
      do i=1,ivis
!
!     ...visualization element accounting for VLEVEL
         call get_vis_elem(vis_on_type(ntype),i,n_vert_offset(iel), nverl)
!
         GEOM_OBJ(k+1) = l
         GEOM_OBJ(k+2:k+1+j) = nverl(1:j)
         k = k + (1+j)
         ice_subd = ice_subd + 1
         icn_subd = icn_subd + nvert(ntype)
      enddo
   enddo
!$OMP END DO
!$OMP END PARALLEL
!
!   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
!   if (RANK .eq. ROOT) write(*,300) end_time - start_time
!
!..Step 4 : Collect on host
!
   if (DISTRIBUTED .and. .not. HOST_MESH) then
      count = 1; IcE=0; IcN=0;
      call MPI_REDUCE(ice_subd,IcE,count,MPI_INTEGER,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(icn_subd,IcN,count,MPI_INTEGER,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      if (RANK .eq. ROOT) then
         count = 3*IcP
         call MPI_REDUCE(MPI_IN_PLACE,GEOM_PTS,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         count = ico
         call MPI_REDUCE(MPI_IN_PLACE,GEOM_OBJ,count,MPI_INTEGER,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      else
         count = 3*IcP
         call MPI_REDUCE(GEOM_PTS,GEOM_PTS,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         count = ico
         call MPI_REDUCE(GEOM_OBJ,GEOM_OBJ,count,MPI_INTEGER,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      endif
   else
      IcE = ice_subd
      IcN = icn_subd
   endif
!
!..Step 5 : Write to file with HDF5
!
!   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
   if (RANK .eq. ROOT) call geometry_write(Sname,len(Sname),Sfile,len(Sfile))
!   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
!   if (RANK .eq. ROOT) write(*,300) end_time - start_time
!
!..Step 6 : Deallocate
!
   call geometry_close()
!
end subroutine geom2vtk

