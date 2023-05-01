!----------------------------------------------------------------------------------------
!> Purpose : write geometry to .h5 file
!!
!> @param[in ] Sname - name/description (e.g., 'Geometry')
!> @param[in ] Sfile - geometry file name (e.g., ../outputs/paraview/geom_00000.h5)
!> @param[out] IcE   - Number of all vis object subelements in the active mesh
!> @param[out] IcN   - Number of all vis object subelement vertices in the active mesh
!> @param[out] IcP   - Number of points (coords) of all vis objects in the active mesh
!!
!> @date Oct 2019
!----------------------------------------------------------------------------------------
!
#include "typedefs.h"
!
subroutine geom2vtk(Sname,Sfile, IcE,IcN,IcP)
!
   use data_structure3D
   use upscale
   use paraview
   use MPI           , only: MPI_COMM_WORLD,MPI_SUM,MPI_INTEGER
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
   real(8),dimension(  MAXbrickH) :: shapH
   real(8),dimension(3,MAXbrickH) :: gradH
   character(len=4) :: etype
   integer :: iv, nV, ico, iel, mdle, nrdofH, i, j, k, l, ivis, ndom
   integer :: nodesl(27), norder(19), nverl(27)
   integer :: norientl(27), nedge_orient(12), nface_orient(6)
!
   real(8) :: xi(3), x(3), xnod(3,MAXbrickH), dxdxi(3,3)
!
!..OpenMP parallelization: auxiliary variables
   integer, dimension(NRELES) :: n_vert_offset,n_obj_offset,n_elem_vert
   integer, dimension((MAXP+1)**3) :: elem_connectivity
   integer, allocatable :: offsets_connectivity(:)
   integer, allocatable :: elem_types(:)
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
      etype = NODES(mdle)%type
      vis_obj = vis_on_type(etype)
      n_vert_offset(iel) = IcP
      n_obj_offset(iel) = ico
      n_elem_vert(iel) = vis_obj%NR_VERT
      IcP = IcP + vis_obj%NR_VERT
      ico = ico + vis_obj%NR_ELEM*(nobj_conf(etype)+1)
   enddo
   ! write(*,*) " I am here"
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
!$OMP PRIVATE(mdle,ndom,etype,iv,xi,nV,subd,    &
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
      etype = NODES(mdle)%type
      nV = n_elem_vert(iel)
!
      do iv=1,nV
!
!     ...visualization point accounting for VLEVEL
         call get_vis_point(vis_on_type(etype),iv-1, xi)
!
!     ...H1 shape functions at visualization point
         call shape3DH(etype,xi,norder,nedge_orient,nface_orient, nrdofH,shapH,gradH)
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
!$OMP PRIVATE(mdle,etype,i,j,k,l,ivis,ndom,  &
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
      etype = NODES(mdle)%type
!
      j = nobj_conf(etype)
      l = ivis_type(etype)
      k = n_obj_offset(iel)
      call get_vis_nrelem(etype, ivis)
      do i=1,ivis
!
!     ...visualization element accounting for VLEVEL
         call get_vis_elem(vis_on_type(etype),i,n_vert_offset(iel), nverl)
!
         GEOM_OBJ(k+1) = l
         GEOM_OBJ(k+2:k+1+j) = nverl(1:j)
         k = k + (1+j)
         ice_subd = ice_subd + 1
         icn_subd = icn_subd + nvert(etype)
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
   if (RANK .eq. ROOT) then
      if(VIS_FORMAT .eq. 0) then
         call geometry_write(Sname,len(Sname),Sfile,len(Sfile))
      else
         nV = size(GEOM_PTS,dim=2)
         write(PARAVIEW_IO,1032)
         write(PARAVIEW_IO,1033)
         write(PARAVIEW_IO,1034)nV,ice
         write(PARAVIEW_IO,1035)
         write(PARAVIEW_IO,1036)

         1032 format("<VTKFile type=","""UnstructuredGrid""",">")
         1033 format("<UnstructuredGrid>")
         1034 format("<Piece NumberOfPoints=""",i0,""""," NumberOfCells=""",i0,""">")
         1035 format("<Points>")
         1036 format("<DataArray type=""","Float32"" Name=""Points"" NumberOfComponents=""3"" format=""ascii"">")

         !writting ascii format data for the coordinates
         
         do count = 1,nV
            write(PARAVIEW_IO,*) GEOM_PTS(1,count),GEOM_PTS(2,count),GEOM_PTS(3,count)
            ! Encode(PARAVIEW_IO,*) GEOM_PTS(1,count),GEOM_PTS(2,count),GEOM_PTS(3,count)

         enddo

         write(PARAVIEW_IO,1037)
         write(PARAVIEW_IO,1038)

         1037 format("</DataArray>")
         1038 format("</Points>")

         write(PARAVIEW_IO,1039)
         write(PARAVIEW_IO,1040)

         1039 format("<Cells>")
         1040 format("<DataArray type=""","Int64"" IdType=""1""  Name=""connectivity"" format=""ascii"">")

         k = 0
         l = 0
         elem_connectivity = ZERO
         allocate(offsets_connectivity(IcE))
         allocate(elem_types(IcE))
         offsets_connectivity = ZERO
         elem_types = ZERO
         
         do count  = 1,IcE

            k = k + 1
            l = GEOM_OBJ(k)
            elem_types(count) = l
            j = nobj_conf_reverse(l)

            do iv = 1,j
               elem_connectivity(iv) = GEOM_OBJ(k+iv) 
            enddo
            write(PARAVIEW_IO,*) elem_connectivity(1:j)

            k = k + j
            if(count .eq. 1) then
            offsets_connectivity(count) = k - 1
            else
            offsets_connectivity(count) = k - count
            endif
         enddo

         write(PARAVIEW_IO,1041)
         write(PARAVIEW_IO,1042)
         

         1041 format("</DataArray>")
         1042 format("<DataArray type=""","Int64"" IdType=""1""  Name=""offsets"" format=""ascii"">")
         
         write(PARAVIEW_IO,*) offsets_connectivity

         write(PARAVIEW_IO,1043)
         write(PARAVIEW_IO,1044)

         1043 format("</DataArray>")
         1044 format("<DataArray type=""","UInt8"" Name=""types"" format=""ascii"">")

         write(PARAVIEW_IO,*) elem_types

         write(PARAVIEW_IO,1045)
         write(PARAVIEW_IO,1046)
         ! write(PARAVIEW_IO,1047)
         ! write(PARAVIEW_IO,1048)
         ! write(PARAVIEW_IO,1049)

         1045 format("</DataArray>")
         1046 format("</Cells>")
         write(PARAVIEW_IO,*) "<PointData Scalars=""","scalars"">"

         deallocate(offsets_connectivity)
         deallocate(elem_types)
         ! close(PARAVIEW_IO)
      endif
   endif
!   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
!   if (RANK .eq. ROOT) write(*,300) end_time - start_time
!
!..Step 6 : Deallocate
!
   call geometry_close()
!
end subroutine geom2vtk

