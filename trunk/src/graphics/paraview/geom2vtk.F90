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
!> @date Mar 2023
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

   if(VIS_VTU) then 
      allocate(VTU_element_type_offset(NRELES))
   endif

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

      if(VIS_VTU) VTU_element_type_offset(iel) = m
      
      IcP = IcP + vis_obj%NR_VERT
!
!  ...second-order vis doesn't yet work with mixed meshes (paraview error)
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
!..Additional Computational for computing the total number of elements
!  (including subelements: vlevel > 0)
!..VTU Format needs this information a-priori as VTU needs headers
!  with this information before appending the data.
   if(VIS_VTU) then
      if(SECOND_ORDER_VIS) then
         allocate(ELEM_TYPES(NRELES))
         ELEM_TYPES = ZERO
      else
         do iel=1,NRELES
            mdle = ELEM_ORDER(iel)

            if (PARAVIEW_DOMAIN.ne.0) then
               call find_domain(mdle, ndom)
               if (ndom.ne.PARAVIEW_DOMAIN) cycle
            endif

            ntype = NODES(mdle)%ntype
            call get_vis_nrelem(ntype, ivis)
            ice_subd = ice_subd + ivis
         enddo
         allocate(ELEM_TYPES(ice_subd))
         ELEM_TYPES = ZERO
      endif
      ice_subd = 0
   endif
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
!
      call get_vis_nrelem(ntype, ivis)
!
      do i=1,ivis
!
!     ...visualization element accounting for VLEVEL
         call get_vis_elem(vis_on_type(ntype),i,n_vert_offset(iel), nverl)
!
         if (SECOND_ORDER_VIS) then
            GEOM_OBJ(k+1:k+j) = nverl(1:j)
            if(VIS_VTU) ELEM_TYPES(VTU_element_type_offset(iel) + i) = l
         else
            GEOM_OBJ(k+1) = l
            GEOM_OBJ(k+2:k+1+j) = nverl(1:j)
            if(VIS_VTU) ELEM_TYPES(VTU_element_type_offset(iel) + i) = l
         endif
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
      call MPI_REDUCE(ice_subd,IcE,count,MPI_INTEGER,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      call MPI_REDUCE(icn_subd,IcN,count,MPI_INTEGER,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      if (RANK .eq. ROOT) then
         count = 3*IcP
         call MPI_REDUCE(MPI_IN_PLACE,GEOM_PTS,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         count = ico
         call MPI_REDUCE(MPI_IN_PLACE,GEOM_OBJ,count,MPI_INTEGER,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)

         if(VIS_VTU) then
            count = size(ELEM_TYPES)
            call MPI_REDUCE(MPI_IN_PLACE,ELEM_TYPES,count,MPI_INTEGER,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         endif
      else
         count = 3*IcP
         call MPI_REDUCE(GEOM_PTS,GEOM_PTS,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         count = ico
         call MPI_REDUCE(GEOM_OBJ,GEOM_OBJ,count,MPI_INTEGER,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
  
         if(VIS_VTU) then
            count = size(ELEM_TYPES)
            call MPI_REDUCE(ELEM_TYPES,ELEM_TYPES,count,MPI_INTEGER,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
         endif
      endif
   else
      IcE = ice_subd
      IcN = icn_subd
   endif
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!
!..Step 5 : Write to file with HDF5
!
   if (RANK .eq. ROOT) then
      if (.not. VIS_VTU) then
         call geometry_write(Sname,len(Sname),Sfile,len(Sfile))
      else
         call write_VTU_headers(IcE)
      endif
   endif
!
   if(VIS_VTU) then
      deallocate(ELEM_TYPES)
      deallocate(VTU_element_type_offset)
   endif
!
!..Step 6 : Deallocate
!
   call geometry_close()
!
end subroutine geom2vtk



!----------------------------------------------------------------------------------------
!> @brief writes the header for VTU file and also appends the co-ordinate, connectivity and element type data in VTU file.
!!
!> @param[in ] IcE - total number of elements for visualization
!!
!> @date Mar 2023
!----------------------------------------------------------------------------------------
subroutine write_VTU_headers(IcE)
!
!
   use data_structure3D
   use upscale
   use paraview
   use physics
!
   implicit none
!
   integer,intent(in)      ::  IcE
!
   integer, allocatable    :: elem_connectivity(:,:)
   integer, allocatable    :: offsets_connectivity(:)
!
!..keeps track of data offset and size when using VTU format
!  (VIS_FORMAT = 1)
   integer                 :: VTU_data_offset
   integer                 :: VTU_data_size
!..Auxiliary variables
   character(len=80)       :: str1, str2, str3, str4
   integer                 :: k,l,iv,j,count, nV
   integer                 :: iphys,iload,icomp
!
!..Main header starts here
!
!..number of vertices or nodes
   nV = size(GEOM_PTS,dim=2)
!
   write(PARAVIEW_IO) ''// '<?xml version="1.0"?>' //char(10)
   write(PARAVIEW_IO) ''// '<VTKFile type="UnstructuredGrid" version="0.1" ' //  &
                           'byte_order="LittleEndian">'                      // char(10)
   write(PARAVIEW_IO) '  '// '<UnstructuredGrid>' // char(10)
   write(str1, '(i0.0)') nV
   write(str2, '(i0.0)') IcE
   write(PARAVIEW_IO) '    ' // '<Piece NumberOfPoints="' // trim(str1) //  &
                     '" NumberOfCells ="' // trim(str2) // '">' // char(10)
!..Main header ends here
!
!..Header for coordinates of vertices/nodes..!
   write(PARAVIEW_IO) '      ' // '<Points>' // char(10)
   write(str1, '(i1)')   0                          ! data_offset
   write(str2, '(i0.0)') 8 * 8                     ! real precision
   write(PARAVIEW_IO) '        ' // '<DataArray type="Float' // trim(str2) // '"' //  &
                     ' NumberOfComponents="3"'                     //  &
                     ' format="appended"'                          //  &
                     ' offset="' // trim(str1) // '">' // char(10)
   write(PARAVIEW_IO) '        ' // '</DataArray>' // char(10)
   write(PARAVIEW_IO) '      ' // '</Points>'    // char(10)
!
!..Header for connectivity data starts here
   VTU_data_offset = 0
   VTU_data_offset =  VTU_data_offset + 4 + nV * 3 * 8
   write(PARAVIEW_IO) '      ' // '<Cells>' // char(10)
!
!..connectivity
   write(str1, '(i0.0)') VTU_data_offset        ! data_offset
   write(str2, '(i0.0)') 4 * 8                  ! integer precision
   write(PARAVIEW_IO) '        ' // '<DataArray type="Int' // trim(str2) // '"' //  &
                     ' Name="connectivity"'                      //  &
                     ' format="appended"'                        //  &
                     ' offset="' // trim(str1) // '">' // char(10)
   write(PARAVIEW_IO) '        ' // '</DataArray>' // char(10)
!..Header for connectivity data ends here
!
!..computations done to compute offsets for each element connectivity
   k = 0
   l = 0
!
   allocate(offsets_connectivity(IcE))
   allocate(elem_connectivity(IcE,(MAXP+1)**3))
   elem_connectivity = ZERO
   offsets_connectivity = ZERO
!
   do count  = 1,IcE
      if(SECOND_ORDER_VIS) then
         l = ELEM_TYPES(count)
         j = nobj_conf_VTU(l)
         do iv = 1,j
            elem_connectivity(count,iv) = GEOM_OBJ(k+iv) 
         enddo
         k = k + j
         offsets_connectivity(count) = k
      else
         k = k + 1
         l = ELEM_TYPES(count)
         j = nobj_conf_VTU(l)

         do iv = 1,j
            elem_connectivity(count,iv) = GEOM_OBJ(k+iv) 
         enddo
         k = k + j
         if(count .eq. 1) then
         offsets_connectivity(count) = k - 1
         else
         offsets_connectivity(count) = k - count
         endif
      endif
   enddo
!
!..Header for connectivity data offset starts here
   VTU_data_offset = VTU_data_offset + 4 + offsets_connectivity(IcE) * 4
   write(str1, '(i0.0)') VTU_data_offset              ! data_offset
   write(str2, '(i0.0)') 4 * 8                        ! integer precision
   write(PARAVIEW_IO) '        ' // '<DataArray type="Int' // trim(str2) // '"' //  &
   ' Name="offsets"'                           //  &
   ' format="appended"'                        //  &
   ' offset="' // trim(str1) // '">' // char(10)
   write(PARAVIEW_IO) '        ' // '</DataArray>' // char(10)
!..Header for connectivity data offset ends here
!
!..Header for element types starts here
   VTU_data_offset = VTU_data_offset + 4 + IcE * 4 
   write(str1, '(i0.0)') VTU_data_offset              ! data_offset
   write(str2, '(i0.0)') 4 * 8                        ! integer precision
   write(PARAVIEW_IO) '        ' // '<DataArray type="Int' // trim(str2) // '"' //  &
   ' Name="types"'                             //  &
   ' format="appended"'                        //  &
   ' offset="' // trim(str1) // '">' // char(10)
   write(PARAVIEW_IO) '        ' // '</DataArray>' // char(10)
   write(PARAVIEW_IO) '      ' // '</Cells>' // char(10)
!..Header for element types ends here
!
!..Headers for attributes starts here
   write(PARAVIEW_IO) '      ' // '<PointData Scalars="scalars">' // char(10)
!..Headers for attributes ends here
!
   VTU_data_offset = VTU_data_offset + 4 + IcE * 4
!
   if (.not. PARAVIEW_DUMP_ATTR) goto 50
!
!..Header for solution data starts here
!
!..loop over NRCOMS
   do iload=1,NRCOMS
!  ...loop over physics variables
      do iphys=1,NR_PHYSA
         if (IPARATTR_VTU(iphys) .eq. 0) cycle
!     ...loop over components
         do icomp=1,NR_COMP(iphys)
!
            if (IPARATTR_VTU(iphys) .lt. icomp) cycle
!
            select case(D_TYPE(iphys))
!
               case(CONTIN)
                  write(str1, '(i0.0)') VTU_data_offset        ! data_offset
                  write(str2, '(i0.0)') 8 * 8                  ! integer precision
                  write(str3, '(i0.0)') iphys
                  write(str4, '(i0.0)') icomp
                  write(PARAVIEW_IO) '        ' // '<DataArray type="Float' // trim(str2) // '"' //  &
                  ' Name='//'"H1_'//trim(str3)//'_'//trim(str4)//'_var"'                     //  &
                  ' format="appended"'                          //  &
                  ' offset="' // trim(str1) // '">' // char(10)
                  write(PARAVIEW_IO) '        ' // '</DataArray>' // char(10)
!
                  VTU_data_offset = VTU_data_offset + 4 + nV * 8
!
               case(TANGEN)
                  write(str1, '(i0.0)') VTU_data_offset        ! data_offset
                  write(str2, '(i0.0)') 8 * 8                  ! integer precision
                  write(str3, '(i0.0)') iphys
                  write(PARAVIEW_IO) '        ' // '<DataArray type="Float' // trim(str2) // '"' //  &
                  ' Name='//'"HCurl_'//trim(str3)//'_var"'//' NumberOfComponents="3"' //  &
                  ' format="appended"'                          //  &
                  ' offset="' // trim(str1) // '">' // char(10)
                  write(PARAVIEW_IO) '        ' // '</DataArray>' // char(10)
!
                  VTU_data_offset = VTU_data_offset + 4 + nV * 3 * 8
!
               case(NORMAL)
                  write(str1, '(i0.0)') VTU_data_offset        ! data_offset
                  write(str2, '(i0.0)') 8 * 8                  ! integer precision
                  write(str3, '(i0.0)') iphys
                  write(PARAVIEW_IO) '        ' // '<DataArray type="Float' // trim(str2) // '"' //  &
                  ' Name='//'"HDiv_'//trim(str3)//'_var"'//' NumberOfComponents="3"' //  &
                  ' format="appended"'                          //  &
                  ' offset="' // trim(str1) // '">' // char(10)
                  write(PARAVIEW_IO) '        ' // '</DataArray>' // char(10)
!
                  VTU_data_offset = VTU_data_offset + 4 + nV * 3 * 8
!
               case(DISCON)
                  write(str1, '(i0.0)') VTU_data_offset        ! data_offset
                  write(str2, '(i0.0)') 8 * 8                  ! integer precision
                  write(str3, '(i0.0)') iphys
                  write(str4, '(i0.0)') icomp
                  write(PARAVIEW_IO) '        ' // '<DataArray type="Float' // trim(str2) // '"' //  &
                  ' Name='//'"L2_'//trim(str3)//'_'//trim(str4)//'_var"'                     //  &
                  ' format="appended"'                          //  &
                  ' offset="' // trim(str1) // '">' // char(10)
                  write(PARAVIEW_IO) '        ' // '</DataArray>' // char(10)
                  VTU_data_offset = VTU_data_offset + 4 + nV * 8
            end select
!     ...end loop over components
         enddo
!  ...end loop over physics variables
      enddo
!..end loop over NRCOMS
   enddo
!
   50 continue
!
!..Header for solution data ends here
!
!..closing data headers
   write(PARAVIEW_IO) '      ' // '</PointData>' // char(10)
   write(PARAVIEW_IO) '    ' // '</Piece>'             // char(10)
   write(PARAVIEW_IO) '  ' // '</UnstructuredGrid>'  // char(10)
!
!..appending data in binary
   write(PARAVIEW_IO) '' // '<AppendedData encoding="raw">' // char(10)
   write(PARAVIEW_IO) '_'
!
!..writing point coordinates
   VTU_data_size =  nV * 3 * 8
   write(PARAVIEW_IO) VTU_data_size ! data_size =  three coordinates for each node * 8
   do count = 1, nV
      write(PARAVIEW_IO) GEOM_PTS(1,count),GEOM_PTS(2,count),GEOM_PTS(3,count)
   enddo
!
!..writing connectivity data
   VTU_data_size = offsets_connectivity(IcE) * 4
   write(PARAVIEW_IO) VTU_data_size
   do count  = 1,IcE
      l = ELEM_TYPES(count)
      j = nobj_conf_VTU(l)
      do iv = 1,j
      write(PARAVIEW_IO) elem_connectivity(count,iv)
      enddo
   enddo
!
!..writing connectivity offsets
   VTU_data_size = IcE * 4
   write(PARAVIEW_IO) VTU_data_size
   do count = 1, IcE
      write(PARAVIEW_IO) offsets_connectivity(count)
   enddo
!
!..writing element types
   VTU_data_size = IcE * 4
   write(PARAVIEW_IO) VTU_data_size
   do count = 1, IcE
      write(PARAVIEW_IO) ELEM_TYPES(count)
   enddo
!
!..deallocation
   deallocate(offsets_connectivity)
   deallocate(elem_connectivity)
!
end subroutine write_VTU_headers
