!----------------------------------------------------------------------------------------
!> @brief writes the header of the VTU file and also appends the coordinates,
!!        connectivity and element type data in VTU file.
!!
!> @param[in] IcE - total number of elements for visualization
!!
!> @date Sep 2023
!----------------------------------------------------------------------------------------
subroutine write_VTU_headers(IcE)
!
   use data_structure3D
   use upscale
   use paraview
   use physics
!
   implicit none
!
   integer,intent(in)      :: IcE
!
   integer, allocatable    :: elem_connectivity(:,:)
   integer, allocatable    :: offsets_connectivity(:)
!
!..keeps track of data offset and size when using VTU format
!  (VIS_FORMAT = 1)
   integer                 :: VTU_data_offset
   integer                 :: VTU_data_size
!..Auxiliary variables
   character(len=16)       :: str1,str2,str3,str4
   integer                 :: k,l,iv,j,count,nV
   integer                 :: iphys,iload,icomp,jcomp
!
   integer                 :: ipart,npart
   character(len=4)        :: suffix
!
!..Export both real- and imaginary-part for complex solutions
#if C_MODE
   npart = 2
#else
   npart = 1
#endif
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
   write(str1, '(i1)')   0                         ! data_offset
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
   elem_connectivity = 0
   offsets_connectivity = 0
!
   do count  = 1,IcE
      if (SECOND_ORDER_VIS) then
         l = VTU_ELEM_TYPES(count)
         j = nobj_conf_VTU(l)
         do iv = 1,j
            elem_connectivity(count,iv) = GEOM_OBJ(k+iv)
         enddo
         k = k + j
         offsets_connectivity(count) = k
      else
         k = k + 1
         l = VTU_ELEM_TYPES(count)
         j = nobj_conf_VTU(l)
         do iv = 1,j
            elem_connectivity(count,iv) = GEOM_OBJ(k+iv)
         enddo
         k = k + j
         if (count .eq. 1) then
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
!..loop over solution copies
   do iload=1,NRCOMS
!
      if (.not. PARAVIEW_LOAD(iload)) cycle
!
      jcomp = 0
!
!  ...loop over physics variables
      do iphys=1,NR_PHYSA
!
         if (.not. PARAVIEW_ATTR(iphys)) then
            jcomp = jcomp + NR_COMP(iphys)
            cycle
         endif
!
!     ...loop over components
         do icomp=1,NR_COMP(iphys)
            jcomp = jcomp+1
!
!        ...loop over real/imag parts
            do ipart = 1,npart
               if (ipart .eq. 1) then
                  if (.not. PARAVIEW_COMP_REAL(jcomp)) cycle
                  suffix = "real"
               else
                  if (.not. PARAVIEW_COMP_IMAG(jcomp)) cycle
                  suffix = "imag"
               endif
!
!           ...specify type of physics variable
               select case(D_TYPE(iphys))
!
!           ...H1-variable
               case(CONTIN)
                  write(str1, '(i0.0)') VTU_data_offset        ! data_offset
                  write(str2, '(i0.0)') 8 * 8                  ! integer precision
                  write(str3, '(i2.2)') icomp
                  write(str4, '(i2.2)') iload
                  write(PARAVIEW_IO) '        '//'<DataArray type="Float'//trim(str2)//'"'//    &
                     ' Name="'//PHYSA(iphys)//'_comp_'//trim(str3)//                             &
                                             '_load_'//trim(str4)//'_'//suffix//'"'//           &
                     ' format="appended"'//                                                     &
                     ' offset="'//trim(str1)//'">'//char(10)
                  write(PARAVIEW_IO) '        '//'</DataArray>'//char(10)
                  VTU_data_offset = VTU_data_offset + 4 + nV * 8
!
!           ...H(curl)-variable
               case(TANGEN)
                  write(str1, '(i0.0)') VTU_data_offset        ! data_offset
                  write(str2, '(i0.0)') 8 * 8                  ! integer precision
                  write(str3, '(i2.2)') icomp
                  write(str4, '(i2.2)') iload
                  write(PARAVIEW_IO) '        '// '<DataArray type="Float'//trim(str2)//'"'//   &
                     ' Name="'//PHYSA(iphys)//'_comp_'//trim(str3)//                             &
                                             '_load_'//trim(str4)//'_'//suffix//'"'//           &
                     ' NumberOfComponents="3"'//                                                &
                     ' format="appended"'//                                                     &
                     ' offset="'//trim(str1)//'">'//char(10)
                  write(PARAVIEW_IO) '        '//'</DataArray>'//char(10)
                  VTU_data_offset = VTU_data_offset + 4 + nV * 3 * 8
!
!           ...H(div)-variable
               case(NORMAL)
                  write(str1, '(i0.0)') VTU_data_offset        ! data_offset
                  write(str2, '(i0.0)') 8 * 8                  ! integer precision
                  write(str3, '(i2.2)') icomp
                  write(str4, '(i2.2)') iload
                  write(PARAVIEW_IO) '        '//'<DataArray type="Float'//trim(str2)//'"'//    &
                     ' Name="'//PHYSA(iphys)//'_comp_'//trim(str3)//                             &
                                             '_load_'//trim(str4)//'_'//suffix//'"'//           &
                     ' NumberOfComponents="3"'//                                                &
                     ' format="appended"'//                                                     &
                     ' offset="'//trim(str1)//'">'//char(10)
                  write(PARAVIEW_IO) '        '//'</DataArray>'//char(10)
                  VTU_data_offset = VTU_data_offset + 4 + nV * 3 * 8
!
!           ...L2-variable
               case(DISCON)
                  write(str1, '(i0.0)') VTU_data_offset        ! data_offset
                  write(str2, '(i0.0)') 8 * 8                  ! integer precision
                  write(str3, '(i2.2)') icomp
                  write(str4, '(i2.2)') iload
                  write(PARAVIEW_IO) '        '//'<DataArray type="Float'//trim(str2)//'"'//    &
                     ' Name="'//PHYSA(iphys)//'_comp_'//trim(str3)//                             &
                                             '_load_'//trim(str4)//'_'//suffix//'"'//           &
                     ' format="appended"'//                                                     &
                     ' offset="'//trim(str1)//'">'//char(10)
                  write(PARAVIEW_IO) '        '//'</DataArray>'//char(10)
                  VTU_data_offset = VTU_data_offset + 4 + nV * 8
!
               end select
!        ...end loop over real/imag parts
            enddo
!     ...end loop over components
         enddo
!  ...end loop over physics variables
      enddo
!..end loop over solution copies
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
   write(PARAVIEW_IO) VTU_data_size ! data_size = 3 coordinates for each node * 8
   do count = 1, nV
      write(PARAVIEW_IO) GEOM_PTS(1,count),GEOM_PTS(2,count),GEOM_PTS(3,count)
   enddo
!
!..writing connectivity data
   VTU_data_size = offsets_connectivity(IcE) * 4
   write(PARAVIEW_IO) VTU_data_size
   do count  = 1,IcE
      l = VTU_ELEM_TYPES(count)
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
      write(PARAVIEW_IO) VTU_ELEM_TYPES(count)
   enddo
!
!..deallocation
   deallocate(offsets_connectivity)
   deallocate(elem_connectivity)
!
end subroutine write_VTU_headers
