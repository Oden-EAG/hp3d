#include "typedefs.h"
!
!-------------------------------------------------------------------------------------------
!> Purpose : paraview driver
!-------------------------------------------------------------------------------------------
!
subroutine my_paraview_driver(IParAttr)
!
   use upscale
   use physics,            only: DTYPE,NR_COMP,NR_PHYSA,PHYSA
   use data_structure3D,   only: NRCOMS
   use environment,        only: PREFIX,QUIET_MODE
   use paraview,           only: PARAVIEW_DUMP_ATTR,FILE_VIS, &
                                 PARAVIEW_DOMAIN,VLEVEL, &
                                 PARAVIEW_DUMP_GEOM
   use mpi_param,          only: RANK,ROOT
   use par_mesh,           only: DISTRIBUTED,HOST_MESH
   use mg_data_structure,  only: NRELES_GRID
!
   implicit none
!
   integer, intent(in) :: IParAttr(NR_PHYSA)
!
   real(8) :: time
   integer :: idx,iphys,iload,icomp
!
   integer, save :: id = -1
   logical, save :: initialized = .false.
!
   character(len=2) :: vis_level
!
!-------------------------------------------------------------------------------------------
!
!..Set true to write out geometry file on every call to this routine
!  Set false to write out geometry file only on first call to this routine
   PARAVIEW_DUMP_GEOM = .true.
!
!..load files for visualization upscale
   if (SECOND_ORDER_VIS) then
      if (.not. initialized) then
         call load_vis(TETR_VIS,trim(FILE_VIS)//'/tetra10',TETR)
         call load_vis(PRIS_VIS,trim(FILE_VIS)//'/prism18',PRIS)
         call load_vis(HEXA_VIS,trim(FILE_VIS)//'/hexa27' ,BRIC)
         initialized = .true.
      endif
   else
      if (.not. initialized) then
         call load_vis(TETR_VIS,trim(FILE_VIS)//'/tetra_'//trim(VLEVEL),TETR)
         call load_vis(PRIS_VIS,trim(FILE_VIS)//'/prism_'//trim(VLEVEL),PRIS)
         call load_vis(HEXA_VIS,trim(FILE_VIS)//'/hexa_'//trim(VLEVEL),BRIC)
         initialized = .true.
      endif
   endif
!
   time=-1.d0
!
!..integer id to append to Fname
   id=id+1
!
!  WRITE GEOMETRY WITH VIS 0 TO FILE (VISUALIZE ADAPTIVE REFS)
   if (RANK .eq. ROOT) then
      call paraview_begin_vis0(id,time)
   endif
!
   if (SECOND_ORDER_VIS) then
!  ...simply call normal routine
      call paraview_geom
   else
      vis_level = VLEVEL; VLEVEL = '0'
      call load_vis(TETR_VIS,trim(FILE_VIS)//'/tetra_'//trim(VLEVEL),TETR)
      call load_vis(PRIS_VIS,trim(FILE_VIS)//'/prism_'//trim(VLEVEL),PRIS)
      call load_vis(HEXA_VIS,trim(FILE_VIS)//'/hexa_'//trim(VLEVEL),BRIC)
!
      call paraview_geom_vis0
   endif
!
!..If multigrid solver active
   if (allocated(NRELES_GRID)) then
      if (DISTRIBUTED .and. (.not. HOST_MESH)) then
!     ...write subdomains
         call paraview_attr_subd(id, 1)
!     ...write number of elements
         call paraview_attr_subd(id, 2)
!     ...write polynomial order
         call paraview_attr_subd(id, 3)
!     ...write padded subdomain overlap
         call paraview_attr_subd(id, 4)
      endif
   endif
!
   if (SECOND_ORDER_VIS) then
!  ...do nothing; second order geometry still initialized
   else
      VLEVEL = vis_level
      call load_vis(TETR_VIS,trim(FILE_VIS)//'/tetra_'//trim(VLEVEL),TETR)
      call load_vis(PRIS_VIS,trim(FILE_VIS)//'/prism_'//trim(VLEVEL),PRIS)
      call load_vis(HEXA_VIS,trim(FILE_VIS)//'/hexa_'//trim(VLEVEL),BRIC)
   endif
   if (RANK .eq. ROOT) then
      call paraview_end_vis0
   endif
!..END ADAPTIVE REFS GEOMETRY WRITING
!
   if (RANK .eq. ROOT) then
      call paraview_begin(id,time) ! [OPENS THE XMF FILE, WRITES HEADER]
   endif
!
!  -- GEOMETRY --
   call paraview_geom
!
!  -- PHYSICAL ATTRIBUTES --
   if (.not. PARAVIEW_DUMP_ATTR) goto 90
!
!..loop over rhs's
   do iload=1,1 !NRCOMS
!
!  ...loop over physics variables
      do iphys=1,NR_PHYSA

         if (IParAttr(iphys) .eq. 0) cycle
!
!     ...loop over components
         do icomp=1,NR_COMP(iphys)
            if (IParAttr(iphys) .ge. icomp) then
!
!           ...encode iload, iphys, icomp into a single attribute's index
               idx = iload*100 + iphys*10 + icomp*1
!
               select case(DTYPE(iphys))
!
!                 -- H1 --
                  case('contin')
                     call paraview_attr_scalar(id,idx)
!
!                 -- H(curl) --
                  case('tangen')
                     call paraview_attr_vector(id,idx)
!
!                 -- H(div) --
                  case('normal')
                     call paraview_attr_vector(id,idx)
!
!                 -- L2 --
                  case('discon')
                     call paraview_attr_scalar(id,idx)
!
               end select
!
               if (.not. QUIET_MODE .and. RANK .eq. ROOT) then
                  write(*,8000) PHYSA(iphys),DTYPE(iphys),icomp,iload
 8000             format(/,1x,a5,7x,a6,12x,i1,5x,i2)
               endif
!
            endif
!     ...end loop over components of physics variable
         enddo
!  ...end loop over physics variables
      enddo
!..end loop over rhs's
   enddo
!
  90 continue
!
   if (RANK .eq. ROOT) then
      call paraview_end ! [CLOSES THE XMF FILE, WRITES FOOTER]
   endif
!
end subroutine my_paraview_driver

subroutine paraview_geom_vis0
!
   use environment , only : PREFIX
   use paraview    , only : PARAVIEW_IO,PARAVIEW_DUMP_GEOM,PARAVIEW_DIR
   use mpi_param   , only : RANK,ROOT
!
   implicit none
   integer,         save :: id=-1
   character(len=5),save :: postfix
   integer,         save :: ice, icn, icp
!
!-------------------------------------------------------------------------------------------
!
!...h5 file is produced on 1st visit, OR as required by the user
   if (PARAVIEW_DUMP_GEOM .or. (id == -1)) then
!
!  ...increment visitation flag
      id=id+1
!
!  ...convert integer to string
      write(postfix,"(I5.5)") id
!
!  ...produce .h5 file
      call geom2vtk("Geometry",trim(PARAVIEW_DIR)//trim(PREFIX)//"geomV0_"//trim(postfix)//".h5", ice,icn,icp)
!
   endif
!
   if (RANK .ne. ROOT) goto 90
!
!..write to .xmf file
   write(PARAVIEW_IO,1012) ice
   write(PARAVIEW_IO,1013) (ice+icn)
   write(PARAVIEW_IO,1014) trim(PREFIX), trim(postfix)
   write(PARAVIEW_IO,1015)
   write(PARAVIEW_IO,1016)
   write(PARAVIEW_IO,1017)
   write(PARAVIEW_IO,1018) icp
   write(PARAVIEW_IO,1019) trim(PREFIX), trim(postfix)
   write(PARAVIEW_IO,1020)
   write(PARAVIEW_IO,1021)
!
 1012 format("      <Topology TopologyType='Mixed' NumberOfElements='",i8,"'>")
 1013 format("        <DataItem Dimensions='",i12,"' NumberType='Int' Precision='4' Format='HDF'>")
 1014 format("        ",a,"geomV0_",a,".h5:/Objects")
 1015 format("        </DataItem>")
 1016 format("      </Topology>")
 1017 format("      <Geometry GeometryType='XYZ'>")
 1018 format("        <DataItem Dimensions='",i10, " 3' NumberType='Float' Precision='4' Format='HDF'>")
 1019 format("        ",a,"geomV0_",a,".h5:/Coords")
 1020 format("        </DataItem>")
 1021 format("      </Geometry>")
!
   90 continue
!
end subroutine paraview_geom_vis0


subroutine paraview_begin_vis0(Id,Time)
!
   use paraview    , only : PARAVIEW_IO,PARAVIEW_DIR,paraview_init
   use environment , only : PREFIX
!
   implicit none
!
   integer, intent(in) :: Id
   real(8), intent(in) :: Time
!
   character(len=5) :: postfix
!
!-------------------------------------------------------------------------------------------
!
   call paraview_init
!
!..convert integer to string
   write(postfix,"(I5.5)") Id
!
!..open .xmf file
   open(unit=PARAVIEW_IO , file=trim(PARAVIEW_DIR)//trim(PREFIX)//"V0_"//trim(postfix)//'.xmf')
!
   write(PARAVIEW_IO, 1001)
   write(PARAVIEW_IO, 1002)
   write(PARAVIEW_IO, 1003)
!
!..non negative time is provided
   if (Time >= 0.d0) then
      write(PARAVIEW_IO, 1004) Time
   end if
!
!..HEADER of .xmf file
 1001 format("<Xdmf xmlns:xi='http://www.w3.org/2003/XInclude' Version='2.1'>")
 1002 format("  <Domain>")
 1003 format("    <Grid Name='Geometry' GridType='Uniform'>")
 1004 format("    <Time Value='", f14.10, "' />")
!
end subroutine paraview_begin_vis0
!
!
!-------------------------------------------------------------------------------------------
!> Purpose : print footer, close file
!!
!> @date Oct 2019
!-------------------------------------------------------------------------------------------
subroutine paraview_end_vis0
!
   use paraview , only : PARAVIEW_IO, paraview_finalize
!
   implicit none
!
!..FOOTER of .xmf file
   write(PARAVIEW_IO, 1004)
   write(PARAVIEW_IO, 1005)
   write(PARAVIEW_IO, 1006)
 1004 format("    </Grid>")
 1005 format("  </Domain>")
 1006 format("</Xdmf>")
!
   close(PARAVIEW_IO)
!
   call paraview_finalize
!
end subroutine paraview_end_vis0





!-------------------------------------------------------------------------------------------
!> Purpose : print subdomain info
!!
!> @date Oct 2019
!-------------------------------------------------------------------------------------------
subroutine paraview_attr_subd(Id, Idx)
!
   use environment , only : PREFIX
   use paraview    , only : PARAVIEW_IO,PARAVIEW_DIR
   use par_mesh    , only : DISTRIBUTED,HOST_MESH
   use mpi_param   , only : RANK,ROOT
!
   implicit none

   integer, intent(in) :: Id
   integer, intent(in) :: Idx
!
   character(len=60) :: fname,nick
   integer           :: ic,iload,iattr,icomp
   character(len=5)  :: postfix
   character(len=2)  :: comp,load
!
!-------------------------------------------------------------------------------------------
!
!
   if ((.not. DISTRIBUTED) .or. HOST_MESH) then
      write(*,*) 'This routines outputs subdomain information (distributed meshes only)'
      stop 1
   endif
!
!..convert integer to string
   write(postfix,"(I5.5)") Id
!
!..file name and nickname for attribute
   select case(idx)
   case(1)
      fname = "subdomains_"
      nick  = "subdomains"
   case(2)
      fname = "nreles_"
      nick  = "nreles"
   case(3)
      fname = "order_"
      nick  = "order"
   case(4)
      fname = "overlap_"
      nick  = "overlap"
   case default
      write(*,*) 'paraview_attr_subd: bad idx - ', idx
      stop 1
   end select
!
!..write to .h5 file
   call scalar2vtk_mg("Scalar",  &
      trim(PARAVIEW_DIR)//trim(PREFIX)//"scalar_"//trim(fname)//postfix//".h5", &
      trim(nick),Idx, ic)
!
   if (RANK .eq. ROOT) then
!
!  ...write to .xmf file
      write(PARAVIEW_IO, 1101) trim(nick)
      write(PARAVIEW_IO, 1102) ic
      write(PARAVIEW_IO, 1103) trim(PREFIX), trim(fname), trim(postfix), trim(nick)
      write(PARAVIEW_IO, 1104)
      write(PARAVIEW_IO, 1105)
!
 1101 format("      <Attribute Name='",a,"' AttributeType='Scalar' Center='Node'>")
 1102 format("        <DataItem Dimensions='",i10, " 1' NumberType='Float' Precision='4' Format='HDF'>")
 1103 format("        ",a,"scalar_",a,a,".h5:",a)
 1104 format("        </DataItem>")
 1105 format("      </Attribute>")
!
   endif
!
end subroutine paraview_attr_subd






!-------------------------------------------------------------------------------------------
!> Purpose : output extra attributes for multigrid solver (subd, nreles_subd, order, overlap)
!!
!> @date Feb 2023
!-------------------------------------------------------------------------------------------
subroutine scalar2vtk_mg(Sname,Sfile,Snick,Idx, Ic)
!
   use data_structure3D
   use element_data
   use physics          , only: DTYPE,ADRES
   use upscale
   use paraview
   use par_mesh         , only: DISTRIBUTED,HOST_MESH
   use MPI              , only: MPI_COMM_WORLD,MPI_SUM,MPI_INTEGER
   use mpi_param        , only: RANK,ROOT,NUM_PROCS
   use mg_data_structure, only: GRID, NODES_MG, CURRENT_GRID
   use index_map
!
   implicit none
!
   character(len=*), intent(in ) :: Sname
   character(len=*), intent(in ) :: Sfile
   character(len=*), intent(in ) :: Snick
   integer,          intent(in ) :: Idx
   integer,          intent(out) :: Ic
!
   type(vis) :: vis_obj
   integer   :: ntype
   integer   :: mdle, mdleC, iel, iv, nV
!
   real(8) :: val
   integer :: ibeg,iattr,icomp,isol,iload,ireal,ndom
   integer :: nord1,nord2,nord3,nord12,norder(19)
!
   real(8), external :: dreal_part,dimag_part
!
!..OpenMP parallelization: auxiliary variables
   integer, dimension(NRELES) :: n_vert_offset, n_elem_vert
!
   type(indmap) :: mdle2iel
!
!..MPI
   integer :: ierr,count,subd,Igrid
!
!----------------------------------------------------------------------------------------
!
!------------------------------------------------------------------
!..Step 1 : Preliminary calculations (offsets, etc.)
!------------------------------------------------------------------
!
   Ic=0
!
!..create list of mdle node numbers,
!  and count number of visualization points
   do iel=1,NRELES
      mdle = ELEM_ORDER(iel)
      if (PARAVIEW_DOMAIN.ne.0) then
         call find_domain(mdle, ndom)
         if (ndom.ne.PARAVIEW_DOMAIN) cycle
      endif
      ntype = NODES(mdle)%ntype
      vis_obj = vis_on_type(ntype)
      n_vert_offset(iel) = Ic
      n_elem_vert(iel) = vis_obj%nr_vert
      Ic = Ic + vis_obj%nr_vert
   enddo
!
   call attr_init(1,Ic); ATTR_VAL(1,1:Ic) = 0.d0
!
!------------------------------------------------------------------
! Step 2 : Write subdomain info
!------------------------------------------------------------------
!
   Igrid = CURRENT_GRID
!
   select case(Idx)
   case(1)
!
!..export subdomain numbers
!
      do iel=1,NRELES
         mdle = ELEM_ORDER(iel)
         call get_subd(mdle, subd)
         if (RANK .ne. subd) cycle
!
!     ...select appropriate visualization object
         ntype = NODES(mdle)%ntype
         nV = n_elem_vert(iel)
!
!     ...loop over vertices of visualization object
         do iv=1,nV
            ATTR_VAL(1,n_vert_offset(iel)+iv) = real(RANK,8)
         enddo
!  ...end loop over elements
      enddo
!
   case(2)
!
!..export number of elements in subdomain
!
   do iel=1,NRELES
      mdle = ELEM_ORDER(iel)
      call get_subd(mdle, subd)
      if (RANK .ne. subd) cycle
!
!  ...select appropriate visualization object
      ntype = NODES(mdle)%ntype
      nV = n_elem_vert(iel)
!
!  ...loop over vertices of visualization object
      do iv=1,nV
         ATTR_VAL(1,n_vert_offset(iel)+iv) = real(GRID(Igrid)%nreles_subd,8)
      enddo
!..end loop over elements
   enddo
!
   case(3)
!
!..export polynomial order
!
   do iel=1,NRELES
      mdle = ELEM_ORDER(iel)
      call get_subd(mdle, subd)
      if (RANK .ne. subd) cycle
!
!  ...select appropriate visualization object
      ntype = NODES(mdle)%ntype
      nV = n_elem_vert(iel)
!
      call decode(NODES(Mdle)%order, nord12,nord3)
      call decode(nord12, nord1,nord2)
!
!  ...loop over vertices of visualization object
      do iv=1,nV
         ATTR_VAL(1,n_vert_offset(iel)+iv) = real(nord3,8)
      enddo
!..end loop over elements
   enddo
!
   case(4)
!
!..export overlap due to ghosting
!
!  ...Initialize map to get iel from mdle
      mdle2iel = index_map_init(Grid(Igrid)%mdlel,Grid(Igrid)%nreles)
!
!  ...loop elements
      do iel = 1,NRELES
         mdle = ELEM_ORDER(iel)
!
!     ...get coarse grid element
         call find_master(Igrid,mdle, mdleC)
!
!     ...if coarse element is in my padded subdomain, mark
         if (map_index(mdle2iel,mdleC) .gt. 0) then
!        ...select appropriate visualization object
            ntype = NODES(mdle)%ntype
            nV = n_elem_vert(iel)
!
!        ...loop vertices of visualization object
            do iv=1,nV
               ATTR_VAL(1,n_vert_offset(iel)+iv) = 1.d0
            enddo
         endif
!
!  ...end loop over elements
      enddo
!
      call index_map_finalize(mdle2iel)
!
   case default
      write(*,*) 'scalar2vtk_mg: bad idx - ', idx
      stop 1
   end select
!
!------------------------------------------------------------------
!  Step 3 : Collect on host
!------------------------------------------------------------------
!
   if (DISTRIBUTED .and. .not. HOST_MESH) then
      count = Ic
      if (RANK .eq. ROOT) then
         call MPI_REDUCE(MPI_IN_PLACE,ATTR_VAL,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      else
         call MPI_REDUCE(ATTR_VAL,ATTR_VAL,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      endif
   endif
!
!..if exporting overlap, adjust
   if (idx .eq. 4 .and. RANK .eq. ROOT) then
      ATTR_VAL = ATTR_VAL - 1.d0
   endif
!
!------------------------------------------------------------------
!  Step 4 : Write to file with HDF5
!------------------------------------------------------------------
!
   if (RANK .eq. ROOT) call attr_write(Sname,len(Sname),Sfile,len(Sfile),Snick,len(Snick))
!
!------------------------------------------------------------------
!  Step 5 : Deallocate
!------------------------------------------------------------------
!
   call attr_close()
!
end subroutine scalar2vtk_mg
