!----------------------------------------------------------------------------------------
!> @brief      write vector attribute to .h5 file
!!
!> @param[in ] Sname   - name/description (e.g., 'Vector')
!> @param[in ] Sfile   - attribute file name (e.g., ../output/paraview/vector_00000.h5)
!> @param[in ] Snick   - attribute nickname to be used in .h5 file
!> @param[in ] Idx     - index identifying vector attribute
!> @param[out] Ic      - number of vertices of visualization object
!!
!> @date       Mar 2023
!----------------------------------------------------------------------------------------
!
#include "typedefs.h"
!
subroutine vector2vtk(Sname,Sfile,Snick,Idx, Ic)
!
   use data_structure3D
   use element_data
   use physics
   use upscale
   use paraview
   use MPI
   use mpi_param        , only: RANK,ROOT,NUM_PROCS
   use par_mesh         , only: DISTRIBUTED,HOST_MESH
!
   implicit none
!
   character(len=*), intent(in ) :: Sname
   character(len=*), intent(in ) :: Sfile
   character(len=*), intent(in ) :: Snick
   integer,          intent(in ) :: Idx
   integer,          intent(out) :: Ic
!
   type(vis)                           :: vis_obj
   integer                             :: ntype
   integer                             :: mdle, nflag, iel, iv, nV
   real(8), dimension(3)               :: xi,x,val
   integer, dimension(12)              :: nedge_orient
   integer, dimension(6)               :: nface_orient
   integer, dimension(19)              :: norder
   real(8), dimension(3,MAXbrickH)     :: xnod
   real(8), dimension(3,3)             :: dxdxi
!
   VTYPE, dimension(MAXEQNH,MAXbrickH) :: zdofH
   VTYPE, dimension(MAXEQNE,MAXbrickE) :: zdofE
   VTYPE, dimension(MAXEQNV,MAXbrickV) :: zdofV
   VTYPE, dimension(MAXEQNQ,MAXbrickQ) :: zdofQ
   VTYPE, dimension(  MAXEQNH  )       :: zsolH
   VTYPE, dimension(  MAXEQNH,3)       :: zgradH
   VTYPE, dimension(3,MAXEQNE  )       :: zsolE
   VTYPE, dimension(3,MAXEQNE  )       :: zcurlE
   VTYPE, dimension(3,MAXEQNV  )       :: zsolV
   VTYPE, dimension(  MAXEQNV  )       :: zdivV
   VTYPE, dimension(  MAXEQNQ  )       :: zsolQ
!
   integer :: ibeg,iattr,icomp,isol,iload,ireal,ndom
   real(8), external :: dreal_part,dimag_part
!
   integer :: VTU_data_size
!
!..OpenMP parallelization: auxiliary variables
   integer, dimension(NRELES) :: n_vert_offset, n_elem_vert
!
!..MPI
   integer :: ierr,count,subd
!
!----------------------------------------------------------------------------------------
!
!..Step 1 : Preliminary calculations (offsets, etc.)
   Ic=0
!
!..decode
   ireal = abs(Idx)/Idx
   iload = abs(Idx)/100
   iattr = abs(Idx) - iload*100 ; iattr=iattr/10
   icomp = abs(Idx) - iload*100 - iattr*10
!
!..address of 1st component for the attribute
   ibeg=ADRES(iattr)
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
   call attr_init(3,Ic); ATTR_VAL(1:3,1:Ic) = 0
!
!..Step 2 : Write attribute of interest
!
   nflag=1
!$OMP PARALLEL DO                                     &
!$OMP PRIVATE(mdle,ndom,ntype,iv,nV,xi,               &
!$OMP         xnod,zdofH,zdofE,zdofV,zdofQ,           &
!$OMP         norder,nedge_orient,nface_orient,       &
!$OMP         x,dxdxi,zsolH,zgradH,zsolE,zcurlE,      &
!$OMP         zsolV,zdivV,zsolQ,isol,val,subd)        &
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
!  ...element info for computing solution
      call nodcor(         mdle, xnod)
      call solelm(         mdle, zdofH,zdofE,zdofV,zdofQ)
      call find_elem_nodes(mdle, norder,nedge_orient,nface_orient)
!
!  ...select appropriate visualization object
      ntype = NODES(mdle)%ntype
      nV = n_elem_vert(iel)
!
!  ...loop over nodes of visualization object
      do iv=1,nV
!
!     ...visualization point accounting for VLEVEL
         call get_vis_point(vis_on_type(ntype),iv-1, xi)
!
!     ...compute element solution
         call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod, &
                      zdofH,zdofE,zdofV,zdofQ,nflag,x,dxdxi,         &
                      zsolH,zgradH,zsolE,zcurlE,zsolV,zdivV,zsolQ)
!
!     ...approximation space
         select case(D_TYPE(iattr))
!
!        -- H^1 -- (FOR GRADIENT INFO)
         case(CONTIN)
            isol = (iload-1)*NRHVAR + ibeg + icomp
!
!        ...REAL part
            if (ireal == 1) then
               val(1) = dreal_part(zgradH(isol,1))
               val(2) = dreal_part(zgradH(isol,2))
               val(3) = dreal_part(zgradH(isol,3))
!        ...IMAGINARY part
            else
               val(1) = dimag_part(zgradH(isol,1))
               val(2) = dimag_part(zgradH(isol,2))
               val(3) = dimag_part(zgradH(isol,3))
            endif
!
!        -- H(curl) --
         case(TANGEN)
            isol = (iload-1)*NREVAR + ibeg + icomp
!
!        ...REAL part
            if (ireal == 1) then
               val(1) = dreal_part(zsolE(1,isol))
               val(2) = dreal_part(zsolE(2,isol))
               val(3) = dreal_part(zsolE(3,isol))
!        ...IMAGINARY part
            else
               val(1) = dimag_part(zsolE(1,isol))
               val(2) = dimag_part(zsolE(2,isol))
               val(3) = dimag_part(zsolE(3,isol))
            endif
!
!        -- H(div) --
         case(NORMAL)
            isol = (iload-1)*NRVVAR + ibeg + icomp
!
!        ...REAL part
            if (ireal == 1) then
               val(1) = dreal_part(zsolV(1,isol))
               val(2) = dreal_part(zsolV(2,isol))
               val(3) = dreal_part(zsolV(3,isol))
!        ...IMAGINARY part
            else
               val(1) = dimag_part(zsolV(1,isol))
               val(2) = dimag_part(zsolV(2,isol))
               val(3) = dimag_part(zsolV(3,isol))
            endif
!
         case default
            write(*,*) 'vector2vtk: unexpected approximation space. stop.'
            stop
         end select
!
!     ...add value to visualization object
         ATTR_VAL(1:3,n_vert_offset(iel)+iv) = val
!
!  ...end loop over vis object vertices
      enddo
!..end loop over elements
   enddo
!$OMP END PARALLEL DO
!
!..Step 3 : Collect on host
!
   if (DISTRIBUTED .and. .not. HOST_MESH) then
      count = 3*Ic
      if (RANK .eq. ROOT) then
         call MPI_REDUCE(MPI_IN_PLACE,ATTR_VAL,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      else
         call MPI_REDUCE(ATTR_VAL,ATTR_VAL,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      endif
   endif
!
!..Step 4 : Write to file with HDF5
!
   if (RANK .eq. ROOT) then
      if (.not. VIS_VTU) then
         call attr_write(Sname,len(Sname),Sfile,len(Sfile),Snick,len(Snick))
      else
!     ...appending the attribute data in VTU file
         nV = size(ATTR_VAL,dim=2)
         VTU_data_size = nV * 3 * 8
         write(PARAVIEW_IO) VTU_data_size
         do iv = 1,nV
            write(PARAVIEW_IO) ATTR_VAL(1,iv),ATTR_VAL(2,iv),ATTR_VAL(3,iv)
         enddo
      endif
   endif
!
!..Step 5 : Deallocate
!
   call attr_close()
!
end subroutine vector2vtk
