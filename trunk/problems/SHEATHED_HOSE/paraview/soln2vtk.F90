!----------------------------------------------------------------------------------------
!> @brief write scalar attribute to .h5 file
!!
!> @param[in ] Sname   - "Scalar"
!> @param[in ] Sfile   - attribute file name (e.g., ../output/paraview/scalar_00000.h5)
!> @param[in ] Snick   - attribute's nickname to be used in .h5 file
!> @param[in ] Scenter - "Node"
!> @param[in ] Scomp   - component of stress
!                         1 = sigma_rr
!                         2 = sigma_rt
!                         3 = sigma_tt
!> @param[out] Ic      - number of vertices of visualization object
!!
!> @date Mar 16
!----------------------------------------------------------------------------------------
!
subroutine soln2vtk(Sname, Sfile, Snick, Scenter, Scomp, Ic)
!
   use data_structure3D
   use element_data
   use upscale
   use paraview
   use mpi_wrapper
   use par_mesh         , only: DISTRIBUTED,HOST_MESH
   use sheathed_isotropic_materials
!
   implicit none
   character(len=*), intent(in ) :: Sname
   character(len=*), intent(in ) :: Sfile
   character(len=*), intent(in ) :: Snick
   character(len=*), intent(in ) :: Scenter
   integer,          intent(in ) :: Scomp
   integer,          intent(out) :: Ic
!
   type(vis) :: vis_obj
   integer   :: etype
   integer                              :: mdle, nflag, iel, iv, ndom, nV
   real(8),  dimension(3)                :: xi,x
   integer, dimension(12)               :: nedge_orient
   integer, dimension(6)                :: nface_orient
   integer, dimension(19)               :: norder
   real(8),  dimension(3,MAXbrickH)      :: xnod
   real(8), dimension(MAXEQNH,MAXbrickH) :: zdofH
   real(8), dimension(MAXEQNE,MAXbrickE) :: zdofE
   real(8), dimension(MAXEQNV,MAXbrickV) :: zdofV
   real(8), dimension(MAXEQNQ,MAXbrickQ) :: zdofQ
   real(8), dimension(3,3)               :: dxdxi
   real(8), dimension(  MAXEQNH  )       :: zsolH
   real(8), dimension(  MAXEQNH,3)       :: zgradH
   real(8), dimension(3,MAXEQNE  )       :: zsolE
   real(8), dimension(3,MAXEQNE  )       :: zcurlE
   real(8), dimension(3,MAXEQNV  )       :: zsolV
   real(8), dimension(  MAXEQNV  )       :: zdivV
   real(8), dimension(  MAXEQNQ  )       :: zsolQ
!
!  ...stiffness tensor
   real(8), dimension(3,3,3,3) :: C
!
!  ...stress and displacement
   real(8), dimension(3,3) :: sigma
!
!  ...miscellaneous
   real(8)  :: val
   real(8)  :: r2
   integer :: i,j,m,n
!
!..OpenMP parallelization: auxiliary variables
   integer, dimension(NRELES) :: n_vert_offset, n_elem_vert
!..MPI
   integer :: ierr,count,subd
!
!----------------------------------------------------------------------------------------
!
!..Step 1 : Preliminary calculations (offsets, etc.)
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
      etype = NODES(mdle)%ntype
      vis_obj = vis_on_type(etype)
      n_vert_offset(iel) = Ic
      n_elem_vert(iel) = vis_obj%nr_vert
      Ic = Ic + vis_obj%nr_vert
   enddo
!
!
   call attr_init(1,Ic); ATTR_VAL(1,1:Ic) = 0
!
!..Step 2 : Write attribute of interest
!
   nflag=1
!$OMP PARALLEL DO                                     &
!$OMP PRIVATE(mdle,ndom,etype,iv,nV,xi,               &
!$OMP         xnod,zdofH,zdofE,zdofV,zdofQ,           &
!$OMP         norder,nedge_orient,nface_orient,       &
!$OMP         x,dxdxi,zsolH,zgradH,zsolE,zcurlE,      &
!$OMP         zsolV,zdivV,zsolQ,val,subd)             &
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
      etype = NODES(mdle)%ntype
      nV = n_elem_vert(iel)
!
!  ...loop over vertices of visualization object
      do iv=1,nV
!
!     ...visualization point accounting for VLEVEL
         call get_vis_point(vis_on_type(etype),iv-1, xi)
!
!     ...compute element solution
         call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod, &
                        zdofH,zdofE,zdofV,zdofQ,nflag,x,dxdxi,         &
                        zsolH,zgradH,zsolE,zcurlE,zsolV,zdivV,zsolQ)
!     ...compute radius^2
         r2 = x(2)**2+x(3)**2
!
!     ...compute stress and displacement in Cartesian coordinates
         if (ndom.eq.1) then
            !  STRESS
            call getC(x,ndom, C)
            sigma=0.d0
            do i=1,3; do j=1,3; do m=1,3; do n=1,3
            sigma(i,j) = sigma(i,j)  &
                       + C(i,j,m,n)*zgradH(m,n)
            enddo; enddo; enddo; enddo
         else
            !  STRESS
            !          ( sigma1   sigma4  sigma5 )
            !  sigma = ( sigma4   sigma2  sigma6 )
            !          ( sigma5   sigma6  sigma3 )
            sigma(1,1) = zsolQ(4)
            sigma(2,2) = zsolQ(5)
            sigma(3,3) = zsolQ(6)
            sigma(1,2) = zsolQ(7); sigma(2,1) = zsolQ(7) 
            sigma(1,3) = zsolQ(8); sigma(3,1) = zsolQ(8)
            sigma(2,3) = zsolQ(9); sigma(3,2) = zsolQ(9)
         endif
!
         select case(Scomp)
!
!        -- sigma_rr --
         case(1)
            val =   sigma(2,2)*x(2)**2  &
                   + 2*sigma(2,3)*x(2)*x(3)  &
                   +   sigma(3,3)*x(3)**2
            val = val/r2
!
!        -- sigma_rt --
         case(2)
            val = (sigma(3,3)-sigma(2,2))*x(2)*x(3)  &
                   + sigma(2,3)*(x(2)**2-x(3)**2)
            val = val/r2
!
!        -- sigma_tt --
         case(3)
            val =   sigma(2,2)*x(3)**2  &
                   - 2*sigma(2,3)*x(2)*x(3)  &
                   +   sigma(3,3)*x(2)**2
            val = val/r2
!
         endselect
!
!     ...add value to visualization object
         ATTR_VAL(1,n_vert_offset(iel)+iv) = val
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
      count = Ic
      if (RANK .eq. ROOT) then
         call MPI_REDUCE(MPI_IN_PLACE,ATTR_VAL,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      else
         call MPI_REDUCE(ATTR_VAL,ATTR_VAL,count,MPI_REAL8,MPI_SUM,ROOT,MPI_COMM_WORLD,ierr)
      endif
   endif
!
!..Step 4 : Write to file with HDF5
!
   if (RANK .eq. ROOT) call attr_write(Sname,len(Sname),Sfile,len(Sfile),Snick,len(Snick))
!
!..Step 5 : Deallocate
!
   call attr_close()
!
endsubroutine soln2vtk
