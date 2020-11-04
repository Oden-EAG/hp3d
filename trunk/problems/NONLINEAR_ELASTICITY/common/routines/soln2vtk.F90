!----------------------------------------------------------------------------------------
!> Purpose : write scalar attribute to .h5 file
!!
!> @param[in ] Sname   - name for attribute's .h5 file
!> @param[in ] Snick   - attribute's nickname to be used in .h5 file
!> @param[in ] Scenter - "Node"
!> @param[in ] Domain    - number of subdomain to be dispalyed
!                         0 = domain 1 and 2
!                         1 = domain 1
!                         2 = domain 2
!> @param[in ] Icomp   - component of DISP_ATTR to display (see common_prob_data)
!> @param[out] Ic      - number of vertices of visualization object
!!
!> @date Aug 2020
!----------------------------------------------------------------------------------------
!
#include "typedefs.h"
!
subroutine soln2vtk(Sname, Sfile, Snick, Icomp, Ic)
!
      use data_structure3D
      use element_data
      use physics          , only : DTYPE,ADRES
      use upscale
      use paraview
      use MPI              , only: MPI_COMM_WORLD,MPI_SUM,MPI_INTEGER
      use mpi_param        , only: RANK,ROOT,NUM_PROCS
      use par_mesh         , only: DISTRIBUTED,HOST_MESH
      use common_prob_data
!
      implicit none
      character(len=*), intent(in ) :: Sname
      character(len=*), intent(in ) :: Sfile
      character(len=*), intent(in ) :: Snick
      integer,          intent(in ) :: Icomp
      integer,          intent(out) :: Ic
!
      type(vis)      :: vis_obj
      character(len=4) :: etype
      integer                             :: mdle, nflag, iel, iv, ndom, nv
      real*8,  dimension(3)               :: xi,x
      integer, dimension(12)              :: nedge_orient
      integer, dimension(6)               :: nface_orient
      integer, dimension(19)              :: norder
      real*8,  dimension(3,MAXbrickH)     :: xnod
      VTYPE, dimension(MAXEQNH,MAXbrickH) :: zdofH
      VTYPE, dimension(MAXEQNE,MAXbrickE) :: zdofE
      VTYPE, dimension(MAXEQNV,MAXbrickV) :: zdofV
      VTYPE, dimension(MAXEQNQ,MAXbrickQ) :: zdofQ
      real*8, dimension(3,3)              :: dxdxi
      VTYPE, dimension(  MAXEQNH  )       :: zsolH
      VTYPE, dimension(  MAXEQNH,3)       :: zgradH
      VTYPE, dimension(3,MAXEQNE  )       :: zsolE
      VTYPE, dimension(3,MAXEQNE  )       :: zcurlE
      VTYPE, dimension(3,MAXEQNV  )       :: zsolV
      VTYPE, dimension(  MAXEQNV  )       :: zdivV
      VTYPE, dimension(  MAXEQNQ  )       :: zsolQ
!
!  ...miscellaneous
      real*8 :: val
!
!..OpenMP parallelization: auxiliary variables
   integer, dimension(NRELES) :: n_vert_offset, n_elem_vert
!..MPI
   integer :: ierr,count,subd
!
!----------------------------------------------------------------------------------------
!
! Step 0 : Clear vis
! if (Stype.eq."Scalar") 
    ! call vis_scalar_clear
! if (Stype.eq."Vector") call vis_vector_clear
!
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
      etype = NODES(mdle)%type
      vis_obj = vis_on_type(etype)
      n_vert_offset(iel) = Ic
      n_elem_vert(iel) = vis_obj%nr_vert
      Ic = Ic + vis_obj%nr_vert
    enddo

    call attr_init(1,Ic); ATTR_VAL(1,1:Ic) = 0

!
      nflag=1
!$OMP PARALLEL DO                                     &
!$OMP PRIVATE(mdle,ndom,etype,iv,nV,xi,               &
!$OMP         xnod,zdofH,zdofE,zdofV,zdofQ,           &
!$OMP         norder,nedge_orient,nface_orient,       &
!$OMP         x,dxdxi,zsolH,zgradH,zsolE,zcurlE,      &
!$OMP         zsolV,zdivV,zsolQ,val,subd)        &
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
!       element info for computing solution
      call nodcor(         mdle, xnod)
      call solelm(         mdle, zdofH,zdofE,zdofV,zdofQ)
      call find_elem_nodes(mdle, norder,nedge_orient,nface_orient)
      call find_domain(    mdle, ndom)

!
!       select appropriate visualization object
      etype = NODES(mdle)%type
      nV = n_elem_vert(iel)
!
!       loop over vertices of visualization object
      do iv=1,nv
!
!         visualization point accounting for VLEVEL
        call get_vis_point(vis_obj, iv-1, xi)
!
!         compute element solution
        call soleval(mdle,xi,nedge_orient,nface_orient,norder,xnod,zdofH,zdofE,zdofV,zdofQ,nflag, &
                     x,dxdxi,zsolH,zgradH,zsolE,zcurlE,zsolV,zdivV,zsolQ                           )
!         set component to evaluate in soldis 
        ICHOOSE_DISP = icomp
!         use soldis with rn = e_3 
!         (since we are evaluating at vertices there is no well defined normal)
        call soldis(mdle,xi,x,(/0.d0,0.d0,1.d0/),zsolH,zgradH,zsolE,zcurlE,zsolV,zdivV,zsolQ, val)
!
!         add value to visualization object
        ! if (Stype.eq."Scalar") 
        ! call vis_scalar_add_value(val)
        ! if (Stype.eq."Vector") call vis_vector_add_value(val)
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
!   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
    if (RANK .eq. ROOT) call attr_write(Sname,len(Sname),Sfile,len(Sfile),Snick,len(Snick))
!   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
!   if (RANK .eq. ROOT) write(*,300) end_time - start_time
!
!..Step 5 : Deallocate
!
   call attr_close()


!     Step 2 : write to file and clear
      ! if (Stype.eq."Scalar") then
      ! call vis_scalar_write(Stype,len(Stype),Sname,  len(Sname), &
      !                       Snick,len(Snick),Scenter,len(Scenter))
      ! call vis_scalar_clear
      ! elseif (Stype.eq."Vector") then
      !   call vis_vector_write(Stype,len(Stype),Sname,  len(Sname), &
      !                         Snick,len(Snick),Scenter,len(Scenter))
      !   call vis_vector_clear
      ! endif
!
!
end subroutine soln2vtk
