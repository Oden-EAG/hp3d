!-----------------------------------------------------------------------
!> @brief Routine performs isotropic refinement on marked elements
!> @param[in]   Nr_elem_ref - total number of marked elements
!> @param[in]   Mdle_ref    - Middle node number of the marked element
!> @param[out]  Flag_pref   - a array which indicates which of the marked elements have under 
!!                            isotropic p-refinement
!> @date May 2024
!-----------------------------------------------------------------------
subroutine Finegrid_padap(Nr_elem_ref,Mdle_ref,Flag_pref)

   use data_structure3D
   use environment     , only: QUIET_MODE
   use mpi_param       , only: ROOT, RANK
   use MPI             , only: MPI_COMM_WORLD,MPI_COMM_WORLD
!
   implicit none
!
   integer,                            intent(in)  :: Nr_elem_ref
   integer, dimension(Nr_elem_ref),    intent(in)  :: Mdle_ref
   integer, dimension(Nr_elem_ref),    intent(out) :: Flag_pref
!
   integer :: nord_new,nord,nordx,nordy,nordz,naux,pord
   integer :: mdle,iel,kref,ierr
   integer :: etype
!
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
   if (RANK .eq. ROOT)  write(*,*) 'Starting hp refining the mesh to compute fine mesh solution'
   if (RANK .eq. ROOT)  write(*,*) 'NRELES = ', NRELES
!
!..p-refinement of marked elements
   do iel = 1,nr_elem_ref
      mdle = mdle_ref(iel)
      etype = NODES(mdle)%ntype
      nord = NODES(mdle)%order
      select case(etype)
         case(MDLB)
            nord_new = nord + 111
            call decode(nord_new,naux,nordz)
            call decode(naux,nordx,nordy)
            pord = MAX(nordx,nordy,nordz)
            if (pord .gt. MAXP) then
                  go to 700
            endif 
         case default
            write(*,*) 'refine_DPG: READING UNEXPECTED ELEMENT TYPE: ',etype
            call pause
      end select
!  ...p-refine only if pord is less than MAXP : thats why we have a goto above in if statement
      call nodmod(mdle,nord_new)  
      flag_pref(iel) = 1
!  ...not performing p-ref to any element which are already at the highest polynomial order.
      700 continue
   enddo
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!..raise order of approximation on non-middle nodes by enforcing minimum rule
   call enforce_min_rule
   call par_verify
   call update_gdof
   call update_Ddof
!
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
   if (RANK .eq. ROOT) write(*,*) "Performing isotropic h-refinement of marked elements"
!..h-refinement of marked elements
   do iel = 1,nr_elem_ref
      mdle = mdle_ref(iel)
      etype = NODES(mdle)%ntype
      nord = NODES(mdle)%order
      select case(etype)
         case(MDLB)
            kref = 111
         case default
            write(*,*) 'refine_DPG: READING UNEXPECTED ELEMENT TYPE: ',etype
            call pause
      end select
      call refine(mdle,kref)
   enddo
!
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
   call close_mesh
   call enforce_min_rule
   call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!
   call par_verify
   call update_gdof
   call update_Ddof
   call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!    
end subroutine Finegrid_padap
