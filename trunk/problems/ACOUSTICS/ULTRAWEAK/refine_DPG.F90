!--------------------------------------------------------------------------
!> @brief      Marks and executes refinemnts based on the current solution
!!             and refinement type given
!!
!> @param[in]  Irefine  - refinement type
!!                            0: no refinement; only display residual/error
!!                            1: uniform refinements
!!                            2: adaptive refinemens (isotropic)
!> @param[in]  Nreflag  - refinement kind
!!                            1: h-refinement
!!                            2: p-refinement
!!                            3: hp-refinement, depending on element size
!> @param[in]  Factor   - threshold for adaptive marking
!> @param[out] Nstop    - flag indicating whether refinement increased DOFs
!!                            0: refinement increased DOFs
!!                            1: refinement did not increase DOFs
!!
!> @date       July 2023
!--------------------------------------------------------------------------
subroutine refine_DPG(Irefine,Nreflag,Factor, Nstop)
!
   use common_prob_data_UW
   use control
   use element_data
   use data_structure3D
   use environment  , only: QUIET_MODE
   use assembly_sc  , only: NRDOF_CON, NRDOF_TOT
   use parametersDPG, only: NORD_ADD
   use par_mesh     , only: DISTRIBUTED,HOST_MESH
   use mpi_param    , only: ROOT,RANK,NUM_PROCS
   use MPI          , only: MPI_COMM_WORLD,MPI_SUM,MPI_COMM_WORLD,   &
                            MPI_REAL8,MPI_INTEGER,MPI_IN_PLACE,MPI_MAX
   use sorts        , only: qsort_duplet
!
   implicit none
!
   integer, intent(in)  :: Irefine
   integer, intent(in)  :: Nreflag
   real(8), intent(in)  :: Factor
   integer, intent(out) :: Nstop
!
   integer :: Nflag(NR_PHYSA)
   integer, parameter :: max_step = 50
   integer, dimension(max_step), save :: nrdof_tot_mesh
   integer, dimension(max_step), save :: nrdof_con_mesh
   real(8), dimension(max_step), save :: residual_mesh
   real(8), dimension(max_step), save :: rate_mesh
   real(8), dimension(max_step), save :: error_mesh
   real(8), dimension(max_step), save :: rel_error_mesh
   real(8), dimension(max_step), save :: rate_error_mesh
!
   real(8) :: elem_resid(NRELES)
   integer :: elem_ref_flag(NRELES)
   integer :: nr_elem_ref
   integer, allocatable :: mdle_ref(:), elem_pref(:), elem_pref_ord(:)
   integer, allocatable :: ref_type(:), ref_kind(:)
!
!..adaptive refinements
   real(8) :: res_sum
!
   integer, save :: istep = 0
   integer, save :: irefineold = 0
!
   real(8) :: resid_subd,resid_tot,elem_resid_max
   real(8) :: errorH,errorE,errorV,errorQ,rnormH,rnormE,rnormV,rnormQ
   real(8) :: error_tot,rnorm_tot,error_subd,rnorm_subd
!
   integer :: i,ic,mdle,iel,kref,subd,count,ierr,icpref
   integer :: iref,nord_comp(3),nordx_new,nordy_new,nordz_new
   integer :: nord,nord_new
   integer :: idec_pref,refine_strategy
   real(8) :: x(3), xnod(3,8)
!
!..element type
   integer :: ntype
!
   real(8) :: MPI_Wtime,start_time,end_time
!
!..printing flag
   integer :: iprint = 1
!
!-----------------------------------------------------------------------
!
!..Set refinement strategy
!     1 : Greedy
!     2 : Doerfler
   refine_strategy = 2
!
!..field variables flag - used inside ndof_DPG and for error calc.
   nflag(1:NR_PHYSA)=(/0,0,1/)
!
!..initialize
   elem_resid(1:NRELES)    = 0.d0
   elem_ref_flag(1:NRELES) = 0
   allocate(mdle_ref(NRELES)); mdle_ref = 0
!
!..increase step if necessary
!..irefineold=0 means no refinement was made in the previous step
!..if first call or if a refinement was made, increase step
   if ((istep.eq.0).or.(irefineold.ne.0)) istep=istep+1
   irefineold=Irefine
!
!..get dof count from solver
   nrdof_tot_mesh(istep) = NRDOF_TOT
   nrdof_con_mesh(istep) = NRDOF_CON
   if (NRDOF_TOT .eq. 0) then
      istep=istep-1
      goto 70
   endif
!
!..initialize residual and error
   resid_subd = 0.d0
   error_subd = 0.d0
   rnorm_subd = 0.d0
   elem_resid_max = 0.d0
!
!..start timer
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!
!..residual/error computation
!
!$OMP PARALLEL DO                                     &
!$OMP PRIVATE(errorH,errorE,errorV,errorQ,            &
!$OMP         rnormH,rnormE,rnormV,rnormQ,            &
!$OMP         mdle,subd)                              &
!$OMP SCHEDULE(DYNAMIC)                               &
!$OMP REDUCTION(+:resid_subd,error_subd,rnorm_subd)   &
!$OMP REDUCTION(MAX:elem_resid_max)
   do iel=1,NRELES
      mdle = ELEM_ORDER(iel)
      call get_subd(mdle, subd)
      if (DISTRIBUTED .and. (RANK .ne. subd)) cycle
!
!  ...compute L^2 norm of element residual
      call elem_residual(mdle, elem_resid(iel),elem_ref_flag(iel))
!
!  ...sum element residual and update max
      resid_subd = resid_subd + elem_resid(iel)
      elem_resid_max = max(elem_resid_max,elem_resid(iel))
!
      if (NEXACT.ge.1) then
!     ...compute L^2 norm of element error
         call element_error(mdle,Nflag, errorH,errorE,errorV,errorQ,  &
                                        rnormH,rnormE,rnormV,rnormQ)
         error_subd = error_subd + errorH + errorE + errorV + errorQ
         rnorm_subd = rnorm_subd + rnormH + rnormE + rnormV + rnormQ
      endif
   enddo
!$OMP END PARALLEL DO
!
!..collect error and residual from subdomains
   resid_tot = 0.d0; error_tot = 0.d0; rnorm_tot = 0.d0
   if (DISTRIBUTED .and. (.not. HOST_MESH)) then
      count = 1
      call MPI_ALLREDUCE(resid_subd,resid_tot,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(error_subd,error_tot,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(rnorm_subd,rnorm_tot,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(MPI_IN_PLACE,elem_resid_max,count,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ierr)
      count = NRELES
      call MPI_ALLREDUCE(MPI_IN_PLACE,elem_resid,count,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
   else
      resid_tot = resid_subd
      error_tot = error_subd
      rnorm_tot = rnorm_subd
   endif
!
!..end timer
   call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
   if (.not. QUIET_MODE) then
      if (RANK .eq. ROOT) write(*,1010) end_time-start_time
 1010 format('elem_residual  : ',f12.5,'  seconds')
   endif
!
   if (RANK .ne. ROOT) goto 70
!
!..update and display convergence history
   residual_mesh(istep) = sqrt(resid_tot)
   if (NEXACT.ge.1) then
      error_mesh(istep)     = sqrt(error_tot)
      rel_error_mesh(istep) = sqrt(error_tot/rnorm_tot)
   endif
!
!..compute decrease rate for the residual and error
   select case(istep)
   case(1)
      rate_mesh(istep) = 0.d0
      if (NEXACT.ge.1) rate_error_mesh(istep) = 0.d0
   case default
      rate_mesh(istep) =  &
      log(residual_mesh(istep-1)/residual_mesh(istep))/  &
      log(float(nrdof_tot_mesh(istep-1))/float(nrdof_tot_mesh(istep)))
      if (NEXACT.ge.1) then
         rate_error_mesh(istep) = &
         log(error_mesh(istep-1)/error_mesh(istep))/  &
         log(float(nrdof_tot_mesh(istep-1))/float(nrdof_tot_mesh(istep)))
      endif
   end select
!
!..print out the history of refinements
   write(*,*)
   write(*,*) 'HISTORY OF REFINEMENTS'
   if (NEXACT.eq.0) write(*,7005)
   if (NEXACT.gt.0) write(*,7006)
 7005  format(' mesh |',' nrdof_tot |',' nrdof_con |','    residual   |','   residual rate  ',/)
 7006  format(' mesh |',' nrdof_tot |',' nrdof_con |','    residual   |','   residual rate  |', &
              ' field error  |','rel field error|','   error rate ',/)
!
   do i=1,istep
      if (NEXACT.eq.0) then
         write(*,7003) i,nrdof_tot_mesh(i),nrdof_con_mesh(i),residual_mesh(i),rate_mesh(i)
 7003    format(2x,i2,'  | ',2(i9,' | '),es12.5,'  |',f7.2)
      else
         write(*,7004) i,nrdof_tot_mesh(i),nrdof_con_mesh(i),residual_mesh(i),rate_mesh(i), &
                       error_mesh(i),rel_error_mesh(i),rate_error_mesh(i)
 7004    format(2x,i2,'  | ',2(i9,' | '),es12.5,'  |',f7.2,'          ', &
                2(' | ',es12.5),'  |',f7.2)
      endif
      if (i .eq. istep) write(*,*)
   enddo
!
   70 continue
!
!..mark elements for adaptive refinement
   if (Irefine .eq. IADAPTIVE) then
      nr_elem_ref = 0
      if (refine_strategy .eq. 1) then
!
!     Strategy 1 (greedy strategy)
!
         do iel=1,NRELES
            mdle = ELEM_ORDER(iel)
            if (elem_resid(iel) > Factor*elem_resid_max) then
               nr_elem_ref = nr_elem_ref + 1
               mdle_ref(nr_elem_ref) = mdle
            endif
 1020       format(A,I2,A,I9,A,es12.5)
         enddo
      else
!
!     Strategy 2 (Doerfler's strategy)
!
!     ...sort elements (iel_ref) according to residual values
         mdle_ref(1:NRELES) = ELEM_ORDER(1:NRELES)
         call qsort_duplet(mdle_ref,elem_resid,NRELES,1,NRELES)
!
!     ...verify the array is sorted
         do iel = 1,NRELES-1
            if (elem_resid(iel) .lt. elem_resid(iel+1)) then
               write(*,1021) char(9), 'elem_resid not sorted: elem_resid(',iel  ,') = ',elem_resid(iel), &
                                         ', elem_resid(',iel+1,') = ',elem_resid(iel+1)
1021           format(A,I9,A,es12.5,A,I9,A,es12.5,/)
               stop
            endif
         enddo
!
         res_sum = 0.d0
         do iel = 1,NRELES
            mdle = mdle_ref(iel)
!
1022        format(A,I2,A,I9,A,es12.5)
            nr_elem_ref = nr_elem_ref + 1
            res_sum = res_sum + elem_resid(iel)
            if (res_sum > Factor*resid_tot) exit
         enddo
!
         if (RANK.eq.ROOT) write(*,1023) char(9), 'nr_elem_ref = ', nr_elem_ref
1023     format(A,A,I7)
      endif
   endif
!
   Nstop = 0
!
!..decide refinement types
   select case(Irefine)
!  ...uniform refinements
      case(IUNIFORM)
         continue
!
!  ...adaptive refinements
      case(IADAPTIVE)
!
!     ...allocate refinement arrays
         allocate(ref_type(nr_elem_ref)); ref_type(:) = 0
         allocate(ref_kind(nr_elem_ref)); ref_kind(:) = 0
!
!     ...get refinement kind
         do iel = 1,nr_elem_ref
            mdle = mdle_ref(iel)
!
!        ...Select_refinements can only be called on SUBD
!           Nreflag=1,2 doesn't need to cycle, but communication is cheap anyway
            call get_subd(mdle, subd)
            if (DISTRIBUTED .and. (RANK .ne. subd)) cycle
!
            ntype = NODES(mdle)%ntype
!
!        ...set potential h-ref
            select case(ntype)
               case(MDLB)
                  kref = 111 ! iso
               case(MDLP)
                  kref = 11  ! iso
               case default
                  write(*,*) char(9), 'refine_DPG: READING UNEXPECTED ELEMENT TYPE: ',ntype
                  stop 1
            end select
!
!        ...set potential p-ref
            nord = NODES(mdle)%order
            select case(ntype)
            case(MDLB)
               call decod(nord,10,3,nord_comp(1:3))
               nordx_new = min(nord_comp(1)+1,MAXP)
               nordy_new = min(nord_comp(2)+1,MAXP)
               nordz_new = min(nord_comp(3)+1,MAXP)
               nord_new = nordx_new*100+nordy_new*10+nordz_new
            case default
               write(*,*) 'refine_DPG: UNFINISHED' ; stop 1
            end select
!
!        ...if element already at max order, h-refine instead
            if (nord_new.eq.nord) then
               idec_pref=0
            else
               idec_pref=1
            endif
!
            select case(Nreflag)
            case(1)
               ref_type(iel) = 1
               ref_kind(iel) = kref
            case(2)
               ref_type(iel) = 2
               ref_kind(iel) = nord_new
            case(3,4)
               call select_refinement(mdle,iref)
               select case(iref)
!           ...h-refinement selected
               case(1)
                  ref_type(iel) = 1
                  ref_kind(iel) = kref
!           ...p-refinement selected
               case(2)
                  if (idec_pref.eq.1) then
                     ref_type(iel) = 2
                     ref_kind(iel) = nord_new
!              ...p-refinement selected but already at max order
                  else
                     ref_type(iel) = 1
                     ref_kind(iel) = kref
                  endif
               end select
            end select
         enddo
!
!     ...communicate refinements
         if (DISTRIBUTED .and. nr_elem_ref .gt. 0) then
            count = nr_elem_ref
            call MPI_ALLREDUCE(MPI_IN_PLACE,ref_type,count,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
            call MPI_ALLREDUCE(MPI_IN_PLACE,ref_kind,count,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
         endif
!
      case default; Nstop = 1
!
   end select
!
!-----------------------------------------------------------------------
!                         REFINE AND UPDATE MESH
!-----------------------------------------------------------------------
!
   Nstop = 0
!..use appropriate strategy depending on refinement type
   select case(Irefine)
!  ...uniform refinements
      case(IUNIFORM)
!     ...isotropic href
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
!
         select case (Nreflag)
         case(1)
            call global_href
         case(2)
            call global_pref
         case default
            call global_href
         end select
!
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
         if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2020) end_time-start_time
         call update_gdof
         call update_Ddof
!
!  ...adaptive refinements
      case(IADAPTIVE)
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
         if (RANK .eq. ROOT)  write(*,*) char(9), 'Starting adaptive refinements...'
         if (RANK .eq. ROOT)  write(*,*) char(9), 'NRELES = ', NRELES
!
!     ...arrays to buffer p-refinements
         allocate (elem_pref(nr_elem_ref))
         allocate (elem_pref_ord(nr_elem_ref))
         icpref = 0
!
!     ...execute h-refinements, buffer prefinements
         do iel = 1,nr_elem_ref
            mdle = mdle_ref(iel)
!
            select case(ref_type(iel))
            case(0)
               cycle
            case(1)
               kref = ref_kind(iel)
               call refine(mdle,kref)
            case(2)
               icpref = icpref + 1
               elem_pref(icpref) = mdle
!           ...not currently used, execute_pref does uniform p ref
               elem_pref_ord(icpref) = nord_new
            case default
               if (RANK.eq.ROOT) write(*,*) 'Elem:', mdle, 'not refined, ref_type = ',ref_type(iel)
            end select
         enddo
!
!     ...deallocate refinement type arrays
         deallocate(ref_type,ref_kind)
!
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
         if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2020) end_time-start_time
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
         call close_mesh
         call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
         if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2025) end_time-start_time
!
         if (icpref.gt.0) then
            call execute_pref(elem_pref(1:icpref), icpref)
            call enforce_min_rule
         endif
!
         deallocate(elem_pref,elem_pref_ord)
!
         call par_verify
         call update_gdof
         call update_Ddof
!
         if (RANK .eq. ROOT) write(*,*) char(9), 'NRELES = ', NRELES
         if (RANK .eq. ROOT) write(*,*) char(9), 'Finished adaptive refinements...'
!
      case default; Nstop = 1
!
   end select
!
  2020 format(' refinement : ',f12.5,'  seconds')
  2025 format(' close mesh : ',f12.5,'  seconds')
  2030 format(A,I8,', ',I9)
!
   if (IBC_PROB.eq.3 .or. IBC_PROB.eq.4 .or. IBC_PROB.eq.6) call propagate_flag(2,3)
!
end subroutine refine_DPG





!--------------------------------------------------------------------------
!> @brief      Select between h- and p- refinements by comparing element
!!             size to wavelength
!!
!> @param[in]  Mdle - middle node
!> @param[out] Iref - refinement kind
!!                      1: h-refinement
!!                      2: p-refinement
!!
!> @date       July 2023
!--------------------------------------------------------------------------
   subroutine select_refinement(Mdle,Iref)
!
      use data_structure3D
      use common_prob_data_UW
      use element_data
!
      implicit none
!
      integer, intent(in)  :: Mdle
      integer, intent(out) :: Iref
!
      real(8) :: hmin, lambda
      integer :: ntype, nrv, nre, loc
      logical :: found
!
!  ...for cavity problem; list of singular edges and vertices
      integer, parameter :: nr_vert = 8
      integer, parameter :: nr_edge = 12
      integer, parameter :: vert_singular(8)  = (/ 562,563,571,570,634,626,627,635 /)
      integer, parameter :: edge_singular(12) = (/ 1587,1568,1590,1591,     &
                                                   1415,1414,1411,1392,     &
                                                   1589,1592,1569,1566 /)
!
!---------------------------------------------------------------------
!
      ntype = NODES(Mdle)%ntype
!
!  ...get element nodes
      call elem_nodes(Mdle, nodesl,norientl)
      nrv = nvert(ntype)
      nre = nedge(ntype)
!
!  ...look for a singular vertex
      found = .false.
!
      if (PROB_KIND.eq.PROB_CAVITY)
!
!     ...check whether any element vertices are singular
         do iv=1,nrv
            call locate(nodesl(iv),vert_singular,nr_vert, loc)
            if (loc.ne.0) then
               found = .true.
               exit
            endif
         enddo

         if (.not. found) then
!        ...look for a singular edge or its father
            do ie = 1, nre
               nod = nodesl(nrv+ie)
!
!           ...find father
               call locate_father(nod,nodp)
!
!           ...if not edge node, skip
               if (NODES(nodp)%ntype .ne. MEDG) cycle
!
!           ...check whether father is singular
               call locate(nodp,edge_singular,nr_edge, loc)
               if (loc.ne.0) then
                  found = .true.
                  exit
               endif
            enddo
         endif
      endif
!
      call find_hmin(Mdle,hmin)
!
      lambda = OMEGA/(2.0d0*PI-1.d-9)
!
!  ...if less than 2 element per wave, or if singular element; h-refine
      if (hmin.gt. 0.5d0/lambda .or. found) then
         Iref=1
      else
         Iref=2
      endif
!
   end subroutine select_refinement





!--------------------------------------------------------------------------
!> @brief      Find initial-mesh ancestor of a node
!!
!> @param[in]  nod   - node
!> @param[out] nodp  - initial mesh ancestor of node
!!
!> @date       July 2023
!--------------------------------------------------------------------------
   subroutine locate_father(nod,nodp)
!
      use data_structure3D
!
      implicit none
!
      integer, intent(in) :: nod
      integer, intent(out):: nodp
      integer             :: nods
!
!-------------------------------------------------------
!
      nodp = nod
      nods = nod
!
      do while (NODES(nods)%father .gt. 0)
         nodp = NODES(nods)%father
         nods  = nodp
      enddo
!
   end subroutine locate_father
   




!--------------------------------------------------------------------------
!> @brief      Find minimum edge length of a hexahedral element
!!
!> @param[in]  Mdle - middle node
!> @param[out] Iref - refinement kind
!!                      1: h-refinement
!!                      2: p-refinement
!!
!> @date       July 2023
!--------------------------------------------------------------------------
   subroutine find_hmin(Mdle,Hmin)
!
      use data_structure3D

      implicit none
!
      integer, intent(in)  :: Mdle
      real*8,  intent(out) :: Hmin
!
!  ...geometry dof
      real(8) :: xnod(3,MAXbrickH), h1, h2, h3
      real(8) :: dx1(3), dx2(3), dx3(3)
      integer :: i
!
!--------------------------------------------------------------------------
!
      if (NODES(mdle)%ntype.ne.MDLB) then
         write(*,*) 'find_hmin: only implemented for hexahedral elements so far'
         stop 1
      endif
!
!  ...determine min element size
      call nodcor(Mdle,xnod)
!
      dx1(:) = xnod(:,2)-xnod(:,1)
      dx2(:) = xnod(:,4)-xnod(:,1)
      dx3(:) = xnod(:,5)-xnod(:,1)
!
      h1 = dsqrt(dx1(1)**2 + dx1(2)**2 + dx1(3)**2)
      h2 = dsqrt(dx2(1)**2 + dx2(2)**2 + dx2(3)**2)
      h3 = dsqrt(dx3(1)**2 + dx3(2)**2 + dx3(3)**2)
!
      Hmin = min(h1,h2,h3)
!
   end subroutine find_hmin
