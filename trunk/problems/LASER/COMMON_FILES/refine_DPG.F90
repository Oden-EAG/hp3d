!--------------------------------------------------------------------
!
!     routine name      - refine_DPG
!
!--------------------------------------------------------------------
!
!     latest revision:  - Sept 17
!
!     purpose:          - refines elements assuming problem has
!                         already been solved. If uniform refinements
!                         it refines everything. Otherwise it follows
!                         greedy strategy based on residual to
!                         determine which elements to refine.
!                         It also displays dof, residuals and
!                         relevant errors at each step.
!
!     arguments:
!
!     in:
!             Irefine   = 0 if no refinement -> only display
!                       = 1 if uniform refinements
!                       = 2 if anisotropic refinements (provide kref)
!                       = 3 if adaptive refinements
!             Nreflag   = 1 for h-refinements
!                       = 2 for p-refinements
!                       = 3 for h- or p-refinement, depending
!                           upon the element size
!             Factor    - if element error \ge Factor*max_error
!                         then the element is refined
!             nflag     - integer vector of length NR_PHYSA to indicate
!                         which PHYS Attribute one should compute the error for
!             physNick  -  nickname for PHYS Attribute:
!                          physNick = HEVQ, for e.g., physNick = 1001 means
!                          1 H and 1 Q variable etc.
!             ired      - logical, .true. if compute residual is true
!             iso_ans   - integer to define how many anisotropic refinement
!     out:
!             Nstop     = 1 if refinement did not increase dof
!
!---------------------------------------------------------------------
!
subroutine refine_DPG(Irefine,Nreflag,Factor,nflag,physNick,ires,iso_ans, Nstop)
!
   use CommonParam
   use control
   use data_structure3D
   use environment, only : QUIET_MODE
   use m_assembly, ONLY: NRDOF_CON, NRDOF_TOT
   use parametersDPG, ONLY : NORD_ADD
!
   implicit none
!
   integer, intent(in)  :: Irefine
   integer, intent(in)  :: Nreflag
   real*8,  intent(in)  :: Factor
   integer, dimension(NR_PHYSA), intent(in) :: nflag
   integer, intent(in)  :: physNick
   integer, intent(in)  :: iso_ans
   logical, intent(in)  :: ires
   integer, intent(out) :: Nstop
!
   integer, parameter :: max_step = 300
   integer, dimension(max_step), save :: nrdof_tot_mesh
   integer, dimension(max_step), save :: nrdof_con_mesh
   real*8,  dimension(max_step), save :: residual_mesh
   real*8,  dimension(max_step), save :: rate_mesh
   real*8,  dimension(max_step), save :: error_mesh
   real*8,  dimension(max_step), save :: rel_error_mesh
   real*8,  dimension(max_step), save :: rate_error_mesh

!
   real*8,  allocatable :: elem_resid(:)
   integer, allocatable :: elem_ref_flag(:)
   integer, allocatable :: list_elem(:)
   integer, allocatable :: elem_pref(:), elem_pref_ord(:)
   integer, allocatable :: list_ref_flags(:)
!
   integer, save :: istep = 0
   integer, save :: irefineold = 0
   integer :: nrdofField, nrdofDPG
   integer :: nr_elem_to_refine
   integer, dimension(NRELES)   :: N_elem
   real*8 :: residual, elem_resid_max, err, rnorm
   real*8 :: errorH,errorE,errorV,errorQ,rnormH,rnormE,rnormV,rnormQ
   integer :: i,iprint,ic,mdle,iel,kref
   integer, dimension(3) :: vref, jref
   integer :: nord, nordx, nordy, nordz, nordyz
   integer :: nordx_new, nordy_new, nordz_new, nord_new
   integer :: idec_pref,idec_href, iref, enforce_flag, icpref
   integer :: ref_xyz,iii
!..auxiliary variables for timing
   real*8 :: start, OMP_get_wtime
!
!-----------------------------------------------------------------------
!
!..initialize
   iprint = 0
   allocate (elem_resid(NRELES),elem_ref_flag(NRELES))
   allocate (list_elem(NRELES),list_ref_flags(NRELES))
!
!..increase step if necessary
!..irefineold=0 means no refinement was made in the previous step
!..if first call or if a refinement was made, increase step
   if ((istep.eq.0).or.(irefineold.ne.0)) istep=istep+1
   irefineold=Irefine
!
!..get dof count from solver
!  call ndof_DPG !calculates nrdofDPG and nrdofField
   nrdof_tot_mesh(istep) = NRDOF_TOT
   nrdof_con_mesh(istep) = NRDOF_CON
!
!..create list of mdle nods numbers
   mdle=0
   do iel=1,NRELES
      call nelcon(mdle, mdle)
      N_elem(iel) = mdle
   enddo
!
!..initialize global residual and error
   residual       = 0.d0
   elem_resid_max = 0.d0
   err   = 0.d0
   rnorm = 0.d0
!
!..start timer
   start = OMP_get_wtime()
!
!..residual computation
!
!$OMP PARALLEL DO                                     &
!$OMP PRIVATE(errorH,errorE,errorV,errorQ,            &
!$OMP         rnormH,rnormE,rnormV,rnormQ)            &
!$OMP SCHEDULE(DYNAMIC)                               &
!$OMP REDUCTION(+:residual,err,rnorm)
   do iel=1,NRELES
      if(ires) then
        call elem_residual(N_elem(iel), elem_resid(iel),elem_ref_flag(iel))
        residual = residual + elem_resid(iel)
      endif
      if (NEXACT.ge.1) then
         call element_error(N_elem(iel),nflag, errorH,errorE,errorV,errorQ,  &
                                               rnormH,rnormE,rnormV,rnormQ)
         select case(physNick)
            case(1)
               err   = err   + errorQ
               rnorm = rnorm + rnormQ
            case(10)
               err   = err   + errorV
               rnorm = rnorm + rnormV
            case(100)
               err   = err   + errorE
               rnorm = rnorm + rnormE
            case(1000)
               err   = err   + errorH
               rnorm = rnorm + rnormH
            case(1001)
               err   = err   + errorH+ errorQ
               rnorm = rnorm + rnormH+ rnormQ
            case default
               err   = err   + errorH + errorE + errorV + errorQ
               rnorm = rnorm + rnormH + rnormE + rnormV + rnormQ
         end select
      endif
   enddo
!$OMP END PARALLEL DO
!
!..end timer
   if (.not. QUIET_MODE) then
      write(*,10) OMP_get_wtime()-start
 10   format(' elem_residual  : ',f12.5,'  seconds')
   endif
!
!..store maximum element residual
   elem_resid_max = maxval(elem_resid)
!
!..update and display convergence history
   residual_mesh(istep) = sqrt(residual)
   if (NEXACT.ge.1) then
      error_mesh(istep)     = sqrt(err)
      rel_error_mesh(istep) = sqrt(err/rnorm)
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
   write(*,*)
   write(*,*) 'HISTORY OF REFINEMENTS'
   if (NEXACT.eq.0) write(*,7005)
   if (NEXACT.gt.0) write(*,7006)
 7005  format(' mesh |',' nrdof_tot|',' nrdof_con|','    residual   |','   residual rate  ')
 7006  format(' mesh |',' nrdof_tot|',' nrdof_con|','    residual   |','   residual rate  |',   &
              ' field error  |','rel field error|','   error rate ')
   write(*,*)
!
   do i=1,istep
      if (NEXACT.eq.0) then
         write(*,*)
         write(*,7003) i,nrdof_tot_mesh(i),nrdof_con_mesh(i), &
                         residual_mesh(i),rate_mesh(i)
 7003 format(2x,i2,'  | ',2(i8,' | '),e12.5,'  |',f7.2)
      else
      write(*,7004) i,nrdof_tot_mesh(i),nrdof_con_mesh(i),       &
                    residual_mesh(i),rate_mesh(i),error_mesh(i), &
                    rel_error_mesh(i),rate_error_mesh(i)
 7004 format(2x,i2,'  | ',2(i8,' | '),e12.5,'  |',f7.2,'          ', &
                 2(' | ',e12.5),'  |',f7.2)
      endif
   enddo
!
!-----------------------------------------------------------------------
!                         REFINE AND UPDATE MESH
!-----------------------------------------------------------------------
!
!..use appropriate strategy depending on refinement type
!  NOTE: if Irefine=0 nothing is done, so only relevant info (the
!        code above this line) is displayed.
   select case(Irefine)
!  ...uniform refinements
      case(IUNIFORM)
         call global_href
         call close_mesh
         call update_gdof
         call update_Ddof
         nr_elem_to_refine = NRELES
         write(*,*)'NRELES is', NRELES
!   ...anisotropic refinements
      case(IANISOTROPIC)
 30      write(*,*)'set direction to refine: 1 - xy direction, 2 - z direction'
         !read(*,*), ref_xyz
         ref_xyz = 2
         if((iso_ans.le.-1).or.(ref_xyz.ne.1).and.(ref_xyz.ne.2)) then
            write(*,*) 'iso_ans,ref_xyz are', iso_ans,ref_xyz
            write(*,*) 'error in anisotropic refinements: enter valid input'
            go to 30
         endif
         if(iso_ans.ne.-1) then
            do iii=1,iso_ans
            call setAnisoRef(ref_xyz)
         enddo
         call close_mesh
         call update_gdof
         call update_Ddof
         else
            do while(1.gt.0)
               call setAnisoRef(ref_xyz)
            enddo
         call close_mesh
         call update_gdof
         call update_Ddof
         endif
         nr_elem_to_refine = NRELES
         write(*,*)'NRELES is', NRELES
!
!  ...adaptive refinements
      case(IADAPTIVE)
!     ...use the greedy strategy to determine which elements to refine
         ic=0
         mdle=0
         do iel=1,NRELES
            call nelcon(mdle, mdle)
            if (elem_resid(iel).ge.Factor*elem_resid_max) then
               ic=ic+1
               list_elem(ic) = mdle
               list_ref_flags(ic) = elem_ref_flag(iel)
               iprint = 1
               if (iprint.eq.1) then
                  write(*,7037) ic,mdle
 7037             format('refine_DPG: ADDING TO THE LIST: ic,mdle = ',2i6)
               endif
               iprint = 0
            endif
         enddo
         nr_elem_to_refine = ic
         allocate (elem_pref(nr_elem_to_refine))
         allocate (elem_pref_ord(nr_elem_to_refine))
         icpref=0
!
         if (nr_elem_to_refine.gt.0) then
!        ...refine the elements from the list
            do iel=1,nr_elem_to_refine
               mdle = list_elem(iel)
               kref = list_ref_flags(iel)
!
               call decode(kref,jref,vref(3))
               call decode(jref,vref(1),vref(2))
!
!           ...restrictions on h-refinements
               idec_href=1
!           ...p-refinements limited by the maximum order supported by the code
               nord = NODES(mdle)%order
               select case(NODES(mdle)%type)
                  case('mdlb')
                     call decode(nord,nordyz,nordz)
                     call decode(nordyz,nordx,nordy)
                     nordx_new = nordx; nordy_new = nordy ; nordz_new = nordz
                     if (vref(1).eq.1) then
                        nordx_new = min(nordx+1,MAXP-NORD_ADD)
                     endif
                     if (vref(2).eq.1) then
                        nordy_new = min(nordy+1,MAXP-NORD_ADD)
                     endif
                     if (vref(3).eq.1) then
                        nordz_new = min(nordz+1,MAXP-NORD_ADD)
                     endif
                     nord_new =nordx_new*100+nordy_new*10+nordz_new
                  case default
                     write(*,*) 'refine_DPG: UNFINISHED' ; stop 1
               end select
               if (nord_new.eq.nord) then
                  idec_pref=0
               else
                  idec_pref=1
               endif
!
               select case(Nreflag)
                  case(1)
                     if (idec_href.eq.1) then
                        call refine(mdle,kref)
                     endif
                  case(2)
                     if (idec_pref.eq.1) then
                        call nodmod(mdle,nord_new)
                        ! icpref=icpref+1
                        ! elem_pref(icpref) = mdle
                        ! elem_pref_ord(icpref) = nord_new
                     endif
                  case(3)
                     call select_refinement(mdle,iref)
                     select case(iref)
!
!              ...h-refinement selected
                  case(1)
                     if (idec_href.eq.1) then
                        call refine(mdle,kref)
                     endif
!
!              ...p-refinement selected
                  case(2)
                     if (idec_pref.eq.1) then
                        call nodmod(mdle,nord_new)
                        ! icpref=icpref+1
                        ! elem_pref(icpref) = mdle
                        ! elem_pref_ord(icpref) = nord_new
                     else
                        if (idec_href.eq.1) then
                           call refine(mdle,kref)
                        endif
                     endif
               end select
            end select
         enddo
!
!     ...close the mesh
         write(*,*) 'Calling close'
         call close_mesh
         write(*,*) 'Calling enforce_max_rule'
         call enforce_max_rule
         call update_gdof
         call update_ddof
         nr_elem_to_refine = NRELES
         write(*,*)'NRELES is', NRELES
!
         ! if (icpref.gt.0) call execute_pref(elem_pref(1:icpref), icpref)
         ! if (icpref.gt.0) &
         ! call perform_pref(elem_pref(1:icpref),elem_pref_ord(1:icpref), icpref)
         ! call enforce_min_rule
!
!     ...update geometry and Dirichlet flux dof after the refinements
         ! write(*,*) 'Calling update_gdof'
         ! call update_gdof
!
!     ...if case of impedance condition there is no need to update Dirichlet dof
!        if (IBC_PROB .ne. BC_IMPEDANCE) then
!           write(*,*) 'Calling update_ddof'
!           call update_Ddof
!        endif
         ! write(*,*) 'Done'
      endif
      deallocate(elem_pref, elem_pref_ord)
!
   end select
!
!-----------------------------------------------------------------------
!                   DETERMINE IF MESH WAS REFINED
!-----------------------------------------------------------------------
!..determine if mesh was refined
   if (Irefine.gt.0) then
      if (nr_elem_to_refine.eq.0) then
         Nstop = 1 ! no elements were refined
      else
         Nstop = 0 ! as expected, elements were refined
      endif
   endif
!
!-----------------------------------------------------------------------
!                              FINALIZE
!-----------------------------------------------------------------------
!
   deallocate(elem_resid)
   deallocate(elem_ref_flag)
   deallocate(list_elem)
   deallocate(list_ref_flags)
!
end subroutine refine_DPG
!
!
!--------------------------------------------------------------------
!
!     routine name      - select_refinement
!
!--------------------------------------------------------------------
!
!     latest revision:  - Aug 17
!
!     purpose:          - criterion set by the user for choosing
!                         between h and p refinement
!
!     arguments:
!
!     in:
!             mdle      - element node number
!     out:
!             iref      = 1 for h-refinement
!                       = 2 for p-refinement
!
!---------------------------------------------------------------------
!
subroutine select_refinement(Mdle,Iref)
!
   use data_structure3D
   use CommonParam
!
   implicit none
!
   integer, intent(in)  :: Mdle
   integer, intent(out) :: Iref
!..geometry dof
   real*8  :: xnod(3,MAXbrickH), hmin
!
!---------------------------------------------------------------------
!
!..determine nodes coordinates
!
   call nodcor(Mdle, xnod)
!
!..determine min element size
   hmin = min(abs(xnod(1,2)-xnod(1,1)),abs(xnod(2,4)-xnod(2,1)),  &
               abs(xnod(3,5)-xnod(3,1)))
!
   if (hmin.gt. 6.0d0/OMEGA) then
      Iref=1
   else
      Iref=2
   endif
!
!
end subroutine select_refinement
