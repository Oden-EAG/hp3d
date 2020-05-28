!--------------------------------------------------------------------
!
!     routine name      - refine_DPG
!
!--------------------------------------------------------------------
!
!     latest revision:  - Apr 2019
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
!             Nflag     - integer vector of length NR_PHYSA to indicate
!                         which PHYS Attribute one should compute the error for
!             PhysNick  - nickname for PHYS Attribute: physNick = HEVQ
!                         e.g., physNick = 1001 means 1H1 and 1L2 variable, etc.
!             Ires      - logical, .true. if compute residual is true
!             Iso_ans   - integer to define how many anisotropic refinement
!     out:
!             Nstop     = 1 if refinement did not increase dof
!
!---------------------------------------------------------------------
!
subroutine refine_DPG(Irefine,Nreflag,Factor,Nflag,PhysNick,Ires,Iso_ans, Nstop)
!
   use commonParam
   use control
   use data_structure3D
   use environment,   only: QUIET_MODE
   use assembly_sc,   only: NRDOF_CON, NRDOF_TOT
   use parametersDPG, only: NORD_ADD
!
   implicit none
!
   integer, intent(in)  :: Irefine
   integer, intent(in)  :: Nreflag
   real(8), intent(in)  :: Factor
   integer, intent(in)  :: Nflag(NR_PHYSA)
   integer, intent(in)  :: PhysNick
   integer, intent(in)  :: Iso_ans
   logical, intent(in)  :: Ires
   integer, intent(out) :: Nstop
!
   integer, parameter :: max_step = 300
   integer, dimension(max_step), save :: nrdof_tot_mesh
   integer, dimension(max_step), save :: nrdof_con_mesh
   real(8), dimension(max_step), save :: residual_mesh
   real(8), dimension(max_step), save :: rate_mesh
   real(8), dimension(max_step), save :: error_mesh
   real(8), dimension(max_step), save :: rel_error_mesh
   real(8), dimension(max_step), save :: rate_error_mesh
!
   real(8), allocatable :: elem_resid(:)
   integer, allocatable :: elem_ref_flag(:)
   integer, allocatable :: list_elem(:)
   integer, allocatable :: elem_pref(:), elem_pref_ord(:)
   integer, allocatable :: list_ref_flags(:)
!
   integer, save :: istep = 0
   integer, save :: irefineold = 0
!
   integer :: nrdofField, nrdofDPG
   integer :: nr_elem_to_refine
   integer :: n_elem(NRELES)
   real(8) :: residual, elem_resid_max, err, rnorm
   real(8) :: errorH,errorE,errorV,errorQ,rnormH,rnormE,rnormV,rnormQ
   integer :: i,ic,mdle,iel,kref,iii
   integer :: vref(3), jref(3)
   integer :: nord, nordx, nordy, nordz, nordyz
   integer :: nordx_new, nordy_new, nordz_new, nord_new
   integer :: idec_pref,idec_href, iref, enforce_flag, icpref
!
!..geometry dof (work space for nodcor)
   real(8) :: xnod(3,MAXbrickH)
   real(8) :: maxz,minz,midz
!
!..element type
   character(len=4) :: etype
!
!..auxiliary variables for timing
   real(8) :: start, OMP_get_wtime
!
!..printing flag
   integer :: iprint = 0
!
!-----------------------------------------------------------------------
!
!..initialize
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
      n_elem(iel) = mdle
   enddo
!
!..initialize global residual and error
   residual       = 0.d0
   elem_resid     = 0.d0
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
!$OMP         rnormH,rnormE,rnormV,rnormQ,            &
!$OMP         mdle,etype,xnod,maxz,minz,midz)         &
!$OMP SCHEDULE(DYNAMIC)                               &
!$OMP REDUCTION(+:residual,err,rnorm)
   do iel=1,NRELES
      mdle = n_elem(iel)
      if (USE_PML) then
         xnod = ZERO
         call nodcor(Mdle, xnod)
         etype = NODES(Mdle)%type
         select case(etype)
            case('mdlb')
               maxz = maxval(xnod(3,1:8))
               minz = minval(xnod(3,1:8))
            case('mdlp')
               maxz = maxval(xnod(3,1:6))
               minz = minval(xnod(3,1:6))
            case default
               write(*,*) 'refine_DPG: invalid etype param. stop.'
               stop
         end select
         !midz = minz + (maxz-minz)/2.d0
         if (maxz .gt. PML_REGION) cycle
      endif
      if(ires) then
         call elem_residual(mdle, elem_resid(iel),elem_ref_flag(iel))
         residual = residual + elem_resid(iel)
      endif
      if (NEXACT.ge.1) then
         call element_error(mdle,nflag, errorH,errorE,errorV,errorQ,  &
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
               err   = err   + errorH + errorQ
               rnorm = rnorm + rnormH + rnormQ
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
   do iel=1,NRELES
      mdle = n_elem(iel)
      call find_domain(mdle, i)
      if(elem_resid(iel)>0.75d0*elem_resid_max) then
         !write(*,2020) 'ndom=',i,', mdle=',mdle,', r =',elem_resid(iel)
      endif
 2020 format(A,I2,A,I5,A,es12.5)
   enddo
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
   write(*,*); write(*,*)
   write(*,*) 'HISTORY OF REFINEMENTS'
   if (NEXACT.eq.0) write(*,7005)
   if (NEXACT.gt.0) write(*,7006)
 7005  format(' mesh |',' nrdof_tot|',' nrdof_con|','    residual   |','   residual rate  ')
 7006  format(' mesh |',' nrdof_tot|',' nrdof_con|','    residual   |','   residual rate  |', &
              ' field error  |','rel field error|','   error rate ',/)
!
   do i=1,istep
      if (NEXACT.eq.0) then
         write(*,7003) i,nrdof_tot_mesh(i),nrdof_con_mesh(i),residual_mesh(i),rate_mesh(i)
 7003    format(2x,i2,'  | ',2(i8,' | '),es12.5,'  |',f7.2)
      else
         write(*,7004) i,nrdof_tot_mesh(i),nrdof_con_mesh(i),residual_mesh(i),rate_mesh(i), &
                       error_mesh(i),rel_error_mesh(i),rate_error_mesh(i)
 7004    format(2x,i2,'  | ',2(i8,' | '),es12.5,'  |',f7.2,'          ', &
                2(' | ',es12.5),'  |',f7.2)
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
         call update_gdof
         call update_Ddof
         nr_elem_to_refine = NRELES
         write(*,*) 'NRELES is', NRELES
!   ...anisotropic refinements
      case(IANISOTROPIC)
!     ...choose direction of refinements
         select case (ANISO_FLAG)
            case (IREFINE_X)
               write(*,*) 'Anisotropic refinement in x.'
            case (IREFINE_Y)
               write(*,*) 'Anisotropic refinement in y.'
            case (IREFINE_Z)
               write(*,*) 'Anisotropic refinement in z.'
            case (IREFINE_XY)
               write(*,*) 'Anisotropic refinement in xy.'
            case (IREFINE_XZ)
               write(*,*) 'Anisotropic refinement in xz.'
            case (IREFINE_YZ)
               write(*,*) 'Anisotropic refinement in yz.'
            case default
               write(*,*) 'refine_DPG: invalid ANISO_FLAG param. stop.'
               stop
         end select
         if (iso_ans.lt.1) then
            write(*,*) 'refine_DPG: invalid iso_ans param. stop.'
            stop
         endif
         do iii=1,iso_ans
            call setAnisoRef(ANISO_FLAG)
         enddo
         call close_mesh
         call update_gdof
         call update_Ddof
!
         nr_elem_to_refine = NRELES
         write(*,*) 'NRELES is', NRELES
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
   ! write(*,*) 'propagate_flag'
   if (IBCFLAG .eq. 3) then
      call propagate_flag(2,9)
      call propagate_flag(3,9)
   endif
   ! write(*,*) 'propagate_flag finished'
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
   use commonParam
!
   implicit none
!
   integer, intent(in)  :: Mdle
   integer, intent(out) :: Iref
!..geometry dof
   real(8) :: xnod(3,MAXbrickH), hmin
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
