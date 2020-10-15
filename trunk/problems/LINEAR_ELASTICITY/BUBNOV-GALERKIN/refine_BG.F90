!--------------------------------------------------------------------
!
!     routine name      - refine_BG
!
!--------------------------------------------------------------------
!
!     latest revision:  - Sep 15
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
!                       = 2 if adaptive refinements
!             Nreflag   = 1 for h-refinements
!                       = 2 for p-refinements
!                       = 3 for h- or p-refinement, depending
!                           upon the element size
!             Factor    - if element error \ge Factor*max_error
!                         then the element is refined
!     out:
!             Nstop     = 1 if refinement did not increase dof
!
!---------------------------------------------------------------------
!
      subroutine refine_BG(Irefine,Nreflag,Factor, Nstop)
!
      use control
      use data_structure3D
      use uhm2
      use common_prob_data
!
      implicit none
      integer, intent(in)  :: Irefine
      integer, intent(in)  :: Nreflag
      real*8,  intent(in)  :: Factor
      integer, intent(out) :: Nstop
!
      integer, parameter :: max_step = 20
      integer, dimension(max_step), save :: nrdofBG_mesh
      real*8,  dimension(max_step), save :: rate_mesh
      real*8,  dimension(max_step), save :: error_mesh
      real*8,  dimension(max_step), save :: rel_error_mesh
      real*8,  dimension(max_step), save :: rate_error_mesh

!
      integer, allocatable, dimension(:) :: elem_ref_flag
      integer, allocatable, dimension(:) :: list_elem
      integer, allocatable, dimension(:) :: list_ref_flags
!
      integer, save :: istep = 0
      integer, save :: irefineold = 0
      integer :: nr_elem_to_refine
      integer, dimension(NR_PHYSA) :: nflag
      real*8 :: err, rnorm, el_err,el_norm
      real*8 :: errorH,errorE,errorV,errorQ,rnormH,rnormE,rnormV,rnormQ
      integer :: i,iprint,ic,mdle,iel,kref,idec_href
!
!-----------------------------------------------------------------------
!                              INITIALIZE
!-----------------------------------------------------------------------
      iprint=0
      if (iprint.eq.1) write(*,*) 'refine_BG: DEBUGGING...'
!
!
      if (iprint.ge.1) write(*,*) 'refine_BG: ALLOCATING...'
      allocate (elem_ref_flag(NRELES),list_elem(NRELES), &
                list_ref_flags(NRELES))
!
!-----------------------------------------------------------------------
!                      INCREASE STEP IF NECESSARY
!-----------------------------------------------------------------------
!  ...irefineold=0 means no refinement was made in the previous step
!  ...if first call or if a refinement was made, increase step
      if ((istep.eq.0).or.(irefineold.ne.0)) istep=istep+1
      irefineold=Irefine
!
!-----------------------------------------------------------------------
!                            CALCULATE DOF
!-----------------------------------------------------------------------
!  ...field variables flag - used inside ndof_BG and for error calc.
      nflag(1:NR_PHYSA)=(/1/)
!  ...calculate dof
      nrdofBG_mesh(istep) = NRDOFSH + NRDOFSE + NRDOFSV + NRDOFSQ
!
!-----------------------------------------------------------------------
!                           CALCULATE ERROR
!-----------------------------------------------------------------------
!
      err = 0.d0; rnorm = 0.d0
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        if (NEXACT.ge.1) then
          ! call element_error(mdle,nflag,errorH,errorE,errorV,errorQ,  &
          !                               rnormH,rnormE,rnormV,rnormQ)
          ! err   = err   + (errorH+errorE+errorV+errorQ) !error norm
          ! rnorm = rnorm + (rnormH+rnormE+rnormV+rnormQ) !exact sol norm
          call custom_error(mdle, el_err,el_norm)
          err   = err + el_err !error norm
          rnorm = rnorm + el_norm !exact sol norm
        endif
      enddo
!
!-----------------------------------------------------------------------
!                UPDATE AND DISPLAY CONVERGENCE HISTORY
!-----------------------------------------------------------------------
!
      if (NEXACT.ge.1) then
         error_mesh(istep) = sqrt(err)
         rel_error_mesh(istep) = sqrt(err/rnorm)
      endif
!
!  ...compute decrease rate for the residual and error
      select case(istep)
      case(1)
        rate_mesh(istep) = 0.d0
        if (NEXACT.ge.1) rate_error_mesh(istep) = 0.d0
      case default
        rate_mesh(istep) =  &
         log(float(nrdofBG_mesh(istep-1))/float(nrdofBG_mesh(istep)))
        if (NEXACT.ge.1) then
          rate_error_mesh(istep) = &
          log(rel_error_mesh(istep-1)/rel_error_mesh(istep))/  &
          log(float(nrdofBG_mesh(istep-1))/float(nrdofBG_mesh(istep)))
        endif
      end select
!
!  ...print out the history of refinements
      write(*,*) 'HISTORY OF REFINEMENTS'
      if (NEXACT.eq.0) then
        write(*,*) 'Error: no exact solution'
        return
      endif
      if (NEXACT.gt.0) write(*,7006)
 7006   format(' mesh |',       &
               ' nrdofBG  |', &
               ' field error  |',  &
               'rel field error|', &
               ' error rate ')
      do i=1,istep
        write(*,7004) i,nrdofBG_mesh(i),error_mesh(i), &
                      rel_error_mesh(i),rate_error_mesh(i)
 7004   format(2x,i2,'  | ',i8,2(' | ',e12.5),'  | ',f7.2)
      enddo
!
!-----------------------------------------------------------------------
!                   CUSTOM DISPLAY OF RELEVANT ERRORS
!-----------------------------------------------------------------------
!                          Comment if desired
!-----------------------------------------------------------------------
!
!  ...CUSTOM - comment if you do not wish to have this
      call custom_display
!
!-----------------------------------------------------------------------
!                         REFINE AND UPDATE MESH
!-----------------------------------------------------------------------
!
!  ...use appropriate strategy depening on refinement type
!     NOTE: if Irefine=0 nothing is done, so only relevant info (the
!           code above this line) is displayed.
      select case(Irefine)
!  ...uniform refinements
      case(IUNIFORM)
        call global_href
        call update_gdof
        call update_ddof
        call verify_orient
        call verify_neig
        nr_elem_to_refine = NRELES
!
!  ...adaptive refinements
      case(IADAPTIVE)
        write(*,*) 'ERROR: No a posteriori error estimator'
!
      end select
!
!-----------------------------------------------------------------------
!                   DETERMINE IF MESH WAS REFINED
!-----------------------------------------------------------------------
!  ...determine if mesh was refined
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
      if (iprint.ge.1) write(*,*) 'refine_BG: DEALLOCATING...'
      deallocate (elem_ref_flag,list_elem,list_ref_flags)
!
      if (iprint.ge.1) write(*,*) 'refine_BG: EXITING...'
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!                         AUXILIARY SUBFUNCTIONS
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      contains
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!                  CUSTOM DISPLAY OF RELEVANT ERRORS
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        subroutine custom_display
!
        integer, parameter :: n_dispvar = 3
        integer, dimension(NR_PHYSA) :: nerrflag
        real*8, dimension(max_step,n_dispvar), save :: other_err
        real*8, dimension(max_step,n_dispvar), save :: other_rel_err
        real*8, dimension(n_dispvar) :: err_more,rnorm_more
!    ...compute other related errors which the user might want to know
        if (NEXACT.ge.1) then
          mdle=0
          err_more=0.d0; rnorm_more=0.d0
          do iel=1,NRELES
            call nelcon(mdle, mdle)
!           Displacement error
            nerrflag(1:NR_PHYSA)=(/1/)
            call element_error(mdle,nerrflag, errorH,errorE,errorV,errorQ, &
                                              rnormH,rnormE,rnormV,rnormQ)
            err_more(1)=err_more(1)+(errorH+errorE+errorV+errorQ)
            rnorm_more(1)=rnorm_more(1)+(rnormH+rnormE+rnormV+rnormQ)
!           Stress error
            nerrflag(1:NR_PHYSA)=(/0/)
            call element_error(mdle,nerrflag, errorH,errorE,errorV,errorQ, &
                                              rnormH,rnormE,rnormV,rnormQ)
            err_more(2)=err_more(2)+(errorH+errorE+errorV+errorQ)
            rnorm_more(2)=rnorm_more(2)+(rnormH+rnormE+rnormV+rnormQ)
!           Combined error=Stress+Displacement
            nerrflag(1:NR_PHYSA)=(/1/)
            call element_error(mdle,nerrflag, errorH,errorE,errorV,errorQ, &
                                              rnormH,rnormE,rnormV,rnormQ)
            err_more(3)=err_more(3)+(errorH+errorE+errorV+errorQ)
            rnorm_more(3)=rnorm_more(3)+(rnormH+rnormE+rnormV+rnormQ)
          enddo
          do i=1,3
            other_err(istep,i)=sqrt(err_more(i))
            other_rel_err(istep,i)=other_err(istep,i)/sqrt(rnorm_more(i))
          enddo
!      ...print extra information in new table
          write(*,*) 'RELEVANT ERRORS'
          write(*,7008)
 7008     format(' mesh |',       &
                 '    DispErr    |', &
                 '  DispRelErr   |', &
                 '   StrsErr     |', &
                 '  StrsRelErr   |', &
                 '   CombErr     |', &
                 '  CombRelErr   ')
          do i=1,istep
            write(*,7009) i,other_err(i,1),other_rel_err(i,1), &
                            other_err(i,2),other_rel_err(i,2), &
                            other_err(i,3),other_rel_err(i,3)
 7009       format(2x,i2,6('  | ',e12.5))
          enddo
        endif
!
        end subroutine custom_display
!
      end subroutine refine_BG