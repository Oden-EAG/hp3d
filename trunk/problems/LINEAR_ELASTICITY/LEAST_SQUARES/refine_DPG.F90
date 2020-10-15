!--------------------------------------------------------------------
!
!     routine name      - refine_DPG
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
      subroutine refine_DPG(Irefine,Nreflag,Factor, Nstop)
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
      integer, dimension(max_step), save :: nrdofDPG_mesh
      integer, dimension(max_step), save :: nrdofField_mesh      
      real*8, dimension(max_step), save :: residual_mesh
      real*8, dimension(max_step), save :: rate_mesh
      real*8, dimension(max_step), save :: error_mesh
      real*8, dimension(max_step), save :: rel_error_mesh
      real*8, dimension(max_step), save :: rate_error_mesh

!
      real*8, allocatable, dimension(:) :: elem_resid
      integer, allocatable, dimension(:) :: elem_ref_flag
      integer, allocatable, dimension(:) :: list_elem
      integer, allocatable, dimension(:) :: list_ref_flags
!
      integer, save :: istep = 0
      integer, save :: irefineold = 0
      integer :: nrdofField, nrdofDPG
      integer :: nr_elem_to_refine
      integer, dimension(NR_PHYSA) :: nflag
      real*8 :: residual, elem_resid_max, err, rnorm
      real*8 :: errorH,errorE,errorV,errorQ,rnormH,rnormE,rnormV,rnormQ
      integer :: i,iprint,ic,mdle,iel,kref,idec_href
!
!-----------------------------------------------------------------------
!                              INITIALIZE
!-----------------------------------------------------------------------
      iprint=0
      if (iprint.eq.1) write(*,*) 'refine_DPG: DEBUGGING...'
!
!
      if (iprint.ge.1) write(*,*) 'refine_DPG: ALLOCATING...'
      allocate (elem_resid(NRELES),elem_ref_flag(NRELES),  &
                list_elem(NRELES),list_ref_flags(NRELES))
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
!  ...field variables flag - used inside ndof_DPG and for error calc.
      nflag(1:NR_PHYSA)=(/1,1/)
!  ...calculate dof
      call ndof_DPG !calculates nrdofDPG and nrdofField
      nrdofDPG_mesh(istep)   = nrdofDPG
      nrdofField_mesh(istep) = nrdofField
!
!-----------------------------------------------------------------------
!                         CALCULATE RESIDUAL
!-----------------------------------------------------------------------
!
      residual = 0.d0; elem_resid_max = 0.d0
      err = 0.d0; rnorm = 0.d0
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        call elem_residual(mdle, elem_resid(iel),elem_ref_flag(iel))
        if (NEXACT.ge.1) then
          call element_error(mdle,nflag,errorH,errorE,errorV,errorQ,  &
                                        rnormH,rnormE,rnormV,rnormQ)
          err   = err   + (errorH+errorE+errorV+errorQ) !error norm
          rnorm = rnorm + (rnormH+rnormE+rnormV+rnormQ) !exact sol norm
        endif
        if (iprint.eq.1) write(*,7010) iel,mdle,elem_resid(iel)
 7010     format('refine_DPG: iel,mdle = ',i5,i6,  &
                  ' residual^2 = ',e12.5)
        elem_resid_max = max(elem_resid_max,elem_resid(iel))
        residual = residual + elem_resid(iel)
      enddo
!
!-----------------------------------------------------------------------
!                UPDATE AND DISPLAY CONVERGENCE HISTORY
!-----------------------------------------------------------------------
!
      residual_mesh(istep) = sqrt(residual)
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
         log(residual_mesh(istep-1)/residual_mesh(istep))/  &
         log(float(nrdofDPG_mesh(istep-1))/float(nrdofDPG_mesh(istep)))
        if (NEXACT.ge.1) then
          rate_error_mesh(istep) = &
          log(rel_error_mesh(istep-1)/rel_error_mesh(istep))/  &
          log(float(nrdofDPG_mesh(istep-1))/float(nrdofDPG_mesh(istep)))
        endif
      end select
!
!  ...print out the history of refinements
      write(*,*) 'HISTORY OF REFINEMENTS'
      if (NEXACT.eq.0) write(*,7005)
      if (NEXACT.gt.0) write(*,7006)
 7005   format(' mesh |',       &
               ' nrdofDPG |', &
               'nrdofField|', &
               '    residual   |', &
               '  rate   ')
 7006   format(' mesh |',       &
               ' nrdofDPG |', &
               'nrdofField|', &
               '    residual   |', &
               '  rate   |',   &
               ' field error  |',  &
               'rel field error|', &
               ' error rate ')
      do i=1,istep
        if (NEXACT.eq.0) then
          write(*,7003) i,nrdofDPG_mesh(i),nrdofField_mesh(i), &
                         residual_mesh(i),rate_mesh(i)
 7003     format(2x,i2,'  | ',2(i8,' | '),e12.5,'  | ',f7.2)
        else
          write(*,7004) i,nrdofDPG_mesh(i),nrdofField_mesh(i), &
                         residual_mesh(i),rate_mesh(i),error_mesh(i), &
                          rel_error_mesh(i),rate_error_mesh(i)
 7004     format(2x,i2,'  | ',2(i8,' | '),e12.5,'  | ',f7.2, &
                  2(' | ',e12.5),'  | ',f7.2)
        endif
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
!    ...use the greedy strategy to determine which elements to refine
        ic=0
        mdle=0
        do iel=1,NRELES
          call nelcon(mdle, mdle)
          if (elem_resid(iel).ge.Factor*elem_resid_max) then
            ic=ic+1
            list_elem(ic) = mdle
            list_ref_flags(ic) = elem_ref_flag(iel)
            if (iprint.eq.1) then
              write(*,7037) ic,mdle
 7037         format('refine_DPG: ADDING TO THE LIST: ic,mdle = ',2i6)
            endif
          endif
        enddo
        nr_elem_to_refine = ic
! 
        if (nr_elem_to_refine.gt.0) then
!      ...refine the elements from the list
          do iel=1,nr_elem_to_refine
            mdle = list_elem(iel)
            kref = list_ref_flags(iel)
!    
!        ...restrictions on h-refinements
            idec_href=1
    !!!        call find_gen(mdle,ngen)
    !!!        select case(NODES(mdle)%type)
    !!!        case('mdlt')
    !!!          if (ngen(1).eq.MAXGENT) idec_href=0
    !!!        case('mdlq')
    !!!         if (ngen(1).eq.MAXGENQ(1)) idec_href=0
    !!!         if (ngen(2).eq.MAXGENQ(2)) idec_href=0
    !!!        end select
    !
    !  .....p-refinements limited by the maximum order supported
    !       by the code
    !        nord = NODES(mdle)%order
    !!!        select case(NODES(mdle)%type)
    !!!        case('mdlt')
    !!!          nord_new = min(nord+1,MAXP-NORD_ADD-1)
    !!!        case('mdlq')
    !!!          call decode(nord, nordh,nordv)
    !!!         nordh_new = nordh; nordv_new = nordv
    !!!          if (kref(1).eq.1) then
    !!!            nordh_new = min(nordh+1,MAXP-NORD_ADD-1)
    !!!          endif
    !!!          if (kref(2).eq.1) then
    !!!            nordv_new = min(nordv+1,MAXP-NORD_ADD-1)
    !!!          endif
    !!!          nord_new = nordh_new*10+nordv_new
    !!!        end select
    !!!        if (nord_new.eq.nord) then
    !!!          idec_pref=0
    !!!        else
    !!!          idec_pref=1
    !!!        endif
  !
            select case(Nreflag)
            case(1)
              if (idec_href.eq.1) then
                call refine(mdle,kref)
                write(*,*) 'mdle,kref = ',mdle,kref
              endif
            case(2)
              write(*,*) 'p UNFINISHED '; stop1
            case(3)
             write(*,*) 'hp UNFINISHED '; stop1
            end select
          enddo
          !
!      ...close the mesh
          call close
!      ...update geometry and Dirichlet flux dof after the refinements
          call update_gdof
          call update_Ddof
        endif
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
      if (iprint.ge.1) write(*,*) 'refine_DPG: DEALLOCATING...'
      deallocate (elem_resid,elem_ref_flag,list_elem,list_ref_flags)
!
      if (iprint.ge.1) write(*,*) 'refine_DPG: EXITING...'
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
!                             CALCULATE DOF 
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!    ...defines nrdofDPG and nrdofField for current mesh
        subroutine ndof_DPG
!
        integer :: iphys,nH,nE,nV,nQ,ndofH,ndofE,ndofV,ndofQ
!
!    ...nflag=field variables flag comes from father routine refine_DPG
!    ...calculate number of components in each energy space
        nH=0;nE=0;nV=0;nQ=0
        do iphys=1,NR_PHYSA
          select case(DTYPE(iphys))
          case('contin'); nH = nH + NR_COMP(iphys)
          case('tangen'); nE = nE + NR_COMP(iphys)
          case('normal'); nV = nV + NR_COMP(iphys)
          case('discon'); nQ = nQ + NR_COMP(iphys)
          end select
        enddo
!    ...calculate nrdofField
        nrdofField = 0
        do iphys=1,NR_PHYSA
!         skip the trace variables
          if (nflag(iphys).eq.0) cycle   
          select case(DTYPE(iphys))
          case('contin'); nrdofField=nrdofField+(NRDOFSH/nH)*NR_COMP(iphys)
          case('tangen'); nrdofField=nrdofField+(NRDOFSE/nE)*NR_COMP(iphys)
          case('normal'); nrdofField=nrdofField+(NRDOFSV/nV)*NR_COMP(iphys)
          case('discon'); nrdofField=nrdofField+(NRDOFSQ/nQ)*NR_COMP(iphys)
          end select
        enddo
!    ...calculation of nrdofDPG - first add all possible dof
        nrdofDPG = NRDOFSH + NRDOFSE + NRDOFSV + NRDOFSQ
!    ...then correct by subtracting bubbles
        mdle=0
        do iel=1,NRELES
          call nelcon(mdle, mdle)
!    .....adjust number of nrdofDPG
          call find_ndof(mdle, ndofH,ndofE,ndofV,ndofQ) !per component
          do iphys=1,NR_PHYSA
!        ...skip if field variable, otherwise remove middle dof and 
!           leave trace dof only
            if (nflag(iphys).eq.1) cycle
            select case(DTYPE(iphys))
            case('contin'); nrdofDPG=nrdofDPG-ndofH*NR_COMP(iphys)
            case('tangen'); nrdofDPG=nrdofDPG-ndofE*NR_COMP(iphys)
            case('normal'); nrdofDPG=nrdofDPG-ndofV*NR_COMP(iphys)
            case('discon'); nrdofDPG=nrdofDPG-ndofQ*NR_COMP(iphys)
            end select
          enddo
        enddo
!
        end subroutine ndof_DPG
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
            nerrflag(1:NR_PHYSA)=(/1,0/)
            call element_error(mdle,nerrflag, errorH,errorE,errorV,errorQ, &
                                              rnormH,rnormE,rnormV,rnormQ)
            err_more(1)=err_more(1)+(errorH+errorE+errorV+errorQ)
            rnorm_more(1)=rnorm_more(1)+(rnormH+rnormE+rnormV+rnormQ)
!           Stress error
            nerrflag(1:NR_PHYSA)=(/0,1/)
            call element_error(mdle,nerrflag, errorH,errorE,errorV,errorQ, &
                                              rnormH,rnormE,rnormV,rnormQ)
            err_more(2)=err_more(2)+(errorH+errorE+errorV+errorQ)
            rnorm_more(2)=rnorm_more(2)+(rnormH+rnormE+rnormV+rnormQ)
!           Combined error=Stress+Displacement
            nerrflag(1:NR_PHYSA)=(/1,1/)
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
      end subroutine refine_DPG