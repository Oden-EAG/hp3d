












subroutine Finegrid_padap(nr_elem_ref,mdle_ref,xnod_ref,flag_pref)

    use common_prob_data
    use control
    use data_structure3D
    use environment     , only: QUIET_MODE
    use assembly_sc     , only: NRDOF_CON, NRDOF_TOT
    use parametersDPG
    use par_mesh        , only: DISTRIBUTED, HOST_MESH
    use mpi_param       , only: ROOT, RANK,NUM_PROCS
    use MPI             , only: MPI_COMM_WORLD, MPI_SUM,MPI_COMM_WORLD, &
                                MPI_REAL8, MPI_INTEGER, MPI_IN_PLACE, MPI_MAX

    implicit none
    integer, parameter :: Irefine = 2
    integer, parameter :: adap_strat = 1
    ! integer, parameter :: Factor = 0.75
    ! logical :: Ires = .true.
    integer, intent(in) :: nr_elem_ref
    integer, intent(out) :: flag_pref(nr_elem_ref)

    integer, dimension(nr_elem_ref),    intent(in)  :: mdle_ref
    real(8), dimension(nr_elem_ref,3,MAXbrickH),    intent(out)  ::  xnod_ref
    
    real(8) :: dummy_xnod(3,MAXbrickH)

    integer :: nord_new,is,nord,nordx,nordy,nordz,naux,pord
    integer :: i,ic,mdle,iel,kref,subd,count,ierr
    integer :: Nreflag
    integer :: Nstop
    real(8) :: MPI_Wtime,start_time,end_time
    !..element type
    character(len=4) :: etype


    select case(Irefine)

    case(IADAPTIVE)
        call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
        if (RANK .eq. ROOT)  write(*,*) 'Starting hp refining the mesh to compute fine mesh solution'
        if (RANK .eq. ROOT)  write(*,*) 'NRELES = ', NRELES        

        if(adap_strat .eq. 0) then
            do iel = 1,nr_elem_ref
             mdle = mdle_ref(iel)
             etype = NODES(mdle)%type
             nord = NODES(mdle)%order
             select case(etype)
             case('mdlb')
                kref = 111
                nord_new = nord + 111
                
                call decode(nord_new,naux,nordz)
                call decode(naux,nordx,nordy)
                pord = MAX(nordx,nordy,nordz)
                if (pord .gt. MAXP) then
                    write(*,1001) 'local pref: mdle,p,MAXP = ',mdle,pord,MAXP,'. stop.'
                    stop
                    1001 format(A,I7,I3,I3,A)
                endif

             case default
                write(*,*) 'refine_DPG: READING UNEXPECTED ELEMENT TYPE: ',etype
                call pause
             end select
            !p refine and then h-refine
                call nodmod(mdle,nord_new)
                call refine(mdle,kref)
               !  call break(mdle,kref)
            enddo
        elseif(adap_strat .eq. 1) then
            ! p refinement of marked elements
            do iel = 1,nr_elem_ref
                mdle = mdle_ref(iel)
                etype = NODES(mdle)%type
                nord = NODES(mdle)%order
                select case(etype)
                   case('mdlb')
                    !   kref = 111 ! iso
                      nord_new = nord + 111
                    !   write(*,*) "the elem is = ", mdle,nord_new
                      call decode(nord_new,naux,nordz)
                      call decode(naux,nordx,nordy)
                      pord = MAX(nordx,nordy,nordz)
                     !  write(*,*) mdle,pord,MAXP
                      if (pord .gt. MAXP) then
                          write(*,1002) 'local pref: mdle,p,MAXP = ',mdle,pord,MAXP,'. stop.'
                          go to 700
                          1002 format(A,I7,I3,I3,A)
                      endif 
                    

                   case default
                      write(*,*) 'refine_DPG: READING UNEXPECTED ELEMENT TYPE: ',etype
                      call pause
                end select
                !p refine only if pord is less than MAXP : thats why we have a goto above in if statement
                call nodmod(mdle,nord_new)  
                flag_pref(iel) = 1
                700 continue !not performing p-ref to any element which are already at the highest polynomial order.

            enddo

            call MPI_BARRIER (MPI_COMM_WORLD, ierr)
            ! end_time = MPI_Wtime()
            ! if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2020) end_time-start_time
            ! call MPI_BARRIER (MPI_COMM_WORLD, ierr)
            ! start_time = MPI_Wtime()
            call enforce_min_rule
            call MPI_BARRIER (MPI_COMM_WORLD, ierr)
            ! end_time = MPI_Wtime()
            ! if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2025) end_time-start_time
            !..raise order of approximation on non-middle nodes by enforcing minimum rule
            
            call par_verify
            call update_gdof
            call update_Ddof

            ! write(*,*)  "here 2 = ",MAXNODS,NRNODS,NPNODS
            call MPI_BARRIER (MPI_COMM_WORLD, ierr)
            ! start_time = MPI_Wtime()
            
            !h-refinement of marked elements

            do iel = 1,nr_elem_ref
                mdle = mdle_ref(iel)
                etype = NODES(mdle)%type
                nord = NODES(mdle)%order
                select case(etype)
                   case('mdlb')
                      kref = 111 ! iso
                    !   nord_new = nord + 111
                    !   write(*,*) "the elem is = ", mdle,nord_new
                    !   call decode(nord_new,naux,nordz)
                    !   call decode(naux,nordx,nordy)
                    !   pord = MAX(nordx,nordy,nordz)  
                   case default
                      write(*,*) 'refine_DPG: READING UNEXPECTED ELEMENT TYPE: ',etype
                      call pause
                end select
                
                call refine(mdle,kref)
            enddo

            ! call enforce_min_rule
            call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
            if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2020) end_time-start_time
            call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
            call close_mesh
            ! call enforce_min_rule
            call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
            if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2025) end_time-start_time
            !..raise order of approximation on non-middle nodes by enforcing minimum rule
            
            call par_verify
            call update_gdof
            call update_Ddof
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        endif


        ! if (RANK .eq. ROOT) write(*,*) 'NRELES = ', NRELES
        ! if (RANK .eq. ROOT) write(*,*) 'Finished hp refining selected elements'


    case default; Nstop = 1
    end select


    ! do iel = 1,nr_elem_ref
    !     mdle = mdle_ref(iel)
    !     etype = NODES(mdle)%type
    !     nord = NODES(mdle)%order
    !     call nodcor(mdle,dummy_xnod)
    !     xnod_ref(iel,1:3,1:MAXbrickH) = dummy_xnod(1:3,1:MAXbrickH)
    ! enddo

    ! select case(Irefine)

    ! case(IADAPTIVE)
    !     call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
    !     if (RANK .eq. ROOT)  write(*,*) 'Starting hp refining the mesh to compute fine mesh solution'
    !     if (RANK .eq. ROOT)  write(*,*) 'NRELES = ', NRELES        

    !     if(adap_strat .eq. 0) then
    !         do iel = 1,nr_elem_ref
    !          mdle = mdle_ref(iel)
    !          etype = NODES(mdle)%type
    !          nord = NODES(mdle)%order
    !          select case(etype)
    !          case('mdlb')
    !             kref = 111

    !          case default
    !             write(*,*) 'refine_DPG: READING UNEXPECTED ELEMENT TYPE: ',etype
    !             call pause
    !          end select
    !          ! h-refine
    !          call refine(mdle,kref)
            
    !         enddo
    !     elseif(adap_strat .eq. 1) then
    !         do iel = 1,nr_elem_ref
    !             mdle = mdle_ref(iel)
    !             etype = NODES(mdle)%type
    !             nord = NODES(mdle)%order
    !             select case(etype)
    !                case('mdlb')
    !                   kref = 111 ! iso
    !                case default
    !                   write(*,*) 'refine_DPG: READING UNEXPECTED ELEMENT TYPE: ',etype
    !                   call pause
    !             end select
                
    !             !h refine
    !             call refine(mdle,kref)
    !         enddo
    !     endif
    !   !   call enforce_min_rule
    !     call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
    !     if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2020) end_time-start_time
    !     call MPI_BARRIER (MPI_COMM_WORLD, ierr); start_time = MPI_Wtime()
    !     call close_mesh
    !     call MPI_BARRIER (MPI_COMM_WORLD, ierr); end_time = MPI_Wtime()
    !     if ((.not. QUIET_MODE) .and. (RANK .eq. ROOT)) write(*,2025) end_time-start_time
        
    !     call par_verify
    !     call update_gdof
    !     call update_Ddof

    !     ! if (RANK .eq. ROOT) write(*,*) 'NRELES = ', NRELES
    !     if (RANK .eq. ROOT) write(*,*) 'Finished h refining selected elements'


    ! case default; Nstop = 1
    ! end select






    2020 format(' fine mesh refinement : ',f12.5,'  seconds')
    2025 format(' close mesh : ',f12.5,'  seconds')
    2030 format(A,I8,', ',I9)



end subroutine Finegrid_padap