!> Purpose : solve a system of problem using UHM
#include "implicit_none.h"
subroutine uhm_solve_problem
  use assembly
  use control
  use data_structure3D
  use physics
  use uhm
  use locker

  implicit none
  integer :: &
       nodm(MAXNODM),ndofmH(MAXNODM),ndofmE(MAXNODM), &
       ndofmV(MAXNODM),ndofmQ(MAXNODM), idx(NRINDEX), &
       mvarH(MAXNODM),mvarE(MAXNODM),mvarV(MAXNODM),mvarQ(MAXNODM)
  integer :: &
       nrdofs(NR_PHYSA), max_dofs(NR_PHYSA)
  integer :: &
       iprint, istat,  &
       i, j, k, mm, iel, ino, iphy, irhs, ivar, jvar, lda, mdle, nod, &
       nvarH, nvarE, nvarV, nvarQ, &
       ndofH, ndofE, ndofV, ndofQ, &
       nrdofm, nrdofc, nrnodm, max_dofm, max_dofc
  real*8 :: t_base, t_after
  VTYPE,pointer :: ZVOID0,ZVOID1

  iprint = 0

  ! Step 1 : Initialization
  !~~~~~~~~~~~~~~~~~~~~~~~~

  ! 1.1 set number of rhs in assembly
  NR_RHS = UHM_N_RHS

  ! 1.2 allocate array for the assembly dirichlet 
  call assembly_begin

  ! 1.3 allocate pointers to UHM same number of NRELES
  call uhm_alloc(NRELES)

  ! 1.4 set ISYM_FLAG = 2 as unsymmetric, UHM does not use compressed form
  ISYM_FLAG      = 2

  ! 1.5 lock datastructure
  call data_lock

  ! Step 2 : Symbolic interface UHM to HP3D
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ! 2.1 initialize the max dof
  MAXDOFM = 0; MAXDOFC =0; MAXDOFS = 0;

  ! 2.2 loop over all elements

  write(*,*) 'uhm_solve_problem: BEGIN OF FEEDING CONNECTIVITY FOR ELTS '
  call uhm_timer        (t_base)
  mdle = 0
  do iel=1, NRELES
     call nelcon_from_locker(iel, mdle)

     ! 2.2.1 celem catch modified degree of freedoms excluding all dirichlet BC
     call celem( &
          mdle, 1, &
          nrdofs,nrdofm,nrdofc, &
          nodm,ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm, &
          ZVOID0,ZVOID1)

     ! 2.2.2 update max dof for each physical attributes
     do i=1, NR_PHYSA
        MAXDOFS(i) = max(MAXDOFS(i), nrdofs(i))
     end do

     MAXDOFM = max(MAXDOFM, nrdofm)
     MAXDOFC = max(MAXDOFC, nrdofc)

     ! 2.2.3 if element has non-zero dofs interface to UHM
     if (nrdofc.gt.0) then

        ! 2.2.4 create UHM element 
        call uhm_add_element(UHM_MESH, UHM_ELTS(iel)%id)

        ! counter
        k  = 0
        mm = 0

        ! 2.2.5 loop over physics :: fill UHM_ELTS
        ! id(1) represent node number
        ! id(2) physics number, different physics is considered independent node
        ! n_dof should handle all # of eqns and # of variables
        !       UHM does is not necessary to know all the details
        !
        do ino=nrnodm, 1, -1
           if (ndofmH(ino).gt.0) then
              k = k + 1
              mm = mm + ndofmH(ino)
              UHM_ELTS(iel)%nods(1,k) = nodm(ino)
              UHM_ELTS(iel)%nods(2,k) = UHM_PHYSICS_CONTIN
              UHM_ELTS(iel)%n_dof(k)  = ndofmH(ino)
              UHM_ELTS(iel)%p(k)      = NODES(nodm(ino))%order
           end if
           
           if (ndofmE(ino).gt.0) then
              k = k + 1
              mm = mm + ndofmE(ino)
              UHM_ELTS(iel)%nods(1,k) = nodm(ino)
              UHM_ELTS(iel)%nods(2,k) = UHM_PHYSICS_TANGEN
              UHM_ELTS(iel)%n_dof(k)  = ndofmE(ino)
              UHM_ELTS(iel)%p(k)      = NODES(nodm(ino))%order
           end if
           
           if ( (ndofmV(ino).gt.0).or.(ndofmQ(ino).gt.0) ) then
              write(*,*) 'uhm_solve_problem: NOT SUPPORT'
              stop 1
           end if
        end do


        ! 2.2.6 store the effective nodes
        UHM_ELTS(iel)%n_nods = k
        UHM_ELTS(iel)%mm     = mm

        ! print statement
        if (iprint.eq.1) then
           write(*,*) 'mdle = ',mdle
           write(*,*) 'nrdofs(1:',NR_PHYSA,')=',nrdofs(1:NR_PHYSA)
           write(*,*) 'nrdofm=',nrdofm
           write(*,*) 'nrdofc=',nrdofc
           write(*,*) 'nodm(1:',nrnodm,')=', nodm(1:nrnodm)
           write(*,*) 'ndofmH(1:',nrnodm,')=', ndofmH(1:nrnodm)
           write(*,*) 'ndofmE(1:',nrnodm,')=', ndofmE(1:nrnodm)

           call uhm_elt_disp(NSTD_OUT, UHM_ELTS(iel))
           call pause
           cycle
        end if

        ! 2.2.7 create nodes and add the nodes to associate element in UHM
        do ino=1, UHM_ELTS(iel)%n_nods
           call uhm_add_node( &
                UHM_MESH, &
                UHM_ELTS(iel)%nods(1,ino), &
                UHM_ELTS(iel)%nods(2,ino), &
                UHM_ELTS(iel)%n_dof(ino), &
                UHM_ELTS(iel)%p(ino), nod )
           call uhm_element_add_node( UHM_ELTS(iel)%id, nod )
        end do

     end if
  end do
  call uhm_timer        (t_after)
  write(*,*) 'uhm_solve_problem: END OF FEEDING CONNECTIVITY FOR ELTS '
  write(*,*) 'uhm_solve_problem: TIME ', t_after - t_base, ' SEC'

  ! Step 3 : Allocate memory for assembly module with max dofs
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  allocate(UHM_ZSTIFF(MAXDOFM*MAXDOFM), UHM_ZXLOAD(MAXDOFM*NR_RHS), stat=istat)
  if (istat.ne.SUCCESS) then
     call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
  end if

  call assembly_alloc

  ! Step 4 : Analyze the problem
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  call uhm_timer        (t_base)
  select case (UHM_METHOD)
  case (UHM_INTERF_PARDISO);       call uhm_analyze_pardiso_1
  case (UHM_INTERF_WSMP);          call uhm_analyze_wsmp_1
  case (UHM_INTERF_MUMPS);         call uhm_analyze_mumps_1
  case default;                    call uhm_analyze_uhm
  end select
  call uhm_timer        (t_after)
  write(*,*) 'uhm_solve_problem: ANALYZE TIME  ', t_after - t_base, ' SEC'
  
  ! Step 5 : Matrix interface UHM to HP3D
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ! 5.1 loop over all elements
  write(*,*) 'uhm_solve_problem: BEGIN OF FEEDING MATRIX FROM CELEM FOR ELTS '
  call uhm_timer        (t_base)
  mdle = 0

  !!!$OMP PARALLEL DO PRIVATE(iel, mdle) 
  do iel=1, NRELES
     call nelcon_from_locker(iel, mdle)
     
     ! 5.1.1 celem to catch modified matrix
     call celem( &
          mdle, 2, &
          nrdofs,nrdofm,nrdofc, &
          nodm,ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm, &
          UHM_ZXLOAD,UHM_ZSTIFF)

     ! 5.1.2 if element has non-zero dofs interface to UHM
     if (nrdofc.gt.0) then

        call uhm_copy_in( &
             UHM_MESH, &
             UHM_ELTS(iel), &
             UHM_DATATYPE, &
             UHM_ELTS(iel)%mm, &
             UHM_ELTS(iel)%mm, &
             UHM_PHYSICS_MULTI, &
             UHM_ELTS(iel)%nods, &
             UHM_LHS, &
             UHM_ZSTIFF)

        call uhm_copy_in( &
             UHM_MESH, &
             UHM_ELTS(iel), &
             UHM_DATATYPE, &
             UHM_ELTS(iel)%mm, &
             UHM_N_RHS, &
             UHM_PHYSICS_MULTI, &
             UHM_ELTS(iel)%nods, &
             UHM_RHS, &
             UHM_ZXLOAD)

     end if
  end do
  !!!$OMP END PARALLEL DO

  call uhm_timer        (t_after)
  write(*,*) 'uhm_solve_problem: TIME ', t_after - t_base, ' SEC'
  write(*,*) 'uhm_solve_problem: END OF FEEDING MATRIX FROM CELEM FOR ELTS '

  ! Step 6 : Decompose the matrix
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  select case (UHM_METHOD)
  case (UHM_INTERF_PARDISO);       call uhm_analyze_pardiso_2
  case (UHM_INTERF_WSMP);          call uhm_analyze_wsmp_2
  case (UHM_INTERF_MUMPS);         call uhm_analyze_mumps_2
  end select

  call uhm_timer        (t_base)
  select case (UHM_METHOD)
  case (UHM_INTERF_PARDISO);       call uhm_pardiso_decompose (UHM_PARDISO)
  case (UHM_INTERF_WSMP);          call uhm_wsmp_decompose    (UHM_WSMP)
  case (UHM_INTERF_MUMPS);         call uhm_mumps_decompose   (UHM_MUMPS)
  case default;                    call uhm_decompose_uhm
  end select
  call uhm_timer        (t_after)
  write(*,*) 'uhm_solve_problem: DECOMPOSE TIME  ', t_after - t_base, ' SEC'


  ! Step 7 : Solve problem
  !~~~~~~~~~~~~~~~~~~~~~~~
  call uhm_timer        (t_base)
  select case (UHM_METHOD)
  case (UHM_INTERF_PARDISO);      
     call uhm_pardiso_solve(UHM_PARDISO)
     call uhm_pardiso_export_matrix_uhm(UHM_PARDISO, UHM_MESH)
  case (UHM_INTERF_WSMP);          
     call uhm_wsmp_solve(UHM_WSMP)
     call uhm_wsmp_refine(UHM_WSMP)
     call uhm_wsmp_export_matrix_uhm(UHM_WSMP, UHM_MESH)
  case (UHM_INTERF_MUMPS);      
     call uhm_mumps_solve(UHM_MUMPS)
     call uhm_mumps_export_matrix_uhm(UHM_MUMPS, UHM_MESH)
  case default;                
     call uhm_solve_uhm
  end select
  call uhm_timer        (t_after)
  write(*,*) 'uhm_solve_problem: SOLVE TIME  ', t_after - t_base, ' SEC'


  ! Step 8 : Check residual
  !~~~~~~~~~~~~~~~~~~~~~~~~
  call uhm_timer        (t_base)
  select case (UHM_METHOD)
  case (UHM_INTERF_PARDISO, UHM_INTERF_WSMP, UHM_INTERF_MUMPS);          
  case default;                
     call uhm_check_residual_uhm
  end select
  call uhm_timer        (t_after)
  write(*,*) 'uhm_solve_problem: CHECK TIME  ', t_after - t_base, ' SEC'

  ! Step 9 : Solution out
  !~~~~~~~~~~~~~~~~~~~~~~~

  ! 9.1 loop over all elements
  write(*,*) 'uhm_solve_problem: BEGIN OF STORING SOLUTIONS FOR ELTS'
  call uhm_timer        (t_base)
  do iel=1, NRELES

     ! counters for H1, H(curl), H(div), L2 components supported by node "nod"
     mvarH(1:MAXNODM)=0 ; mvarE(1:MAXNODM)=0 ; mvarV(1:MAXNODM)=0 ; mvarQ(1:MAXNODM)=0

     ! 9.1.2 if element has non-zero dofs interface to UHM
     if (UHM_ELTS(iel)%mm.gt.0) then

        lda  = UHM_ELTS(iel)%mm

        ! 9.1.3 catch solution from UHM
        UHM_ZXLOAD = ZERO
        call uhm_copy_out( &
             UHM_MESH, &
             UHM_ELTS(iel), &
             UHM_DATATYPE, &
             UHM_ELTS(iel)%mm, &
             UHM_N_RHS, &
             UHM_PHYSICS_MULTI, &
             UHM_ELTS(iel)%nods, &
             UHM_RHS, &
             UHM_ZXLOAD)

        ! 9.1.4 
        do irhs=1, UHM_N_RHS
           i = 0
           
           ! 9.1.4.1 loop over nodes of modified element
           do ino=1, UHM_ELTS(iel)%n_nods
              nod  = UHM_ELTS(iel)%nods(1,ino)
              iphy = UHM_ELTS(iel)%nods(2,ino)
              
              call find_ndof(nod, ndofH, ndofE, ndofV, ndofQ)
              call get_index(nod, idx)

              if (iprint.eq.1) then
                 write(*,*) 'uhm_solve_problem: nod, ndofH, nvarE, ndofV, ndofQ ', &
                      nod, ndofH, ndofE, ndofV, ndofQ
              endif

              select case(iphy)
              case(UHM_PHYSICS_CONTIN) 

                 ! loop over H1 dofs associated to node
                 do j=1, ndofH
 
                    ! number of H1 components stored so far
                    ivar=mvarH(ino)

                    ! loop over H1, H(curl), H(div), L2 components 
                    do k=1,NRINDEX
                       select case(idx(k))

                       ! H1 component with Dirchlet BC
                       case(1)
                          ! increment index, do no store     
                          ivar=ivar+1

                       ! free H1 component
                       case(2) 
                          ! increment index     
                          ivar=ivar+1 

                          ! store
                          i=i+1
                          NODES(nod)%zdofH(ivar,j) = UHM_ZXLOAD(lda*(irhs-1)+i)
                       end select
                    end do

                 end do

                 !  update number of stored H1 components
                 mvarH(ino)=ivar

              case(UHM_PHYSICS_TANGEN)
                 do j=1, ndofE
                    do ivar=1, nvarE
                       i = i + 1
                       NODES(nod)%zdofE(ivar,j) = UHM_ZXLOAD(lda*(irhs-1)+i)
                    end do
                 end do
              case default
                 write(*,*) 'uhm_solve_problem: NOT SUPPORT (iphy = ', iphy, ')'
                 stop 1
              end select
           end do
        end do
     end if
  end do
  call uhm_timer        (t_after)
  write(*,*) 'uhm_solve_problem: END OF STORING SOLUTIONS FOR ELTS'
  write(*,*) 'uhm_solve_problem: SOLUTION OUT TIME  ', t_after - t_base, ' SEC'

  ! Step 10 : Finalize
  !~~~~~~~~~~~~~~~~~~~~~~~~
  ISYM_FLAG      = UHM_ISYM_FLAG

  call uhm_timer        (t_base)
  deallocate(UHM_ZSTIFF, stat=istat)
  deallocate(UHM_ZXLOAD, stat=istat)

  call uhm_dealloc
  call assembly_dealloc
  call assembly_end

  call data_unlock
  call uhm_timer        (t_after)
  write(*,*) 'uhm_solve_problem: FINALIZATION TIME  ', t_after - t_base, ' SEC'

end subroutine uhm_solve_problem
