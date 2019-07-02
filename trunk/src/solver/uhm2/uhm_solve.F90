!> Purpose : solve a system of problem using UHM
#include "implicit_none.h"
subroutine uhm_solve
  use environment
  use assembly
  use control
  use data_structure3D
  use physics
  use solvermod
  use uhm2

  implicit none
  integer :: &
       nodm(MAXNODM),ndofmH(MAXNODM),ndofmE(MAXNODM), &
       ndofmV(MAXNODM),ndofmQ(MAXNODM),               &      
       mvarH(MAXNODM),mvarE(MAXNODM),mvarV(MAXNODM),mvarQ(MAXNODM)
  integer :: &
       nrdofs(NR_PHYSA), max_dofs(NR_PHYSA), idx(NRINDEX)
  integer :: &
       ime, istat, ndom, iper, &
       i, j, k, mm, iel, ino, iph, irhs, ivar,  jvar, mdle, nod, &
       nvarH, nvarE, nvarV, nvarQ, &
       ndofH, ndofE, ndofV, ndofQ, &
       nrdofm, nrdofc, nrnodm, max_dofm, max_dofc, n_dofs
  real*8 :: t, residual, dper
  logical iflag

  ! Step 1 : Initial Setup
  !~~~~~~~~~~~~~~~~~~~~~~~

  ! 1.1 set number of rhs in assembly
  call uhm_solver_get_rhs(UHM_SOLVER_PTR, NR_RHS)

  ! 1.2 allocate array for the assembly dirichlet 
  call assembly_begin

  ! 1.3 initialize a temporary work array for UHM communication
  UHM_N_ELTS = 0
  call uhm_alloc(NRELES)

  ! 1.4 set ISYM_FLAG = 2 (row) or 3 (column) as unsymmetric
  UHM_ISYM_TEMP = ISYM_FLAG     
  ISYM_FLAG     = 3 

  ! Step 2 : Dumping Connectivity
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  ! 2.1 initialize the max dof
  MAXDOFM = 0; 
  MAXDOFC = 0; 
  MAXDOFS = 0;

  ! 2.2 loop over all elements
  write(*,*) 'uhm_solve: BEGIN OF FEEDING CONNECTIVITY '
  call uhm_time_in

  mdle = 0; ime = 0
  do iel=1, NRELES
     call nelcon(mdle, mdle)

     ! 2.2.0 if this element is not in the compute domain, skip it.
     call find_domain(mdle, ndom)
     if (.not.is_compute_domain(ndom)) then
        ! write(*,*) ' Mdle is cycled' , mdle, ndom
        cycle
     end if

     ! 2.2.1 celem catch modified degree of freedoms excluding all dirichlet BC
     call celem( &
          mdle, 1, &
          nrdofs,nrdofm,nrdofc, &
          nodm,ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm, &
          ZVOID,ZVOID)

     ! 2.2.2 update max dof for each physical attributes
     do i=1, NR_PHYSA
        MAXDOFS(i) = max0(MAXDOFS(i), nrdofs(i))
     end do

     MAXDOFM = max0(MAXDOFM, nrdofm)
     MAXDOFC = max0(MAXDOFC, nrdofc)

!      if (nrdofc.gt.UHM_MAX_ELT_NODES) then
!         write(*,*) 'Increase UHM_MAX_ELT_NODES = ', UHM_MAX_ELT_NODES 
!         write(*,*) '# Modified Element Nodes   = ', nrdofc
!      end if

     ime = ime + 1 

     ! 2.2.3 if element has non-zero dofs interface to UHM
     if (nrdofc.gt.0) then

        ! counters
        k = 0; mm = 0

        ! 2.2.5 loop over physics :: fill UHM_ELTS
        do ino=nrnodm, 1, -1
           if (ndofmH(ino).gt.0) then
              k = k + 1
              mm = mm + ndofmH(ino)
              UHM_ELTS(ime)%nodes(1,k)  = nodm(ino)
              UHM_ELTS(ime)%nodes(2,k)  = UHM_PHYSICS_CONTIN
              UHM_ELTS(ime)%idofs(k)    = ndofmH(ino)
              UHM_ELTS(ime)%ikinds(k)   = 0
              UHM_ELTS(ime)%iweights(k) = 0 
           end if
        end do
        do ino=nrnodm, 1, -1
           if (ndofmE(ino).gt.0) then
              k = k + 1
              mm = mm + ndofmE(ino)
              UHM_ELTS(ime)%nodes(1, k)  = nodm(ino)
              UHM_ELTS(ime)%nodes(2, k)  = UHM_PHYSICS_TANGEN
              UHM_ELTS(ime)%idofs(k)     = ndofmE(ino)
              UHM_ELTS(ime)%ikinds(k)    = 0
              UHM_ELTS(ime)%iweights(k)  = 0
           end if
        end do
        do ino=nrnodm, 1, -1
           if (ndofmV(ino).gt.0) then
              k = k + 1
              mm = mm + ndofmV(ino)
              UHM_ELTS(ime)%nodes(1, k)  = nodm(ino)
              UHM_ELTS(ime)%nodes(2, k)  = UHM_PHYSICS_NORMAL
              UHM_ELTS(ime)%idofs(k)     = ndofmV(ino)
              UHM_ELTS(ime)%ikinds(k)    = 0
              UHM_ELTS(ime)%iweights(k)  = 0 
           end if
        end do
        do ino=nrnodm, 1, -1
           if (ndofmQ(ino).gt.0) then
              k = k + 1
              mm = mm + ndofmQ(ino)
              UHM_ELTS(ime)%nodes(1, k)  = nodm(ino)
              UHM_ELTS(ime)%nodes(2, k)  = UHM_PHYSICS_DISCON
              UHM_ELTS(ime)%idofs(k)     = ndofmQ(ino)
              UHM_ELTS(ime)%ikinds(k)    = 0
              UHM_ELTS(ime)%iweights(k)  = 0 
           end if
        end do

        ! 2.2.6 store the effective nodes
        UHM_ELTS(ime)%id      = mdle
        UHM_ELTS(ime)%n_nodes = k
        UHM_ELTS(ime)%n_dofs  = mm

	if (mm.ne.nrdofc) then	
	   write(*,*) ' mm, nrdofc = ', mm, nrdofc
	   call pause
	end if     
        
        if (UHM_SOLVER_UPDATE_ENABLED) then
           call uhm_solver_copy_in_test( &
                UHM_SOLVER_PTR, &
                UHM_ELTS(ime)%nodes, UHM_NIDS, UHM_ELTS(ime)%n_nodes, &
                UHM_ELTS(ime)%itest)
        else
           UHM_ELTS(ime)%itest = .FALSE.
        end if
        
        ! print statement
        if (UHM_VERBOSE) then
           write(*,*) 'mdle = ',mdle
           write(*,*) 'nrdofs(1:',NR_PHYSA,')=',nrdofs(1:NR_PHYSA)
           write(*,*) 'nrdofm=',nrdofm
           write(*,*) 'nrdofc=',nrdofc
           write(*,*) 'nodm(1:',nrnodm,')=', nodm(1:nrnodm)
           write(*,*) 'ndofmH(1:',nrnodm,')=', ndofmH(1:nrnodm)
           write(*,*) 'ndofmE(1:',nrnodm,')=', ndofmE(1:nrnodm)

           call uhm_element_disp(NSTD_OUT, UHM_ELTS(ime))
           call uhm_hold
           cycle
        end if

     end if
  end do
  UHM_N_ELTS = ime
  call uhm_time_out(t)
  write(*,*) 'uhm_solve: END OF FEEDING CONNECTIVITY , ', UHM_N_ELTS
  if (UHM_VERBOSE) then
     call uhm_hold
  end if

  write(*,*) 'uhm_solve: FEEDING CONNECTIVITY TIME ', t, ' SEC'
  write(*,*) ' '

  ! Step 3 : Allocate Storage for Assembly with MAX DOFS
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  allocate(UHM_ZSTIFF(MAXDOFM*MAXDOFM), UHM_ZXLOAD(MAXDOFM*NR_RHS), stat=istat)
  if (istat.ne.SUCCESS) then
     call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
  end if

  call assembly_alloc

  ! Step 4 : Dump Matrices to UHM 
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  iper = 1
  dper = UHM_N_ELTS/25.d0

  ! 4.1 loop over all elements
  write(*,*) 'uhm_solve: BEGIN OF FEEDING MATRICES  '
  call uhm_time_in
  call uhm_solver_copy_in_begin(UHM_SOLVER_PTR)
  do iel=1, UHM_N_ELTS
     mdle = UHM_ELTS(iel)%id

     if (UHM_ELTS(iel)%itest) then
        cycle
     end if

     if (iel.gt.iper*dper) then
        write(*,*) 'uhm_solve: Feeding matrices ', iper*4, '% done '
        iper = iper + 1
     end if
     
     ! 4.1.1 celem to catch modified matrix
     ! write(*,*) 'iel, mdle = ', iel,mdle, 'UHM remembers = ',   UHM_ELTS(iel)%id
     call celem( &
          mdle, 2, &
          nrdofs,nrdofm,nrdofc, &
          nodm,ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm, &
          UHM_ZXLOAD,UHM_ZSTIFF)

     ! 4.1.2 if element has non-zero dofs interface to UHM
     if (UHM_VERBOSE) then
        write(*,*) 'mdle = ', mdle, 'UHM remembers = ',   UHM_ELTS(iel)%id
     end if
     
     n_dofs = UHM_ELTS(iel)%n_dofs
     if (n_dofs.eq.nrdofc) then
        call uhm_solver_copy_in_ab( &
             UHM_SOLVER_PTR, &
             UHM_ELTS(iel)%id,    UHM_EIDS, &
             UHM_ELTS(iel)%nodes, UHM_NIDS, &
             UHM_ELTS(iel)%idofs, &
             UHM_ELTS(iel)%iweights, &
             UHM_ELTS(iel)%ikinds, &
             UHM_ELTS(iel)%n_nodes, &
             n_dofs, UHM_ZSTIFF, &
             n_dofs, UHM_ZXLOAD) 
     else
        write(*,*) 'uhm_solve: sanity check fail, mdle = ', mdle
        write(*,*) '           n_dofs.ne.nrdofc = ', n_dofs, nrdofc
        stop 1
     end if

  end do
  call uhm_solver_copy_in_end(UHM_SOLVER_PTR)
  call uhm_time_out(t)
  write(*,*) 'uhm_solve: END OF FEEDING MATRICES '
  if (UHM_VERBOSE) then
     call uhm_solver_content(UHM_SOLVER_PTR)
     call uhm_hold
  end if
  write(*,*) 'uhm_solve: FEEDING MATRICES TIME ', t, ' SEC'
  write(*,*) ' '
 
  ! Step 5 : Analysis
  !~~~~~~~~~~~~~~~~~~~
  call uhm_time_in
  call uhm_solver_create(UHM_SOLVER_PTR)
  call uhm_time_out(t)
  if (UHM_VERBOSE) then
     call uhm_hold
  end if
  write(*,*) 'uhm_solve: ANALYSIS AND CREATE WORKSPACE TIME  ', t, ' SEC'
  write(*,*) ' '
  if (UHM_SOLVER_WRITE_ASSEMBLED_MATRIX_ONLY) then
     call uhm_solver_dump_assembled_matrix(UHM_SOLVER_PTR, trim(PREFIX)//'mat.mm')
  else
    ! Step 6 : Decompose
     !~~~~~~~~~~~~~~~~~~~
     if (UHM_SOLVER_UPDATE_ENABLED) then
        call uhm_time_in
        call uhm_solver_updecompose(UHM_SOLVER_PTR)
        call uhm_time_out(t)
        write(*,*) 'uhm_solve: UPDATE DECOMPOSE TIME  ', t, ' SEC'
        write(*,*) ' '
     else 
        call uhm_time_in
        call uhm_solver_decompose(UHM_SOLVER_PTR)
        call uhm_time_out(t)
        write(*,*) 'uhm_solve: DECOMPOSE TIME  ', t, ' SEC'
        write(*,*) ' '
     end if
     if (UHM_VERBOSE) then
        call uhm_hold
     end if
     
     ! Step 7 : Solve
     !~~~~~~~~~~~~~~~
     call uhm_time_in
     call uhm_solver_solve(UHM_SOLVER_PTR)
     call uhm_time_out(t)
     if (UHM_VERBOSE) then
        call uhm_hold
     end if
     write(*,*) 'uhm_solve: SOLVE TIME  ', t, ' SEC'
     write(*,*) ' '
     
     ! Step 8 : Check Residual
     !~~~~~~~~~~~~~~~~~~~~~~~
     call uhm_time_in
     call uhm_solver_check(UHM_SOLVER_PTR, residual)
     call uhm_time_out(t)
     write(*,*) ' '
     write(*,*) 'uhm_solve: RESIDUAL = ', residual
     if (residual.gt.UHM_RESIDUAL_THRES) then
        write(*,*) 'uhm_solve: WARNING !!!!!!! RESIDUAL IS GREATER THAN THRESHOLD'
        call uhm_hold
     end if
     if (UHM_VERBOSE) then
        call uhm_hold
     end if
     
     write(*,*) 'uhm_solve: CHECK TIME  ', t, ' SEC'
     write(*,*) ' '
     
     
     ! Step 9 : Solution out
     !~~~~~~~~~~~~~~~~~~~~~~~
     
     ! 9.1 loop over all elements
     write(*,*) 'uhm_solve: BEGIN OF STORING SOLUTIONS'
     call uhm_time_in
     do iel=1, UHM_N_ELTS
        
        ! counters for H1, Hcurl, Hdiv, L2 components supported by "node"
        mvarH(1:MAXNODM)=0 ; mvarE(1:MAXNODM)=0 ; mvarV(1:MAXNODM)=0 ; mvarQ(1:MAXNODM)=0
        
        ! 9.1.2 if element has non-zero dofs interface to UHM
        n_dofs = UHM_ELTS(iel)%n_dofs
        
        ! 9.1.3 catch solution from UHM
        UHM_ZXLOAD = ZERO
        call uhm_solver_copy_out_nodes( &
             UHM_SOLVER_PTR, &
             UHM_ELTS(iel)%nodes, UHM_NIDS, UHM_ELTS(iel)%idofs, &
             UHM_ELTS(iel)%n_nodes, &
             n_dofs, UHM_ZXLOAD)
        
        if (UHM_VERBOSE) then
           call uhm_element_disp(NSTD_OUT, UHM_ELTS(iel))        
           do i=1,NR_RHS
              write(*,*) 'solution, size, col=', n_dofs,i
              write(*,*) UHM_ZXLOAD((i-1)*n_dofs+1:i*n_dofs)
           end do
           call uhm_hold
        end if
        
        ! 9.1.4 
        do irhs=1, NR_RHS
           i = 0
           
           do ino=1, UHM_ELTS(iel)%n_nodes
              nod  = UHM_ELTS(iel)%nodes(1, ino)
              iph  = UHM_ELTS(iel)%nodes(2, ino)
              
              call find_ndof(nod, ndofH, ndofE, ndofV, ndofQ)
              call get_index(nod, idx)
              
              if (UHM_VERBOSE) then
                 write(*,*) 'uhm_solve: nod =', nod, ' ndofH,E,V,Q = ', &
                      ndofH, ndofE, ndofV, ndofQ, '  idx = ', idx
              endif

              
              select case(iph)
              case(UHM_PHYSICS_CONTIN)
                 do j=1, ndofH
                    ivar=mvarH(ino)
                    
                    ! loop over phys
                    do k=1,NRINDEX
                       select case(idx(k))
                       case(1) ! H1 component with Dirchlet BC  
                          ! increment index, do no store
                          ivar=ivar+1 
                          
                       case(2)! H1 component 
                          ivar=ivar+1
                          
                          ! store
                          i=i+1
                          NODES(nod)%zdofH(ivar,j) = UHM_ZXLOAD(n_dofs*(irhs-1)+i)  
                       end select
                    end do
                 enddo
                 mvarH(ino)=ivar
                 
              case(UHM_PHYSICS_TANGEN)
                 do j=1, ndofE
                    ivar=mvarE(ino)
                    
                    ! loop over phys
                    do k=1,NRINDEX
                       select case(idx(k))
                       case(3) ! Hcurl component with Dirchlet BC  
                          ! increment index, do no store
                          ivar=ivar+1 
                          
                       case(4)! Hcurl component 
                          ivar=ivar+1
                          
                          ! store
                          i=i+1
                          NODES(nod)%zdofE(ivar,j) = UHM_ZXLOAD(n_dofs*(irhs-1)+i)  
                       end select
                    end do
                 enddo
                 mvarE(ino)=ivar
                 
              case(UHM_PHYSICS_NORMAL)
                 do j=1, ndofV
                    ivar=mvarV(ino)
                    
                    ! loop over phys
                    do k=1,NRINDEX
                       select case(idx(k))
                       case(5) ! Hdiv component with Dirchlet BC  
                          ! increment index, do no store
                          ivar=ivar+1 
                          
                       case(6)! Hdiv component 
                          ivar=ivar+1
                          
                          ! store
                          i=i+1
                          NODES(nod)%zdofV(ivar,j) = UHM_ZXLOAD(n_dofs*(irhs-1)+i)  
                       end select
                    end do
                 enddo
                 mvarV(ino)=ivar
                 
              case(UHM_PHYSICS_DISCON)
                 do j=1, ndofQ
                    ivar=mvarQ(ino)
                    
                    ! loop over phys
                    do k=1,NRINDEX
                       select case(idx(k))
                       case(7) ! L2 component with Dirchlet BC  
                          ! increment index, do no store
                          ivar=ivar+1 
                          
                       case(8)! L2 component 
                          ivar=ivar+1
                          

                          ! store
                          i=i+1
                          NODES(nod)%zdofQ(ivar,j) = UHM_ZXLOAD(n_dofs*(irhs-1)+i)  

                       end select
                    end do
                 enddo
                 mvarQ(ino)=ivar
                 
              case default
                 write(*,*) 'uhm_solve: NOT SUPPORT, ', &
                      ' node = ', nod, ' phy  = ', iph
                 stop 1
              end select
           end do
        end do
     end do
     call uhm_time_out(t)
     write(*,*) 'uhm_solve: END OF STORING SOLUTIONS FOR ELTS'
     write(*,*) 'uhm_solve: SOLUTION OUT TIME  ', t, ' SEC'
     write(*,*) ' '
  end if

  ! Step 9 : Finalize
  !~~~~~~~~~~~~~~~~~~~~~~~~
  if (UHM_SOLVER_REPORT_SHOW) then
     call uhm_solver_report(UHM_SOLVER_PTR)
  end if
  if (UHM_SOLVER_WRITE_ASSEMBLED_MATRIX) then
     call uhm_solver_dump_assembled_matrix(UHM_SOLVER_PTR, 'mat.ijv')
  end if
  call uhm_solver_clear(UHM_SOLVER_PTR)

  ISYM_FLAG = UHM_ISYM_TEMP

  deallocate(UHM_ZSTIFF, UHM_ZXLOAD, stat=istat)
  if (istat.ne.SUCCESS) then
     call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
  end if

  call uhm_dealloc
  call assembly_dealloc
  call assembly_end

end subroutine uhm_solve
