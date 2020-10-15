

subroutine nl_elast_solve(N_incr,Step_size,Abs_tol,Rel_tol,Maxiter,Nupdate,Lsflag,Bfgs,pvflag)
use nl_solver_module
use hyperelasticity
use common_prob_data
use physics
implicit none
integer, intent(in) :: N_incr,Maxiter,Nupdate
real*8 , intent(in) :: Step_size,Abs_tol,Rel_tol
logical, intent(in) :: Lsflag,Bfgs
logical, intent(in),optional :: pvflag

integer :: iParAttr(NR_PHYSA)
logical :: icustom(NR_DISP), pvfloc
integer :: n,iter,i,m
real*8  :: tol,res_norm,res_norm0,lstol, update_norm

! check if paraview plot flag is given. If not, default is false
if (present(pvflag)) then
   pvfloc = pvflag
else
   pvfloc=.false.
endif

iParAttr = NR_COMP
! ! do not dump trace vaiables
! do i=1,NR_PHYSA
!    if ( PHYSAi(i) ) iParAttr(i) = 0
! enddo
! Variables in DISP_ATTR have the following order
!(/'u_co1','u_co2','u_co3','u_nrm', &
! 't_co1','t_co2','t_co3','t_nrm', &
! 'det_F','press','sigm1','sigm2', &
! 'sigm3','lmda1','lmda2','lmda3'/)
! icustom(i) is true if we want to dump to paraview the i-th variable of DISP_ATTR
icustom = .false.
icustom(9:16) = .true.

write(*,*) 'nl_elast_solve: STARTING...'
write(*,*) ''

lstol = 0.5d0

call set_nl_solver_params(lstol,lsflag,bfgs)


! print parameters
write(*,*) 'nl_elast_solve: N_incr          =',N_incr
write(*,*) 'nl_elast_solve: Step_size       =',Step_size
write(*,*) 'nl_elast_solve: Abs_tol         =',Abs_tol
write(*,*) 'nl_elast_solve: Rel_tol         =',Rel_tol
write(*,*) 'nl_elast_solve: MAXITER         =',Maxiter
write(*,*) 'nl_elast_solve: NUPDATE         =',Nupdate
write(*,*) 'nl_elast_solve: LINESEARCH_FLAG =',LINESEARCH_FLAG
if (LINESEARCH_FLAG) write(*,*) 'nl_elast_solve: LINESEARCH_TOL  =',LINESEARCH_TOL
write(*,*) 'nl_elast_solve: BFGS_FLAG       =',BFGS_FLAG

!  set ground load_factor
if (.not.SAVE_LOAD_FACTOR) then
   LOAD_FACTOR_PREV = 0.d0
else
   LOAD_FACTOR_PREV = LOAD_FACTOR
endif

! prepare connectivity and array sizing (using celem, etc.) - THIS ASSUMES A CONSTANT MESH ALL ALONG
call prep_mumps('G')

! allocate memory for solution and residual
allocate(RES_GLOBAL(NRDOF_CON))


! allocate arrays for BFGS
if ( BFGS_FLAG ) call alloc_BFGS(NRDOF_CON,Nupdate)


! set initial guess(stored in NEW_SOL) for first iteration
if (FIRST_NLSOLVE) then
!  if it's the first call, the guess is zero
   ! NEW_SOL = 0.d0
   ! call map_dofs_to_local(NEW_SOL)
   call reset_coms
   write(*,*) 'nl_elast_solve: dof values have been reset'
   if (pvfloc) call my_paraview_driver(iParAttr,icustom)
! else
!  if not the first call, the guess is assigned to be the DOF values in memory
   ! call map_dofs_to_global(NEW_SOL)
endif

do n=1,N_incr
!
   write(*,*)'::::::::::::::::::::::::::::::::::::::::::::::::::::::::'
   write(*,4002) n,N_incr
 4002 format('nl_elast_solve: beginning load step =',i4,' out of ',i4,/)
! 
!  given the increment, a change in a factor of the load and b.c. functions is made
   call load_increment(n,Step_size)

!  update the Dirichlet dof in both coms
   call update_Ddof
   write(*,*) 'nl_elast_solve: Dirichlet dof have been updated with new load factor'

   if (NRDOF_CON.lt.1) then
      write(*,*) 'nl_elast_solve: there is no DOF to solve'
      cycle
   endif
! 
!  begin itetration count
   iter = 0
!  First call to assembly to obtain tangent A and global residual RES_GLOBAL
   call dealloc_schur
   call assembly_mumps(2, NRDOF_CON, RES_GLOBAL)
   write(*,*) 'nl_elast_solve: initial global residual assembled'
   write(*,*) 'nl_elast_solve: initial global matrix assembled'



   if (BFGS_FLAG) RES_PREV = RES_GLOBAL

 !   write(*,*) 'MUMPS_PAR%NNZ=',MUMPS_PAR%NNZ

 !   ! write(*,*) 'MUMPS_PAR%ICNTL(6)=',MUMPS_PAR%ICNTL(6)
 !   write(*,*) 'IRN JCN A'
 !   do i=1,MUMPS_PAR%NNZ
 !      write(*,3993) MUMPS_PAR%IRN(i) , MUMPS_PAR%JCN(i) , MUMPS_PAR%A(i)
 ! 3993 format(i4,i4,es22.15)
 !   enddo







!  get LU factorization of tangent A through MUMPS
   call factorize_tangent_mumps
   write(*,*) 'nl_elast_solve: initial global matrix factorized'

!  get norm of initial residual vector
   call resid_norm(RES_GLOBAL,NRDOF_CON,res_norm0)

   write(*,4000) ' nl_elast_solve: initial residual = ',res_norm0
 4000 format(A,es10.3)
! 
   res_norm = res_norm0
   tol = Rel_tol!**(0.5d0)
   if (n.eq.N_incr) tol = Rel_tol
   do while(    iter.lt.MAXITER       .and. &
            res_norm.gt.tol*res_norm0 .and. &
            res_norm.gt.Abs_tol            )
!     count iteration
      iter = iter + 1

      write(*,*) ''
      write(*,*) 'nl_elast_solve: beginning iteration ',iter
      write(*,*) ''
!  ...MODIF NEWTON METHOD: only factorize tangent A if iter is multiple of NUPDATE
      m = mod(iter-1,NUPDATE)+1
      if ( iter.gt.1 .and. m.eq.1 ) then
         call dealloc_schur
         call assembly_mumps(2,NRDOF_CON,RES_GLOBAL)
         write(*,*) ''
         write(*,*) 'nl_elast_solve:   global matrix assembled for iteration ',iter
!        get LU factorization of tangent A through MUMPS
         call factorize_tangent_mumps
         write(*,*) 'nl_elast_solve:  global matrix factorized for iteration ',iter
         write(*,*) ''
      endif

   ! write(*,*) 'RES_GLOBAL'
   ! write(*,*) RES_GLOBAL(1:18)
!
! 
!     compute BFGS updates
      if ( BFGS_FLAG .and. m.gt.1 ) then
!        update current residual
         RES_CURR = RES_GLOBAL
!        compute and store vectors w and w
         call BFGS_UPDATES(iter,Nupdate,NRDOF_CON,MUMPS_PAR%RHS,RES_PREV,RES_CURR)
      endif
!
!
!     apply BFGS updates to residual
      if ( BFGS_FLAG .and. m.gt.1 )                                              &
      call BFGS_RES(iter,Nupdate,NRDOF_CON,RES_CURR,RES_GLOBAL)
!
!
!     since the tangent A is already factorized, we ask mumps to solve A*delta = res
!     solution delta is stored in MUMPS_PAR%RHS
      call solve_delta_mumps(RES_GLOBAL)


   ! write(*,*) 'MUMPS_PAR%RHS'
   ! write(*,*) MUMPS_PAR%RHS(1:18)
!
!
!     apply BFGS updates to solution delta
      if ( BFGS_FLAG .and. m.gt.1 ) then
         call BFGS_SOL(iter,Nupdate,NRDOF_CON,MUMPS_PAR%RHS)
         ! retrieve the original residual
         RES_GLOBAL = RES_CURR
      endif
!
!
!     store the delta dofs in first slot (and perform the backward static condensation)
      call map_dofs_to_local(MUMPS_PAR%RHS)
      write(*,*) 'nl_elast_solve:              delta solved for iteration ',iter
!     we compute factor s for line search (if not enabled, this sets LINESEARCH_FACTOR=1)
      call get_delta_factor
      ! call get_step_illinois
      write(*,*) 'nl_elast_solve: line s. factor determined for iteration ',iter
!     update solution: slot 2 += ls_factor*delta and reset first slot
      call update_sol(LINESEARCH_FACTOR)
      write(*,*) 'nl_elast_solve:       new solution stored for iteration ',iter
!
!     obtain new RES_GLOBAL - if we're using LineSearch we already have the residual
      if (.not.LINESEARCH_FLAG)  then
         call dealloc_schur
         call assembly_mumps(1,NRDOF_CON, RES_GLOBAL)
      endif

      write(*,*) 'nl_elast_solve:    new residual assembled for iteration ',iter
! 
!
!     for BFGS store copy of current residual
      if( BFGS_FLAG .and. m.gt.1 )  RES_PREV = RES_CURR
!
! 
!     get norm of residual vector
      call resid_norm(RES_GLOBAL,NRDOF_CON,res_norm)
      write(*,4003) res_norm,iter
 4003 format(   ' nl_elast_solve: residual norm =',es10.3,' for iteration ',i12)
!
      call resid_norm(MUMPS_PAR%RHS,NRDOF_CON,update_norm)
      write(*,4004) LINESEARCH_FACTOR * update_norm , iter
 4004 format(   ' nl_elast_solve:sol.update norm=',es10.3,' for iteration ',i12)
! 
      write(*,*) 'nl_elast_solve:                     completed iteration ',iter
! 
   enddo
   write(*,*) '::::::::::::::::::::::::::::::::::::::::::::::'
   write(*,*) ''
   write(*,4000) 'nl_elast_solve: final residual norm  = ',res_norm
   write(*,4000) 'nl_elast_solve: final residual ratio = ',res_norm/res_norm0

   if (res_norm.le.tol*res_norm0.or.res_norm.le.Abs_tol) then
      write(*,4001) n,iter
 4001 format(/,'nl_elast_solve: convergence attained in step ',i4,' in ',i5,' iterations.',/)

      if (pvfloc) call my_paraview_driver(iParAttr,icustom)

      if (NR_PHYSA.gt.1) then
         write(*,*) 'computing DPG residual...'
         call residual
      endif

   else

      write(*,*) '!!!'
      write(*,*) 'nl_elast_solve: SOLUTION HAS DIVERGED IN STEP ',n
      write(*,*) '!!!'
      goto 99

   endif

enddo

  99  FIRST_NLSOLVE=.false.

  if(.not.SAVE_LOAD_FACTOR) FIRST_NLSOLVE=.true.

call close_solver

if (BFGS_FLAG) call dealloc_BFGS

deallocate(RES_GLOBAL)
write(*,*) ''
write(*,*) 'nl_elast_solve: DONE'
write(*,*) '::::::::::::::::::::'

end subroutine
