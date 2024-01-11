!----------------------------------------------------------------------
! exec_case
!----------------------------------------------------------------------
   subroutine exec_case(idec)
!
      use commonParam
      use control
      use data_structure3D
      use par_mesh
      use paraview      , only: VLEVEL,paraview_select_attr
      use zoltan_wrapper, only: zoltan_w_partition,zoltan_w_eval
      use mpi_param
      use MPI           , only: MPI_COMM_WORLD,MPI_CHARACTER,MPI_INTEGER,MPI_LOGICAL
!
      implicit none
!
      integer, intent(in) :: idec
!
      integer :: nstop
      logical :: solved
!
      integer :: flag(2)
      integer :: physNick
!
      logical :: iPvAttr(2)
      character(len=2) :: vis_level
!
      integer :: fld,numPts,i,mdle,kref,refs
      real(8) :: res
!
      integer :: src,count,ierr
!
!----------------------------------------------------------------------
!
      physNick = 1
      flag=0; flag(2)=1
!
      solved = .false.
!
      select case(idec)
!
!     ...Paraview graphics
         case(3)
!
            iPvAttr = (/.false.,.true./) ! write field output only
            if (RANK .eq. ROOT) then
               write(*,300) ' paraview output: select fields (T/F)' , &
                            '  - 2 H(curl) (E and H flux)'          , &
                            '  - 6 L2      (E and H field)'
               read (*,*) iPvAttr(1), iPvAttr(2)
           300 format(A,/,A,/,A)
            endif
!
            count = 2; src = ROOT
            call MPI_BCAST (iPvAttr,count,MPI_LOGICAL,src,MPI_COMM_WORLD,ierr)
!
            if (RANK .eq. ROOT) then
               write(*,*) 'paraview output: select VLEVEL (0-4)...'
               read (*,*) vis_level
               select case(vis_level)
                  case('0','1','2','3','4')
                     VLEVEL = vis_level
                  case default
                     write(*,*) ' invalid VLEVEL. setting VLEVEL=3 (default).'
                     VLEVEL = '3'
               end select
            endif
!
            count = len(VLEVEL); src = ROOT
            call MPI_BCAST (VLEVEL,count,MPI_CHARACTER,src,MPI_COMM_WORLD,ierr)
!
            call paraview_select_attr(iPvAttr)
            call paraview_driver
            call MPI_BARRIER (MPI_COMM_WORLD, ierr)
!
!     ...print data structure (interactive)
         case(10); call result
!
!     ...print general data structure info
         case(11)
            write(*,110) NRELIS,NRELES,NRNODS
 110        format(' NRELIS,NRELES,NRNODS            = ',3I10)
            write(*,111) NRDOFSH,NRDOFSE,NRDOFSV,NRDOFSQ
 111        format(' NRDOFSH,NRDOFSE,NRDOFSV,NRDOFSQ = ',4I10)
            write(*,112) MAXNODS,NPNODS
 112        format(' MAXNODS,NPNODS                  = ',2I10)
!
!     ...print current partition (elems)
         case(15)
            write(*,*) 'printing current partition (elems)...'
            call print_partition
!
!     ...print current subdomains (nodes)
         case(16)
            write(*,*) 'printing current subdomains (nodes)...'
            call print_subd
!
!     ...print current partition coordinates
         case(17)
            write(*,*) 'printing current partition coordinates...'
            call print_coord
!
!     ...single uniform h-refinement
         case(20)
            write(*,*) 'global h-refinement...'
            call global_href
            if (IBCFLAG .eq. 3) then
               call propagate_flag(3,3)
            endif
            call update_gdof
            call update_Ddof
!
!     ...single uniform p-refinement
         case(21)
            write(*,*) 'global p-refinement...'
            call global_pref
            call update_gdof
            call update_Ddof
!
!     ...distribute mesh
         case(30)
            write(*,*) 'distribute mesh...'
            call distr_mesh
!
!     ...collect dofs on ROOT processor
         case(31)
            write(*,*) 'collecting dofs on ROOT...'
            call collect_dofs
!
!     ...evaluate current partition
         case(33)
            if (DISTRIBUTED) then
               write(*,*) 'evaluating current partition...'
               call zoltan_w_eval
            else
               write(*,*) 'distribute mesh first to use Zoltan...'
            endif
!
!     ...run mesh verification routines
         case(35)
            write(*,*) 'verify distributed mesh consistency...'
            call par_verify
!
!     ...solve problem with par_nested (nested MPI MUMPS)
         case(40)
            write(*,*) 'calling nested MUMPS (MPI) solver...'
            call par_nested('H')
!
!     ...solve problem with par_mumps_sc (MPI MUMPS)
         case(41)
            write(*,*) 'calling MUMPS (MPI) solver...'
            call par_mumps_sc('H')
!
!     ...solve problem with mumps_sc (OpenMP MUMPS)
         case(42)
            write(*,*) 'calling MUMPS (OpenMP) solver...'
            call mumps_sc('H')
!
!     ...solve problem with pardiso_sc (OpenMP Pardiso)
         case(43)
            write(*,*) 'calling Pardiso (OpenMP) solver...'
            call pardiso_sc('H')
!
!     ...solve problem with Frontal solver (sequential)
         case(44)
            write(*,*) 'calling Frontal (Seq) solver...'
            call solve1(1)
!
!     ...solve problem with PETSc solver (MPI)
         case(45)
            write(*,*) 'calling PETSc (MPI) solver...'
            call petsc_solve('H')
!
!     ...compute error for problem with known solution
         case(50)
            if (NEXACT .eq. 0) then
               write(*,*) 'NEXACT=0. Returning...'
            else
               write(*,*) 'computing error...'
               call exact_error(flag,physNick)
            endif
!
!     ...compute the residual
         case(51)
            write(*,*) 'computing residual...'
            call residual(res)
!
!     ...debugging routines
         case(70)
            if (RANK.eq.ROOT) then
               write(*,*) 'Select on mdle node from the list: '
               do i=1,NRELES
                  write(*,2610) ELEM_ORDER(i)
             2610 format(I6)
               enddo
               read(*,*) mdle
            endif
            if (NUM_PROCS .gt. 1) then
               count = 1; src = ROOT
               call MPI_BCAST (mdle,count,MPI_INTEGER,src,MPI_COMM_WORLD,ierr)
            endif
            select case(NODES(mdle)%ntype)
               case(MDLB); kref = 110
               case(MDLP); kref = 10
            end select
            call refine(mdle,kref)
            call close_mesh
            call update_gdof
            call update_Ddof
         case(71)
            do i = 1,13
               write(*,100) 'Iteration ', i
               mdle = ELEM_ORDER(NRELES/2)
               select case(NODES(mdle)%ntype)
                  case(MDLB); kref = 110
                  case(MDLP); kref = 10
               end select
               call refine(mdle,kref)
               call close_mesh
               call update_gdof
               call update_Ddof
            enddo
        100 format(/,'/////////////////////////////////////////////////////////////', &
                   /,'             ',A,I4,/)
!
         case default
            write(*,*) 'exec_case: unknown case...'
      end select
!
   end subroutine exec_case

