!--------------------------------------------------------
!> Purpose : contains the environmental variables for UHM
module uhm

  use kinds
  use error

  !------------------------------------------------------
  ! ** type def
  type uhm_elt
     integer(kind=C_PTR) :: id
     integer :: n_nods, mm
     integer, dimension(2,60) :: nods
     integer, dimension(60)   :: n_dof, p
  end type uhm_elt

  ! ** Two ways include fortran interface from UHM
  !    1. Copy paste all parameters from uhm/wrapper/fort.f90
  !    2. Use include statement for all subroutines.

  ! ** From UHM wrapper
  integer, parameter :: UHM_REAL               = 101
  integer, parameter :: UHM_COMPLEX            = 103

  integer, parameter :: UHM_FULL_MATRIX        = 100
  integer, parameter :: UHM_LOWER_TRIANGULAR   = 200
  integer, parameter :: UHM_UPPER_TRIANGULAR   = 300

  integer, parameter :: UHM_LHS                = 1
  integer, parameter :: UHM_RHS                = 2

  integer, parameter :: UHM_LEAF_TO_ROOT       = 1
  integer, parameter :: UHM_ROOT_TO_LEAF       = 0

  integer, parameter :: UHM_NODE_KIND_DEFAULT  = 0
  integer, parameter :: UHM_NODE_KIND_BOUNDARY = 1

  integer, parameter :: UHM_UNASSEMBLED        = 0
  integer, parameter :: UHM_ASSEMBLED          = 1

  integer, parameter :: UHM_PHYSICS_SINGLE     = 1
  integer, parameter :: UHM_PHYSICS_MULTI      = 2

  integer, parameter :: UHM_CHOL               = 1
  integer, parameter :: UHM_CHOL_PIV           = 2
  integer, parameter :: UHM_CHOL_INCPIV        = 3

  integer, parameter :: UHM_LDL_NOPIV          = 4
  integer, parameter :: UHM_LDL_PIV            = 5
  integer, parameter :: UHM_LDL_INCPIV         = 6

  integer, parameter :: UHM_LU_NOPIV           = 7
  integer, parameter :: UHM_LU_PIV             = 8
  integer, parameter :: UHM_LU_INCPIV          = 9

  integer, parameter :: UHM_QR                 = 10

  integer, parameter :: UHM_SYMMETRY   = 0
  integer, parameter :: UHM_UNSYMMETRY = 1

  character(len=9), dimension(10), parameter :: UHM_METHOD_NAME = &
       (/'CHOL     ','CHOL_PIV ','CHOL_IPIV', &
         'LDL_NOPIV','LDL_PIV  ','LDL_INPIV', &
         'LU_NOPIV ','LU_piv   ','LU_INCPIV', &
         'QR       '/)

  ! ** parameter for physics variables
  integer, parameter :: UHM_PHYSICS_CONTIN = 1
  integer, parameter :: UHM_PHYSICS_TANGEN = 2
  integer, parameter :: UHM_PHYSICS_NORMAL = 3
  integer, parameter :: UHM_PHYSICS_DISCON = 4

  integer, parameter :: UHM_INTERF_PARDISO = 11
  integer, parameter :: UHM_INTERF_WSMP    = 12
  integer, parameter :: UHM_INTERF_MUMPS   = 13

  integer, dimension(2), parameter :: UHM_SYM_MUMPS   = (/2, 0/)

  !------------------------------------------------------
  ! ** variables
  integer :: UHM_METHOD, UHM_N_RHS, UHM_DATATYPE, UHM_ISYM_FLAG

  ! pointer to mesh object
  integer(kind=C_PTR) :: UHM_MESH, UHM_PARDISO, UHM_WSMP, UHM_MUMPS

  ! relavent uhm objects to FE
  type (uhm_elt), target, allocatable :: UHM_ELTS(:)

  !
#if C_MODE
  complex*16, allocatable :: UHM_ZSTIFF(:),UHM_ZXLOAD(:)
#else
  real*8,     allocatable :: UHM_ZSTIFF(:),UHM_ZXLOAD(:)
#endif

  !$OMP THREADPRIVATE (UHM_ZSTIFF, UHM_ZXLOAD)

contains


  !--------------------------------------------------------
  !> Purpose : display all possible methods
  !! @param[in] Nout - display device name
  subroutine uhm_method_disp_all(Nout)
    implicit none
    integer, intent(in) :: Nout
    integer :: i
    do i=1,10
       call uhm_method_disp(Nout, i)
    enddo
  end subroutine uhm_method_disp_all

  !--------------------------------------------------------
  !> Purpose : display the method name for given method number
  !! @param[in] Nout    - display device name
  !! @param[in] Nmethod - method to be shown
  subroutine uhm_method_disp(Nout, Nmethod)
    implicit none
    integer, intent(in) :: Nout, Nmethod
    write(Nout,*) 'Method : ', Nmethod, UHM_METHOD_NAME(Nmethod)
  end subroutine uhm_method_disp

  !--------------------------------------------------------
  !> Purpose : display the element
  !! @param[in] Elt  - type 'uhm_elt' to display
  !! @param[in] Nout - display device name
  subroutine uhm_elt_disp(Nout, Elt)
    implicit none
    integer, intent(in) :: Nout
    type (uhm_elt), intent(in) :: Elt
    write(Nout,*) ' - uhm elt - '
    write(Nout,*) ' id : ', Elt%id
    write(Nout,*) ' n_nods :', Elt%n_nods, ' mm : ', Elt%mm
    write(Nout,*) ' id(1)  : ', Elt%nods(1,1:Elt%n_nods)
    write(Nout,*) ' id(2)  : ', Elt%nods(2,1:Elt%n_nods)
    write(Nout,*) ' n_dof  : ', Elt%n_dof(1:Elt%n_nods)
    write(Nout,*) ' '
  end subroutine uhm_elt_disp


  !--------------------------------------------------------
  !> Purpose : initialize flame and create uhm mesh object
  !! @param[in] N_thread - number of threads to use
  !! @param[in] N_blocksize - internal subblock matrix size
  !! @param[in] N_method - factorization method( UHM_CHOL, UHM_LU_PIV, UHM_QR and so on )
  subroutine uhm_initialize(N_thread, N_blocksize, N_method)
    use control, only:ISYM_FLAG
    implicit none
    integer, intent(in) :: N_thread, N_blocksize, N_method
    integer :: n_datatype

#if C_MODE
    n_datatype = UHM_COMPLEX
#else
    n_datatype = UHM_REAL
#endif

    call uhm_initialize_fla

    call uhm_mesh_create(UHM_MESH)
    call uhm_unlock(UHM_MESH)

    call uhm_set_num_threads    (N_thread)
    call uhm_set_hier_block_size(N_blocksize)

    ! back up symmetry flag
    UHM_ISYM_FLAG = ISYM_FLAG

    ! check valid method
    select case (N_method)
    case (1,2,3,4,5,6)
       call uhm_set_symmetry(UHM_MESH, UHM_SYMMETRY)
    case (7,8,9,10)
       call uhm_set_symmetry(UHM_MESH, UHM_UNSYMMETRY)
    case (UHM_INTERF_PARDISO);
    case (UHM_INTERF_WSMP);
    case (UHM_INTERF_MUMPS);
    case default
       write(*,*) 'uhm:: method = ', N_method
       call logic_error(ERR_INVALID_VALUE,__FILE__,__LINE__)
    end select

    UHM_METHOD   = N_method
    UHM_DATATYPE = n_datatype

    select case (UHM_METHOD)
    case (UHM_INTERF_PARDISO);  call uhm_pardiso_create(UHM_PARDISO, UHM_DATATYPE)
    case (UHM_INTERF_WSMP);     call uhm_wsmp_create   (UHM_WSMP,    UHM_DATATYPE)
    case (UHM_INTERF_MUMPS);    call uhm_mumps_create  (UHM_MUMPS,   UHM_DATATYPE)
    end select

  end subroutine uhm_initialize

  !--------------------------------------------------------
  !> Purpose : delete uhm mesh object and finalize flame
  subroutine uhm_finalize
    implicit none
    real*8 :: buffer_used, buffer_max_used

    select case (UHM_METHOD)
    case (UHM_INTERF_PARDISO);
       call uhm_pardiso_finalize(UHM_PARDISO)
       call uhm_pardiso_delete(UHM_PARDISO)
    case (UHM_INTERF_WSMP);
       call uhm_wsmp_finalize(UHM_WSMP)
       call uhm_wsmp_delete(UHM_WSMP)
    case (UHM_INTERF_MUMPS);
       call uhm_mumps_finalize(UHM_MUMPS)
       call uhm_mumps_delete(UHM_MUMPS)
    end select

    call uhm_unlock(UHM_MESH)
    call uhm_mesh_delete(UHM_MESH)
    call uhm_finalize_fla

  end subroutine uhm_finalize

  !--------------------------------------------------------
  !> Purpose : allocate workspace
  !! @param[in] N_elts - number of elements, same as NRELES in datastructure
  subroutine uhm_alloc(N_elts)
    implicit none
    integer, intent(in) :: N_elts
    integer :: istat

    allocate(UHM_ELTS(N_elts), stat=istat)
    if (istat.ne.SUCCESS) then
       call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
    endif
  end subroutine uhm_alloc

  !--------------------------------------------------------
  !> Purpose : deallocate workspace
  subroutine uhm_dealloc
    implicit none
    integer :: istat
    deallocate(UHM_ELTS, stat=istat)
    if (istat.ne.SUCCESS) then
       call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
    endif
  end subroutine uhm_dealloc

  !--------------------------------------------------------
  !> Purpose : set the number of RHS
  !! @param[in] N_rhs - number of RHS, should be same as NR_RHS in assembly
  subroutine uhm_set_nr_rhs(N_rhs)
    implicit none
    integer, intent(in) :: N_rhs
    UHM_N_RHS = N_rhs
  end subroutine uhm_set_nr_rhs

  !--------------------------------------------------------
  !> Purpose : get the number of RHS
  !! @param[out] N_rhs - number of RHS
  subroutine uhm_get_nr_rhs(N_rhs)
    implicit none
    integer, intent(out) :: N_rhs
    N_rhs = UHM_N_RHS
  end subroutine uhm_get_nr_rhs

  !--------------------------------------------------------
  !> Purpose : analyze the sparse system
  subroutine uhm_analyze_uhm
    implicit none

    integer :: &
         n_dof, n_nonzero_factor, i_hier_flag, i_mult_flag, &
         n_blocksize, n_thread;
    real*8  :: &
         t_base, t_after, flop_decompose, flop_solve, buffer

    ! perform multi-level bisection
    call uhm_timer        (t_base)
    call uhm_build_tree   (UHM_MESH)
    call uhm_timer        (t_after)
    t_after = t_after - t_base;

    ! create scheduler and does not allow the change of connectivity
    call uhm_lock         (UHM_MESH)

    !!call uhm_mesh_disp    (UHM_MESH)

    call uhm_create_matrix_without_buffer(UHM_MESH, UHM_DATATYPE, UHM_N_RHS)
    call uhm_create_element_matrix_buffer(UHM_MESH, 0)

    call uhm_get_n_dof    (UHM_MESH,n_dof)
    call uhm_get_n_nonzero_factor(UHM_MESH, UHM_METHOD.eq.UHM_CHOL, n_nonzero_factor)

    call uhm_estimate_cost( UHM_MESH, &
         UHM_METHOD, UHM_DATATYPE, UHM_N_RHS, &
         flop_decompose, flop_solve, buffer)

    call uhm_is_multithreading_enable(i_mult_flag)
    call uhm_is_hier_matrix_enable   (i_hier_flag)

    call uhm_get_num_threads     (n_thread)
    call uhm_get_hier_block_size (n_blocksize)

    write(*,*) '-- Analysis Report --'
    if (i_mult_flag.gt.0) then
       write(*,*) 'Multithreading is enabled'
    else
       write(*,*) 'Multithreading is disabled'
    endif

    if (i_hier_flag.gt.0) then
       write(*,*) 'Hierarchical matrices are enabled'
    else
       write(*,*) 'Hierarchical matrices are disabled'
    endif

    if (UHM_DATATYPE.eq.UHM_REAL) then
       write(*,*) 'REAL datatype'
    else
       write(*,*) 'COMPLEX datatype'
    end if

    write(*,*) 'Number of threads    = ', n_thread
    write(*,*) 'Method               = ', UHM_METHOD_NAME(UHM_METHOD)
    write(*,*) 'Blocksize            = ', n_blocksize
    write(*,*) 'NDOF                 = ', n_dof
    write(*,*) 'Non-zero entries     = ', n_nonzero_factor
    write(*,7001) flop_decompose/1.0e9
    write(*,7002) flop_solve/1.0e9
    ! UHM estimate the buffer amount much smaller. it can misguide
    !!write(*,7003) buffer/1.0e6
    !!write(*,*) 'In addition to above buffer, UHM need a little more workspace'
    write(*,*) '---------------------'
    write(*,7004) t_after
    write(*,*) '---------------------'
    write(*,*) ' '


7001 format('Est. FLOP decompose  = ', f12.5, ' GFLOP')
7002 format('Est. FLOP solve      = ', f12.5, ' GFLOP')
7003 format('Est. buffer used     = ', f12.5, ' MB')
7004 format('Time analysis        = ', f12.5, ' SEC')

  end subroutine uhm_analyze_uhm

  !--------------------------------------------------------
  !> Purpose : decompose the sparse matrix according to a given method
  subroutine uhm_decompose_uhm
    implicit none
    real*8  :: &
         t_base, t_after, buffer_max_used, buffer_used

    call uhm_timer        (t_base)
    select case (UHM_METHOD)
    case(1,2,3,4,5,6,7,8,9,10)
       call uhm_decompose_with_free(UHM_MESH, UHM_METHOD)
    case default
       call logic_error(ERR_INVALID_VALUE, __FILE__,__LINE__)
    end select
    call uhm_timer        (t_after)
    t_after = t_after - t_base

    call uhm_get_buffer_used     (buffer_used)
    call uhm_get_max_buffer_used (buffer_max_used)

    write(*,*) '-- Decomposition Report --'
    write(*,*) 'Method             = ', UHM_METHOD_NAME(UHM_METHOD)

    write(*,7005) buffer_used/1.0e6
    write(*,7006) buffer_max_used/1.0e6
    write(*,*) '---------------------'
    write(*,7007) t_after
    write(*,*) '---------------------'
    write(*,*) ' '

7005 format('Buffer used          = ', f12.5, ' MB')
7006 format('Buffer MAX used      = ', f12.5, ' MB')
7007 format('Time decomposition   = ', f12.5, ' SEC')

  end subroutine uhm_decompose_uhm

  !--------------------------------------------------------
  !> Purpose : perform backward and forward substritution for the factored matrix
  subroutine uhm_solve_uhm
    implicit none
    real*8  :: &
         t_base, t_after

    call uhm_timer        (t_base)
    select case (UHM_METHOD)
    case(1,2,3,4,5,6,7,8,9,10)
       call uhm_solve(UHM_MESH, UHM_METHOD)
    case default
       call logic_error(ERR_INVALID_VALUE, __FILE__,__LINE__)
    end select
    call uhm_timer        (t_after)
    t_after = t_after - t_base

    write(*,*) '-- Solution Report --'
    write(*,*) 'Method             = ', UHM_METHOD_NAME(UHM_METHOD)
    write(*,*) '---------------------'
    write(*,7008) t_after
    write(*,*) '---------------------'
    write(*,*) ' '

7008 format('Time solution        = ', f12.5)

  end subroutine uhm_solve_uhm

  !--------------------------------------------------------
  !> Purpose : calculate residual | Ax - b |
  subroutine uhm_check_residual_uhm
    implicit none
    real*8 :: res
    real*8 :: &
         t_base, t_after

    call uhm_timer        (t_base)
    select case (UHM_METHOD)
    case(1,2,3,4,5,6,7,8,9,10)
       call uhm_check(UHM_MESH, UHM_METHOD)
    case default
       call logic_error(ERR_INVALID_VALUE, __FILE__,__LINE__)
    end select

    res = 0.d0
    call uhm_get_residual(UHM_MESH, res)

    call uhm_timer        (t_after)
    t_after = t_after - t_base

    write(*,*) '-- Residual Report --'
    write(*,*) 'Method              = ', UHM_METHOD_NAME(UHM_METHOD)
    write(*,*) '---------------------'
    write(*,7009) res
    write(*,*) '---------------------'
    write(*,*) ' '

7009 format('Residual             = ', e12.5)

  end subroutine uhm_check_residual_uhm


end module uhm
