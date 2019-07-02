!--------------------------------------------------------
!> Purpose : contains the environmental variables for UHM
module uhm2
  use iso_c_binding, only: c_ptr
  use error 

  ! Copied from following location to reduce stupid dependency 
  ! include"uhm/wrapper/fort/parameters.f90"

  !! Mat_Base_::Type_
  INTEGER, PARAMETER :: UHM_FLOAT               =  100
  INTEGER, PARAMETER :: UHM_DOUBLE              =  101
  INTEGER, PARAMETER :: UHM_FCMPLX              =  102
  INTEGER, PARAMETER :: UHM_DCMPLX              =  103
  INTEGER, PARAMETER :: UHM_INT                 =  104

  INTEGER, PARAMETER :: UHM_FULL_MATRIX         =  601
  INTEGER, PARAMETER :: UHM_LOWER_TRIANGULAR    =  300
  INTEGER, PARAMETER :: UHM_UPPER_TRIANGULAR    =  301

  !! Mesh_Base_::Node_Type_
  INTEGER, PARAMETER :: UHM_NODE_DEFAULT        =   0
  INTEGER, PARAMETER :: UHM_NODE_BOUNDARY       =   1 

  !! Time Line
  INTEGER, PARAMETER :: UHM_TIME_LEVEL_USER     =   1
  INTEGER, PARAMETER :: UHM_TIME_LEVEL_FRONT    =   2
  INTEGER, PARAMETER :: UHM_TIME_LEVEL_INTERNAL =   3

  !! Mesh_Base_::Element_Marker_
  INTEGER, PARAMETER :: UHM_ELEMENT_MARKER_ZERO =   0
  INTEGER, PARAMETER :: UHM_ELEMENT_REUSED      =   1

  !! Solver type
  INTEGER, PARAMETER :: UHM_DIRECT_CHOL_NOPIV     =   11
  INTEGER, PARAMETER :: UHM_DIRECT_CHOL_INCPIV    =   12
  INTEGER, PARAMETER :: UHM_DIRECT_CHOL_PIV       =   13
  INTEGER, PARAMETER :: UHM_DIRECT_LDL_NOPIV      =   21
  INTEGER, PARAMETER :: UHM_DIRECT_LDL_INCPIV     =   22
  INTEGER, PARAMETER :: UHM_DIRECT_LDL_PIV        =   23
  INTEGER, PARAMETER :: UHM_DIRECT_LDL_2X2_INC    =   32
  INTEGER, PARAMETER :: UHM_DIRECT_LDL_2X2_PIV    =   33
  INTEGER, PARAMETER :: UHM_DIRECT_LU_NOPIV       =   41
  INTEGER, PARAMETER :: UHM_DIRECT_LU_INCPIV      =   42
  INTEGER, PARAMETER :: UHM_DIRECT_LU_PIV         =   43

  integer, parameter :: UHM_PHYSICS_CONTIN = 1
  integer, parameter :: UHM_PHYSICS_TANGEN = 2
  integer, parameter :: UHM_PHYSICS_NORMAL = 3
  integer, parameter :: UHM_PHYSICS_DISCON = 4

  integer, parameter :: UHM_MAX_ELT_NODES = 256
  integer, parameter :: UHM_NIDS = 2
  integer, parameter :: UHM_EIDS = 1

  !! Pointer to the solver object 
  type (c_ptr), dimension(2) :: UHM_SOLVER_PTR

  logical :: &
       UHM_VERBOSE                            = .FALSE., &
       UHM_SOLVER_UPDATE_ENABLED              = .FALSE., &
       UHM_SOLVER_REPORT_SHOW                 = .FALSE., &
       UHM_SOLVER_WRITE_ASSEMBLED_MATRIX      = .FALSE., &
       UHM_SOLVER_WRITE_ASSEMBLED_MATRIX_ONLY = .FALSE.
  type uhm_element
     logical :: itest
     integer :: id, n_nodes, n_dofs
     integer, dimension(UHM_NIDS, UHM_MAX_ELT_NODES) :: nodes
     integer, dimension(UHM_MAX_ELT_NODES)   :: idofs, ikinds, iweights
  end type uhm_element

  integer UHM_N_ELTS 
  type (uhm_element), target, allocatable :: UHM_ELTS(:)

  ! stck to 1D array structure... it is confusing in communicating with celem.
#if C_MODE
  complex*16, allocatable :: UHM_ZSTIFF(:),UHM_ZXLOAD(:)
  complex*16              :: ZVOID
#else
  real*8,     allocatable :: UHM_ZSTIFF(:),UHM_ZXLOAD(:)
  real*8                  :: ZVOID
#endif

  integer :: UHM_ISYM_TEMP 
  real*8 :: UHM_RESIDUAL_THRES = 1.0e-8

contains
  !--------------------------------------------------------
  subroutine uhm_is_petsc_interfaced(Iflag)
    implicit none
    logical, intent(out) :: Iflag
    integer :: i
    i = transfer(UHM_SOLVER_PTR(2), 1)
    if (i.gt.100) then
       Iflag = .TRUE.
    else
       Iflag = .FALSE.
    end if
  end subroutine uhm_is_petsc_interfaced
  
  subroutine uhm_element_disp(Nout, Elt)
    implicit none
    integer, intent(in) :: Nout
    type (uhm_element), intent(in) :: Elt
    write(Nout,*) ' - uhm element - '
    write(Nout,*) ' id, n_nodes, n_dofs, itest : ', Elt%id,Elt%n_nodes,Elt%n_dofs, Elt%itest
    write(Nout,*) ' nodes  : ', Elt%nodes(1,1:Elt%n_nodes)
    write(Nout,*) ' phys   : ', Elt%nodes(2,1:Elt%n_nodes)
    write(Nout,*) ' idofs  : ', Elt%idofs(1:Elt%n_nodes)
    write(Nout,*) ' kinds  : ', Elt%ikinds(1:Elt%n_nodes)
    write(Nout,*) ' weights: ', Elt%iweights(1:Elt%n_nodes)
    write(Nout,*) ' ' 
  end subroutine uhm_element_disp

  subroutine uhm_alloc(Nelts)
    implicit none
    integer, intent(in) :: Nelts
    integer :: istat

    allocate(UHM_ELTS(Nelts), stat=istat)
    if (istat.ne.SUCCESS) then
       call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
    endif
  end subroutine uhm_alloc

  subroutine uhm_dealloc
    implicit none
    integer :: istat
    deallocate(UHM_ELTS, stat=istat)
    if (istat.ne.SUCCESS) then
       call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
    endif
  end subroutine uhm_dealloc

end module uhm2
