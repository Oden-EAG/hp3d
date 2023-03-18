!> @brief Defines multiphysics info
!> @date Mar 2023
module physics

  ! ...number of different physics attributes
  integer, save :: NR_PHYSA

  ! ...list of physics attributes
  character(len=5), save, allocatable :: PHYSA(:)

  ! ...enabled/disabled physics attributes flags
  logical, save, allocatable :: PHYSAm(:)

  ! ...the corresponding number of components
  integer, save, allocatable :: NR_COMP(:)

  ! ...the corresponding discretization type
  integer, save, allocatable :: D_TYPE(:)
  integer, parameter :: CONTIN = 1
  integer, parameter :: TANGEN = 2
  integer, parameter :: NORMAL = 3
  integer, parameter :: DISCON = 4

  !  ...interface variable flag
  logical, save, allocatable :: PHYSAi(:)
  
  !  ...homogeneous Dirichlet BC variable flag
  logical, save, allocatable :: PHYSAd(:)
  !
  !  ...location of the first component for the attribute
  integer, save, allocatable :: ADRES(:)
  !
  !  ...total number of H1,H(curl),H(div) and L2 variables
  integer, save :: NRHVAR,NREVAR,NRVVAR,NRQVAR
  !
  !  ...number of entries in index
  integer, save :: NRINDEX
  !
  !----------------------------------------------------------------------
  !
  !  ...number of different node cases
  !  ...ex: physical attributes = {a,b,c}
  !  .......cases               = {a,b,c; ab,ac,bc; abc}
  integer, save :: NRCASES
  !
  !     number of H1 unknowns per 'case' node
  integer, save, allocatable :: NREQNH(:)
  !
  !     number of H(curl) unknowns per 'case' node
  integer, save, allocatable :: NREQNE(:)
  !
  !     number of H(div) unknowns per 'case' node
  integer, save, allocatable :: NREQNV(:)
  !
  !     number of L^2 unknowns per 'case' node
  integer, save, allocatable :: NREQNQ(:)
  !
  interface dumpout_physics
     module procedure dumpout_physics_to_file
     module procedure dumpout_physics_to_default
  end interface

  interface dumpin_physics
     module procedure dumpin_physics_from_file
     module procedure dumpin_physics_from_default
  end interface
  !
contains
!
!-----------------------------------------------------------------------
!
!> @brief Return string representation of a discretization type
!> @date Mar 2023
      function S_DType(IDtype)
         Integer IDtype
         character(6) S_DType
         select case(IDtype)
            case(CONTIN); S_DType = 'contin'
            case(TANGEN); S_DType = 'tangen'
            case(NORMAL); S_DType = 'normal'
            case(DISCON); S_DType = 'discon'
            case default
               write(*,*) 'S_DType: IDtype = ', IDtype
               stop
         end select
      end function S_DType
!
!-----------------------------------------------------------------------
!
!> @brief Return integer representation of a discretization type
!> @date Mar 2023
      function I_DType(SDtype)
         character(6) SDtype
         integer I_DType
         select case(SDtype)
            case('contin'); I_DType = CONTIN
            case('tangen'); I_DType = TANGEN
            case('normal'); I_DType = NORMAL
            case('discon'); I_DType = DISCON
            case default
               write(*,*) 'I_DType: SDtype = ', SDtype
               stop
         end select
      end function I_DType

  !> Purpose : allocate data structure for multiphysics
  subroutine alloc_physics
    !
    if (allocated(PHYSA)) then
       deallocate(PHYSA,PHYSAm,NR_COMP,D_TYPE,PHYSAi,PHYSAd,ADRES, &
                  NREQNH,NREQNE,NREQNV,NREQNQ)
    endif
    !
    allocate(PHYSA(NR_PHYSA))
    allocate(PHYSAm(NR_PHYSA))
    allocate(NR_COMP(NR_PHYSA))
    allocate(D_TYPE(NR_PHYSA))
    allocate(PHYSAi(NR_PHYSA))
    allocate(PHYSAd(NR_PHYSA))
    allocate(ADRES(NR_PHYSA))
    !
    allocate(NREQNH(NRCASES))
    allocate(NREQNE(NRCASES))
    allocate(NREQNV(NRCASES))
    allocate(NREQNQ(NRCASES))
    !
  end subroutine alloc_physics

  !> Purpose : deallocate data structure for multiphysics
  subroutine dealloc_physics
    if (allocated(PHYSA)) then
       deallocate(PHYSA,PHYSAm,NR_COMP,D_TYPE,PHYSAi,PHYSAd,ADRES, &
                  NREQNH,NREQNE,NREQNV,NREQNQ)
    endif
  end subroutine dealloc_physics
  !
  !> Purpose : routine dumps out multiphysics data structure
  !! @param fp file to dumpout
  subroutine dumpout_physics_to_file(fp)
    implicit none
    character(len=*), intent(in) :: fp
    integer, parameter           :: ndump = 31
    integer                      :: i
    open(unit=ndump,file=fp, &
         form='formatted',access='sequential',status='unknown')
    !
    write(ndump,*) NR_PHYSA
    write(ndump,*) (PHYSA(i),'  ',i=1,NR_PHYSA)
    write(ndump,*) (PHYSAm(i),'  ',i=1,NR_PHYSA)
    write(ndump,*) NR_COMP
    write(ndump,*) (D_TYPE(i),'  ',i=1,NR_PHYSA)
    write(ndump,*) (PHYSAi(i),'  ',i=1,NR_PHYSA)
    write(ndump,*) (PHYSAd(i),'  ',i=1,NR_PHYSA)
    write(ndump,*) ADRES
    write(ndump,*) NRHVAR,NREVAR,NRVVAR,NRQVAR
    write(ndump,*) NRINDEX
    write(ndump,*) NRCASES
    write(ndump,*) NREQNH
    write(ndump,*) NREQNE
    write(ndump,*) NREQNV
    write(ndump,*) NREQNQ
    !
    close(ndump)
  end subroutine dumpout_physics_to_file

  !> Purpose : dumpout to default location "files/dumpPHYS".
  subroutine dumpout_physics_to_default
    call dumpout_physics_to_file('files/dumpPHYS')
  end subroutine dumpout_physics_to_default
  !
  !> Purpose : routine dumps in multiphysics data structure
  !! @param fp file to dumpout
  subroutine dumpin_physics_from_file(fp)
    implicit none
    character(len=*), intent(in) :: fp
    integer, parameter           :: ndump = 31
    open(unit=ndump,file=fp, &
         form='formatted',access='sequential',status='unknown')
    !
    read(ndump,*) NR_PHYSA
    allocate(PHYSA(NR_PHYSA))
    read(ndump,*) PHYSA
    allocate(PHYSAm(NR_PHYSA))
    read(ndump,*) PHYSAm
    allocate(NR_COMP(NR_PHYSA))
    read(ndump,*) NR_COMP
    allocate(D_TYPE(NR_PHYSA))
    read(ndump,*) D_TYPE
    allocate(PHYSAi(NR_PHYSA))
    read(ndump,*) PHYSAi
    allocate(PHYSAd(NR_PHYSA))
    read(ndump,*) PHYSAd
    allocate(ADRES(NR_PHYSA))
    read(ndump,*) ADRES
    read(ndump,*) NRHVAR,NREVAR,NRVVAR,NRQVAR
    read(ndump,*) NRINDEX
    read(ndump,*) NRCASES
    allocate(NREQNH(NRCASES))
    read(ndump,*) NREQNH
    allocate(NREQNE(NRCASES))
    read(ndump,*) NREQNE
    allocate(NREQNV(NRCASES))
    read(ndump,*) NREQNV
    allocate(NREQNQ(NRCASES))
    read(ndump,*) NREQNQ
    !
    close(ndump)
  end subroutine dumpin_physics_from_file

  !> Purpose : dumpin to default location "files/dumpPHYS".
  subroutine dumpin_physics_from_default
    call dumpin_physics_from_file('files/dumpPHYS')
  end subroutine dumpin_physics_from_default

  !> @brief Return discretization type for a physics attribute
  !> @date Mar 2023
  function Discretization_type(Phys)
  character(len=5) :: Phys
  character(len=6) :: Discretization_type
  integer :: n
  call locate_char(Phys,PHYSA,NR_PHYSA, n)
  select case(n)
    case(0); write(*,*) 'Discretization_type: ATTRIBUTE = ',Phys; stop
    case default; Discretization_type = S_DType(D_TYPE(n))
  end select
  end function Discretization_type

  !> @brief Compute the number of interface variables
  !> @date Mar 2023
  subroutine get_num_interf_vars(InrHvar,InrEvar,InrVvar)
  InrHvar=0; InrEvar=0; InrVvar=0
  do i=1,NR_PHYSA
     if (PHYSAi(i)) then
       select case(D_TYPE(i))
          !  .....H^1
       case(CONTIN)
          InrHvar = InrHvar+1
          !  .....H(curl)
       case(TANGEN)
          InrEvar = InrEvar+1
          !  .....H(div)
       case(NORMAL)
          InrVvar = InrVvar+1
       case(DISCON)
          write(*,*) 'get_num_interf_vars: L2 INTERFACE VARIABLE ??'
          stop 1
       end select
     endif
  enddo
  end subroutine get_num_interf_vars

end module physics
