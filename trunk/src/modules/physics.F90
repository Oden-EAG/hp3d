!> Purpose :
!!   module defines multiphysics info
!!   @revision Aug 18
module physics

  ! number of different physics attributes
  integer, save :: NR_PHYSA

  ! list of physics attributes
  character(len=5), save, allocatable :: PHYSA(:)

  ! active/disactive physics attributes flags
  logical, save, allocatable :: PHYSAm(:)

  !  the corresponding number of components
  integer, save, allocatable :: NR_COMP(:)

  !  ...the corresponding discretization type
  character(len=6), save, allocatable :: DTYPE(:)

  !  ...interface variable flag
  logical, save, allocatable :: PHYSAi(:)
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

  !> Purpose : allocate data structure for multiphysics
  subroutine alloc_physics
    !
    if (allocated(PHYSA)) then
       deallocate(PHYSA,PHYSAm,NR_COMP,DTYPE,PHYSAi,ADRES, &
                  NREQNH,NREQNE,NREQNV,NREQNQ)
    endif
    !
    allocate(PHYSA(NR_PHYSA))
    allocate(PHYSAm(NR_PHYSA))
    allocate(NR_COMP(NR_PHYSA))
    allocate(DTYPE(NR_PHYSA))
    allocate(PHYSAi(NR_PHYSA))
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
       deallocate(PHYSA,PHYSAm,NR_COMP,DTYPE,PHYSAi,ADRES, &
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
    write(ndump,*) (DTYPE(i),'  ',i=1,NR_PHYSA)
    write(ndump,*) (PHYSAi(i),'  ',i=1,NR_PHYSA)
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
    integer                      :: i
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
    allocate(DTYPE(NR_PHYSA))
    read(ndump,*) DTYPE
    allocate(PHYSAi(NR_PHYSA))
    read(ndump,*) PHYSAi
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

  !> Purpose : return discretization type for a physics attribute
  function Discretization_type(Phys)
  character(len=5) :: Phys
  character(len=6) :: Discretization_type
  call locate_char(Phys,PHYSA,NR_PHYSA, n)
  select case(n)
  case(0); write(*,*) 'Discretization_type: ATTRIBUTE = ',Phys,' DOES NOT EXIST'; stop 1
  case default; Discretization_type = DTYPE(n)
  end select
  end function Discretization_type

  !> Purpose : compute the number of interface variables
  subroutine get_num_interf_vars(InrHvar,InrEvar,InrVvar)
  InrHvar=0; InrEvar=0; InrVvar=0
  do i=1,NR_PHYSA
     if (PHYSAi(i)) then
       select case(DTYPE(i))
          !  .....H^1
       case('contin')
          InrHvar = InrHvar+1
          !  .....H(curl)
       case('tangen')
          InrEvar = InrEvar+1
          !  .....H(div)
       case('normal')
          InrVvar = InrVvar+1
       case('discon')
          write(*,*) 'get_num_interf_vars: L2 INTERFACE VARIABLE ??'
          stop 1
       end select
     endif
  enddo
  end subroutine get_num_interf_vars




end module physics


