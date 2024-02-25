!> Purpose : read input file and set physics and data structure
subroutine read_input(Fp)
  use data_structure3D , only: MAXNODS
  use environment      , only: QUIET_MODE
  use parameters       , only: MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ,NRRHS
  use physics
  implicit none
  !----------------------------------------------------------------------
  ! input arguments
  character(len=*), intent(in) :: Fp

  ! local variables
  integer      :: narray(10)
  integer      :: i, j
  character(6) :: dtype
  integer, parameter :: nin = 103
  
#if DEBUG_MODE
  integer :: iprint
  iprint=0
#endif
  !----------------------------------------------------------------------
  !!write(*,*) 'QUIET_MODE = ',QUIET_MODE; call pause
  !----------------------------------------------------------------------
  ! file open
  open(unit=nin,file=Fp, &
       form='formatted',access='sequential',status='old',action='read')
  !
  !  ...read in the maximum number of nodes
  if (MAXNODS .gt. 0) then
    read(nin,*) i ! use preset user value for MAXNODS
  else
    read(nin,*) MAXNODS
  endif
  !
  !  ...read number of physics variables
  read(nin,*) NR_PHYSA
  !
  !  ...set number of cases
  NRCASES = 2**NR_PHYSA-1
  !
  !  ...allocate the data structure arrays
  call alloc_physics
  !
  !  ...read in the physics attributes, the corresponding discretization
  !     type and number of components, compute the number of entries in
  !     index and the total number of different components
  NRHVAR=0; NREVAR=0; NRVVAR=0; NRQVAR=0
  NRINDEX=0
  do i=1,NR_PHYSA
     read(nin,*) PHYSA(i),dtype,NR_COMP(i)
     D_TYPE(i) = I_DType(dtype)
!!!     write(*,7001) PHYSA(i),S_DType(D_TYPE(i)),NR_COMP(i)
!!!7001 format('read_input: PHYSICS ATTRIBUTE, D_TYPE, NR OF COMP = ', &
!!!          a5,2x,a6,i3)
  !
  !  default setting for each attribute is `active'
     PHYSAm(i) = .true.
  !
  !  default setting for each attribute is not an interface variable
     PHYSAi(i) = .false.
  !
  !  default setting for each attribute is not a homogeneous Dirichlet variable
     PHYSAd(i) = .false.
  !
     select case(D_TYPE(i))
        !  .....H^1
     case(CONTIN)
        ADRES(i) = NRHVAR
        NRHVAR = NRHVAR + NR_COMP(i)
        !  .....H(curl)
     case(TANGEN)
        ADRES(i) = NREVAR
        NREVAR = NREVAR + NR_COMP(i)
        !  .....H(div)
     case(NORMAL)
        ADRES(i) = NRVVAR
        NRVVAR = NRVVAR + NR_COMP(i)
        !  .....L^2
     case(DISCON)
        ADRES(i) = NRQVAR
        NRQVAR = NRQVAR + NR_COMP(i)
     end select
  enddo
  NRINDEX = NRHVAR + NREVAR + NRVVAR + NRQVAR
  NRINDEX_HEV = NRHVAR + NREVAR + NRVVAR
  !
  if (NRINDEX_HEV > MAX_NRINDEX_HEV) then
    write(*,*) 'NRINDEX_HEV, MAX_NRINDEX_HEV = ',NRINDEX_HEV,MAX_NRINDEX_HEV
    stop
  endif
  !
  !  set parameter module variables used by solelm
  MAXEQNH = NRHVAR * NRRHS
  MAXEQNE = NREVAR * NRRHS
  MAXEQNV = NRVVAR * NRRHS
  MAXEQNQ = NRQVAR * NRRHS
  !
  close(nin)
  !----------------------------------------------------------------------
  ! printing
if (.not. QUIET_MODE) write(*,*) '-- Physics --'
if (.not. QUIET_MODE) write(*,9999) MAXNODS
9999 format(' MAXNODS = ',i12)
  do i=1, NR_PHYSA
if (.not. QUIET_MODE) write(*,1003) PHYSA(i),S_DType(D_TYPE(i)),NR_COMP(i)
1003 format(' PHYSA,D_TYPE,NR_COMP = ',a5,2x,a6,2x,i1)
  enddo
if (.not. QUIET_MODE) write(*,*)''

  !----------------------------------------------------------------------
  !  check the order in which the attributes have been defined
  do i=1,NR_PHYSA-1
     if (D_TYPE(i+1).lt.D_TYPE(i)) then
        write(*,6001); stop 1
     endif
  enddo

6001 format('read_input: WRONG ORDER OF PHYSICS ATTRIBUTES')
  !----------------------------------------------------------------------
  !  ...generate info on number of H1, H(curl), H(div) and L2 components
  !     for each case
  do i=1,NRCASES
     NREQNH(i) = 0; NREQNE(i) = 0; NREQNV(i) = 0; NREQNQ(i) = 0
     narray = 0
     call decod(i,2,NR_PHYSA, narray)
     do j=1,NR_PHYSA
        if (narray(j).eq.1) then
           select case(D_TYPE(j))
           case(CONTIN)
              NREQNH(i) = NREQNH(i) + NR_COMP(j)
           case(TANGEN)
              NREQNE(i) = NREQNE(i) + NR_COMP(j)
           case(NORMAL)
              NREQNV(i) = NREQNV(i) + NR_COMP(j)
           case(DISCON)
              NREQNQ(i) = NREQNQ(i) + NR_COMP(j)
           end select
        endif
     enddo
#if DEBUG_MODE
     if (iprint.eq.1) then
       write(*,9998)i,NREQNH(i),NREQNE(i),NREQNV(i),NREQNQ(i)
9998   format(' Icase,NREQNH,NREQNE,NREQNV,NREQNQ = ',5(i2,2x))
     endif
#endif
  enddo

end subroutine read_input
