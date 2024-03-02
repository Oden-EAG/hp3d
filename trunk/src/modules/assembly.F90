!> @brief define matrices for FE assembly
module assembly
   use parameters
   use physics
   use error
   implicit none

#include "typedefs.h"

! Type definition
!~~~~~~~~~~~~~~~~
   type super_vector
      integer :: nrow
      VTYPE, allocatable :: vect(:)
   end type super_vector
   !
   type super_array
      integer :: ncol, nrow
      VTYPE, allocatable :: array(:,:)
   end type super_array

   save

! Member variables
!~~~~~~~~~~~~~~~~~
!..number of right-hand sides (load vectors)
   integer :: NR_RHS

!..maximum number of element local dof for a physics attribute
   integer, allocatable :: MAXDOFS(:)

!..maximum number of modified element dof for all variables and compression mode
   integer :: MAXDOFM, MAXDOFC

!..element local load vectors stiffness matrices
   type(super_array), allocatable :: BLOC(:), ALOC(:,:), AAUX(:)
!$OMP THREADPRIVATE (BLOC, ALOC, AAUX)

!..modified element load vector and stiffness matrix
   VTYPE, allocatable :: ZBMOD(:,:), ZAMOD(:,:)
!$OMP THREADPRIVATE (ZBMOD, ZAMOD)

!..dof extraction vector and dirichlet data flag
   integer, allocatable :: NEXTRACT(:), IDBC(:)
!$OMP THREADPRIVATE (NEXTRACT, IDBC)

!..Dirichlet data
   VTYPE, allocatable :: ZDOFD(:,:)
!$OMP THREADPRIVATE (ZDOFD)

contains

!> @brief create workspace for the celem call
   subroutine assembly_begin
      integer :: istat

!  ...set the possible maximum and allocate workspace
      MAXDOFM = &
         MAXbrickH*NRHVAR + MAXbrickE*NREVAR + &
         MAXbrickV*NRVVAR + MAXbrickQ*NRQVAR
      allocate( &
         NEXTRACT(MAXDOFM), IDBC(MAXDOFM), &
         ZDOFD(MAXDOFM,NR_RHS), MAXDOFS(NR_PHYSA), &
         stat=istat)

      if (istat.ne.SUCCESS) then
         call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
      endif

   end subroutine assembly_begin

   subroutine assembly_begin_par


!  ...set the possible maximum and allocate workspace
      MAXDOFM = &
         MAXbrickH*NRHVAR + MAXbrickE*NREVAR + &
         MAXbrickV*NRVVAR + MAXbrickQ*NRQVAR
      allocate( &
         NEXTRACT(MAXDOFM), IDBC(MAXDOFM), &
         ZDOFD(MAXDOFM,NR_RHS))

   end subroutine assembly_begin_par

!> @brief free workspace
   subroutine assembly_end
      integer :: istat
      deallocate( &
         NEXTRACT, IDBC, ZDOFD, MAXDOFS, &
         stat=istat)

      if (istat.ne.SUCCESS) then
         call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
      endif

   end subroutine assembly_end

   subroutine assembly_end_par
      deallocate(NEXTRACT, IDBC, ZDOFD)
   end subroutine assembly_end_par

!> @brief allocate arrays
   subroutine assembly_alloc
      integer :: i, j, istat
      allocate( &
         BLOC(NR_PHYSA), AAUX(NR_PHYSA), ALOC(NR_PHYSA,NR_PHYSA), &
         stat=istat)

      if (istat.ne.SUCCESS) then
         call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
      endif

      do i=1,NR_PHYSA
         BLOC(i)%nrow = MAXDOFS(i)
         BLOC(i)%ncol = NR_RHS
         allocate(BLOC(i)%array(MAXDOFS(i),NR_RHS), stat=istat)

         if (istat.ne.SUCCESS) then
            call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
         endif

         do j=1,NR_PHYSA
            ALOC(i,j)%nrow = MAXDOFS(i)
            ALOC(i,j)%ncol = MAXDOFS(j)
            allocate(ALOC(i,j)%array(MAXDOFS(i),MAXDOFS(j)), stat=istat)

            if (istat.ne.SUCCESS) then
               call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
            endif

         enddo
         AAUX(i)%nrow = MAXDOFM
         AAUX(i)%ncol = MAXDOFS(i)
         allocate(AAUX(i)%array(MAXDOFM,MAXDOFS(i)), stat=istat)

         if (istat.ne.SUCCESS) then
            call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
         endif
      enddo
      allocate(ZBMOD(MAXDOFM,NR_RHS), ZAMOD(MAXDOFM,MAXDOFM), stat=istat)

      if (istat.ne.SUCCESS) then
         call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
      endif

   end subroutine assembly_alloc

!> @brief deallocate arrays
   subroutine assembly_dealloc
      integer :: i, j, istat
      do i=1,NR_PHYSA
         if (allocated(BLOC(i)%array)) deallocate(BLOC(i)%array, stat=istat)
         if (istat.ne.SUCCESS) then
            call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
         endif

         do j=1,NR_PHYSA
            if (allocated(ALOC(i,j)%array)) deallocate(ALOC(i,j)%array, stat=istat)
            if (istat.ne.SUCCESS) then
               call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
            endif
         enddo

         if (allocated(AAUX(i)%array)) deallocate(AAUX(i)%array, stat=istat)
         if (istat.ne.SUCCESS) then
            call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
         endif
      enddo

      deallocate(BLOC,AAUX,ALOC,ZBMOD,ZAMOD, stat=istat)
      if (istat.ne.SUCCESS) then
         call logic_error(ERR_ALLOC_FAILURE,__FILE__,__LINE__)
      endif

   end subroutine assembly_dealloc

end module assembly
