#include "typedefs.h"
!----------------------------------------------------------------------
!> @brief   module stores information about h-refined elements
!> @date    Oct 2023
!----------------------------------------------------------------------
!
   module refinements_history
!
      use node_types
      use parameters, only: NDIMEN,MAXEQNH,MAXEQNQ,MAXEQNV,MAXEQNQ, &
                            MAXbrickH,MAXbrickE,MAXbrickV,MAXbrickQ
!
      implicit none
!
!  ...anticipated max number of refined elements
      integer :: MAX_ELEMS_REF
!
!  ...number of refined elements
      integer :: NR_ELEMS_REF
!
!----------------------------------------------------------------------
!  REFINED ELEMENT                                                    |
!----------------------------------------------------------------------
      type refined_element
!
!  ......element middle node
         integer :: mdle
!
!  ......element dof
         real*8, dimension(:,:), pointer :: xnod
         VTYPE,  dimension(:,:), pointer :: zdofH
         VTYPE,  dimension(:,:), pointer :: zdofE
         VTYPE,  dimension(:,:), pointer :: zdofV
         VTYPE,  dimension(:,:), pointer :: zdofQ
      endtype refined_element
!
!-----------------------------------------------------------------------
!
!  ...data structure arrays
      type(refined_element), allocatable, save :: ELEMS_REF(:)
!
!-----------------------------------------------------------------------
!
   contains
!
!-----------------------------------------------------------------------
!
!..allocate memory for data structure
   subroutine allocref
!
      integer :: nel
!
      if (allocated(ELEMS_REF)) then
         write(*,*) 'allocref: WARNING !! ELEMS_REFINED', &
                    ' HAS NOT BEEN DEALLOCATED'
         call deallocref
      endif
!
      allocate(ELEMS_REF(MAX_ELEMS_REF))
      do nel=1,MAX_ELEMS_REF
         ELEMS_REF(nel)%mdle = -1
         nullify (ELEMS_REF(nel)%xnod)
         nullify (ELEMS_REF(nel)%zdofH)
         nullify (ELEMS_REF(nel)%zdofE)
         nullify (ELEMS_REF(nel)%zdofV)
         nullify (ELEMS_REF(nel)%zdofQ)
      enddo
!
      NR_ELEMS_REF = 0
!
   end subroutine allocref
!
!-----------------------------------------------------------------------
!
!..deallocate ELEMS_REF
   subroutine deallocref
!
      integer :: nel
!
      do nel=1,NR_ELEMS_REF
        if (associated(ELEMS_REF(nel)%xnod))  deallocate(ELEMS_REF(nel)%xnod)
        if (associated(ELEMS_REF(nel)%zdofH)) deallocate(ELEMS_REF(nel)%zdofH)
        if (associated(ELEMS_REF(nel)%zdofE)) deallocate(ELEMS_REF(nel)%zdofE)
        if (associated(ELEMS_REF(nel)%zdofV)) deallocate(ELEMS_REF(nel)%zdofV)
        if (associated(ELEMS_REF(nel)%zdofQ)) deallocate(ELEMS_REF(nel)%zdofQ)
      enddo
      deallocate(ELEMS_REF)
      NR_ELEMS_REF=0
!
!
   end subroutine deallocref
!
!-----------------------------------------------------------------------
!
!..increase MAX_ELEMS_REF
   subroutine increase_MAXELEMS_REF()
!
      type(refined_element), allocatable :: ELEMS_REF_NEW(:)
      integer :: max_ELEMS_REF_new,nel
!
      if (.not. allocated(ELEMS_REF)) then
         write(*,*) 'increase_MAXELEMS_REF: ELEMS_REF not allocated. returning...'
         return
      endif
!
!  ...determine size of new ELEMS_REF array
      max_ELEMS_REF_new = 2*MAX_ELEMS_REF
!
!  ...allocate new ELEMS_REF array twice the size of the old
      allocate(ELEMS_REF_NEW(max_ELEMS_REF_new))
!
!  ...copy data from old ELEMS_REF array into the new one
      ELEMS_REF_NEW(1:MAX_ELEMS_REF) = ELEMS_REF(1:MAX_ELEMS_REF)
!
      !$OMP PARALLEL DO
      do nel=MAX_ELEMS_REF+1,max_ELEMS_REF_new
        ELEMS_REF_NEW(nel)%mdle = -1
        nullify (ELEMS_REF_NEW(nel)%xnod)
        nullify (ELEMS_REF_NEW(nel)%zdofH)
        nullify (ELEMS_REF_NEW(nel)%zdofE)
        nullify (ELEMS_REF_NEW(nel)%zdofV)
        nullify (ELEMS_REF_NEW(nel)%zdofQ)
      enddo
      !$OMP END PARALLEL DO
!
!  ...move ELEMS_REF pointer to the new array,
!     and deallocate the old array
      call move_alloc(ELEMS_REF_NEW, ELEMS_REF)
      MAX_ELEMS_REF = max_ELEMS_REF_new
!
      end subroutine increase_MAXELEMS_REF
!
      subroutine save_element(Mdle,Xnod,ZdofH,ZdofE,ZdofV,ZdofQ)
      use data_structure3D
!
      integer :: Mdle
      real*8, dimension(NDIMEN,MAXbrickH)  ::  Xnod
      VTYPE,  dimension(MAXEQNH,MAXbrickH) ::  ZdofH
      VTYPE,  dimension(MAXEQNE,MAXbrickE) ::  ZdofE
      VTYPE,  dimension(MAXEQNV,MAXbrickV) ::  ZdofV
      VTYPE,  dimension(MAXEQNQ,MAXbrickQ) ::  ZdofQ
!
      integer :: norder(19),nel,nrdofH,nrdofE,nrdofV,nrdofQ
!
      if (.not.allocated(ELEMS_REF)) then
        if (MAX_ELEMS_REF.le.0) then
          write(*,*) 'break: SET MAX_ELEMS_REF'
          stop 1
        endif
        allocate(ELEMS_REF(MAX_ELEMS_REF))
      endif
!
      if (NR_ELEMS_REF.eq.MAX_ELEMS_REF) call increase_MAXELEMS_REF
!
      nel = NR_ELEMS_REF + 1
      NR_ELEMS_REF = nel
!
      ELEMS_REF(nel)%mdle = Mdle
      call find_order(Mdle, norder)
      call celndof(NODES(Mdle)%ntype,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
!
      if (nrdofH.gt.0) then
        allocate(ELEMS_REF(nel)%xnod (NDIMEN, nrdofH))
        allocate(ELEMS_REF(nel)%zdofH(MAXEQNH,nrdofH))
        ELEMS_REF(nel)%xnod (1:NDIMEN, 1:nrdofH) = Xnod (1:NDIMEN, 1:nrdofH)
        ELEMS_REF(nel)%zdofH(1:MAXEQNH,1:nrdofH) = ZdofH(1:MAXEQNH,1:nrdofH)
      endif
      if (nrdofE.gt.0) then
        allocate(ELEMS_REF(nel)%zdofE(MAXEQNE,nrdofE))
        ELEMS_REF(nel)%zdofE(1:MAXEQNE,1:nrdofE) = ZdofE(1:MAXEQNE,1:nrdofE)
      endif
      if (nrdofV.gt.0) then
        allocate(ELEMS_REF(nel)%zdofV(MAXEQNE,nrdofV))
        ELEMS_REF(nel)%zdofV(1:MAXEQNV,1:nrdofV) = ZdofV(1:MAXEQNV,1:nrdofV)
      endif
      if (nrdofQ.gt.0) then
        allocate(ELEMS_REF(nel)%zdofQ(MAXEQNE,nrdofQ))
        ELEMS_REF(nel)%zdofQ(1:MAXEQNQ,1:nrdofQ) = ZdofQ(1:MAXEQNQ,1:nrdofQ)
      endif
!
   end subroutine save_element
!
!
   end module refinements_history
