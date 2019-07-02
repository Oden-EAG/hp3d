!=====================================================================
module SPLIT_SURF                                                    !
!---------------------------------------------------------------------
  use kinds                                                          !
  use GMP                                                            !
  use control                                                        !
!---------------------------------------------------------------------
! ** SHARED VARIABLES **                                             | 
!---------------------------------------------------------------------
  IMPLICIT NONE                                                      !
  SAVE                                                               !
! splitting surface number                                           |
  integer                              :: Nsplit_MOD                 !
! number of bounding surfaces                                        |
  integer                              :: Nr_bound_MOD               !
! number of conforming surfaces                                      |
  integer                              :: Nr_confm_MOD               !
! bounding surfaces numbers                                          |
  integer, dimension(:), allocatable   :: Ns_bound_MOD               !
! conforming sufaces numbers                                         |
  integer, dimension(:), allocatable   :: Ns_confm_MOD               !
! shifts                                                             |
  real                                 :: Dh1_MOD,Dh2_MOD            !
! original GMP parameters                                            |
  integer :: nrcurve_orig, nrtrian_orig, nrrecta_orig, &             !
             nrprism_orig, nrtetra_orig, nrpyram_orig                !
! number of figures to split                                         |
  integer                              :: nr_figs_to_split           !
! info about new figs, points, curves                                |
  integer, dimension(:,:), allocatable :: figs_split                 !
  integer, dimension(:),   allocatable :: new_point                  !
  integer, dimension(:,:), allocatable :: new_curves                 !
! normal to split plane                                              |
  real, dimension(3)                   :: NORMAL                     !
! anticipated maximum number of new curves                           |
  integer, parameter                   :: MAX_NEW_CURVES = 1000      !
  contains                                                           !
!---------------------------------------------------------------------
! ** SHARED SUBROUTINES **                                           |
!---------------------------------------------------------------------
! allocate_Ns_bound_MOD                                              |
! deallocate_Ns_bound_MOD                                            |
! allocate_Ns_confm_MOD                                              |
! deallocate_Ns_confm_MOD                                            |
! deallocate_FIGS_SPLIT                                              |
! deallocate_NEW_POINT                                               |
! allocate_NEW_CURVES                                                |
! deallocate_NEW_CURVES                                              |
!=====================================================================
!
!
!
!-------------------------------------------------------------------------------
subroutine allocate_Ns_bound_MOD
!-------------------------------------------------------------------------------
  integer :: status
!-------------------------------------------------------------------------------
! printing flag (0,1)
#define I_PRINT 0
!
#if I_PRINT >= 1
    write(*,*)'allocate_Ns_bound_MOD: allocating Ns_bound_MOD...'
#endif    
    allocate(Ns_bound_MOD(Nr_bound_MOD), STAT = status)
    if (status .ne. 0) then
      write(*,*)'allocate_Ns_bound_MOD: Ns_bound_MOD not allocated.'
      stop
    endif
#if I_PRINT >= 1
    write(*,*)'allocate_Ns_bound_MOD: Ns_bound_MOD allocated successfully.'
#endif    
end subroutine allocate_Ns_bound_MOD
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
subroutine deallocate_Ns_bound_MOD
!-------------------------------------------------------------------------------
  integer :: status
!-------------------------------------------------------------------------------
! printing flag (0,1)
#define I_PRINT 0
!
#if I_PRINT >= 1
    write(*,*)'deallocate_Ns_bound_MOD: deallocating Ns_bound_MOD...'
#endif    
    deallocate(Ns_bound_MOD, STAT = status)
    if (status .ne. 0) then
      write(*,*)'deallocate_Ns_bound_MOD: Ns_bound_MOD not deallocated.'
      stop
    endif
#if I_PRINT >= 1
    write(*,*)'deallocate_Ns_bound_MOD: Ns_bound_MOD deallocated successfully.'
#endif    
end subroutine deallocate_Ns_bound_MOD
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
subroutine allocate_Ns_confm_MOD
!-------------------------------------------------------------------------------
  integer :: status
!-------------------------------------------------------------------------------
! printing flag (0,1)
#define I_PRINT 0
!
#if I_PRINT >= 1
    write(*,*)'allocate_Ns_confm_MOD: allocating Ns_confm_MOD...'
#endif    
    allocate(Ns_confm_MOD(Nr_confm_MOD), STAT = status)
    if (status .ne. 0) then
      write(*,*)'allocate_Ns_confm_MOD: Ns_confm_MOD not allocated.'
      stop
    endif
#if I_PRINT >= 1
    write(*,*)'allocate_Ns_confm_MOD: Ns_confm_MOD allocated successfully.'
#endif    
end subroutine allocate_Ns_confm_MOD
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
subroutine deallocate_Ns_confm_MOD
!-------------------------------------------------------------------------------
  integer :: status
!-------------------------------------------------------------------------------
! printing flag (0,1)
#define I_PRINT 0
!
#if I_PRINT >= 1
    write(*,*)'deallocate_Ns_confm_MOD: deallocating Ns_confm_MOD...'
#endif    
    deallocate(Ns_confm_MOD, STAT = status)
    if (status .ne. 0) then
      write(*,*)'deallocate_Ns_confm_MOD: Ns_confm_MOD not deallocated.'
      stop
    endif
#if I_PRINT >= 1
    write(*,*)'deallocate_Ns_confm_MOD: Ns_confm_MOD deallocated successfully.'
#endif    
end subroutine deallocate_Ns_confm_MOD
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
subroutine deallocate_FIGS_SPLIT
!-------------------------------------------------------------------------------
  integer :: status
!-------------------------------------------------------------------------------
! printing flag (0,1)  
#define I_PRINT 0
!
    deallocate(figs_split, STAT = status)
    if (status .ne. 0) then
      write(*,*)'deallocate_FIGS_SPLIT: figs_split not deallocated.'
      stop
    endif
#if I_PRINT >= 1
    write(*,*)'deallocate_FIGS_SPLIT: figs_split deallocated successfully.'
#endif    
end subroutine deallocate_FIGS_SPLIT
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
subroutine deallocate_NEW_POINT
!-------------------------------------------------------------------------------
  integer :: status
!-------------------------------------------------------------------------------
! printing flag (0,1)  
#define I_PRINT 0
! 
    deallocate(new_point, STAT = status)
    if (status .ne. 0) then
      write(*,*)'deallocate_NEW_POINT: new_point not deallocated.'
      stop
    endif
#if I_PRINT >= 1
    write(*,*)'deallocate_NEW_POINT: new_point deallocated successfully.'
#endif    
end subroutine deallocate_NEW_POINT
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
subroutine allocate_NEW_CURVES
!-------------------------------------------------------------------------------
  integer :: status
!-------------------------------------------------------------------------------
! printing flag (0,1)  
#define I_PRINT 0
!
    allocate(new_curves(2,MAX_NEW_CURVES), STAT = status)
    if (status .ne. 0 ) then
      write(*,*)'allocate_NEW_CURVES: new_curves not allocated.'
      stop
    endif
#if I_PRINT >= 1
    write(*,*)'allocate_NEW_CURVES: new_curves allocate successfully.'    
#endif
end subroutine allocate_NEW_CURVES
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
subroutine deallocate_NEW_CURVES
!-------------------------------------------------------------------------------
  integer :: status
!-------------------------------------------------------------------------------
! printing flag (0,1)  
#define I_PRINT 0
!
    deallocate(new_curves, STAT = status)
    if (status .ne. 0) then
      write(*,*)'deallocate_NEW: new_curves not deallocated.'
      stop
    endif
#if I_PRINT >= 1
    write(*,*)'deallocate_NEW: new_curves deallocated successfully.'
#endif        
end subroutine deallocate_NEW_CURVES
!-------------------------------------------------------------------------------
!
!
!===============================================================================
end module SPLIT_SURF                                                          !
!===============================================================================
