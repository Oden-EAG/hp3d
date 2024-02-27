!
#include "typedefs.h"
!
#if HP3D_USE_INTEL_MKL
!
!----------------------------------------------------------------------
!
!   module name        - pardiso_data
!
!----------------------------------------------------------------------
!
!   latest revision    - May 2018
!
!   purpose            - module stores data for PARDISO coarse solve
!
!----------------------------------------------------------------------
!
      module pardiso_data
!
      implicit none
!
!  ...problem workspace
      integer                 :: PRDS_N,PRDS_NZ,PRDS_NRHS
      integer,    allocatable :: PRDS_IA(:), PRDS_JA(:)
      VTYPE,      allocatable :: PRDS_A(:),PRDS_RHS(:),PRDS_XSOL(:)
!
!  ...pardiso workspace
!  ...internal solver memory pointer for 64-bit architectures
      integer(8)              :: PRDS_PT(64)
      integer                 :: PRDS_IPARM(64)
      integer                 :: PRDS_MTYPE, PRDS_MAXFCT, PRDS_MNUM
      integer                 :: PRDS_PHASE, PRDS_ERROR, PRDS_MSGLVL
      integer, allocatable    :: PRDS_PERM(:)
      character(1)            :: PRDS_TYPE
!
      contains
!
!----------------------------------------------------------------------
!
      subroutine start_pardiso
!
      implicit none
!
#if HP3D_COMPLEX
      select case(PRDS_TYPE)
      case('S')
         PRDS_MTYPE   = 6  ! complex symmetric
      case('H')
         PRDS_MTYPE   = 4  ! complex hermitian positive definite
      case('I')
         PRDS_MTYPE   = -4 ! complex hermitian indefinite
      case default
         PRDS_MTYPE   = 13  ! complex non-symmetric
      end select
#else
      select case(PRDS_TYPE)
      case('S')
         PRDS_MTYPE   = -2  ! real symmetric indefinite
      case('H')
         PRDS_MTYPE   = 2  ! real symmetric positive definite
      case default
         PRDS_MTYPE   = 11  ! real non-symmetric
      end select
#endif
      PRDS_IPARM(8) = 0
!
      call pardisoinit(PRDS_PT, PRDS_MTYPE, PRDS_IPARM)
!
      PRDS_MNUM    = 1
      PRDS_MAXFCT  = 1
      PRDS_MSGLVL  = 0       ! with statistical no information
!
      end subroutine start_pardiso
!
!----------------------------------------------------------------------
!
      subroutine finalize_pardiso
!
      implicit none
!
      PRDS_PHASE     =  0   ! release internal memory
      call pardiso(PRDS_PT,PRDS_MAXFCT,PRDS_MNUM,PRDS_MTYPE,PRDS_PHASE,    &
                   PRDS_N,PRDS_A,PRDS_IA,PRDS_JA,PRDS_PERM,PRDS_NRHS,      &
                   PRDS_IPARM,PRDS_MSGLVL,PRDS_RHS,PRDS_XSOL,PRDS_ERROR)
!
      PRDS_PHASE     = -1    ! release internal memory
      call pardiso(PRDS_PT,PRDS_MAXFCT,PRDS_MNUM,PRDS_MTYPE,PRDS_PHASE,    &
                   PRDS_N,PRDS_A,PRDS_IA,PRDS_JA,PRDS_PERM,PRDS_NRHS,      &
                   PRDS_IPARM,PRDS_MSGLVL,PRDS_RHS,PRDS_XSOL,PRDS_ERROR)
!
      deallocate(PRDS_A, PRDS_IA,PRDS_JA,PRDS_RHS,PRDS_XSOL,PRDS_PERM)
!
      end subroutine finalize_pardiso
!
!----------------------------------------------------------------------
!
      end module pardiso_data
!

#endif
