!
#include "typedefs.h"
!
!----------------------------------------------------------------------
!> @brief       Wrapper for element stiffness matrix assembly
!!
!> @param[in]   mdle    - middle node number
!> @param[out]  Itest   - index for assembly
!!                         (trivial since not multiphysics application)
!> @param[out]  Itrial  - index for assembly
!!                         (trivial since not multiphysics application)
!!
!> @date        July 2023
!----------------------------------------------------------------------
   subroutine elem(Mdle, Itest,Itrial)
!
      use assembly
      use commonParam
      use data_structure3D
      use parametersDPG
      use physics, only: NR_PHYSA
!
      implicit none
!
!  ...declare input/output variables
      integer, intent(in)  :: Mdle
      integer, intent(out) :: Itest(NR_PHYSA)
      integer, intent(out) :: Itrial(NR_PHYSA)
!
!  ...number of test and trial degrees of freedom
      integer :: nrdofH,   nrdofE,   nrdofV,   nrdofQ
      integer :: nrdofHH,  nrdofEE,  nrdofVV,  nrdofQQ
      integer :: ndofHmdl, ndofEmdl, ndofVmdl, ndofQmdl
      integer :: nrdofEi, nrdofVi
      integer :: nrTest,  nrTrial
!
      integer :: nrv, nre, nrf
!
!  ...element order and enriched order
      integer :: norder(19), norderP(19), nordP
!
!  ...middle node type
      integer :: ntype
!
!  ...timer
      real(8) :: start_time, end_time
!
!----------------------------------------------------------------------
!
      Itest (1:NR_PHYSA) = 0
      Itrial(1:NR_PHYSA) = 0
      norder (1:19) = 0
      norderP(1:19) = 0
!
      ntype = NODES(Mdle)%ntype
      nrv = nvert(ntype); nre = nedge(ntype); nrf = nface(ntype)
!
!  ...determine order of approximation
      call find_order(Mdle, norder)
!
!  ...set the enriched order of approximation
      select case(ntype)
         case(MDLB)
            nordP = NODES(Mdle)%order + NORD_ADD*111
         case(MDLP)
            nordP = NODES(Mdle)%order + NORD_ADD*11
         case(MDLN,MDLD)
            nordP = NODES(Mdle)%order + NORD_ADD
         case default
            write(*,*) 'elem: invalid ntype param. stop.'
            stop
      end select
!
!  ...note: compute_enriched_order works only for hexa and prism currently
      call compute_enriched_order(ntype,nordP, norderP)
!  ...compute nrdof for trial
      call celndof(ntype,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
!  ...compute nrdof for test
      call celndof(ntype,norderP, nrdofHH,nrdofEE,nrdofVV,nrdofQQ)
!  ...compute number of bubble dofs
      call ndof_nod(ntype,norder(nre+nrf+1), ndofHmdl,ndofEmdl,ndofVmdl,ndofQmdl)
!
!  ...node supports all physical attributes
!  ...2 physical attributes: case = 2^2-1 = 3
      if (NODES(Mdle)%case .ne. 3) then
         write(*,*) 'elem: Mdle,NODES(Mdle)%case = ', &
              Mdle,NODES(Mdle)%case, '. stop.'
         call logic_error(ERR_INVALID_VALUE, __FILE__,__LINE__)
         stop
      endif
!
!  ...calculate number of interface dofs
      nrdofEi = nrdofE-ndofEmdl
!
!  ...calculate number of trial and test dofs
      nrTest  = 2*nrdofEE
      nrTrial = 2*nrdofEi + 6*nrdofQ
!
!  ...signal case
      Itest(1:2)=1; Itrial(1:2)=1
!
      call elem_maxwell(Mdle,nrTest,nrTrial,                   &
         nrdofEE,nrdofH,nrdofE,nrdofQ,nrdofEi,                 &
         BLOC(1)%nrow,BLOC(2)%nrow,                            &
         BLOC(1)%array,ALOC(1,1)%array,ALOC(1,2)%array,        &
         BLOC(2)%array,ALOC(2,1)%array,ALOC(2,2)%array)
!
   end subroutine elem







!----------------------------------------------------------------------
!> @brief       Compute enriched order for nodes of element
!!
!> @param[in]   ntype   - middle node type
!> @param[in]   Nord    - scalar element order
!> @param[out]  Norder  - vector element order
!!
!> @date        July 2023
!----------------------------------------------------------------------
   subroutine compute_enriched_order(ntype,Nord, Norder)
!
      use parameters, only : MODORDER
      use node_types
!
      implicit none
!
      integer, intent(in)  :: ntype
      integer, intent(in)  :: Nord
      integer, intent(out) :: Norder(19)
!
      integer :: temp(2)
      integer :: nordF(3),nordB(3)
!
!----------------------------------------------------------------------
!
!  ...see implementation of BrokenExactSequence in shape functions
      select case(ntype)
         case(MDLB)
            call decod(Nord,MODORDER,2, temp) !xy face, z edge
            nordF(1) = temp(1); nordB(3) = temp(2)
            call decod(nordF(1),MODORDER,2, nordB(1:2)) !x edge, y edge
            call encod((/nordB(1),nordB(3)/),MODORDER,2, nordF(2)) !xz face
            call encod(nordB(2:3),MODORDER,2, nordF(3)) !yz face
!        ...edges
            Norder(1:4)   = (/nordB(1),nordB(2),nordB(1),nordB(2)/) !x,y,x,y edges
            Norder(5:8)   = (/nordB(1),nordB(2),nordB(1),nordB(2)/) !x,y,x,y edges
            Norder(9:12)  = nordB(3) !z edges
!        ...faces
            Norder(13:14) = nordF(1) !xy faces
            Norder(15:18) = (/nordF(2),nordF(3),nordF(2),nordF(3)/) !xz,yz,xz,yz faces
!        ...element interior
            Norder(19)    = Nord
         case(MDLP)
            call decod(Nord,MODORDER,2, nordB(1:2))
!        ...edges (xy)
            norder(1:6)   = nordB(1)
!        ...edges (z)
            norder(7:9)   = nordB(2)
!        ...triangular faces
            norder(10:11) = nordB(1)
!        ...quadrilateral faces (z)
            norder(12:14) = Nord
!        ...element interior
            norder(15)    = Nord
         case(MDLN,MDLD)
            norder(1:11) = Nord
         case default
            write(*,*) 'compute_enriched_order: incorrect element type. stop.'
            stop
      end select
!
   end subroutine compute_enriched_order
!
