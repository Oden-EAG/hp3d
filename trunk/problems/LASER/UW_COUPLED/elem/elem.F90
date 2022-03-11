!--------------------------------------------------------------------
!
!     routine name      - elem
!
!--------------------------------------------------------------------
!
!     latest revision:  - Oct 2019
!
!     purpose:          - routine returns unconstrained (ordinary)
!                         stiffness matrix and load vector
!                         for either:
!                         1) single step of transient heat problem
!                         2) multiple steps of transient heat problem
!                         3) steady state Maxwell problem (signal)
!                         4) steady state Maxwell problem (pump)
!
!     arguments:
!        in:
!              Mdle     - element middle node number
!
!        out:
!              Itest    - index for assembly
!              Itrial   - index for assembly
!
!
!-----------------------------------------------------------------------
!
#include "typedefs.h"
!
subroutine elem(Mdle, Itest,Itrial)
!
   use assembly
   use commonParam
   use data_structure3D
   use parametersDPG
   use physics   , only: NR_PHYSA
   use mpi_param , only: RANK,ROOT
   use laserParam, only: PLANE_PUMP
!
   implicit none
!
!----------------------------------------------------------------------
!
!..declare input/output variables
   integer,                    intent(in)  :: Mdle
   integer,dimension(NR_PHYSA),intent(out) :: Itest,Itrial
!
!..number of test and trial degrees of freedom
   integer :: nrdofH ,nrdofE ,nrdofV ,nrdofQ
   integer :: nrdofHH,nrdofEE,nrdofVV,nrdofQQ
   integer :: nrdofEi,nrdofVi,nrTest,nrTrial
!
   integer :: ndofHmdl,ndofEmdl,ndofVmdl,ndofQmdl
   integer :: nrv,nre,nrf
!
!..element order and enriched order
   integer :: norder(19),norderP(19),nordP
!
!..element type
   character(len=4) :: etype
!
!..fld_flag refers to either pump (0) or signal (1) field
   integer :: fld_flag
!
!..timer
   real(8) :: MPI_Wtime,start_time,end_time
!
!----------------------------------------------------------------------
!
!..start timer
!   start_time = MPI_Wtime()
!
   Itest (1:NR_PHYSA) = 0
   Itrial(1:NR_PHYSA) = 0
   norder (1:19) = 0
   norderP(1:19) = 0
!
   etype = NODES(Mdle)%type
   nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
!..determine order of approximation
   call find_order(Mdle, norder)
!
!..set the enriched order of approximation
   select case(etype)
      case('mdlb')
         nordP = NODES(Mdle)%order+NORD_ADD*111
      case('mdln','mdld')
         nordP = NODES(Mdle)%order+NORD_ADD
      case('mdlp')
         nordP = NODES(Mdle)%order+NORD_ADD*11
      case default
         write(*,*) 'elem: invalid etype param. stop.'
         stop
   end select
!..note: compute_enriched_order works only for hexa and prism currently
   call compute_enriched_order(etype,nordP, norderP)
!..compute nrdof for trial
   call celndof(etype,norder, nrdofH,nrdofE,nrdofV,nrdofQ)
!..compute nrdof for test
   call celndof(etype,norderP, nrdofHH,nrdofEE,nrdofVV,nrdofQQ)
!..compute number of bubble dofs
   call ndof_nod(etype,norder(nre+nrf+1), ndofHmdl,ndofEmdl,ndofVmdl,ndofQmdl)
!
!..node supports all physical attributes
!..6 physical attributes: case = 2^6-1 = 63
   if (NODES(Mdle)%case .ne. 63) then
      write(*,*) 'elem: Mdle,NODES(Mdle)%case = ', &
           Mdle,NODES(Mdle)%case, '. stop.'
      call logic_error(ERR_INVALID_VALUE, __FILE__,__LINE__)
      stop
   endif
!
!..NO_PROBLEM: (1) - single step of transient heat equation
!  ............(2) - multiple steps of transient heat equation
!  ............(3) - time harmonic Maxwell - for signal EH-field
!  ............(4) - time harmonic Maxwell - for pump EH-field
   select case(NO_PROBLEM)
!  ...Heat cases
      case(1,2)
!     ...calculate number of interface dofs
         nrdofVi = nrdofV-ndofVmdl
!     ...calculate number of trial and test dofs
         nrTest  = nrdofHH
         nrTrial = nrdofH + nrdofVi
         Itest(1)=1; Itrial(1)=1
         Itest(4)=1; Itrial(4)=1
         call elem_heat(Mdle,nrTest,nrTrial,                            &
            nrdofHH,nrdofH,nrdofV,nrdofVi,                              &
            BLOC(1)%nrow,BLOC(4)%nrow,                                  &
            BLOC(1)%array,ALOC(1,1)%array,ALOC(1,4)%array,              &
            BLOC(4)%array,ALOC(4,1)%array,ALOC(4,4)%array)
!
!  ...Maxwell cases - signal or pump EH-field
      case(3,4)
!     ...calculate number of interface dofs
         nrdofEi = nrdofE-ndofEmdl
!     ...calculate number of trial and test dofs
         nrTest  = 2*nrdofEE
         nrTrial = 2*nrdofEi + 6*nrdofQ
!
         if (NO_PROBLEM .eq. 3) then
!        ...signal case
            Itest(2)=1; Itrial(2)=1
            Itest(5)=1; Itrial(5)=1
            fld_flag = 1;
!
            if (FAST_INT .eq. 1 .and. etype .eq. 'mdlb') then
               call elem_maxwell_fi_hexa(Mdle,fld_flag,nrTest,nrTrial,  &
                  nrdofEE,nrdofH,nrdofE,nrdofQ,nrdofEi,                 &
                  BLOC(2)%nrow,BLOC(5)%nrow,                            &
                  BLOC(2)%array,ALOC(2,2)%array,ALOC(2,5)%array,        &
                  BLOC(5)%array,ALOC(5,2)%array,ALOC(5,5)%array)
            elseif (FAST_INT .eq. 1 .and. etype .eq. 'mdlp') then
               call elem_maxwell_fi_pris(Mdle,fld_flag,nrTest,nrTrial,  &
                  nrdofEE,nrdofH,nrdofE,nrdofQ,nrdofEi,                 &
                  BLOC(2)%nrow,BLOC(5)%nrow,                            &
                  BLOC(2)%array,ALOC(2,2)%array,ALOC(2,5)%array,        &
                  BLOC(5)%array,ALOC(5,2)%array,ALOC(5,5)%array)
            else
               call elem_maxwell(Mdle,fld_flag,nrTest,nrTrial,          &
                  nrdofEE,nrdofH,nrdofE,nrdofQ,nrdofEi,                 &
                  BLOC(2)%nrow,BLOC(5)%nrow,                            &
                  BLOC(2)%array,ALOC(2,2)%array,ALOC(2,5)%array,        &
                  BLOC(5)%array,ALOC(5,2)%array,ALOC(5,5)%array)
            endif
!
         else
!        ...notify user that PLANE_PUMP=0 is required for pump Maxwell solve
            if (PLANE_PUMP.ne.0) then
               if (RANK.eq.ROOT) then
                  write(*,*) 'elem: pump solve via Maxwell requires PLANE_PUMP=0. stop.'
               endif
               stop
            endif
!        ...pump case
            Itest(3)=1; Itrial(3)=1
            Itest(6)=1; Itrial(6)=1
            fld_flag = 0;
            if (FAST_INT .eq. 1 .and. etype .eq. 'mdlb') then
               call elem_maxwell_fi_hexa(Mdle,fld_flag,nrTest,nrTrial,  &
                  nrdofEE,nrdofH,nrdofE,nrdofQ,nrdofEi,                 &
                  BLOC(3)%nrow,BLOC(6)%nrow,                            &
                  BLOC(3)%array,ALOC(3,3)%array,ALOC(3,6)%array,        &
                  BLOC(6)%array,ALOC(6,3)%array,ALOC(6,6)%array)
            elseif (FAST_INT .eq. 1 .and. etype .eq. 'mdlp') then
               call elem_maxwell_fi_pris(Mdle,fld_flag,nrTest,nrTrial,  &
                  nrdofEE,nrdofH,nrdofE,nrdofQ,nrdofEi,                 &
                  BLOC(3)%nrow,BLOC(6)%nrow,                            &
                  BLOC(3)%array,ALOC(3,3)%array,ALOC(3,6)%array,        &
                  BLOC(6)%array,ALOC(6,3)%array,ALOC(6,6)%array)
            else
               call elem_maxwell(Mdle,fld_flag,nrTest,nrTrial,          &
                  nrdofEE,nrdofH,nrdofE,nrdofQ,nrdofEi,                 &
                  BLOC(3)%nrow,BLOC(6)%nrow,                            &
                  BLOC(3)%array,ALOC(3,3)%array,ALOC(3,6)%array,        &
                  BLOC(6)%array,ALOC(6,3)%array,ALOC(6,6)%array)
            endif
         endif
!..end select for select case(NO_PROBLEM)
   end select
!
!..end timer
!   end_time = MPI_Wtime()
!      !$OMP CRITICAL
!      write(*,10) etype, end_time-start_time
!      !write(*,11) end_time-start_time
! 10   format(A,' elem : ',f12.5,'  seconds')
! 11   format(f12.5)
!      !$OMP END CRITICAL
!
end subroutine elem
!
!
!----------------------------------------------------------------------
!  routine: compute_enriched_order
!----------------------------------------------------------------------
!  purpose: - compute enriched order vector based on input Nord
!             (enriched order based on mdle order and Nord only)
!----------------------------------------------------------------------
subroutine compute_enriched_order(EType,Nord, Norder)
!
   use parameters, only : MODORDER
!
   implicit none
!
   character(len=4), intent(in)  :: Etype
   integer         , intent(in)  :: Nord
   integer         , intent(out) :: Norder(19)
!
   integer :: temp(2)
   integer :: nordF(3),nordB(3)
!
!----------------------------------------------------------------------
!
!..see implementation of BrokenExactSequence in shape functions
   select case(Etype)
      case('mdlb')
         call decod(Nord,MODORDER,2, temp) !xy face, z edge
         nordF(1) = temp(1); nordB(3) = temp(2)
         call decod(nordF(1),MODORDER,2, nordB(1:2)) !x edge, y edge
         call encod((/nordB(1),nordB(3)/),MODORDER,2, nordF(2)) !xz face
         call encod(nordB(2:3),MODORDER,2, nordF(3)) !yz face
!     ...edges
         Norder(1:4)   = (/nordB(1),nordB(2),nordB(1),nordB(2)/) !x,y,x,y edges
         Norder(5:8)   = (/nordB(1),nordB(2),nordB(1),nordB(2)/) !x,y,x,y edges
         Norder(9:12)  = nordB(3) !z edges
!     ...faces
         Norder(13:14) = nordF(1) !xy faces
         Norder(15:18) = (/nordF(2),nordF(3),nordF(2),nordF(3)/) !xz,yz,xz,yz faces
!     ...element interior
         Norder(19)    = Nord
      case('mdlp')
         call decod(Nord,MODORDER,2, nordB(1:2))
!     ...edges (xy)
         norder(1:6)   = nordB(1)
!     ...edges (z)
         norder(7:9)   = nordB(2)
!     ...triangular faces
         norder(10:11) = nordB(1)
!     ...quadrilateral faces (z)
         norder(12:14) = Nord
!     ...element interior
         norder(15)    = Nord
      case default
         write(*,*) 'compute_enriched_order: only available for hexa and prism. stop.'
         stop
   end select
!
end subroutine compute_enriched_order
!
