!----------------------------------------------------------------------
!> @brief       Calculate modified stiffness matrix and load vector
!!              This is an advanced routine; it is not recommended that
!!              all but the most advanced users modify this routine.
!!
!> @param[in]   mdle    - middle node number
!> @param[in]   idec    - flag whether to compute matrices
!!                         1: only compute constrained elem info
!!                         2: also compute stiffness matrix and load
!> @param[out]  Nrdofs  - # of element local dof for each physics
!> @param[out]  Nrdofm  - # of modified element dof in the expanded
!> @param[out]  Nrdofc  - # of modified element dof after compression
!> @param[out]  Nodm    - unconstrained nodes
!> @param[out]  NdofmH  - the number of H^1 dof
!> @param[out]  NdofmE  - the number of H(curl) dof
!> @param[out]  NdofmV  - the number of H(div) dof
!> @param[out]  NdofmQ  - the number of L^2 dof
!> @param[out]  Nrnodm  - number of the modified element nodes
!> @param[out]  Bload   - 1D array containing the modified load vector
!> @param[out]  Astif   - 1D array containing the modified stiffness matrix
!!                        in packed form
!!
!> @date       July 2023
!----------------------------------------------------------------------
   subroutine celem(Mdle,Idec,Nrdofs,Nrdofm,Nrdofc,Nodm, &
                    NdofmH,NdofmE,NdofmV,NdofmQ,Nrnodm,Bload,Astif)
   use physics
   use data_structure3D
!
!--------------------------------------------------------------------------
! 
  implicit none
!  
   integer,                      intent(in)  :: Mdle,Idec
   integer, dimension(NR_PHYSA), intent(out) :: Nrdofs
   integer,                      intent(out) :: Nrdofm,Nrdofc
   integer, dimension(MAXNODM),  intent(out) :: Nodm
   integer, dimension(MAXNODM),  intent(out) :: NdofmH,NdofmE,NdofmV,NdofmQ
   integer,                      intent(out) :: Nrnodm
   complex(8),                   intent(out) :: Bload(*),Astif(*)
   integer, dimension(NRINDEX)               :: nbcond
!   
!--------------------------------------------------------------------------
!
!..This is a hack to eliminate the bubbles in the trace physics variables,
!  so that only the boundary dof are included for the trace variables.
!  This is done by saying that the value is known (to be zero) for the
!  bubble (middle/interior) dof by using the BC flag (1 if known, 0 if unknown)
!  The 'known' dof are eliminated by static condensation locally
   nbcond = 0; nbcond(1:2) = 1; ! removes the bubbles
   call encod(nbcond,2,NRINDEX, NODES(Mdle)%bcond)
!
!..redirect to the system routine
   call celem_system(Mdle,Idec, &
                     Nrdofs,Nrdofm,Nrdofc, &
                     Nodm,NdofmH,NdofmE,NdofmV,NdofmQ,Nrnodm, &
                     Bload,Astif)
!
!..reset the BC flags back to zero (unknown)
   NODES(Mdle)%bcond = 0
!
!
   end subroutine celem
