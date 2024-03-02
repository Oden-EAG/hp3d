!--------------------------------------------------------------------
!
!     routine name      - celem
!
!--------------------------------------------------------------------
!
!     latest revision:  - July 2019
!
!> @brief         - calculate modified stiffness matrix and load vector
!
!     arguments:
!        in:
!             Mdle      - middle node of an element
!             Idec      - 1 only info for nodes, 2 matrices
!        out:
!             Nrdofs    - # of element local dof for each physics
!             Nrdofm    - # of modified element dof in the expanded
!             Nrdofc    - # of modified element dof after compression
!             Nodm      - actual nodes in the order
!             NdofmH    - the number of H1 dof
!             NdofmE    - the number of H(curl) dof
!             NdofmV    - the number of H(div) dof
!             NdofmQ    - the number of L2 dof
!             Nrnodm    - number of the modified element nodes
!             Bload     - 1D array containing the modified load vector
!             Astif     - 1D array containing the modified stiffness matrix
!---------------------------------------------------------------------
subroutine celem(Mdle,Idec,Nrdofs,Nrdofm,Nrdofc,Nodm, &
                 NdofmH,NdofmE,NdofmV,NdofmQ,Nrnodm,Bload,Astif)
!
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
!
!--------------------------------------------------------------------------
!
!..redirect to the system routine
   call celem_system(Mdle,Idec, &
                     Nrdofs,Nrdofm,Nrdofc, &
                     Nodm,NdofmH,NdofmE,NdofmV,NdofmQ,Nrnodm, &
                     Bload,Astif)
!
!
end subroutine celem
