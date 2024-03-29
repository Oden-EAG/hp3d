!-----------------------------------------------------------------------
!
!   routine name       - coarse2macro
!
!-----------------------------------------------------------------------
!
!   latest revision    - Feb 2018
!
!   purpose            - routine applies the prolongation operator to
!                        a coarse-grid solution vector using the
!                        coefficients from logic_macro
!                        (including Dirichlet dof)
!
!   arguments :
!     in:
!             Igrid    -
!             Iel      - natural order of elements
!             Zu_c     - solution vector for the modified coarse element
!             nrdof_c  - length of Zu_c
!     out:
!             Zu_m     - solution vector for the modified macro-element
!             nrdof_m  - length of Zu_m
!
!----------------------------------------------------------------------
!
!..TODO accommodate usage of PHYSAm (loop through physics variables)
!..     see solout_mg or compute_sol_dof for example
!
#include "typedefs.h"
!
   subroutine coarse2macro(Igrid,Iel,Zu_c,Nrdof_c, Zu_m,Nrdof_m)
!
   use data_structure3D,  only: NODES, get_index
   use mg_data_structure, only: GRID
   use physics,           only: NRHVAR, NREVAR, NRVVAR, NRINDEX
   use parameters,        only: NRCOMS, MAXEQNH, MAXEQNE, MAXEQNV,      &
                                MAXbrickH, MAXbrickE, MAXbrickV, ZERO
   use assembly,          only: NR_RHS
!
   implicit none
!
!
   integer, intent(in)  :: Igrid,Iel, Nrdof_c, Nrdof_m
!
   VTYPE, intent(in)  :: Zu_c(nrdof_c)
   VTYPE, intent(out) :: Zu_m(nrdof_m)
!
!..locals
   VTYPE           :: zdofH_c(MAXEQNH,MAXbrickH),     &
                      zdofE_c(MAXEQNE,MAXbrickE),     &
                      zdofV_c(MAXEQNV,MAXbrickV)
   VTYPE           :: zdofH_m(MAXEQNH,GRID(Igrid)%constr(Iel)%nrdofH_macro),   &
                      zdofE_m(MAXEQNE,GRID(Igrid)%constr(Iel)%nrdofE_macro),   &
                      zdofV_m(MAXEQNV,GRID(Igrid)%constr(Iel)%nrdofV_macro)
!
!..component counters for the nodes (used in case of multiple loads)
   integer, allocatable :: mvarH(:),mvarE(:), mvarV(:)
!
   integer              :: iH(MAXEQNH), iE(MAXEQNE), iV(MAXEQNV)
   integer              :: nrnod_macro, nrdofH_macro, nrdofE_macro, nrdofV_macro
   integer              :: nrnodm, nvarH, nvarE, nvarV, ivar
   integer              :: naH, naE, naV,i,k,load,nn,nod,j
   integer              :: kH, kE, kV, kp,lH, lE,lV
!
!..decoded index for a node
   integer              :: index(NRINDEX)
!
!--------------------------------------------------------------------------
!
!..macro-element number of nodes
   nrnod_macro = GRID(Igrid)%constr(Iel)%nrnod_macro
!
!..modified coarse element number of nodes
   nrnodm   = GRID(Igrid)%constr(Iel)%nrnodm
!
!--------------------------------------------------------------------------
!
!..compute total number of dof for the macro-element
   naH = 0; naE = 0; naV = 0
   do i=1,nrnod_macro
      naH = naH + GRID(Igrid)%constr(Iel)%ndofH_macro(i)
      naE = naE + GRID(Igrid)%constr(Iel)%ndofE_macro(i)
      naV = naV + GRID(Igrid)%constr(Iel)%ndofV_macro(i)
   enddo

   nrdofH_macro = naH; nrdofE_macro = naE; nrdofV_macro = naV

!..initiate dof counter
   nn = 0;
!..initialize the number of H1,H(curl),H(div),L2 copied so far
   ivar=0

!..initialize component counters
   allocate(mvarH(nrnodm)); mvarH = 0
   allocate(mvarE(nrnodm)); mvarE = 0
   allocate(mvarV(nrnodm)); mvarV = 0

!..initialize dof matrices
   zdofH_c = ZERO; zdofE_c = ZERO; zdofV_c = ZERO
!
!..loop through loads
   do load=1,NR_RHS
!
!  ...H1 dof
      if (NRHVAR.eq.0) goto 100
!
!  ...initialize H1 column indices
      iH = 0
!
!  ...loop through coarse element nodes
      do i = 1, nrnodm
!
!     ...pick up the node number from the list
         nod = GRID(Igrid)%constr(Iel)%nodm(i)
!
!     ...compute the number of active H1 variables for the node
         call get_index(nod, index)
!
!     ...loop through the nodal dof (potentially NONE)
         do j=1,GRID(Igrid)%constr(Iel)%ndofmH(i)
!        ...loop through the components
            ivar=mvarH(i)
            do k=1,NRINDEX
               select case(index(k))
               case(1)
                  ivar = ivar + 1
                  iH(ivar) = iH(ivar) + 1
               case(2)
                  ivar = ivar + 1
                  iH(ivar) = iH(ivar) + 1
                  nn = nn+1
                  zdofH_c(ivar,iH(ivar)) = Zu_c(nn)
               end select
!
!        ...end of loop through NRINDEX
            enddo
!     ...end of loop through the nodal dof
         enddo
!
!     ...update the number of components stored so far
         mvarH(i) = ivar

!  ...end of loop through coarse element nodes
      enddo
!
  100 continue
!
!  ...H(curl) dof
      if (NREVAR.eq.0) goto 200
      ivar = 0
!
!  ...initialize H(curl) column indices
      iE = 0
!
!  ...loop through coarse element nodes
      do i = 1, nrnodm
!     ...pick up the node number from the list
         nod = GRID(Igrid)%constr(Iel)%nodm(i)
!
!     ...compute the number of active H(curl) variables for the node
         call get_index(nod, index)
!
         do j=1,GRID(Igrid)%constr(Iel)%ndofmE(i)
!
!        ...loop through the components
            ivar=mvarE(i)
            do k=1,NRINDEX
               select case(index(k))
               case(3)
                  ivar=ivar+1
                  iE(ivar) = iE(ivar) + 1
               case(4)
                  ivar=ivar+1
                  iE(ivar) = iE(ivar) + 1
                  nn=nn+1
                  zdofE_c(ivar,iE(ivar)) = Zu_c(nn)
               end select
!
!        ...end of loop through NRINDEX
            enddo
!
!     ...end of loop through the nodal dof
         enddo
!
!     ...update the number of components stored so far
         mvarE(i) = ivar

!  ...end of loop through coarse element nodes
      enddo
!
  200 continue
!
!  ...H(div) dof
      if (NRVVAR.eq.0) goto 300
!
!  ...initialize H(div) column indices
      iV = 0
!
!  ...loop through coarse element nodes
      do i = 1, nrnodm
!     ...pick up the node number from the list
         nod = GRID(Igrid)%constr(Iel)%nodm(i)
!
!     ...compute the number of active H(curl) variables for the node
         call get_index(nod, index)
!
         do j=1,GRID(Igrid)%constr(Iel)%ndofmV(i)!/nvarV
!
!        ...loop through the components
            ivar=mvarV(i)
            do k=1,NRINDEX
               select case(index(k))
               case(5)
                  ivar=ivar+1
                  iV(ivar) = iV(ivar) + 1
               case(6)
                  ivar=ivar+1
                  iV(ivar) = iV(ivar) + 1
                  nn=nn+1
                  zdofV_c(ivar,iV(ivar)) = Zu_c(nn)
               end select
!
!        ...end of loop through NRINDEX
            enddo
!
!     ...end of loop through the nodal dof
         enddo
!
!     ...update the number of components stored so far
         mvarV(i) = ivar

!  ...end of loop through coarse element nodes
      enddo

  300 continue

!..end of loop through loads
   enddo
!
   deallocate(mvarH, mvarE, mvarV)
!
!--------------------------------------------------------------------
!
!...calculate the macro element dof
!
!..initialize the macro dof matrices
   ZdofH_m = ZERO; ZdofE_m = ZERO; ZdofV_m = ZERO
!
!..loop through the macro-element H1 dof
   do kH=1,nrdofH_macro
!
!  ...accumulate for the values
      do kp=1,GRID(Igrid)%constr(Iel)%nrconH_macro(kH)
         lH = GRID(Igrid)%constr(Iel)%nacH_macro(kp,kH)
         do ivar=1,NRHVAR*NRCOMS
            ZdofH_m(ivar,kH) = ZdofH_m(ivar,kH)     &
                             + GRID(Igrid)%constr(Iel)%constrH_macro(kp,kH)*zdofH_c(ivar,lH)
         enddo
      enddo
!...end of loop through macro-element H1 dof
   enddo
!
!..loop through the macro-element H(curl) dof
   do kE=1,nrdofE_macro
!
!  ...accumulate for the values
      do kp=1,GRID(Igrid)%constr(Iel)%nrconE_macro(kE)
         lE = GRID(Igrid)%constr(Iel)%nacE_macro(kp,kE)
         do ivar=1,NREVAR*NRCOMS
            ZdofE_m(ivar,kE) = ZdofE_m(ivar,kE)    &
                             + GRID(Igrid)%constr(Iel)%ConstrE_macro(kp,kE)*zdofE_c(ivar,lE)
         enddo
      enddo
!...end of loop through macro-element H(curl) dof
   enddo
!
!..loop through the macro-element H(curl) dof
   do kV=1,nrdofV_macro
!
!  ...accumulate for the values
      do kp=1,GRID(Igrid)%constr(Iel)%nrconV_macro(kV)
         lV = GRID(Igrid)%constr(Iel)%nacV_macro(kp,kV)
         do ivar=1,NRVVAR*NRCOMS
            ZdofV_m(ivar,kV) = ZdofV_m(ivar,kV)      &
                             + GRID(Igrid)%constr(Iel)%constrV_macro(kp,kV)*zdofV_c(ivar,lV)
         enddo
      enddo
!...end of loop through macro-element H(curl) dof
   enddo
!
!..reconstruct the macro element solution vector
!
!..initiate dof counter
   nn = 0;
!..initialize the number of H1,H(curl),H(div),L2 copied so far
   ivar=0

!..initiate component counters
   allocate(mvarH(nrnod_macro)); mvarH = 0
   allocate(mvarE(nrnod_macro)); mvarE = 0
   allocate(mvarV(nrnod_macro)); mvarV = 0
!
!..loop through loads
   do load=1,NR_RHS
!
!  ...H1 dof
      if (NRHVAR.eq.0) goto 400
!
!  ...initialize H1 column indices
      iH = 0
!
!  ...loop through macro element nodes
      do i = 1, nrnod_macro
!
!     ...pick up the node number from the list
         nod = GRID(Igrid)%constr(Iel)%nod_macro(i)
!
!     ...compute the number of active H1 variables for the node
         call get_index(nod, index)
!
!     ...loop through the nodal dof (potentially NONE)
         do j=1,GRID(Igrid)%constr(Iel)%ndofH_macro(i)
!        ...loop through the components
            ivar=mvarH(i)
            do k=1,NRINDEX
               select case(index(k))
               case(1)
                  ivar = ivar + 1
                  iH(ivar) = iH(ivar) + 1
               case(2)
                  ivar = ivar + 1
                  iH(ivar) = iH(ivar) + 1
                  nn = nn+1
                  Zu_m(nn) = zdofH_m(ivar,iH(ivar))
               end select
!
!        ...end of loop through NRINDEX
            enddo
!     ...end of loop through the nodal dof
         enddo
!
!     ...update the number of components stored so far
         mvarH(i) = ivar

!  ...end of loop through macro element nodes
      enddo
!
  400 continue
!
!  ...H(curl) dof
      if (NREVAR.eq.0) goto 500
!
!  ...initialize H(curl) column indices
      iE = 0
!
!  ...loop through macro element nodes
      do i = 1, nrnod_macro
!     ...pick up the node number from the list
         nod = GRID(Igrid)%constr(Iel)%nod_macro(i)
!
!     ...compute the number of active H(curl) variables for the node
         call get_index(nod, index)
!
         do j=1,GRID(Igrid)%constr(Iel)%ndofE_macro(i)!/nvarE
!
!        ...loop through the components
            ivar=mvarE(i)
            do k=1,NRINDEX
               select case(index(k))
               case(3)
                  ivar=ivar+1
                  iE(ivar) = iE(ivar) + 1
               case(4)
                  ivar=ivar+1
                  iE(ivar) = iE(ivar) + 1
                  nn=nn+1
                  Zu_m(nn) = zdofE_m(ivar,iE(ivar))
               end select
!
!        ...end of loop through NRINDEX
            enddo
!
!     ...end of loop through the nodal dof
         enddo
!
!     ...update the number of components stored so far
         mvarE(i) = ivar

!  ...end of loop through macro element nodes
      enddo

  500 continue
!
!  ...H(div) dof
      if (NRVVAR.eq.0) goto 600
!
!  ...initialize H(div) column indices
      iV = 0
!
!  ...loop through macro element nodes
      do i = 1, nrnod_macro

!     ...pick up the node number from the list
         nod = GRID(Igrid)%constr(Iel)%nod_macro(i)
!
!     ...compute the number of active H(curl) variables for the node
         call get_index(nod, index)
!
         do j=1,GRID(Igrid)%constr(Iel)%ndofV_macro(i)!/nvarV
!
!        ...loop through the components
            ivar=mvarV(i)
            do k=1,NRINDEX
               select case(index(k))
               case(5)
                  ivar=ivar+1
                  iV(ivar) = iV(ivar) + 1
               case(6)
                  ivar=ivar+1
                  iV(ivar) = iV(ivar) + 1
                  nn=nn+1
                  Zu_m(nn) = zdofV_m(ivar,iV(ivar))
               end select
!
!        ...end of loop through NRINDEX
            enddo
!
!     ...end of loop through the nodal dof
         enddo
!
!     ...update the number of components stored so far
         mvarV(i) = ivar

!  ...end of loop through macro element nodes
      enddo
!
  600 continue
!
!..end of loop through loads
   enddo
!
   deallocate(mvarH, mvarE, mvarV)
!
!
   end subroutine coarse2macro
