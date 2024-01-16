!
#include "typedefs.h"
!
!-------------------------------------------------------------------------------
!
!   routine name       - propagate_flag
!
!-------------------------------------------------------------------------------
!
!   latest revision    - June 2021
!
!   purpose            - Propagates Nflag (customized BC flag 2-9) from
!                        element faces to element edges and vertices.
!                        The flag is passed provided ALL adjacent faces
!                        share the flag (or have a Dirichlet flag).
!                        Since a node stores only a 0/1 BC flag per component,
!                        propagating Nflag corresponds to setting BC flag = 1
!                        on the corresponding component (Icomp).
!                        Background:
!                        Impedance boundary condition (BC) should not be
!                        inherited by edges and vertices from a face
!                        during refinement, unless the edge/vertex is only
!                        adjacent to impedance and dirichlet faces.
!
!   arguments:
!     in:
!              Icomp   - Physics attribute component number (1,..,NRINDEX_HEV)
!              Nflag   - A custom BC flag (2-9); e.g., impedance BC flag
!
!-------------------------------------------------------------------------------
!
! TODO: IF POSSIBLE, THE PROPAGATE_FLAG ROUTINE SHOULD BE AUTOMATED --
!       IT IS TOO COMPLICATED TO BE SUPPLIED BY THE USER.
!       With new BC implementation, to propagate impendance flags correctly,
!       this routine must check faces for the impedance flag (Nflag=3) on the
!       \hat H component (Icomp=3 [signal], Icomp=5 [pump]).
!       - For nodes that are not adjacent to ONLY impedance/Dirichlet faces,
!         the \hat H component (Icomp) BC flag must be set to 0 (instead of 1).
!       - For nodes that are adjacent to ONLY impedance/Dirichlet faces, the BC
!         flag of \hat H must be set to 1 (DOFs treated like Dirichlet DOFs).
!
subroutine propagate_flag(Icomp,Nflag)
!
   use data_structure3D
   use commonParam, only: IBCFLAG
!
   implicit none
!
   integer, intent(in) :: Icomp,Nflag
!
   integer :: etype
   integer :: iel,mdle,ifc,nrfn,i,j,nod
!
!..element nodes and orientations, face nodes
   integer :: nodesl(27),norientl(27),nface_nodes(9)
!
!..element face BC flags, decoded BC flag for a node
   integer :: ibc(6,NRINDEX),nodflag(NRINDEX_HEV)
!
!-------------------------------------------------------------------------------
!
   if (IBCFLAG .ne. 3) then
      write(*,*) 'propagate_flag called for IBCFLAG.ne.3, returning...'
      return
   endif
!
   if ((Icomp.lt.1) .or. (Icomp.gt.NRINDEX_HEV)) then
      write(*,*) 'propagate_flag: invalid Icomp = ', Icomp
      return
   endif
!
!..loop through active elements
!$OMP PARALLEL                                     &
!$OMP PRIVATE(etype,mdle,ifc,nrfn,i,j,nod,nodesl,  &
!$OMP         norientl,nface_nodes,ibc,nodflag)
!
!$OMP DO
   do iel=1,NRELES
      mdle = ELEM_ORDER(iel)
      etype = NODES(mdle)%ntype
!
!  ...determine element nodes
      call elem_nodes(mdle, nodesl,norientl)
!
!  ...get the element boundary conditions flags
      call find_bc(mdle, ibc)
!
!  ...loop through element faces
      do ifc=1,nface(etype)
!
!     ...if face has a Dirichlet BC flag on this component,
!        then neither propagate Nflag from this face to its edges/vertices,
!        nor prohibit another face from passing Nflag to the edges/vertices.
         if (ibc(ifc,Icomp).eq.1) cycle
!
!     ...determine face node numbers
         call face_nodes(etype,ifc, nface_nodes,nrfn)
!
!     ...loop through the face nodes
!$OMP CRITICAL
         do i=1,nrfn !-1
            j = nface_nodes(i)
            nod = nodesl(j)
!
!        ...if node belongs to a face that has impedance BC (Nflag),
!           then propagate the flag unless prohibited by another adjacent face
            if (ibc(ifc,Icomp).eq.Nflag) then
               if (NODES(nod)%visit.ne.-Nflag) then
                  NODES(nod)%visit = Nflag
               endif
!        ...prohibit the flag to be passed to the node
!           (if node belongs to a face that has no impedance or Dirichlet BC)
            else
               NODES(nod)%visit = -Nflag
            endif
         enddo
!$OMP END CRITICAL
      enddo
   enddo
!$OMP END DO
!
!..change -Nflag to zero
!$OMP DO
   do nod=1,NRNODS
      if (NODES(nod)%visit.eq.0) cycle
      call decod(NODES(nod)%bcond,2,NRINDEX_HEV, nodflag)
      if (NODES(nod)%visit.eq.-Nflag) then
         nodflag(Icomp) = 0
      elseif (NODES(nod)%visit.eq.Nflag) then
         nodflag(Icomp) = 1
      endif
      call encod(nodflag,2,NRINDEX_HEV, NODES(nod)%bcond)
!
   enddo
!$OMP END DO
!$OMP END PARALLEL
!
   call reset_visit
!
end subroutine propagate_flag
!
!
!-------------------------------------------------------------------------------
!..just to display the current size of the
!  data structure NODES (in bytes)
!-------------------------------------------------------------------------------
subroutine my_sizetest
!
   use data_structure3D
!
   implicit none
!
   integer :: nod, nH,nE,nV,nQ
   integer(8) :: size_coord, size_h1, size_hcurl, size_hdiv, size_l2, size_tot
   size_coord = 0
   size_h1 = 0
   size_hcurl = 0
   size_hdiv = 0
   size_l2 = 0
!
   do nod=1,NRNODS
       if (.not. associated(NODES(nod)%dof)) cycle
       if (associated(NODES(nod)%dof%coord)) then
         size_coord = size_coord + STORAGE_SIZE(NODES(nod)%dof%coord)
       endif
       if (associated(NODES(nod)%dof%zdofH)) then
         size_h1 = size_h1 + STORAGE_SIZE(NODES(nod)%dof%zdofH)
       endif
       if (associated(NODES(nod)%dof%zdofE)) then
         size_hcurl = size_hcurl + STORAGE_SIZE(NODES(nod)%dof%zdofE)
       endif
       if (associated(NODES(nod)%dof%zdofV)) then
         size_hdiv = size_hdiv + STORAGE_SIZE(NODES(nod)%dof%zdofV)
       endif
       if (associated(NODES(nod)%dof%zdofQ)) then
         size_l2 = size_l2 + STORAGE_SIZE(NODES(nod)%dof%zdofQ)
       endif
   size_tot = size_coord + size_h1 + size_hcurl + size_hdiv + size_l2
   enddo
!
   write(*,*) '1: total DOF size is: ', size_tot
!
   nH=0; nE=0; nV=0; nQ=0
!
   write(*,*) 'my_sizetest: NRHVAR, NREVAR, NRVVAR, NRQVAR = ', &
                            NRHVAR, NREVAR, NRVVAR, NRQVAR
   do nod = 1, NRNODS
      if (Is_inactive(nod)) cycle
      select case(NODES(nod)%ntype)
         case(VERT)
            nH = nH + 1
         case(MEDG)
            nH = nH + MAXP-1
            nE = nE + NREVAR*MAXP
         case(MDLQ)
            nH = nH + NRHVAR*MAXmdlqH
            nE = nE + NREVAR*MAXmdlqE
            nV = nV + NRVVAR*MAXmdlqV
         case(MDLB)
            nH = nH + NRHVAR*MAXmdlbH
            nE = nE + NREVAR*MAXmdlbE
            nV = nV + NRVVAR*MAXmdlbV
            nQ = nQ + NRQVAR*MAXmdlbQ
      end select
   enddo
!
   size_tot = (nH + nE + nV + nQ)*16
   write(*,*) '2: total DOF size is: ', size_tot
   call pause
!
end subroutine my_sizetest
!
!
!-------------------------------------------------------------------------------
!  routine:          set_PML
!  purpose:          sets the PML data
!  last modified:    Jan 2018
!-------------------------------------------------------------------------------
subroutine set_PML
   use commonParam
   implicit none
   PML_REGION=ZL-PML_FRAC*ZL
end subroutine set_PML
!
!
!-------------------------------------------------------------------------------
!   routine:            get_Beta
!
!   last modified:      Oct 2018
!
!   purpose:            returns the PML function at a point Xp
!
!   arguments:
!       in:             Xp       - coordinates (x,y,z) at a physical point
!                       Fld_flag - 1 (signal) / 0 (pump)
!       out:            Zbeta    - PML stretch function that stretches only in z-direction
!                       Zdbeta   - z-derivative of PML stretch function
!                       Zd2beta  - second z-derivative of PML stretch function
!
!-------------------------------------------------------------------------------
subroutine get_Beta(Xp,Fld_flag, Zbeta,Zdbeta,Zd2beta)
!
   use commonParam
   use laserParam
!
   implicit none
!
   real(8), intent(in)  :: Xp(3)
   integer, intent(in)  :: Fld_flag
   VTYPE,   intent(out) :: Zbeta,Zdbeta,Zd2beta
!
   real(8) :: z,a,b,c,L,n,rho,drho,d2rho
   real(8) :: pml_left
!
!..OMEGA_RATIO_SIGNAL or OMEGA_RATIO_PUMP
   real(8) :: OMEGA_RATIO_FLD
!
!-------------------------------------------------------------------------------
!
!..set OMEGA_RATIO_FLD
   select case(Fld_flag)
      case(0)
         OMEGA_RATIO_FLD = OMEGA_RATIO_PUMP
      case(1)
         OMEGA_RATIO_FLD = OMEGA_RATIO_SIGNAL ! 1.0d0
      case default
      write(*,*) 'get_Beta: invalid Fld_flag param. stop.'
         stop
   end select
!
!..TODO: automatically tune PML for envelope wavenumber
   !if (ENVELOPE) OMEGA_RATIO_FLD = OMEGA_RATIO_FLD / 100.d0
   !if (ENVELOPE) OMEGA_RATIO_FLD = OMEGA_RATIO_FLD / 400.d0
   if (ENVELOPE) OMEGA_RATIO_FLD = OMEGA_RATIO_FLD / 800.d0
!
   z = Xp(3)
   b = PML_REGION
   L = ZL
   zbeta = z
   zdbeta = ZONE
!..check if the wave is exponential growth (test case)
   if(.false. .and. SIGMA.ne.ZERO) then ! deactivated (not used currently)
      a = EXP_COEFF
      c = 5.d0
      n = 9.d0
      if(z.gt.b) then
         rho = a*c*(z-b)**n*exp(z/L)
         drho = a*c*(z-b)**(n-1.d0)*dexp(z/L)*(n+(z-b)/L)
         d2rho = a*c*dexp(z/L)*(z-b)**(n-2.d0)*(((z-b)/L)**2+2.d0*n*((z-b)/L)+n*(n-1))
         zbeta = z - ZI*rho/OMEGA
         zdbeta = 1.d0 - ZI*drho/OMEGA
         zd2beta = -ZI*d2rho/OMEGA
      endif
   else
      c = 25.d0
      n = 3.d0
!  ...check for co-pumping: PML on same side for both signal and pump
      if(COPUMP.eq.1) then
         if(z.gt.PML_REGION) then
            rho = c*((z-b)/(ZL-PML_REGION))**n
            drho = c*n*((z-b)/(ZL-PML_REGION))**(n-1.d0)*(1.d0/(ZL-PML_REGION))
            d2rho = c*n*(n-1.d0)*((z-b)/(ZL-PML_REGION))**(n-2.d0)*(1.d0/(ZL-PML_REGION)**2)
            if((rho.le.0.d0).or.(drho.le.0.d0).or.(d2rho.le.0)) then
               write(*,*) ' get_Beta: rho, drho,d2rho are negative. stop.'
               stop
            endif
            zbeta = z - ZI*rho/(OMEGA*OMEGA_RATIO_FLD)
            zdbeta = 1.d0 - ZI*drho/(OMEGA*OMEGA_RATIO_FLD)
            zd2beta = -ZI*d2rho/(OMEGA*OMEGA_RATIO_FLD)
         endif
!  ...next check for counter-pumping: PML on opposite side for signal and pump
!  ...PML @ z.gt.PML_REGION for signal and @ z.lt.PML_FRAC*ZL for pump
      elseif(COPUMP.eq.0) then
         pml_left = PML_FRAC*ZL;
!     ...signal
         if((Fld_flag.eq.1) .and. (z.gt.PML_REGION)) then
            rho = c*((z-b)/(ZL-PML_REGION))**n
            drho = c*n*((z-b)/(ZL-PML_REGION))**(n-1.d0)*(1.d0/(ZL-PML_REGION))
            d2rho = c*n*(n-1.d0)*((z-b)/(ZL-PML_REGION))**(n-2.d0)*(1.d0/(ZL-PML_REGION)**2)
            if((rho.le.0.d0).or.(drho.le.0.d0).or.(d2rho.le.0)) then
               write(*,*) ' get_Beta: rho, drho,d2rho are negative. stop.'
               stop
            endif
            zbeta = z - ZI*rho/(OMEGA*OMEGA_RATIO_FLD)
            zdbeta = 1.d0 - ZI*drho/(OMEGA*OMEGA_RATIO_FLD)
            zd2beta = -ZI*d2rho/(OMEGA*OMEGA_RATIO_FLD)
!     ...pump
         else if((Fld_flag.eq.0) .and. (z.lt.pml_left)) then
            rho = c*((pml_left-z)/pml_left)**n
            drho = c*n*((pml_left-z)/pml_left)**(n-1.d0)*(1.d0/(pml_left))
            d2rho = c*n*(n-1.d0)*((pml_left-z)/pml_left)**(n-2.d0)*(1.d0/(pml_left)**2)
            if((rho.le.0.d0).or.(drho.le.0.d0).or.(d2rho.le.0)) then
               write(*,*) ' get_Beta: rho, drho,d2rho are negative. stop.'
               stop
            endif
            zbeta = z - ZI*rho/(OMEGA*OMEGA_RATIO_FLD)
            zdbeta = 1.d0 - ZI*drho/(OMEGA*OMEGA_RATIO_FLD)
            zd2beta = -ZI*d2rho/(OMEGA*OMEGA_RATIO_FLD)
         endif
!  ...endif for COPUMP check
      endif
!..endif for the wave is exponential growth
   endif
!
end subroutine get_Beta
