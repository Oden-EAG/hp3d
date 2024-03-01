c----------------------------------------------------------------------
c
c   routine name       - solin2
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 2023
c
c   purpose            - an interface routine for the frontal solver
c                        returning for an element its stiffness matrix
c                        and load vector
c
c   arguments :
c     in:
c             Iel      - sequential element number
c             Mdum
c             Ifg
c             Nrhs0
c             Mdest
c     out:
c             Zalhs    - stiffness matrix
c             Zarhs    - load vector
c
c----------------------------------------------------------------------
#include "typedefs.h"
      subroutine solin2(Iel,Mdum,Ifg,Nrhs0,Mdest, Zalhs,Zarhs)
c
      use control         , only: ISTC_FLAG
      use physics         , only: NR_PHYSA
      use data_structure3D, only: MAXNODM
      use frsolmod
      implicit none
c
      integer, intent(in)  :: Iel,Mdum,Ifg,Nrhs0,Mdest
      VTYPE  , intent(out) :: Zalhs(*),Zarhs(*)
c
c  ...nodes for a modified element and the corresponding number
c     of H1,H(curl),H(div) and L2 dof
      integer :: nodm(MAXNODM),ndofmH(MAXNODM),ndofmE(MAXNODM),
     .                         ndofmV(MAXNODM),ndofmQ(MAXNODM)
c
c  ...number of variables for each physics attribute for an element
      integer :: nrdofs(NR_PHYSA)
c
      integer :: mdle,i,nrdofm,nrdofc,nrnodm
c
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
c
c---------------------------------------------------------------------
c
      if (REORDER) then
c
c  .....find the element number
        mdle = NEW_ELEM_ORDER(Iel)
      else
c
c  .....find the element number using the standard ordering of elements
        mdle=0
        do i=1,Iel
          call nelcon(mdle, mdle)
        enddo
      endif
c
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,7001) Iel,mdle
 7001   format('solin2: Iel,mdle = ',i8,i10)
      endif
#endif
c
c  ...evaluate the element load vector and stiffness matrix
      if (ISTC_FLAG) then
        call celem_systemI(Iel,mdle,2,
     .           nrdofs,nrdofm,nrdofc,
     .           nodm,ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,
     .           Zarhs,Zalhs)
      else
        call celem(mdle,2,
     .           nrdofs,nrdofm,nrdofc,
     .           nodm,ndofmH,ndofmE,ndofmV,ndofmQ,nrnodm,
     .           Zarhs,Zalhs)
      endif
c
#if HP3D_DEBUG
      if (iprint.eq.1) then
         write(*,*)'Element',mdle,'RHS'
         do i=1,nrdofc
           if(abs(Zarhs(i)).gt.1.d-8)then
           endif
         enddo
      endif
#endif
      end subroutine solin2
