#include "typedefs.h"

!> @brief   copy local dof into datastructure
!! @param[in] Nod                     - a node number
!! @param[in] Ncoms                   - solution component set: 1,...,NRCOMS
!! @param[in] ZvalH,ZvalE,ZvalV,ZvalQ - H1,H(curl),H(div) and L2 dof
!!                                      for the node in the expanded mode
!> @date    Sep 2023
subroutine dof_in(Nod,Ncoms,ZvalH,ZvalE,ZvalV,ZvalQ)

  use data_structure3D

  implicit none

  ! ** Arguments
  integer, intent(in) :: Nod,Ncoms
  VTYPE,   intent(in) :: ZvalH(MAXEQNH,*), ZvalE(MAXEQNE,*), &
                         ZvalV(MAXEQNE,*), ZvalQ(MAXEQNQ,*)
  !
  ! ** Locals
  ! node physics attributes flags
  integer, dimension(NR_PHYSA) :: ncase
  integer :: i, j, k, ivarH, ivarE, ivarV, ivarQ, icase
  integer :: ndofH, ndofE, ndofV, ndofQ, nvar


  ! decode the physical attributes of the node
  icase = NODES(Nod)%case
  call decod(icase,2,NR_PHYSA, ncase)

  ! determine number of dof for the node
  call find_ndof(Nod, ndofH,ndofE,ndofV,ndofQ)

  ! initiate the node component counters
  ivarH=0; ivarE=0; ivarV=0; ivarQ=0

  if (.not. associated(NODES(Nod)%dof)) then
    write(*,*) 'dof_in: nod dof not associated. returning.'
    return
  endif

  ! loop through multiple loads
  do j=1,NRRHS

     ! loop through physical attributes of the node
     do i=1,NR_PHYSA
        if (ncase(i).ne.1) cycle

        ! loop through the components of the attribute
        do k=1,NR_COMP(i)
           select case(D_TYPE(i))
              case(CONTIN)
                 nvar = (j-1)*NRHVAR+ADRES(i)+k
                 ivarH = ivarH+1
                 if (ndofH.gt.0) then
                    NODES(Nod)%dof%zdofH(ivarH,1:ndofH,Ncoms) = ZvalH(nvar,1:ndofH)
                 endif
              case(TANGEN)
                 nvar = (j-1)*NREVAR+ADRES(i)+k
                 ivarE = ivarE+1
                 if (ndofE.gt.0) then
                    NODES(Nod)%dof%zdofE(ivarE,1:ndofE,Ncoms) = ZvalE(nvar,1:ndofE)
                 endif
              case(NORMAL)
                 nvar = (j-1)*NRVVAR+ADRES(i)+k
                 ivarV = ivarV+1
                 if (ndofV.gt.0) then
                    NODES(Nod)%dof%zdofV(ivarV,1:ndofV,Ncoms) = ZvalV(nvar,1:ndofV)
                 endif
              case(DISCON)
                 nvar = (j-1)*NRQVAR+ADRES(i)+k
                 ivarQ = ivarQ+1
                 if (ndofQ.gt.0) then
                    NODES(Nod)%dof%zdofQ(ivarQ,1:ndofQ,Ncoms) = ZvalQ(nvar,1:ndofQ)
                 endif
              case default
                 write(*,*) 'dofin: D_TYPE = ', S_DType(D_TYPE(i))
                 call logic_error(ERR_INVALID_VALUE,__FILE__,__LINE__)
           end select
        ! loop through the components of the attribute
        enddo
     ! loop through physical attributes of the node
     enddo
  ! loop through multiple loads
  enddo

end subroutine dof_in
