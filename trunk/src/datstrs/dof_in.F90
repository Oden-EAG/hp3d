!> Purpose : copy local dof into datastructure
!! @param[in] Nod         - a node number
!! @param[in] ZvalH,ZvalE,ZvalV,ZvalQ - H1,H(curl),H(div) and L2 dof
!!                        for the node in the expanded mode

#include "typedefs.h"
subroutine dof_in(Nod,ZvalH,ZvalE,ZvalV,ZvalQ)
  use data_structure3D
  implicit none
  !
  ! ** Arguements
  integer, intent(in) :: Nod
  VTYPE,   intent(in) :: &
       ZvalH(MAXEQNH,*), ZvalE(MAXEQNE,*), &
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

  ! loop through multiple copies of variables
  do j=1,NRCOMS

     ! loop through physical attributes of the node
     do i=1,NR_PHYSA
        if (ncase(i).eq.1) then
           ! loop through the components of the attribute
           do k=1,NR_COMP(i)
              select case(DTYPE(i))
              case('contin')
                 nvar = (j-1)*NRHVAR+ADRES(i)+k
                 ivarH = ivarH+1
                 if (ndofH.gt.0) then
                    NODES(Nod)%dof%zdofH(ivarH,1:ndofH) = ZvalH(nvar,1:ndofH)
                 endif
              case('tangen')
                 nvar = (j-1)*NREVAR+ADRES(i)+k
                 ivarE = ivarE+1
                 if (ndofE.gt.0) then
                    NODES(Nod)%dof%zdofE(ivarE,1:ndofE) = ZvalE(nvar,1:ndofE)
                 endif
              case('normal')
                 nvar = (j-1)*NRVVAR+ADRES(i)+k
                 ivarV = ivarV+1
                 if (ndofV.gt.0) then
                    NODES(Nod)%dof%zdofV(ivarV,1:ndofV) = ZvalV(nvar,1:ndofV)
                 endif
              case('discon')
                 nvar = (j-1)*NRQVAR+ADRES(i)+k
                 ivarQ = ivarQ+1
                 if (ndofQ.gt.0) then
                    NODES(Nod)%dof%zdofQ(ivarQ,1:ndofQ) = ZvalQ(nvar,1:ndofQ)
                 endif
              case default
                 write(*,*) 'dofin: NOT SUPPORT DTYPE ', DTYPE(i)
                 call logic_error(ERR_INVALID_VALUE,__FILE__,__LINE__)
              end select
           enddo
        endif
     enddo
  enddo
end subroutine dof_in




