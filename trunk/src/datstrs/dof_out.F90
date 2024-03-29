#include "typedefs.h"

!> @brief   copy dof for a node from data structure into local arrays
!> @param[in]     Nod                     - a node number
!> @param[in]     Ncoms                   - solution component set: 1,...,NRCOMS
!> @param[in,out] KdofH,KdofE,KdofV,KdofQ - dof counters for the local arrays
!> @param[out]    ZvalH,ZvalE,ZvalV,ZvalQ - H1,H(curl),H(div) and L2 dof from the
!!                                          data structure in the expanded mode
!> @date    Sep 2023
subroutine dof_out( Nod,Ncoms,               &
                    KdofH,KdofE,KdofV,KdofQ, &
                    ZvalH,ZvalE,ZvalV,ZvalQ  )

  use data_structure3D

  implicit none

  ! ** Arguments
  integer, intent(in)    :: Nod, Ncoms
  integer, intent(inout) :: KdofH, KdofE, KdofV, KdofQ
  VTYPE,   intent(out)   :: ZvalH(MAXEQNH,*), ZvalE(MAXEQNE,*), &
                            ZvalV(MAXEQNV,*), ZvalQ(MAXEQNQ,*)

  ! ** Locals
  integer :: ncase(NR_PHYSA)
  integer :: i, j, k, ivarH, ivarE, ivarV, ivarQ
  integer :: ndofH, ndofE, ndofV, ndofQ, nvar

#if HP3D_DEBUG
  integer :: iprint
  iprint=0
#endif

  ! decode the physical attributes of the node
  call decod(NODES(Nod)%case,2,NR_PHYSA, ncase)

  ! determine number of dof for the node
  call find_ndof(Nod, ndofH,ndofE,ndofV,ndofQ)

  ! initiate the node component counters
  ivarH=0; ivarE=0; ivarV=0; ivarQ=0

  if (.not. associated(NODES(Nod)%dof)) then
    write(*,*) 'dof_out: nod dof not associated. returning.'
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

                    if (.not.associated(NODES(Nod)%dof%zdofH)) then
                      write(*,*)'zdofH not associated for Nod = ',Nod
                      call pause
                      return
                    endif

                    ZvalH(nvar,KdofH+1:KdofH+ndofH) = &
                         NODES(Nod)%dof%zdofH(ivarH,1:ndofH,Ncoms)
#if HP3D_DEBUG
                    if (iprint.eq.1) then
                       write(*,7001) ivarH,nvar,KdofH+1,KdofH+ndofH
                       write(*,*) 'ZvalH = ',ZvalH(nvar,KdofH+1:KdofH+ndofH)
                    endif
#endif
                 endif

              case(TANGEN)
                 nvar = (j-1)*NREVAR+ADRES(i)+k
                 ivarE = ivarE+1
                 if (ndofE.gt.0) then
                    ZvalE(nvar,KdofE+1:KdofE+ndofE) = &
                         NODES(Nod)%dof%zdofE(ivarE,1:ndofE,Ncoms)
#if HP3D_DEBUG
                    if (iprint.eq.2) then
                       write(*,7001) ivarE,nvar,KdofE+1,KdofE+ndofE
                       write(*,*) 'ZvalE = ', ZvalE(nvar,KdofE+1:KdofE+ndofE)
                    endif
#endif
                 endif

              case(NORMAL)
                 nvar = (j-1)*NRVVAR+ADRES(i)+k
                 ivarV = ivarV+1
                 if (ndofV.gt.0) then
                    ZvalV(nvar,KdofV+1:KdofV+ndofV) = &
                         NODES(Nod)%dof%zdofV(ivarV,1:ndofV,Ncoms)
#if HP3D_DEBUG
                    if (iprint.eq.3) then
                       write(*,7001) ivarV,nvar,KdofV+1,KdofV+ndofV
                       write(*,*) 'ZvalV = ', ZvalV(nvar,KdofV+1:KdofV+ndofV)
                    endif
#endif
                 endif

              case(DISCON)
!                 write(*,*) 'dof_out: ndofQ = ',ndofQ
                 nvar = (j-1)*NRQVAR+ADRES(i)+k
                 ivarQ = ivarQ+1
                 if (ndofQ.gt.0) then
                    ZvalQ(nvar,KdofQ+1:KdofQ+ndofQ) = &
                         NODES(Nod)%dof%zdofQ(ivarQ,1:ndofQ,Ncoms)
                 endif
!                 if (ndofQ.gt.0) then
!                   write(*,*) 'nvar,ivarQ,KdofQ = ',nvar,ivarQ,KdofQ
!                   write(*,*) 'ZvalQ(nvar,KdofQ+1:KdofQ+ndofQ) = ', &
!                               ZvalQ(nvar,KdofQ+1:KdofQ+ndofQ)
!                   write(*,*) 'NODES(Nod)%dof%zdofQ(ivarQ,1:ndofQ,Ncoms) = ', &
!                               NODES(Nod)%dof%zdofQ(ivarQ,1:ndofQ,Ncoms)
!                 endif

              case default
                 write(*,*) 'dofout: D_TYPE = ', S_DType(D_TYPE(i))
           end select
7001       format('dof_out: ivar = ', i3,2x,'(',i3,',',i3,':',i3,') = ')
        ! loop through the components of the attribute
        enddo
     ! loop through physical attributes of the node
     enddo
  ! loop through multiple loads
  enddo

  ! update the local dof counters
  KdofH = KdofH + ndofH
  KdofE = KdofE + ndofE
  KdofV = KdofV + ndofV
  KdofQ = KdofQ + ndofQ

end subroutine dof_out
