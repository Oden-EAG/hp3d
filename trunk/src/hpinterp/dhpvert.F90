!> Purpose :  routine determines Dirichlet dof for a vertex
!! @param[in]  Iflag - a flag specifying which of the objects the vertex
!! @param[in]  No    - number of a specific object
!! @param[in]  Xi    - reference coordinates of the point
!! @param[in]  Icase - node case
!!
!! @param[out] ZdofH - updated dirichlet bc
#include "implicit_none.h"
subroutine dhpvert(Mdle,Iflag,No,Xi,Icase, ZdofH)
  use data_structure3D
  implicit none
  ! ** Argumenet
  !---------------------------------------------------------------------
  integer,                                 intent(in)  :: Iflag,No,Icase,Mdle
  real*8, dimension(3),                    intent(in)  :: Xi
  VTYPE,  dimension(NRCOMS*NREQNH(Icase)), intent(out) :: ZdofH

  ! ** Locals
  !---------------------------------------------------------------------
  ! Dirichlet BC data
  VTYPE :: &
       zvalH(MAXEQNH), &
       zdvalH(MAXEQNH,3), zdvalHdxi(MAXEQNH,3),zdvalHdt(MAXEQNH), &
       zvalE(3,MAXEQNE), zdvalE(3,MAXEQNE,3), &
       zvalV(3,MAXEQNV), zdvalV(3,MAXEQNV,3)

  real*8,  dimension(3)   :: x
  real*8,  dimension(3,3) :: void

  integer, dimension(NR_PHYSA) :: ncase
  integer :: iprint, ivarH, nvarH, iattr,iload,icomp
  !---------------------------------------------------------------------
      iprint=0
!
!  ...determine coordinates of the point
      select case(Iflag)
      case(5)                 ; call prism(No,Xi, x,void)
      case(6)                 ; call  hexa(No,Xi, x,void)
      case(7)                 ; call tetra(No,Xi, x,void)
      case(8)                 ; call pyram(No,Xi, x,void)
      case default
        write(*,1000)iflag
1000    format(' dhpvert: unknown iflag = ',i7)
        stop
      end select
      if (iprint.eq.1) then
        write(*,7010) No, Xi, x
 7010   format('dhpvert: No = ',i3,' Xi = ',3f8.3,' x = ',3f8.3)
      endif
!
!  ...dirichlet value in compact mode
      call dirichlet(Mdle,x,Icase, zvalH,zdvalH,zvalE,zdvalE,zvalV,zdvalV)
!
!  ...shift the data skipping irrelevant entries
      call decod(Icase,2,NR_PHYSA, ncase)
!
      ivarH=0 ; nvarH=0

!  ...loop through multiple copies of variables
      do iload=1,NRCOMS
!
!  ......loop through physical attributes
         do iattr=1,NR_PHYSA
!
!  .........loop through components
            do icomp=1,NR_COMP(iattr)
               select case(DTYPE(iattr))
               case('contin')
                  ivarH=ivarH+1
                  if (ncase(iattr) == 1) then
                     nvarH=nvarH+1
                     ZdofH(nvarH) = zvalH(ivarH)
                  endif
               endselect
!
!  .........loop through components
            enddo
!  ......loop through physical attributes
         enddo
!  ...loop through multiple copies of variables
      enddo
!
!
endsubroutine dhpvert
