!
#include "typedefs.h"
!
!-----------------------------------------------------------------------
!> Purpose :  routine determines Dirichlet dof for a vertex
!! @param[in]  Iflag - a flag specifying which of the objects the vertex
!! @param[in]  No    - number of a specific object
!! @param[in]  Xi    - reference coordinates of the point
!! @param[in]  Icase - node case
!! @param[in]  Bcond - node BC flags
!!
!! @param[out] ZdofH - updated dirichlet BC dof
!-----------------------------------------------------------------------
subroutine dhpvert(Mdle,Iflag,No,Xi,Icase,Bcond, ZdofH)
!
   use data_structure3D
   implicit none
!
!..Arguments
   integer, intent(in)  :: Iflag,No,Icase,Bcond,Mdle
   real(8), intent(in)  :: Xi(3)
   VTYPE  , intent(out) :: ZdofH(NRCOMS*NREQNH(Icase))
!
!..Locals
!..work space for dirichlet
   VTYPE :: zvalH(  MAXEQNH),zdvalH(  MAXEQNH,3), &
            zvalE(3,MAXEQNE),zdvalE(3,MAXEQNE,3), &
            zvalV(3,MAXEQNV),zdvalV(3,MAXEQNV,3)
!
   real(8), dimension(3)   :: x
   real(8), dimension(3,3) :: void
!
!..decimal representation of Icase
   integer :: ncase(NR_PHYSA)
!
!..decimal representation of Bcond
   integer :: ibcnd(NRINDEX)
!
   integer :: ivarH,nvarH,iphys,iload,icomp,ic
!
#if DEBUG_MODE
   integer :: iprint
   iprint=0
#endif
!
!---------------------------------------------------------------------
!
!..determine coordinates of the point
   select case(Iflag)
      case(5)                 ; call prism(No,Xi, x,void)
      case(6)                 ; call  hexa(No,Xi, x,void)
      case(7)                 ; call tetra(No,Xi, x,void)
      case(8)                 ; call pyram(No,Xi, x,void)
      case default
         write(*,7010) Iflag
 7010    format(' dhpvert: unknown Iflag = ',i7)
         stop
   end select
!
#if DEBUG_MODE
   if (iprint.eq.1) then
      write(*,7020) No, Xi, x
 7020 format('dhpvert: No = ',i3,' Xi = ',3f8.3,' x = ',3f8.3)
   endif
#endif
!
!..dirichlet value in compact mode
   call dirichlet(Mdle,x,Icase, zvalH,zdvalH,zvalE,zdvalE,zvalV,zdvalV)
!
!..shift the data skipping irrelevant entries
   call decod(Icase,2,NR_PHYSA, ncase)
   call decod(Bcond,2,NRINDEX,  ibcnd)
!
   ivarH=0; nvarH=0

!..loop through multiple copies of variables
   do iload=1,NRCOMS
!
!  ...initiate the BC component counter
      ic=0
!
!  ...loop through physical attributes
      do iphys=1,NR_PHYSA
!
!     ...loop through components
         do icomp=1,NR_COMP(iphys)
!
!        ...if the variable is supported by the node, update the BC component counter
            if (ncase(iphys).eq.1) ic=ic+1
!
!        ...select the discretization type
            select case(DTYPE(iphys))
               case('contin')
                  ivarH=ivarH+1
!
!              ...if the variable is supported by the node
                  if (ncase(iphys).eq.1) then
!
!                 ...update the H1 variable counter
                     nvarH=nvarH+1
!
!                 ...do not write dof if physics attribute is deactivated
                     if (.not. PHYSAm(iphys)) exit
!
!                 ...store Dirichlet dof
                     if (ibcnd(ic).eq.1) ZdofH(nvarH) = zvalH(ivarH)
                  endif
            end select
!
!     ...loop through components
         enddo
!
!  ...loop through physical attributes
      enddo
!
!..loop through multiple copies of variables
   enddo
!
#if DEBUG_MODE
   if (iprint.eq.1) then
      write(*,7030) ZdofH(1:NRCOMS*NREQNH(Icase))
 7030 format('dhpvert: ZdofH = ',10e12.5)
   endif
#endif
!
end subroutine dhpvert
