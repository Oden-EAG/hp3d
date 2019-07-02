!--------------------------------------------------------------------
!> Purpose - routine performs an element p-refinement (possibly a
!!           derefinement.)
!!
!! @param[in] Mdle - element middle node 
!! @param[in] Ip   - increment or decrement to x1-direction.
!! @param[in] Iq   - increment or decrement to x2-direction.
!! @param[in] Ir   - increment or decrement to x3-direction.
!!
!> @date Dec 14
!--------------------------------------------------------------------
! Remark: after calling this routine (possibly multiple times) you
! should call the following routines:
!
!   call enforce_min_rule
!   call update_gdof
!   call update_ddof
!--------------------------------------------------------------------
!
subroutine p_refine(Mdle,Ip,Iq,Ir)
!
      use data_structure3D , only : NODES,MAXP
!      
      implicit none
      integer, intent(in) :: Mdle,Ip,Iq,Ir
!
      character(len=4) :: ntype
      integer :: nord, nordh, jp, jq, jr
!
!--------------------------------------------------------------------
!
!     initialize (needed since ALL orders are incremented)
      jp=0 ; jq=0 ; jr=0
!
!     node type and order of approximation
      ntype=NODES(Mdle)%type ; nord=NODES(Mdle)%order
!
!     decode node current order of approximation
      select case(ntype)
      case('medg','mdlt','mdln','mdld') ; jp = nord
      case('mdlq','mdlp')               ; call decode(nord, jp,jq)
      case('mdlb')                      ; call decode(nord, nordh,jr)
                                          call decode(nordh, jp,jq)
      endselect
!
!     increment
      jp = jp + ip ; jq = jq + iq ; jr = jr + ir
!
!     check upper bound
      if (jp.gt.MAXP) jp = MAXP
      if (jq.gt.MAXP) jq = MAXP
      if (jr.gt.MAXP) jr = MAXP
!    
!     check lower bound
      if (jp.lt.1   ) jp = 1
      if (jq.lt.1   ) jq = 1
      if (jr.lt.1   ) jr = 1
!
!     encode new order of approximation
      select case (ntype)
      case('medg','mdlt','mdln','mdld'); nord = jp
      case('mdlq','mdlp');               nord = 10*jp+jq
      case('mdlb');                      nord = 100*jp+10*jq+jr
      end select
!
!     modify node order of approximation (and everything related)
      call nodmod(Mdle, nord)
!
!
endsubroutine p_refine
