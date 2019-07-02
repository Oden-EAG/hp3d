!--------------------------------------------------------------------
!> Purpose : routine performs a global p-enrichment
!!
!> @date Aug 17
!--------------------------------------------------------------------
!

!--------------------------------------------------------------------
!
   subroutine global_pref

   use data_structure3D , only : NRELES
   use data_structure3D , only : NRNODS,NODES
!      
   implicit none
   integer :: mdle,iel,nord
!
!  ...loop over active elements
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        nord = NODES(mdle)%order
        select case(NODES(mdle)%type)
        case('mdln','mdld'); nord = nord+1
        case('mdlp'); nord = nord+11
        case('mdlb'); nord = nord+111
        end select
        call nodmod(mdle, nord)
      enddo

   end subroutine global_pref



   subroutine global_punref

   use data_structure3D , only : NRELES
   use data_structure3D , only : NRNODS,NODES
!      
   implicit none
   integer :: mdle,iel,nord
!
!  ...loop over active elements
      mdle=0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        nord = NODES(mdle)%order
        select case(NODES(mdle)%type)
        case('mdln','mdld'); nord = nord-1
        case('mdlp'); nord = nord-11
        case('mdlb'); nord = nord-111
        end select
        call nodmod(mdle, nord)
      enddo


   end subroutine global_punref