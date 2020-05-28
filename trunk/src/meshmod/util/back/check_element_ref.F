c----------------------------------------------------------------------
c
c   routine name       - check_ref
c
c----------------------------------------------------------------------
c
c   latest revision    - May 10
c
c   purpose            - given two refinement flags for a middle node
c                        routine checks if the second refinement
c                        accomodates the first one
c
c
c   arguments :
c     in:
c             Type     - middle node type
c             Kref1    - first refinement flag
c             Kref2    - second refinement flag, possible upgrade
c                        of Kref1
c     out:
c             Ipass    = 0 if the update is impossible
c                      = 1 if the update is possible
c
c----------------------------------------------------------------------
c
      subroutine check_ref(Type,Kref1,Kref2, Ipass)
c
#include "syscom.blk"
c
      character(len=4) :: Type
c
      iprint=0
c
      Ipass=1
      select case(Type)
      case('mdln')
        if (Kref1.ne.0) then
          if (Kref1.ne.Kref2) Ipass=0
        endif
      case('mdlp','mdld')
        call decode(Kref1, krefh1,krefv1)
        call decode(Kref2, krefh2,krefv2)
        if (krefh1.gt.krefh2) Ipass=0
        if (krefv1.gt.krefv2) Ipass=0
      case('mdlb')
        call decode(Kref1, krefh1,krefz1)
        call decode(krefh1, krefx1,krefy1)
        call decode(Kref2, krefh2,krefz2)
        call decode(krefh2, krefx2,krefy2)
        if (krefx1.gt.krefx2) Ipass=0
        if (krefy1.gt.krefy2) Ipass=0
        if (krefz1.gt.krefz2) Ipass=0
      end select
      if (iprint.eq.1) then
        write(*,7001) Kreff1,Kreff2,Ipass
 7001   format('check_ref: Kreff1,Kreff2,Ipass = ',3i3)
        call pause
      endif
c
c
      end subroutine









