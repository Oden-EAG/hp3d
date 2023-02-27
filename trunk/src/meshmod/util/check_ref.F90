!----------------------------------------------------------------------------
!> @brief      check compatability two refinements
!!
!> @param[in]  Type  - face type
!> @param[in]  Kref1 - existing refinement
!> @param[in]  Kref2 - trial refinement
!> @param[out] Ipass - compatability flag ( 0 - cannot upgrade )
!!
!> @date       Feb 2023
!----------------------------------------------------------------------------
subroutine check_ref(Ntype,Kref1,Kref2, Ipass)
!
      use node_types
      implicit none
      integer, intent(in)  :: Ntype
      integer, intent(in)  :: Kref1,Kref2
      integer, intent(out) :: Ipass
      integer :: krefh1,krefv1,krefh2,krefv2
      integer :: krefx1,krefx2,krefy1,krefy2,krefz1,krefz2
!----------------------------------------------------------------------------
!
!  ...initialize to compatible refinements
      Ipass=1
!
      select case(Ntype)
!
!============================================================================
!  Triangle                                                                 |
!============================================================================
      case(MDLT)
        if ( (Kref1.ne.0).and.(Kref1.ne.Kref2) ) then
          Ipass=0
        endif
!
!============================================================================
!  Quad                                                                     |
!============================================================================
      case(MDLQ)
        call decode(Kref1, krefh1,krefv1)
        call decode(Kref2, krefh2,krefv2)
        if (krefh1.gt.krefh2) Ipass=0
        if (krefv1.gt.krefv2) Ipass=0
!
!============================================================================
!  Tet                                                                      |
!============================================================================
      case(MDLN)
        if (Kref1.ne.0) then
          if (Kref1.ne.Kref2) Ipass=0
        endif
!
!============================================================================
!  Prism / Pyramid                                                          |
!============================================================================
!
      case(MDLP,MDLD)
        call decode(Kref1, krefh1,krefv1)
        call decode(Kref2, krefh2,krefv2)
        if (krefh1.gt.krefh2) Ipass=0
        if (krefv1.gt.krefv2) Ipass=0
!
!============================================================================
!  Hexa                                                                     |
!============================================================================
      case(MDLB)
        call decode(Kref1,  krefh1,krefz1)
        call decode(krefh1, krefx1,krefy1)
!
        call decode(Kref2,  krefh2,krefz2)
        call decode(krefh2, krefx2,krefy2)
        if (krefx1.gt.krefx2) Ipass=0
        if (krefy1.gt.krefy2) Ipass=0
        if (krefz1.gt.krefz2) Ipass=0
!
      case default
!
      endselect
!
!
end subroutine check_ref
