c----------------------------------------------------------------------
c
c   routine name       - find_element_ref_flag
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 2023
c
c   purpose            - given a middle node type and refinement flag
c                        'Kref', and refinement flags 'Kreff' for the
c                        element faces, routine checks the consistency
c                        of the refinements and updates the middle
c                        node refinement flag to accomodate face
c                        refinements
c
c
c   arguments :
c     in:
c             Ntype    - middle node type
c             Kref     - refinement flag for the middle node of
c                        the element
c             Kreff    - refinement flags for the faces
c     out:
c             Krefm    - modified refinement flag for the middle node
c
c----------------------------------------------------------------------
c
      subroutine find_element_ref_flag(Ntype,Kref,Kreff, Krefm)
c
      use element_data
c
      implicit none
c
      integer, intent(in)  :: Ntype,Kref
      integer, intent(in)  :: Kreff(6)
      integer, intent(out) :: Krefm
c
      integer :: kreff_loc(6)
c
c  ...possible element refinement flags
      integer :: mdln_ref_flags(1:13) = (/11,12,13,
     .                                   21,22,23,24,25,26,
     .                                   31,32,33,34/)
      integer :: mdlp_ref_flags(1:3) = (/11,10,01/)
      integer :: mdld_ref_flags(1:2) = (/10,01/)
      integer :: mdlb_ref_flags(1:7) = (/111,110,101,011,100,010,001/)
c
      integer :: i,ifc,ipass,kref_loc
c
#if DEBUG_MODE
      integer :: iprint
      iprint=0
#endif
c
      Krefm=-1
      select case(Ntype)
      case(MDLB)
c  .....loop through possible refinements
        do i=7,1,-1
          kref_loc = mdlb_ref_flags(i)
          call check_ref(Ntype,Kref,kref_loc, ipass)
          if (ipass.eq.0) cycle
          call find_face_ref_flags(Ntype,kref_loc, kreff_loc)
          do ifc=1,6
            call check_ref(MDLQ,Kreff(ifc),kreff_loc(ifc), ipass)
          enddo
          if (ipass.eq.1) then
            Krefm = kref_loc; exit
          endif
        enddo
      case(MDLN)
        do i=1,13
          kref_loc = mdln_ref_flags(i)
          call check_ref(Ntype,Kref,kref_loc, ipass)
          if (ipass.eq.0) cycle
          call find_face_ref_flags(Ntype,kref_loc, kreff_loc)
          do ifc=1,4
            call check_ref(MDLT,Kreff(ifc),kreff_loc(ifc), ipass)
          enddo
          if (ipass.eq.1) then
            Krefm = kref_loc; exit
          endif
        enddo
      case(MDLP)
        do i=3,1,-1
          kref_loc = mdlp_ref_flags(i)
          call check_ref(Ntype,Kref,kref_loc, ipass)
          if (ipass.eq.0) cycle
          call find_face_ref_flags(Ntype,kref_loc, kreff_loc)
          do ifc=1,2
            call check_ref(MDLT,Kreff(ifc),kreff_loc(ifc), ipass)
#if DEBUG_MODE
            if (iprint.eq.1) then
              write(*,8001) ifc,Kreff(ifc),kreff_loc(ifc),ipass
 8001         format('ifc,Kreff(ifc),kreff_loc(ifc),ipass = ',4i4)
            endif
#endif
          enddo
          do ifc=3,5
            call check_ref(MDLQ,Kreff(ifc),kreff_loc(ifc), ipass)
#if DEBUG_MODE
            if (iprint.eq.1) then
              write(*,8001) ifc,Kreff(ifc),kreff_loc(ifc),ipass
            endif
#endif
          enddo
          if (ipass.eq.1) then
            Krefm = kref_loc; exit
          endif
        enddo
      case(MDLD)
        do i=1,2
          kref_loc = mdld_ref_flags(i)
          call check_ref(Ntype,Kref,kref_loc, ipass)
          if (ipass.eq.0) cycle
          call find_face_ref_flags(Ntype,kref_loc, kreff_loc)
          call check_ref(MDLQ,Kreff(1),kreff_loc(1), ipass)
          do ifc=2,5
            call check_ref(MDLT,Kreff(ifc),kreff_loc(ifc), ipass)
          enddo
          if (ipass.eq.1) then
            Krefm = kref_loc; exit
          endif
        enddo
      end select
c
#if DEBUG_MODE
      if (iprint.eq.1) then
        write(*,7001) S_Type(Ntype),Kref,Krefm,Kreff(1:nface(Ntype))
 7001   format('find_element_ref_flag: Type,Kref,Krefm,Kreff = ',
     .          a5,2x,i2,2x,i2,6x,6(2x,i2))
        call pause
      endif
#endif
c
      end subroutine find_element_ref_flag
c
c----------------------------------------------------------------------
c
c   routine name       - find_element_ref_flag2
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 2023
c
c   purpose            - given an element  middle node type, refinement
c                        flags 'Kreff' for the element faces, and
c                        refinement flags 'Krefe' for the element edges,
c                        routine determines a minimum refinement to
c                        accomodate the requested refinements of faces
c                        and edges
c
c
c   arguments :
c     in:
c             Ntype    - middle node type
c             Kreff    - refinement flags for faces
c             Krefe    - refinement flags for edges
c     out:
c             Kref     - refinement flag for the middle node
c
c----------------------------------------------------------------------
c
      subroutine find_element_ref_flag2(Ntype,Kreff,Krefe, Kref)
c
      use element_data
c
      implicit none
c
      integer, intent(in)  :: Ntype
      integer, intent(in)  :: Kreff(6),Krefe(12)
      integer, intent(out) :: Kref
c
      integer :: kreff_loc(6),krefe_loc(12)
c
c  ...possible element refinement flags
      integer :: mdln_ref_flags(1:13) = (/11,12,13,
     .                                   21,22,23,24,25,26,
     .                                   31,32,33,34/)
      integer :: mdlp_ref_flags(1:3) = (/11,10,01/)
      integer :: mdld_ref_flags(1:2) = (/10,01/)
      integer :: mdlb_ref_flags(1:7) = (/111,110,101,011,100,010,001/)
c
      integer :: i,ifc,ipass,kref_loc,nr_ref
c
#if DEBUG_MODE
      integer :: iprint
      iprint=0
#endif
c
c  ...determine number of refinements to check
      select case(Ntype)
        case(MDLB); nr_ref=7
        case(MDLN); nr_ref=13
        case(MDLP); nr_ref=3
        case(MDLD); nr_ref=2
      end select
c
c  ...loop through possible refinements
      do i=nr_ref,1,-1
        ipass=1
        select case(Ntype)
          case(MDLB); kref_loc = mdlb_ref_flags(i)
          case(MDLN); kref_loc = mdln_ref_flags(i)
          case(MDLP); kref_loc = mdlp_ref_flags(i)
          case(MDLD); kref_loc = mdld_ref_flags(i)
        end select
c
        call find_face_ref_flags(Ntype,kref_loc, kreff_loc)
        call find_edge_ref_flags(Ntype,kref_loc, krefe_loc)
        do ifc=1,nface(Ntype)
          call check_ref(MDLT,Kreff(ifc),kreff_loc(ifc), ipass)
        enddo
        do ie=1,nedge(Ntype)
          if (krefe_loc(ie).lt.Krefe(ie)) ipass=0
        enddo
        if (ipass.eq.1) then
          Kref = kref_loc; exit
        endif
      enddo
c
#if DEBUG_MODE
      if (iprint.eq.1) then
        write(*,7001) S_Type(Ntype)
 7001   format('find_element_ref_flag2: Type = ',a5)
        write(*,7002) Kreff(1:nface(Ntype))
 7002   format('                        Face flags = ',6i3)
        write(*,7003) Krefe(1:nedge(Ntype))
 7003   format('                        Edge flags = ',6i3)
        write(*,7004) Kref
 7004   format('                        Element flag = ',i3)
      endif
#endif
c
      end subroutine find_element_ref_flag2
