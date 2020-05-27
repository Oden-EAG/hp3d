c----------------------------------------------------------------------
c
c   routine name       - find_element_ref_flag
c
c----------------------------------------------------------------------
c
c   latest revision    - May 10
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
c             Type     - middle node type
c             Kref     - refinement flag for the middle node of
c                        the element
c             Kreff    - refinement flags for the faces
c     out:
c             Krefm    - modified refinement flag for the middle node
c
c----------------------------------------------------------------------
c
      subroutine find_element_ref_flag(Type,Kref,Kreff, Krefm)
c
      use element_data
#include "syscom.blk"
c
      character(len=4) :: Type
c
      dimension Kreff(6)
      dimension kreff_loc(6)
c
c  ...possible element refinement flags
      integer :: mdln_ref_flags(1:13) = (/11,12,13,
     .                                   21,22,23,24,25,26,
     .                                   31,32,33,34/)
      integer :: mdlp_ref_flags(1:3) = (/11,10,01/)
      integer :: mdld_ref_flags(1:2) = (/10,01/)
      integer :: mdlb_ref_flags(1:7) = (/111,110,101,011,100,010,001/)
c
      iprint=0
c
      Krefm=-1
      select case(Type)
      case('mdln')
c
c  .....loop through possible refinements
        do i=1,13
          kref_loc = mdln_ref_flags(i)
          call check_ref(Type,Kref,kref_loc, ipass)
          if (ipass.eq.0) cycle
          call find_face_ref_flags(Type,kref_loc, kreff_loc)
          do if=1,4
            call check_ref('mdlt',Kreff(if),kreff_loc(if), ipass)
          enddo
          if (ipass.eq.1) then
            Krefm = kref_loc; exit
          endif
        enddo
      case('mdlp')
        do i=3,1,-1
          kref_loc = mdlp_ref_flags(i)
          call check_ref(Type,Kref,kref_loc, ipass)
          if (ipass.eq.0) cycle
          call find_face_ref_flags(Type,kref_loc, kreff_loc)
          do if=1,2
            call check_ref('mdlt',Kreff(if),kreff_loc(if), ipass)
            if (iprint.eq.1) then
              write(*,8001) if,Kreff(if),kreff_loc(if),ipass
 8001         format('if,Kreff(if),kreff_loc(if),ipass = ',4i4)
            endif
          enddo
          do if=3,5
            call check_ref('mdlq',Kreff(if),kreff_loc(if), ipass)
            if (iprint.eq.1) then
              write(*,8001) if,Kreff(if),kreff_loc(if),ipass
            endif
          enddo
          if (ipass.eq.1) then
            Krefm = kref_loc; exit
          endif
        enddo
      case('mdld')
        do i=1,2
          kref_loc = mdld_ref_flags(i)
          call check_ref(Type,Kref,kref_loc, ipass)
          if (ipass.eq.0) cycle
          call find_face_ref_flags(Type,kref_loc, kreff_loc)
          call check_ref('mdlq',Kreff(1),kreff_loc(1), ipass)
          do if=2,5
            call check_ref('mdlt',Kreff(if),kreff_loc(if), ipass)
          enddo
          if (ipass.eq.1) then
            Krefm = kref_loc; exit
          endif
        enddo
      case('mdlb')
        do i=7,1,-1
          kref_loc = mdlb_ref_flags(i)
          call check_ref(Type,Kref,kref_loc, ipass)
          if (ipass.eq.0) cycle
          call find_face_ref_flags(Type,kref_loc, kreff_loc)
          do if=1,6
            call check_ref('mdlq',Kreff(if),kreff_loc(if), ipass)
          enddo
          if (ipass.eq.1) then
            Krefm = kref_loc; exit
          endif
        enddo
      end select
c
      if (iprint.eq.1) then
        write(*,7001) Type,Kref,Krefm,Kreff(1:nface(Type))
 7001   format('find_element_ref_flag: Type,Kref,Krefm,Kreff = ',
     .          a5,2x,i2,2x,i2,6x,6(2x,i2))
        call pause
      endif
c
c
      end subroutine
c
c----------------------------------------------------------------------
c
c   routine name       - find_element_ref_flag2
c
c----------------------------------------------------------------------
c
c   latest revision    - May 10
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
c             Type     - middle node type
c             Kreff    - refinement flags for faces
c             Krefe    - refinement flags for edges
c     out:
c             Kref     - refinement flag for the middle node
c
c----------------------------------------------------------------------
c
      subroutine find_element_ref_flag2(Type,Kreff,Krefe, Kref)
c
      use element_data
#include "syscom.blk"
c
      character(len=4) :: Type
c
      dimension Kreff(6),Krefe(12)
      dimension kreff_loc(6),krefe_loc(12)
c
c  ...possible element refinement flags
      integer :: mdln_ref_flags(1:13) = (/11,12,13,
     .                                   21,22,23,24,25,26,
     .                                   31,32,33,34/)
      integer :: mdlp_ref_flags(1:3) = (/11,10,01/)
      integer :: mdld_ref_flags(1:2) = (/10,01/)
      integer :: mdlb_ref_flags(1:7) = (/111,110,101,011,100,010,001/)
c
      iprint=0
c
c  ...determine number of refinements to check
      select case(Type)
      case('mdln'); nr_ref=13
      case('mdlp'); nr_ref=3
      case('mdld'); nr_ref=2
      case('mdlb'); nr_ref=7
      end select
c
c  ...loop through possible refinements
      do i=nr_ref,1,-1
        ipass=1
        select case(Type)
        case('mdln'); kref_loc = mdln_ref_flags(i)
        case('mdlp'); kref_loc = mdlp_ref_flags(i)
        case('mdld'); kref_loc = mdld_ref_flags(i)
        case('mdlb'); kref_loc = mdlb_ref_flags(i)
        end select
c
        call find_face_ref_flags(Type,kref_loc, kreff_loc)
        call find_edge_ref_flags(Type,kref_loc, krefe_loc)
        do if=1,nface(Type)
          call check_ref('mdlt',Kreff(if),kreff_loc(if), ipass)
        enddo
        do ie=1,nedge(Type)
          if (krefe_loc(ie).lt.Krefe(ie)) ipass=0
        enddo
        if (ipass.eq.1) then
          Kref = kref_loc; exit
        endif
      enddo
c
      if (iprint.eq.1) then
        write(*,7001) Type
 7001   format('find_element_ref_flag2: Type = ',a5)
        write(*,7002) Kreff(1:nface(Type))
 7002   format('                        Face flags = ',6i3)
        write(*,7003) Krefe(1:nedge(Type))
 7003   format('                        Edge flags = ',6i3)
        write(*,7004) Kref
 7004   format('                        Element flag = ',i3)
      endif
c
c
      end subroutine










