!---------------------------------------------------------------------------------
subroutine update_connectivities_ORIG_ENTITIES
!---------------------------------------------------------------------------------
! LATEST REVISION: Jul 09
!
! PURPOSE: routine updates connectivities to points for the
!          original mesh entities
!---------------------------------------------------------------------------------
  use SPLIT_SURF
!---------------------------------------------------------------------------------
  IMPLICIT NONE
!---------------------------------------------------------------------------------
  integer :: nc,ileft,iright,ion,iv,np,nt,nr,npri,ntet,npyr
!---------------------------------------------------------------------------------
! printing flag (0,1,2)
#define I_PRINT 0
!
#if I_PRINT >= 1
    write(*,*)'update_connectivities_ORIG_ENTITIES: updating connectivities...'
#endif
!---------------------------------------------------------------------------------
!  STEP 1: update CURVES --> POINTS connectivities
!
! ..loop through original curves
    do nc = 1, nrcurve_orig
      ileft = 0;  iright = 0;  ion = 0
      do iv = 1, 2
        np = CURVES(nc)%EndPoNo(iv)
! ......a twin point exists
        if (new_point(np) .gt. 0) then
          ion = ion + 1
! ......point on NEG side of split surface
        elseif (new_point(np) .eq. -1) then
          ileft = ileft + 1
! ......point on POS side of split surface
        elseif (new_point(np) .eq. -2) then
          iright = iright + 1
        endif
      enddo
      if ((ion .eq. 1) .and. (iright .eq. 1)) then
! ......loop over curve endpoints
        do iv = 1, 2
          np = CURVES(nc)%EndPoNo(iv)
! ........if endpoint has a twin point
          if (new_point(np) .gt. 0)  CURVES(nc)%EndPoNo(iv) = new_point(np)
        enddo
      endif
! ..end of loop through curves
    enddo
#if I_PRINT >= 2
    write(*,*)'update_connectivities: CURVES connectivities updated.'
#endif
!
!---------------------------------------------------------------------------------
!  STEP 2: update TRIANGLES --> POINTS connectivities
!
! ..loop through original triangles
    do nt = 1, nrtrian_orig
      ileft = 0;  iright = 0;  ion = 0
! ....loop through vertices
      do iv = 1, 3
        np = TRIANGLES(nt)%VertNo(iv)
! ......vertex was duplicated (INT)
        if ((new_point(np) .gt. 0) .and. (new_point(np) .ne. np)) then
          ion = ion + 1
! ......vertex on NEG_SIDE of split surface
        elseif (new_point(np) .eq. -1) then
          ileft = ileft + 1
! ......vertex on POS_SIDE of split surface
        elseif (new_point(np) .eq. -2) then
          iright = iright + 1
        endif
! ....end loop through vertices
      enddo
      if ((ion .gt. 0) .and. (iright .gt. 0) .and. (ileft .eq. 0)) then
        do iv = 1, 3
          np = TRIANGLES(nt)%VertNo(iv)
          if (new_point(np) .gt. 0)  TRIANGLES(nt)%VertNo(iv) = new_point(np)
        enddo
      endif
      if ((ion .gt. 0) .and. (iright .gt. 0) .and. (ileft .gt. 0)) then
        write(*,*)   'update_connectivities: triangle crossing split surface!'
        write(*,*)   '*******************  nt = ',nt
        do iv = 1, 3
          np = TRIANGLES(nt)%VertNo(iv)
          write(*,*) '**  iv,np,new_point(np) = ',iv,np,new_point(np)
        enddo
        call print_GMP
        stop
      endif
! ..end of loop through original triangles
    enddo
!
! ..loop through rectangles
      do nr=1,nrrecta_orig
        ileft=0; iright=0; ion=0
        do iv=1,4
          np = RECTANGLES(nr)%VertNo(iv)
          if ((new_point(np).gt.0).and.(new_point(np).ne.np)) then
            ion=ion+1
          elseif (new_point(np).eq.-1) then
            ileft=ileft+1
          elseif (new_point(np).eq.-2) then
            iright=iright+1
          endif
        enddo
        if ((ion.gt.0).and.(iright.gt.0).and.(ileft.eq.0)) then
          do iv=1,4
            np = RECTANGLES(nr)%VertNo(iv)
            if (new_point(np).gt.0)  RECTANGLES(nr)%VertNo(iv) = new_point(np)
          enddo
        endif
        if ((ion.gt.0).and.(iright.gt.0).and.(ileft.gt.0)) then
          write(*,*) 'split_surface: RECTANGLE CROSSING SPLIT SURFACE'
          write(*,*) 'nr = ',nr
          do iv=1,3
            np = RECTANGLES(nr)%VertNo(iv)
            write(*,*) 'iv,np,new_point(np) = ',iv,np,new_point(np)
          enddo
          call print_GMP
          stop
        endif
      enddo
!
!  ...loop through prisms
      do npri=1,nrprism_orig
#if I_PRINT >= 2
          write(*,*) 'split_surface: RECONNECTING PRISM ',npri
#endif
        ileft=0; iright=0; ion=0
        do iv=1,6
          np = PRISMS(npri)%VertNo(iv)
          if ((new_point(np).gt.0).and.(new_point(np).ne.np)) then
            ion=ion+1
          elseif (new_point(np).eq.-1) then
            ileft=ileft+1
          elseif (new_point(np).eq.-2) then
            iright=iright+1
          endif
        enddo
        if ((ion.gt.0).and.(iright.gt.0).and.(ileft.eq.0)) then
          do iv=1,6
            np = PRISMS(npri)%VertNo(iv)
            if (new_point(np).gt.0)  PRISMS(npri)%VertNo(iv) = new_point(np)
          enddo
        endif
        if ((ion.gt.0).and.(iright.gt.0).and.(ileft.gt.0)) then
          write(*,*) 'split_surface: PRISM CROSSING SPLIT SURFACE'
          stop
        endif
      enddo
!
!  ...loop through tets
      do ntet=1,nrtetra_orig
        ileft=0; iright=0; ion=0
        do iv=1,4
          np = TETRAS(ntet)%VertNo(iv)
          if ((new_point(np).gt.0).and.(new_point(np).ne.np)) then
            ion=ion+1
          elseif (new_point(np).eq.-1) then
            ileft=ileft+1
          elseif (new_point(np).eq.-2) then
            iright=iright+1
          endif
        enddo
        if ((ion.gt.0).and.(iright.gt.0).and.(ileft.eq.0)) then
          do iv=1,4
            np = TETRAS(ntet)%VertNo(iv)
            if (new_point(np).gt.0)  TETRAS(ntet)%VertNo(iv) = new_point(np)
          enddo
        endif
        if ((ion.gt.0).and.(iright.gt.0).and.(ileft.gt.0)) then
          write(*,*) 'split_surface: TET CROSSING SPLIT SURFACE'
          stop
        endif
      enddo
!
!  ...loop through pyramids
      do npyr=1,nrpyram_orig
        ileft=0; iright=0; ion=0
        do iv=1,5
          np = PYRAMIDS(npyr)%VertNo(iv)
          if ((new_point(np).gt.0).and.(new_point(np).ne.np)) then
            ion=ion+1
          elseif (new_point(np).eq.-1) then
            ileft=ileft+1
          elseif (new_point(np).eq.-2) then
            iright=iright+1
          endif
        enddo
        if ((ion.gt.0).and.(iright.gt.0).and.(ileft.eq.0)) then
          do iv=1,5
            np = PYRAMIDS(npyr)%VertNo(iv)
            if (new_point(np).gt.0) PYRAMIDS(npyr)%VertNo(iv) = new_point(np)
          enddo
        endif
        if ((ion.gt.0).and.(iright.gt.0).and.(ileft.gt.0)) then
          write(*,*) 'split_surface: PYRAMID CROSSING SPLIT SURFACE'
          stop
        endif
      enddo
#if I_PRINT >= 1
    write(*,*)'update_connectivities: done!'
#endif
!
end subroutine update_connectivities_ORIG_ENTITIES
