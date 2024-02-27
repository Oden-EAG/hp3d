!----------------------------------------------------------------------
!
!   routine name       - find_curve_to_block
!
!----------------------------------------------------------------------
!
!   latest revision    - Mar 2023
!
!   purpose            - routine finds blocks adjacent to a curve
!
!   arguments :
!     in:
!               Nc     - a curve
!               Maxbl  - dimension of Neigbl
!                        (anticipated max number of adjacent blocks)
!     out:
!               Nrbl   - number of adjacent blocks
!               Neigbl - list of the blocks nicknames
!
!----------------------------------------------------------------------
!
subroutine find_curve_to_block(Nc,Maxbl, Nrbl,Neigbl)
!
      use GMP
!
      implicit none
!
      integer :: Nc,Maxbl,Nrbl
      integer :: Neigbl(Maxbl)
!
      integer :: if,is,lab,nf,num
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
!
      if (iprint.eq.1) write(*,7001) Nc
 7001 format('find_curve_to_block: Nc = ',i8)
#endif
!
!  ...initiate number of neighboring blocks
      Nrbl=0
!
!  ...loop through figures adjacent to the curve
      do if=1,CURVES(Nc)%NrFig
        call decode(abs(CURVES(Nc)%FigNo(if)), nf,lab)
#if HP3D_DEBUG
        if (iprint.eq.1) then
          write(*,7002) if,CURVES(Nc)%FigNo(if),nf,lab
 7002     format('if,CURVES(Nc)%FigNo(if),nf,lab = ',i2,i8,2x,i8,i2)
        endif
#endif
        select case(lab)
!
!  .....triangle
        case(1)
!
!  .......loop through blocks adjacent to the triangle
          do is=1,2
            if (TRIANGLES(nf)%BlockNo(is).eq.0) cycle
            call locate(TRIANGLES(nf)%BlockNo(is),Neigbl,Nrbl, num)
            if (num.eq.0) then
              Nrbl=Nrbl+1
              if (Nrbl.gt.Maxbl) then
                write(*,7010)
 7010           format('find_curve_to_block: ', &
                       'NUMBER OF NEIGHBORS EXCEEDED')
                stop 1
              endif
              Neigbl(Nrbl) = TRIANGLES(nf)%BlockNo(is)
            endif
          enddo
!
!  .....rectangle
        case(2)
!
!  .......loop through blocks adjacent to the rectangle
          do is=1,2
            if (RECTANGLES(nf)%BlockNo(is).eq.0) cycle
            call locate(RECTANGLES(nf)%BlockNo(is),Neigbl,Nrbl, num)
            if (num.eq.0) then
              Nrbl=Nrbl+1
              if (Nrbl.gt.Maxbl) then
                write(*,7010)
                stop 1
              endif
              Neigbl(Nrbl) = RECTANGLES(nf)%BlockNo(is)
            endif
          enddo
        end select
      enddo
#if HP3D_DEBUG
      if (iprint.eq.1) call pause
#endif
!
!
end subroutine find_curve_to_block
