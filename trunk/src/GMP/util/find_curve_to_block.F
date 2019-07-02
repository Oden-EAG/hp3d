c----------------------------------------------------------------------
c
c   routine name       - find_curve_to_block
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Jan 08
c
c   purpose            - routine finds blocks adjacent to a curve
c
c   arguments :
c     in:
c               Nc     - a curve
c               Maxbl  - dimension of Neigbl
c                        (anticipated max number of adjacent blocks)
c     out:
c               Nrbl   - number of adjacent blocks
c               Neigbl - list of the blocks nicknames
c
c----------------------------------------------------------------------
c
      subroutine find_curve_to_block(Nc,Maxbl, Nrbl,Neigbl)
c
      use GMP
c
#include "syscom.blk"
c
      dimension Neigbl(Maxbl)
c
      iprint=0
      if (iprint.eq.1) write(*,7001) Nc
 7001 format('find_curve_to_block: Nc = ',i8)
c
c  ...initiate number of neighboring blocks
      Nrbl=0
c
c  ...loop through figures adjacent to the curve
      do if=1,CURVES(Nc)%NrFig
        call decode(iabs(CURVES(Nc)%FigNo(if)), nf,lab)
        if (iprint.eq.1) then
          write(*,7002) if,CURVES(Nc)%FigNo(if),nf,lab
 7002     format('if,CURVES(Nc)%FigNo(if),nf,lab = ',i2,i8,2x,i8,i2)
        endif
c
        select case(lab)
c
c  .....triangle
        case(1)
c
c  .......loop through blocks adjacent to the triangle
          do is=1,2
            if (TRIANGLES(nf)%BlockNo(is).eq.0) cycle
            call locate(TRIANGLES(nf)%BlockNo(is),Neigbl,Nrbl, num)
            if (num.eq.0) then
              Nrbl=Nrbl+1
              if (Nrbl.gt.Maxbl) then
                write(*,7010)
 7010           format('find_curve_to_block: ',
     .                 'NUMBER OF NEIGHBORS EXCEEDED')
                stop 1
              endif
              Neigbl(Nrbl) = TRIANGLES(nf)%BlockNo(is)
            endif
          enddo
c
c  .....rectangle
        case(2)
c
c  .......loop through blocks adjacent to the rectangle
          do is=1,2
            if (RECTANGLES(nf)%BlockNo(is).eq.0) cycle
            call locate(RECTANGLES(nf)%BlockNo(is),Neigbl,Nrbl, num)
            if (num.eq.0) then
              Nrbl=Nrbl+1
              if (Nrbl.gt.Maxbl) then
                write(*,7001)
                stop 1
              endif
              Neigbl(Nrbl) = RECTANGLES(nf)%BlockNo(is)
            endif
          enddo
        end select
      enddo
      if (iprint.eq.1) call pause
c
c
      end
