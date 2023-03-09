c----------------------------------------------------------------------
c
c   routine name       - find_point_to_block
c
c----------------------------------------------------------------------
c
c   latest revision    - Mar 2023
c
c   purpose            - routine finds blocks adjacent to a point
c
c   arguments :
c     in:
c               Np     - a point
c               Maxbl  - dimension of Neigbl
c                        (anticipated max number of adjacent blocks)
c     out:
c               Nrbl   - number of adjacent blocks
c               Neigbl - list of the blocks nicknames
c
c----------------------------------------------------------------------
c
      subroutine find_point_to_block(Np,Maxbl, Nrbl,Neigbl)
c
      use GMP
      implicit none
c
      integer :: Np,Maxbl,Nrbl
      integer :: Neigbl(Maxbl)
c
      integer :: ic,if,is,nc,nf,num,lab
c
      integer :: iprint=0
c
c  ...initiate number of neighboring blocks
      Nrbl=0
c
c  ...loop through curves connected to the point
      do ic=1,POINTS(Np)%NrCurv
        nc = POINTS(Np)%CurvNo(ic)
c
c  .....loop through figures adjacent to the curve
        do if=1,CURVES(nc)%NrFig
          call decode(abs(CURVES(nc)%FigNo(if)), nf,lab)
c
          select case(lab)
c
c  .......triangle
          case(1)
c
c  .........loop through blocks adjacent to the triangle
            do is=1,2
              if (TRIANGLES(nf)%BlockNo(is).eq.0) cycle
              call locate(TRIANGLES(nf)%BlockNo(is),Neigbl,Nrbl, num)
              if (num.eq.0) then
                Nrbl=Nrbl+1
                if (Nrbl.gt.Maxbl) then
                  write(*,7001)
 7001             format('find_point_to_block: ',
     .                   'NUMBER OF NEIGHBORS EXCEEDED')
                  stop 1
                endif
                Neigbl(Nrbl) = TRIANGLES(nf)%BlockNo(is)
              endif
            enddo
c
c  .......rectangle
          case(2)
c
c  .........loop through blocks adjacent to the rectangle
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
          endselect
        enddo
      enddo
c
      if (iprint.eq.1) then
        write(*,7002) Np, Neigbl(1:Nrbl)
 7002   format('find_point_to_block: Np, Neigbl = ',i5,10i6)
        call pause
      endif
c
c
      end subroutine find_point_to_block
