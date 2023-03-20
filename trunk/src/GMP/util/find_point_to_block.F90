!----------------------------------------------------------------------
!
!   routine name       - find_point_to_block
!
!----------------------------------------------------------------------
!
!   latest revision    - Mar 2023
!
!   purpose            - routine finds blocks adjacent to a point
!
!   arguments :
!     in:
!               Np     - a point
!               Maxbl  - dimension of Neigbl
!                        (anticipated max number of adjacent blocks)
!     out:
!               Nrbl   - number of adjacent blocks
!               Neigbl - list of the blocks nicknames
!
!----------------------------------------------------------------------
!
subroutine find_point_to_block(Np,Maxbl, Nrbl,Neigbl)
!
      use GMP
!
      implicit none
!
      integer :: Np,Maxbl,Nrbl
      integer :: Neigbl(Maxbl)
!
      integer :: ic,if,is,nc,nf,num,lab
!
#if DEBUG_MODE
      integer :: iprint
      iprint=0
#endif
!
!  ...initiate number of neighboring blocks
      Nrbl=0
!
!  ...loop through curves connected to the point
      do ic=1,POINTS(Np)%NrCurv
        nc = POINTS(Np)%CurvNo(ic)
!
!  .....loop through figures adjacent to the curve
        do if=1,CURVES(nc)%NrFig
          call decode(abs(CURVES(nc)%FigNo(if)), nf,lab)
!
          select case(lab)
!
!  .......triangle
          case(1)
!
!  .........loop through blocks adjacent to the triangle
            do is=1,2
              if (TRIANGLES(nf)%BlockNo(is).eq.0) cycle
              call locate(TRIANGLES(nf)%BlockNo(is),Neigbl,Nrbl, num)
              if (num.eq.0) then
                Nrbl=Nrbl+1
                if (Nrbl.gt.Maxbl) then
                  write(*,7001)
 7001             format('find_point_to_block: ', &
                         'NUMBER OF NEIGHBORS EXCEEDED')
                  stop 1
                endif
                Neigbl(Nrbl) = TRIANGLES(nf)%BlockNo(is)
              endif
            enddo
!
!  .......rectangle
          case(2)
!
!  .........loop through blocks adjacent to the rectangle
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
      enddo
!
#if DEBUG_MODE
      if (iprint.eq.1) then
        write(*,7002) Np, Neigbl(1:Nrbl)
 7002   format('find_point_to_block: Np, Neigbl = ',i5,10i6)
        call pause
      endif
#endif
!
!
end subroutine find_point_to_block
