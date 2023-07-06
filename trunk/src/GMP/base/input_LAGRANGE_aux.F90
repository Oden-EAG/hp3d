!----------------------------------------------------------------------
!> Purpose : routine creates a list of blocks adjacent to a curve 
!            specified by its endpoints
!  
!! @param[in ] Np1,Np2       - point numbers
!! @param[in ] Np1_blk       - list of blocks connected to Np1
!! @param[in ] Nr1           - number of blocks connected to Np1
!! @param[in ] Np2_blk       - list of blocks connected to Np1
!! @param[in ] Nr2           - number of blocks connected to Np2
!! @param[in ] Nr            - dimension of List_blk
!
!! @param[out] List_blk      - list of blocks connected to curve
!                              Np1 - Np2 and the corresponding
!                              edge numbers
!! @param[out] Nr_blk        - number of blocks connected to the
!                              curve
!
!! @revision Jul 23
!---------------------------------------------------------------------
!
      subroutine find_curve_blocks(Np1,Nr1,Np1_blk, Np2,Nr2,Np2_blk, &
                                   List_blk,Nr,Nr_blk)
!
      use GMP
      use element_data
      implicit none
!
      integer,  intent(in)  :: Np1,Np2,Nr1,Nr2,Nr
      integer,  intent(in)  :: Np1_blk(Nr1),Np2_blk(Nr2)
      integer,  intent(out) :: List_blk(2,Nr)
      integer,  intent(out) :: Nr_blk
!
      integer :: iprint,ibl,k1,k2,nbl,lab,ntype,ie,iv(1:2),mp(1:2),loc1,loc2,j
!
      iprint=0  
      ibl=0
!
!  ...loop through the blocks connected to the first endpoint 
      do k1=1,Nr1
!
!  .....loop through the blocks connected to the second endpoint 
        do k2=1,Nr2
!
!  .......if the block is on both point lists then add it to the curve list
          if (Np1_blk(k1).eq.Np2_blk(k2)) then
            ibl=ibl+1
            if (ibl.gt.Nr) then
              write(*,*) 'find_curve_blocks: INCREASE Nr'
              stop 1
            endif
            List_blk(1,ibl) = Np1_blk(k1)
            call decode(Np1_blk(k1), nbl,lab)
            select case(lab)
            case(1); ntype = PRIS
            case(2); ntype = BRIC
            case(3); ntype = TETR
            case(4); ntype = PYRA
            end select
            do ie=1,nedge(ntype)
              call edge_to_vert(ntype,ie, iv(1),iv(2))
              select case(ntype)
              case(PRIS); mp(1:2) = PRISMS(nbl)%VertNo(iv(1:2))
              case(BRIC); mp(1:2) = HEXAS(nbl)%VertNo(iv(1:2))
              case(TETR); mp(1:2) = TETRAS(nbl)%VertNo(iv(1:2))
              case(PYRA); mp(1:2) = PYRAMIDS(nbl)%VertNo(iv(1:2))
              end select
              call locate(Np1,mp,2, loc1)
              if (loc1.ne.0) then
                call locate(Np2,mp,2, loc2)
                if (loc2.ne.0) then
                  List_blk(2,ibl) = ie
                endif
              endif
            enddo
          endif
!
!  .....end of loop through blocks connected to the second point
        enddo
!
!  ...end of loop through blocks connected to the first point
      enddo
      Nr_blk = ibl
      if (iprint.eq.1) then
        write(*,7010) Np1,Np2
 7010   format('find_curve_blocks: List_blk for Np1,Np2 = ',2i5)
        write(*,9010) (List_blk(1:2,j),j=1,Nr_blk)
 9010   format(10(i8,1x,i2,3x))
        call pause
      endif
!
!
      end subroutine find_curve_blocks
!
!----------------------------------------------------------------------
!> Purpose : routine creates a list of blocks adjacent to a figure specified
!            by its endpoints
!  
!! @param[in ] Np1,Np2,Np3,Np4  - point numbers (Np4=Np1 for a triangle)
!! @param[in ] Np1_blk          - list of blocks connected to Np1
!! @param[in ] Nr1              - number of blocks connected to Np1
!!             Same for the remaining points
!
!! @param[out] List_blk         - list of blocks connected to curve
!                                 Np1-Np2-Np3(-Np4) and the corresponding
!                                 face numbers
!! @param[out] Nr_blk           - number of blocks connected to the
!                                 figure
!
!! @revision Jul 23
!---------------------------------------------------------------------
!
      subroutine find_figure_blocks(Ftype, &
                                    Np1,Nr1,Np1_blk, Np2,Nr2,Np2_blk, &
                                    Np3,Nr3,Np3_blk, &
                                    List_blk,Nr_blk)
!
      use GMP
      use element_data
      implicit none
!
      integer,  intent(in)  :: Ftype
      integer,  intent(in)  :: Np1,Np2,Np3,Nr1,Nr2,Nr3
      integer,  intent(in)  :: Np1_blk(Nr1),Np2_blk(Nr2),Np3_blk(Nr3)
      integer,  intent(out) :: List_blk(2,2)
      integer,  intent(out) :: Nr_blk
!
      integer :: iprint,ibl,k1,k2,nbl1,loc,nbl,lab,ntype,jf,&
                 iv(4),mp(4),loc1,loc2,loc3,j
!
      iprint=0  
      ibl=0
!
!  ...loop through blocks connected to Np1
      do k1=1,Nr1
        nbl1 = Np1_blk(k1)
!
!  .....loop through blocks connected to Np2
        do k2=1,Nr2
!
!  .......for a block connected to Np1 and Np2
          if (Np2_blk(k2).eq.nbl1) then
!
!  .........look for the common block on the list of blocks connected to Np3
            call locate(nbl1,Np3_blk,Nr3, loc)
            if (loc.ne.0) then
!
!  ...........add the block to the list
              ibl=ibl+1
              List_blk(1,ibl) = nbl1
!
!  ...........find number of faces
              call decode(nbl1, nbl,lab)
              select case(lab)
              case(1); ntype = PRIS
              case(2); ntype = BRIC
              case(3); ntype = TETR
              case(4); ntype = PYRA
              end select
!
!  ...........loop through faces
              do jf=1,nface(ntype)
!
!  .............find the face points 
                call face_to_vert(ntype,jf,iv(1),iv(2),iv(3),iv(4))
                select case(lab)
                case(1); mp(1:4) = PRISMS(nbl)%VertNo(iv(1:4))
                case(2); mp(1:4) = HEXAS(nbl)%VertNo(iv(1:4))
                case(3); mp(1:4) = TETRAS(nbl)%VertNo(iv(1:4))
                case(4); mp(1:4) = PYRAMIDS(nbl)%VertNo(iv(1:4))
                end select
                call locate(Np1,mp,4, loc1)
                if (loc1.ne.0) then
                  call locate(Np2,mp,4, loc2)
                  if (loc2.ne.0) then
                    call locate(Np3,mp,4, loc3)
                    if (loc3.ne.0) then
!
!  ...................save the face number
                      List_blk(2,ibl) = jf
                    endif
                  endif
                endif
!
!  ...........end of loop through faces
              enddo
            endif
          endif  
!
!  .....end of loop through blocks attached to Np2
        enddo
!
!  ...end of loop through blocks attached to Np1
      enddo
!
      Nr_blk = ibl
      if (iprint.eq.1) then
        write(*,7010) Np1,Np2,Np3
 7010   format('find_figure_blocks: List_blk for Np1,Np2,Np3 = ',3i5)
        write(*,9010) (List_blk(1:2,j),j=1,Nr_blk)
 9010   format(10(i8,1x,i2,3x))
        call pause
      endif
!
      end subroutine find_figure_blocks

        
