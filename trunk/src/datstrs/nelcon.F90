!----------------------------------------------------------------------
!
!   routine name       - nelcon
!
!----------------------------------------------------------------------
!
!   latest revision    - Feb 2023
!
!   purpose            - routine enables the natural order of elements
!
!   arguments
!     in:
!                Mdle0 - (middle node of) element
!     out:
!                Mdle1 - (middle node of) element next in the
!                        natural order of elements
!
!----------------------------------------------------------------------
!
subroutine nelcon(Mdle0, Mdle1)
!
      use data_structure3D
      use refinements
!
      implicit none
!
      integer, intent(in)  :: Mdle0
      integer, intent(out) :: Mdle1
!
      integer :: mdle,nfath,nrbros,noson
!
#if DEBUG_MODE
      integer :: iprint
      iprint=0
#endif
!
      mdle = Mdle0
!
!  ...first element in the mesh
      if (mdle.eq.0) then
        mdle = 1
        go to 20
      endif
!
!  ...Step 1: move horizontally, if you can, otherwise move vertically up
   10 continue
      nfath = NODES(mdle)%father
!
!  ...initial mesh element
      if (nfath.lt.0) then
!
!  .....move to the next initial mesh element
        mdle = mdle + 1
        go to 20
      else
!
!  .....find the son number in the family
        call nr_mdle_sons(NODES(nfath)%ntype,NODES(nfath)%ref_kind,nrbros)
!        call locate(mdle,NODES(nfath)%sons,nrbros, noson)
        noson = mdle - NODES(nfath)%first_son + 1
!        if (noson<0 .or. noson>nrbros) call pause
!
#if DEBUG_MODE
        if (iprint.eq.1) then
          write(*,7002) mdle,nfath,nrbros,noson
 7002     format('nelcon: mdle,nfath,nrbros,noson = ',2i7,2i3)
        endif
#endif
!
!  .....if mdle is not the last son in the family, go to the next brother
        if (noson.lt.nrbros) then
!          mdle = NODES(nfath)%sons(noson+1)
          mdle = Son(nfath,noson+1)
          go to 20
        else
!
!  .......move up on the tree
          mdle = nfath
          go to 10
        endif
      endif
!
!  ...Step 2 move vertically down
   20 continue
      do while (NODES(mdle)%ref_kind.gt.0)
!        mdle = NODES(mdle)%sons(1)
        mdle = Son(mdle,1)
      enddo
      Mdle1 = mdle
!
!
#if DEBUG_MODE
      if (iprint.eq.1) then
        write(*,7010) Mdle0,Mdle1
 7010   format('nelcon: Mdle0,Mdle1 = ',2i7)
      endif
#endif
!
end subroutine nelcon


!! Alternate routine
!      subroutine nelcon(Mdle, Mdle_next)
!      integer, intent(in)  :: Mdle
!      integer, intent(out) :: Mdle_next
!
!! Step 0 : if Mdle = 0, set Mdle as the first element. Then start from Step 2.
!      iflag = 1;
!      if (Mdle.eq.0) then
!        Mdle = 1; iflag = 0;
!      end if
!
!! Step 1 : move right on the brothers' list until it meet a node that is not
!      a middle node.
!      do while ( true.and.iflag )
!        father = Mdle->father
!        if (father.lt.0) then   ! if Mdle is initial mesh element
!          Mdle = Mdle + 1       ! move to the right
!          exit                  ! loop exit
!        else
!          ison  = locate(Mdle, father->sons) ! otherwise, find next brother
!          nbros = nr_midle_sons(Mdle)
!          if (ison.lt.nbros) then ! if next brother is middle node
!            Mdle = father->sons( ison + 1 ) ! Mdle = next brother
!            exit                ! loop exit
!          end if
!        end if
!        Mdle = father           ! move up on the tree
!      end if
!
!! Step 2 : move down to first son until it reach the leaf level middle node.
!      do while ( .not.Is_leaf(Mdle) )
!        Mdle = Mdle->sons( 1 )
!      end do
!
!! Assign next middle node
!      Mdle_next = Mdle
!
!      end subroutine nelcon
