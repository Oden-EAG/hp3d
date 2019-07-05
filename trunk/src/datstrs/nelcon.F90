!----------------------------------------------------------------------
!
!   routine name       - nelcon
!
!----------------------------------------------------------------------
!
!   latest revision    - Dec 07
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
      iprint=0
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
        call nr_mdle_sons(NODES(nfath)%type,NODES(nfath)%ref_kind,nrbros)
        call locate(mdle,NODES(nfath)%sons,nrbros, noson)
        if (iprint.eq.1) then
          write(*,7002) mdle,nfath,nrbros,noson
 7002     format('nelcon: mdle,nfath,nrbros,noson = ',2i7,2i3)
        endif
!
!  .....if mdle is not the last son in the family, go to the next brother
        if (noson.lt.nrbros) then
          mdle = NODES(nfath)%sons(noson+1)
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
        mdle = NODES(mdle)%sons(1)
      enddo
      Mdle1 = mdle
!
!
      if (iprint.eq.1) then
        write(*,7010) Mdle0,Mdle1
 7010   format('nelcon: Mdle0,Mdle1 = ',2i7)
      endif
!
end subroutine nelcon


!cccc      I like this routine more
!C$$$      subroutine nelcon(Mdle, Mdle_next)
!C$$$      integer, intent(in)  :: Mdle
!C$$$      integer, intent(out) :: Mdle_next
!
!C$$$! Step 0 : if Mdle = 0, set Mdle as the first element. Then start from Step 2.
!C$$$      iflag = 1;
!C$$$      if (Mdle.eq.0) then
!C$$$        Mdle = 1; iflag = 0;
!C$$$      end if
!
!C$$$! Step 1 : move right on the brothers' list until it meet a node that is not
!C$$$      a middle node.
!C$$$      do while ( true.and.iflag )
!C$$$        father = Mdle->father
!C$$$        if (father.lt.0) then   ! if Mdle is initial mesh element
!C$$$          Mdle = Mdle + 1       ! move to the right
!C$$$          exit                  ! loop exit
!C$$$        else
!C$$$          ison  = locate(Mdle, father->sons) ! otherwise, find next brother
!C$$$          nbros = nr_midle_sons(Mdle)
!C$$$          if (ison.lt.nbros) then ! if next brother is middle node
!C$$$            Mdle = father->sons( ison + 1 ) ! Mdle = next brother
!C$$$            exit                ! loop exit
!C$$$          end if
!C$$$        end if
!C$$$        Mdle = father           ! move up on the tree
!C$$$      end if
!
!C$$$! Step 2 : move down to first son until it reach the leaf level middle node.
!C$$$      do while ( .not.Is_leaf(Mdle) )
!C$$$        Mdle = Mdle->sons( 1 )
!C$$$      end do
!
!C$$$! Assign next middle node
!C$$$      Mdle_next = Mdle
!
!C$$$      end subroutine nelcon


