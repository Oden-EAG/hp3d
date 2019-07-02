!-----------------------------------------------------------------------------
!> Purpose : routine determines the ``natural order of elements for a grid
!            consiting of elements of generation level Ngen OR LOWER
!!
!! @param[in]  Ngen  - a generation level
!! @param[in]  Mdle0 - (midle node of an) element of generation level
!                      Ngen or lower IF THE NODE IS ACTIVE
!! @param[out] Mdle1 - (middle node of the ) element next in the
!                      natural order of elements
!!
!! @revision Nov 17
!-----------------------------------------------------------------------------
!
      subroutine nelconMG(Ngen,Mdle0, Mdle1)
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
        call nr_mdle_sons(NODES(nfath)%type,NODES(nfath)%ref_kind, nrbros)
        call locate(mdle,NODES(nfath)%sons,nrbros, noson)
        if (iprint.eq.1) then
          write(*,7002) mdle,nfath,nrbros,noson
 7002     format('nelconMG: mdle,nfath,nrbros,noson = ',2i7,2i3)
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
!
!  ...THIS LINE is the ONLY departure from the standard nelcon
      do while ((NODES(mdle)%ref_kind.gt.0).and.(Gen_lev(mdle).le.Ngen))
        mdle = NODES(mdle)%sons(1)
      enddo
      Mdle1 = mdle
!
!
      if (iprint.eq.1) then
        write(*,7010) Mdle0,Mdle1
 7010   format('nelconMG: Mdle0,Mdle1 = ',2i7)
      endif
!
      end subroutine nelconMG


