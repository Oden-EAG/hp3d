! ----------------------------------------------------------------------
! 
!    routine name       - nelcon_macro
! 
! ----------------------------------------------------------------------
! 
!    latest revision    - Jan 2018
! 
!    purpose            - routine returns the next element in the fine
!                         mesh corresponding to a coarse element
!                         in the natural order of elements
! 
! 
!    arguments :
!      in:
!               Mdlec   - middle node of a coarse mesh element
!               Mdle0   - middle node of an element identified with
!                         an element
!      out:
!               Mdle1   - middle node of the element next in the
!                         natural ordering of elements
! 
! ----------------------------------------------------------------------
! 
   subroutine nelcon_macro(Mdlec,Mdle0, Mdle1)
! 
   use data_structure3D, ONLY: NODES
   use refinements,      ONLY: nr_mdle_sons
!
   IMPLICIT NONE
!
! ----------------------------------------------------------------------
!
   integer    :: iprint, mdle, mdle0, mdlec, mdle1, nfath, nrbros, noson

! ----------------------------------------------------------------------
   iprint=0
! 
   mdle = Mdle0
! 
! ..first element in the mesh
   if (mdle.eq.0) then
      mdle = Mdlec
      go to 20
   endif
! 
! ...Step 1: move horizontally, if you can, otherwise move vertically up
 10 continue
! 
! ...quit if you have come back to MdleC
   if (mdle.eq.Mdlec) then
      Mdle1=0
      return
   endif
! 
   nfath = NODES(mdle)%father
! 
! ...initial mesh element
   if (nfath.lt.0) then
! ...move to the next initial mesh element
      mdle = mdle + 1
      go to 20
   else
! 
! ...find the son number in the family
      call nr_mdle_sons(NODES(nfath)%type,NODES(nfath)%ref_kind,nrbros)
      call locate(mdle,NODES(nfath)%sons,nrbros, noson)
      if (iprint.eq.1) then
         write(*,7002) mdle,nfath,nrbros,noson
 7002 format('nelcon: mdle,nfath,nrbros,noson = ',2i7,2i3)
      endif
!                    
! ...if mdle is not the last son in the family, go to the next brother
      if (noson.lt.nrbros) then
         mdle = NODES(nfath)%sons(noson+1)
         go to 20
      else
!                 
! ...move up on the tree
         mdle = nfath
         go to 10
      endif
   endif
!  
! ...Step 2 move vertically down                
   20 continue
   do while (NODES(mdle)%ref_kind.gt.0)
      mdle = NODES(mdle)%sons(1)
   enddo
   Mdle1 = mdle
!
   if (iprint.eq.1) then
      write(*,7010) Mdle0,Mdle1
 7010 format('nelcon: Mdle0,Mdle1 = ',2i7)
   endif
! 
   end subroutine nelcon_macro
