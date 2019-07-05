!----------------------------------------------------------------------
!
!   routine name       - find_case
!
!----------------------------------------------------------------------
!
!   computer           - machine independent
!
!   latest revision    - Dec 07
!
!   purpose            - routine identifies the node case number
!                        based on physics info
!
!   arguments :
!     in:
!          Number      - number of physics attributes for a specific
!                        node
!          Phys        - the corresponding list of the attributes
!     out:
!          Icase       - node case number
!
!   required  routines - none
!
!----------------------------------------------------------------------
!
subroutine find_case(Number,Phys, Icase)
!
      use physics
!
      implicit none
!
      integer         , intent(in)  :: Number
      character(len=5), intent(in)  :: Phys(Number)
      integer         , intent(out) :: Icase
!
      integer :: narray(10)
      integer :: iprint, ii, j
!
!----------------------------------------------------------------------
!
      iprint=0
!
      narray(1:10)=0
!
      if (iprint.eq.1) then
        write(*,7001) Number, Phys(1:Number)
 7001   format(' find_case: Phys(1:',i1,') = ',5(a5,' ; '))
      endif
!
!  ...loop through physics attributes
      do j=1,NR_PHYSA
        call locate_char(PHYSA(j),Phys,Number, ii)
        if (ii.ne.0) narray(j)=1
      enddo
!
!  ...encode the case number
      call encod(narray,2,NR_PHYSA, Icase)
      if (iprint.eq.1) then
        write(*,7002) Icase
 7002   format(' find_case: Icase = ',i2)
        call pause
      endif
!
!
end subroutine find_case
