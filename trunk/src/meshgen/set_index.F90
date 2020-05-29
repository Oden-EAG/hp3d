!---------------------------------------------------------------------
!   latest revision    - May 2020
!
!   purpose            - set index for a node using the nodal case
!                        and boundary condition flags
!
!   arguments
!     in:
!          Icase       - node case
!          Iflag       - BC flag
!     out:
!          Index       - Vector indicating presence and kind of
!                        particular variables encoded decimally
!                        into a single (LONG) integer
!
!          Explanation of index in the expanded (or decimal) mode
!          indexd. For the i-th component:
!
!          indexd(i) = 0  component does not exist
!                    = 1  H1 component with Dirichlet BC flag
!                    = 2  free H1 component
!                    = 3  H(curl) component with Dirichlet BC flag
!                    = 4  free H(curl) component
!                    = 5  H(div) component with Dirichlet BC flag
!                    = 6  free H(div) component
!                    = 7  L2 component with Dirichlet BC flag
!                    = 8  free L2 component
!---------------------------------------------------------------------
!
subroutine set_index(Icase,Iflag, Index)
!
      use physics
!
      implicit none

      integer,    intent(in)  :: Icase,Iflag
      integer(8), intent(out) :: Index
!
!  ...local variables
!  ...index in the decimal form
      integer,dimension(NRINDEX)  :: indexd
!  ...binary version of Icase
      integer,dimension(NR_PHYSA) :: ncase
!  ...decimal version of the BC flag
      integer,dimension(NR_PHYSA) :: ibcd
!  ...others
      integer :: i,j,ic
!
#if DEBUG_MODE
      integer :: iprint = 0
#endif
!
!-----------------------------------------------------------------------
!
!  ...decode the Icase flag, decode the BC flag
      call decod(Icase, 2,NR_PHYSA, ncase)
      call decod(Iflag,10,NR_PHYSA, ibcd )
!
!  ...initiate index counter
      ic=0
!
!  ...loop through the physics attributes
      do i=1,NR_PHYSA
!
!  .....if physical attribute is absent, skip its components
        if (ncase(i).eq.0) then
          do j=1,NR_COMP(i)
            ic=ic+1 ; indexd(ic)=0
          enddo
!
!  .....physical attribute is present
        else
          select case(DTYPE(i))
!
!  .......H1 variable
          case('contin')
!
!  .........loop through components
            do j=1,NR_COMP(i)
              ic=ic+1
!
!  ...........free H1 component
              indexd(ic)=2
!
!  ...........Dirichlet BC on ALL components
              if (ibcd(i).eq.1) indexd(ic)=1
!
!  ...........Dirichlet BC on 2nd and 3rd components
              if ((ibcd(i).eq.3).and.((j.eq.2).or.(j.eq.3))) then
                indexd(ic)=1
              endif
!
!  ...........Dirichlet BC on 1st and 3rd components
              if ((ibcd(i).eq.4).and.((j.eq.1).or.(j.eq.3))) then
                indexd(ic)=1
              endif
!
!  ...........Dirichlet BC on 1st and 2nd components
              if ((ibcd(i).eq.5).and.((j.eq.1).or.(j.eq.2))) then
                indexd(ic)=1
              endif
!
!  ...........Dirichlet BC on 1st component
              if ((ibcd(i).eq.6).and.(j.eq.1)) indexd(ic)=1
!
!  ...........Dirichlet BC on 2nd component
              if ((ibcd(i).eq.7).and.(j.eq.2)) indexd(ic)=1
!
!  ...........Dirichlet BC on 3rd component
              if ((ibcd(i).eq.8).and.(j.eq.3)) indexd(ic)=1
            enddo
!
!  .......H(curl) variable
          case('tangen')
!
!  .........loop through components
            do j=1,NR_COMP(i)
              ic=ic+1
!
!  ...........free H(curl) component
              indexd(ic)=4
!
!  ...........Dirichlet BC
              if (ibcd(i).eq.1) indexd(ic)=3
!
!  ...........Dirichlet BC on 2nd and 3rd components
              if ((ibcd(i).eq.3).and.((j.eq.2).or.(j.eq.3))) then
                indexd(ic)=3
              endif
!
!  ...........Dirichlet BC on 1st and 3rd components
              if ((ibcd(i).eq.4).and.((j.eq.1).or.(j.eq.3))) then
                indexd(ic)=3
              endif
!
!  ...........Dirichlet BC on 1st and 2nd components
              if ((ibcd(i).eq.5).and.((j.eq.1).or.(j.eq.2))) then
                indexd(ic)=3
              endif
!
!  ...........Dirichlet BC on 1st component
              if ((ibcd(i).eq.6).and.(j.eq.1)) indexd(ic)=3
!
!  ...........Dirichlet BC on 2nd component
              if ((ibcd(i).eq.7).and.(j.eq.2)) indexd(ic)=3
!
!  ...........Dirichlet BC on 3rd component
              if ((ibcd(i).eq.8).and.(j.eq.3)) indexd(ic)=3
!
!          ...specific for EM impedance BC
!  ...........Impedance BC for Primal Maxwell
              !if ((ibcd(i).eq.8).and.(j.eq.1)) indexd(ic)=3
!  ...........Impedance BC for UW Maxwell
              !if ((ibcd(i).eq.9).and.(j.eq.2)) indexd(ic)=3
            enddo
!
!  .......H(div) variable
          case('normal')
!
!  .........loop through components
            do j=1,NR_COMP(i)
              ic=ic+1
!
!  ...........free H(div) component
              indexd(ic)=6
!
!  ...........Dirichlet BC
              if (ibcd(i).eq.1) then
                indexd(ic)=5
              endif
!
!  ...........Dirichlet BC on 1st component
              if ((ibcd(i).eq.3).and.(j.eq.1)) then
                indexd(ic)=5
              endif
!
!  ...........Dirichlet BC on 2nd component
              if ((ibcd(i).eq.4).and.(j.eq.2)) then
                indexd(ic)=5
              endif
!
!  ...........Dirichlet BC on 3rd component
              if ((ibcd(i).eq.5).and.(j.eq.3)) then
                indexd(ic)=5
              endif
!
!  ...........Dirichlet BC on 2nd and 3rd components
              if ((ibcd(i).eq.6).and.((j.eq.2).or.(j.eq.3))) then
                indexd(ic)=5
              endif
!
!  ...........Dirichlet BC on 1st and 3rd components
              if ((ibcd(i).eq.7).and.((j.eq.1).or.(j.eq.3))) then
                indexd(ic)=5
              endif
!
!  ...........Dirichlet BC on 1st and 2nd components
              if ((ibcd(i).eq.8).and.((j.eq.1).or.(j.eq.2))) then
                indexd(ic)=5
              endif
            enddo
!
!  ...........specific for acoustic impedance BC:
              !if (ibcd(i).eq.9) then
!             ...Eliminate H(div) dof (trick to avoid singular ZalocVV)
              !  indexd(ic)=5
              !endif
!
!  .......L2 variable
          case('discon')
!
!  .........loop through components
            do j=1,NR_COMP(i)
              ic=ic+1
!
!  ...........free L2 component
              indexd(ic)=8
!
!  ...........Dirichlet BC
              if (ibcd(i).eq.1) indexd(ic)=7
            enddo
          end select
        endif
!
!  ...end of loop through physics attributes
      enddo
      if (ic.ne.NRINDEX) then
        write(*,*) 'set_index: INCONSISTENCY.'
        stop
      endif
!
      call encodLong(indexd,10,NRINDEX, Index)
!
!
#if DEBUG_MODE
      if (iprint.eq.1) then
        write(*,7001) Icase,Iflag,indexd
 7001   format('set_index: Icase,Iflag,indexd = ',i3,i3,3x,10i1)
        write(*,7002) Index
 7002   format('           Index = ',i10)
        call pause
      endif
#endif
!
!
end subroutine
