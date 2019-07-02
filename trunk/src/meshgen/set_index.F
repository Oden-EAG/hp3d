c---------------------------------------------------------------------
c   latest revision    - Jun 15
c
c   purpose            - set index for a node using the nodal case
c                        and boundary condition flags
c
c   arguments
c     in:
c          Icase       - node case
c          Iflag       - BC flag
c     out:
c          Index       - Vector indicating presence and kind of
c                        particular variables encoded decimally
c                        into a single (LONG) integer
c
c          Explanation of index in the expanded (or decimal) mode
c          indexd. For the i-th component:
c
c          indexd(i) = 0  component does not exist
c                    = 1  H1 component with Dirichlet BC flag
c                    = 2  free H1 component
c                    = 3  H(curl) component with Dirichlet BC flag
c                    = 4  free H(curl) component
c                    = 5  H(div) component with Dirichlet BC flag
c                    = 6  free H(div) component
c                    = 7  L2 component with Dirichlet BC flag
c                    = 8  free L2 component
c---------------------------------------------------------------------
c
      subroutine set_index(Icase,Iflag, Index)
c
      use physics
c
      implicit none
c
      integer,   intent(in)  :: Icase,Iflag
      integer*8, intent(out) :: Index
c
c  ...local variables
c  ...index in the decimal form
      integer,dimension(NRINDEX) :: indexd
c  ...binary version of Icase
      integer,dimension(NR_PHYSA) :: ncase
c  ...decimal version of the BC flag
      integer,dimension(NR_PHYSA) :: ibcd
c  ...others
      integer :: i,j,ic,iprint
c
c-----------------------------------------------------------------------
c
      iprint=0
c
c  ...decode the Icase flag, decode the BC flag
      call decod(Icase, 2,NR_PHYSA, ncase)
      call decod(Iflag,10,NR_PHYSA, ibcd )
c
c  ...initiate index counter
      ic=0
c
c  ...loop through the physics attributes
      do i=1,NR_PHYSA
c
c  .....if physical attribute is absent, skip its components
        if (ncase(i).eq.0) then
          do j=1,NR_COMP(i)
            ic=ic+1 ; indexd(ic)=0
          enddo
c
c  .....physical attribute is present
        else
          select case(DTYPE(i))
c
c  .......H1 variable
          case('contin')
c
c  .........loop through components
            do j=1,NR_COMP(i)
              ic=ic+1
c
c  ...........free H1 component
              indexd(ic)=2
c
c  ...........Dirichlet BC on ALL components
              if (ibcd(i).eq.1) indexd(ic)=1
c
c  ...........Dirichlet BC on 2nd and 3rd components
              if ((ibcd(i).eq.3).and.((j.eq.2).or.(j.eq.3))) then
                indexd(ic)=1
              endif
c
c  ...........Dirichlet BC on 1st and 3rd components
              if ((ibcd(i).eq.4).and.((j.eq.1).or.(j.eq.3))) then
                indexd(ic)=1
              endif
c
c  ...........Dirichlet BC on 1st and 2nd components
              if ((ibcd(i).eq.5).and.((j.eq.1).or.(j.eq.2))) then
                indexd(ic)=1
              endif
c
c  ...........Dirichlet BC on 1st component
              if ((ibcd(i).eq.6).and.(j.eq.1)) indexd(ic)=1
c
c  ...........Dirichlet BC on 2nd component
              if ((ibcd(i).eq.7).and.(j.eq.2)) indexd(ic)=1
c
c  ...........Dirichlet BC on 3rd component
              if ((ibcd(i).eq.8).and.(j.eq.3)) indexd(ic)=1
            enddo
c
c  .......H(curl) variable
          case('tangen')
c
c  .........loop through components
            do j=1,NR_COMP(i)
              ic=ic+1
c
c  ...........free H(curl) component
              indexd(ic)=4
c
c  ...........Dirichlet BC
              if (ibcd(i).eq.1) indexd(ic)=3
c
c  ...........Dirichlet BC on 2nd and 3rd components
              if ((ibcd(i).eq.3).and.((j.eq.2).or.(j.eq.3))) then
                indexd(ic)=3
              endif
c
c  ...........Dirichlet BC on 1st and 3rd components
              if ((ibcd(i).eq.4).and.((j.eq.1).or.(j.eq.3))) then
                indexd(ic)=3
              endif
c
c  ...........Dirichlet BC on 1st and 2nd components
              if ((ibcd(i).eq.5).and.((j.eq.1).or.(j.eq.2))) then
                indexd(ic)=3
              endif
c
c  ...........Dirichlet BC on 1st component
              if ((ibcd(i).eq.6).and.(j.eq.1)) indexd(ic)=3
c
c  ...........Dirichlet BC on 2nd component
              if ((ibcd(i).eq.7).and.(j.eq.2)) indexd(ic)=3
c
c  ...........Dirichlet BC on 3rd component
              if ((ibcd(i).eq.8).and.(j.eq.3)) indexd(ic)=3
c              
c          ...specific for impedance    
c  ...........Dirichlet BC on 2rd component
              if ((ibcd(i).eq.9).and.(j.eq.2)) indexd(ic)=3              
            enddo
c
c  .......H(div) variable
          case('normal')
c
c  .........loop through components
            do j=1,NR_COMP(i)
              ic=ic+1
c
c  ...........free H(div) component
              indexd(ic)=6
c
c  ...........Dirichlet BC
              if (ibcd(i).eq.1) then
                indexd(ic)=5
              endif
c  ...........impedance bc: 
              if (ibcd(i).eq.9) then
c             ...Eliminate H(div) dof (trick to avoid singular ZalocVV)
                indexd(ic)=5
              endif              
c
c  ...........Dirichlet BC on 1st component
              if ((ibcd(i).eq.3).and.(j.eq.1)) then
                indexd(ic)=5
              endif
c
c  ...........Dirichlet BC on 2nd component
              if ((ibcd(i).eq.4).and.(j.eq.2)) then
                indexd(ic)=5
              endif
c
c  ...........Dirichlet BC on 3rd component
              if ((ibcd(i).eq.5).and.(j.eq.3)) then
                indexd(ic)=5
              endif
c
c  ...........Dirichlet BC on 2nd and 3rd components
              if ((ibcd(i).eq.6).and.((j.eq.2).or.(j.eq.3))) then
                indexd(ic)=5
              endif
c
c  ...........Dirichlet BC on 1st and 3rd components
              if ((ibcd(i).eq.7).and.((j.eq.1).or.(j.eq.3))) then
                indexd(ic)=5
              endif
c
c  ...........Dirichlet BC on 1st and 2nd components
              if ((ibcd(i).eq.8).and.((j.eq.1).or.(j.eq.2))) then
                indexd(ic)=5
              endif
            enddo
c
c  .......L2 variable
          case('discon')
c
c  .........loop through components
            do j=1,NR_COMP(i)
              ic=ic+1
c
c  ...........free L2 component
              indexd(ic)=8
c
c  ...........Dirichlet BC
              if (ibcd(i).eq.1) indexd(ic)=7
            enddo
          end select
        endif
c
c  ...end of loop through physics attributes
      enddo
      if (ic.ne.NRINDEX) then
        write(*,*) 'set_index: INCONSISTENCY '
        stop
      endif
c
      call encodLong(indexd,10,NRINDEX, Index)
c
c
      if (iprint.eq.1) then
        write(*,7001) Icase,Iflag,indexd
 7001   format('set_index: Icase,Iflag,indexd = ',i3,i3,3x,10i1)
        write(*,7002) Index
 7002   format('           Index = ',i10)
        call pause
      endif
c
c
      end
