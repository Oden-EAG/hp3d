c----------------------------------------------------------------------
c
c   routine name       - findno
c
c   latest revision    - Mar 2023
c
c   purpose            - given a point on the projected plane,
c                        routine identifies a GMP block
c                        that contains the point and that is
c                        closest to the point of view
c
c   arguments :
c     in:
c              Xy_in   - coordinates of a projected point
c     out:
c              Number  = Number of the figure, if the search
c                        was successful
c                      = 0 otherwise
c----------------------------------------------------------------------
c
      subroutine findno(Xy_in, Number)
c
      use GMP
      use graphmod
c
      implicit none
c
      real(8) :: Xy_in(2)
      integer :: Number
c
c  ...reference and physical coordinates
      real(8) :: eta(2),x(3)
c
c  ...list of figures
      integer, parameter :: maxfig=10
      integer :: list_figs(maxfig)
      real(8) :: xyz_list(1:3,maxfig), xyz_min(1:3)
c
c  ...blocks adjacent to a figure
      integer :: nbla(2)
c
      real(8) :: d,d_min,fact
      integer :: i,i_min,if,ifg,ii,is,ivis,j,lab
      integer :: nbl,nfl,norient,nr,nr_fig,nh,nick,nvoid
c
      integer :: iprint=0
c
      nr_fig=0
c
c.....loop through hexahedra
      do nh = 1,NRHEXAS
        nick=10*nh+2
c
c  .....check visibility of the hexahedron
        call locate (nick,IGINV,NRINVBL, ii)
        if (ii.ne.0) cycle
c
c  .....loop through the hexahedron faces
        do if = 1,6
          call decode(abs(HEXAS(nh)%FigNo(if)), nr,norient)
c
c  .......adjacent blocks
          nbla(1:2)=abs(RECTANGLES(nr)%BlockNo(1:2))
          ivis=0
          do is=1,2
c
c  .........boundary face
            if (nbla(is).eq.0) then
              ivis=1
            elseif (nbla(is).ne.nick) then
              call locate (nbla(is),IGINV,NRINVBL, ii)
              if (ii.ne.0) ivis=1
            endif
          enddo
          if (ivis.eq.0) cycle
c
c  .......determine the face coordinates
          call findcoord(nr,2,Xy_in, eta,x,nfl)
c
c  .......if the inverse map has converged...
          if (nfl.eq.1) then
            if ( (eta(1).ge.0.d0).and.(eta(1).le.1.d0).and.
     .           (eta(2).ge.0.d0).and.(eta(2).le.1.d0) ) then
c
c     .........store on the list
               nr_fig = nr_fig+1
               if (nr_fig.gt.maxfig) then
                 write(*,*) 'findno: maxfig = ',maxfig
                 stop 1
               endif
               list_figs(nr_fig) = nr*10+2
               xyz_list(1:3,nr_fig) = x(1:3)
             endif
           endif
c
c  ......end of loop through faces
         enddo
c
c  ....end of loop through hexas
       enddo
c
c
      if (nr_fig.eq.0) then
c
c  .....no has been found, return with zero flag
         Number=0
      else
c
c  .....select the point that is closest to the point of view
        fact = 10.d0
        i_min=0
        d_min = 1.d30
        do i=1,nr_fig
          d = 0.d0
          do j=1,3
            d = d + (fact*RN(j) - xyz_list(j,i))**2
          enddo
          d = sqrt(d)
          if (d.lt.d_min) then
            d_min = d
            i_min = list_figs(i)
            xyz_min(1:3) = xyz_list(1:3,i)
          endif
        enddo
        call decode(i_min,ifg,lab)
        select case(lab)
        case(1)
          write(*,*) 'findno: UNFINISHED'
        case(2)
          do is=1,2
            nbl = RECTANGLES(ifg)%BlockNo(is)
            if (nbl.ne.0) then
              call locate (nbl,IGINV,NRINVBL, ii)
c
c  ...........visible block
              if (ii.eq.0) then
                call decode(nbl, Number, nvoid)
              endif
            endif
          enddo
        end select
      endif
c
      write(*,*)'findno: NUMBER = ',NUMBER
c
      end subroutine findno
