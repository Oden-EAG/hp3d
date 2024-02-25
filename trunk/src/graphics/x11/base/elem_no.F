#if HP3D_USE_X11

c----------------------------------------------------------------------
c
c   routine name       - elem_no
c
c----------------------------------------------------------------------
c
c   latest revision    - Feb 2023
c
c   purpose            - given a point on the projected plane,
c                        routine identifies an element's face
c                        that contains the point and that is
c                        closest to the point of view
c
c   arguments :
c     in:
c                Xy_in - coordinates of a projected point
c     out:
c                Mdle  = the (middle node of) element, if the search
c                        has been succesful
c                      = 0 otherwise
c
c----------------------------------------------------------------------
c
      subroutine elem_no(Xy_in, Mdle)
c
      use data_structure3D
      use graphmod
      use node_types
c
      implicit none
c
      real(8), intent(in)  :: Xy_in(2)
      integer, intent(out) :: Mdle
c
c  ...master face coordinates
      real(8) :: t(2)
c
c  ...physical coordinates
      real(8) :: xyz(3)
c
c  ...neighbors of an element
      integer :: neig(4,6)
c
c  ...list of elements intersecting with the projecting line and
c     physical coordinates of the corresponding intersection points
      integer, parameter :: maxlist = 10
      integer :: list_elem(maxlist)
      real(8) :: xyz_list(3,maxlist)
c
      real(8) :: d,d_min,fact
      integer :: iel,i_min,nr_elem,loc,ivis,is,ifc,locn,nfl
      integer :: i,j,ifinside,ndom,ndomn,mdlen
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
c
c  ...initiate number of elements found
      nr_elem=0
c
c  ...loop through elements
      Mdle = 0
      do iel=1,NRELES
c
        call nelcon(Mdle, Mdle)
c
c  .....check if on the list of invisible elements
        call locate(Mdle,IGINV,NRINVBL, loc)
        if (loc.gt.0) cycle
        call find_domain(Mdle, ndom)
        if (NDOMAIN(ndom).eq.0) cycle
c
c  .....get neighbors
        call find_neig(Mdle, neig)
c
c  .....loop through element sides
        do ifc=1,nface(NODES(Mdle)%ntype)
c
c  .......check visibility of the face...
          ivis=0
c
c  .......loop through the face neighbors
          do is=1,4
            mdlen = neig(is,ifc)
            if (mdlen.eq.0) then
              ivis=1
            else
              call locate(mdlen,IGINV,NRINVBL, locn)
              if (locn.gt.0) ivis=1
              call find_domain(mdlen, ndomn)
              if (NDOMAIN(ndomn).eq.0) ivis=1
            endif
          enddo
          if (ivis.eq.0) cycle
          call invmap_face(Mdle,ifc,Xy_in, t,xyz,nfl)
#if HP3D_DEBUG
          if (iprint.eq.1) then
            write(*,7002) Mdle,ifc,t,xyz,nfl
 7002       format('elem_no: Mdle,ifc,t,xyz,nfl = ',
     .                       i4,i2,2x,2f6.2,2x,3f6.2,i2)
          endif
#endif
c
c  .......if the inverse map has converged...
          if (nfl.eq.0) then
c
c  .........check if within the face
            ifinside=0
            select case(face_type(NODES(Mdle)%ntype,ifc))
            case(TRIA)
              if ((t(1).gt.0.d0).and.(t(2).gt.0.d0).and.
     .            (t(1)+t(2).lt.1.d0)) ifinside=1
            case(RECT)
              if ((t(1).gt.0.d0).and.(t(2).gt.0.d0).and.
     .            (t(1).lt.1.d0).and.(t(2).lt.1.d0)) ifinside=1
            end select
            if (ifinside.eq.1) then
c
c  ...........store on the list
              nr_elem = nr_elem+1
              if (nr_elem.gt.maxlist) then
                write(*,*) 'elem_no: INCREASE maxlist'
                stop 1
              endif
              list_elem(nr_elem) = Mdle
              xyz_list(1:3,nr_elem) = xyz(1:3)
            endif
          endif
c
c  .....end of loop through the element faces
        enddo
c
c  ...end of loop through elements
      enddo
c
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,*) 'elem_no: list_elem,xyz_list  = '
        do i=1,nr_elem
          write(*,7004) list_elem(i),xyz_list(1:3,i)
 7004     format(i4,2x,3f8.3)
        enddo
        call pause
      endif
#endif
c
      if (nr_elem.eq.0) then
c
c  .....no element has been found, return with zero flag
        Mdle=0
      else
c
c  .....select the element that is closest to the point of view
        fact = 10.d0
        i_min=0
        d_min = 1.d30
        do i=1,nr_elem
          d = 0.d0
          do j=1,3
            d = d + (fact*RN(j) - xyz_list(j,i))**2
          enddo
          d = sqrt(d)
          if (d.lt.d_min) then
            d_min = d
            i_min = i
          endif
        enddo
        Mdle = list_elem(i_min)
      endif
#if HP3D_DEBUG
      if (iprint.eq.1) then
        write(*,*)'elem_no: FINAL Mdle = ',Mdle
        call pause
      endif
#endif
c
c
      end subroutine elem_no

#endif
