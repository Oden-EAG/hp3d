!---------------------------------------------------------------------------------
!> @brief return isotropic refinement flag
!!
!! @param[in ] Nod  - a node number
!! @param[out] Kref - isotropic refinement flag
!!
!> @date Feb 2023
!---------------------------------------------------------------------------------
subroutine get_isoref(Nod, Kref)
!
      use data_structure3D
!
      implicit none
      integer, intent(in)     :: Nod
      integer, intent(out)    :: Kref
!
      real(8), dimension(3,8) :: xsub
      real(8), dimension(3)   :: dist, xi
      real(8), dimension(3,2) :: x
!
      integer, dimension(27)  :: nodesl,norientl
      integer, dimension(2)   :: iv
      integer, dimension(2,3), parameter :: ie = &
            reshape( (/1,6, 3,5, 4,2/), (/2,3/) )
      integer :: iflag, no, i,j,k, loc
!
#if HP3D_DEBUG
      integer :: iprint
      iprint=0
#endif
!---------------------------------------------------------------------------------
!
      Kref = 0
!
!  ...select refinement based on node type
      select case(NODES(Nod)%ntype)
!
!     EDGE
      case(MEDG) ; Kref=1
!
!     TRIANGLE
      case(MDLT) ; Kref=1
!
!     QUAD
      case(MDLQ) ; Kref=11
!
!     PRISM
      case(MDLP) ; Kref=11
!
!     BRICK
      case(MDLB) ; Kref=111
!
!     TET
      case(MDLN)
        call refel     (Nod, iflag,no,xsub)
        call elem_nodes(Nod, nodesl,norientl)
!
!  .....measure 3 diagonals
        dist(1:3)=0.d0
        do i=1,3
          do j=1,2
            do k=1,2
              iv(k) = (TETRA_EDGE_TO_VERT(k, ie(j,i)))
            enddo
!  .........compute the middle point between two end points
            xi = 0.d0
            do k=1,2
              xi = xi + xsub(1:3, iv(k))/2.d0
            enddo
!  .........compute physical coordinate for xi
            call hpvert(iflag,no,xi, x(1:3,j))
          enddo
!
!  .......distance
          do j=1,3
            dist(i) = dist(i) + (x(j,1) - x(j,2))**2
          enddo
        enddo
!
!  .....select longest diagonal
        loc=1
        do i=1,3
          if (dist(loc) > dist(i))  loc=i
        enddo
!
!  .....selec appropriate refinement kind
        Kref = 10 + loc
!
#if HP3D_DEBUG
        if (iprint.eq.1) then
          write(*,7000) Kref, dist
7000      format(' get_isoref : Kref = ',i3,' dist = ',3f8.3)
        endif
#endif
!
      case(MDLD)
        write(*,*) 'get_isoref: no isotropic refinement for pyramid.'
        stop
      endselect
!
end subroutine get_isoref
