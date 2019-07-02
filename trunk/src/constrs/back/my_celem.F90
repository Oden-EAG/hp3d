      subroutine my_celem(Vec)

      implicit none
      real*8, dimension(:), pointer, intent(out) :: Vec

      nullify(Vec)

      end subroutine my_celem
