!> @brief Saves common variables used in frontal solver
!> @date Mar 2024
module czdirio
!
   implicit none
!
   save
!
#if HP3D_COMPLEX
   integer    :: nbuf(9), lenf(9), irsave(9)
   complex(8) :: storage(100)
#endif
!
end module czdirio
