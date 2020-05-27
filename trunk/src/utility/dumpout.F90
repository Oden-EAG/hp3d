subroutine dumpout
!      
      use GMP              , only : dumpout_GMP
      use data_structure3D , only : dumpout_physics,dumpout_hp3d
!      
      implicit none
      write(*,*)'Dumping out GMP...'     ; call dumpout_GMP
      write(*,*)'Done!'
      write(*,*)'Dumping out physics...' ; call dumpout_physics
      write(*,*)'Done!'
      write(*,*)'Dumping out HP3D...'    ; call dumpout_hp3d('./files/dumpc3D')
      write(*,*)'Done!'
!
!
end subroutine dumpout
