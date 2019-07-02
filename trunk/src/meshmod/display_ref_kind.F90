subroutine display_ref_kind
!
  write(*,*)'================================================================='
  write(*,*)'NODE TYPE      // ISOTROPIC REFINEMENT // ANISOTROPIC REFINEMENT '
  write(*,*)'mdlp (prism)      11                      10,01                  '
  write(*,*)'mdlb (hexa)       111                     110,101,011,100,010,001'
  write(*,*)'mdln (tetra)      11,12,13                24,32                  '
  write(*,*)'mdld (pyramid)                            10                     '
  write(*,*)'================================================================='
!
end subroutine 
