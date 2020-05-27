subroutine get_ddet(A, Adet)
  implicit none
  real(8), dimension(3,3), intent(in)  :: A
  real(8),                 intent(out) :: Adet
  Adet = A(1,1)*A(2,2)*A(3,3) &
       + A(2,1)*A(3,2)*A(1,3) &
       + A(3,1)*A(1,2)*A(2,3) &
       - A(3,1)*A(2,2)*A(1,3) &
       - A(1,1)*A(3,2)*A(2,3) &
       - A(2,1)*A(1,2)*A(3,3)
end subroutine get_ddet

subroutine get_zdet(A, Adet)
  implicit none
  complex(8), dimension(3,3), intent(in)  :: A
  complex(8),                 intent(out) :: Adet
  Adet = A(1,1)*A(2,2)*A(3,3) &
       + A(2,1)*A(3,2)*A(1,3) &
       + A(3,1)*A(1,2)*A(2,3) &
       - A(3,1)*A(2,2)*A(1,3) &
       - A(1,1)*A(3,2)*A(2,3) &
       - A(2,1)*A(1,2)*A(3,3)
end subroutine get_zdet

subroutine dinvert(A, Ainv, Adet)
  implicit none
  real(8), dimension(3,3), intent(in)  :: A
  real(8), dimension(3,3), intent(out) :: Ainv
  real(8),                 intent(out) :: Adet
  real(8) :: det(3)

  call get_ddet(A, Adet)

  det(1) =   A(2,2)*A(3,3) - A(3,2)*A(2,3)
  det(2) = - A(2,1)*A(3,3) + A(3,1)*A(2,3)
  det(3) =   A(2,1)*A(3,2) - A(3,1)*A(2,2)

  Ainv(1,1) = det(1)/Adet
  Ainv(2,1) = det(2)/Adet
  Ainv(3,1) = det(3)/Adet

  det(1) =   A(3,2)*A(1,3) - A(1,2)*A(3,3)
  det(2) =   A(1,1)*A(3,3) - A(3,1)*A(1,3)
  det(3) = - A(1,1)*A(3,2) + A(3,1)*A(1,2)

  Ainv(1,2) = det(1)/Adet
  Ainv(2,2) = det(2)/Adet
  Ainv(3,2) = det(3)/Adet

  det(1) =   A(1,2)*A(2,3) - A(2,2)*A(1,3)
  det(2) = - A(1,1)*A(2,3) + A(2,1)*A(1,3)
  det(3) =   A(1,1)*A(2,2) - A(2,1)*A(1,2)

  Ainv(1,3) = det(1)/Adet
  Ainv(2,3) = det(2)/Adet
  Ainv(3,3) = det(3)/Adet

end subroutine dinvert


subroutine zinvert(A, Ainv, Adet)
  implicit none
  complex(8), dimension(3,3), intent(in)  :: A
  complex(8), dimension(3,3), intent(out) :: Ainv
  complex(8),                 intent(out) :: Adet
  complex(8) :: det(3)

  call get_zdet(A, Adet)

  det(1) =   A(2,2)*A(3,3) - A(3,2)*A(2,3)
  det(2) = - A(2,1)*A(3,3) + A(3,1)*A(2,3)
  det(3) =   A(2,1)*A(3,2) - A(3,1)*A(2,2)

  Ainv(1,1) = det(1)/Adet
  Ainv(2,1) = det(2)/Adet
  Ainv(3,1) = det(3)/Adet

  det(1) =   A(3,2)*A(1,3) - A(1,2)*A(3,3)
  det(2) =   A(1,1)*A(3,3) - A(3,1)*A(1,3)
  det(3) = - A(1,1)*A(3,2) + A(3,1)*A(1,2)

  Ainv(1,2) = det(1)/Adet
  Ainv(2,2) = det(2)/Adet
  Ainv(3,2) = det(3)/Adet

  det(1) =   A(1,2)*A(2,3) - A(2,2)*A(1,3)
  det(2) = - A(1,1)*A(2,3) + A(2,1)*A(1,3)
  det(3) =   A(1,1)*A(2,2) - A(2,1)*A(1,2)

  Ainv(1,3) = det(1)/Adet
  Ainv(2,3) = det(2)/Adet
  Ainv(3,3) = det(3)/Adet

end subroutine zinvert


