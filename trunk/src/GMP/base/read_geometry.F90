subroutine read_geometry(Fp)
!
   use GMP
   use control
!
   implicit none
!
   character(len=*),intent(in) :: Fp
!
   select case(INPUT_FILE)
      case(DUMPIN_)        ; call dumpin_GMP       ('files/dumpGMP')
         write(*,*)'dumping in GMP...'
         write(*,*)''
      case(DEFAULT_)       ; call input_DEFAULT    (Fp)
!
!  ...LEGACY: DO NOT USE!
!     case(LEGACY_)        ; call input_geometry   (Fp)
!     case(COMPACT_)       ; call input_GMPdata    (Fp)
!     case(RECONSTRUCT_)   ; call input_reconstruct(Fp)
!     case(NETGEN_)        ; call input_NETGEN     (Fp)
      case default
         write(*,*) 'read_geometry: INPUT_FILE'; stop
   end select
!
!
end subroutine read_geometry
