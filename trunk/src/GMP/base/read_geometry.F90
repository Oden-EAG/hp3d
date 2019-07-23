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
      case(DUMPIN_LEGACY_) ; call dumpin_GMP_LEGACY('files/dumpGMP')      
         write(*,*)'dumping in GMP (LEGACY)...'
         write(*,*)''
      case(DUMPIN_)        ; call dumpin_GMP       ('files/dumpGMP')      
         write(*,*)'dumping in GMP...'
         write(*,*)''
!  ...LEGACY: DO NOT USE!
      case(LEGACY_)        ; call input_geometry   (Fp)
      case(DEFAULT_)       ; call input_DEFAULT    (Fp)
!  ...LEGACY: DO NOT USE!      
      case(COMPACT_)       ; call input_GMPdata    (Fp)
      case(RECONSTRUCT_)   ; call input_reconstruct(Fp)
      case(NETGEN_)        ; call input_NETGEN     (Fp)
      end select    
!
!      
endsubroutine read_geometry
