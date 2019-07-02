!
!------------------------------------------------------------------------
!
!    function name      - Facetype
!
!------------------------------------------------------------------------
!
!    latest revision   - Jan 2018
!
!    purpose           - returns the face type
!                        (temporary function)
!
!   arguments :
!     in:
!         Type         - element type
!           if         - face number
!    out:
!         Facetype     - face type
!
!-----------------------------------------------------------------------
!
   function Facetype(Type,If)
!
   character(len=4) :: Facetype,Type
!
   select case(Type)
   case('mdlp')
      select case(If)
      case(1,2); 
         Facetype = 'mdlt'
      case(3,4,5)
         Facetype = 'mdlq'
      end select
   case('mdlb')
      Facetype = 'mdlq'
   case('mdln')
      Facetype = 'mdlt'
   case('mdld')
      select case(If)
      case(1)
         Facetype = 'mdlq'
      case(2,3,4,5)
         Facetype = 'mdlt' 
      end select
   end select
!   
   end function Facetype