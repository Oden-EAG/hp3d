c----------------------------------------------------------------------
c
c   routine name       - shape3Q
c
c----------------------------------------------------------------------
c
c   latest revision    - Oct 2012
c
c   purpose            - routine evaluates L2 shape functions for 3D
c                        elements of various types
c
c   arguments :
c     in:
c         Type         - element type
c         Xi           - master element coordinates [0,1]x[0,1]x[0,1]
c         Norder       - order of approximation for the nodes (in the H1 sense)
c
c     out:
c         Nrdof        - number of the element shape functions      
c         ShapQ        - values of shape functions
c----------------------------------------------------------------------
c
      subroutine shape3Q_0(Type,Xi,Norder,Nrdof,ShapQ)
c      
      use parameters

      implicit none
c 
c  ...subroutine arguments
      character(len=4)                     , intent(in)  :: Type 
      double precision, dimension(3)       , intent(in)  :: Xi
      integer , dimension(19)              , intent(in)  :: Norder
      integer ,                              intent(out) :: Nrdof
      double precision, dimension(MAXBrickQ),intent(out) :: ShapQ
c
c  ...local variables
      integer :: iprint
c  
c----------------------------------------------------------------------
c      
      iprint=0
c      
      if (iprint.eq.1) then
        write(*,7001)Type,Xi(1:3)
 7001   format(' shape3Q: Type,Xi = ',a5,2x,3(e12.5,2x))
      endif
c      
      select case(Type)
      case('pris','mdlp')
        write(*,*) 'shape3Q: L2 functions for prism not implemented.'
        write(*,*) '         Switch to SHAPE=1.'
        stop
c        
      case('bric','mdlb')
        call shapeQb(Xi,Norder,Nrdof,ShapQ)
c        
      case('tetr','mdln')
        write(*,*) 'shape3Q: L2 functions for tet not implemented'
        write(*,*) '         Switch to SHAPE=1.'
        stop
c        
      case('pyra','mdld')
        write(*,*) 'shape3Q: L2 functions for pyramid not implemented'
        write(*,*) '         Switch to SHAPE=1.'
        stop
c        
      case default
        write(*,*) 'shape3Q: error, Type = ', Type
        stop
c        
      endselect
c
c
      endsubroutine shape3Q_0
