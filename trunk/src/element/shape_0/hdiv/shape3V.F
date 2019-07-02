c----------------------------------------------------------------------
c
c   routine name       - shape3V
c
c----------------------------------------------------------------------
c
c   latest revision    - Nov 13
c
c   purpose            - routine evaluates H(div) shape functions for 
c                        3D elements of various types
c
c   arguments :
c     in:
c         Type         - element type
c         Xi           - master element coordinates
c         Norder       - order of approximation for the nodes
c         Nface_orient - face orientations
c     out:
c         Nrdof        - number of the element shape functions      
c         ShapV        - shape functions
c         DivV         - div of shape functions
c
c----------------------------------------------------------------------
c
      subroutine shape3V_0(Type,Xi,Norder,Nface_orient, 
     .                     Nrdof,ShapV,DivV)
c      
      use parameters
#include "syscom.blk"
c 
      character (len=*),intent(in) :: Type 
      dimension Nface_orient(6),Norder(19)
      dimension Xi(3),ShapV(3,MAXbrickV),DivV(MAXbrickV)
c 
c----------------------------------------------------------------------
c      
      write(*,*)'OLD CONSTRUCTION OF SHAPE FUNCTIONS DID'
      write(*,*)'NOT INCLUDE H(div)! PLEASE SWITCH TO NEW'
      write(*,*)'CONSTRUCTION.'
      stop
c
c
      endsubroutine shape3V_0
