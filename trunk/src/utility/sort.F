c-----------------------------------------------------------------------
c
c   routine name       - sort
c
c-----------------------------------------------------------------------
c
c   latest revision    - jul,1991
c
c   purpose            - routine sorts out elements in the decreasing
c                        values in Val order      
c
c   usage              - call sort(List,Val,N)
c
c   arguments :
c     in:
c             List     - a list of elements
c             Val      - a list of element values
c             N        - number of elements on the list      
c     out:
c             List     - sorted out elements
c             Val      - sorted out values      
c
c-----------------------------------------------------------------------
c
      subroutine sort(List,Val,N)
c      
      implicit none
c
      integer              , intent(in   ) :: N
      integer, dimension(N), intent(inout) :: List
      real(8), dimension(N), intent(inout) :: Val
c      
      real(8) :: valmax,valaux
      integer :: iprint,i,imax,j,laux
c      
c-----------------------------------------------------------------------
c
      iprint=0
c      
c  ...printing
      if (iprint.eq.1) then
        write(*,*) 'ORIGINAL LISTS, N=', N
        do i=1,N
          write(*,7001)i,List(i),Val(i)
        enddo
      endif
c      
c  ...DO NOT USE UNINITIALIZED VARIABLES
c     ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
c     
c  ...loop through elements
      do i=1,N-1
        valmax= Val(i)
        imax  = i
c
c  .....loop through remaining elements          
        do j=i,N
          if (Val(j).gt.valmax) then
            valmax=Val(j) ; imax=j
          endif
        enddo
c
c  .....switch places between the i-th and the maximum value elements
        valaux   =Val(i)    ; laux      =List(i)
        Val(i)   =Val(imax) ; List(i)   =List(imax)
        Val(imax)=valaux    ; List(imax)=laux
c        
c  ...end loop through elements        
      enddo
c
c  ...printing      
      if (iprint.eq.1) then
        write(*,*) 'SORTED OUT, N=', N
        do i=1,N
          write(*,7001) i,List(i),Val(i)
 7001     format(1x,'i,List,Val = ',i4,i6,e12.5)
        enddo
      endif
c     
c
      endsubroutine sort
