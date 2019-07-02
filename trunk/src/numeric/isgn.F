c----------------------------------------------------------------------
c
c   routine name       - sgn
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Nov 09
c
c   purpose            - function return the sign of expression
c
c   arguments :
c     in:
c     out:
c----------------------------------------------------------------------
c
      integer function isgn(I)
c
      integer I
c
      if (I.eq.0) isgn=  0
      if (I.gt.0) isgn=  1
      if (I.lt.0) isgn= -1
c
      return
      end


     
        
      

