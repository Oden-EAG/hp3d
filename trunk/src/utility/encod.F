c----------------------------------------------------------------------
c
c   routine name       - encod
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Aug 07
c
c   purpose            - routine encodes a sequence of digits
c                        into a nickname from left to right
c
c   arguments :
c     in:
c               Narray - digits to be encoded
c               Mod    - system number
c               N      - number of the digits stored in the nickname
c     out:
c               Nick   - the nickname
c
c   required  routines -
c
c----------------------------------------------------------------------
c
      subroutine encod(Narray,Mod,N, Nick)
c
      dimension Narray(N)
c
      Nick=Narray(1)
      do i=2,N
        Nick = Nick*Mod + Narray(i)
      enddo
c
c
      end
c
c
c
c----------------------------------------------------------------------
c
c   routine name       - encodLong
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Jun 15
c
c   purpose            - routine encodes a sequence of digits
c                        into a LONG nickname from left to right
c
c   arguments :
c     in:
c               Narray - digits to be encoded
c               Mod    - system number
c               N      - number of the digits stored in the nickname
c     out:
c               Nick   - the LONG (integer) nickname
c
c   required  routines -
c
c----------------------------------------------------------------------
c
      subroutine encodLong(Narray,Mod,N, Nick)
c
      implicit none
c
      integer,             intent(in)    :: Mod,N
      integer,dimension(N),intent(in)    :: Narray
      integer(8),           intent(out)  :: Nick
c
c  ...local variables
      integer    :: i
      integer(8) :: longmod
c  ...be careful with integer type - need to convert!
      longmod = int(Mod,8)
      Nick = int(Narray(1),8)
      do i=2,N
c    ...be careful with integer type - need to convert!
        Nick = Nick*longmod + int(Narray(i),8)
      enddo
c
c
      end subroutine encodLong
