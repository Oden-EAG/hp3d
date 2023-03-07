c Routines:
c  - checkorder
c----------------------------------------------------------------------
c
c     routine name      - checkorder
c
c----------------------------------------------------------------------
c
c     latest revision:  - Feb 2023
c
c     purpose:          - routine checks whether polynomial orders of
c                         edges, faces and bubbles are within
c                         appropriate bounds. Returns values that help
c                         to more efficiently size temporary arrays in
c                         other routines in order to save memory.
c
c     arguments:
c
c     in:
c       Ntype           - element type
c       Dtype           - discretization type (H1,H(curl),H(div),L2)
c       Norder          - polynomial order for the nodes (H1 sense)
c       MaxOrd          - maximum polynomial order allowed
c
c     out:
c       Nsize           - highest polynomial order in Norder and other
c                         values that help appropriately size arrays
c
c----------------------------------------------------------------------
c
      subroutine checkorder(Ntype,Dtype,Norder,MaxOrd, Nsize)
c
      use parameters , only : MODORDER
      use node_types
c
      implicit none
      integer, intent(in)  :: Ntype
      character(len=6), intent(in)  :: Dtype
      integer, intent(in)  :: Norder(19)
      integer, intent(in)  :: MaxOrd
      integer, intent(out) :: Nsize(2)
      integer :: minp,maxp,nps(27),m,f,nordF(2),nordB(3),iprint
      logical :: confident
c
c  ...intialize printing flag
      iprint = 0
c
c  ...initialize output variable
      Nsize = 0
c
c  ...The flag confident determines whether you are confident enough to
c     NOT make a check to Norder. In that case, it is assumed all orders
c     are less than MaxOrd, and the size of the arrays assume the worst
c     case scenario (with MaxOrd).
c  ...If you are not confident enough (default), this routine will check
c     whether the values in Norder lie below MaxOrd, and will determine
c     the maximum order in Norder and size temporary arrays accordingly
c     in the shape element routines.
      confident=.FALSE.
      if (confident) then
        minp=1
        maxp=MaxOrd
      else
        select case(Ntype)
        case(BRIC,MDLB)
          m=12
          nps(1:m)=Norder(1:m)
          do f=1,6
            call decod(Norder(12+f),MODORDER,2, nordF)
            nps(m+1)=nordF(1)
            nps(m+2)=nordF(2)
            m=m+2
          enddo
          call decod(Norder(19),MODORDER,3, nordB)
          nps(m+1)=nordB(1)
          nps(m+2)=nordB(2)
          nps(m+3)=nordB(3)
          m=m+3
        case(TETR,MDLN)
          m=11
          nps(1:m)=Norder(1:m)
        case(PRIS,MDLP)
          m=11
          nps(1:m)=Norder(1:m)
          do f=1,3
            call decod(Norder(11+f),MODORDER,2, nordF)
            nps(m+1)=nordF(1)
            nps(m+2)=nordF(2)
            m=m+2
          enddo
          call decod(Norder(15),MODORDER,2, nordB(1:2))
          nps(m+1)=nordB(1)
          nps(m+2)=nordB(2)
          m=m+2
        case(PYRA,MDLD)
          m=8
          nps(1:m)=Norder(1:m)
          call decod(Norder(9),MODORDER,2, nordF)
          nps(m+1)=nordF(1)
          nps(m+2)=nordF(2)
          m=m+2
          nps(m+1:m+5)=Norder(10:14)
          m=m+5
        case(QUAD,MDLQ,RECT)
          m=4
          nps(1:m)=Norder(1:m)
          call decod(Norder(5),MODORDER,2, nordF)
          nps(m+1)=nordF(1)
          nps(m+2)=nordF(2)
          m=m+2
        case(TRIA,MDLT)
          m=4
          nps(1:m)=Norder(1:m)
        case(SEGM)
          m=1
          nps(1:m)=Norder(1:m)
        end select
c    ...print this when debugging
        if (iprint.ge.1) then
          write(*,7006) S_Type(Ntype),Dtype
 7006     format('checkorder: Type, Dtype = ',1A6,2x,1A6)
          write(*,7007) nps(1:m)
 7007     format('checkorder: nps = ',27i2)
        endif
c    ...Find the values of minp and maxp inside nps.
        minp=minval(nps(1:m))
        maxp=maxval(nps(1:m))
c    ...Determine if minp and maxp are within the bounds.
c       Otherwise stop the code as order is too low or too high.
        if (minp.lt.1) then
          write(*,7001) minp
          write(*,7002)
 7001     format('checkorder: Polynomial order = ',i3,' is less than 1')
 7002     format('            Order must be at least 1')
          stop 1
        else if (maxp.gt.MaxOrd) then
          write(*,7003,advance="no") maxp
          write(*,7004) MaxOrd
          write(*,7005) MaxOrd
 7003     format('checkorder: Polynomial order = ',i3)
 7004     format(                        ' is more than MaxOrd = ',i3)
 7005     format('            Order must be at most MaxOrd = ',i3)
          stop 1
        endif
      endif
c
c  ...The sizing of the arrays is based on the value of maxp and is
c     dependent on element type and discretization type (H1,Hcurl,...).
c     Compare these values with module parameters.
      Nsize(1)=maxp
      select case(NType)
      case(BRIC,MDLB)
        select case(Dtype)
        case('contin');Nsize(2)=(maxp+1)**3
        case('tangen');Nsize(2)=3*maxp*(maxp+1)**2
        case('normal');Nsize(2)=3*maxp**2*(maxp+1)
        case('discon');Nsize(2)=maxp**3
        end select
      case(TETR,MDLN)
        select case(Dtype)
        case('contin');Nsize(2)=(maxp+1)*(maxp+2)*(maxp+3)/6
        case('tangen');Nsize(2)=maxp*(maxp+2)*(maxp+3)/2
        case('normal');Nsize(2)=maxp*(maxp+1)*(maxp+3)/2
        case('discon');Nsize(2)=maxp*(maxp+1)*(maxp+2)/6
        end select
      case(PRIS,MDLP)
        select case(Dtype)
        case('contin');Nsize(2)=(maxp+1)*(maxp+2)*(maxp+1)/2
        case('tangen');Nsize(2)=maxp*(maxp+2)*(maxp+1)
     .                          +(maxp+1)*(maxp+2)*maxp/2
        case('normal');Nsize(2)=maxp*(maxp+2)*maxp
     .                          +maxp*(maxp+1)*(maxp+1)/2
        case('discon');Nsize(2)=maxp*(maxp+1)*maxp/2
        end select
      case(PYRA,MDLD)
        select case(Dtype)
        case('contin');Nsize(2)=5+8*(maxp-1)+(maxp-1)**2
     .                          +2*(maxp-2)*(maxp-1)+(maxp-1)**3
        case('tangen');Nsize(2)=8*maxp+2*maxp*(maxp-1)
     .                          +4*maxp*(maxp-1)+3*(maxp-1)**2*maxp
        case('normal');Nsize(2)=maxp**2+2*maxp*(maxp+1)
     .                          +3*(maxp-1)*maxp**2
        case('discon');Nsize(2)=maxp**3
        end select
      case(QUAD,MDLQ,RECT)
        select case(Dtype)
        case('contin');         Nsize(2)=(maxp+1)**2
        case('tangen','normal');Nsize(2)=2*maxp*(maxp+1)
        case('discon');         Nsize(2)=maxp**2
        end select
      case(TRIA,MDLT)
        select case(Dtype)
        case('contin');         Nsize(2)=(maxp+1)*(maxp+2)/2
        case('tangen','normal');Nsize(2)=maxp*(maxp+2)
        case('discon');         Nsize(2)=maxp*(maxp+1)/2
        end select
      case(SEGM)
        select case(Dtype)
        case('contin');                  Nsize(2)=maxp+1
        case('tangen','normal','discon');Nsize(2)=maxp
        end select
      end select
c
c  ...print this when debugging
      if (iprint.ge.1) then
        write(*,7008) Nsize(1),Nsize(2)
 7008   format('checkorder: Nsize(1:2) = ',i2,i5)
      endif
c
c
      end subroutine checkorder
