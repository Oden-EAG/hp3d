c----------------------------------------------------------------------
c
c   routine name       - GMP_err
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - Oct, 2003
c
c   purpose            - routine evaluates H01-error of the geometry
c                        interpolation error
c
c   usage              - call GMP_err(Error)
c
c   arguments :
c     out:
c             Error    - the error
c
c----------------------------------------------------------------------
c
      subroutine GMP_err(Error)
c
      use data_structure
c
#include "syscom.blk"
#include "cinout.blk"
#include "cpmax.blk"
c
c
      dimension xsub(3,8),zeta(3,8)
      dimension norder(19),xnod(3,(NPMAX+1)**3)
c
c  explanation of variables:
c     xsub   - reference coordinates of the element wrt the
c              block, the element belongs to
c     zeta   - master element coordinates of the subelement
c     norder - order of approximation for an element
c     xnod   - contains geometry dof for an element
c
c
c----------------------------------------------------------------------
c
c
      iprint=1
c
c  ...set the tollerance constant
      eps = 1.d-4
c
c  ...initiate Data
      nel=0
      kind=0
c
      Error = 0.d0
c

c
c  ...loop through the body elements
      mdle = 0
      do  iel=1,NRELEB
c
c     ...initiate the error
         derror = 0.d0
c
c
         call nelconb(mdle,mdle)
c
c
c     ...determine order of approximation
         call findapb(Mdle,norder)
c
c     ...determine nodes coordinates
         call nodcorb(Mdle,xnod)
c
c     ...determine the geometrical block for the father element and the
c        corresponding reference coordinates
         call refelb(Mdle,6,nick,xsub)
c
c
c     ...initiate number of subdivisions
         nsub = 1
c
c     ...increment number of subdivisions
 10      nsub = nsub*2
         if (nsub.gt.10) go to 110
c
         d = 1.d0/nsub
c     ...initiate the new element error
         derrnew = 0.d0
c
c     ...loop through subelements...
         do k=1,nsub
            do j=1,nsub
               do i= 1,nsub
c            ....the first vertex
                  zeta(1,1) =  (i-1)*d
                  zeta(2,1) =(j-1)*d
c             ....the second vertex
                  zeta(1,2) = i*d
                  zeta(2,2) = (j-1)*d
c             ....the third vertex
                  zeta(1,3) = i*d
                  zeta(2,3) = j*d
c             ....the fourth vertex
                  zeta(1,4) = (i-1)*d
                  zeta(2,4) = j*d
c
                  zeta(3,1:4) = (k-1)*d
c             ....from fifth to eighth vertex
                  do ip = 1,4
                     zeta(1:2,ip+4) = zeta(1:2,ip)
                  enddo
                  zeta(3,5:8) = k*d
c
c
c        .........accumulate for the error
                  call eleb_err(nick,norder,xnod,xsub,zeta,derr)
                  derrnew = derrnew + derr
               enddo
            enddo
         enddo
c
c     ...compare the errors
         diff = abs(derrnew - Derror)
         derror = derrnew
         if (diff.lt.eps) then
            go to 100
         else
            go to 10
         endif
c
c     ...issue a warning...
 110     write(*,*) 'NO CONVERGENCE IN GMP_ERR !!'
         write(*,*) 'nsub,error = ',nsub,derror
         stop
c
 100     continue
         write(*,7001) mdle, derror
 7001    format('GMP_err: mdle,derror = ',i5,e12.5)
         Error = Error + derror
      enddo
c
      Error = log(error)
c
      return
      end

