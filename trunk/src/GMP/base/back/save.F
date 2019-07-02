c----------------------------------------------------------------------
c
c   routine name       - check_reconstruct
c
c----------------------------------------------------------------------
c
c   latest revision    - May 07
c
c   purpose            - Routine verifies orientation and determinants
c                        for a hexahedral mesh with some points on
c                        an oriented surface
c
c   arguments            none
c
c---------------------------------------------------------------------
c
      subroutine check_reconstruct
c
      use control
      use GMP
#include "syscom.blk"
#include "cinout.blk"
c
c  ...vertex coordinates, normals, surface flags
      dimension xv(3,8),rn(3,8),nflag(8)
c
c  ...normals for a face
      dimension rnf(3,4)
c
c  ...order of approximation (linear)
      dimension norder(19)
c
c  ...shape functions and their derivatives
      dimension vshap(8),dvshap(3,8)
c
c  ...master element coordinates of vertices
      dimension xi(1:3,8)
c
c  ...derivatives wrt reference coordinates
      dimension dxdxi(3,3,8)
c
c  ...vertex points on a face, face coordinates
      dimension nface_vert(4,6), nface_axes(2,6)
c
c  ...a flag indicating whether a point has been visited
      integer, dimension(:), pointer:: nvisit
c
      data norder /1,1,1,1,1,1,1,1,1,1,1,1,
     .             11,11,11,11,11,11, 111/
c
      data nface_vert /1,2,3,4, 5,6,7,8, 1,2,6,5,
     .                 2,3,7,6, 3,4,8,7, 1,4,8,5/
c
      data nface_axes /1,2, 1,2, 1,3, 2,3, 1,3, 2,3/
c
      data xi /0.d0,0.d0,0.d0, 1.d0,0.d0,0.d0,
     .         1.d0,1.d0,0.d0, 0.d0,1.d0,0.d0,
     .         0.d0,0.d0,1.d0, 1.d0,0.d0,1.d0,
     .         1.d0,1.d0,1.d0, 0.d0,1.d0,1.d0/
c
      iprint=0
      allocate(nvisit(NRPOINT))
      nvisit(1:NRPOINT) = 0
 10   continue
c
c  ...loop through the hexas
      do ih=1,NRHEXAS
c
c  .....determine vertex coordinates
        do iv=1,8
          np = HEXAS(ih)%VertNo(iv)
          xv(1:3,iv) = POINTS(np)%Rdata(1:3)
          if (POINTS(np)%Type.eq.'CoorNrm') then
            nflag(iv)=1
            rn(1:3,iv) = POINTS(np)%Rdata(4:6)
          else
            nflag(iv)=0
            rn(1:3,iv) = 0.d0
          endif
        enddo
c
c  .....check determinants at vertices
        do iv=1,8
          call shape3b(xi(1:3,iv),norder, nrdof,vshap,dvshap)
          do i=1,3
            dxdxi(1:3,i,iv) = 0.d0
            do k=1,nrdof
              dxdxi(1:3,i,iv) = dxdxi(1:3,i,iv)
     .                        + xv(1:3,k)*dvshap(i,k)
            enddo
          enddo
          call mixed_product(dxdxi(1:3,1,iv),dxdxi(1:3,2,iv),
     .                       dxdxi(1:3,3,iv), rjac)
          if (iprint.eq.1) then
            write(*,7004) ih,iv,rjac
 7004       format('check_reconstruct: ih,iv,rjac = ',
     .                  i5,i2,e12.5)
          endif
          if (rjac.lt.0.d0) then
            write(*,7001) ih,iv
 7001       format('check_reconstruct: NEGATIVE JACOBIAN at ih,iv = ',
     .             i5,i2)
          endif
        enddo
c
c  .....check normals for faces
        do if=1,6
          idec=1
          do ifv=1,4
            iv = nface_vert(ifv,if)
            if (nflag(iv).eq.0) idec=0
          enddo
          if (idec.eq.1) then
            do ifv=1,4
              iv = nface_vert(ifv,if)
              call mixed_product(dxdxi(1:3,nface_axes(1,if),iv),
     .                           dxdxi(1:3,nface_axes(2,if),iv),
     .                           rn(1:3,iv), det)
              if (ifv.eq.1) rsign = sign(1.d0,det)
              if (iprint.eq.1) then
                write(*,7003) ih,if,ifv,rsign*det
 7003           format('check_reconstruct: ih,if,ifv,rsign*det = ',
     .                  i5,i2,i2,e12.5)
              endif
              if (rsign*det.lt.0.d0) then
                write(*,7002) ih,if,ifv
 7002           format('check_reconstruct: INCONSISTENT NORMAL FOR ',
     .                 'ih,if,ifv = ',i5,2i2)
c
c  .............first vertex point
                iv1 = nface_vert(1,if); np1 = HEXAS(ih)%VertNo(iv1)
                iv = nface_vert(ifv,if); np = HEXAS(ih)%VertNo(iv)
                if (nvisit(np1).eq.0) then
                  write(*,7005) np1
 7005             format('check_reconstruct: MODIFYING NORMAL FOR np =',
     .                    i5)
                  POINTS(np1)%Rdata(4:6) = - POINTS(np1)%Rdata(4:6)
                  nvisit(np1) = 1
                  go to 10
                elseif (nvisit(np).eq.0) then
                  write(*,7005) np
                  POINTS(np)%Rdata(4:6) = - POINTS(np)%Rdata(4:6)
                  nvisit(np) = 1
                  go to 10
                else
                  write(*,7006) rsign*det
 7006             format('check_reconstruct: CANNOT RECOVER...',e12.5)
                  stop 1
                endif
              endif
            enddo
          endif
        enddo
c
c  ...end of loop through hexas
      enddo
      deallocate(nvisit)
c
c
      end
















c
c  ...........check if the normal is pointing outside of the hexa
              step = 1.d0; inside=0
              do iter=1,10
                temp(1:3) = xv(1:3,iv) -


              call mixed_product(dxdxi(1:3,nface_axes(1,if),iv),
     .                           dxdxi(1:3,nface_axes(2,if),iv),
     .                           rn(1:3,iv), det)
              rsign = rsign + det
            enddo
            if (rsign.gt.0.d0) then
              sgn = 1.d0
            elseif (rsign.lt.0.d0) then
              sgn = -1.d0
            else
              write(*,*) 'check_reconstruct: SIGN PROBLEM'
              stop 1
            endif
            do ifv=1,4
              iv = nface_vert(ifv,if)
              rnf(1:3) = sgn*rn(1:3,iv)
c
c  ...........project onto the face normal
              call cross_product(dxdxi(1:3,nface_axes(1,if),iv),
     .                           dxdxi(1:3,nface_axes(2,if),iv),
     .                           temp)
              call scalar_product(temp,temp, s)
              s = sqrt(s)
              temp(1:3) = temp(1:3)/s
              call scalar_product(rnf,temp, s)
              if (s.lt.0.d0) then
                np = HEXAS(ih)%VertNo(iv)
                write(*,7005) np,ih,if
 7005           format('check_reconstruct: MODIFYING NORMAL FOR ',
     .                 'np,ih,if =',i5,i5,i2)
c
c  .............in plane component
                rnf(1:3) = rnf(1:3) - s*temp(1:3)
                call scalar_product(rnf(1:3),rnf(1:3), s)
                s = sqrt(s)
                rnf(1:3) = rnf(1:3)/s
c
c  .............combine the in plane component with the new
c               normal component
                POINTS(np)%Rdata(4:6) = (cos(theta_min)*rnf(1:3)
     .                                + sin(theta_min)*temp(1:3))*sgn
                nvisit(np) = nvisit(np)+1
                if (nvisit(np).gt.10) then
                  write(*,7006)
 7006             format('check_reconstruct: CANNOT RECOVER...')
                  cycle
                endif
                go to 10
              endif
            enddo
          endif
        enddo



c----------------------------------------------------------------------
c
c   routine name       - correct_input
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - May 07
c
c   purpose            - Routine corrects Jesicca's input file by
c                        adding the points on the lower section to
c                        the list of surface points
c
c   arguments            none
c
c   required  routines -
c
c---------------------------------------------------------------------
c
      subroutine correct_input
c
      use control
      use GMP
#include "syscom.blk"
#include "cinout.blk"
c
      dimension temp(3)
c
      nr_added=0
c
c  ...loop through all rectangles
      do nr=1,NRRECTA
c
c  .....skip if not on the surface
        if (RECTANGLES(nr)%BlockNo(2).ne.0) cycle
c
c  .....loop through the rectangle vertices
        do iv=1,4
          np = RECTANGLES(nr)%VertNo(iv)
          if (POINTS(np)%Type.ne.'CoorNrm') then
            nr_added = nr_added + 1
            POINTS(np)%Type = 'CoorNrm'
            temp(1:3) = POINTS(np)%Rdata(1:3)
            deallocate(POINTS(np)%Rdata)
            allocate(POINTS(np)%Rdata(1:6))
            POINTS(np)%Rdata(1:3) = temp(1:3)
c
c  .........set normal to the unit vector of the position vector
            call normalize(temp)
            POINTS(np)%Rdata(4:6) = temp(1:3)
          endif
        enddo
c
c  ...end of loop through rectangles
      enddo
c
      write(*,7001) nr_added
 7001 format('correct_input: ',i4,' SURFACE POINTS ADDED')
      call pause
c
      end
