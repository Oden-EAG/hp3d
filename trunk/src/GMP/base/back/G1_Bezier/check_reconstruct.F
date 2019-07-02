c----------------------------------------------------------------------
c
c   routine name       - check_reconstruct
c
c----------------------------------------------------------------------
c
c   computer           - machine independent
c
c   latest revision    - May 07
c
c   purpose            - Routine verifies orientation and determinants
c                        for a hexahedral mesh with some points on
c                        an oriented surface
c
c   arguments            none
c
c   required  routines -
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
c  ...normal for a face vertex
      dimension rnf(3)
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
c  ...neigboring vertices for a vertex in the cube
      dimension neig_vert(3,8)
c
c  ...a flag indicating whether a point has been visited
      integer, dimension(:), pointer:: nvisit
c
c  ...a work space
      dimension temp(3),temp1(3),iface(6),
     .          aa(3,3),bb(3),uu(3),ipiv(3)
c
      data norder /1,1,1,1,1,1,1,1,1,1,1,1,
     .             11,11,11,11,11,11, 111/
c
      data neig_vert/2,4,5, 3,1,6, 4,2,7, 1,3,8,
     .               8,6,1, 5,7,2, 6,8,3, 7,5,4/
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
      data pi /3.14159265358979312d0/
c
      theta_max = pi/2.d0*7.d0/9.d0
c
c
      iprint=0
      allocate(nvisit(NRPOINT))
      nvisit(1:NRPOINT) = 0
c
      nrsplit=0
c  ...loop over tets
      do ih=1,NRTETRA
c  .....loop over tet's vertices
        do iv=1,4
          np = TETRAS(ih)%VertNo(iv)
          if (POINTS(np)%Type.eq.'CoorNrm') then
            nflag(iv)=1
          else
            nflag(iv)=0
          endif
c  .....end loop over vertices
        enddo
c
c  .....check if only one face is on the reconstructed surface
        ii=0
c  .....loop over tet's faces
        do if=1,4
          call decode(TETRAS(ih)%FigNo(if), nr,lab)
          idec=1
          do ifv=1,4
            iv = nface_vert(ifv,if)
            if (nflag(iv).eq.0) idec=0
          enddo
c
c  .......check for neighbors, there should be only one
          if ((TRIANGLES(nr)%BlockNo(1).ne.0).and.
     .        (TRIANGLES(nr)%BlockNo(2).ne.0)) idec=0
          if (idec.eq.1) then
            ii=ii+1
            iface(ii) = if
          endif
        enddo
        if (ii.gt.1) then
          write(*,8001) ih,0,ii
 8001     format(i6,2x,6i3)
          nrsplit=nrsplit+1
        endif
      enddo
      if (nrsplit.gt.0) then
        write(*,8002) nrsplit
 8002   format('check_reconstruct: NUMBER OF HEXAS TO SPLIT = ',i4)
        stop 1
      endif
c
c  ...loop through the hexas
      do ih=1,NRTETRA
        if (iprint.eq.5) then
          write(*,5038) ih, TETRAS(ih)%VertNo(1:8)
 5038     format('check_reconstruct: ih, VERTICES = ',
     .           i5,3x,8i5)
        endif
c
c  .....determine vertex coordinates
        do iv=1,8
          np = TETRAS(ih)%VertNo(iv)
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
          call shapeHb(xi(1:3,iv),norder, nrdof,vshap,dvshap)
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
            call pause
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
c
c  .........face on the reconstructed surface
c
c  .........loop through the face vertices
            do ifv=1,4
              iv = nface_vert(ifv,if)
              np = TETRAS(ih)%VertNo(iv)
c
c  ...........compute components of the normal wrt local vertex
c             coordinates
              do k=1,3
                ivk = neig_vert(k,iv)
                aa(1:3,k) = xv(1:3,ivk) - xv(1:3,iv)
              enddo
              bb(1:3) = POINTS(np)%Rdata(4:6)
              if (iprint.eq.5) then
                do i=1,3
                  write(*,7027) i,aa(i,1:3),bb(i)
 7027             format('i,aa,bb = ',i2,2x,3f8.3,2x,f8.3)
                enddo
              endif
              call decomp(3,3,aa,ipiv,iflag)
              if (iflag.eq.1) then
                write(*,*) 'check_reconstruct: SINGULAR MATRIX'
                stop 1
              endif
              call gauss2(3,3,aa,ipiv,bb, uu)
              if (iprint.eq.5) then
                write(*,7007) np,ih,iv,if,ifv,uu
 7007           format('check_reconstruct: np,ih,iv,if,ifv,uu = ',
     .                 i5,i5,3i2,2x,3f8.3)
                call pause
              endif
c
              if ((uu(1).ge.0.d0).and.(uu(2).ge.0.d0).and.
     .          (uu(3).ge.0.d0)) then
                if (nvisit(np).eq.0) then
                  write(*,7005) np,ih,if,ifv
 7005             format('check_reconstruct: MODIFYING NORMAL FOR ',
     .                   'np,ih,if,ifv =',i5,i5,2i2)
                  nvisit(np)=1
                  POINTS(np)%Rdata(4:6) = -POINTS(np)%Rdata(4:6)
                else
                  write(*,7006) np,ih,if,ifv
 7006             format('check_reconstruct: CANNOT REVERSE NORMAL ',
     .                   'FOR np,ih,if,ifv =',i5,i5,2i2)
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
ccc      return
c
c-----------------------------------------------------------------------
c
c  ...geometrical corrections...........................................
 10   continue
c
c  ...number of corrected normals
      normals_corrected=0
c
c  ...loop through the hexas
      do ih=1,NRTETRA
c
c  .....determine vertex coordinates
        do iv=1,8
          np = TETRAS(ih)%VertNo(iv)
          xv(1:3,iv) = POINTS(np)%Rdata(1:3)
          if (POINTS(np)%Type.eq.'CoorNrm') then
            nflag(iv)=1
          else
            nflag(iv)=0
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
c
c  .........face on the reconstructed surface
c
c  .........loop through the face vertices
            do ifv=1,4
              iv = nface_vert(ifv,if)
              np = TETRAS(ih)%VertNo(iv)
c
c  ...........enforce minimum angle condition
c
c  ...........look for the vertex connected to the face vertex
c             that it is NOT on the face
              do k=1,3
                ivk = neig_vert(k,iv)
                call locate(ivk,nface_vert(1:4,if),4, ii)
                if (ii.eq.0) exit
              enddo
c
c  ...........edge vector
              temp(1:3) = xv(1:3,iv) - xv(1:3,ivk)
              call normalize(temp)
              call scalar_product(temp,POINTS(np)%Rdata(4:6), s)
              if (s.lt.cos(theta_max)) then
c
c  .............modify the component along the edge vector
                temp1(1:3) = POINTS(np)%Rdata(4:6) - s*temp(1:3)
                call normalize(temp1)
                POINTS(np)%Rdata(4:6) = sin(theta_max*.99d0)*temp1(1:3)
     .                                + cos(theta_max*.99d0)*temp(1:3)
                write(*,7047) np,ih,POINTS(np)%Rdata(1:3)
 7047           format('input_reconstruct: CORRECTED NORMAL AT np,ih =',
     .                  1x,2i5,3f8.3)
                normals_corrected = normals_corrected+1
              endif
            enddo
          endif
        enddo
c
c  ...end of loop through hexas
      enddo
      if (normals_corrected.gt.0) then
        write(*,*) 'normals_corrected = ',normals_corrected
        call pause
        go to 10
      endif
c
c
      end


