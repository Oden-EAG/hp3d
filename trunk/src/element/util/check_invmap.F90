subroutine check_invmap
!
!--------------------------------------------------------------------------
      use data_structure3D
      use control , only : EXGEOM
!--------------------------------------------------------------------------
!
      implicit none
!
      integer :: mdle,nel,ifig,nint,l,k,i,idec,nrdofH
      real*8  :: s
      character (len=4) :: ftype
      real*8,dimension(3)           :: xi,xi_aux,x,temp
      real*8,dimension(3,MAXbrickH) :: xnod
      real*8,dimension(2)           :: t
      real*8,dimension(3,3)         :: dxdxi
      real*8,dimension(3,2)         :: dxdt,dxidt
!
      integer,dimension(12) :: nedge_orient
      integer,dimension(6)  :: nface_orient
      integer,dimension(19) :: norder
      integer,dimension(5)  :: nordf
!
      real*8,dimension(  MAXbrickH) :: vshapH
      real*8,dimension(3,MAXbrickH) :: dvshapH
!
!  ...integration points      
      real*8,dimension(3,MAX_NINT3) :: xiloc
      real*8,dimension(  MAX_NINT3) :: wxi
      real*8,dimension(2,MAXquadH)  :: tloc
      real*8,dimension(  MAXquadH)  :: wt
!      
      real*8, parameter :: eps=1.d-10
      integer :: iprint_invmap,iprint
      common /cinvmap/ iprint_invmap
!
!--------------------------------------------------------------------------
!
      iprint=1
      write(*,*)'checking inverse map... '
!      
!  ...loop over active elements      
      do nel=1,NRELES
        mdle = ELEM_ORDER(nel)
!
        if (iprint.eq.1) then
          write(*,9999)mdle,NODES(mdle)%type
 9999     format(' mdle,type = ',i7,2x,a4)          
        endif
!
!  .....order, orientations, gdof's        
        call find_order( mdle, norder)
        call find_orient(mdle, nedge_orient,nface_orient)
        if (EXGEOM.eq.0) call nodcor(mdle, xnod)
!
!--------------------------------------------------------------------------        
!  E L E M E N T   I N T E R I O R                                        |
!--------------------------------------------------------------------------        
        call set_3Dint(NODES(Mdle)%type,norder, nint,xiloc,wxi)
!
!  .....loop over interior integration points
        do l=1,nint
!
!  .......Gauss point        
          xi(1:3)=xiloc(1:3,l)
!
          select case(EXGEOM)
!  .......parametric element          
          case(0)
!  .........shape functions
            call shape3H(NODES(Mdle)%type,xi,norder,nedge_orient, &
                         nface_orient, nrdofH,vshapH,dvshapH)
!
!  .........accumulate
            x(1:3)=0.d0 ; dxdxi(1:3,1:3)=0.d0 
            do k=1,nrdofH
              x(1:3) = x(1:3) + xnod(1:3,k)*vshapH(k)
              do i=1,3
                dxdxi(1:3,i) = dxdxi(1:3,i) + xnod(1:3,k)*dvshapH(i,k)
              enddo
            enddo
!  .......exact geometry element            
          case(1)
            call exact_geom(mdle,xi, x,dxdxi)
          endselect
!          
!  .......compute inverse map          
          call invmap(NODES(mdle)%type,x,norder,nedge_orient, &
                      nface_orient,xnod ,idec,xi_aux)
!          
!  .......check                        
          temp(1:3)=xi_aux(1:3)-xi(1:3)
          call norm(temp, s)
          if (s.gt.eps) then
            write(*,*)'Interior point FAIL'
            write(*,6665)mdle,NODES(mdle)%type
 6665       format(' mdle,type = ',i7,2x,a4)           
            write(*,6666)xi(    1:3)
 6666       format(' xi     = ',3(e12.5,2x))
            write(*,6667)xi_aux(1:3)
 6667       format(' xi_aux = ',3(e12.5,2x))
          endif
!
!  .....end loop over interior integration points          
        enddo
!
!--------------------------------------------------------------------------        
!  E L E M E N T   F A C E S                                              |        
!--------------------------------------------------------------------------        
!  .....loop over element faces
        do ifig=1,nface(NODES(mdle)%type)
!
!  .......set up the element quadrature
          ftype=face_type(NODES(mdle)%type,ifig)
!
!  .......determine order for the face
          call face_order(NODES(mdle)%type,ifig,norder, nordf)
          call set_2Dint(ftype,nordf, nint,tloc,wt)
!
!  .......loop through face integration points
          do l=1,nint
!
!  .........Gauss point
            t(1:2)=tloc(1:2,l)
!
!  .........determine the master element coordinates
            call face_param(NODES(mdle)%type,ifig,t, xi,dxidt)
!
            select case(EXGEOM)
!  .........parametric element            
            case(0)
!  ...........derivatives and values of the shape functions
              call shape3H(NODES(mdle)%type,xi,norder,nedge_orient, &
                           nface_orient, nrdofH,vshapH,dvshapH)
!
!  ...........accumulate
              x(1:3)=0.d0 ; dxdxi(1:3,1:3)=0.d0 
              do k=1,nrdofH
                x(1:3) = x(1:3) + xnod(1:3,k)*vshapH(k)
                do i=1,3
                  dxdxi(1:3,i) = dxdxi(1:3,i) + xnod(1:3,k)*dvshapH(i,k)
                enddo
              enddo
!  .........exact geometry element              
            case(1)
              call exact_geom(mdle,xi, x,dxdxi)
            endselect
!
!  .........compute inverse map          
            call invmap(NODES(mdle)%type,x,norder,nedge_orient, &
                        nface_orient,xnod ,idec,xi_aux)
!            
            temp(1:3)=xi_aux(1:3)-xi(1:3)
            call norm(temp, s)
            if (s.gt.eps) then
              write(*,6668)ifig
 6668         format(' Face point FAIL, ifig = ',i1)              
              write(*,6665)mdle,NODES(mdle)%type
              write(*,6666)xi(    1:3)
              write(*,6667)xi_aux(1:3)
            endif
!
!  .......end of loop over integration points          
          enddo
!          
!  .....end of loop over element faces
        enddo
!
      enddo
!
      write(*,*)'DONE!'      
!
!      
end subroutine check_invmap
