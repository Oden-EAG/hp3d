c----------------------------------------------------------------------
c                                                                     
c     routine name      - set_3Dint
c                                                                     
c----------------------------------------------------------------------
c                                                                     
c     latest revision:  - Sep 17
c                                                                     
c     purpose:          - routine sets up quadrature data for a standard
c                         3D element, accounting for different element
c                         types and orders of approximation
c                                                                    
c     arguments:                                                     
c                                                                     
c     in:              
c             Type      - element type
c             Norder    - order of approximation
c        
c     out:    
c             Nint      - number of integration points
c             Xiloc     - integration points
c             Waloc     - weights
c             
c----------------------------------------------------------------------
c
      subroutine set_3Dint(Type,Norder, Nint,Xiloc,Waloc)
c
      use parameters, only : MAXP,MAX_NINT3
#include "syscom.blk"
c
      character(len=4),intent(in) :: Type
      dimension Norder(19),Xiloc(3,MAX_NINT3),Waloc(MAX_NINT3)
c
      call set_3Dint_aux(Type,Norder,MAXP,MAX_NINT3, Nint,Xiloc,Waloc)
c
      end subroutine set_3Dint
c
c----------------------------------------------------------------------
c                                                                     
c     routine name      - set_3Dint_DPG
c                                                                     
c----------------------------------------------------------------------
c                                                                     
c     latest revision:  - Sep 17
c                                                                     
c     purpose:          - routine sets up quadrature data for a DPG
c                         3D element, accouting for different element
c                         types and orders of approximation
c                                                                    
c     arguments:                                                     
c                                                                     
c     in:              
c             Type      - element type
c             Norder    - order of approximation
c        
c     out:    
c             Nint      - number of integration points
c             Xiloc     - integration points
c             Waloc     - weights
c             
c----------------------------------------------------------------------
c
      subroutine set_3Dint_DPG(Type,Norder, Nint,Xiloc,Waloc)
c
      use parametersDPG, only : MAXPP,MAXNINT3ADD
#include "syscom.blk"
c
      character(len=4),intent(in) :: Type
      dimension Norder(19),Xiloc(3,MAXNINT3ADD),Waloc(MAXNINT3ADD)
c
      call set_3Dint_aux(Type,Norder,MAXPP,MAXNINT3ADD, 
     .                   Nint,Xiloc,Waloc)
c
      end subroutine set_3Dint_DPG
c
c----------------------------------------------------------------------
c                                                                     
c     routine name      - set_3Dint_aux
c                                                                     
c----------------------------------------------------------------------
c                                                                     
c     latest revision:  - Sep 17
c                                                                     
c     purpose:          - routine sets up quadrature data for a 3D
c                         element, accouting for different element
c                         types and orders of approximation
c                                                                    
c     arguments:                                                     
c                                                                     
c     in:              
c             Type      - element type
c             Norder    - order of approximation
c             Maxp      - maximum p
c             Max_nint3 - maximum number of integration points
c        
c     out:    
c             Nint      - number of integration points
c             Xiloc     - integration points
c             Waloc     - weights
c             
c----------------------------------------------------------------------
c
      subroutine set_3Dint_aux(Type,Norder,Maxp,Max_nint3, 
     .                         Nint,Xiloc,Waloc)
c
      use parameters,        only : MODORDER
      use control,           only : INTEGRATION
      use gauss_quadrature,  only : INITIALIZED,
     .                              NSELECT     ,NRGAUPO,
     .                              NSELECT_TETS,NRGAUPO_TETS,
     .                              XIGAUSS     ,WAGAUSS,
     .                              XIGAUS1     ,WAGAUS1,
     .                              XIGAUSS_TETS,WAGAUSS_TETS
c      
#include "syscom.blk"
c
      character(len=4),intent(in) :: Type
      dimension Norder(19),Xiloc(3,Max_nint3),Waloc(Max_nint3)
      integer :: nordhv(2), nordxyz(3)
c
c----------------------------------------------------------------------

c  ...initialize if needed
      if (.NOT. INITIALIZED) call init_gauss_quadrature
c      
      iprint=0
      if (iprint.eq.1) then
        write(*,7001) Type,Norder
 7001   format('set_3Dint: Type, Norder = ',a4,2x,19i4)
      endif
c
      select case(Type)
c
c======================================================================      
c  PRISM                                                              |
c======================================================================      
      case('mdlp','pris')
c
c  .....determine order of approximation              
        nordh=0 ; nordz=0
        do i=1,15
          select case(i)
c         horizontal edges & faces
          case(1,2,3,4,5,6, 10,11) ; nordh=max(nordh,Norder(i))
c         vertical edges             
          case(7,8,9)              ; nordz=max(nordz,Norder(i))
c         vertical faces & interior
          case(12,13,14, 15)
            call decod(Norder(i),MODORDER,2, nordhv)
            nordh=max(nordh,nordhv(1)) ; nordz=max(nordz,nordhv(2))
          endselect
        enddo
c
c  .....account for over-integration
        nordh=min(nordh+INTEGRATION,Maxp)
        nordz=min(nordz+INTEGRATION,Maxp)
c        
c  .....compute number of integration points         
        kint =NSELECT(nordh)
        nintx=NRGAUPO(kint)
        nintz=nordz+1
        Nint =nintx*nintz
c
c  .....compute integration points and weights        
        l=0
        do l2=1,nintz ; do l1=1,nintx
          l=l+1
          Xiloc(1:2,l)=XIGAUSS(2:3,l1,kint)
          Xiloc(  3,l)=XIGAUS1(l2,nintz)
          Waloc(    l)=WAGAUSS(l1,kint)/2.d0*WAGAUS1(l2,nintz)
        enddo ; enddo
c
c======================================================================      
c  BRICK                                                              |
c======================================================================      
      case('mdlb','bric')
c
c  .....determine order of approximation
        nordx=0 ; nordy=0 ; nordz=0
        do i=1,19
          select case(i)
          case(1,3,5,7)
            nordx = max(nordx,Norder(i))
          case(2,4,6,8)
            nordy = max(nordy,Norder(i))
          case(9,10,11,12)
            nordz = max(nordz,Norder(i))
          case(13,14)
            call decod(Norder(i),MODORDER,2, nordhv)
            nordx = max(nordx,nordhv(1))
            nordy = max(nordy,nordhv(2))
          case(15,17)
            call decod(Norder(i),MODORDER,2, nordhv)
            nordx = max(nordx,nordhv(1))
            nordz = max(nordz,nordhv(2))
          case(16,18)
            call decod(Norder(i),MODORDER,2, nordhv)
            nordy = max(nordy,nordhv(1))
            nordz = max(nordz,nordhv(2))
          case(19)
            call decod(Norder(i),MODORDER,3, nordxyz)
            nordx = max(nordx,nordxyz(1))
            nordy = max(nordy,nordxyz(2))
            nordz = max(nordz,nordxyz(3))
          end select
        enddo
c
c  .....account for overintegration
        nordx=min(nordx+INTEGRATION,Maxp)
        nordy=min(nordy+INTEGRATION,Maxp)
        nordz=min(nordz+INTEGRATION,Maxp)
c
c  .....compute number of integration points        
        nintx=nordx+1
        ninty=nordy+1
        nintz=nordz+1
        Nint=nintx*ninty*nintz
c
c  .....compute integration points and weights        
        l=0
        do l3=1,nintz ; do l2=1,ninty ; do l1=1,nintx
          l=l+1
          Xiloc(1,l)=XIGAUS1(l1,nintx)
          Xiloc(2,l)=XIGAUS1(l2,ninty)
          Xiloc(3,l)=XIGAUS1(l3,nintz)
          Waloc(  l)=WAGAUS1(l1,nintx)*WAGAUS1(l2,ninty)*
     .               WAGAUS1(l3,nintz)
        enddo ; enddo ; enddo
c
c======================================================================      
c  TETRAHEDRON                                                        |
c======================================================================      
      case('mdln','tetr')
c
c  .....determine order of approximation
        nord=0
        do i=1,11 ; nord=max(nord,Norder(i)) ; enddo
c
c  .....account for overintegration
        nord=min(nord+INTEGRATION,Maxp)
c
c  .....retrive number of integration points        
        kint=NSELECT_TETS(nord)
        Nint=NRGAUPO_TETS(kint)
c
c  .....retrive integration points and weights        
        Xiloc(1:3,1:Nint)=XIGAUSS_TETS(1:3,1:Nint,kint)
        Waloc(    1:Nint)=WAGAUSS_TETS(    1:Nint,kint)
c
c======================================================================      
c  PYRAMID                                                            |
c======================================================================      
      case('mdld','pyra')
c              
c  .....determine order of approximation
        nordx=0 ; nordy=0 ; nordz=0
        do i=1,14
          select case(i)
          case(1,3)
            nordx = max(nordx,Norder(i))
          case(2,4)
            nordy = max(nordy,Norder(i))
          case(5,6,7,8)
            nordz = max(nordz,Norder(i))
          case(9)
            call decod(Norder(i),MODORDER,2, nordhv)
            nordx = max(nordx,nordhv(1))
            nordy = max(nordy,nordhv(2))
          case(10,12)
            nordx = max(nordx,Norder(i))
            nordz = max(nordz,Norder(i))
          case(11,13)
            nordy = max(nordy,Norder(i))
            nordz = max(nordz,Norder(i))
          case(14)
            nordx = max(nordx,Norder(i))
            nordy = max(nordy,Norder(i))
            nordz = max(nordz,Norder(i))
          end select
        enddo
c
c  .....account for overintegration
        nordx=min(nordx+INTEGRATION,Maxp) 
        nordy=min(nordy+INTEGRATION,Maxp)
        nordz=min(nordz+INTEGRATION,Maxp)
c
c  .....compute number of integration points        
        nintx=nordx+1
        ninty=nordy+1
        nintz=nordz+2
        Nint=nintx*ninty*nintz
c
c  .....compute integration points and weights        
        l=0
        do l3=1,nintz ; do l2=1,ninty ; do l1=1,nintx
          l=l+1
          Xiloc(3,l)=XIGAUS1(l3,nintz)
          factor    =1.d0-Xiloc(3,l)
          Xiloc(1,l)=XIGAUS1(l1,nintx)*factor
          Xiloc(2,l)=XIGAUS1(l2,ninty)*factor
          Waloc(  l)=WAGAUS1(l1,nintx)*WAGAUS1(l2,ninty)
     .              *WAGAUS1(l3,nintz)*factor**2
        enddo ; enddo ; enddo
c
      case default
        write(*,*) 'set_3Dint:Mdle,Type = ',Mdle,Type
        stop
      endselect
c
c
      end subroutine set_3Dint_aux
c
c
c      
c--------------------------------------------------------------------
c                                                                     
c     routine name      - test_set_3Dint
c                                                                     
c-------------------------------------------------------------------- 
c                                                                     
c     latest revision:  - Dec 08
c                                                                     
c     purpose:          - routine test quadrature data for a 3D
c                         element
c                                                                    
c     arguments         - none
c                                                                     
c-----------------------------------------------------------------------
c
      subroutine test_set_3Dint
c
      use parameters
      use parametersDPG
#include "syscom.blk"
c
      dimension Norder(19),Xiloc(3,MAX_NINT3),Waloc(MAX_NINT3)
c
c  ...force type
      real*8 :: z      
c      
c-----------------------------------------------------------------------
c      
      character(len=4) :: type
      type = 'mdld'
c              
      iprint=1
c
c  ...loop through orders of approximation
      do np=1,8
        Norder(1:14)=np; Norder(9)=np*10+np
        if (iprint.eq.1) then
          write(*,7001) Type,Norder
 7001     format('test_set_3Dint: Type, Norder = ',a4,2x,19i4)
        endif
        call set_3Dint(type,norder, nint,xiloc,waloc)
c
c  .....loop through monomials
        do iz=0,2*np
        do iy=0,2*np
        do ix=0,2*np
          if (ix+iy+iz.gt.2*np) cycle
          s=0.d0
          do l=1,nint
            x=xiloc(1,l) ; y=xiloc(2,l) ; z=xiloc(3,l)
            wa=waloc(l)
            s = s + (x**ix*y**iy*(1.d0-z)**iz)*wa
          enddo
          sexact = 1.d0/((ix+1.d0)*(iy+1.d0)*(ix+iy+iz+3.d0))
          if (iprint.eq.1) then
            write(*,7002) ix,iy,iz,s,sexact
 7002       format('test_set_3Dint: ix,iy,iz,s,sexact = ',
     .              3i2,2x,2e25.15)
          endif
          if (abs(s-sexact).gt.1.d-14) call pause
        enddo
        enddo
        enddo
      enddo
c
c
      end
