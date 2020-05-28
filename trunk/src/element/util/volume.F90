!-----------------------------------------------------------------------
subroutine volume(Vol)
!
      use data_structure3D
      use element_data
!
      implicit none
!
!  ...dummy arguments
      real(8), intent(out) :: Vol
!
!  ...local variables
      integer, dimension(19) :: norder
      character(len=4) :: type
      real(8) :: xiloc(3,MAX_NINT3),wxi(MAX_NINT3)
      real(8) :: xi(3), x(3)
      real(8) :: dxdxi(3,3), dxidx(3,3)
      real(8) :: wa, weight, rjac
      integer :: iprint, mdle, nint, l, iflag, i
!-----------------------------------------------------------------------
!
      iprint=0
!
!  ...initialize
      Vol=0.d0
!
!  ...use 5th order to perform integration
      norder(1:19)=5
!
!  ...loop over initial mesh elements
      do mdle=1,NRELIS
!
!  .....integration points
        type = NODES(mdle)%type
        call set_3Dint(type,norder, nint,xiloc,wxi)
!
!  .....loop over integration points
        do l=1,nint
          xi(1:3)=xiloc(1:3,l) ; wa=wxi(l)
!
!  .......evaluate exact geometry map
          call exact_geom(mdle,xi, x,dxdxi)
!
!  .......evaluate the inverse derivatives and jacobian
          call geom(dxdxi, dxidx,rjac,iflag)
          if (iflag.ne.0) then
            write(*,*)'volume: NEGATIVE JACOBIAN Mdle = ',Mdle
            write(*,*) '        rjac = ',rjac
            call pause
          endif
!
!  ........total weight
           weight = wa*rjac
!
!  ........accumulate
           Vol = Vol + weight
!
!  ......end of loop over integration points
         enddo
!
!  ...end of loop over initial mesh elements
      enddo
!
!
endsubroutine volume
