!>@brief This module replaces common file "wvglob.blk"
!>@date  Feb 2024
module wvglob

   integer :: IOPEN,LWTYPE(10),LWSIZE(10,4),LSCALE(10), &
              IGPSP(10),ILEPSP,ILESCR,ILEDF,ICURRWIN
   real(8) :: XGSCAL(10,4)

!  explanation of variables:
!
!       IOPEN           .ne.0 if windows initialized at all
!       LSCALE          0/1 for scaling in windows
!       XGSCAL(4)       - xscale,yscale,xoff,yoff
!                       x <- Xworld*xscale+xoff
!       IGPSP           0/1 for postscript echo
!       ILEPSP          no of pspactive windows
!       ILESCR          no of screen windows
!       ILEDF           no of dump file windows

end module wvglob
