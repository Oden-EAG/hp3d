c-----------------------------------------------------------------------
c   routine name       - customize_geometry_toy_sphere
c-----------------------------------------------------------------------
c
c   latest revision    - Mar 09
c
c   purpose            - routine customizes geometry for toy_sphere
c                        problem
c
c-----------------------------------------------------------------------
      subroutine customize_geometry_toy_sphere
c
      use control
      use GMP
      use AF
#include "syscom.blk"
c
c  ...layer thickness
      dimension thickness(10), no_dom(10), xp(3)
c
c  ...set number of truncating sphere
      nts = 1
c
c  ...set the number of layers to add and their thickness
      n_lay = 1
      thickness(1:n_lay) = (/ 0.2d0 /)
c
c  ...set subdomain numbers....
      do i = 1, n_lay
        no_dom(i) = 2
      enddo
c
c  ...rescale sphere radius
      r = 0.8d0
      SURFACES(nts)%Rdata(4) = r
      do i = 1, NRPOINT
        POINTS(i)%Rdata(1:3) = r*POINTS(i)%Rdata(1:3)
      end do
c
c  ...add layers
      write(*,1000)
 1000 format(' === customize_geometry_toy_sphere  ====================')
 1007 format(' =======================================================')
      write(*,1008)r
 1008 format(' sphere radius set to r = ', F4.2)
      write(*,1001)
 1001 format(' ADDING PRISMATIC LAYER')
      write(*,1002) nts
 1002 format(' - truncating surface = ', I2)
      write(*,1003) n_lay
 1003 format(' - number of layers   = ', I2)
      write(*,1004) thickness(1:n_lay)
 1004 format(' - layers thickness   = ', 10F4.2)
      write(*,1005) no_dom(1:n_lay)
 1005 format(' - domains number     = ', 5I2)
c
      call add_prism_layer(nts,n_lay,thickness,no_dom)
      write(*,1006)
 1006 format(' PRISMATIC LAYER ADDED')
      write(*,1007)
c
      end subroutine customize_geometry_toy_sphere
