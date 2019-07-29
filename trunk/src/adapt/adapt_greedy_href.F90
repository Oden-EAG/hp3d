#if DEBUG_MODE

!-----------------------------------------------------------------------
!> Purpose : routine employs the greedy strategy to performe adaptive
!!           h-refinements
!!
!> @param[in] Iphy  -
!> @param[in] Imode -
!> @param[in] Eps   -
!> @param[in] Ath   -
!!
!> rev@Dec 12
!-----------------------------------------------------------------------
subroutine adapt_greedy_href(Iphy,Imode,Eps,Ath)
      use data_structure3D
!
      implicit none
      integer, intent(in) :: Iphy, Imode
      real*8 , intent(in) :: Eps, Ath
!
      real*8 :: eta, err, vol_mdle, vol_old, vol_new
      integer, allocatable :: nlist(:)
      real*8 , allocatable :: rlist(:)
!       
      integer :: iprint, mdle, iel, i, kref, istat, j, nreles_old
!-----------------------------------------------------------------------
!
      iprint=1
    
      ! list of elements to be refined  
      allocate(nlist(NRELES),rlist(NRELES), STAT=istat)
      if (istat.ne.0) then
         write(*,*)'adapt_greedy: nlist, rlist not allocated!'
         stop 1
      endif
    
      ! loop over active elements
      mdle=0 ; vol_old=0.d0
      do iel=1,NRELES
         call nelcon(mdle, mdle)
         call volume_hp_mdle(mdle, vol_mdle)
         vol_old = vol_old + vol_mdle
    
         nlist(iel) = mdle
         
         rlist(iel) = sqrt(sum(NODES(mdle)%error(0:3,Iphy))) 
         err = err + rlist(iel)
      enddo
    
      ! check absolute criterion
      if (err.gt.Ath) then
         ! sort list in descending order
         call sort(nlist,rlist,NRELES)
         write(*,*) 'adapt_greedy: min, max = ', rlist(1), rlist(NRELES)
         
         eta = (1.d0 - Eps)*rlist(1)
         nreles_old=NRELES
         
         ! refine elements from the list   
         j=0
         do i=1,nreles_old
            mdle = nlist(i)
            err  = rlist(i)
            if (err.lt.eta) exit
            if (NODES(mdle)%ref_kind.ne.0  ) cycle
            !
            j=j+1
            
            !  ...refine element
            select case (Imode) 
            case (0); call get_isoref(mdle, kref)
            case (1); call get_anisoref(mdle, NODES(mdle)%error(1:3,Iphy), kref)
            end select
            
            call refine(mdle, kref)
    
         enddo
         call close
         call update_gdof
         call update_ddof
    
         call volume_hp(vol_new)
     
     write(*,7004) j,float(j)/float(nreles_old)
7004 format(' --- number of refined elements,ratio = ',i7,2x,f5.3 )
     write(*,7005)vol_old,vol_new
7005 format(' --- voume old,new = ',2(e12.5,2x))
  else
     write(*,7006) err, Ath
7006 format(' --- err, abs threshold = ',2(e12.5,2x))
  end if

  !   deallocate lists
  deallocate(nlist,rlist, STAT=istat)
  if (istat.ne.0) then
     write(*,*)'adapt_geom: nlist, rlist not deallocated!'
     stop
  endif


endsubroutine adapt_greedy_href

#endif
