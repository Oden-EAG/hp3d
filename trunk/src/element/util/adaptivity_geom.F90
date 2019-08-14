#if DEBUG_MODE

!------------------------------------------------------------------------
!> Purpose : performs adaptive h-refinements of the geometry
!
!> rev@Jan 13
!------------------------------------------------------------------------
subroutine adaptivity_geom(Eps, Nref,Ratio)
!
      use data_structure3D , only : NRELES , NODES
      implicit none
! 
      real*8, intent(in ) :: Eps
      integer,intent(out) :: Nref
      real*8, intent(out) :: Ratio
! 
      integer, allocatable :: nlist(:)
      real*8 , allocatable :: rlist(:)
!
      real*8 :: eta,vol
      integer :: iprint,mdle,i,kref,istat,j,nreles_save
!      
      integer,                save :: ivis=0
      real*8, dimension(20,2),save :: rwork
      integer,dimension(20,1),save :: iwork
! 
!-----------------------------------------------------------------------
!
      iprint=0
!  
!  ...list of elements to be refined  
      allocate(nlist(NRELES), STAT=istat)
      if (istat.ne.0) then
         write(*,*)'adaptivity_geom: nlist not allocated!'
         stop
      endif
      allocate(rlist(NRELES), STAT=istat)
      if (istat.ne.0) then
         write(*,*)'adaptivity_geom: rlist not allocated!'
         stop
      endif
!  
!  ...loop over active elements and record error
      do i=1,NRELES
         mdle = ELEM_ORDER(i)
         nlist(i)=mdle ; rlist(i)=NODES(mdle)%error(0,0)
      enddo
!  
!  ...sort list in descending order
      call sort(nlist,rlist,NRELES)
      eta=Eps*rlist(1) ; nreles_save=NRELES
!    
!  ...refine elements from the list   
      j=0
      do i=1,nreles_save
         if (rlist(i)                .lt.eta) exit
         if (NODES(nlist(i))%ref_kind.ne.0  ) cycle
         j=j+1
!      
!  ......printing
         if (iprint.eq.1) then
            write(*,7003)j,nlist(i),NODES(nlist(i))%Type,rlist(i)
 7003    format(' j,mdle,type,err = ',i4,2x,i5,2x,a4,2x,e12.5)
         endif
!  
!  ......refine element
         call get_isoref(nlist(i), kref)
         call refine(    nlist(i), kref)
      enddo
!      
      call close
      call update_gdof
!  
!  ...number of refined elements, ratio
      Nref=j ; Ratio=float(j)/float(nreles_save)
      call volume_hp(vol)
!  
!  ...report
      ivis=ivis+1
      rwork(ivis,1)=vol
      rwork(ivis,2)=100.d0*Ratio
      iwork(ivis,1)=Nref
!
      write(*,*)'-- Geometry Adapt. Report --'
      do i=1,ivis
        write(*,9999)i,iwork(i,1),rwork(i,2),rwork(i,1)
 9999   format(' i,n_ref,perc_ref,volume = ',i2,2x,i4,2x,f6.2,2x,e12.5)       
      enddo
      write(*,*)''

!  ...deallocate lists
      deallocate(nlist, STAT=istat)
      if (istat.ne.0) then
         write(*,*)'adaptivity_geom: nlist not deallocated!'
         stop
      endif
      deallocate(rlist, STAT=istat)
      if (istat.ne.0) then
         write(*,*)'adaptivity_geom: rlist not deallocated!'
         stop
      endif
!
!
endsubroutine adaptivity_geom

#endif
