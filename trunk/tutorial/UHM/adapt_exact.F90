#if DEBUG_MODE
!
!------------------------------------------------------------------------
!> Purpose : performs adaptive h-refinements of the geometry
!------------------------------------------------------------------------
subroutine adapt_geom
    !
    use data_structure3D
    implicit none
    ! 
    !   ...local variables
    real*8 :: eta, vol_old, vol_new
    ! 
    integer, allocatable :: nlist(:)
    real*8 , allocatable :: rlist(:)
    !
    integer :: iprint, mdle, iel, i, kref, istat, j, nreles_save
    ! 
    !   ...percentage of elements to be refined  
    real*8 , parameter :: eps    = 0.5d0
!-----------------------------------------------------------------------
    iprint=0
!  
!   list of elements to be refined  
    allocate(nlist(NRELES), STAT=istat)
    if (istat.ne.0) then
      write(*,*)'adapt_geom: nlist not allocated!'
      stop
    endif        
    allocate(rlist(NRELES), STAT=istat)
    if (istat.ne.0) then
      write(*,*)'adapt_geom: rlist not allocated!'
      stop
    endif 
!  
!   loop over active elements
    mdle=0
    do iel=1,NRELES
      call nelcon(mdle, mdle)
      nlist(iel) = mdle ; rlist(iel) = NODES(mdle)%error
    enddo
!  
!   compute volume
    call volume_hp(vol_old)
!  
!   sort list in descending order
    call sort(nlist,rlist,NRELES)
    eta = eps*rlist(1)
    nreles_save=NRELES
!    
!   refine elements from the list   
    do i=1,nreles_save
      if (rlist(i)                .lt.eta) exit
      if (NODES(nlist(i))%ref_kind.ne.0  ) cycle
      j=j+1
!      
!  ...printing
      if (iprint.eq.1) then
        write(*,7003)j,nlist(i),NODES(nlist(i))%Type,rlist(i)
7003    format(' j,mdle,type,err = ',i4,2x,i5,2x,a4,2x,e12.5)
      endif
!  
!  ...refine element
      call get_isoref(nlist(i), kref)
      call refine(    nlist(i), kref)
      call close
    enddo
    call update_gdof
    call volume_hp(vol_new)
!  
    write(*,7004)j,float(j)/float(nreles_save)
7004 format(' number of processed elements,ratio = ',i7,2x,f5.3 )
    write(*,7005)vol_old,vol_new
7005 format(' voume old,new = ',2(e12.5,2x))
!  
!   deallocate lists
    deallocate(nlist, STAT=istat)
    if (istat.ne.0) then
      write(*,*)'adapt_geom: nlist not deallocated!'
      stop
    endif        
    deallocate(rlist, STAT=istat)
    if (istat.ne.0) then
      write(*,*)'adapt_geom: rlist not deallocated!'
      stop
    endif 
!  
!
end subroutine adapt_geom
!
!
!
!
!------------------------------------------------------------------------
!> Purpose : performs adaptive h-refinements of the geometry
!------------------------------------------------------------------------
subroutine adapt_h_refinement_back
    use data_structure3D
    implicit none
    ! 
    !   ...local variables
    real*8 :: Err, Rnorm , derr, dnorm
    ! 
    integer, allocatable :: nlist(:)
    real*8 , allocatable :: rlist(:)
    !
    real*8  :: eta ,err_rate, vol
    integer :: iprint, mdle, iel, ierr, i, kref, istat, j, nreles_save
    ! 
    !   ...percentage of elements to be refined  
    real*8 , parameter :: eps    = 0.5d0
    ! 
    !   ...variables saved between subsequent calls 
    integer, save :: ivis        = 0
    integer, save :: nrgdof_save
    real*8 , save ::    err_save,vol_save
!-----------------------------------------------------------------------
  iprint=0
!
! infinite loop
  do
!  
!   list of elements to be refined  
    allocate(nlist(NRELES), STAT=istat)
    if (istat.ne.0) then
      write(*,*)'adapt_h_refinement: nlist not allocated!'
      stop
    endif        
    allocate(rlist(NRELES), STAT=istat)
    if (istat.ne.0) then
      write(*,*)'adapt_h_refinement: rlist not allocated!'
      stop
    endif 
!    
!   initialize global quantities  
    Err = 0.d0 ; Rnorm = 0.d0
!  
    write(*,*)'L2 error only? 1-No, 0-Yes' 
    read(*,*)ierr
!!!    ierr=1
!  
!   loop over active elements
    mdle=0
    do iel=1,NRELES
      call nelcon(mdle, mdle)
      call geometry_error_element(mdle,ierr, derr,dnorm)
      nlist(iel) = mdle ; rlist(iel) = derr
!  
!     accumulate        
      Err   = Err   + derr
      Rnorm = Rnorm + dnorm
    enddo
!  
    Err   = sqrt(Err)
    Rnorm = sqrt(Rnorm)
!  
!   compute volume
    call volume_hp(vol)
!      
!   if not first visit, compute rate
    if (ivis.gt.0) then
      write(*,7002)nrgdof_save,NRGDOF
7002  format(' adapt_h_refinement: NRGDOF old,new = ',2(i7,2x))
!
      write(*,7000)vol_save,vol
7000  format(' adapt_h_refinement: volume old,new = ',2(e12.5,2x))
!
      err_rate = log(err_save/Err)/log(float(nrgdof_save)/NRGDOF)
      write(*,7006)Err,Rnorm,Err/Rnorm,err_rate
7006  format(' adapt_h_refinement: Err,Rnorm,Err/Rnorm,Rate = ',3(e12.5,2x),f9.6)   
!   else rate is meaningless
    else
      write(*,7005)Err,Rnorm,Err/Rnorm
7005  format(' adapt_h_refinement: Err,Rnorm,Err/Rnorm      = ',3(e12.5,2x)     )   
    endif
!
!   save error and number of gdofs
    err_save=Err ; nrgdof_save=NRGDOF ; vol_save=vol
!
    write(*,*)'adaptive h-refinements? 1-Yes, 0-No' 
    read(*,*)ierr
!!!    ierr=1
    if (ierr.eq.0) return
!  
!   sort list in descending order
    call sort(nlist,rlist,NRELES)
    eta = eps*rlist(1)
!  
!   save number of active elements, and initialize counter of refined elements
    nreles_save=NRELES ; j=0
!    
!   refine elements from the list   
    do i=1,nreles_save
      if (rlist(i)                .lt.eta) exit
      if (NODES(nlist(i))%ref_kind.ne.0  ) cycle
      j=j+1
!      
!   ..printing
      if (iprint.eq.1) then
        write(*,7003)j,nlist(i),NODES(nlist(i))%Type,rlist(i)
7003    format(' j,mdle,type,err = ',i4,2x,i5,2x,a4,2x,e12.5)
      endif
!  
!   ..refine element
      call get_isoref(nlist(i), kref)
      call refine(    nlist(i), kref)
      call close
    enddo
!  
    write(*,7004)j,float(j)/float(nreles_save)
7004 format(' number of processed elements,ratio = ',i7,2x,f5.3 )
!    
    call update_gdof
!  
!   deallocate lists
    deallocate(nlist, STAT=istat)
    if (istat.ne.0) then
      write(*,*)'adapt_h_refinement: nlist not deallocated!'
      stop
    endif        
    deallocate(rlist, STAT=istat)
    if (istat.ne.0) then
      write(*,*)'adapt_h_refinement: rlist not deallocated!'
      stop
    endif 
!  
!   raise visitation flag
    ivis=1
!    
  enddo
!
!
end subroutine adapt_h_refinement_back

#endif
