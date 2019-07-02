module poorman
  implicit none
  
  

  private :: partition_real, partition_real_index
  public  :: qsort_real, qsort_real_index



  interface qsort
     module procedure qsort_real, qsort_real_index
  end interface

contains

end module quick_sort
