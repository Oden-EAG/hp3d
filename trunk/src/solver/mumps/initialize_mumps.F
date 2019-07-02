       subroutine initialize_mumps
c
       use mumps_data
       use parallelism
c
       call start_mumps
       call alloc_mumps
c
       return
       end
