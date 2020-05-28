!
! system interfaces for RPN_MPI functions
!
interface
  function get_host_id() result(h) bind(C,name='gethostid')
    import :: C_LONG
    integer(C_LONG) :: h
  end function get_host_id
  function numa_node(cpu) result(n) bind(C,name='numa_node_of_cpu')
    import :: C_INT
    integer(C_INT), intent(IN), value :: cpu
    integer(C_INT) :: n
  end function 
  function sched_get_my_cpu() result(c) bind(C,name='sched_getcpu')
    import :: C_INT
    integer(C_INT) :: c
  end function sched_get_my_cpu
  function c_alarm(seconds) result(seconds_since) BIND(C,name='alarm')
    use ISO_C_BINDING
    implicit none
    integer(C_INT), intent(IN), value :: seconds
    integer(C_INT) :: seconds_since
  end function c_alarm
end interface
