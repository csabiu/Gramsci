program test_pq
  use kdtree2_precision_module
  use kdtree2_priority_queue_module
  implicit none
  type(kdtree2_result), target, allocatable :: results(:)
  type(pq), target :: q_storage
  type(pq), pointer :: q
  type(kdtree2_result) :: e
  real(kdkind) :: pri

  allocate(results(10))
  q_storage = pq_create(results)
  q => q_storage

  call assert(q%heap_size == 0, "initial heap size")
  pri = pq_insert(q, 5.0_kdkind, 1)
  call assert(q%heap_size == 1, "heap size after first insert")
  call assert(abs(pri-5.0_kdkind) < 1.0e-12_kdkind, "max priority after first insert")

  pri = pq_insert(q, 10.0_kdkind, 2)
  call assert(q%heap_size == 2, "heap size after second insert")
  call assert(abs(pri-10.0_kdkind) < 1.0e-12_kdkind, "max priority after second insert")

  call pq_extract_max(q, e)
  call assert(e%idx == 2, "extract max returned wrong idx")
  call assert(q%heap_size == 1, "heap size after extract")

  call pq_extract_max(q, e)
  call assert(e%idx == 1, "extract max returned wrong idx second time")
  call assert(q%heap_size == 0, "heap size after second extract")

  print *, "All tests passed"
contains
  subroutine assert(cond, msg)
    logical, intent(in) :: cond
    character(*), intent(in) :: msg
    if (.not. cond) then
       print *, 'Test failed:', msg
       stop 1
    end if
  end subroutine assert
end program test_pq
