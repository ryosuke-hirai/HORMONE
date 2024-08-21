module mpi_utils
#ifdef MPI
 use mpi
#endif
  implicit none

  interface allreduce_mpi
    module procedure allreduce_mpi_real8
    module procedure allreduce_mpi_real8_array
    module procedure allreduce_mpi_int4
    module procedure allreduce_mpi_int4_array
  end interface allreduce_mpi

  public :: allreduce_mpi

  integer :: myrank, nprocs
  integer :: mpi_subarray_default, mpi_subarray_spc, mpi_subarray_gravity, mpi_type_sink_prop, mpi_type_sink_prop_array, mpi_subarray_extgrtv

  contains

  subroutine init_mpi()
#ifdef MPI
    integer :: ierr
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
#else
    myrank = 0
    nprocs = 1
#endif
  end subroutine init_mpi

  subroutine finalize_mpi()
#ifdef MPI
    integer :: ierr
    call MPI_FINALIZE(ierr)
#endif
  end subroutine finalize_mpi

  subroutine allreduce_mpi_real8(op_str, value)
    character(len=*), intent(in) :: op_str
    real(8), intent(inout) :: value

#ifdef MPI
    real(8) :: value_reduced
    integer :: op, ierr

    select case (op_str)
    case ('sum')
      op = MPI_SUM
    case ('max')
      op = MPI_MAX
    case ('min')
      op = MPI_MIN
    case default
      error stop "Invalid operation specified"
    end select

    call MPI_ALLREDUCE(value, value_reduced, 1, MPI_REAL8, op, MPI_COMM_WORLD, ierr)
    value = value_reduced
#else
    ! Check if operation is valid, but do nothing
    select case (op_str)
    case ('sum', 'max', 'min')
      value = value
    case default
      error stop "Invalid operation specified"
    end select
#endif

  end subroutine allreduce_mpi_real8

  subroutine allreduce_mpi_real8_array(op_str, value)
    character(len=*), intent(in) :: op_str
    real(8), intent(inout) :: value(:)

#ifdef MPI
    real(8), allocatable :: value_reduced(:)
    integer :: op, ierr

    allocate(value_reduced,mold=value)

    select case (op_str)
    case ('sum')
      op = MPI_SUM
    case ('max')
      op = MPI_MAX
    case ('min')
      op = MPI_MIN
    case default
      error stop "Invalid operation specified"
    end select

    call MPI_ALLREDUCE(value, value_reduced, size(value), MPI_REAL8, op, MPI_COMM_WORLD, ierr)
    value = value_reduced

    deallocate(value_reduced)
#else
    ! Check if operation is valid, but do nothing
    select case (op_str)
    case ('sum', 'max', 'min')
      value = value
    case default
      error stop "Invalid operation specified"
    end select
#endif

  end subroutine allreduce_mpi_real8_array

  subroutine allreduce_mpi_int4(op_str, value)
    character(len=*), intent(in) :: op_str
    integer, intent(inout) :: value

#ifdef MPI
    integer :: value_reduced
    integer :: op, ierr

    select case (op_str)
    case ('sum')
      op = MPI_SUM
    case ('max')
      op = MPI_MAX
    case ('min')
      op = MPI_MIN
    case default
      error stop "Invalid operation specified"
    end select

    call MPI_ALLREDUCE(value, value_reduced, 1, MPI_INTEGER, op, MPI_COMM_WORLD, ierr)
    value = value_reduced
#else
    ! Check if operation is valid, but do nothing
    select case (op_str)
    case ('sum', 'max', 'min')
      value = value
    case default
      error stop "Invalid operation specified"
    end select
#endif

  end subroutine allreduce_mpi_int4

  subroutine allreduce_mpi_int4_array(op_str, value)
    character(len=*), intent(in) :: op_str
    integer, intent(inout) :: value(:)

#ifdef MPI
    integer, allocatable :: value_reduced(:)
    integer :: op, ierr

    allocate(value_reduced,mold=value)

    select case (op_str)
    case ('sum')
      op = MPI_SUM
    case ('max')
      op = MPI_MAX
    case ('min')
      op = MPI_MIN
    case default
      error stop "Invalid operation specified"
    end select

    call MPI_ALLREDUCE(value, value_reduced, size(value), MPI_INTEGER, op, MPI_COMM_WORLD, ierr)
    value = value_reduced

    deallocate(value_reduced)
#else
    ! Check if operation is valid, but do nothing
    select case (op_str)
    case ('sum', 'max', 'min')
      value = value
    case default
      error stop "Invalid operation specified"
    end select
#endif

  end subroutine allreduce_mpi_int4_array

  subroutine barrier_mpi
#ifdef MPI
    integer :: ierr
    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
  end subroutine

  subroutine stop_mpi(code)
    integer, intent(in) :: code
#ifdef MPI
    call finalize_mpi
#endif
    if (myrank==0 .and. code/=0) then
      error stop code
    endif

    stop
  end subroutine

end module mpi_utils
